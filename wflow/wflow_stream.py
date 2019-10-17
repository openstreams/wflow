#!/usr/bin/python

"""
Run the wflow_stream hydrological model..


usage

::


    wflow_stream [-C casename][-R runId][-c configfile]
                [-T last_step][-S first_step][-s seconds]
                [-p inputparameter multiplication][-l loglevel]

    -C: set the name  of the case (directory) to run

    -R: set the name runId within the current case

    -c name of the config file (in the case directory)

    -T: Set end time of the run: yyyy-mm-dd hh:mm:ss

    -S: Set start time of the run: yyyy-mm-dd hh:mm:ss

    -s: Set the model timesteps in seconds

    -P: set parameter change string (e.g: -P "self.FC = self.FC * 1.6") for non-dynamic variables

    -p: set parameter change string (e.g: -P "self.Precipitation = self.Precipitation * 1.11") for
        dynamic variables

    -l: loglevel (most be one of DEBUG, WARNING, ERROR)

"""

import getopt
import pdb
import pcraster as pcr
import numpy as np

from wflow.wf_DynamicFramework import *
from wflow.wflow_adapt import *

# import scipy


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args:
        print(msg)
    print(__doc__)
    sys.exit(0)


class WflowModel(pcraster.framework.DynamicModel):
    """
  The user defined model class. This is your work!
  """

    def __init__(self, cloneMap, Dir, RunDir, configfile):
        """
      *Required*
      
      The init function **must** contain what is shown below. Other functionality
      may be added by you if needed.
      
      """
        pcraster.framework.DynamicModel.__init__(self)

        self.caseName = os.path.abspath(Dir)
        self.clonemappath = os.path.join(os.path.abspath(Dir), "staticmaps", cloneMap)
        pcr.setclone(self.clonemappath)
        self.runId = RunDir
        self.Dir = os.path.abspath(Dir)
        self.configfile = configfile
        self.SaveDir = os.path.join(self.Dir, self.runId)


    def parameters(self):
        """
      List all the parameters (both static and forcing here). Use the wf_updateparameters()
      function to update them in the initial section (static) and the dynamic section for
      dynamic parameters and forcing date.

      Possible parameter types are:

      + staticmap: Read at startup from map
      + statictbl: Read at startup from tbl, fallback to map (need Landuse, Soil and TopoId (subcatch) maps!
      + timeseries: read map for each timestep
      + monthlyclim: read a map corresponding to the current month (12 maps in total)
      + dailyclim: read a map corresponding to the current day of the year
      + hourlyclim: read a map corresponding to the current hour of the day (24 in total)


      :return: List of modelparameters
      """
        modelparameters = []

        # creating mapstacks
        self.P_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "precipitation", "/inmaps/Precip"
        )  # timeseries for rainfall
        self.PET_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "PET", "/inmaps/PM"
        )  # timeseries for potential evapotranspiration
        self.TEMP_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "temperature", "/inmaps/temperature"
        )  # timeseries for air temperature
        
        # Meteo and other forcing
        modelparameters.append(
            self.ParamType(
                name="precipitation",
                stack=self.P_mapstack,
                type="timeseries",
                default=0.0,
                verbose=True,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="PET",
                stack=self.PET_mapstack,
                type="timeseries",
                default=0.0,
                verbose=True,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="temperature",
                stack=self.TEMP_mapstack,
                type="timeseries",
                default=10.0,
                verbose=True,
                lookupmaps=[],
            )
        )

        # Static model parameters
        modelparameters.append(
            self.ParamType(
                name="C",
                stack="staticmaps/C.map",
                type="staticmap",
                default=0.0,
                verbose=False,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="cropf",
                stack="staticmaps/cropf.map",
                type="staticmap",
                default=1.0,
                verbose=False,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="whc",
                stack="staticmaps/whc.map",
                type="staticmap",
                default=100.0,
                verbose=False,
                lookupmaps=[],
            )
        )

        modelparameters.append(
            self.ParamType(
                name="meltf",
                stack="staticmaps/meltf.map",
                type="staticmap",
                default=0.3,
                verbose=False,
                lookupmaps=[],
            )
        )
        
        modelparameters.append(
            self.ParamType(
                name="togwf",
                stack="staticmaps/togwf.map",
                type="staticmap",
                default=0.4,
                verbose=False,
                lookupmaps=[],
            )
        )
        
        self.wf_multparameters()

        return modelparameters

    def stateVariables(self):

        states = [  "Ground_water",
                    "Available_water",
                    "snow",
                    ]

        return states

    def supplyCurrentTime(self):
        """
      *Optional*
      
      Supplies the current time in seconds after the start of the run
      This function is optional. If it is not set the framework assumes
      the model runs with daily timesteps.
      
      Output:
      
          - time in seconds since the start of the model run
          
      """

        return self.currentTimeStep() * int(
            configget(self.config, "model", "timestepsecs", "86400")
        )

    def suspend(self):
        """
      *Required*
      
      Suspends the model to disk. All variables needed to restart the model
      are saved to disk as pcraster maps. Use resume() to re-read them
      
      This function is required. 
      
    """

        self.logger.info("Saving initial conditions...")
        #: It is advised to use the wf_suspend() function
        #: here which will suspend the variables that are given by stateVariables
        #: function.
        self.wf_suspend(self.SaveDir + "/outstate/")

    def initial(self):

        """
    *Required*
    
    Initial part of the model, executed only once. It reads all static model
    information (parameters) and sets-up the variables used in modelling.
    
    This function is required. The contents is free. However, in order to
    easily connect to other models it is advised to adhere to the directory
    structure used in the other models.
    
    """
        #: pcraster option to calculate with units or cells. Not really an issue
        #: in this model but always good to keep in mind.
        pcr.setglobaloption("unittrue")

        self.timestepsecs = int(
            configget(self.config, "model", "timestepsecs", "86400")
        )
        self.basetimestep = 86400
        # Reads all parameter from disk
        self.wf_updateparameters()
        self.wf_multparameters()    # needed so parameters can be altered in the inifile
        self.logger.info("Starting Dynamic run...")

    def resume(self):
        """ 
    *Required*

    This function is required. Read initial state maps (they are output of a 
    previous call to suspend()). The implementation shown here is the most basic
    setup needed.
    
    """
        self.logger.info("Reading initial conditions...")
        #: It is advised to use the wf_resume() function
        #: here which pick up the variable save by a call to wf_suspend()

        if self.reinit:
            self.logger.warning("Setting initial states to default")
            for s in self.stateVariables():
                exec("self." + s + " = cover(1.0)")
        else:
            try:
                self.wf_resume(self.Dir + "/instate/")
            except:
                self.logger.warning("Cannot load initial states, setting to default")
                for s in self.stateVariables():
                    exec("self." + s + " = cover(1.0)")

    def default_summarymaps(self):
        """
      *Optional*

      Return a default list of variables to report as summary maps in the outsum dir.
      The ini file has more options, including average and sum
      """
        return ["self.precipitation"]

    def dynamic(self):
        """
      *Required*
      
      This is where all the time dependent functions are executed. Time dependent
      output should also be saved here.
      """
        self.wf_updateparameters()  # read the temperature map for each step (see parameters())

        self.wf_multparameters()    # needed so parameters can be altered in the inifile under [variable_change_once]

        # snow routine
        snowfall = self.precipitation
        snowfall = pcr.scalar(self.temperature < 3.0) * snowfall
        self.snow = self.snow + snowfall
        melt = (self.meltf * self.temperature)
        melt = pcr.scalar(self.temperature > 3.0) * melt
        melt = pcr.max(0.0, pcr.min(self.snow, melt))
        self.snow = self.snow - melt
        self.precipitation = self.precipitation - snowfall + melt

        # soil storage
        self.PotenEvap = self.PET * self.cropf
        Peff = self.precipitation - self.PotenEvap
        Aw_1 = self.Available_water

        # Soil wetting below capacity
        self.whc = pcr.max(0.0, self.whc)
        self.below_cap = (Aw_1 + pcr.max(0.0, Peff)) <= self.whc
        self.available_water_below_cap = (self.Available_water * pcr.scalar(self.below_cap)) + (pcr.max(0.0, Peff) * pcr.scalar(self.below_cap))

        # Soil wetting above capacity
        self.above_cap = (Aw_1 + Peff) > self.whc
        self.excess = (Aw_1 * pcr.scalar(self.above_cap)) + (Peff * pcr.scalar(self.above_cap)) - (self.whc * pcr.scalar(self.above_cap))
        self.available_water_above_cap = (self.whc * pcr.scalar(self.above_cap))
        self.Available_water = self.available_water_below_cap + self.available_water_above_cap
        
        
        # soil drying
        self.drying = Peff <= 0.0
        self.exp_content = (Peff)/(self.whc)
        self.available_water_drying = (Aw_1 * pcr.scalar(self.drying)) * pcr.exp(self.exp_content)

        # soil that is not drying
        self.not_drying = Peff > 0.0
        self.Available_water = self.Available_water * pcr.scalar(self.not_drying)
        self.Available_water = self.Available_water + self.available_water_drying
        self.excess = self.excess * pcr.scalar(self.not_drying)

        # seepage to groundwater
        self.runoff = self.togwf * self.excess
        self.togw = self.excess - self.runoff

        # adding water to groundwater and taking water from groundwater
        self.Ground_water = self.Ground_water + self.togw
        self.sloflo = (self.Ground_water/self.C)
        self.Ground_water = self.Ground_water - self.sloflo

        # adding water from groundwater to runoff
        self.runoff = self.runoff + self.sloflo

# The main function is used to run the program from the command line

def main(argv=None):
    """
    *Optional but needed it you want to run the model from the command line*
    
    Perform command line execution of the model. This example uses the getopt
    module to parse the command line options.
    
    The user can set the caseName, the runDir, the timestep and the configfile.
    """
    global multpars
    caseName = "default"
    runId = "run_default"
    configfile = "wflow_stream.ini"
    _lastTimeStep = 0
    _firstTimeStep = 0
    timestepsecs = 86400
    wflow_cloneMap = "whc.map"

    # This allows us to use the model both on the command line and to call
    # the model usinge main function from another python script.

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return

    opts, args = getopt.getopt(argv, "C:S:T:c:s:R:l:P:p")

    for o, a in opts:
        if o == "-C":
            caseName = a
        if o == "-R":
            runId = a
        if o == "-c":
            configfile = a
        if o == "-s":
            timestepsecs = int(a)
        if o == "-T":
            _lastTimeStep = int(a)
        if o == "-S":
            _firstTimeStep = int(a)
        if o == "-l":
            exec("loglevel = logging." + a)
        if o == "-P":
            left = a.split("=")[0]
            right = a.split("=")[1]
            configset(
                myModel.config, "variable_change_once", left, right, overwrite=True
            )
        if o == "-p":
            left = a.split("=")[0]
            right = a.split("=")[1]
            configset(
                myModel.config, "variable_change_timestep", left, right, overwrite=True
            )

    if len(opts) < 1:
        usage()

    myModel = WflowModel(wflow_cloneMap, caseName, runId, configfile)
    dynModelFw = wf_DynamicFramework(
        myModel, _lastTimeStep, firstTimestep=_firstTimeStep
    )
    dynModelFw.createRunId(NoOverWrite=False, level=logging.DEBUG)
    dynModelFw._runInitial()
    dynModelFw._runResume()
    # dynModelFw._runDynamic(0,0)
    dynModelFw._runDynamic(_firstTimeStep, _lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()

if __name__ == "__main__":
    main()
