#!/usr/bin/python

"""
Definition of the wflow_floodmap model.
---------------------------------------

Runs flood mapping (very basic) based on an existing
wflow_sbm|hbv model run. 

The wflow_sbm|hbv model must have saved mapstacks for
water level and discharge for each timestep (run*****.*** and lev*****.***).
If the name of you Q and/or H maps are different specify these in the
[inputmapstacks] section, e.g:
    
::

    [inputmapstacks]
    Q = runDyn
    H = levDyn
    

If a wflow_bankfull.map map is present in the staticmaps directory that map will
be used to determine if the river is flooding, otherwise bankfull
is determined using: Bankful = RiverWidth/60 (the RiverWidth map is taken
from the runid/outsum directory)

Ini file settings:

::
    
    [model]
    # Maximum distance between a cell to be flooded cell and a river cell
    # or already flooded cel. Functions as a max flooding velocity
    # Never set lower that the length of one cell.
    maxflooddist= 0.3
    


Usage:
wflow_floodmap  -C case -R Runid -c inifile -h -I

    -C: set the name  of the case (directory) to run
    
    -I: generate initial conditions from scratch
    
    -R: set the name runId within the current case
    
    -c name of the config file (in the case directory)
    
    -h displays help information

    -l: loglevel (most be one of DEBUG, WARNING, ERROR)
    
$Author: schelle $
$Id: wflow_floodmap.py 916 2014-02-11 14:49:35Z schelle $
$Rev: 916 $
"""

# TODO: update to update framework

import os.path

import pcraster.framework
from wflow.wf_DynamicFramework import *
from wflow.wflow_adapt import *


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
      Initialize the object
      
      """
        pcraster.framework.DynamicModel.__init__(self)

        self.caseName = os.path.abspath(Dir)
        self.clonemappath = os.path.join(os.path.abspath(Dir), "staticmaps", cloneMap)
        pcr.setclone(self.clonemappath)
        self.runId = RunDir
        self.Dir = os.path.abspath(Dir)
        self.configfile = configfile
        self.SaveDir = os.path.join(self.Dir, self.runId)

    def stateVariables(self):
        """ 
      *Required*
      
      Returns a list of state variables that are essential to the model. 
      This list is essential for the resume and suspend functions to work.
      
      This function is specific for each model and **must** be present. This is
      where you specify the state variables of you model. If your model is stateless
      this function must return and empty array (states = [])
      

      
      :var FloodExtent.map: Current FloodExtent
      
      """
        states = ["FloodExtent"]

        return states

    def supplyCurrentTime(self):
        """
      *Optional*
      
      Supplies the current time in seconds after the start of the run
      This function is optional. If it is not set the framework assumes
      the model runs with daily timesteps.
      
      Ouput:
      
          - time in seconds since the start of the model run
          
      """

        return self.currentTimeStep() * int(
            configget(self.config, "model", "timestepsecs", "86400")
        )

    def suspend(self):
        """      
      Suspends the model to disk. All variables needed to restart the model
      are saved to disk as pcraster maps. Use resume() to re-read them
      
    """

        # self.logger.info("Saving initial conditions...")
        #: It is advised to use the wf_suspend() function
        #: here which will suspend the variables that are given by stateVariables
        #: function.
        self.logger.info("Saving initial conditions...")
        self.wf_suspend(os.path.join(self.SaveDir, "outstate"))

        pcr.report(
            pcr.ifthen(self.MaxDepth > 0.0, self.MaxDepth),
            os.path.join(self.SaveDir, "outsum", "MaxDepth.map"),
        )
        pcr.report(
            pcr.ifthen(pcr.scalar(self.MaxExt) > 0.0, self.MaxExt),
            os.path.join(self.SaveDir, "outsum", "MaxExt.map"),
        )

    def parameters(self):
        """
    Define all model parameters here that the framework should handle for the model
    See wf_updateparameters and the parameters section of the ini file
    If you use this make sure to all wf_updateparameters at the start of the dynamic section
    and at the start/end of the initial section
    """
        modelparameters = []

        # Static model parameters e.g.
        # modelparameters.append(self.ParamType(name="RunoffGeneratingGWPerc",stack="intbl/RunoffGeneratingGWPerc.tbl",type="static",default=0.1))
        # 3: Input time series ###################################################
        self.WL_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "WaterLevel", "/inmaps/H"
        )  # timeseries for level

        modelparameters.append(
            self.ParamType(
                name="WL",
                stack=self.WL_mapstack,
                type="timeseries",
                default=0.0,
                verbose=True,
                lookupmaps=[],
            )
        )

        return modelparameters

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

        #: Note the use of the configget functione below. This way you sepcify a default
        #: for a parameter but it can be overwritten by the uses in the ini file.
        self.timestepsecs = int(
            configget(self.config, "model", "timestepsecs", "86400")
        )
        self.maxdist = float(configget(self.config, "model", "maxflooddist", "1E31"))
        self.reinit = int(configget(self.config, "run", "reinit", "0"))

        self.wf_updateparameters()

        self.basetimestep = 86400
        self.SaveMapDir = self.Dir + "/" + self.runId + "/outmaps"
        self.Altitude = pcr.readmap(self.Dir + "/staticmaps/wflow_dem")
        self.River = pcr.readmap(self.Dir + "/staticmaps/wflow_river")
        self.Ldd = pcr.readmap(self.Dir + "/staticmaps/wflow_ldd")
        self.BankFull = pcrut.readmapSave(self.Dir + "/staticmaps/wflow_bankfull", 0.0)
        self.RiverWidth = pcr.readmap(
            self.Dir + "/" + self.runId + "/outsum/RiverWidth.map"
        )
        self.BankFull = pcr.ifthenelse(
            self.BankFull == 0.0, self.RiverWidth / 60.0, self.BankFull
        )

        self.FloodDepth = pcr.scalar(pcr.cover(0.0))
        self.MaxExt = pcr.boolean(pcr.cover(0.0))
        self.MaxDepth = pcr.cover(0.0)
        self.logger.info("End of initial...")

    def resume(self):
        """ 
    *Required*

    This function is required. Read initial state maps (they are output of a 
    previous call to suspend()). The implementation showns her is the most basic 
    setup needed.
    
    """
        # self.logger.info("Reading initial conditions...")
        #: It is advised to use the wf_resume() function
        #: here which pick upt the variable save by a call to wf_suspend()
        if self.reinit == 1:
            self.logger.info("Setting initial conditions to default (zero)")
            self.FloodExtent = pcr.cover(pcr.boolean(0))
        else:
            self.wf_resume(os.path.join(self.Dir, "instate"))

    def dynamic(self):
        """
      *Required*
      
      This is where all the time dependent functions are executed. Time dependent
      output should also be saved here.
      
      :var self.FLOOD: Actual flood level [m]
      :var self.FloodExtent: Actual flood extent [-]
      
      """
        self.FloodDepth = pcr.scalar(self.FloodDepth) * 0.0

        self.wf_updateparameters()

        self.WLatRiver = pcr.ifthenelse(
            pcr.scalar(self.River) > 0, self.WL + self.Altitude, pcr.scalar(0.0)
        )
        # WL surface if level > bankfull. For the eventual surface substract bankfull as measure for river depth
        self.water_surf = pcr.cover(
            pcr.ifthen(
                self.WLatRiver > (self.BankFull + self.Altitude),
                self.WLatRiver - self.BankFull,
            ),
            0.0,
        )
        self.water_surf_id = pcr.ordinal(pcr.uniqueid(pcr.boolean(self.water_surf)))

        # Check how many points over bankfull we have
        tmp = pcr.pcr2numpy(pcr.mapmaximum(self.water_surf_id), 0)
        fld_points_a = tmp[0, 0]

        self.logger.info(
            "Step: "
            + str(self.currentStep)
            + ". River cells over bankfull: "
            + str(fld_points_a)
        )
        # Only do mapping of the number of points is larger then 0

        if fld_points_a < 1:
            self.FloodDepth = pcr.ifthen(self.FloodDepth > 1e31, self.FloodDepth)
            self.distfromriv = self.FloodDepth
            self.spread = self.FloodDepth
            self.FloodExtent = pcr.cover(pcr.boolean(0))
            self.FloodDepth = pcr.scalar(self.FloodExtent)
        else:
            # find zones connect to a rivercell > bankfull
            self.RiverCellZones = pcr.subcatchment(self.Ldd, self.water_surf_id)
            self.spreadRivLev = pcr.areaaverage(
                pcr.ifthen(self.water_surf > 0, self.water_surf), self.RiverCellZones
            )
            self.spreadRivDemLev = pcr.areaaverage(
                pcr.ifthen(self.water_surf > 0, self.Altitude), self.RiverCellZones
            )

            # add the new first estimate to the old extent
            self.FloodExtent = pcr.cover(
                pcr.boolean(self.FloodExtent), pcr.boolean(self.water_surf_id)
            )
            # determine the distance to the nearest already flooded celll
            self.distfromriv = ldddist(
                self.Ldd, self.FloodExtent, 1
            )  # is in units of model (degree here)

            # a cell is flooded if the bottomlevel is lower than the waterlevel of the nearest river cell.
            self.FloodDepth = pcr.ifthenelse(
                self.spreadRivLev - self.Altitude >= 0.0,
                self.spreadRivLev - self.Altitude,
                0.0,
            )
            # Exclude points too far away()
            self.FloodDepth = pcr.ifthenelse(
                self.distfromriv > self.maxdist, 0.0, self.FloodDepth
            )
            self.FloodDepth = pcr.ifthen(self.FloodDepth > 0.0, self.FloodDepth)
            self.FloodExtent = pcr.ifthenelse(
                self.FloodDepth > 0.0, pcr.boolean(1), pcr.boolean(0)
            )

            # Keep track of af depth and extent
            self.MaxDepth = pcr.max(self.MaxDepth, pcr.cover(self.FloodDepth, 0))
            self.MaxExt = pcr.max(
                pcr.scalar(self.MaxExt), pcr.scalar(pcr.cover(self.FloodExtent, 0))
            )

        # reporting of maps is done by the framework (see ini file)


# The main function is used to run the program from the command line


def main(argv=None):
    """
    *Optional*
    
    Perform command line execution of the model. This example uses the getopt
    module to parse the command line options.
    
    The user can set the caseName, the runDir, the timestep and the configfile.
    """
    global multpars
    caseName = "default"
    runId = "run_default"
    configfile = "wflow_floodmap.ini"
    _lastTimeStep = 0
    _firstTimeStep = 0
    timestepsecs = 86400
    wflow_cloneMap = "wflow_subcatch.map"
    runinfoFile = "runinfo.xml"
    loglevel = logging.DEBUG

    # This allows us to use the model both on the command line and to call
    # the model usinge main function from another python script.

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return

    opts, args = getopt.getopt(argv, "C:hS:T:c:s:R:fIs:l:")

    for o, a in opts:
        if o == "-C":
            caseName = a
        if o == "-R":
            runId = a
        if o == "-c":
            configfile = a
        if o == "-l":
            exec("loglevel = logging." + a)
        if o == "-s":
            timestepsecs = int(a)
        if o == "-T":
            _lastTimeStep = int(a)
        if o == "-S":
            _firstTimeStep = int(a)
        if o == "-h":
            usage()

    if _lastTimeStep < _firstTimeStep:
        print(
            "The starttimestep ("
            + str(_firstTimeStep)
            + ") is cmaller than the last timestep ("
            + str(_lastTimeStep)
            + ")"
        )
        usage()

    myModel = WflowModel(wflow_cloneMap, caseName, runId, configfile)
    dynModelFw = wf_DynamicFramework(
        myModel, _lastTimeStep, firstTimestep=_firstTimeStep
    )
    dynModelFw.createRunId(NoOverWrite=False, level=loglevel)
    for o, a in opts:
        if o == "-I":
            configset(myModel.config, "model", "reinit", "1", overwrite=True)
        if o == "-s":
            configset(myModel.config, "model", "timestepsecs", a, overwrite=True)

    dynModelFw._runInitial()
    dynModelFw._runResume()
    # dynModelFw._runDynamic(0,0)
    dynModelFw._runDynamic(_firstTimeStep, _lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()


if __name__ == "__main__":
    main()
