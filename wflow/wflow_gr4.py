#!/usr/bin/python

"""
Definition of the wflow_gr4 model.
----------------------------------


Usage:
wflow_gr4  [-l loglevel][-c configfile][-f][-h] -C case -R Runid -

    -C: set the name  of the case (directory) to run
    
    -R: set the name runId within the current case
    
    -c name of the config file (in the case directory)
    
    -f: Force overwrite of existing results

    -h: print usage information
    
    -l: loglevel (most be one of DEBUG, WARNING, ERROR)
    
$Author: schelle $
$Id: wflow_gr4.py 923 2014-03-13 13:48:37Z schelle $
$Rev: 923 $


NOTES
-----

- The max length of the arrays is determined  by the X4 parameter (int(X4))
- The X4 parameter is always uniform over that catchment. Howvere, the state
  of the UH is determined per grid cell.

"""

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


def pcr_tanh(x):
    """
    define tanh for pcraster objects
    
    """
    return (pcr.exp(x) - pcr.exp(-x)) / (pcr.exp(x) + pcr.exp(-x))


def initUH1(X4, D):
    """
    Initialize the UH1 unit hydrograph
    
    Input:   
        - X4
        - D
    
    Returns:   
        - UH1, SH1
    """
    NH = int(numpy.ceil(X4))

    t = np.arange(1, NH + 1)
    SH1 = numpy.minimum(1.0, (t / X4) ** D)

    # Use numpy.diff to get the UH, insert value at zero to complete
    UH1 = numpy.diff(SH1, axis=0)
    UH1 = numpy.insert(UH1, 0, SH1[0])
    return UH1, SH1


def initUH2(X4, D):
    """
    Initialize the UH2 unit hydrograph
    
    Input:
    
        - X4
        - D
    
    Returns:
    
        - UH2, SH2
    """
    NH = int(numpy.ceil(X4))

    t1 = np.arange(1, NH)
    t2 = np.arange(NH, 2 * NH + 1)

    SH2_1 = 0.5 * (t1 / X4) ** D
    SH2_2 = 1 - 0.5 * (numpy.maximum(0, 2 - t2 / X4)) ** D

    SH2 = numpy.minimum(1.0, numpy.hstack((SH2_1, SH2_2)))

    # Use numpy.diff to get the UH, insert value at zero to complete
    UH2 = numpy.diff(SH2, axis=0)
    UH2 = numpy.insert(UH2, 0, SH2[0])
    return UH2, SH2


def mk_qres(N):
    """
    Returns an array (or ayyar of maps) to store the 
    delayed flow in
    
    Input:
    
        - N nr op steps
        
    Ouput:
    
        - nr of steps elemenst initialized with zeros's
    """

    uhq = []

    for i in range(0, N):
        uhq.append(pcr.cover(0.0))

    return uhq


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

    def stateVariables(self):
        """ 
      returns a list of state variables that are essential to the model. 
      This list is essential for the resume and suspend functions to work.
      
      This function is specific for each model and **must** be present.


      :var self.S_X1: production reservoir content at the beginning of the time step (divided by X1) [mm]
      :var self.R_X3: routing reservoir content at the beginning of the time step (divided by X3) [mm]
      
      .. todo::
      
          add routing state vars
          
      """
        states = ["S_X1", "R_X3", "QUH1", "QUH2"]

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
            configget(self.config, "model", "timestepsecs", "3600")
        )

    def suspend(self):
        """
      *Required*
      
      Suspends the model to disk. All variables needed to restart the model
      are saved to disk as pcraster maps. Use resume() to re-read them
      
      This function is required. 
      
    """

        self.logger.info("Saving initial conditions...")
        self.wf_suspend(os.path.join(self.SaveDir, "outstate"))

        if self.OverWriteInit:
            self.logger.info("Saving initial conditions over start conditions...")
            self.wf_suspend(os.path.join(self.SaveDir, "/instate"))

    def initial(self):

        """
    Initial part of the gr4 model, executed only once. Reads all static model
    information (parameters) and sets-up the variables used in modelling.
    
    :var dt.tbl: time step (1) [hour]
    :var B.tbl: routing ratio (0.9) [-]
    :var NH: UH dimension (number) taken from ini file [-]
    :var D.tbl: variable for hourly time steps (1.25) [-]
    :var C.tbl: variable (number) [hour]
    
    *Parameters*
    
    :var X1.tbl: capacity of the production store, accounts for soil moisture (number) [mm]
    :var X2.tbl: water exchange coefficient (number) [mm]
    :var X3.tbl: capacity of the routing store (number) [mm]
    :var X4 (in ini): time base of the unit hydrograph (number) [hour]
    
    
    """
        #: pcraster option to calculate with units or cells. Not really an issue
        #: in this model but always good to keep in mind.
        pcr.setglobaloption("unittrue")
        self.thestep = pcr.scalar(0)
        self.ZeroMap = pcr.cover(0.0)

        self.timestepsecs = int(configget(self.config, "model", "timestepsecs", "3600"))
        self.basetimestep = 3600
        self.reinit = int(configget(self.config, "run", "reinit", "0"))
        self.OverWriteInit = int(configget(self.config, "model", "OverWriteInit", "0"))
        self.SaveMapDir = self.Dir + "/" + self.runId + "/outmaps"
        self.TEMP_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Temperature", "/inmaps/TEMP"
        )
        self.intbl = configget(self.config, "model", "intbl", "intbl")
        self.Altitude = pcr.readmap(self.Dir + "/staticmaps/wflow_dem")

        wflow_subcatch = configget(
            self.config, "model", "wflow_subcatch", "/staticmaps/wflow_subcatch.map"
        )
        wflow_landuse = configget(
            self.config, "model", "wflow_landuse", "/staticmaps/wflow_landuse.map"
        )
        wflow_soil = configget(
            self.config, "model", "wflow_soil", "/staticmaps/wflow_soil.map"
        )
        self.P_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Precipitation", "/inmaps/P"
        )  # timeseries for rainfall
        self.PET_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "EvapoTranspiration", "/inmaps/PET"
        )  # timeseries for rainfall"/inmaps/PET"          # potential evapotranspiration
        sizeinmetres = int(configget(self.config, "layout", "sizeinmetres", "0"))
        subcatch = pcr.ordinal(
            pcr.readmap(self.Dir + wflow_subcatch)
        )  # Determines the area of calculations (all cells > 0)
        subcatch = pcr.ifthen(subcatch > 0, subcatch)
        self.xl, self.yl, self.reallength = pcrut.detRealCellLength(
            self.ZeroMap, sizeinmetres
        )
        self.ToCubic = (
            self.reallength * self.reallength * 0.001
        ) / self.timestepsecs  # m3/s

        self.LandUse = pcr.readmap(
            self.Dir + wflow_landuse
        )  #: Map with lan-use/cover classes
        self.LandUse = pcr.cover(self.LandUse, pcr.nominal(pcr.ordinal(subcatch) > 0))
        self.Soil = pcr.readmap(self.Dir + wflow_soil)  #: Map with soil classes
        self.Soil = pcr.cover(self.Soil, pcr.nominal(pcr.ordinal(subcatch) > 0))
        self.OutputId = pcr.readmap(
            self.Dir + wflow_subcatch
        )  # location of subcatchment

        # hourly time step
        self.dt = int(configget(self.config, "gr4", "dt", "1"))
        # routing ratio found in criteria validation file, first line
        self.B = float(configget(self.config, "gr4", "B", "0.9"))
        # hourly time-steps
        self.D = float(configget(self.config, "gr4", "D", "1.25"))

        # The following parameters are spatial (apart from X4)
        # capacity of the production store, accounts for soil moisture (mm) (>=0)
        self.X1 = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/X1.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            285.72,
        )
        # water exchange coefficient
        self.X2 = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/X2.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            -0.42,
        )
        # capacity of the routing store (mm)
        self.X3 = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/X3.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            169.02,
        )
        # time base of the unit hydrograph (hr)
        # self.X4=self.readtblDefault(self.Dir + "/" + self.intbl + "/X4.tbl",self.LandUse,subcatch,self.Soil,32.85)
        self.X4 = float(configget(self.config, "gr4", "X4", "32.85"))
        # Set static initial values here #########################################
        # Number of UH units
        self.NH = int(numpy.ceil(self.X4))

        self.UH1, self.SH1 = initUH1(self.X4, self.D)
        self.UH2, self.SH2 = initUH2(self.X4, self.D)

        self.QUH1 = mk_qres(self.NH)
        self.QUH2 = mk_qres(self.NH * 2)

        self.logger.info("End of initial section...")

    def resume(self):
        """ 
    *Required*

    This function is required. Read initial state maps (they are output of a 
    previous call to suspend()). The implementation showns her is the most basic 
    setup needed.
    
    """
        self.logger.info("Reading initial conditions...")
        #: It is advised to use the wf_resume() function
        #: here which pick upt the variable save by a call to wf_suspend()
        if self.reinit == 1:
            # STATES
            self.S_X1 = 245.4900 / self.X1  # STATE(1),level in production store
            self.R_X3 = 43.9031 / self.X3  # STATE(2),level in routing store
            self.QUH1 = mk_qres(self.NH)
            self.QUH2 = mk_qres(self.NH * 2)
        else:
            self.wf_resume(os.path.join(self.Dir, "instate"))

    def dynamic(self):
        """
      *Required*
      
      :var self.Pn: net precipitation [mm]
      :var self.En: net evapotranspiration [mm]
      :var self.Ps: part of Pn that feeds the production reservoir [mm]
      :var self.Es: evaporation quantity substracted from the production reservoir [mm]
      """

        self.logger.debug(
            "Step: "
            + str(int(self.thestep + self._d_firstTimeStep))
            + "/"
            + str(int(self._d_nrTimeSteps))
        )
        self.thestep = self.thestep + 1

        self.Precipitation = pcr.cover(self.wf_readmap(self.P_mapstack, 0.0), 0.0)
        self.PotEvaporation = pcr.cover(self.wf_readmap(self.PET_mapstack, 0.0), 0.0)

        # ROUTING WATER AND PRODUCTION RESERVOIR PERCOLATION ========================================================

        self.Pn = pcr.ifthenelse(
            self.Precipitation >= self.PotEvaporation,
            self.Precipitation - self.PotEvaporation,
            pcr.scalar(0.0),
        )
        self.En = pcr.ifthenelse(
            self.Precipitation >= self.PotEvaporation,
            pcr.scalar(0.0),
            self.PotEvaporation - self.Precipitation,
        )
        self.Ps = (self.X1 * (1 - (self.S_X1) ** 2) * pcr_tanh(self.Pn / self.X1)) / (
            1 + self.S_X1 * pcr_tanh(self.Pn / self.X1)
        )
        self.Es = (
            self.S_X1 * self.X1 * (2 - self.S_X1) * pcr_tanh(self.En / self.X1)
        ) / (1 + (1 - self.S_X1) * pcr_tanh(self.En / self.X1))
        self.Ps = pcr.ifthenelse(
            self.Precipitation >= self.PotEvaporation, self.Ps, pcr.scalar(0.0)
        )
        self.Es = pcr.ifthenelse(
            self.Precipitation >= self.PotEvaporation, pcr.scalar(0.0), self.Es
        )

        self.Sprim_X1 = (
            self.S_X1 + ((self.Ps - self.Es) * self.dt) / self.X1
        )  # reservoir new content
        # Filter out value < 0 in self.Sprim_X1
        self.Sprim_X1 = pcr.max(0.0, self.Sprim_X1)

        self.Perc = (
            self.Sprim_X1 * self.X1 * (1 - (1 + (self.Sprim_X1 / 5.25) ** 4) ** -0.25)
        )  # percolation
        self.S_X1 = (
            self.Sprim_X1 - (self.Perc * self.dt) / self.X1
        )  # reservoir new content

        self.Pr = self.Perc + (self.Pn - self.Ps)  # quantity to routing

        # ACTUAL ROUTING =====================================================
        # UH1 has a memory of int(X4) steps

        # ouput of UH1 =========================================================

        for j in range(0, self.NH):  # UH1 output for each time step
            self.QUH1[j] = self.QUH1[j] + float(self.UH1[j]) * self.Pr

        self.Q9 = self.B * self.QUH1[0]

        # Add the current Q to the UH res
        for j in range(0, 2 * self.NH):  # UH2 output for each time step
            self.QUH2[j] = self.QUH2[j] + float(self.UH2[j]) * self.Pr

        self.Q1prim = self.QUH2[0]
        # Get final runoff
        self.Q1 = (1 - self.B) * self.Q1prim
        self.F = self.X2 * (self.R_X3) ** 3.5  # water subterranean exchange
        self.Rprim_X3 = (
            self.R_X3 + (self.Q9 + self.F) / self.X3
        )  # new routing reservoir level
        self.Qr = (
            self.Rprim_X3 * self.X3 * (1.0 - (1.0 + (self.Rprim_X3) ** 4) ** -0.25)
        )  # routing output
        self.R_X3 = self.Rprim_X3 - self.Qr / self.X3  # new routing reservoir level
        self.Qd = pcr.max(0.0, self.Q1 + self.F)  # flow component Qd
        self.Q = self.Qr + self.Qd  # total flow Q in mm/hr
        # Updated this line to get total Q per basin
        self.SurfaceRunoff = pcr.areatotal(self.Q * self.ToCubic, self.OutputId)

        # Remove first item from the UH stacks and add a new empty one at the end
        self.QUH1 = delete(self.QUH1, 0)
        self.QUH1 = append(self.QUH1, pcr.cover(0.0))
        self.QUH2 = delete(self.QUH2, 0)
        self.QUH2 = append(self.QUH2, pcr.cover(0.0))


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
    configfile = "wflow_gr4.ini"
    _lastTimeStep = 0
    _firstTimeStep = 0
    timestepsecs = 3600
    wflow_cloneMap = "wflow_subcatch.map"
    NoOverWrite = True
    loglevel = logging.DEBUG

    # This allows us to use the model both on the command line and to call
    # the model usinge main function from another python script.

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return

    opts, args = getopt.getopt(argv, "C:S:T:c:s:R:fhIXi:l:")

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
        if o == "-f":
            NoOverWrite = 0
        if o == "-h":
            usage()
        if o == "-l":
            exec("loglevel = logging." + a)

    if len(opts) <= 1:
        usage()

    if _lastTimeStep < _firstTimeStep:
        print(
            "The starttimestep ("
            + str(_firstTimeStep)
            + ") is smaller than the last timestep ("
            + str(_lastTimeStep)
            + ")"
        )
        usage()

    myModel = WflowModel(wflow_cloneMap, caseName, runId, configfile)
    dynModelFw = wf_DynamicFramework(
        myModel, _lastTimeStep, firstTimestep=_firstTimeStep
    )
    dynModelFw.createRunId(NoOverWrite=NoOverWrite, level=loglevel)

    for o, a in opts:
        if o == "-X":
            configset(myModel.config, "model", "OverWriteInit", "1", overwrite=True)
        if o == "-I":
            configset(myModel.config, "model", "reinit", "1", overwrite=True)
        if o == "-i":
            configset(myModel.config, "model", "intbl", a, overwrite=True)
        if o == "-s":
            configset(myModel.config, "model", "timestepsecs", a, overwrite=True)
        if o == "-c":
            configset(myModel.config, "model", "configfile", a, overwrite=True)

    dynModelFw._runInitial()
    dynModelFw._runResume()
    # dynModelFw._runDynamic(0,0)
    dynModelFw._runDynamic(_firstTimeStep, _lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()


if __name__ == "__main__":
    main()
