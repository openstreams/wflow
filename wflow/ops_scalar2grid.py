#!/usr/bin/python

"""
Definition of the ops_scalar2grid program.
---------------------------------------

Converts sclar timeseries (in tss format) to grids.

Usage:
ops_scalar2grid  -C case -R Runid -c inifile

    -C: set the name  of the case (directory) to run
    
    -R: set the name runId within the current case
    
    -c name of the config file (in the case directory)

"""

from wflow.wf_DynamicFramework import *
from wflow.wflow_adapt import *


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args:
        print(msg)
    print(__doc__)
    sys.exit(0)


class WflowModel(DynamicModel):
    """
  The user defined model class. This is your work!
  """

    def __init__(self, cloneMap, Dir, RunDir, configfile):
        """
      *Required*
      
      The init function **must** contain what is shown below. Other functionality
      may be added by you if needed.
      
      """
        DynamicModel.__init__(self)
        setclone(Dir + "/staticmaps/" + cloneMap)
        self.runId = RunDir
        self.caseName = Dir
        self.Dir = Dir
        self.configfile = configfile

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

        # Static model parameters
        modelparameters.append(
            self.ParamType(
                name="Stations",
                stack="staticmaps/ops_scalar2grid_stations.map",
                type="staticmap",
                default=0.0,
                verbose=False,
                lookupmaps=[],
            )
        )

        return modelparameters

    def stateVariables(self):
        """ 
      *Required*
      
      Returns a list of state variables that are essential to the model. 
      This list is essential for the resume and suspend functions to work.
      
      This function is specific for each model and **must** be present. This is
      where you specify the state variables of you model. If your model is stateless
      this function must return and empty array (states = [])
      
      In the simple example here the TSoil variable is a state 
      for the model.
      
      :var TSoil: Temperature of the soil [oC]
      """
        states = []

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
        self.wf_suspend(self.Dir + "/outstate/")

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
        setglobaloption("unittrue")

        self.timestepsecs = int(
            configget(self.config, "model", "timestepsecs", "86400")
        )
        self.basetimestep = 86400
        # Reads all parameter from disk
        self.wf_updateparameters()

        self.interpolationmethod = configget(
            self.config, "model", "interpolationmethod", "inverse"
        )
        self.inversepower = int(configget(self.config, "model", "inversepower", "3"))
        self.ToInterpolate = configsection(self.config, "interpolate")
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
        return ["self.Stations"]

    def dynamic(self):
        """
      *Required*
      
      This is where all the time dependent functions are executed. Time dependent
      output should also be saved here.
      """

        self.wf_updateparameters()  # read the temperature map for each step (see parameters())

        self.Stations = ordinal(self.Stations)

        for var in self.ToInterpolate:
            tss = configget(self.config, "interpolate", var, None)
            tmp = timeinputscalar(self.Dir + "/" + tss, self.Stations)

            if self.interpolationmethod == "thiessen":
                Unq = uniqueid(boolean(abs(tmp) + 1.0))
                GaugeArea = spreadzone(ordinal(cover(Unq, 0)), 0, 1)
                exec("self." + var + " = areaaverage(tmp,GaugeArea)")
            elif self.interpolationmethod == "inverse":
                exec(
                    "self."
                    + var
                    + "=inversedistance(1,tmp,"
                    + str(self.inversepower)
                    + ",0,0)"
                )
            else:
                print("not implemented:" + self.interpolationmethod)


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
    configfile = "ops_scalar2grid.ini"
    _lastTimeStep = 10
    _firstTimeStep = 1
    timestepsecs = 86400
    wflow_cloneMap = "ops_scalar2grid_stations.map"

    # This allows us to use the model both on the command line and to call
    # the model usinge main function from another python script.

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return

    opts, args = getopt.getopt(argv, "C:S:T:c:s:R:")

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

    if len(opts) <= 1:
        usage()

    myModel = WflowModel(wflow_cloneMap, caseName, runId, configfile)
    dynModelFw = wf_DynamicFramework(
        myModel, _lastTimeStep, firstTimestep=_firstTimeStep
    )
    dynModelFw.createRunId(NoOverWrite=False, level=logging.DEBUG)
    dynModelFw._runInitial()
    dynModelFw._runResume()
    dynModelFw._runDynamic(0, 0)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()


if __name__ == "__main__":
    main()
