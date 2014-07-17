#!/usr/bin/python

"""
Definition of the wflow_W3RA model.
---------------------------------------

The model is modified from the Australian Water Resources Assessment
Landscape (AWRA-L) model version 0.5

W3RA is documented in
van Dijk et al. (2013), Water Resour. Res., 49, 2729-2746, doi:10.1002/wrcr.20251
URL: http://onlinelibrary.wiley.com/doi/10.1002/wrcr.20251/abstract
More comprehensive documentation of AWRA-L version 0.5 can be found in:

Van Dijk, A.I.J.M. (2010) The Australian water resources assessment system
(version 0.5), 3.0.5.Technical description of the landscape hydrology model
(AWRA-L). WIRADA Technical Report, CSIRO Water for a Healthy Country
Flagship, Canberra.
URL: http://www.clw.csiro.au/publications/waterforahealthycountry/2010/wfhc-aus-water-resources-assessment-system.pdf
The section references below refer to the sections in the AWRA-L report.
Changes compared to that code are indicated, e.g. by commenting out
redundant code.

Further question please contact albert.vandijk@anu.edu.au
Port to Python/PCRaster: Deltares


Usage:
wflow_W3RA  -C case -R Runid -c inifile

    -C: set the name  of the case (directory) to run
    
    -R: set the name runId within the current case
    
    -c name of the config file (in the case directory)
    
$Author: schelle $
$Id: wflow_sceleton.py 898 2014-01-09 14:47:06Z schelle $
$Rev: 898 $
"""

import numpy
import os
import os.path
import shutil, glob
import getopt
try:
    from  wflow.wf_DynamicFramework import *
except ImportError:
    from  wf_DynamicFramework import *
    
try:
    from  wflow.wflow_adapt  import *
except ImportError:
    from  wflow_adapt import *  
import scipy



def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)

class WflowModel(DynamicModel):  
  """
  The user defined model class. T
  """
  
  def __init__(self, cloneMap,Dir,RunDir,configfile):
      """
      *Required*
      
      The init function **must** contain what is shown below. O
      ther functionality
      may be added by you if needed.
      
      """
      DynamicModel.__init__(self)   
      setclone(Dir + "/staticmaps/" + cloneMap)
      self.runId=RunDir      
      self.caseName=Dir
      self.Dir = Dir
      self.configfile = configfile
     


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
      states = ['S0','Ss','Sd','Sg','Sr','Mleaf','FreeWater','DrySnow']
      
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
      
      return self.currentTimeStep() * int(configget(self.config,'model','timestepsecs','86400'))
  
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


    self.timestepsecs = int(configget(self.config,'model','timestepsecs','86400'))
    
    self.basetimestep=86400
    self.SaveMapDir = self.Dir + "/" + self.runId + "/outmaps"

    # Define here the W3RA mapstacks (best to read these via netcdf)

    self.TEMP_mapstack=self.Dir + configget(self.config,"inputmapstacks","Temperature","/inmaps/TEMP")

    self.Pg_mapstack          =   self.Dir + configget(self.config,"inputmapstacks","Pg","/inmaps/Pg")
    self.Rg_mapstack          =   self.Dir + configget(self.config,"inputmapstacks","Rg","/inmaps/Rg")
    self.Ta_mapstack          =   self.Dir + configget(self.config,"inputmapstacks","Ta","/inmaps/Ta")
    self.T24_mapstack         =   self.Dir + configget(self.config,"inputmapstacks","T24","/inmaps/T24")
    self.pe_mapstack          =   self.Dir + configget(self.config,"inputmapstacks","pe","/inmaps/pe")
    self.pair_mapstack        =   self.Dir + configget(self.config,"inputmapstacks","pair","/inmaps/pair")
    self.u2_mapstack          =   self.Dir + configget(self.config,"inputmapstacks","u2","/inmaps/u2")
    self.fday_mapstack        =   self.Dir + configget(self.config,"inputmapstacks","fday","/inmaps/fday")
    self.ns_alb_mapstack      =   self.Dir + configget(self.config,"inputmapstacks","ns_alb","/inmaps/ns_alb")

    self.Altitude=readmap(self.Dir + "/staticmaps/wflow_dem")

    # Add reading of parameters here

    self.Nhru        = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.Fhru        = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.SLA         = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.LAIref      = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.Sgref       = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.S0FC        = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.SsFC        = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.SdFC        = readmap(self.Dir + "/staticmaps/Nhru.map")
    # fday
    self.Vc          = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.alb_dry     = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.alb_wet     = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.w0ref_alb   = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.Gfrac_max   = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.fvegref_G   = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.hveg        = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.Us0         = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.Ud0         = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.wslimU      = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.wdlimU      = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.cGsmax      = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.FsoilEmax   = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.w0limE      = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.FwaterE     = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.S_sls       = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.ER_frac_ref = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.InitLoss    = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.PrefR       = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.FdrainFC    = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.beta        = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.Fgw_conn    = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.K_gw        = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.K_rout      = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.LAImax      = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.Tgrow       = readmap(self.Dir + "/staticmaps/Nhru.map")
    self.Tsenc       = readmap(self.Dir + "/staticmaps/Nhru.map")


    self.logger.info("Starting Dynamic run...")


  def resume(self):
    """ 
    *Required*

    This function is required. Read initial state maps (they are output of a 
    previous call to suspend()). The implementation showns here is the most basic 
    setup needed.
    
    """
    self.logger.info("Reading initial conditions...")
    #: It is advised to use the wf_resume() function 
    #: here which pick up the variable save by a call to wf_suspend()
    try:
        self.wf_resume(self.Dir + "/instate/")
    except:
        self.logger.warn("Cannot load initial states, setting to default")
        for s in self.stateVariables():
            exec "self." + s + " = cover(1.0)"


  def default_summarymaps(self):
      """
      *Optional*

      Return a default list of variables to report as summary maps in the outsum dir.
      """
      return ['self.Altitude']

  def dynamic(self):
      """
      *Required*
      
      This is where all the time dependent functions are executed. Time dependent
      output should also be saved here.
      """

      #Put the W3RA here. Stuff from W3RA_timestep_model.m
      


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
    configfile="wflow_W3RA.ini"
    _lastTimeStep = 10
    _firstTimeStep = 1
    timestepsecs=86400
    wflow_cloneMap = 'wflow_subcatch.map'
    
    # This allows us to use the model both on the command line and to call 
    # the model usinge main function from another python script.
    
    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return     

    opts, args = getopt.getopt(argv, 'C:S:T:c:s:R:')
    
    for o, a in opts:
        if o == '-C': caseName = a
        if o == '-R': runId = a
        if o == '-c': configfile = a
        if o == '-s': timestepsecs = int(a)
        if o == '-T': _lastTimeStep=int(a)
        if o == '-S': _firstTimeStep=int(a)
        
    if (len(opts) <=1):
        usage()
        
    myModel = WflowModel(wflow_cloneMap, caseName,runId,configfile)
    dynModelFw = wf_DynamicFramework(myModel, _lastTimeStep,firstTimestep=_firstTimeStep)
    dynModelFw.createRunId(NoOverWrite=False,level=logging.DEBUG)    
    dynModelFw._runInitial()
    dynModelFw._runResume()
    dynModelFw._runDynamic(_firstTimeStep,_lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()
    

if __name__ == "__main__":
    main()
