#!/usr/bin/python

"""
Definition of the wflow_sceleton model.
---------------------------------------

This simple model calculates soil temperature using
air temperature as a forcing.

Usage:
wflow_sceleton  -C case -R Runid -c inifile

    -C: set the name  of the case (directory) to run
    
    -R: set the name runId within the current case
    
    -c name of the config file (in the case directory)
    
$Author: schelle $
$Id: wflow_sceleton.py 898 2014-01-09 14:47:06Z schelle $
$Rev: 898 $
"""

import reservoir_Si
# import reservoir_Sa
import reservoir_Su
import reservoir_Sf
import reservoir_Ss
import reservoir_Sus
import JarvisCoefficients

import pdb
import numpy
import os
import os.path
import shutil, glob
import getopt

from wflow.wf_DynamicFramework import *
from wflow.wflow_adapt import *
#import scipy
from copy import copy as copylist

# TODO: see below
"""
Inlezen tijdseries (grids)
Nieuwe lezen parameters
Reservoir nul een doorgeefreservoir maken
Multiplication with cell surface aanpassen
Verwijderen IRURFR_L statements?
Documentatie updaten!
"""


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)

class WflowModel(DynamicModel):  
  """
  The user defined model class. This is your work!
  """
  
  def __init__(self, cloneMap,Dir,RunDir,configfile):
      """
      *Required*
      
      The init function **must** contain what is shown below. Other functionality
      may be added by you if needed.
      
      """
      DynamicModel.__init__(self)   
      setclone(os.path.join(Dir, 'staticmaps', cloneMap))
      self.runId = RunDir
      self.caseName = os.path.abspath(Dir)
      self.Dir = os.path.abspath(Dir)
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

      #Static model parameters
      modelparameters.append(self.ParamType(name="Altitude",stack="staticmaps/wflow_dem.map",type="staticmap",default=0.0,verbose=False,lookupmaps=[]))

      # Meteo and other forcing
      modelparameters.append(self.ParamType(name="Temperature",stack="inmaps/TEMP",type="timeseries",default=10.0,verbose=False,lookupmaps=[]))

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
      states = ['Si','Su','Sus','Sf','Ss']
                       
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
    [report(self.Si[i], self.SaveDir + "/outmaps/Si" + self.NamesClasses[i] +".map") for i in self.Classes]        
    [report(self.Su[i], self.SaveDir + "/outmaps/Su" + self.NamesClasses[i] +".map") for i in self.Classes]        
    [report(self.Sus[i], self.SaveDir + "/outmaps/Sus" + self.NamesClasses[i] +".map") for i in self.Classes]        
    [report(self.Sf[i], self.SaveDir + "/outmaps/Sf" + self.NamesClasses[i] +".map") for i in self.Classes]        
    [report(self.Sr[i], self.SaveDir + "/outmaps/Sr" + self.NamesClasses[i] +".map") for i in self.Classes]        
    report(self.Ss, self.SaveDir + "/outmaps/Ss.map")
    
    [report(self.percent[i], self.SaveDir + "/outmaps/percent" + self.NamesClasses[i] +".map") for i in self.Classes]        
    report(self.percentArea,self.SaveDir + "/outmaps/percentArea.map")
    report(self.surfaceArea,self.SaveDir + "/outmaps/surfaceArea.map")
        
    report(self.sumprecip,self.SaveDir + "/outsum/sumprecip.map")
    report(self.sumevap,self.SaveDir + "/outsum/sumevap.map")
    report(self.sumpotevap,self.SaveDir + "/outsum/sumpotevap.map")
    report(self.sumtemp,self.SaveDir + "/outsum/sumtemp.map")
    report(self.sumrunoff,self.SaveDir + "/outsum/sumrunoff.map")
    report(self.sumwb,self.SaveDir + "/outsum/sumwb.map")


      
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

    self.teller=0
    
    self.thestep = scalar(0)
    #: files to be used in case of timesries (scalar) input to the model
    #files for forcing data
    self.precipTss = os.path.join(self.Dir,
                                  configget(self.config, "model", "Pfile_1", ""))
    self.evapTss = os.path.join(self.Dir,
                                configget(self.config, "model", "Efile_1", ""))
    self.tempTss = os.path.join(self.Dir,
                                configget(self.config, "model", "Tfile_1", ""))
    self.precipTss2 = os.path.join(self.Dir,
                                   configget(self.config, "model", "Pfile_2", "")) 
    self.evapTss2 = os.path.join(self.Dir,
                                 configget(self.config, "model", "Efile_2", ""))
    self.tempDMTss = os.path.join(self.Dir,
                                  configget(self.config, "model", "TDMfile_2", ""))
    self.radnTss = os.path.join(self.Dir,
                                configget(self.config, "model", "RNfile_2", ""))
    self.radsTss = os.path.join(self.Dir,
                                configget(self.config, "model", "RSfile_2", ""))
    self.sgammaTss = os.path.join(self.Dir,
                                  configget(self.config, "model", "SGfile_2", ""))
    self.vpdTss = os.path.join(self.Dir,
                               configget(self.config, "model", "VPDfile_2", ""))
    self.windTss = os.path.join(self.Dir,
                                configget(self.config, "model", "Wfile_2", ""))
    self.daySTss = os.path.join(self.Dir,
                                configget(self.config, "model", "DSfile_2", ""))
    self.dayETss = os.path.join(self.Dir,
                                configget(self.config, "model", "DEfile_2", ""))
#    self.laiTss = configget(self.config,"model","LAIfile_2","")
         
  
    self.logger.info("running for " + str(self.nrTimeSteps()) + " timesteps")           #keeping track of number of timesteps
    
        
    # Set and get defaults from ConfigFile here ###################################
    self.timestepsecs = int(configget(self.config,
                                      "model", "timestepsecs", "3600"))  # number of seconds in a timestep   
    self.scalarInput = int(configget(self.config,
                                     "model", "ScalarInput", "1"))  # forcing data in maps (0) or timeseries (1)
    self.InputSeries = int(configget(self.config,
                                     "model", "InputSeries", "1"))  # forcing data in maps (0) or timeseries (1)
    self.reinit = int(configget(self.config,
                                "model", "reinit", "0"))    
    self.maxTransit = float(configget(self.config,
                                    "model", "maxTransitTime", "34"))  # maximum Transit time in cacthment
    self.distForcing = int(configget(self.config,
                                     "model", "DistForcing", "10"))             #number of different forcing inputs (eg. number of rainfall stations)
    self.maxGaugeId = int(configget(self.config,
                                    "model", "maxGaugeId", "10"))  # highest index of all used meteo stations
    self.IRURFR_L = int(configget(self.config,
                                  "model", "L_IRURFR", "0"))  # combination of reservoirs that are distributed (1: all these reservoirs are distributed)        
    self.URFR_L = int(configget(self.config,
                                "model", "L_URFR", "0"))  # combination of reservoirs that are distributed (1: all these reservoirs are distributed)            
    self.FR_L = int(configget(self.config,
                              "model", "L_FR", "0"))  # combination of reservoirs that are distributed (1: all these reservoirs are distributed)        
    self.Ctime = int(configget(self.config,
                               "model", "spinUp_time", "7775"))  # number of timesteps for which no data needs to be recorded
    self.NamesClasses = eval(str(configget(self.config,
                                           "model", "classes", "['W','H','P']")))  # classes used in model
    self.Classes = [x for x in range(len(self.NamesClasses))]  # numbering of classes
    
    # selection of reservoir conceputalisatie - codes are described in reservoir files
    self.selectSi = configget(self.config, "model",
                             "selectSi", "0, 0, 0").replace(
                             ' ', '').replace('[', '').replace(
                             ']', '').replace(
                             'None', '').split(',')
    self.selectSu = configget(self.config, "model",
                              "selectSu", "0, 0, 0").replace(
                             ' ', '').replace('[', '').replace(
                             ']', '').replace(
                             'None', '').split(',')
    self.selectSus = configget(self.config, "model",
                               "selectSus", "0, 0, 0").replace(
                             ' ', '').replace('[', '').replace(
                             ']', '').replace(
                             'None', '').split(',')
    self.selectSf = configget(self.config, "model",
                              "selectSf", "0, 0, 0").replace(
                             ' ', '').replace('[', '').replace(
                             ']', '').replace(
                             'None', '').split(',')
    self.selectSs = configget(self.config, "model", "selectSs", "groundWaterCombined3")
    self.selectSr = configget(self.config, "model",
                              "selectSr", "0, 0, 0").replace(
                             ' ', '').replace('[', '').replace(
                             ']', '').replace(
                             'None', '').split(',')
    # static maps to use (normally default)
    wflow_subcatch = configget(self.config,
                               "model", "wflow_subcatch", "staticmaps/wflow_subcatch.map")
    wflow_dem  = configget(self.config,
                           "model", "wflow_dem", "staticmaps/wflow_dem.map")
    wflow_ldd = configget(self.config,
                          "model", "wflow_ldd", "staticmaps/wflow_ldd.map")
    wflow_gauges  = configget(self.config,
                              "model", "wflow_gauges", "staticmaps/wflow_gauges.map")
    wflow_mgauges  = configget(self.config,
                               "model", "wflow_mgauges", "staticmaps/wflow_mgauges.map")
    wflow_surfaceArea = configget(self.config,
                                  "model", "wflow_surfaceArea", "staticmaps/wflow_surfaceArea.map")
    wflow_transit = configget(self.config,
                              "model", "wflow_transit", "staticmaps/wflow_transit.map")  
    wflow_percent = [configget(self.config,
                               "model", "wflow_percent_" + str(self.Classes[i]),
                               "staticmaps/wflow_percent" + str(self.Classes[i]) +".map") for i in self.Classes]
    self.rst_laiTss = [configget(self.config,
                                 "model","rst_lai_" + str(self.Classes[i]),
                                 "staticmaps/rst_lai_" + str(self.Classes[i]) +".map") for i in self.Classes]
    
    # 2: Input base maps ########################################################  
    subcatch = ordinal(readmap(os.path.join(self.Dir, wflow_subcatch))) # Determines the area of calculations (all cells > 0)
    subcatch = ifthen(subcatch > 0, subcatch)
    
    self.Altitude = readmap(os.path.join(self.Dir, wflow_dem)) * scalar(defined(subcatch)) #: The digital elevation map (DEM)
    self.TopoLdd = readmap(os.path.join(self.Dir, wflow_ldd))        #: The local drinage definition map (ldd)
    self.TopoId = readmap(os.path.join(self.Dir, wflow_subcatch))        #: Map define the area over which the calculations are done (mask)
    self.TopoId = ifthen(scalar(self.TopoId) > 0,self.TopoId)
    self.surfaceArea = scalar(readmap(os.path.join(self.Dir, wflow_surfaceArea)))        #: Map with surface area per cell
    self.totalArea = areatotal(self.surfaceArea, nominal(self.TopoId))
    self.percentArea = self.surfaceArea/self.totalArea    
    self.Transit = scalar(readmap(os.path.join(self.Dir, wflow_transit)))        #: Map with surface area per cell
    self.gaugesR = nominal(readmap(os.path.join(self.Dir, wflow_gauges)))    
    self.percent = []
    for i in self.Classes:
        self.percent.append(readmap(os.path.join(self.Dir, wflow_percent[i])))

    #MODEL PARAMETERS
    self.sumax = eval(str(configget(self.config, "model", "sumax", "[0]")))
    self.sumin = eval(str(configget(self.config, "model", "sumin", "[0]")))
    self.samax = eval(str(configget(self.config, "model", "samax", "[0]")))
    self.susmax1 = eval(str(configget(self.config, "model", "susmax1", "[0]")))
    self.susmax2 = eval(str(configget(self.config, "model", "susmax2", "[0]")))
    self.susmax3 = eval(str(configget(self.config, "model", "susmax3", "[0]")))
    self.srmax = eval(str(configget(self.config, "model", "sumax", "[0]")))
    self.beta = eval(str(configget(self.config, "model", "beta", "[0]")))
    self.famax = eval(str(configget(self.config, "model", "famax", "[0]")))
    self.Ce = eval(str(configget(self.config, "model", "Ce", "[0]")))
    self.Co = eval(str(configget(self.config, "model", "Co", "[0]")))
    self.D = eval(str(configget(self.config, "model", "D", "[0]")))
    self.Kf = eval(str(configget(self.config, "model", "Kf", "[0]")))
    self.Tf = eval(str(configget(self.config, "model", "Tf", "[0]")))
    self.imax = eval(str(configget(self.config, "model", "imax", "[0]")))
    self.perc = eval(str(configget(self.config, "model", "perc", "[0]")))
    self.cap = eval(str(configget(self.config, "model", "cap", "[0]")))
    self.Kd = eval(str(configget(self.config, "model", "Kd", "[0]")))
    self.Kr = eval(str(configget(self.config, "model", "Kr", "[0]")))
    self.LP = eval(str(configget(self.config, "model", "LP", "[0]")))
    self.Ks = eval(str(configget(self.config, "model", "Ks", "[0]")))
    #Jarvis stressfunctions
    self.JC_Topt = eval(str(configget(self.config, "model", "JC_Topt", "[0]")))
    self.JC_D05 = eval(str(configget(self.config, "model", "JC_D05", "[0]")))
    self.JC_cd1 = eval(str(configget(self.config, "model", "JC_cd1", "[0]")))
    self.JC_cd2 = eval(str(configget(self.config, "model", "JC_cd2", "[0]")))
    self.JC_cr = eval(str(configget(self.config, "model", "JC_cr", "[0]")))
    self.JC_cuz = eval(str(configget(self.config, "model", "JC_cuz", "[0]")))
    self.SuFC = eval(str(configget(self.config, "model", "SuFC", "[0]")))
    self.SuWP = eval(str(configget(self.config, "model", "SuWP", "[0]")))
    self.JC_rstmin = eval(str(configget(self.config, "model", "JC_rstmin", "[0]")))
    self.gamma = eval(str(configget(self.config, "model", "gamma", "[0]")))
    self.Cp = eval(str(configget(self.config, "model", "Cp", "[0]")))
    self.rhoA = eval(str(configget(self.config, "model", "rhoA", "[0]")))
    self.rhoW = eval(str(configget(self.config, "model", "rhoW", "[0]")))
    self.lamda = eval(str(configget(self.config, "model", "lamda", "[0]")))
    
    # initialise list for routing
    self.trackQ = [0*scalar(self.TopoId)] * int(self.maxTransit)
    
    # initialise list for lag function
    self.convQu = [[0*scalar(self.TopoId)] * self.Tf[i] for i in self.Classes]
       

    if self.scalarInput:
        self.gaugesMap = nominal(readmap(os.path.join(self.Dir, wflow_mgauges)))  #: Map with locations of rainfall/evap/temp gauge(s). Only needed if the input to the model is not in maps
    self.OutputId = readmap(os.path.join(self.Dir, wflow_subcatch))  # location of subcatchment
    self.OutputIdRunoff = boolean(ifthenelse(self.gaugesR == 1,1*scalar(self.TopoId),0*scalar(self.TopoId)))       # location of subcatchment  

    self.ZeroMap = 0.0*scalar(subcatch)                    #map with only zero's
  
  
    # For in memory override:
    self.P = self.ZeroMap
    self.PET = self.ZeroMap
    self.TEMP = self.ZeroMap
    
   
    self.logger.info("Linking parameters to landuse, catchment and soil...")

   
    # Initializing of variables

    self.logger.info("Initializing of model variables..")
    self.TopoLdd = lddmask(self.TopoLdd,boolean(self.TopoId))   
    catchmentcells = maptotal(scalar(self.TopoId))
 
    self.sumprecip = self.ZeroMap  # accumulated rainfall for water balance
    self.sumevap = self.ZeroMap  # accumulated evaporation for water balance
    self.sumrunoff = self.ZeroMap  # accumulated runoff for water balance (weigthted for upstream area)
    self.sumpotevap = self.ZeroMap  # accumulated runoff for water balance
    self.sumtemp = self.ZeroMap                          #accumulated runoff for water balance
    self.Q = self.ZeroMap
    self.sumwb = self.ZeroMap
       
    # Define timeseries outputs There seems to be a bug and the .tss files are 
    # saved in the current dir...
    
    self.SaveDir = os.path.join(self.Dir, self.runId)
    self.logger.info("Starting Dynamic run...")

  def resume(self):
    """ 
    *Required*

    This function is required. Read initial state maps (they are output of a 
    previous call to suspend()). The implementation shown here is the most basic
    setup needed.
    
    """
    if self.reinit == 1:
        #self.logger.info("Setting initial conditions to default (zero!)")
        self.logger.info("Setting initial conditions to preset values in main script!!")
        self.Si = [self.ZeroMap] * len(self.Classes)
        self.Su = [self.ZeroMap] * len(self.Classes) 
        self.Sa = [self.ZeroMap] * len(self.Classes) 
        self.Sus = [self.ZeroMap] * len(self.Classes)
        self.Sf = [self.ZeroMap] * len(self.Classes)
        self.Sr = [self.ZeroMap] * len(self.Classes)
        #self.Ss = [self.ZeroMap] * len(self.Classes)       # for separate gw reservoir per class
        self.Ss = self.ZeroMap                              # for combined gw reservoir 
        
        # set initial storage values
        self.Sa = [x + y for (x,y) in zip(self.Su, [130 * scalar(self.TopoId)] * len(self.Classes))]
        self.Su = [x + y for (x,y) in zip(self.Su, [130 * scalar(self.TopoId)] * len(self.Classes))]
        self.Sus = [x + y for (x,y) in zip(self.Sus, [30 * scalar(self.TopoId)] * len(self.Classes))]
        self.Ss = self.Ss + 30 * scalar(self.TopoId)           # for combined gw reservoir 
    
    else:
        self.wf_resume(self.Dir + "/instate/")  
    
    self.wbSi_ = [self.ZeroMap] * len(self.Classes)    
    self.wbSu_ = [self.ZeroMap] * len(self.Classes)    
    self.wbSa_ = [self.ZeroMap] * len(self.Classes)
    self.wbSus_ = [self.ZeroMap] * len(self.Classes)    
    self.wbSr_ = [self.ZeroMap] * len(self.Classes)    
    self.wbSf_ = [self.ZeroMap] * len(self.Classes)    
    self.wbSfrout = self.ZeroMap
    self.wbSs = self.ZeroMap
    
    self.Ei_ = [self.ZeroMap] * len(self.Classes)
    self.Pe_ = [self.ZeroMap] * len(self.Classes)
    self.Si_ = [self.ZeroMap] * len(self.Classes)
    self.Eu_ = [self.ZeroMap] * len(self.Classes)
    self.Er_ = [self.ZeroMap] * len(self.Classes)
    self.Qu_ = [self.ZeroMap] * len(self.Classes)
    self.Qd_ = [self.ZeroMap] * len(self.Classes)
    self.Qo_ = [self.ZeroMap] * len(self.Classes)
    self.Qr_ = [self.ZeroMap] * len(self.Classes)
    self.Cap_ = [self.ZeroMap] * len(self.Classes)
    self.Perc_ = [self.ZeroMap] * len(self.Classes)
    self.Fa_ = [self.ZeroMap] * len(self.Classes)
    self.Qf_ = [self.ZeroMap] * len(self.Classes)
    #self.Qs_ = [self.ZeroMap] * len(self.Classes)       # for separate gw reservoir per class
    self.Qs_ = self.ZeroMap                             # for combined gw reservoir 
    self.Qflag_ = [self.ZeroMap] * len(self.Classes)
    self.Qfcub_ = [self.ZeroMap] * len(self.Classes)
    self.Ep_ = [self.ZeroMap] * len(self.Classes)
    self.EpD_ = [self.ZeroMap] * len(self.Classes)
    
    self.JC_temp_ = [self.ZeroMap] * len(self.Classes)
    self.JC_vpd_ = [self.ZeroMap] * len(self.Classes)
    self.JC_rad_ = [self.ZeroMap] * len(self.Classes)
    self.JC_sm_ = [self.ZeroMap] * len(self.Classes)
    self.JC_k_ = [self.ZeroMap] * len(self.Classes)


  def default_summarymaps(self):
      """
      *Optional*

      Return a default list of variables to report as summary maps in the outsum dir.
      The ini file has more options, including average and sum
      """
      return ['self.Altitude']

  def dynamic(self):
        """
        *Required*
        
        This is where all the time dependent functions are executed. Time dependent
        output should also be saved here.
        """
        # TODO: change rainfall .tss files into grids
        # self.wf_updateparameters() # read the temperature map for each step (see parameters())
        self.teller=self.teller+1
        
        #self.logger.debug("Step: "+str(int(self.thestep + self._d_firstTimeStep))+"/"+str(int(self._d_nrTimeSteps)))
        self.thestep = self.thestep + 1
        
        self.Si_t = copylist(self.Si)
        self.Su_t = copylist(self.Su)
        self.Sa_t = copylist(self.Sa)
        self.Sus_t = copylist(self.Sus)
        self.Sf_t = copylist(self.Sf)
        self.Sr_t = copylist(self.Sr)
        self.Ss_t = self.Ss
        self.trackQ_t = copylist(self.trackQ)    #copylist(self.trackQ)
        self.convQu_t = [copylist(self.convQu[i]) for i in self.Classes]    #copylist(self.convQu)
        
        if self.scalarInput:
            if self.InputSeries == 1:
                self.Precipitation = timeinputscalar(self.precipTss, self.gaugesMap)
                self.PotEvaporation = timeinputscalar(self.evapTss, self.gaugesMap)
                self.Temperature = timeinputscalar(self.tempTss, self.gaugesMap)
            elif self.InputSeries == 2:        
                self.Precipitation = timeinputscalar(self.precipTss2, self.gaugesMap)            
                self.EpDay = timeinputscalar(self.evapTss2, self.gaugesMap)            
                self.Tmean = timeinputscalar(self.tempDMTss, self.gaugesMap)
                self.Rn = timeinputscalar(self.radnTss, self.gaugesMap)
                self.rad_si = timeinputscalar(self.radsTss, self.gaugesMap)
                self.sgamma = timeinputscalar(self.sgammaTss, self.gaugesMap)
                self.vpd = timeinputscalar(self.vpdTss, self.gaugesMap)
                self.wind2m = timeinputscalar(self.windTss, self.gaugesMap)
                self.DS = timeinputscalar(self.daySTss, self.gaugesMap)
                self.DE = timeinputscalar(self.dayETss, self.gaugesMap)
    #            self.LAI = timeinputscalar(self.laiTss,self.gaugesMap)
                self.rst_lai = [timeinputscalar(self.rst_laiTss[i], self.gaugesMap) for i in self.Classes]
    
        else:
           self.Precipitation=cover(self.wf_readmap(self.P_mapstack, 0.0), 0.0)
           self.PotEvaporation=cover(self.wf_readmap(self.PET_mapstack, 0.0), 0.0)
           self.Inflow=pcrut.readmapSave(self.Inflow_mapstack, 0.0)
           if self.ExternalQbase:
               self.Seepage = cover(self.wf_readmap(self.Seepage_mapstack, 0.0), 0.0)
           else:
               self.Seepage=cover(0.0)
           self.Temperature=self.wf_readmap(self.TEMP_mapstack, 0.0)
           
        if self.IRURFR_L:
            self.PotEvaporation = areatotal(self.PotEvaporation * self.percentArea, nominal(self.TopoId))
            self.Precipitation = areatotal(self.Precipitation * self.percentArea, nominal(self.TopoId))
            self.Temperature = areaaverage(self.Temperature * self.percentArea, nominal(self.TopoId))
                
       
    
        for k in self.Classes:   
            
        #SNOW =================================================================================================
        
        #nu nog even niet gecodeerd
        
        
        #INTERCEPTION =========================================================================================
            if self.selectSi[k]:
                eval_str = 'reservoir_Si.{:s}(self, k)'.format(self.selectSi[k])
            else:
                eval_str = 'reservoir_Si.interception_no_reservoir(self, k)'
            eval(eval_str)
        #UNSATURATED ZONE ======================================================================================
            if self.selectSu[k]:
                eval_str = 'reservoir_Su.{:s}(self, k)'.format(self.selectSu[k])
            else:
                eval_str = 'reservoir_Si.unsatZone_no_reservoir(self, k)'
            eval(eval_str)
            
        #COMBINED SATURATED AND UNSATURATED ZONE ========================================================================        
            if self.selectSus[k]:
                eval_str = 'reservoir_Sus.{:s}(self, k)'.format(self.selectSus[k])
                eval(eval_str)
            
        #FAST RUNOFF RESERVOIR ===================================================================================
            if self.selectSf[k]:
                eval_str = 'reservoir_Sf.{:s}(self, k)'.format(self.selectSf[k])
            else:
                eval_str = 'reservoir_Si.fastRunoff_no_reservoir(self, k)'
            eval(eval_str)
                
         #RIPARIAN ZONE RESERVOIR ==================================================================================
            if self.selectSr[k]:
                eval_str = 'reservoir_Sr.{:s}(self, k)'.format(self.selectSr[k])
                eval(eval_str)
        
        #SLOW RUNOFF RESERVOIR ===========================================================================
                       
        
        #TOTAL RUNOFF =============================================================================================
    #    self.Qfcub = (sum([x*y for x,y in zip(self.Qf_,self.percent)]) + sum([x*y for x,y in zip(self.Qo_,self.percent)]) + sum([x*y for x,y in zip(self.Qd_,self.percent)]) + sum([x*y for x,y in zip(self.Qr_,self.percent)]))/ 1000 * self.surfaceArea
        self.Qfcub = (sum([x*y for x,y in zip(self.Qf_,self.percent)]))/ 1000 * self.surfaceArea
        reservoir_Sf.routingQf_combined(self)
        
        if self.selectSs:
            eval_str = 'reservoir_Ss.{:s}(self)'.format(self.selectSs)
        else:
            eval_str = 'reservoir_Ss.groundWater_no_reservoir(self)'
        eval(eval_str)
        # pdb.set_trace()
        # for separate gw reservoir per class    
        # self.Qtlag = self.Qflag_/ self.timestepsecs + sum([x*y for x,y in zip(self.Qs_,self.percent)])/ 1000 * self.surfaceArea/ self.timestepsecs
        # for combinzed gw reservoir    
        self.Qtlag = self.Qflag_/ self.timestepsecs + self.Qs_/ 1000 * self.surfaceArea/ self.timestepsecs
        self.QLagTot = areatotal(self.Qtlag, nominal(self.TopoId))                    #catchment total runoff with looptijd
        
      
        # WATER BALANCE (per reservoir, per cell) ========================================================================================   
        self.QtlagWB = (self.Qtlag / self.surfaceArea) * 1000 * self.timestepsecs
        self.convQuWB = [sum(self.convQu[i]) for i in self.Classes]
        self.convQuWB_t = [sum(self.convQu_t[i]) for i in self.Classes]
        self.trackQWB = (sum(self.trackQ) / self.surfaceArea) * 1000
        self.trackQWB_t = (sum(self.trackQ_t) / self.surfaceArea) * 1000
        self.WB = self.Precipitation - sum(multiply(self.Ei_,self.percent)) - sum(multiply(self.Eu_,self.percent)) - sum(multiply(self.Er_,self.percent)) - self.QtlagWB - sum(multiply(self.Si,self.percent)) + sum(multiply(self.Si_t,self.percent)) - sum(multiply(self.Su,self.percent)) + sum(multiply(self.Su_t,self.percent)) - sum(multiply(self.Sus,self.percent)) + sum(multiply(self.Sus_t,self.percent)) - sum(multiply(self.Sf,self.percent)) + sum(multiply(self.Sf_t,self.percent)) - sum(multiply(self.Sr,self.percent)) + sum(multiply(self.Sr_t,self.percent)) - sum(multiply(self.Ss,self.percent)) + sum(multiply(self.Ss_t,self.percent)) - self.trackQWB + self.trackQWB_t - sum(multiply(self.convQuWB,self.percent)) + sum(multiply(self.convQuWB_t,self.percent))
        
    #    #fuxes and states in m3/h
        self.P = areatotal(self.Precipitation / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Ei = areatotal(sum(multiply(self.Ei_,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))  
        self.Eu = areatotal(sum(multiply(self.Eu_,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Er = areatotal(sum(multiply(self.Er_,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.QtotnoRout = areatotal(self.Qfcub + self.Qs_/ 1000 * self.surfaceArea,nominal(self.TopoId))    
        self.Qtot = self.QLagTot * self.timestepsecs
        self.SiWB = areatotal(sum(multiply(self.Si,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Si_WB = areatotal(sum(multiply(self.Si_t,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.SuWB = areatotal(sum(multiply(self.Su,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Su_WB = areatotal(sum(multiply(self.Su_t,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.SaWB = areatotal(sum(multiply(self.Sa,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Sa_WB = areatotal(sum(multiply(self.Sa_t,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.SusWB = areatotal(sum(multiply(self.Sus,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Sus_WB = areatotal(sum(multiply(self.Sus_t,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.SfWB = areatotal(sum(multiply(self.Sf,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Sf_WB = areatotal(sum(multiply(self.Sf_t,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.SrWB = areatotal(sum(multiply(self.Sr,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Sr_WB = areatotal(sum(multiply(self.Sr_t,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.SsWB = areatotal(self.Ss / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Ss_WB = areatotal(self.Ss_t / 1000 * self.surfaceArea,nominal(self.TopoId))    
        self.convQuWB = areatotal(sum(multiply([sum(self.convQu[i]) for i in self.Classes],self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.convQu_WB = areatotal(sum(multiply([sum(self.convQu_t[i]) for i in self.Classes],self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.trackQWB = areatotal(sum(self.trackQ),nominal(self.TopoId))
        self.trackQ_WB = areatotal(sum(self.trackQ_t),nominal(self.TopoId)) 
        
        #WBtot in m3/s
        self.WBtot = (self.P - self.Ei - self.Eu - self.Er - self.Qtot - self.SiWB + self.Si_WB - self.SuWB + self.Su_WB - self.SaWB + self.Sa_WB - self.SusWB + self.Sus_WB - self.SfWB + self.Sf_WB - self.SrWB + self.Sr_WB - self.SsWB + self.Ss_WB - self.convQuWB +self.convQu_WB -self.trackQWB + self.trackQ_WB) / self.timestepsecs     
        
        
        # SUMMED FLUXES ======================================================================================
        self.sumprecip=self.sumprecip  +  self.Precipitation                     #accumulated rainfall for water balance (m/h)
        self.sumevap=self.sumevap + sum(multiply(self.Ei_,self.percent)) + sum(multiply(self.Eu_,self.percent)) + sum(multiply(self.Er_,self.percent))                         #accumulated evaporation for water balance (m/h)
        try:    
            self.sumpotevap=self.sumpotevap + self.PotEvaporation         #accumulated potential evaporation (m/h)
        except:
            self.sumpotevap=self.EpHour
        self.sumrunoff=self.sumrunoff  + self.Qtlag * 1000 * self.timestepsecs / self.surfaceArea                        #accumulated runoff for water balance (m/h)
        self.sumwb=self.sumwb + self.WB
    
        self.sumE = sum(multiply(self.Ei_,self.percent)) + sum(multiply(self.Eu_,self.percent)) + sum(multiply(self.Er_,self.percent))

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
    configfile="wflow_topoflex.ini"
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

    opts, args = getopt.getopt(argv, 'C:S:T:Ic:s:R:')
    
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

    for o, a in opts:
        if o == '-I': configset(myModel.config, 'model', 'reinit', '1', overwrite=True) 

    dynModelFw._runInitial()
    dynModelFw._runResume()
    dynModelFw._runDynamic(_firstTimeStep,_lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()
    

if __name__ == "__main__":
    main()
