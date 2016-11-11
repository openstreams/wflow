#!/usr/bin/python
# Wflow is Free software, see below:
# 
# Copyright (c) Hylke Beck (JRC) J. Schellekens  2005-2013
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#TODO: split off routing

"""
Run the wflow_hbvl (hbv light) hydrological model..

usage: 
wflow_hbv::
    
      [-h][-v level][-F runinfofile][-L logfile][-C casename][-R runId]
      [-c configfile][-T timesteps][-s seconds][-W][-E][-N][-U discharge]
      [-P parameter multiplication][-X][-l loglevel]
      
-F: if set wflow is expected to be run by FEWS. It will determine
    the timesteps from the runinfo.xml file and save the output initial
    conditions to an alternate location. Also set fewsrun=1 in the .ini file!
    
-f: Force overwrite of existing results    

-T: Set the number of timesteps to run

-N: No lateral flow, use runoff response function to generate fast runoff

-s: Set the model timesteps in seconds

-I: re-initialize the initial model conditions with default

-i: Set input table directory (default is intbl)

-x: run for subcatchment only (e.g. -x 1)

-C: set the name  of the case (directory) to run

-R: set the name runId within the current case

-L: set the logfile

-c: name of wflow the configuration file (default: Casename/wflow_hbv.ini). 

-h: print usage information

-U: The argument to this option should be a .tss file with measured discharge in
    [m^3/s] which the program will use to update the internal state to match
    the measured flow. The number of columns in this file should match the 
    number of gauges in the wflow\_gauges.map file.
    
-u: list of gauges/columns to use in update. Format:
    -u [1 , 4 ,13]
    The above example uses column 1, 4 and 13
    
-P: set parameter change string (e.g: -P 'self.FC = self.FC * 1.6') for non-dynamic variables
    
-p: set parameter change string (e.g: -P 'self.Precipitation = self.Precipitation * 1.11') for
    dynamic variables

-l: loglevel (most be one of DEBUG, WARNING, ERROR)

-X overwrites the initial values at the end of each timestep


"""

import numpy
import os
import os.path
import shutil, glob
import getopt

from wflow.wf_DynamicFramework import *
from wflow.wf_DynamicFramework import *
from wflow.wflow_adapt import *
from wflow_adapt import *    
        
#import scipy
#import pcrut



wflow = "wflow_hbv"


#: columns used in updating
updateCols = [] #: columns used in updating
""" Column used in updating """


def usage(*args):
    """
    Print usage information
    
    -  *args: command line arguments given
    """
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)

class WflowModel(DynamicModel):

  """
  The user defined model class.

  """
  

  def __init__(self, cloneMap,Dir,RunDir,configfile):
      DynamicModel.__init__(self)
      self.caseName = os.path.abspath(Dir)
      self.clonemappath = os.path.join(os.path.abspath(Dir),"staticmaps",cloneMap)
      setclone(self.clonemappath)
      self.runId = RunDir
      self.Dir = os.path.abspath(Dir)
      self.configfile = configfile
      self.SaveDir = os.path.join(self.Dir,self.runId)


  def stateVariables(self):
      """ 
      returns a list of state variables that are essential to the model. 
      This list is essential for the resume and suspend functions to work.
      
      This function is specific for each model and **must** be present.

     :var self.DrySnow: Snow pack [mm]
     :var self.FreeWater:  Available free water [mm]
     :var self.UpperZoneStorage: Water in the upper zone [mm]
     :var self.LowerZoneStorage: Water in the lower zone [mm]
     :var self.SoilMoisture: Soil moisture [mm]

      """
      states = ['FreeWater', 'SoilMoisture',
                 'UpperZoneStorage',
                 'LowerZoneStorage',
                 'DrySnow']
      
      return states
      
      
    # The following are made to better connect to deltashell/openmi
  def supplyCurrentTime(self):
      """
      gets the current time in seconds after the start of the run
      
      Ouput:
          - time in seconds since the start of the model run
      """
      return self.currentTimeStep() * int(configget(self.config,'model','timestepsecs','86400'))
  
  def parameters(self):
    """
    Define all model parameters here that the framework should handle for the model
    See wf_updateparameters and the parameters section of the ini file
    If you use this make sure to all wf_updateparameters at the start of the dynamic section
    and at the start/end of the initial section
    """
    modelparameters = []

    #Static model parameters e.g.
    #modelparameters.append(self.ParamType(name="RunoffGeneratingGWPerc",stack="intbl/RunoffGeneratingGWPerc.tbl",type="static",default=0.1))

    # Meteo and other forcing

    modelparameters.append(self.ParamType(name="Precipitation",stack="inmaps/P",type="timeseries",default=0.0,verbose=False,lookupmaps=[]))
    modelparameters.append(self.ParamType(name="PotEvaporation",stack="inmaps/PET",type="timeseries",default=0.0,verbose=False,lookupmaps=[]))
    modelparameters.append(self.ParamType(name="Temperature",stack="inmaps/TEMP",type="timeseries",default=10.0,verbose=False,lookupmaps=[]))




    return modelparameters


  def suspend(self):
    """
      Suspends the model to disk. All variables needed to restart the model
      are saved to disk as pcraster maps. Use resume() to re-read them
    """
    
    
    self.logger.info("Saving initial conditions...")
    self.wf_suspend(os.path.join(self.SaveDir,"outstate"))
   
    if self.OverWriteInit:            
        self.logger.info("Saving initial conditions over start conditions...")
        self.wf_suspend(os.path.join(self.SaveDir,"instate"))


    if self.fewsrun:
        self.logger.info("Saving initial conditions for FEWS...")
        self.wf_suspend(os.path.join(self.Dir, "outstate"))
        

  def initial(self):
      
    """
    Initial part of the model, executed only once. Reads all static model
    information (parameters) and sets-up the variables used in modelling.
    
    *HBV Soil*
    
    :var FC.tbl: Field Capacity (260.0) [mm]  
    :var BETA.tbl: exponent in soil runoff generation equation (1.8)  [-]
    :var LP.tbl: fraction of Fieldcapacity below which actual evaporation=potential evaporation (0.53000)
    :var K2.tbl: Recession constant baseflow (0.02307)
  
    *If SetKquickFlow is set to 1*
 
    :var K1.tbl: (0.09880)
    :var SUZ.tbl: Level over wich K0 is used (100.0) 
    :var K0.tbl: (0.3)
    

    :var PERC.tbl: Percolation from Upper to Lowerzone (0.4000)  [mm/day]
    :var CFR.tbl: Refreezing efficiency constant in refreezing of freewater in snow (0.05000)
    :var PCORR.tbl: Correction factor for precipitation (1.0)
    :var RFCF.tbl: Correction factor for rainfall (1.0)      
    :var SFCF.tbl: Correction factor for snowfall(1.0)     
    :var CEVPF.tbl: Correction factor for potential evaporation (1.0)
    :var EPF.tbl: Exponent of correction factor for evaporation on days with precipitation(0.0)
    :var ECORR.tbl: Evap correction (1.0)
   
    
    *Snow modelling parameters*
    
    :var TTI.tbl: critical temperature for snowmelt and refreezing  (1.000) [oC]
    :var TT.tbl: defines interval in which precipitation falls as rainfall and snowfall (-1.41934) [oC]
    :var CFMAX.tbl: meltconstant in temperature-index ( 3.75653) [-]
    :var WHC.tbl: fraction of Snowvolume that can store water (0.1) [-]
 
    
    """
    global statistics
    global multpars
    global updateCols    
    
    setglobaloption("unittrue")
    
    
    self.thestep = scalar(0)

    #: files to be used in case of timesries (scalar) input to the model
    
    #: name of the tss file with precipitation data ("../intss/P.tss")


    self.logger.info("running for " + str(self.nrTimeSteps()) + " timesteps") 

    
        # Set and get defaults from ConfigFile here ###################################

    self.interpolMethod = configget(self.config,"model","InterpolationMethod","inv")
    self.reinit = int(configget(self.config,"run","reinit","0"))
    self.fewsrun = int(configget(self.config,"run","fewsrun","0"))
    self.OverWriteInit = int(configget(self.config,"model","OverWriteInit","0"))

    self.intbl = configget(self.config,"model","intbl","intbl")
    self.timestepsecs = int(configget(self.config,"model","timestepsecs","86400"))
    self.P_style = int(configget(self.config,"model","P_style","1"))
    self.PET_style = int(configget(self.config,"model","PET_style","1"))
    self.TEMP_style = int(configget(self.config,"model","TEMP_style","1"))

    sizeinmetres = int(configget(self.config,"layout","sizeinmetres","0"))

    # static maps to use (normally default)
    wflow_subcatch = configget(self.config,"model","wflow_subcatch","staticmaps/wflow_subcatch.map")
    wflow_dem  = configget(self.config,"model","wflow_dem","staticmaps/wflow_dem.map")
    wflow_landuse  = configget(self.config,"model","wflow_landuse","staticmaps/wflow_landuse.map")
    wflow_soil  = configget(self.config,"model","wflow_soil","staticmaps/wflow_soil.map")
    wflow_gauges  = configget(self.config,"model","wflow_gauges","staticmaps/wflow_gauges.map")

    # 2: Input base maps ########################################################  
    subcatch = ordinal(self.wf_readmap(os.path.join(self.Dir,wflow_subcatch),0.0,fail=True))  # Determines the area of calculations (all cells > 0)
    subcatch = ifthen(subcatch > 0, subcatch)

    self.Altitude=self.wf_readmap(os.path.join(self.Dir,wflow_dem),0.0,fail=True) * scalar(defined(subcatch)) #: The digital elevation map (DEM)
    self.TopoId=self.wf_readmap(os.path.join(self.Dir, wflow_subcatch),0.0,fail=True)        #: Map define the area over which the calculations are done (mask)

    # read landuse and soilmap and make sure there are no missing points related to the
    # subcatchment map. Currently sets the lu and soil type  type to 1
    self.LandUse=self.wf_readmap(os.path.join(self.Dir , wflow_landuse),0.0,fail=True)#: Map with lan-use/cover classes
    self.LandUse=cover(self.LandUse,nominal(ordinal(subcatch) > 0))
    self.Soil=self.wf_readmap(os.path.join(self.Dir , wflow_soil),0.0,fail=True)#: Map with soil classes
    self.Soil=cover(self.Soil,nominal(ordinal(subcatch) > 0))
    self.OutputLoc=self.wf_readmap(os.path.join(self.Dir , wflow_gauges),0.0,fail=True)  #: Map with locations of output gauge(s)

    
     # Temperature correction per cell to add
    self.TempCor=self.wf_readmap(os.path.join(self.Dir , configget(self.config,"model","TemperatureCorrectionMap","staticmap/swflow_tempcor.map")),0.0)
    self.OutputId=self.wf_readmap(os.path.join(self.Dir , wflow_subcatch),0.0,fail=True)       # location of subcatchment
  
    self.ZeroMap=0.0*scalar(defined(self.Altitude))                    #map with only zero's
  
    # 3: Input time series ###################################################
    self.P_mapstack=self.Dir + configget(self.config,"inputmapstacks","Precipitation","/inmaps/P") # timeseries for rainfall
    self.PET_mapstack=self.Dir + configget(self.config,"inputmapstacks","EvapoTranspiration","/inmaps/PET") # timeseries for rainfall"/inmaps/PET"          # potential evapotranspiration
    self.TEMP_mapstack=self.Dir + configget(self.config,"inputmapstacks","Temperature","/inmaps/TEMP") # timeseries for rainfall "/inmaps/TEMP"          # global radiation
    # For in memory override:
    self.P = self.ZeroMap
    self.PET = self.ZeroMap
    self.TEMP = self.ZeroMap
    # Set static initial values here #########################################

    self.Latitude  =  ycoordinate(boolean(self.Altitude))
    self.Longitude =  xcoordinate(boolean(self.Altitude))
  
    self.logger.info("Linking parameters to landuse, catchment and soil...")

    # TODO: Set default properly
    # TODO: make unit test, running model
    #HBV Soil params
    # + BETA.tif
    # + CFMAX.tif
    # + CFR.tif
    # + CWH.tif -> WHC.tif
    # + FC.tif
    # + K0.tif
    # + K1.tif
    # + K2.tif
    # + LP.tif
    # MAXBAS.tif
    # + PCORR.tif
    # + PERC.tif
    # + SFCF.tif
    # + TT.tif
    # + UZL.tif 

    self.FC = self.readtblDefault(self.Dir + "/" + self.intbl + "/FC.tbl",self.LandUse,subcatch,self.Soil,260.0)

    self.BETA= self.readtblDefault(self.Dir + "/" + self.intbl + "/BETA.tbl",self.LandUse,subcatch,self.Soil,1.8)  # exponent in soil runoff generation equation
    self.K0= self.readtblDefault(self.Dir + "/" + self.intbl + "/K0.tbl",self.LandUse,subcatch,self.Soil, 0.02307) # Recession constant baseflow   #K4=0.07; BASEFLOW:LINEARRESERVOIR
    self.K1= self.readtblDefault(self.Dir + "/" + self.intbl + "/K2.tbl",self.LandUse,subcatch,self.Soil, 0.02307) # Recession constant baseflow   #K4=0.07; BASEFLOW:LINEARRESERVOIR
    self.K2= self.readtblDefault(self.Dir + "/" + self.intbl + "/K2.tbl",self.LandUse,subcatch,self.Soil, 0.02307) # Recession constant baseflow   #K4=0.07; BASEFLOW:LINEARRESERVOIR
    self.LP= self.readtblDefault(self.Dir + "/" + self.intbl + "/LP.tbl",self.LandUse,subcatch,self.Soil, 0.4000) # percolation from Upper to Lowerzone (mm/day)
    self.UZL= self.readtblDefault(self.Dir + "/" + self.intbl + "/UZL.tbl",self.LandUse,subcatch,self.Soil, 0.4000) # percolation from Upper to Lowerzone (mm/day)
    self.PERC= self.readtblDefault(self.Dir + "/" + self.intbl + "/PERC.tbl",self.LandUse,subcatch,self.Soil, 0.4000) # percolation from Upper to Lowerzone (mm/day)
    self.CFR=self.readtblDefault(self.Dir + "/" + self.intbl + "/CFR.tbl",self.LandUse,subcatch,self.Soil, 0.05000)        # refreezing efficiency constant in refreezing of freewater in snow
    self.PCORR=self.readtblDefault(self.Dir + "/" + self.intbl + "/PCORR.tbl",self.LandUse,subcatch,self.Soil, 1.0)      # correction factor for precipitation
    self.SFCF=self.readtblDefault(self.Dir + "/" + self.intbl + "/SFCF.tbl",self.LandUse,subcatch,self.Soil, 1.0)     # correction factor for snowfall
    self.CFMAX= self.readtblDefault(self.Dir + "/" + self.intbl + "/CFMAX.tbl",self.LandUse,subcatch,self.Soil, 2.0)        # maximum capillary rise from runoff response routine to soil moisture routine
    self.WHC= self.readtblDefault(self.Dir + "/" + self.intbl + "/WHC.tbl",self.LandUse,subcatch,self.Soil, 2.0)        # maximum capillary rise from runoff response routine to soil moisture routine
    self.TTI=self.readtblDefault(self.Dir + "/" + self.intbl + "/TTI.tbl" ,self.LandUse,subcatch,self.Soil,1.0)
    self.TT=self.readtblDefault(self.Dir + "/" + self.intbl + "/TT.tbl" ,self.LandUse,subcatch,self.Soil,-1.41934)
    #Cfmax = 3.75653 # meltconstant in temperature-index
    self.RFCF=self.readtblDefault(self.Dir + "/" + self.intbl + "/RFCF.tbl",self.LandUse,subcatch,self.Soil,1.0)      # correction factor for rainfall
    self.CEVPF= self.readtblDefault(self.Dir + "/" + self.intbl + "/CEVPF.tbl",self.LandUse,subcatch,self.Soil, 1.0)   # correction factor for potential evaporation (1.15 in in forested areas )
    self.EPF= self.readtblDefault(self.Dir + "/" + self.intbl + "/EPF.tbl",self.LandUse,subcatch,self.Soil, 0.0)    # exponent of correction factor for evaporation on days with precipitation
    self.ECORR= self.readtblDefault(self.Dir + "/" + self.intbl + "/ECORR.tbl",self.LandUse,subcatch,self.Soil, 1.0)    # evap correction

    # Determine real slope and cell length
    self.xl,self.yl,self.reallength = pcrut.detRealCellLength(self.ZeroMap,sizeinmetres)

    # Multiply parameters with a factor (for calibration etc) -P option in command line
    self.wf_multparameters()

    self.SnowWater = self.ZeroMap

    # Initializing of variables
    self.logger.info("Initializing of model variables..")
    self.QMMConv = self.timestepsecs/(self.reallength * self.reallength * 0.001) #m3/s --> mm
    self.ToCubic = (self.reallength * self.reallength * 0.001) / self.timestepsecs # m3/s

    self.FieldCapacity=self.FC                    #: total water holding capacity of the soil
    self.Treshold=self.LP*self.FieldCapacity      # Threshold soilwaterstorage above which AE=PE

    self.logger.info("End of initial section.")


  def default_summarymaps(self):
      """
      Returns a list of default summary-maps at the end of a run.
      This is model specific. You can also add them to the [summary]section of the ini file but stuff
      you think is crucial to the model should be listed here

       Example:

      """
      lst = ['self.csize','self.upsize','self.TTI','self.TT','self.WHC',
             'self.Slope','self.N','self.xl','self.yl','self.reallength','self.DCL','self.Bw',]

      return lst

  def resume(self):
    """ read initial state maps (they are output of a previous call to suspend()) """
    
    if self.reinit == 1:
        self.logger.info("Setting initial conditions to default (zero!)")
        self.FreeWater =  cover(0.0) #: Water on surface (state variable [mm])
        self.SoilMoisture =  self.FC #: Soil moisture (state variable [mm])
        self.UpperZoneStorage = 0.2 * self.FC #: Storage in Upper Zone (state variable [mm])
        self.LowerZoneStorage = 1.0/(3.0 * self.K2) #: Storage in Uppe Zone (state variable [mm])
        self.DrySnow=cover(0.0) #: Snow amount (state variable [mm])
    else:
        self.wf_resume(os.path.join(self.Dir, "instate"))

    self.initstorage=self.FreeWater + self.DrySnow + self.SoilMoisture + self.UpperZoneStorage + self.LowerZoneStorage

    
    
  def dynamic(self):
    
    """    
    Below a list of variables that can be save to disk as maps or as 
    timeseries (see ini file for syntax):
        
    *Dynamic variables*
    
    :var self.Snow: Snow depth [mm]
    :var self.SnowWater: water content of the snow [mm]
    :var self.LowerZoneStorage: water content of the lower zone [mm]
    :var self.UpperZoneStorage: water content of the Upper zone [mm]
    :var self.Q2: Specific runoff (baseflow part) per cell [mm]
    :var self.Percolation: actual percolation to the lower zone [mm]
    :var self.SoilMoisture: actual soil moisture [mm]
    :var self.Q1: specific runoff (quickflow part) [mm]
    :var self.Q0: specific runoff (quickflow), If K upper zone is precalculated [mm]

    
    *Static variables*
        
    :var self.Altitude: The altitude of each cell [m]
    :var self.ToCubic: Mutiplier to convert mm to m^3/s for fluxes
    """

    self.logger.debug("Step: " + str(int(self.currentStep)) + "/" + str(int(self._d_nrTimeSteps)))
    self.thestep = self.thestep + 1

    self.wf_updateparameters() # read forcing an dynamic parameters


    # Apply correction factor to precipitation
    self.Precipitation =  self.PCORR * self.Precipitation
    self.Temperature=cover(self.wf_readmap(self.TEMP_mapstack,10.0),10.0)
    self.Temperature = self.Temperature + self.TempCor

    # Multiply input parameters with a factor (for calibration etc) -p option in command line (no also in ini)

    self.wf_multparameters()

    RainFrac=ifthenelse(1.0*self.TTI == 0.0,ifthenelse(self.Temperature <= self.TT,scalar(0.0),scalar(1.0)),min((self.Temperature-(self.TT-self.TTI/2.0))/self.TTI,scalar(1.0)))
    RainFrac=max(RainFrac,scalar(0.0))               #fraction of precipitation which falls as rain
    SnowFrac=1.0-RainFrac                    #fraction of self.Precipitation which falls as snow
    self.Precipitation=self.SFCF*SnowFrac*self.Precipitation+self.RFCF*RainFrac*self.Precipitation # different correction for rainfall and snowfall

    self.PotEvaporation=exp(-self.EPF*self.Precipitation)*self.ECORR * self.PotEvaporation  # correction for potential evaporation on wet days
    self.PotEvaporation=self.CEVPF*self.PotEvaporation  # Correct per landuse

    SnowFall=SnowFrac*self.Precipitation  #: snowfall depth
    RainFall=RainFrac*self.Precipitation  #: rainfall depth
    PotSnowMelt=ifthenelse(self.Temperature > self.TT,self.CFMAX*(self.Temperature-self.TT),scalar(0.0)) #Potential snow melt, based on temperature
    PotRefreezing=ifthenelse(self.Temperature < self.TT, self.CFMAX*self.CFR*(self.TT-self.Temperature),0.0)    #Potential refreezing, based on temperature

    Refreezing=ifthenelse(self.Temperature < self.TT,min(PotRefreezing,self.FreeWater),0.0)   	#actual refreezing
    self.SnowMelt=min(PotSnowMelt,self.DrySnow)          #actual snow melt
    self.DrySnow=self.DrySnow+SnowFall+Refreezing-self.SnowMelt     #dry snow content
    self.FreeWater=self.FreeWater-Refreezing               #free water content in snow
    MaxFreeWater=self.DrySnow*self.WHC
    self.FreeWater=self.FreeWater+self.SnowMelt+RainFall
    InSoil = max(self.FreeWater-MaxFreeWater,0.0)   #abundant water in snow pack which goes into soil
    self.FreeWater=self.FreeWater-InSoil


    # Soil and evaporation
    soil_wetness = (self.SoilMoisture/self.FC) ** self.BETA
    soil_wetness = max(min(soil_wetness, 1.0),0.0)
    recharge = (self.Precipitation+InSoil) * soil_wetness
    self.SoilMoisture = self.SoilMoisture+self.Precipitation+InSoil-recharge
    excess = self.SoilMoisture-self.FC
    excess = max(excess,0.0)
    self.SoilMoisture = self.SoilMoisture-excess
    evapfactor = self.SoilMoisture / (self.LP*self.FC)
    evapfactor = min(max(evapfactor,0.0), 1.0)
    #----------------
    self.ActEvap = self.PotEvaporation*evapfactor
    self.ActEvap = min(self.SoilMoisture, self.ActEvap)
    self.SoilMoisture = self.SoilMoisture-self.ActEvap

    # Groundwater boxes
    self.UpperZoneStorage = self.UpperZoneStorage+recharge+excess
    self.actPERC = min(self.UpperZoneStorage, self.PERC)
    self.UpperZoneStorage = self.UpperZoneStorage-self.actPERC
    self.Q0 = self.K0 * max(self.UpperZoneStorage-self.UZL, 0.0)
    self.UpperZoneStorage = self.UpperZoneStorage-self.Q0
    self.Q1 = self.K1*self.UpperZoneStorage
    self.UpperZoneStorage = self.UpperZoneStorage-self.Q1
    self.LowerZoneStorage = self.LowerZoneStorage+self.actPERC
    self.Q2 = self.K2*self.LowerZoneStorage
    self.LowerZoneStorage = self.LowerZoneStorage-self.Q2

    DirectRunoffStorage= self.Q0 + self.Q1 + self.Q2

    self.InwaterMM=max(0.0,DirectRunoffStorage)
    self.Inwater=self.InwaterMM * self.ToCubic
    self.QuickFlowCubic = (self.Q0 + self.Q1) * self.ToCubic
    self.BaseFlowCubic = self.Q2 * self.ToCubic
    




# The main function is used to run the program from the command line

def main(argv=None):  
    """
    Perform command line execution of the model.
    """      
    global multpars
    global updateCols
    caseName = "default_hbv"
    runId = "run_default"
    configfile="wflow_hbvl.ini"
    LogFileName="wflow.log"    
    _lastTimeStep = 0
    _firstTimeStep = 1
    fewsrun=False
    runinfoFile="runinfo.xml"
    timestepsecs=86400
    wflow_cloneMap = 'wflow_subcatch.map'
    NoOverWrite=1
    loglevel = logging.DEBUG
    
    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return     
    
    ## Main model starts here
    ########################################################################
    try:
        opts, args = getopt.getopt(argv, 'c:QXS:F:hC:Ii:T:R:u:s:P:p:Xx:U:fl:L:')
    except getopt.error, msg:
        pcrut.usage(msg)
    
    for o, a in opts:
        if o == '-F': 
            runinfoFile = a
            fewsrun = True

        if o == '-C': caseName = a
        if o == '-R': runId = a
        if o == '-L': LogFileName = a 
        if o == '-l': exec "loglevel = logging." + a           
        if o == '-c': configfile = a
        if o == '-s': timestepsecs = int(a)
        if o == '-T': _lastTimeStep=int(a)
        if o == '-S': _firstTimeStep=int(a)
        if o == '-h': usage()
        if o == '-f': NoOverWrite = 0
        

     
    if fewsrun: 
        ts = getTimeStepsfromRuninfo(runinfoFile,timestepsecs)
        starttime = getStartTimefromRuninfo(runinfoFile)
        if (ts):
            _lastTimeStep =  ts# * 86400/timestepsecs
            _firstTimeStep = 1 
        else:
            print "Failed to get timesteps from runinfo file: " + runinfoFile
            exit(2)
    else:
        starttime = dt.datetime(1990,01,01)
       
    if _lastTimeStep < _firstTimeStep:
        print "The starttimestep (" + str(_firstTimeStep) +") is smaller than the last timestep (" + str(_lastTimeStep) + ")"
        usage()
 
    myModel = WflowModel(wflow_cloneMap, caseName,runId,configfile)
    dynModelFw = wf_DynamicFramework(myModel, _lastTimeStep,firstTimestep=_firstTimeStep,datetimestart=starttime)
    dynModelFw.createRunId(NoOverWrite=NoOverWrite,logfname=LogFileName,level=loglevel,doSetupFramework=False)

    for o, a in opts:
        if o == '-P':
            left = a.split('=')[0]
            right = a.split('=')[1]
            configset(myModel.config,'variable_change_once',left,right,overwrite=True)
        if o == '-p':
            left = a.split('=')[0]
            right = a.split('=')[1]
            configset(myModel.config,'variable_change_timestep',left,right,overwrite=True)
        if o == '-X': configset(myModel.config,'model','OverWriteInit','1',overwrite=True)
        if o == '-I': configset(myModel.config,'model','reinit','1',overwrite=True)
        if o == '-i': configset(myModel.config,'model','intbl',a,overwrite=True)
        if o == '-s': configset(myModel.config,'model','timestepsecs',a,overwrite=True)
        if o == '-x': configset(myModel.config,'model','sCatch',a,overwrite=True)
        if o == '-c': configset(myModel.config,'model','configfile', a,overwrite=True)
        if o == '-M': configset(myModel.config,'model','MassWasting',"0",overwrite=True)
        if o == '-Q': configset(myModel.config,'model','ExternalQbase','1',overwrite=True)
        if o == '-U':
            configset(myModel.config,'model','updateFile',a,overwrite=True)
            configset(myModel.config,'model','updating',"1",overwrite=True)
        if o == '-u':
            exec "zz =" +  a
            updateCols = zz

    dynModelFw.setupFramework()
    dynModelFw._runInitial()
    dynModelFw._runResume()
    dynModelFw._runDynamic(_firstTimeStep,_lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()
    
    
    
    os.chdir("../../")


if __name__ == "__main__":
    main()
