#!/usr/bin/python
# Wflow is Free software, see below:
# 
# Copyright (c) J. Schellekens/Deltares 2005-2013
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

#TODO: remove dots in dynamic phase (default pcraster progress (how?)
#TODO: formal test runs against SMHI model
#TODO: split off routing

# $Rev:: 904           $:  Revision of last commit
# $Author:: schelle    $:  Author of last commit
# $Date:: 2014-01-13 1#$:  Date of last commit
"""
Run the wflow_hbv hydrological model..

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
    [m^3/s] which the progam will use to update the internal state to match 
    the measured flow. The number of columns in this file should match the 
    number of gauges in the wflow\_gauges.map file.
    
-u: list of gauges/columns to use in update. Format:
    -u [1 , 4 ,13]
    The above example uses column 1, 4 and 13
    
-P: set parameter multiply dictionary (e.g: -P {'self.FirstZoneDepth' : 1.2}
    to increase self.FirstZoneDepth by 20%, multiply with 1.2)
    
-p: set input parameter (dynamic, e.g. precip) multiply dictionary 
    (e.g: -p {'self.Precipitation' : 1.2} to increase Precipitation 
    by 20%, multiply with 1.2)    
    
-l: loglevel (most be one of DEBUG, WARNING, ERROR)

-X overwrites the initial values at the end of each timestep

$Author: schelle $
$Id: wflow_hbv.py 904 2014-01-13 14:39:24Z schelle $
$Rev: 904 $
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
#import pcrut



wflow = "wflow_hbv: "
wflowVersion = "$Revision: 904 $  $Date: 2014-01-13 15:39:24 +0100 (Mon, 13 Jan 2014) $" 
"""revision of the model"""


#: columns used in updating
updateCols = [] #: columns used in updating
""" Column sused in updating """


multpars = {} #: Dictionary with parameters and multipliers (static) (used in calibration)
multdynapars = {} #: Dictionary with parameters and multipliers (dynamic) (used in calibration)



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
  The user defined model class. All maps are defined here for documentation 
  purposes
    
  
  """
  
  
  def __init__(self, cloneMap,Dir,RunDir,configfile):
      DynamicModel.__init__(self)   
      setclone(Dir + "/staticmaps/" + cloneMap)
      self.runId=RunDir      
      self.caseName=Dir
      self.Dir = Dir + "/"
      self.configfile = configfile
      self.SaveDir = self.Dir + "/" + self.runId + "/"

      
        
  
      
  def updateRunOff(self):
      """
      Updates the kinematic wave reservoir
      """
        
      self.WaterLevel=(self.Alpha*pow(self.SurfaceRunoff,self.Beta))/self.Bw      
      # wetted perimeter (m)
      P=self.Bw+(2*self.WaterLevel)
      # Alpha  
      self.Alpha=self.AlpTerm*pow(P,self.AlpPow)
      self.OldKinWaveVolume = self.KinWaveVolume
      self.KinWaveVolume = self.WaterLevel * self.Bw * self.DCL


  def stateVariables(self):
      """ 
      returns a list of state variables that are essential to the model. 
      This list is essential for the resume and suspend functions to work.
      
      This function is specific for each model and **must** be present.

     :var self.SurfaceRunoff: Surface runoff in the kin-wave resrvoir [m^3/s]
     :var self.WaterLevel: Water level in the kin-wave resrvoir [m]
     :var self.DrySnow: Snow pack [mm]
     :var self.FreeWater:  Available free water [mm]
     :var self.UpperZoneStorage: Water in the upper zone [mm]
     :var self.LowerZoneStorage: Water in the lower zone [mm]
     :var self.SoilMoisture: Soil moisture [mm]
     :var self.InterceptionStorage: Amount of water on the Canopy [mm]
      
      """
      states = ['FreeWater', 'SoilMoisture',
                 'UpperZoneStorage',
                 'LowerZoneStorage',
                 'InterceptionStorage',
                 'SurfaceRunoff',
                 'WaterLevel',
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
  


  def suspend(self):
    """
      Suspends the model to disk. All variables needed to restart the model
      are saved to disk as pcraster maps. Use resume() to re-read them
    """
    
    
    self.logger.info("Saving initial conditions...")
    self.wf_suspend(self.SaveDir + "/outstate/")
   
    if self.OverWriteInit:            
        self.logger.info("Saving initial conditions over start conditions...")
        self.wf_suspend(self.SaveDir + "/instate/")


    if self.fewsrun:
        self.logger.info("Saving initial conditions for FEWS...")
        self.wf_suspend(self.Dir + "/outstate/")
        
    report(self.sumprecip,self.SaveDir + "/outsum/sumprecip.map")
    report(self.sumevap,self.SaveDir + "/outsum/sumevap.map")
    report(self.sumpotevap,self.SaveDir + "/outsum/sumpotevap.map")
    report(self.sumtemp,self.SaveDir + "/outsum/sumtemp.map")
    report(self.sumlevel,self.SaveDir + "/outsum/sumlevel.map")
    report(self.sumrunoff/catchmenttotal(1,self.TopoLdd),self.SaveDir + "/outsum/sumrunoff.map")

    report(self.suminflow,self.SaveDir + "/outsum/suminflow.map")
      
  def initial(self):
      
    """
    Initial part of the model, executed only once. Reads all static model
    information (parameters) and sets-up the variables used in modelling.
    
    *HBV Soil*
    
    :var FC.tbl: Field Capacity (260.0) [mm]  
    :var BetaSeepage.tbl: exponent in soil runoff generation equation (1.8)  [-] 
    :var LP.tbl: fraction of Fieldcapacity below which actual evaporation=potential evaporation (0.53000)
    :var K4.tbl: Recession constant baseflow (0.02307)
  
    *If SetKquickFlow is set to 1*
 
    :var KQuickFlow.tbl: (0.09880)
    :var SUZ.tbl: Level over wich K0 is used (100.0) 
    :var K0.tbl: (0.3)
    
    *If SetKquickFlow is set to 0*
    
    :var KHQ.tbl: recession rate at flow HQ (0.09880)
    :var HQ.tbl: high flow rate HQ for which recession rate of upper reservoir is known (3.27000)
    :var AlphaNL.tbl: measure of non-linearity of upper reservoir (1.1)
        
    :var PERC.tbl: Percolation from Upper to Lowerzone (0.4000)  [mm/day]
    :var CFR.tbl: Refreezing efficiency constant in refreezing of freewater in snow (0.05000)
    :var Pcorr.tbl: Correction factor for precipitation (1.0)      
    :var RFCF.tbl: Correction factor for rainfall (1.0)      
    :var SFCF.tbl: Correction factor for snowfall(1.0)     
    :var Cflux.tbl: Maximum capillary rise from runoff response routine to soil moisture routine     (2.0)        
    :var ICF.tbl: Maximum interception storage (in forested AND non-forested areas) (2.0)
    :var CEVPF.tbl: Correction factor for potential evaporation (1.0)
    :var EPF.tbl: Exponent of correction factor for evaporation on days with precipitation(0.0)
    :var ECORR.tbl: Evap correction (1.0)
   
    
    *Snow modelling parameters*
    
    :var TTI.tbl: critical temperature for snowmelt and refreezing  (1.000) [oC]
    :var TT.tbl: defines interval in which precipitation falls as rainfall and snowfall (-1.41934) [oC]
    :var Cfmax.tbl: meltconstant in temperature-index ( 3.75653) [-]
    :var WHC.tbl: fraction of Snowvolume that can store water (0.1) [-]
 
    
    """
    global statistics
    global multpars
    global updateCols    
    
    setglobaloption("unittrue")
    
    
    self.thestep = scalar(0)

    #: files to be used in case of timesries (scalar) input to the model
    
    #: name of the tss file with precipitation data ("../intss/P.tss")
    self.precipTss = "../intss/P.tss" 
    self.evapTss="../intss/PET.tss" #: name of the tss file with potential evap data ("../intss/PET.tss")
    self.tempTss="../intss/T.tss" #: name of the tss file with temperature  data ("../intss/T.tss")
    self.inflowTss="../intss/Inflow.tss" #: NOT TESTED name of the tss file with inflow data ("../intss/Inflow.tss")
    self.SeepageTss="../intss/Seepage.tss" #: NOT TESTED name of the tss file with seepage data ("../intss/Seepage.tss")"    
        

    self.logger.info("running for " + str(self.nrTimeSteps()) + " timesteps") 

    
        # Set and get defaults from ConfigFile here ###################################
    self.scalarInput = int(configget(self.config,"model","ScalarInput","0"))
    self.Tslice = int(configget(self.config,"model","Tslice","1"))
    self.interpolMethod = configget(self.config,"model","InterpolationMethod","inv")
    self.reinit = int(configget(self.config,"model","reinit","0"))
    self.fewsrun = int(configget(self.config,"model","fewsrun","0"))
    self.OverWriteInit = int(configget(self.config,"model","OverWriteInit","0"))
    self.updating = int(configget(self.config,"model","updating","0"))
    self.updateFile = configget(self.config,"model","updateFile","no_set")

    self.sCatch = int(configget(self.config,"model","sCatch","0"))
    self.intbl = configget(self.config,"model","intbl","intbl")
    self.timestepsecs = int(configget(self.config,"model","timestepsecs","86400"))
    self.P_style = int(configget(self.config,"model","P_style","1"))
    self.PET_style = int(configget(self.config,"model","PET_style","1"))
    self.TEMP_style = int(configget(self.config,"model","TEMP_style","1"))
    
    
    self.modelSnow = int(configget(self.config,"model","ModelSnow","1"))
    sizeinmetres = int(configget(self.config,"layout","sizeinmetres","0"))
    alf = float(configget(self.config,"model","Alpha","60"))
    Qmax = float(configget(self.config,"model","AnnualDischarge","300"))
    self.UpdMaxDist =float(configget(self.config,"model","UpdMaxDist","100"))
    self.MaxUpdMult =float(configget(self.config,"model","MaxUpdMult","1.3"))
    self.MinUpdMult =float(configget(self.config,"model","MinUpdMult","0.7"))    
    self.UpFrac =float(configget(self.config,"model","UpFrac","0.8"))    
    self.ExternalQbase=int(configget(self.config,'model','ExternalQbase','0'))
    self.SetKquickFlow=int(configget(self.config,'model','SetKquickFlow','0'))

    # static maps to use (normally default)
    wflow_subcatch = configget(self.config,"model","wflow_subcatch","/staticmaps/wflow_subcatch.map")
    wflow_dem  = configget(self.config,"model","wflow_dem","/staticmaps/wflow_dem.map")
    wflow_ldd = configget(self.config,"model","wflow_ldd","/staticmaps/wflow_ldd.map")
    wflow_river  = configget(self.config,"model","wflow_river","/staticmaps/wflow_river.map")
    wflow_riverlength  = configget(self.config,"model","wflow_riverlength","/staticmaps/wflow_riverlength.map")
    wflow_riverlength_fact  = configget(self.config,"model","wflow_riverlength_fact","/staticmaps/wflow_riverlength_fact.map")
    wflow_landuse  = configget(self.config,"model","wflow_landuse","/staticmaps/wflow_landuse.map")
    wflow_soil  = configget(self.config,"model","wflow_soil","/staticmaps/wflow_soil.map")
    wflow_gauges  = configget(self.config,"model","wflow_gauges","/staticmaps/wflow_gauges.map")
    wflow_inflow  = configget(self.config,"model","wflow_inflow","/staticmaps/wflow_inflow.map")
    wflow_mgauges  = configget(self.config,"model","wflow_mgauges","/staticmaps/wflow_mgauges.map")
    wflow_riverwidth = configget(self.config,"model","wflow_riverwidth","/staticmaps/wflow_riverwidth.map")
    
  
    # 2: Input base maps ########################################################  
    subcatch=ordinal(readmap(self.Dir + wflow_subcatch)) # Determines the area of calculations (all cells > 0)
    subcatch = ifthen(subcatch > 0, subcatch)
    if self.sCatch > 0:
        subcatch = ifthen(subcatch == sCatch,subcatch)
    
    self.Altitude=readmap(self.Dir + wflow_dem) * scalar(defined(subcatch)) #: The digital elevation map (DEM)
    self.TopoLdd=readmap(self.Dir + wflow_ldd)        #: The local drinage definition map (ldd)
    self.TopoId=readmap(self.Dir + wflow_subcatch)        #: Map define the area over which the calculations are done (mask)
    self.River=cover(boolean(readmap(self.Dir + wflow_river)),0) #: river network map. Fro those cell that belong to a river a specific width is used in the kinematic wave caulations
    self.RiverLength=pcrut.readmapSave(self.Dir + wflow_riverlength,0.0)
    # Factor to multiply riverlength with (defaults to 1.0)    
    self.RiverLengthFac=pcrut.readmapSave(self.Dir + wflow_riverlength_fact,1.0)

    # read landuse and soilmap and make sure there are no missing points related to the
    # subcatchment map. Currently sets the lu and soil type  type to 1
    self.LandUse=readmap(self.Dir + wflow_landuse)#: Map with lan-use/cover classes
    self.LandUse=cover(self.LandUse,nominal(ordinal(subcatch) > 0))
    self.Soil=readmap(self.Dir + wflow_soil)#: Map with soil classes
    self.Soil=cover(self.Soil,nominal(ordinal(subcatch) > 0))
    self.OutputLoc=readmap(self.Dir + wflow_gauges)  #: Map with locations of output gauge(s)
    self.InflowLoc=nominal(pcrut.readmapSave(self.Dir + wflow_inflow,0.0))  #: Map with location of abstractions/inflows.
    self.SeepageLoc=pcrut.readmapSave(self.Dir + wflow_inflow,0.0)  #: Seapage from external model (if configured)
    RiverWidth=pcrut.readmapSave(self.Dir + wflow_riverwidth,0.0)
    
    
     # Temperature correction per cell to add
    self.TempCor=pcrut.readmapSave(self.Dir + configget(self.config,"model","TemperatureCorrectionMap","staticmapswflow_tempcor.map"),0.0)
 
                      
    if self.scalarInput:
        self.gaugesMap=readmap(self.Dir + wflow_mgauges) #: Map with locations of rainfall/evap/temp gauge(s). Only needed if the input to the model is not in maps
    self.OutputId=readmap(self.Dir + wflow_subcatch)       # location of subcatchment
  
    self.ZeroMap=0.0*scalar(defined(self.Altitude))                    #map with only zero's
  
    # 3: Input time series ###################################################
    self.P_mapstack=self.Dir + configget(self.config,"inputmapstacks","Precipitation","/inmaps/P") # timeseries for rainfall
    self.PET_mapstack=self.Dir + configget(self.config,"inputmapstacks","EvapoTranspiration","/inmaps/PET") # timeseries for rainfall"/inmaps/PET"          # potential evapotranspiration
    self.TEMP_mapstack=self.Dir + configget(self.config,"inputmapstacks","Temperature","/inmaps/TEMP") # timeseries for rainfall "/inmaps/TEMP"          # global radiation
    self.Inflow_mapstack=self.Dir + configget(self.config,"inputmapstacks","Inflow","/inmaps/IF") # timeseries for rainfall "/inmaps/IF" # in/outflow locations (abstractions)
    self.Seepage_mapstack=self.Dir + configget(self.config,"inputmapstacks","Seepage","/inmaps/SE") # timeseries for rainfall "/inmaps/SE" # in/outflow locations (abstractions)
    # For in memory override:
    self.P = self.ZeroMap
    self.PET = self.ZeroMap
    self.TEMP = self.ZeroMap
    # Set static initial values here #########################################
   
    
    self.Latitude  =  ycoordinate(boolean(self.Altitude))
    self.Longitude =  xcoordinate(boolean(self.Altitude))
  
    self.logger.info("Linking parameters to landuse, catchment and soil...")

    self.Beta = scalar(0.6) # For sheetflow   
    #self.M=lookupscalar(self.Dir + "/" + modelEnv['intbl'] + "/M.tbl" ,self.LandUse,subcatch,self.Soil) # Decay parameter in Topog_sbm    
    self.N=lookupscalar(self.Dir + "/" + self.intbl + "/N.tbl",self.LandUse,subcatch,self.Soil)  # Manning overland flow
    """ *Parameter:* Manning's N for all non-river cells """
    self.NRiver=lookupscalar(self.Dir + "/" + self.intbl + "/N_River.tbl",self.LandUse,subcatch,self.Soil)  # Manning river   
    """ Manning's N for all cells that are marked as a river """
    
    
    #HBV Soil params
    self.FC=self.readtblDefault(self.Dir + "/" + self.intbl + "/FC.tbl",self.LandUse,subcatch,self.Soil,260.0)  
    self.BetaSeepage= self.readtblDefault(self.Dir + "/" + self.intbl + "/BetaSeepage.tbl",self.LandUse,subcatch,self.Soil,1.8)  # exponent in soil runoff generation equation
    self.LP= self.readtblDefault(self.Dir + "/" + self.intbl + "/LP.tbl",self.LandUse,subcatch,self.Soil, 0.53000)  # fraction of Fieldcapacity below which actual evaporation=potential evaporation (LP) 
    self.K4= self.readtblDefault(self.Dir + "/" + self.intbl + "/K4.tbl",self.LandUse,subcatch,self.Soil, 0.02307) # Recession constant baseflow   #K4=0.07; BASEFLOW:LINEARRESERVOIR 
    if self.SetKquickFlow:
        self.KQuickFlow= self.readtblDefault(self.Dir + "/" + self.intbl + "/KQuickFlow.tbl",self.LandUse,subcatch,self.Soil, 0.09880) # recession rate at flow HQ     #KHQ=0.2; OUTFLOWUPPERZONE_NONLINEARRESERVOIR
        self.SUZ= self.readtblDefault(self.Dir + "/" + self.intbl + "/SUZ.tbl",self.LandUse,subcatch,self.Soil, 100.0) # Level over wich K0 is used
        self.K0= self.readtblDefault(self.Dir + "/" + self.intbl + "/K0.tbl",self.LandUse,subcatch,self.Soil, 0.3) # K0
    else:
        self.KHQ= self.readtblDefault(self.Dir + "/" + self.intbl + "/KHQ.tbl",self.LandUse,subcatch,self.Soil, 0.09880) # recession rate at flow HQ     #KHQ=0.2; OUTFLOWUPPERZONE_NONLINEARRESERVOIR
        self.HQ= self.readtblDefault(self.Dir + "/" + self.intbl + "/HQ.tbl",self.LandUse,subcatch,self.Soil, 3.27000) # high flow rate HQ for which recession rate of upper reservoir is known   #HQ=3.76;    
        self.AlphaNL= self.readtblDefault(self.Dir + "/" + self.intbl + "/AlphaNL.tbl",self.LandUse,subcatch,self.Soil, 1.1) # measure of non-linearity of upper reservoir  #Alpha=1.6;
        
    self.PERC= self.readtblDefault(self.Dir + "/" + self.intbl + "/PERC.tbl",self.LandUse,subcatch,self.Soil, 0.4000) # percolation from Upper to Lowerzone (mm/day)
    self.CFR=self.readtblDefault(self.Dir + "/" + self.intbl + "/CFR.tbl",self.LandUse,subcatch,self.Soil, 0.05000)        # refreezing efficiency constant in refreezing of freewater in snow 
    #self.FoCfmax=self.readtblDefault(self.Dir + "/" + modelEnv['intbl'] + "/FoCfmax.tbl",self.LandUse,subcatch,self.Soil, 0.6000)  # correcton factor for snow melt/refreezing in forested and non-forested areas
    self.Pcorr=self.readtblDefault(self.Dir + "/" + self.intbl + "/Pcorr.tbl",self.LandUse,subcatch,self.Soil, 1.0)      # correction factor for precipitation
    self.RFCF=self.readtblDefault(self.Dir + "/" + self.intbl + "/RFCF.tbl",self.LandUse,subcatch,self.Soil,1.0)      # correction factor for rainfall
    self.SFCF=self.readtblDefault(self.Dir + "/" + self.intbl + "/SFCF.tbl",self.LandUse,subcatch,self.Soil, 1.0)     # correction factor for snowfall
    self.Cflux= self.readtblDefault(self.Dir + "/" + self.intbl + "/Cflux.tbl",self.LandUse,subcatch,self.Soil, 2.0)        # maximum capillary rise from runoff response routine to soil moisture routine    
    self.ICF=  self.readtblDefault(self.Dir + "/" + self.intbl + "/ICF.tbl",self.LandUse,subcatch,self.Soil, 2.0)       # maximum interception storage (in forested AND non-forested areas)
    self.CEVPF= self.readtblDefault(self.Dir + "/" + self.intbl + "/CEVPF.tbl",self.LandUse,subcatch,self.Soil, 1.0)   # correction factor for potential evaporation (1.15 in in forested areas )
    self.EPF= self.readtblDefault(self.Dir + "/" + self.intbl + "/EPF.tbl",self.LandUse,subcatch,self.Soil, 0.0)    # exponent of correction factor for evaporation on days with precipitation
    self.ECORR= self.readtblDefault(self.Dir + "/" + self.intbl + "/ECORR.tbl",self.LandUse,subcatch,self.Soil, 1.0)    # evap correction
    # Soil Moisture  parameters 
    self.ECALT= self.ZeroMap+0.00000     # evaporation lapse per 100m  
    #self.Ecorr=self.ZeroMap+1            # correction factor for evaporation
   
   
    # HBV Snow parameters    
    # critical temperature for snowmelt and refreezing:  TTI= 1.000
    self.TTI=self.readtblDefault(self.Dir + "/" + self.intbl + "/TTI.tbl" ,self.LandUse,subcatch,self.Soil,1.0)        
    # TT = -1.41934 # defines interval in which precipitation falls as rainfall and snowfall
    self.TT=self.readtblDefault(self.Dir + "/" + self.intbl + "/TT.tbl" ,self.LandUse,subcatch,self.Soil,-1.41934)
    #Cfmax = 3.75653 # meltconstant in temperature-index
    self.Cfmax=self.readtblDefault(self.Dir + "/" + self.intbl + "/Cfmax.tbl" ,self.LandUse,subcatch,self.Soil,3.75653)
    # WHC= 0.10000        # fraction of Snowvolume that can store water
    self.WHC=self.readtblDefault(self.Dir + "/" + self.intbl + "/WHC.tbl" ,self.LandUse,subcatch,self.Soil,0.1)
    
    # Determine real slope and cell length
    self.xl,self.yl,self.reallength = pcrut.detRealCellLength(self.ZeroMap,sizeinmetres)
    self.Slope= slope(self.Altitude)
    self.Slope=ifthen(boolean(self.TopoId),max(0.001,self.Slope*celllength()/self.reallength))
    Terrain_angle=scalar(atan(self.Slope))
   
    # Multiply parameters with a factor (for calibration etc) -P option in command line
    for k, v in multpars.iteritems():
        if self.sCatch > 0:
            estr = k + "= ifthenelse(self.TopoId == self.sCatch," +  k + "*" + str(v)+ "," + k + ")"
        else:
            estr = k + "=" + k + "*" + str(v)
        self.logger.info("Parameter multiplication: " +  estr)
        exec estr
    
    self.N=ifthenelse(self.River, self.NRiver, self.N)
    

    # Determine river width from DEM, upstream area and yearly average discharge
    # Scale yearly average Q at outlet with upstream are to get Q over whole catchment
    # Alf ranges from 5 to > 60. 5 for hardrock. large values for sediments
    # "Noah J. Finnegan et al 2005 Controls on the channel width of rivers: 
    # Implications for modeling fluvial incision of bedrock"

    upstr = catchmenttotal(1, self.TopoLdd)
    Qscale = upstr/mapmaximum(upstr) * Qmax
    W = (alf * (alf + 2.0)**(0.6666666667))**(0.375) * Qscale**(0.375) * (max(0.0001,windowaverage(self.Slope,celllength() * 4.0)))**(-0.1875) * self.N **(0.375)
    # Use supplied riverwidth if possible, else calulate
    RiverWidth = ifthenelse(RiverWidth <=0.0, W, RiverWidth)
    
    self.SnowWater = self.ZeroMap
    

    # Which columns/gauges to use/ignore in kinematic wave updating   
    self.UpdateMap = self.ZeroMap 
    
    if self.updating: 
        _tmp =pcr2numpy(self.OutputLoc,0.0)
        gaugear= _tmp
        touse = numpy.zeros(gaugear.shape,dtype='int')  

        for thecol in updateCols:
            idx = (gaugear == thecol).nonzero()
            touse[idx] = thecol
            
        self.UpdateMap = numpy2pcr(Nominal,touse,0.0)
        # Calculate distance to updating points (upstream) annd use to scale the correction
        # ldddist returns zero for cell at the gauges so add 1.0 tp result
        self.DistToUpdPt = cover(min(ldddist(self.TopoLdd,boolean(cover(self.UpdateMap,0)),1) * self.reallength/celllength(),self.UpdMaxDist),self.UpdMaxDist)
        #self.DistToUpdPt = ldddist(self.TopoLdd,boolean(cover(self.OutputId,0.0)),1)
        #* self.reallength/celllength()


    # Initializing of variables
    self.logger.info("Initializing of model variables..")
    self.TopoLdd=lddmask(self.TopoLdd,boolean(self.TopoId))   
    catchmentcells=maptotal(scalar(self.TopoId))
 
    # Used to seperate output per LandUse/management classes
    #OutZones = self.LandUse
    #report(self.reallength,"rl.map")
    #report(catchmentcells,"kk.map")
    self.QMMConv = self.timestepsecs/(self.reallength * self.reallength * 0.001) #m3/s --> mm
    self.ToCubic = (self.reallength * self.reallength * 0.001) / self.timestepsecs # m3/s
    self.sumprecip=self.ZeroMap #: accumulated rainfall for water balance
    self.sumevap=self.ZeroMap   #: accumulated evaporation for water balance
    self.sumrunoff=self.ZeroMap                          #: accumulated runoff for water balance (weigthted for upstream area)
    self.sumlevel=self.ZeroMap                          #: accumulated level for water balance
    self.sumpotevap=self.ZeroMap                          #accumulated runoff for water balance
    self.sumtemp=self.ZeroMap                          #accumulated runoff for water balance
    self.ForecQ_qmec=self.ZeroMap  # Extra inflow to kinematic wave reservoir for forcing in m^/sec
    self.KinWaveVolume=self.ZeroMap
    self.OldKinWaveVolume=self.ZeroMap
    self.Qvolume=self.ZeroMap
    self.Q=self.ZeroMap
    self.suminflow=self.ZeroMap
    # cntd
    self.FieldCapacity=self.FC                               #: total water holding capacity of the soil
    self.Treshold=self.LP*self.FieldCapacity                      # Threshold soilwaterstorage above which AE=PE
    #CatSurface=maptotal(scalar(ifthen(scalar(self.TopoId)>scalar(0.0),scalar(1.0))))                   # catchment surface (in  km2) 

   
    self.Aspect=scalar(aspect(self.Altitude))# aspect [deg]
    self.Aspect  = ifthenelse(self.Aspect <= 0.0 , scalar(0.001),self.Aspect)
    # On Flat areas the Aspect function fails, fill in with average...
    self.Aspect = ifthenelse (defined(self.Aspect), self.Aspect, areaaverage(self.Aspect,self.TopoId))

    

    # Set DCL to riverlength if that is longer that the basic length calculated from grid  
    drainlength = detdrainlength(self.TopoLdd,self.xl,self.yl)
    
    self.DCL=max(drainlength,self.RiverLength) # m
    # Multiply with Factor (taken from upscaling operation, defaults to 1.0 if no map is supplied
    self.DCL = self.DCL * max(1.0,self.RiverLengthFac)
    
    # water depth (m) 
    # set width for kinematic wave to cell width for all cells
    self.Bw=detdrainwidth(self.TopoLdd,self.xl,self.yl)
    # However, in the main river we have real flow so set the width to the 
    # width of the river
    
    self.Bw=ifthenelse(self.River, RiverWidth, self.Bw)
    
    # term for Alpha                             
    self.AlpTerm=pow((self.N/(sqrt(self.Slope))),self.Beta)
    # power for Alpha
    self.AlpPow=(2.0/3.0)*self.Beta
    # initial approximation for Alpha
    
    # calculate catchmentsize
    self.upsize=catchmenttotal(self.xl * self.yl,self.TopoLdd)
    self.csize=areamaximum(self.upsize,self.TopoId)


    self.logger.info("End of initial section.")


  def default_summarymaps(self):
      """
      Returns a list of default summary-maps at the end of a run.
      This is model specific. You can also add them to the [summary]section of the ini file but stuff
      you think is crucial to the model should be listed here

       Example:

      """
      lst = ['self.Cfmax','self.csize','self.upsize','self.TTI','self.TT','self.WHC',
             'self.Slope','self.N','self.xl','self.yl','self.reallength','self.DCL','self.Bw',]

      return lst

  def resume(self):
    """ read initial state maps (they are output of a previous call to suspend()) """
    
    if self.reinit == 1:
        self.logger.info("Setting initial conditions to default (zero!)")
        self.FreeWater =  cover(0.0) #: Water on surface (state variable [mm])
        self.SoilMoisture =  self.FC #: Soil moisture (state variable [mm])
        self.UpperZoneStorage = 0.2 * self.FC #: Storage in Upper Zone (state variable [mm])
        self.LowerZoneStorage = 1.0/(3.0 * self.K4) #: Storage in Uppe Zone (state variable [mm])
        self.InterceptionStorage = cover(0.0) #: Interception Storage (state variable [mm])
        self.SurfaceRunoff = cover(0.0) #: Discharge in kinimatic wave (state variable [m^3/s])
        self.WaterLevel = cover(0.0) #: Water level in kinimatic wave (state variable [m])
        self.DrySnow=cover(0.0) #: Snow amount (state variable [mm])
    else:
        self.wf_resume(self.Dir + "/instate/")

    P=self.Bw+(2.0*self.WaterLevel)
    self.Alpha=self.AlpTerm*pow(P,self.AlpPow)

    self.OldSurfaceRunoff = self.SurfaceRunoff
    
    self.SurfaceRunoffMM=self.SurfaceRunoff * self.QMMConv
        # Determine initial kinematic wave volume
    self.KinWaveVolume = self.WaterLevel * self.Bw * self.DCL
    self.OldKinWaveVolume = self.KinWaveVolume
    
    if not self.SetKquickFlow:
        self.KQuickFlow=(self.KHQ**(1.0+self.AlphaNL))*(self.HQ**-self.AlphaNL)   # recession rate of the upper reservoir, KHQ*UHQ=HQ=kquickflow*(UHQ**alpha)
 
    
    
  def dynamic(self):
    
    """    
    Below a list of variables that can be save to disk as maps or as 
    timeseries (see ini file for syntax):
        
    *Dynamic variables*
    
    :var self.SurfaceRunoff: Surface runoff in the kinematic wave [m^3/s]
    :var self.WaterLevel: Water level in the kinematic wave [m] (above the bottom)
    :var self.InterceptionStorage: actual interception storage [mm]
    :var self.Snow: Snow depth [mm]
    :var self.SnowWater: water content of the snow [mm]
    :var self.LowerZoneStorage: water content of the lower zone [mm]
    :var self.UpperZoneStorage: water content of the Upper zone [mm]
    :var self.BaseFlow: Specific runoff (baseflow part) per cell [mm]
    :var self.Percolation: actual percolation to the lower zone [mm]
    :var self.SoilMoisture: actual soil moisture [mm]
    :var self.QuickFlow: specific runoff (quickflow part) [mm]
    :var self.RealQuickFlow: specific runoff (quickflow), If K upper zone is precalculated [mm]
    :var self.CapFlux: capilary rise [mm]
    :var self.SurfaceRunoffMM: SurfaceRunoff in mm
    :var self.KinWaveVolume: Volume in the kinematic wave reservoir
        
    
    *Static variables*
        
    :var self.Altitude: The altitude of each cell [m]
    :var self.Bw: Width of the river [m]
    :var self.River: booolean map indicating the presence of a river [-]
    :var self.DLC: length of the river within a cell [m]
    :var self.ToCubic: Mutiplier to convert mm to m^3/s for fluxes
    """

    self.logger.debug("Step: "+str(int(self.thestep + self._d_firstTimeStep))+"/"+str(int(self._d_nrTimeSteps)))
    self.thestep = self.thestep + 1

    
    if self.scalarInput:
            self.Precipitation = timeinputscalar(self.precipTss,self.gaugesMap) * self.Pcorr
            self.Inflow=cover(timeinputscalar(self.inflowTss,self.InflowLoc),0.0)
            #self.Seepage = cover(timeinputscalar(self.SeepageTss,self.SeepageLoc),0)
            self.Precipitation = pcrut.interpolategauges(self.Precipitation,self.interpolMethod)
            self.PotEvaporation=timeinputscalar(self.evapTss,self.gaugesMap)
            self.PotEvaporation = pcrut.interpolategauges(self.PotEvaporation,self.interpolMethod)
            #self.report(self.PotEvaporation,'p')
            self.Temperature=timeinputscalar(self.tempTss,self.gaugesMap)
            self.Temperature = pcrut.interpolategauges(self.Temperature,self.interpolMethod)
            self.Temperature = self.Temperature + self.TempCor
    else:
            self.Precipitation=cover(self.wf_readmap(self.P_mapstack,0.0),0.0) * self.Pcorr
            self.PotEvaporation=cover(self.wf_readmap(self.PET_mapstack,0.0),0.0)
            self.Inflow=cover(self.wf_readmap(self.Inflow_mapstack,0.0),0.0)
            # These ar ALWAYS 0 at present!!!
            #self.Inflow=pcrut.readmapSave(self.Inflow_mapstack,0.0)
            if self.ExternalQbase:
                self.Seepage = cover(self.wf_readmap(self.Seepage_mapstack,0.0),0.0)
            else:
                self.Seepage=cover(0.0)
            self.Temperature=cover(self.wf_readmap(self.TEMP_mapstack,10.0),10.0)
            self.Temperature = self.Temperature + self.TempCor

    # Multiply input parameters with a factor (for calibration etc) -p option in command line
    for k, v in multdynapars.iteritems():
        estr = k + "=" + k + "*" + str(v)
        self.logger.debug("Dynamic Parameter multiplication: " +  estr)
        exec estr    
    
    RainFrac=ifthenelse(1.0*self.TTI == 0.0,ifthenelse(self.Temperature <= self.TT,scalar(0.0),scalar(1.0)),min((self.Temperature-(self.TT-self.TTI/2.0))/self.TTI,scalar(1.0)))  
    RainFrac=max(RainFrac,scalar(0.0))               #fraction of precipitation which falls as rain 
    SnowFrac=1.0-RainFrac                    #fraction of self.Precipitation which falls as snow
    self.Precipitation=self.SFCF*SnowFrac*self.Precipitation+self.RFCF*RainFrac*self.Precipitation # different correction for rainfall and snowfall
    
    self.PotEvaporation=exp(-self.EPF*self.Precipitation)*self.ECORR * self.PotEvaporation  # correction for potential evaporation on wet days
    self.PotEvaporation=self.CEVPF*self.PotEvaporation  # Correct per landuse

    SnowFall=SnowFrac*self.Precipitation  #: snowfall depth 
    RainFall=RainFrac*self.Precipitation  #: rainfall depth
    PotSnowMelt=ifthenelse(self.Temperature > self.TT,self.Cfmax*(self.Temperature-self.TT),scalar(0.0)) #Potential snow melt, based on temperature
    PotRefreezing=ifthenelse(self.Temperature < self.TT, self.Cfmax*self.CFR*(self.TT-self.Temperature),0.0)    #Potential refreezing, based on temperature
    
    
    #PotSnowMelt=self.FoCfmax*PotSnowMelt     	#correction for forest zones 0.6)
    #PotRefreezing=self.FoCfmax*PotRefreezing
    Refreezing=ifthenelse(self.Temperature < self.TT,min(PotRefreezing,self.FreeWater),0.0)   	#actual refreezing    
    self.SnowMelt=min(PotSnowMelt,self.DrySnow)          #actual snow melt
    self.DrySnow=self.DrySnow+SnowFall+Refreezing-self.SnowMelt     #dry snow content 
    self.FreeWater=self.FreeWater-Refreezing               #free water content in snow
    MaxFreeWater=self.DrySnow*self.WHC                     
    self.FreeWater=self.FreeWater+self.SnowMelt+RainFall
    InSoil = max(self.FreeWater-MaxFreeWater,0.0)   #abundant water in snow pack which goes into soil
    self.FreeWater=self.FreeWater-InSoil 

    
    #first part of precipitation is intercepted
    Interception=min(InSoil,self.ICF-self.InterceptionStorage)#: Interception in mm/timestep
    self.InterceptionStorage=self.InterceptionStorage+Interception #: Current interception storage
    NetInSoil=InSoil-Interception   
                  
    self.SoilMoisture=self.SoilMoisture+NetInSoil   
    DirectRunoff=max(self.SoilMoisture-self.FieldCapacity,0.0)    	#if soil is filled to capacity: abundant water runs of directly
    self.SoilMoisture=self.SoilMoisture-DirectRunoff            
    NetInSoil=NetInSoil-DirectRunoff                  		#net water which infiltrates into soil

    IntEvap=min(self.InterceptionStorage,self.PotEvaporation) 	 #: Evaporation from interception storage
    self.InterceptionStorage=self.InterceptionStorage-IntEvap
    
    # I nthe origal HBV code
    RestEvap = max(0.0,self.PotEvaporation-IntEvap)   
        
    SoilEvap=ifthenelse(self.SoilMoisture > self.Treshold,RestEvap,min(RestEvap,self.PotEvaporation*(self.SoilMoisture/self.Treshold))) #: soil evapotranspiration
    self.SoilMoisture=self.SoilMoisture-SoilEvap           #evaporation from soil moisture storage
    
   
    ActEvap=IntEvap+SoilEvap           #: Sum of evaporation components (IntEvap+SoilEvap)
    HBVSeepage=((self.SoilMoisture/self.FieldCapacity)**self.BetaSeepage)*NetInSoil		#runoff water from soil
    self.SoilMoisture=self.SoilMoisture-HBVSeepage        

    Backtosoil=min(self.FieldCapacity-self.SoilMoisture,DirectRunoff) 		#correction for extremely wet periods: soil is filled to capacity
    DirectRunoff=DirectRunoff-Backtosoil
    self.SoilMoisture=self.SoilMoisture+Backtosoil     
    InUpperZone=DirectRunoff+HBVSeepage                         		# total water available for runoff
	
    # Steps is always 1 at the moment
    # calculations for Upper zone
    self.UpperZoneStorage=self.UpperZoneStorage+InUpperZone                 		#incoming water from soil 
    self.Percolation=min(self.PERC,self.UpperZoneStorage)                        		#Percolation
    self.UpperZoneStorage=self.UpperZoneStorage-self.Percolation                  
    self.CapFlux=self.Cflux*(((self.FieldCapacity-self.SoilMoisture)/self.FieldCapacity))   #: Capillary flux flowing back to soil
    self.CapFlux=min(self.UpperZoneStorage,self.CapFlux)
    self.CapFlux=min(self.FieldCapacity-self.SoilMoisture,self.CapFlux)
    self.UpperZoneStorage=self.UpperZoneStorage-self.CapFlux
    self.SoilMoisture=self.SoilMoisture+self.CapFlux

    if not self.SetKquickFlow:
        self.QuickFlow = max(0,self.KQuickFlow*(self.UpperZoneStorage**(1.0+self.AlphaNL)))
        self.RealQuickFlow = self.ZeroMap
    else:
        self.QuickFlow = self.KQuickFlow*self.UpperZoneStorage
        self.RealQuickFlow = max(0,self.K0*(self.UpperZoneStorage - self.SUZ))
    
    """Quickflow volume in mm/timestep"""
    self.UpperZoneStorage=self.UpperZoneStorage-self.QuickFlow-self.RealQuickFlow
 			      
    # calculations for Lower zone
    self.LowerZoneStorage=self.LowerZoneStorage+self.Percolation			
    self.BaseFlow=self.K4*self.LowerZoneStorage #: Baseflow in mm/timestep
    self.LowerZoneStorage=self.LowerZoneStorage-self.BaseFlow
    # Direct runoff generation
    if self.ExternalQbase:
        DirectRunoffStorage=self.QuickFlow+self.Seepage+self.RealQuickFlow
    else:
        DirectRunoffStorage=self.QuickFlow+self.BaseFlow+self.RealQuickFlow

    self.InwaterMM=max(0.0,DirectRunoffStorage)
    self.Inwater=self.InwaterMM * self.ToCubic
    self.QuickFlowCubic = (self.QuickFlow + self.RealQuickFlow) * self.ToCubic
    self.BaseFlowCubic = self.BaseFlow * self.ToCubic
    self.Inwater=self.Inwater + self.Inflow # Add abstractions/inflows in m^3/sec
    
    ##########################################################################
    # Runoff calculation via Kinematic wave ##################################
    ##########################################################################
    # per distance along stream
    q=self.Inwater/self.DCL + self.ForecQ_qmec/self.DCL
    self.OldSurfaceRunoff=self.SurfaceRunoff
    
    self.SurfaceRunoff = kinematic(self.TopoLdd, self.SurfaceRunoff,q,self.Alpha, self.Beta,self.Tslice,self.timestepsecs,self.DCL) # m3/s
    self.SurfaceRunoffMM=self.SurfaceRunoff*self.QMMConv # SurfaceRunoffMM (mm) from SurfaceRunoff (m3/s)
      
    
    self.updateRunOff()
    InflowKinWaveCell=upstream(self.TopoLdd,self.SurfaceRunoff)
    self.MassBalKinWave = (self.KinWaveVolume - self.OldKinWaveVolume)/self.timestepsecs  + InflowKinWaveCell + self.Inwater - self.SurfaceRunoff
    Runoff=self.SurfaceRunoff

    # Updating
    # --------
    # Assume a tss file with as many columns as outpulocs. Start updating for each non-missing value and start with the
    # first column (nr 1). Assumes that outputloc and columns match!

    if self.updating:
        QM = timeinputscalar(self.updateFile, self.UpdateMap) * self.QMMConv
        
        # Now update the state. Just add to the Ustore
        # self.UStoreDepth =  result
        # No determine multiplication ratio for each gauge influence area.
        # For missing gauges 1.0 is assumed (no change).
        # UpDiff = areamaximum(QM,  self.UpdateMap) - areamaximum(self.SurfaceRunoffMM, self.UpdateMap)
        UpRatio = areamaximum(QM,  self.UpdateMap)/areamaximum(self.SurfaceRunoffMM, self.UpdateMap)
    
        UpRatio = cover(areaaverage(UpRatio,self.TopoId),1.0)
        # Now split between Soil and Kyn  wave
        self.UpRatioKyn = min(self.MaxUpdMult,max(self.MinUpdMult,(UpRatio - 1.0) * self.UpFrac + 1.0))
        UpRatioSoil = min(self.MaxUpdMult,max(self.MinUpdMult,(UpRatio - 1.0) * (1.0 - self.UpFrac) + 1.0))

        # update/nudge self.UStoreDepth for the whole upstream area, 
        # not sure how much this helps or worsens things
        UpdSoil = True
        if UpdSoil:      
            toadd = (self.UpperZoneStorage * UpRatioSoil) - self.UpperZoneStorage
            self.UpperZoneStorage = self.UpperZoneStorage + toadd
        
        # Update the kinematic wave reservoir up to a maximum upstream distance
        # TODO:  add (much smaller) downstream updating also?
        MM = (1.0 - self.UpRatioKyn)/self.UpdMaxDist
        self.UpRatioKyn = MM * self.DistToUpdPt + self.UpRatioKyn
        
        self.SurfaceRunoff = self.SurfaceRunoff *  self.UpRatioKyn
        self.SurfaceRunoffMM=self.SurfaceRunoff*self.QMMConv # SurfaceRunoffMM (mm) from SurfaceRunoff (m3/s)
        self.updateRunOff()
 
        Runoff=self.SurfaceRunoff
        

    
    self.sumprecip=self.sumprecip  +  self.Precipitation                     #accumulated rainfall for water balance
    self.sumevap=self.sumevap + ActEvap                           #accumulated evaporation for water balance
    self.sumpotevap=self.sumpotevap + self.PotEvaporation 
    self.sumtemp=self.sumtemp + self.Temperature
    self.sumrunoff=self.sumrunoff  + self.SurfaceRunoffMM                        #accumulated runoff for water balance
    self.sumlevel=self.sumlevel  + self.WaterLevel
    self.suminflow=self.suminflow  + self.Inflow
    
    



# The main function is used to run the program from the command line

def main(argv=None):  
    """
    Perform command line execution of the model.
    """      
    global multpars
    global updateCols
    caseName = "default_hbv"
    runId = "run_default"
    configfile="wflow_hbv.ini"
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
        opts, args = getopt.getopt(argv, 'c:QXS:F:hC:Ii:T:NR:u:s:P:p:Xx:U:fl:')
    except getopt.error, msg:
        pcrut.usage(msg)
    
    for o, a in opts:
        if o == '-F': 
            runinfoFile = a
            fewsrun = True
        if o == '-P': 
            exec ("multpars =" + a,globals(), globals())
        if o == '-p': 
            exec "multdynapars =" + a
            exec ("multdynapars =" + a,globals(), globals())
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
        if (ts):
            _lastTimeStep =  ts# * 86400/timestepsecs
            _firstTimeStep = 1 
        else:
            print "Failed to get timesteps from runinfo file: " + runinfoFile
            exit(2)
       
    if _lastTimeStep < _firstTimeStep:
        print "The starttimestep (" + str(_firstTimeStep) +") is smaller than the last timestep (" + str(_lastTimeStep) + ")"
        usage()
 
    myModel = WflowModel(wflow_cloneMap, caseName,runId,configfile)
    dynModelFw = wf_DynamicFramework(myModel, _lastTimeStep,firstTimestep=_firstTimeStep)
    dynModelFw.createRunId(NoOverWrite=NoOverWrite,logfname=LogFileName,level=loglevel)    
    
    for o, a in opts:
        if o == '-X': configset(myModel.config,'model','OverWriteInit','1',overwrite=True) 
        if o == '-I': configset(myModel.config,'model','reinit','1',overwrite=True) 
        if o == '-i': configset(myModel.config,'model','intbl',a,overwrite=True)
        if o == '-s': configset(myModel.config,'model','timestepsecs',a,overwrite=True)
        if o == '-x': configset(myModel.config,'model','sCatch',a,overwrite=True)
        if o == '-c': configset(myModel.config,'model','configfile', a,overwrite=True)
        if o == '-M': configset(myModel.config,'model','MassWasting',"1",overwrite=True)
        if o == '-N': configset(myModel.config,'model','nolateral','1',overwrite=True) 
        if o == '-Q': configset(myModel.config,'model','ExternalQbase','1',overwrite=True)
        if o == '-U': 
            configset(myModel.config,'model','updateFile',a,overwrite=True)
            configset(myModel.config,'model','updating',"1",overwrite=True)
        if o == '-u': 
            exec "zz =" +  a
            updateCols = zz


    dynModelFw._runInitial()
    dynModelFw._runResume()
    dynModelFw._runDynamic(_firstTimeStep,_lastTimeStep)
    dynModelFw._runSuspend()
    
    
    fp = open(caseName + "/" + runId + "/runinfo/configofrun.ini",'wb')
    myModel.config.write(fp)
    dynModelFw._wf_shutdown()
    
    
    
    os.chdir("../../")


if __name__ == "__main__":
    main()
