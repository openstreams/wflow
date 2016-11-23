#!/usr/bin/python

"""
Definition of the wflow_topoflex model.
---------------------------------------

This simple model calculates soil temperature using
air temperature as a forcing.

Usage:
wflow_topoflex  -C case -R Runid -c inifile

    -C: set the name  of the case (directory) to run
    
    -R: set the name runId within the current case
    
    -c name of the config file (in the case directory)
    

"""

import wflow.reservoir_Si as reservoir_Si
import wflow.reservoir_Sa as reservoir_Sa
import wflow.reservoir_Sw as reservoir_Sw
import wflow.reservoir_Su as reservoir_Su
import wflow.reservoir_Sf as reservoir_Sf
import wflow.reservoir_Ss as reservoir_Ss
import wflow.JarvisCoefficients as JarvisCoefficients

import numpy
import os
import os.path
import shutil, glob
import getopt
import time
import pdb

from wflow.wf_DynamicFramework import *
from wflow.wflow_adapt import *
from copy import deepcopy as copylist

# TODO: see below
"""
Verwijderen IRURFR_L statements?
Documentatie updaten!
Wegschrijven per class, afkortingen van classes gebruiken (zie outputtss_0) Jaap!
Routing functies in apart file onderbrengen, aanroepen, configureerbaar maken welke gebruiken Hessel!
logging toevoegen, ervoor zorgen dat het ook 1 per x aantal stappen weggeschreven kan worden
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

    def __init__(self, cloneMap, Dir, RunDir, configfile):
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

        # Static model parameters
        modelparameters.append(
            self.ParamType(name="Altitude", stack="staticmaps/wflow_dem.map", type="staticmap", default=0.0,
                           verbose=False, lookupmaps=[]))

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
        states = ['Si', 'Su', 'Sf', 'Ss', 'Sw', 'Sa', 'Sfa', 'Qstate']

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

        return self.currentTimeStep() * int(configget(self.config, 'model', 'timestepsecs', '86400'))


    def readtblDefault2(self,pathtotbl,landuse,subcatch,soil, default):
        """
        First check if a prepared  maps of the same name is present
        in the staticmaps directory. next try to
        read a tbl file to match a landuse, catchment and soil map. Returns 
        the default value if the tbl file is not found.
    
        Input:
            -  pathtotbl: full path to table file
            -  landuse: landuse map
            -  subcatch: subcatchment map
            -  soil: soil map
            -  default: default value
    
        Output: 
            - map constructed from tbl file or map with default value
        """
    
        mapname = os.path.dirname(pathtotbl) + "/../staticmaps/" + os.path.splitext(os.path.basename(pathtotbl))[0]+".map"
        if os.path.exists(mapname):
            self.logger.info("reading map parameter file: " + mapname)
            rest = cover(readmap(mapname),default)
        else:
            if os.path.isfile(pathtotbl):
                rest=cover(lookupscalar(pathtotbl,landuse,subcatch,soil), default)
                self.logger.info("Creating map from table: " + pathtotbl)
            else:
                self.logger.warn("tbl file not found (" + pathtotbl + ") returning default value: " + str(default))
                rest = scalar(default)
        
        return rest 
        
        
    def suspend(self):
        """
      *Required*
      
      Suspends the model to disk. All variables needed to restart the model
      are saved to disk as pcraster maps. Use resume() to re-read them
      
      This function is required. 
      
    """

        self.logger.info("Saving initial conditions...")
        
        if self.fewsrun:
            self.logger.info("Saving initial conditions for FEWS...")
#            self.wf_suspend(os.path.join(self.Dir, "outstate"))
            [report(self.Si[i], self.Dir + "/outstate/Si" + self.NamesClasses[i] + ".map") for i in self.Classes if self.selectSi[i]]
            [report(self.Su[i], self.Dir + "/outstate/Su" + self.NamesClasses[i] + ".map") for i in self.Classes if self.selectSu[i]]
            [report(self.Sa[i], self.Dir + "/outstate/Sa" + self.NamesClasses[i] + ".map") for i in self.Classes if self.selectSa[i]]
            [report(self.Sf[i], self.Dir + "/outstate/Sf" + self.NamesClasses[i] + ".map") for i in self.Classes if self.selectSf[i]]
            [report(self.Sfa[i], self.Dir + "/outstate/Sfa" + self.NamesClasses[i] + ".map") for i in self.Classes if self.selectSfa[i]]
            [report(self.Sw[i], self.Dir + "/outstate/Sw" + self.NamesClasses[i] + ".map") for i in self.Classes if self.selectSw[i]]
            report(self.Ss, self.Dir + "/outstate/Ss.map")
            report(self.Qstate, self.Dir + "/outstate/Qstate.map")
        
        #: It is advised to use the wf_suspend() function
        #: here which will suspend the variables that are given by stateVariables
        #: function.
        [report(self.Si[i], self.SaveDir + "/outstate/Si" + self.NamesClasses[i] + ".map") for i in self.Classes if self.selectSi[i]]
        [report(self.Su[i], self.SaveDir + "/outstate/Su" + self.NamesClasses[i] + ".map") for i in self.Classes if self.selectSu[i]]
        [report(self.Sa[i], self.SaveDir + "/outstate/Sa" + self.NamesClasses[i] + ".map") for i in self.Classes if self.selectSa[i]]
        [report(self.Sf[i], self.SaveDir + "/outstate/Sf" + self.NamesClasses[i] + ".map") for i in self.Classes if self.selectSf[i]]
        [report(self.Sfa[i], self.SaveDir + "/outstate/Sfa" + self.NamesClasses[i] + ".map") for i in self.Classes if self.selectSfa[i]]
        [report(self.Sw[i], self.SaveDir + "/outstate/Sw" + self.NamesClasses[i] + ".map") for i in self.Classes if self.selectSw[i]]
        report(self.Ss, self.SaveDir + "/outstate/Ss.map")
        report(self.Qstate, self.SaveDir + "/outstate/Qstate.map")

        [report(self.percent[i], self.SaveDir + "/outmaps/percent" + self.NamesClasses[i] + ".map") for i in
         self.Classes]
        report(self.percentArea, self.SaveDir + "/outmaps/percentArea.map")
        report(self.surfaceArea, self.SaveDir + "/outmaps/surfaceArea.map")

        report(self.sumprecip, self.SaveDir + "/outsum/sumprecip.map")
        report(self.sumevap, self.SaveDir + "/outsum/sumevap.map")
        report(self.sumpotevap, self.SaveDir + "/outsum/sumpotevap.map")
        report(self.sumtemp, self.SaveDir + "/outsum/sumtemp.map")
        report(self.sumrunoff, self.SaveDir + "/outsum/sumrunoff.map")
        report(self.sumwb, self.SaveDir + "/outsum/sumwb.map")

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

        self.thestep = scalar(0)
        #: files to be used in case of timesries (scalar) input to the model
        # files for forcing data
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


        self.logger.info(
            "running for " + str(self.nrTimeSteps()) + " timesteps")  # keeping track of number of timesteps

        self.fewsrun = int(configget(self.config,"model","fewsrun","0"))
        
        # Set and get defaults from ConfigFile here ###################################
        self.timestepsecs = int(configget(self.config,
                                          "model", "timestepsecs", "3600"))  # number of seconds in a timestep
        self.scalarInput = int(configget(self.config,
                                         "model", "ScalarInput", "1"))  # forcing data in maps (0) or timeseries (1)
        self.InputSeries = int(configget(self.config,
                                         "model", "InputSeries", "1"))  # forcing data in maps (0) or timeseries (1)
        self.reinit = int(configget(self.config,
                                    "run", "reinit", "0"))
        
        self.intbl = configget(self.config,
                                    "model","intbl","intbl")
                                    
        self.maxTransit = float(configget(self.config,
                                          "model", "maxTransitTime", "34"))  # maximum Transit time in cacthment
        self.distForcing = int(configget(self.config,
                                         "model", "DistForcing",
                                         "10"))  # number of different forcing inputs (eg. number of rainfall stations)
        self.maxGaugeId = int(configget(self.config,
                                        "model", "maxGaugeId", "10"))  # highest index of all used meteo stations
        self.IRURFR_L = int(configget(self.config,
                                      "model", "L_IRURFR",
                                      "0"))  # combination of reservoirs that are distributed (1: all these reservoirs are distributed)
        self.URFR_L = int(configget(self.config,
                                    "model", "L_URFR",
                                    "0"))  # combination of reservoirs that are distributed (1: all these reservoirs are distributed)
        self.FR_L = int(configget(self.config,
                                  "model", "L_FR",
                                  "0"))  # combination of reservoirs that are distributed (1: all these reservoirs are distributed)
        self.Ctime = int(configget(self.config,
                                   "model", "spinUp_time",
                                   "7775"))  # number of timesteps for which no data needs to be recorded
        self.NamesClasses = eval(str(configget(self.config,
                                               "model", "classes", "['W','H','P']")))  # classes used in model
        self.Classes = [x for x in range(len(self.NamesClasses))]  # numbering of classes

        # selection of reservoir conceputalisatie - codes are described in reservoir files
        self.selectSw = configget(self.config, "model",
                                  "selectSw", "0, 0, 0").replace(
            ' ', '').replace('[', '').replace(
            ']', '').replace(
            'None', '').split(',')
        self.selectSi = configget(self.config, "model",
                                  "selectSi", "0, 0, 0").replace(
            ' ', '').replace('[', '').replace(
            ']', '').replace(
            'None', '').split(',')
        self.selectSa = configget(self.config, "model",
                                  "selectSa", "0, 0, 0").replace(
            ' ', '').replace('[', '').replace(
            ']', '').replace(
            'None', '').split(',')        
        self.selectSu = configget(self.config, "model",
                                  "selectSu", "0, 0, 0").replace(
            ' ', '').replace('[', '').replace(
            ']', '').replace(
            'None', '').split(',')
        self.selectSf = configget(self.config, "model",
                                  "selectSf", "0, 0, 0").replace(
            ' ', '').replace('[', '').replace(
            ']', '').replace(
            'None', '').split(',')
        self.selectSfa = configget(self.config, "model",
                                  "selectSfa", "0, 0, 0").replace(
            ' ', '').replace('[', '').replace(
            ']', '').replace(
            'None', '').split(',')
        self.selectSs = configget(self.config, "model", "selectSs", "groundWaterCombined3")

            
        self.selectRout = configget(self.config, "model", "selectRout", " ")

        # static maps to use (normally default)
        wflow_subcatch = configget(self.config,
                                   "model", "wflow_subcatch", "staticmaps/wflow_subcatch.map")
        wflow_catchArea = configget(self.config,
                                   "model", "wflow_subcatch", "staticmaps/wflow_catchmentAreas.map")
        wflow_dem = configget(self.config,
                              "model", "wflow_dem", "staticmaps/wflow_dem.map")
        wflow_ldd = configget(self.config,
                              "model", "wflow_ldd", "staticmaps/wflow_ldd.map")
        wflow_landuse  = configget(self.config,
                              "model","wflow_landuse","staticmaps/wflow_landuse.map")
        wflow_soil  = configget(self.config,
                                "model","wflow_soil","staticmaps/wflow_soil.map")    
        wflow_gauges = configget(self.config,
                                 "model", "wflow_gauges", "staticmaps/wflow_gauges.map")
        wflow_mgauges = configget(self.config,
                                  "model", "wflow_mgauges", "staticmaps/wflow_mgauges.map")
        wflow_surfaceArea = configget(self.config,
                                      "model", "wflow_surfaceArea", "staticmaps/wflow_surfaceArea.map")
        wflow_transit = configget(self.config,
                                  "model", "wflow_transit", "staticmaps/wflow_transit.map")
        wflow_velocity = configget(self.config,
                                  "model", "wflow_velocity", "staticmaps/wflow_velocity.map")
        wflow_percent = [configget(self.config,
                                   "model", "wflow_percent_" + str(self.Classes[i]),
                                   "staticmaps/wflow_percent" + str(self.Classes[i]) + ".map") for i in self.Classes]
        self.rst_laiTss = [configget(self.config,
                                     "model", "rst_lai_" + str(self.Classes[i]),
                                     "staticmaps/rst_lai_" + str(self.Classes[i]) + ".map") for i in self.Classes]

        # 2: Input base maps ########################################################
        subcatch = ordinal(
            readmap(os.path.join(self.Dir, wflow_subcatch)))  # Determines the area of calculations (all cells > 0)
        subcatch = ifthen(subcatch > 0, subcatch)

        self.Altitude = readmap(os.path.join(self.Dir, wflow_dem)) * scalar(
            defined(subcatch))  #: The digital elevation map (DEM)
        self.TopoLdd = readmap(os.path.join(self.Dir, wflow_ldd))  #: The local drinage definition map (ldd)
        self.TopoId = readmap(
            os.path.join(self.Dir, wflow_subcatch))  #: Map define the area over which the calculations are done (mask)
        self.catchArea = scalar(ifthen(self.TopoId > 0, scalar(1)))
        self.LandUse=ordinal(self.wf_readmap(os.path.join(self.Dir , wflow_landuse),0.0,fail=True))#: Map with lan-use/cover classes
        self.LandUse=cover(self.LandUse,ordinal(ordinal(subcatch) > 0))
        self.Soil=ordinal(self.wf_readmap(os.path.join(self.Dir , wflow_soil),0.0,fail=True))#: Map with soil classes
        self.Soil=cover(self.Soil,ordinal(ordinal(subcatch) > 0)) 
        self.TopoId = ifthen(scalar(self.TopoId) > 0, self.TopoId)
        self.surfaceArea = scalar(readmap(os.path.join(self.Dir, wflow_surfaceArea)))  #: Map with surface area per cell
        self.totalArea = areatotal(self.surfaceArea, nominal(self.TopoId))
        self.percentArea = self.surfaceArea / self.totalArea
        self.Transit = scalar(readmap(os.path.join(self.Dir, wflow_transit)))  #: Map with surface area per cell
        self.velocity = scalar(readmap(os.path.join(self.Dir, wflow_velocity)))  #: Map with surface area per cell
        self.gaugesR = nominal(readmap(os.path.join(self.Dir, wflow_gauges)))
        self.percent = []
        for i in self.Classes:
            self.percent.append(readmap(os.path.join(self.Dir, wflow_percent[i])))

        # MODEL PARAMETERS - VALUES PER CLASS
        self.D = eval(str(configget(self.config, "model", "D", "[0]")))
        self.Tf = eval(str(configget(self.config, "model", "Tf", "[0]")))
        self.Tfa = eval(str(configget(self.config, "model", "Tfa", "[0]")))
        
        # MODEL PARAMETERS - BASED ON TABLES
        self.imax = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/imax" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,1.5) for i in self.Classes]            
        self.sumax = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/sumax" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,70) for i in self.Classes]        
        self.samax = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/samax" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,50) for i in self.Classes]
        self.samin = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/samin" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,0.1) for i in self.Classes]
        self.beta = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/beta" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,0.4) for i in self.Classes]
        self.betaA = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/betaA" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,0.2) for i in self.Classes]        
        self.Kf = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/Kf" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,0.005) for i in self.Classes]
        self.Kfa = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/Kfa" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,0.05) for i in self.Classes]
        self.perc = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/perc" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,0.0035) for i in self.Classes]
        self.cap = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/cap" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,0.028) for i in self.Classes]
        self.LP = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/LP" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,0.15) for i in self.Classes]
        self.Ks = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/Ks" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,0.0004) for i in self.Classes]
        self.Fmax = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/Fmax" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,1) for i in self.Classes]
        self.Fmin = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/Fmin" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,0) for i in self.Classes]
        self.decF = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/decF" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,0.5) for i in self.Classes]
        self.dayDeg = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/dayDeg" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,1) for i in self.Classes]
        self.FrDur0 = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/FrDur0" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,-5) for i in self.Classes]
        self.FrDur1 = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/FrDur1" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,0) for i in self.Classes]
        self.ratFT = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/ratFT" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,1) for i in self.Classes]
        self.Tt = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/Tt" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,1) for i in self.Classes]
        self.Tm = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/Tm" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,2) for i in self.Classes]
        self.Fm = [self.readtblDefault2(self.Dir + "/" + self.intbl + "/Fm" + self.NamesClasses[i] + ".tbl",self.LandUse,subcatch,self.Soil,0.2) for i in self.Classes]
        
        # Jarvis stressfunctions
        self.lamda = eval(str(configget(self.config, "model", "lamda", "[0]")))
        self.lamdaS = eval(str(configget(self.config, "model", "lamdaS", "[0]")))

        # initialise list for routing
        self.trackQ = [0 * scalar(self.catchArea)] * int(self.maxTransit)

        # initialise list for lag function
        self.convQu = [[0 * scalar(self.catchArea)] * self.Tf[i] for i in self.Classes]
        self.convQa = [[0 * scalar(self.catchArea)] * self.Tfa[i] for i in self.Classes]

        if self.scalarInput:
            self.gaugesMap = nominal(readmap(os.path.join(self.Dir,
                                                          wflow_mgauges)))  #: Map with locations of rainfall/evap/temp gauge(s). Only needed if the input to the model is not in maps
        self.OutputId = readmap(os.path.join(self.Dir, wflow_subcatch))  # location of subcatchment
        self.OutputIdRunoff = boolean(ifthenelse(self.gaugesR == 1, 1 * scalar(self.TopoId),
                                                 0 * scalar(self.TopoId)))  # location of subcatchment

        self.ZeroMap = 0.0 * scalar(subcatch)  # map with only zero's

        self.wf_multparameters()


        # For in memory override:
        self.P = self.ZeroMap
        self.PET = self.ZeroMap
        self.TEMP = self.ZeroMap

        self.logger.info("Linking parameters to landuse, catchment and soil...")


        # Initializing of variables

        self.logger.info("Initializing of model variables..")
        self.TopoLdd = lddmask(self.TopoLdd, boolean(self.TopoId))
        catchmentcells = maptotal(scalar(self.TopoId))

        self.sumprecip = self.ZeroMap  # accumulated rainfall for water balance
        self.sumevap = self.ZeroMap  # accumulated evaporation for water balance
        self.sumrunoff = self.ZeroMap  # accumulated runoff for water balance (weigthted for upstream area)
        self.sumpotevap = self.ZeroMap  # accumulated runoff for water balance
        self.sumtemp = self.ZeroMap  # accumulated runoff for water balance
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
            # self.logger.info("Setting initial conditions to default (zero!)")
            self.logger.info("Setting initial conditions to preset values in main script!!")
            self.Si = [self.ZeroMap] * len(self.Classes)
            self.Sw = [self.ZeroMap] * len(self.Classes)
            self.Su = [self.ZeroMap] * len(self.Classes)
            self.Sa = [self.ZeroMap] * len(self.Classes)
            self.Sf = [self.ZeroMap] * len(self.Classes)
            self.Sfa = [self.ZeroMap] * len(self.Classes)
            self.Ss = self.ZeroMap  # for combined gw reservoir
            self.Qstate = self.catchArea * 0  # for combined gw reservoir
            self.Qstate_t = self.catchArea * 0

            # set initial storage values
#            pdb.set_trace()
            self.Sa = [0.05 * self.samax[i] * scalar(self.catchArea) for i in self.Classes]
            self.Su = [self.sumax[i] * scalar(self.catchArea) for i in self.Classes]
            self.Ss = self.Ss + 30 * scalar(self.catchArea)  # for combined gw reservoir

        else:
#            self.wf_resume(self.Dir + "/instate/")
            
            self.Si = []
            for i in self.Classes:
                if self.selectSi[i]:
                    self.Si.append(readmap(os.path.join(self.Dir, 'instate', 'Si' + self.NamesClasses[i] + '.map')))
                else:
                    self.Si.append(self.ZeroMap)
            self.Sw = []
            for i in self.Classes:
                if self.selectSw[i]:
                    self.Sw.append(readmap(os.path.join(self.Dir, 'instate', 'Sw' + self.NamesClasses[i] + '.map')))
                else:
                    self.Sw.append(self.ZeroMap)
            self.Sa = []
            for i in self.Classes:
                if self.selectSa[i]:
                    self.Sa.append(readmap(os.path.join(self.Dir, 'instate', 'Sa' + self.NamesClasses[i] + '.map')))
                else:
                    self.Sa.append(self.ZeroMap)
            self.Su = []
            for i in self.Classes:
                if self.selectSu[i]:
                    self.Su.append(readmap(os.path.join(self.Dir, 'instate', 'Su' + self.NamesClasses[i] + '.map')))
                else:
                    self.Su.append(self.ZeroMap)
            self.Sf = []
            for i in self.Classes:
                if self.selectSf[i]:
                    self.Sf.append(readmap(os.path.join(self.Dir, 'instate', 'Sf' + self.NamesClasses[i] + '.map')))
                else:
                    self.Sf.append(self.ZeroMap)
            self.Sfa = []
            for i in self.Classes:
                if self.selectSfa[i]:                
                    self.Sfa.append(readmap(os.path.join(self.Dir, 'instate', 'Sfa' + self.NamesClasses[i] + '.map')))
                else:
                    self.Sfa.append(self.ZeroMap)
            self.Ss = readmap(os.path.join(self.Dir, 'instate', 'Ss.map'))
            self.Qstate = readmap(os.path.join(self.Dir, 'instate', 'Qstate.map'))

        self.wbSi_ = [self.ZeroMap] * len(self.Classes)
        self.wbSu_ = [self.ZeroMap] * len(self.Classes)
        self.wbSa_ = [self.ZeroMap] * len(self.Classes)
        self.wbSw_ = [self.ZeroMap] * len(self.Classes)
        self.wbSf_ = [self.ZeroMap] * len(self.Classes)
        self.wbSfa_ = [self.ZeroMap] * len(self.Classes)
        self.wbSfrout = self.ZeroMap
        self.wbSs = self.ZeroMap

        self.Ei_ = [self.ZeroMap] * len(self.Classes)
        self.Pe_ = [self.ZeroMap] * len(self.Classes)
        self.Si_ = [self.ZeroMap] * len(self.Classes)
        self.Eu_ = [self.ZeroMap] * len(self.Classes)
        self.Ea_ = [self.ZeroMap] * len(self.Classes)
        self.Ew_ = [self.ZeroMap] * len(self.Classes)
        self.Qu_ = [self.ZeroMap] * len(self.Classes)
        self.Qw_ = [self.ZeroMap] * len(self.Classes)
        self.Qa_ = [self.ZeroMap] * len(self.Classes)
        self.Cap_ = [self.ZeroMap] * len(self.Classes)
        self.Perc_ = [self.ZeroMap] * len(self.Classes)
        self.Fa_ = [self.ZeroMap] * len(self.Classes)
        self.Qf_ = [self.ZeroMap] * len(self.Classes)
        self.Qfa_ = [self.ZeroMap] * len(self.Classes)
        self.Qs_ = self.ZeroMap  # for combined gw reservoir
        self.Qflag_ = [self.ZeroMap] * len(self.Classes)
        self.Qfcub_ = [self.ZeroMap] * len(self.Classes)
        self.Ep_ = [self.ZeroMap] * len(self.Classes)
        self.EpD_ = [self.ZeroMap] * len(self.Classes)
        self.FrDur = [self.ZeroMap] * len(self.Classes)
        self.Ft_ = [self.ZeroMap] * len(self.Classes)

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
        :type self: object
        """

        # TODO: change rainfall .tss files into grids
        self.wf_updateparameters()  # read the temperature map for each step (see parameters())
        # self.logger.debug("Step: "+str(int(self.thestep + self._d_firstTimeStep))+"/"+str(int(self._d_nrTimeSteps)))
        self.thestep = self.thestep + 1

	#if self.thestep == 26:
#		pdb.set_trace()

        self.Si_t = copylist(self.Si)
        self.Sw_t = copylist(self.Sw)
        self.Su_t = copylist(self.Su)
        self.Sa_t = copylist(self.Sa)
        self.Sf_t = copylist(self.Sf)
        self.Sfa_t = copylist(self.Sfa)
        self.Ss_t = self.Ss
        self.trackQ_t = copylist(self.trackQ)  # copylist(self.trackQ)
        self.convQu_t = [copylist(self.convQu[i]) for i in self.Classes]  # copylist(self.convQu)
        self.convQa_t = [copylist(self.convQa[i]) for i in self.Classes]
        
        if self.IRURFR_L:
            self.PotEvaporation = areatotal(self.PotEvaporation * self.percentArea, nominal(self.TopoId))
            self.Precipitation = areatotal(self.Precipitation * self.percentArea, nominal(self.TopoId))
            self.Temperature = areaaverage(self.Temperature * self.percentArea, nominal(self.TopoId))
        
        self.PrecipTotal = self.Precipitation
        if self.selectSw[0] > 0:
            self.Precipitation = ifthenelse(self.Temperature >= self.Tt[0], self.PrecipTotal,0)
            self.PrecipitationSnow = ifthenelse(self.Temperature < self.Tt[0], self.PrecipTotal,0)

        self.EpDay2 = self.EpDay
        self.EpDaySnow2 = self.EpDaySnow
        
        #if self.thestep >= 45:
        	#pdb.set_trace()
        
        for k in self.Classes:

            # SNOW =================================================================================================
            if self.selectSw[k]:
                eval_str = 'reservoir_Sw.{:s}(self, k)'.format(self.selectSw[k])
            else:
                eval_str = 'reservoir_Sw.snow_no_reservoir(self, k)'
            eval(eval_str)


            # INTERCEPTION =========================================================================================
            if self.selectSi[k]:
                eval_str = 'reservoir_Si.{:s}(self, k)'.format(self.selectSi[k])
            else:
                eval_str = 'reservoir_Si.interception_no_reservoir(self, k)'
            eval(eval_str)
            
            # AGRICULTURE ZONE ======================================================================================           
            if self.selectSa[k]:
                eval_str = 'reservoir_Sa.{:s}(self, k)'.format(self.selectSa[k])
            else:
                eval_str = 'reservoir_Sa.agriZone_no_reservoir(self, k)'
            eval(eval_str)            
            
            # UNSATURATED ZONE ======================================================================================
            if self.selectSu[k]:
                eval_str = 'reservoir_Su.{:s}(self, k)'.format(self.selectSu[k])
            else:
                eval_str = 'reservoir_Su.unsatZone_no_reservoir(self, k)'
            eval(eval_str)


            # FAST RUNOFF RESERVOIR ===================================================================================
            if self.selectSf[k]:
                eval_str = 'reservoir_Sf.{:s}(self, k)'.format(self.selectSf[k])
            else:
                eval_str = 'reservoir_Sf.fastRunoff_no_reservoir(self, k)'
            eval(eval_str)

            #FAST AGRICULTURE DITCHES RUNOFF RESERVOIR ===================================================================================
            if self.selectSfa[k]:
                eval_str = 'reservoir_Sf.{:s}(self, k)'.format(self.selectSfa[k])
            else:
                eval_str = 'reservoir_Sf.fastAgriRunoff_no_reservoir(self, k)'
            eval(eval_str)

        
        # TOTAL RUNOFF =============================================================================================
        self.Qftotal = sum([x * y for x, y in zip(self.Qf_, self.percent)]) + sum([x*y for x,y in zip(self.Qfa_,self.percent)])

        # SLOW RUNOFF RESERVOIR ===========================================================================
        if self.selectSs:
            eval_str = 'reservoir_Ss.{:s}(self)'.format(self.selectSs)
        else:
            eval_str = 'reservoir_Ss.groundWater_no_reservoir(self)'
        eval(eval_str)
        
        # ROUTING
        if self.selectRout:
            eval_str = 'reservoir_Sf.{:s}(self)'.format(self.selectRout)
        else:
            eval_str = 'reservoir_Sf.noRouting(self)'
        eval(eval_str)        


        # WATER BALANCE (per reservoir, per cell) ========================================================================================
        self.QtlagWB = (self.Qtlag / self.surfaceArea) * 1000 * self.timestepsecs
        self.convQuWB = [sum(self.convQu[i]) for i in self.Classes]
        self.convQuWB_t = [sum(self.convQu_t[i]) for i in self.Classes]
        self.convQaWB = [sum(self.convQa[i]) for i in self.Classes]
        self.convQaWB_t = [sum(self.convQa_t[i]) for i in self.Classes]
        self.trackQWB = (sum(self.trackQ) / self.surfaceArea) * 1000
        self.trackQWB_t = (sum(self.trackQ_t) / self.surfaceArea) * 1000
        self.WB = self.Precipitation - sum(multiply(self.Ei_, self.percent)) - sum(
            multiply(self.Eu_, self.percent)) - self.QtlagWB - sum(
            multiply(self.Si, self.percent)) + sum(multiply(self.Si_t, self.percent)) - sum(
            multiply(self.Su, self.percent)) + sum(multiply(self.Su_t, self.percent)) - sum(
            multiply(self.Sf, self.percent)) + sum(multiply(self.Sf_t, self.percent)) - sum(
            multiply(self.Ss, self.percent)) + sum(
            multiply(self.Ss_t, self.percent)) - self.trackQWB + self.trackQWB_t - sum(
            multiply(self.convQuWB, self.percent)) + sum(multiply(self.convQuWB_t, self.percent))

        #    #fuxes and states in m3/h
        self.P = areatotal(self.PrecipTotal / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Ei = areatotal(sum(multiply(self.Ei_,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Ea = areatotal(sum(multiply(self.Ea_,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Eu = areatotal(sum(multiply(self.Eu_,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Ew = areatotal(sum(multiply(self.Ew_,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.EwiCorr = areatotal(sum(multiply(multiply(self.Ew_, self.lamdaS / self.lamda), self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Qtot = self.QLagTot * self.timestepsecs
        self.SiWB = areatotal(sum(multiply(self.Si,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Si_WB = areatotal(sum(multiply(self.Si_t,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.SuWB = areatotal(sum(multiply(self.Su,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Su_WB = areatotal(sum(multiply(self.Su_t,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.SaWB = areatotal(sum(multiply(self.Sa,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Sa_WB = areatotal(sum(multiply(self.Sa_t,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.SfWB = areatotal(sum(multiply(self.Sf,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Sf_WB = areatotal(sum(multiply(self.Sf_t,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.SfaWB = areatotal(sum(multiply(self.Sfa,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Sfa_WB = areatotal(sum(multiply(self.Sfa_t,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.SwWB = areatotal(sum(multiply(self.Sw,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Sw_WB = areatotal(sum(multiply(self.Sw_t,self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.SsWB = areatotal(self.Ss / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.Ss_WB = areatotal(self.Ss_t / 1000 * self.surfaceArea,nominal(self.TopoId))    
        self.convQuWB = areatotal(sum(multiply([sum(self.convQu[i]) for i in self.Classes],self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.convQu_WB = areatotal(sum(multiply([sum(self.convQu_t[i]) for i in self.Classes],self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.convQaWB = areatotal(sum(multiply([sum(self.convQa[i]) for i in self.Classes],self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.convQa_WB = areatotal(sum(multiply([sum(self.convQa_t[i]) for i in self.Classes],self.percent)) / 1000 * self.surfaceArea,nominal(self.TopoId))
        self.trackQWB = areatotal(sum(self.trackQ),nominal(self.TopoId))
        self.trackQ_WB = areatotal(sum(self.trackQ_t),nominal(self.TopoId)) 
        self.QstateWB = areatotal(sum(self.Qstate) * self.timestepsecs, nominal(self.TopoId))
        self.Qstate_WB = areatotal(sum(self.Qstate_t) * self.timestepsecs, nominal(self.TopoId))
#        self.QstateWB = areatotal(sum(self.Qstate) * 0.0405, nominal(self.TopoId))
#        self.Qstate_WB = areatotal(sum(self.Qstate_t) * 0.0405, nominal(self.TopoId))
#        self.QstateWB = areatotal(self.Qstate, nominal(self.TopoId))
#        self.Qstate_WB = areatotal(self.Qstate_t, nominal(self.TopoId))
#        
        #WBtot in m3/s
        self.WBtot = (self.P - self.Ei + self.EwiCorr - self.Ew - self.Ea - self.Eu - self.Qtot - self.SiWB + self.Si_WB - self.SuWB + self.Su_WB - self.SaWB + self.Sa_WB - self.SwWB + self.Sw_WB - self.SfWB + self.Sf_WB - self.SfaWB + self.Sfa_WB - self.SsWB + self.Ss_WB - self.convQuWB +self.convQu_WB - self.convQaWB +self.convQa_WB - self.trackQWB + self.trackQ_WB - self.QstateWB + self.Qstate_WB) / self.timestepsecs     

        # SUMMED FLUXES ======================================================================================
        self.sumprecip = self.sumprecip + self.Precipitation  # accumulated rainfall for water balance (m/h)
        self.sumevap = self.sumevap + sum(multiply(self.Ei_, self.percent)) + sum(
            multiply(self.Eu_, self.percent)) + sum(
            multiply(self.Ea_, self.percent)) + sum(
            multiply(self.Ew_, self.percent)) # accumulated evaporation for water balance (m/h)
        try:
            self.sumpotevap = self.sumpotevap + self.PotEvaporation  # accumulated potential evaporation (m/h)
        except:
            self.sumpotevap = self.EpHour
        self.sumrunoff = self.sumrunoff + self.Qtlag * 1000 * self.timestepsecs / self.surfaceArea  # accumulated runoff for water balance (m/h)
        self.sumwb = self.sumwb + self.WB

        self.sumE = sum(multiply(self.Ei_, self.percent)) + sum(multiply(self.Eu_, self.percent))


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
    configfile = "wflow_topoflex.ini"
    LogFileName="wflow.log" 
    _lastTimeStep = 10
    _firstTimeStep = 1
    fewsrun=False
    runinfoFile="runinfo.xml"
    timestepsecs = 86400
    wflow_cloneMap = 'wflow_subcatch.map'
    NoOverWrite=1
    loglevel = logging.DEBUG
	

    # This allows us to use the model both on the command line and to call 
    # the model usinge main function from another python script.

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return

    opts, args = getopt.getopt(argv, 'C:S:T:Ic:s:R:F:fl:L:P:p:')

    for o, a in opts:
        if o == '-F': 
            runinfoFile = a
            fewsrun = True
        if o == '-C': caseName = a
        if o == '-R': runId = a
        if o == '-c': configfile = a
        if o == '-s': timestepsecs = int(a)
        if o == '-T': _lastTimeStep = int(a)
        if o == '-S': _firstTimeStep = int(a)
        if o == '-f': NoOverWrite = 0
        if o == '-L': LogFileName = a 
        if o == '-l': exec "loglevel = logging." + a
    if (len(opts) <= 1):
        usage()
        
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
    dynModelFw.createRunId(NoOverWrite=NoOverWrite,logfname=LogFileName,level=loglevel,model='wflow_topoflex',doSetupFramework=False)


    for o, a in opts:
        if o == '-P':
            left = a.split('=')[0]
            right = a.split('=')[1]
            configset(myModel.config,'variable_change_once',left,right,overwrite=True)
        if o == '-p':
            left = a.split('=')[0]
            right = a.split('=')[1]
            configset(myModel.config,'variable_change_timestep',left,right,overwrite=True)
        if o == '-X': configset(myModel.config, 'model', 'OverWriteInit', '1', overwrite=True)
        if o == '-I': configset(myModel.config, 'model', 'reinit', '1', overwrite=True)
        if o == '-i': configset(myModel.config, 'model', 'intbl', a, overwrite=True)
        if o == '-s': configset(myModel.config, 'model', 'timestepsecs', a, overwrite=True)
        if o == '-x': configset(myModel.config, 'model', 'sCatch', a, overwrite=True)
        if o == '-c': configset(myModel.config, 'model', 'configfile', a, overwrite=True)
        if o == '-M': configset(myModel.config, 'model', 'MassWasting', "0", overwrite=True)
        if o == '-Q': configset(myModel.config, 'model', 'ExternalQbase', '1', overwrite=True)
        if o == '-U':
            configset(myModel.config, 'model', 'updateFile', a, overwrite=True)
            configset(myModel.config, 'model', 'updating', "1", overwrite=True)
        if o == '-u':
            zz = []
            exec "zz =" + a
            updateCols = zz
        if o == '-E': configset(myModel.config, 'model', 'reInfilt', '1', overwrite=True)
        if o == '-R': runId = a
        if o == '-W': configset(myModel.config, 'model', 'waterdem', '1', overwrite=True)

    dynModelFw.setupFramework()
    dynModelFw._runInitial()
    dynModelFw._runResume()
    dynModelFw._runDynamic(_firstTimeStep, _lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()


if __name__ == "__main__":
    main()
