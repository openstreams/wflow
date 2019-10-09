#!/usr/bin/python

# Wflow is Free software, see below:
#
# Copyright (c) J. Schellekens 2005-2011
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


"""
Run the wflow_cqf hydrological model..

wflow_cqf is a model develope oftare the cqflow model and is developed 
specifically for the rio chiquito catchment

usage

::
    
    wflow_cqf [-h][-v level][-F runinfofile][-L logfile][-C casename][-R runId]
          [-c configfile][-T last_step][-S first_step][-s seconds][-W][-E][-N][-U discharge]
          [-P parameter multiplication][-X][-f][-I][-i tbl_dir][-x subcatchId][-u updatecols]
          [-p inputparameter multiplication]
          

    -X: save state at the end of the run over the initial conditions at the start        
    
    -f: Force overwrite of existing results    
    
    -T: Set last timestep
    
    -S: Set the start timestep (default = 1)
    
    -N: No lateral flow, use runoff response function to generate fast runoff
    
    -s: Set the model timesteps in seconds
    
    -I: re-initialize the initial model conditions with default
    
    -i: Set input table directory (default is intbl)
    
    -x: run for subcatchment only (e.g. -x 1)
    
    -C: set the name  of the case (directory) to run
    
    -R: set the name runId within the current case
    
    -L: set the logfile
    
    -E: Switch on reinfiltration of overland flow
    
    -c: name of wflow the configuration file (default: Casename/wflow_cqf.ini). 
    
    -h: print usage information
    
    -W: If set, this flag indicates that an ldd is created for the water level
        for each timestep. If not the water is assumed to flow according to the 
        DEM. Wflow will run a lot slower with this option. Most of the time
        (shallow soil, steep topography) you do not need this option. Also, if you 
        need it you migth actually need another model.
        
    -U: The argument to this option should be a .tss file with measured discharge in
        [m^3/s] which the program will use to update the internal state to match 
        the measured flow. The number of columns in this file should match the 
        number of gauges.
        
    -u: list of gauges/columns to use in update. Format:
        -u [1 , 4 ,13]
        The above example uses column 1, 4 and 13
        Note that this also sets the order in which the updating takes place! In
        general specify downstream gauges first.
        
    -P: set parameter multiply dictionary (e.g: -P {'self.FirstZoneDepth' : 1.2}
        to increase self.FirstZoneDepth by 20%, multiply with 1.2)
        
    -p: set input parameter (dynamic, e.g. precip) multiply dictionary 
        (e.g: -p {'self.Precipitation' : 1.2} to increase Precipitation 
        by 20%, multiply with 1.2)    
        
    -v: set verbosity level


$Author: schelle $
$Id: wflow_cqf.py 909 2014-01-16 15:23:32Z schelle $
$Rev: 909 $
"""

import os.path

import pcraster.framework
from wflow.wf_DynamicFramework import *
from wflow.wflow_adapt import *
from wflow.wflow_funcs import *

wflow = "wflow_cqf: "
wflowVersion = "$Revision: 909 $  $Date: 2014-01-16 16:23:32 +0100 (Thu, 16 Jan 2014) $"

updateCols = []


# Dictionary with parameters and multipliers (used in calibration)
multpars = {}
multdynapars = {}


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args:
        print(msg)
    print(__doc__)
    sys.exit(0)


def actEvap_SBM(RootingDepth, WTable, UStoreDepth, FirstZoneDepth, PotTrans, smoothpar):
    """
    Actual evaporation function:
        
	- first try to get demand from the saturated zone, using the rootingdepth as a limiting factor
	- secondly try to get the remaining water from the unsaturated store

    Input: 
        - RootingDepth,WTable, UStoreDepth,FirstZoneDepth, PotTrans
        
    Output: 
        - ActEvap,  FirstZoneDepth,  UStoreDepth ActEvapUStore
        
    .. todo::
        
        add option to take length of roots in saturated zone into account
        
    """

    # Step 1 from saturated zone, use rootingDepth as a limiting factor
    # rootsinWater = WTable < RootingDepth
    # ActEvapSat = pcr.ifthenelse(rootsinWater,pcr.min(PotTrans,FirstZoneDepth),0.0)
    # new method:
    # use sCurve to determine if the roots are wet.At the moment this ise set
    # to be a 0-1 curve
    wetroots = sCurve(WTable, a=RootingDepth, c=smoothpar)
    ActEvapSat = pcr.min(PotTrans * wetroots, FirstZoneDepth)

    FirstZoneDepth = FirstZoneDepth - ActEvapSat
    RestPotEvap = PotTrans - ActEvapSat

    # now try unsat store
    AvailCap = pcr.min(
        1.0, pcr.max(0.0, (WTable - RootingDepth) / (RootingDepth + 1.0))
    )

    # AvailCap = pcr.max(0.0,pcr.ifthenelse(WTable < RootingDepth,  WTable/RootingDepth,  RootingDepth/WTable))
    MaxExtr = AvailCap * UStoreDepth
    ActEvapUStore = pcr.min(MaxExtr, RestPotEvap, UStoreDepth)
    UStoreDepth = UStoreDepth - ActEvapUStore

    ActEvap = ActEvapSat + ActEvapUStore

    return ActEvap, FirstZoneDepth, UStoreDepth, ActEvapUStore


class WflowModel(pcraster.framework.DynamicModel):
    def __init__(self, cloneMap, Dir, RunDir, configfile):
        pcraster.framework.DynamicModel.__init__(self)
        pcr.setclone(Dir + "/staticmaps/" + cloneMap)
        self.runId = RunDir
        self.caseName = Dir
        self.Dir = Dir + "/"
        self.configfile = configfile

    def updateRunOff(self):
        """
      Updates the kinematic wave reservoir. Should be run after updates to Q
      """
        self.WaterLevel = (self.Alpha * pow(self.SurfaceRunoff, self.Beta)) / self.Bw
        # wetted perimeter (m)
        P = self.Bw + (2 * self.WaterLevel)
        # Alpha
        self.Alpha = self.AlpTerm * pow(P, self.AlpPow)
        self.OldKinWaveVolume = self.KinWaveVolume
        self.KinWaveVolume = self.WaterLevel * self.Bw * self.DCL

    def stateVariables(self):
        """ 
      returns a list of state variables that are essential to the model. 
      This list is essential for the resume and suspend functions to work.
      
      This function is specific for each model and **must** be present.
      
      - CanopyStorage is any needed for subdaily steps
      """
        states = [
            "SurfaceRunoff",
            "WaterLevel",
            "FirstZoneDepth",
            "UStoreDepth",
            "CanopyStorage",
        ]

        return states

    def supplyCurrentTime(self):
        """
      gets the current time in seconds after the start of the run
      """
        return self.currentTimeStep() * self.timestepsecs

    def readtblDefault(self, pathtotbl, landuse, subcatch, soil, default):
        """
    First check if a prepared map of the same name is present
    in the staticmaps directory. next try to
    read a tbl file to match a landuse, catchment and soil map. Returns 
    the default value if the tbl file is not found.
    
    Input: 
        - pathtotbl,landuse,subcatch,soil, default
        
    Output: 
        - map constructed from tbl file or map with default value
    """

        mapname = (
            os.path.dirname(pathtotbl)
            + "/../staticmaps/"
            + os.path.splitext(os.path.basename(pathtotbl))[0]
            + ".map"
        )
        if os.path.exists(mapname):
            self.logger.info("reading map parameter file: " + mapname)
            rest = pcr.cover(pcr.readmap(mapname), default)
        else:
            if os.path.isfile(pathtotbl):
                rest = pcr.lookupscalar(pathtotbl, landuse, subcatch, soil)  #
                self.logger.info("Creating map from table: " + pathtotbl)
            else:
                self.logger.warning(
                    "tbl file not found ("
                    + pathtotbl
                    + ") returning default value: "
                    + str(default)
                )
                rest = pcr.scalar(default)

        return rest

    def suspend(self):

        self.logger.info("Saving initial conditions...")
        self.wf_suspend(self.SaveDir + "/outstate/")

        if self.OverWriteInit:
            self.logger.info("Saving initial conditions over start conditions...")
            self.wf_suspend(self.SaveDir + "/instate/")

        pcr.report(self.CumInwaterMM, self.SaveDir + "/outsum/CumInwaterMM.map")
        pcr.report(self.CumReinfilt, self.SaveDir + "/outsum/CumReinfilt.map")
        pcr.report(self.CumPrec, self.SaveDir + "/outsum/CumPrec.map")
        pcr.report(self.CumEvap, self.SaveDir + "/outsum/CumEvap.map")
        pcr.report(self.CumInt, self.SaveDir + "/outsum/CumInt.map")
        pcr.report(self.CumLeakage, self.SaveDir + "/outsum/CumLeakage.map")
        pcr.report(self.CumPotenEvap, self.SaveDir + "/outsum/CumPotenEvap.map")
        pcr.report(self.CumExfiltWater, self.SaveDir + "/outsum/CumExfiltWater.map")
        pcr.report(self.watbal, self.SaveDir + "/outsum/watbal.map")

    def initial(self):

        """Initial part of the model, executed only once """
        global statistics
        global multpars

        self.thestep = pcr.scalar(0)
        self.basetimestep = 86400
        self.SSSF = False
        pcr.setglobaloption("unittrue")
        intbl = "intbl"
        self.precipTss = "/intss/P.tss"
        self.evapTss = "/intss/PET.tss"
        self.tempTss = "/intss/T.tss"
        self.inflowTss = "/intss/Inflow.tss"
        self.SeepageTss = "/intss/Seepage.tss"

        self.logger.info("running for " + str(self.nrTimeSteps()) + " timesteps")

        # Set and get defaults from ConfigFile here ###################################
        self.scalarInput = int(configget(self.config, "model", "ScalarInput", "0"))
        self.Tslice = int(configget(self.config, "model", "Tslice", "1"))
        self.interpolMethod = configget(
            self.config, "model", "InterpolationMethod", "inv"
        )
        self.reinit = int(configget(self.config, "run", "reinit", "0"))
        self.OverWriteInit = int(configget(self.config, "model", "OverWriteInit", "0"))
        self.updating = int(configget(self.config, "model", "updating", "0"))
        self.updateFile = configget(self.config, "model", "updateFile", "no_set")

        self.sCatch = int(configget(self.config, "model", "sCatch", "0"))
        self.intbl = configget(self.config, "model", "intbl", "intbl")
        self.timestepsecs = int(
            configget(self.config, "model", "timestepsecs", "86400")
        )
        self.modelSnow = int(configget(self.config, "model", "ModelSnow", "1"))
        sizeinmetres = int(configget(self.config, "layout", "sizeinmetres", "0"))
        alf = float(configget(self.config, "model", "Alpha", "60"))
        # TODO: make this into a list for all gauges or a map
        Qmax = float(configget(self.config, "model", "AnnualDischarge", "300"))
        self.UpdMaxDist = float(configget(self.config, "model", "UpdMaxDist", "100"))
        self.ExternalQbase = int(configget(self.config, "model", "ExternalQbase", "0"))
        self.waterdem = int(configget(self.config, "model", "waterdem", "0"))
        WIMaxScale = float(configget(self.config, "model", "WIMaxScale", "0.8"))
        self.reInfilt = int(configget(self.config, "model", "reInfilt", "0"))

        # static maps to use (normally default)
        wflow_subcatch = configget(
            self.config, "model", "wflow_subcatch", "staticmaps/wflow_subcatch.map"
        )
        wflow_dem = configget(
            self.config, "model", "wflow_dem", "staticmaps/wflow_dem.map"
        )
        wflow_ldd = configget(
            self.config, "model", "wflow_ldd", "staticmaps/wflow_ldd.map"
        )
        wflow_river = configget(
            self.config, "model", "wflow_river", "staticmaps/wflow_river.map"
        )
        wflow_riverlength = configget(
            self.config,
            "model",
            "wflow_riverlength",
            "staticmaps/wflow_riverlength.map",
        )
        wflow_riverlength_fact = configget(
            self.config,
            "model",
            "wflow_riverlength_fact",
            "staticmaps/wflow_riverlength_fact.map",
        )
        wflow_landuse = configget(
            self.config, "model", "wflow_landuse", "staticmaps/wflow_landuse.map"
        )
        wflow_soil = configget(
            self.config, "model", "wflow_soil", "staticmaps/wflow_soil.map"
        )
        wflow_gauges = configget(
            self.config, "model", "wflow_gauges", "staticmaps/wflow_gauges.map"
        )
        wflow_inflow = configget(
            self.config, "model", "wflow_inflow", "staticmaps/wflow_inflow.map"
        )
        wflow_mgauges = configget(
            self.config, "model", "wflow_mgauges", "staticmaps/wflow_mgauges.map"
        )

        # 2: Input base maps ########################################################
        subcatch = pcr.ordinal(
            pcr.readmap(self.Dir + wflow_subcatch)
        )  # Determines the area of calculations (all cells > 0)
        subcatch = pcr.ifthen(subcatch > 0, subcatch)
        if self.sCatch > 0:
            subcatch = pcr.ifthen(subcatch == sCatch, subcatch)

        self.Altitude = pcr.readmap(
            self.Dir + wflow_dem
        )  # * pcr.scalar(pcr.defined(subcatch)) # DEM
        self.TopoLdd = pcr.readmap(self.Dir + wflow_ldd)  # Local
        self.TopoId = pcr.ordinal(pcr.readmap(self.Dir + wflow_subcatch))  # area map
        self.River = pcr.cover(pcr.boolean(pcr.readmap(self.Dir + wflow_river)), 0)
        self.RiverLength = pcrut.readmapSave(self.Dir + wflow_riverlength, 0.0)
        # Factor to multiply riverlength with (defaults to 1.0)
        self.RiverLengthFac = pcrut.readmapSave(self.Dir + wflow_riverlength_fact, 1.0)

        # read landuse and soilmap and make sure there are no missing points related to the
        # subcatchment map. Currently sets the lu and soil type  type to 1
        self.LandUse = pcr.readmap(self.Dir + wflow_landuse)
        self.LandUse = pcr.cover(self.LandUse, pcr.nominal(pcr.ordinal(subcatch) > 0))
        self.Soil = pcr.readmap(self.Dir + wflow_soil)
        self.Soil = pcr.cover(self.Soil, pcr.nominal(pcr.ordinal(subcatch) > 0))
        self.OutputLoc = pcr.ordinal(
            pcr.readmap(self.Dir + wflow_gauges)
        )  # location of output gauge(s)
        self.InflowLoc = pcrut.readmapSave(
            self.Dir + wflow_inflow, 0.0
        )  # location abstractions/inflows.
        self.SeepageLoc = pcrut.readmapSave(
            self.Dir + wflow_inflow, 0.0
        )  # location abstractions/inflows.

        # Experimental
        self.RunoffGenSigmaFunction = int(
            configget(self.config, "model", "RunoffGenSigmaFunction", "0")
        )
        self.RunoffGeneratingGWPerc = float(
            configget(self.config, "defaultfortbl", "RunoffGeneratingGWPerc", "0.1")
        )
        self.RunoffGeneratingThickness = float(
            configget(self.config, "defaultfortbl", "RunoffGeneratingThickness", "0.0")
        )

        if self.scalarInput:
            self.gaugesMap = pcr.readmap(
                self.Dir + wflow_mgauges
            )  # location of rainfall/evap/temp gauge(s)
        self.OutputId = pcr.ordinal(
            pcr.readmap(self.Dir + wflow_subcatch)
        )  # location of subcatchment
        # Temperature correction poer cell to add

        self.TempCor = pcrut.readmapSave(
            self.Dir
            + configget(
                self.config,
                "model",
                "TemperatureCorrectionMap",
                "staticmaps/wflow_tempcor.map",
            ),
            0.0,
        )

        self.ZeroMap = 0.0 * pcr.scalar(subcatch)  # map with only zero's

        # 3: Input time series ###################################################
        self.P_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Precipitation", "/inmaps/P"
        )  # timeseries for rainfall
        self.PET_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "EvapoTranspiration", "/inmaps/PET"
        )  # timeseries for rainfall"/inmaps/PET"          # potential evapotranspiration
        self.TEMP_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Temperature", "/inmaps/TEMP"
        )  # timeseries for rainfall "/inmaps/TEMP"          # global radiation
        self.Inflow_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Inflow", "/inmaps/IF"
        )  # timeseries for rainfall "/inmaps/IF" # in/outflow locations (abstractions)
        self.RH_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "RH", "/inmaps/RH"
        )  # timeseries for rainfall "/inmaps/SE" # in/outflow locations (abstractions)
        self.WindSpeed_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "WindSpeed", "/inmaps/wins"
        )  # timeseries for rainfall "/inmaps/SE" # in/outflow locations (abstractions)
        self.RAD_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Radiation", "/inmaps/RAD"
        )  # timeseries for rainfall "/inmaps/SE" # in/outflow locations (abstractions)
        self.Seepage_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Seepage", "/inmaps/SE"
        )  # timeseries for rainfall "/inmaps/SE" # in/outflow locations (abstractions)
        self.HP_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "HP", "/inmaps/HP"
        )  # timeseries for rainfall "/inmaps/SE" # in/outflow locations (abstractions)
        self.WaterCatch_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "WaterCatch", "/inmaps/WC"
        )  # timeseries for rainfall "/inmaps/SE" # in/outflow locations (abstractions)
        # 3: Input time series ###################################################

        # Set static initial values here #########################################
        self.SoilAlbedo = 0.1  # Not used at the moment
        self.pi = 3.1416
        self.e = 2.7183
        self.SScale = 100.0

        self.Latitude = pcr.ycoordinate(pcr.boolean(self.Altitude))
        self.Longitude = pcr.xcoordinate(pcr.boolean(self.Altitude))

        self.logger.info("Linking parameters to landuse, catchment and soil...")
        self.RunoffGeneratingGWPerc = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/RunoffGeneratingGWPerc.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            self.RunoffGeneratingGWPerc,
        )
        self.RunoffGeneratingThickness = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/RunoffGeneratingThickness.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            self.RunoffGeneratingThickness,
        )
        self.Cmax = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/MaxCanopyStorage.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            1.0,
        )
        self.EoverR = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/EoverR.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.1,
        )
        # self.Albedo=pcr.lookupscalar(self.Dir + "\intbl\Albedo.tbl",self.LandUse) # Not used anymore
        self.CanopyGapFraction = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/CanopyGapFraction.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.1,
        )
        self.RootingDepth = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/RootingDepth.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            750.0,
        )  # rooting depth
        #: rootdistpar determien how roots are linked to water table.The number shoudl be negative. A high number means that all roots are wet if
        #: the water table is above the lowest part of the roots. A lower number smooths this.
        self.rootdistpar = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/rootdistpar.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            -80000.0,
        )  # rrootdistpar

        # Soil parameters
        # infiltration capacity if the soil [mm/day]
        self.InfiltCapSoil = (
            self.readtblDefault(
                self.Dir + "/" + self.intbl + "/InfiltCapSoil.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                100.0,
            )
            * self.timestepsecs
            / self.basetimestep
        )
        self.CapScale = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/CapScale.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            100.0,
        )  #
        # infiltration capacity of the compacted
        self.InfiltCapPath = (
            self.readtblDefault(
                self.Dir + "/" + self.intbl + "/InfiltCapPath.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                10.0,
            )
            * self.timestepsecs
            / self.basetimestep
        )
        self.MaxLeakage = (
            self.readtblDefault(
                self.Dir + "/" + self.intbl + "/MaxLeakage.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                0.0,
            )
            * self.timestepsecs
            / self.basetimestep
        )
        # areas (paths) in [mm/day]
        # Fraction area with compacted soil (Paths etc.)
        self.PathFrac = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/PathFrac.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.01,
        )
        # thickness of the soil
        self.FirstZoneThickness = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/FirstZoneCapacity.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            2000.0,
        )
        self.thetaR = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/thetaR.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.01,
        )
        self.thetaS = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/thetaS.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.6,
        )
        # minimum thickness of soild
        self.FirstZoneMinCapacity = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/FirstZoneMinCapacity.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            500.0,
        )

        # FirstZoneKsatVer = $2\inmaps\FirstZoneKsatVer.map
        self.FirstZoneKsatVer = (
            self.readtblDefault(
                self.Dir + "/" + self.intbl + "/FirstZoneKsatVer.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                3000.0,
            )
            * self.timestepsecs
            / self.basetimestep
        )
        self.Beta = pcr.scalar(0.6)  # For sheetflow

        self.M = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/M.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            300.0,
        )  # Decay parameter in Topog_cqf
        self.N = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/N.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.072,
        )  # Manning overland flow
        self.NRiver = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/N_River.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.036,
        )  # Manning river
        self.WaterFrac = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/WaterFrac.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.0,
        )  # Fraction Open water

        # cqflow specific stuff
        self.Albedo = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/Albedo.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.18,
        )  #
        self.LeafAreaIndex = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/LeafAreaIndex.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            4.0,
        )  #
        self.WindSpeedHeigth = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/WindSpeedHeigth.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            4.0,
        )  #
        self.VegetationHeigth = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/VegetationHeigth.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            4.0,
        )  #

        self.xl, self.yl, self.reallength = pcrut.detRealCellLength(
            self.ZeroMap, sizeinmetres
        )
        self.Slope = pcr.slope(self.Altitude)
        # self.Slope=pcr.ifthen(pcr.boolean(self.TopoId),pcr.max(0.001,self.Slope*celllength()/self.reallength))
        self.Slope = pcr.max(0.001, self.Slope * pcr.celllength() / self.reallength)
        Terrain_angle = pcr.scalar(pcr.atan(self.Slope))

        # Multiply parameters with a factor (for calibration etc) -P option in command line
        for k, v in multpars.items():
            estr = k + "=" + k + "*" + str(v)
            self.logger.info("Parameter multiplication: " + estr)
            exec(estr)

        self.N = pcr.ifthenelse(self.River, self.NRiver, self.N)

        # Determine river width from DEM, upstream area and yearly average discharge
        # Scale yearly average Q at outlet with upstream are to get Q over whole catchment
        # Alf ranges from 5 to > 60. 5 for hardrock. large values for sediments
        # "Noah J. Finnegan et al 2005 Controls on the channel width of rivers:
        # Implications for modeling fluvial incision of bedrock"

        upstr = pcr.catchmenttotal(1, self.TopoLdd)
        Qscale = upstr / pcr.mapmaximum(upstr) * Qmax
        W = (
            (alf * (alf + 2.0) ** (0.6666666667)) ** (0.375)
            * Qscale ** (0.375)
            * (pcr.max(0.0001, pcr.windowaverage(self.Slope, pcr.celllength() * 4.0)))
            ** (-0.1875)
            * self.N ** (0.375)
        )
        RiverWidth = W

        # soil thickness based on topographical index (see Environmental modelling: finding simplicity in complexity)
        # 1: calculate wetness index
        # 2: Scale the capacity (now actually a max capacity) based on the index, also apply a minmum capacity
        WI = pcr.ln(
            pcr.accuflux(self.TopoLdd, 1) / self.Slope
        )  # Topographical wetnesss. Scale WI by zone/subcatchment assuming these ara also geological units
        WIMax = pcr.areamaximum(WI, self.TopoId) * WIMaxScale
        self.FirstZoneThickness = pcr.max(
            pcr.min(self.FirstZoneThickness, (WI / WIMax) * self.FirstZoneThickness),
            self.FirstZoneMinCapacity,
        )

        self.FirstZoneCapacity = self.FirstZoneThickness * (self.thetaS - self.thetaR)

        # limit roots to top 99% of first zone
        self.RootingDepth = pcr.min(self.FirstZoneThickness * 0.99, self.RootingDepth)

        # subgrid runoff generation
        self.DemMax = pcr.readmap(self.Dir + "/staticmaps/wflow_demmax")
        self.DrainageBase = pcr.readmap(self.Dir + "/staticmaps/wflow_demmin")
        self.CC = pcr.min(
            100.0,
            -log(1.0 / 0.1 - 1) / pcr.min(-0.1, self.DrainageBase - self.Altitude),
        )

        # if pcr.maptotal(self.RunoffGeneratingThickness <= 0.0):
        self.GWScale = (
            (self.DemMax - self.DrainageBase)
            / self.FirstZoneThickness
            / self.RunoffGeneratingGWPerc
        )
        # else:
        #    self.GWScale = (self.DemMax-self.DrainageBase)/min(self.RunoffGeneratingThickness, self.FirstZoneThickness)

        # Which columns/gauges to use/ignore in updating
        self.UpdateMap = self.ZeroMap

        if self.updating:
            touse = numpy.zeros(gaugear.shape, dtype="int")

            for thecol in updateCols:
                idx = (gaugear == thecol).nonzero()
                touse[idx] = thecol

            self.UpdateMap = pcr.numpy2pcr(pcr.Nominal, touse, 0.0)
            # Calulate distance to updating points (upstream) annd use to scale the correction
            # ldddist returns zero for cell at the gauges so add 1.0 tp result
            self.DistToUpdPt = pcr.cover(
                pcr.min(
                    ldddist(self.TopoLdd, pcr.boolean(pcr.cover(self.UpdateMap, 0)), 1)
                    * self.reallength
                    / pcr.celllength(),
                    self.UpdMaxDist,
                ),
                self.UpdMaxDist,
            )

        # Initializing of variables
        self.logger.info("Initializing of model variables..")
        self.TopoLdd = pcr.lddmask(self.TopoLdd, pcr.boolean(self.TopoId))
        catchmentcells = pcr.maptotal(pcr.scalar(self.TopoId))

        # Used to seperate output per LandUse/management classes
        OutZones = self.LandUse

        self.QMMConv = self.timestepsecs / (
            self.reallength * self.reallength * 0.001
        )  # m3/s --> mm
        self.ToCubic = (
            self.reallength * self.reallength * 0.001
        ) / self.timestepsecs  # m3/s
        self.KinWaveVolume = self.ZeroMap
        self.OldKinWaveVolume = self.ZeroMap
        self.sumprecip = self.ZeroMap  # accumulated rainfall for water balance
        self.sumevap = self.ZeroMap  # accumulated evaporation for water balance
        self.sumrunoff = self.ZeroMap  # accumulated runoff for water balance
        self.sumint = self.ZeroMap  # accumulated interception for water balance
        self.sumleakage = self.ZeroMap
        self.CumReinfilt = self.ZeroMap
        self.sumoutflow = self.ZeroMap
        self.sumsnowmelt = self.ZeroMap
        self.CumRad = self.ZeroMap
        self.SnowMelt = self.ZeroMap
        self.CumPrec = self.ZeroMap
        self.CumInwaterMM = self.ZeroMap
        self.CumInfiltExcess = self.ZeroMap
        self.CumExfiltWater = self.ZeroMap
        self.CumSurfaceWater = self.ZeroMap
        self.CumEvap = self.ZeroMap
        self.CumPotenEvap = self.ZeroMap
        self.CumInt = self.ZeroMap
        self.CumRad = self.ZeroMap
        self.CumLeakage = self.ZeroMap
        self.CumPrecPol = self.ZeroMap
        self.FirstZoneFlux = self.ZeroMap
        self.FreeWaterDepth = self.ZeroMap
        self.SumCellWatBal = self.ZeroMap
        self.PathInfiltExceeded = self.ZeroMap
        self.SoilInfiltExceeded = self.ZeroMap
        self.CumOutFlow = self.ZeroMap
        self.CumCellInFlow = self.ZeroMap
        self.CumIF = self.ZeroMap
        self.CumSeepage = self.ZeroMap
        self.CumActInfilt = self.ZeroMap
        self.Aspect = pcr.scalar(pcr.aspect(self.Altitude))  # aspect [deg]
        self.Aspect = pcr.ifthenelse(self.Aspect <= 0.0, pcr.scalar(0.001), self.Aspect)
        # On Flat areas the Aspect function fails, fill in with average...
        self.Aspect = pcr.ifthenelse(
            pcr.defined(self.Aspect),
            self.Aspect,
            pcr.areaaverage(self.Aspect, self.TopoId),
        )

        # Set DCL to riverlength if that is longer that the basic length calculated from grid
        drainlength = detdrainlength(self.TopoLdd, self.xl, self.yl)

        self.DCL = pcr.max(drainlength, self.RiverLength)  # m
        # Multiply with Factor (taken from upscaling operation, defaults to 1.0 if no map is supplied
        self.DCL = self.DCL * pcr.max(1.0, self.RiverLengthFac)

        # water depth (m)
        # set width for kinematic wave to cell width for all cells
        self.Bw = detdrainwidth(self.TopoLdd, self.xl, self.yl)
        # However, in the main river we have real flow so set the width to the
        # width of the river

        self.Bw = pcr.ifthenelse(self.River, RiverWidth, self.Bw)

        # Add rivers to the WaterFrac, but check with waterfrac map
        self.RiverFrac = pcr.min(
            1.0,
            pcr.ifthenelse(
                self.River, (RiverWidth * self.DCL) / (self.xl * self.yl), 0
            ),
        )
        self.WaterFrac = self.WaterFrac - pcr.ifthenelse(
            (self.RiverFrac + self.WaterFrac) > 1.0,
            self.RiverFrac + self.WaterFrac - 1.0,
            0.0,
        )

        # term for Alpha
        self.AlpTerm = pow((self.N / (pcr.sqrt(self.Slope))), self.Beta)
        # power for Alpha
        self.AlpPow = (2.0 / 3.0) * self.Beta
        # initial approximation for Alpha

        # self.initstorage=pcr.areaaverage(self.FirstZoneDepth,self.TopoId)+areaaverage(self.UStoreDepth,self.TopoId)#+areaaverage(self.Snow,self.TopoId)
        # calculate catchmentsize
        self.upsize = pcr.catchmenttotal(self.xl * self.yl, self.TopoLdd)
        self.csize = pcr.areamaximum(self.upsize, self.TopoId)
        # Save some summary maps
        self.logger.info("Saving summary maps...")
        if self.modelSnow:
            pcr.report(self.Cfmax, self.Dir + "/" + self.runId + "/outsum/Cfmax.map")
            pcr.report(self.TTI, self.Dir + "/" + self.runId + "/outsum/TTI.map")
            pcr.report(self.TT, self.Dir + "/" + self.runId + "/outsum/TT.map")
            pcr.report(self.WHC, self.Dir + "/" + self.runId + "/outsum/WHC.map")

        pcr.report(self.Cmax, self.Dir + "/" + self.runId + "/outsum/Cmax.map")
        pcr.report(
            self.csize, self.Dir + "/" + self.runId + "/outsum/CatchmentSize.map"
        )
        pcr.report(
            self.upsize, self.Dir + "/" + self.runId + "/outsum/UpstreamSize.map"
        )
        pcr.report(self.EoverR, self.Dir + "/" + self.runId + "/outsum/EoverR.map")
        pcr.report(
            self.RootingDepth, self.Dir + "/" + self.runId + "/outsum/RootingDepth.map"
        )
        pcr.report(
            self.CanopyGapFraction,
            self.Dir + "/" + self.runId + "/outsum/CanopyGapFraction.map",
        )
        pcr.report(
            self.InfiltCapSoil,
            self.Dir + "/" + self.runId + "/outsum/InfiltCapSoil.map",
        )
        pcr.report(
            self.InfiltCapPath,
            self.Dir + "/" + self.runId + "/outsum/InfiltCapPath.map",
        )
        pcr.report(self.PathFrac, self.Dir + "/" + self.runId + "/outsum/PathFrac.map")
        pcr.report(self.thetaR, self.Dir + "/" + self.runId + "/outsum/thetaR.map")
        pcr.report(self.thetaS, self.Dir + "/" + self.runId + "/outsum/thetaS.map")
        pcr.report(
            self.FirstZoneMinCapacity,
            self.Dir + "/" + self.runId + "/outsum/FirstZoneMinCapacity.map",
        )
        pcr.report(
            self.FirstZoneKsatVer,
            self.Dir + "/" + self.runId + "/outsum/FirstZoneKsatVer.map",
        )
        pcr.report(self.M, self.Dir + "/" + self.runId + "/outsum/M.map")
        pcr.report(
            self.FirstZoneCapacity,
            self.Dir + "/" + self.runId + "/outsum/FirstZoneCapacity.map",
        )
        pcr.report(Terrain_angle, self.Dir + "/" + self.runId + "/outsum/angle.map")
        pcr.report(self.Slope, self.Dir + "/" + self.runId + "/outsum/slope.map")
        pcr.report(WI, self.Dir + "/" + self.runId + "/outsum/WI.map")
        pcr.report(self.CC, self.Dir + "/" + self.runId + "/outsum/CC.map")
        pcr.report(self.N, self.Dir + "/" + self.runId + "/outsum/N.map")
        pcr.report(
            self.RiverFrac, self.Dir + "/" + self.runId + "/outsum/RiverFrac.map"
        )

        pcr.report(self.xl, self.Dir + "/" + self.runId + "/outsum/xl.map")
        pcr.report(self.yl, self.Dir + "/" + self.runId + "/outsum/yl.map")
        pcr.report(self.reallength, self.Dir + "/" + self.runId + "/outsum/rl.map")
        pcr.report(self.DCL, self.Dir + "/" + self.runId + "/outsum/DCL.map")
        pcr.report(self.Bw, self.Dir + "/" + self.runId + "/outsum/Bw.map")
        pcr.report(
            pcr.ifthen(self.River, self.Bw),
            self.Dir + "/" + self.runId + "/outsum/RiverWidth.map",
        )
        if self.updating:
            pcr.report(
                self.DistToUpdPt,
                self.Dir + "/" + self.runId + "/outsum/DistToUpdPt.map",
            )

        self.SaveDir = self.Dir + "/" + self.runId + "/"
        self.logger.info("Starting Dynamic run...")

    def resume(self):

        if self.reinit == 1:
            self.logger.info("Setting initial conditions to default")
            self.FirstZoneDepth = self.FirstZoneCapacity * 0.85
            self.UStoreDepth = self.FirstZoneCapacity * 0.0
            self.WaterLevel = self.ZeroMap
            self.SurfaceRunoff = self.ZeroMap
            self.Snow = self.ZeroMap
            self.SnowWater = self.ZeroMap
            self.TSoil = self.ZeroMap + 10.0
            self.CanopyStorage = self.ZeroMap

        else:
            self.logger.info("Setting initial conditions from state files")
            self.wf_resume(self.Dir + "/instate/")

        P = self.Bw + (2.0 * self.WaterLevel)
        self.Alpha = self.AlpTerm * pow(P, self.AlpPow)
        self.OldSurfaceRunoff = self.SurfaceRunoff

        self.SurfaceRunoffMM = self.SurfaceRunoff * self.QMMConv
        # Determine initial kinematic wave volume
        self.KinWaveVolume = self.WaterLevel * self.Bw * self.DCL
        self.OldKinWaveVolume = self.KinWaveVolume

        self.SurfaceRunoffMM = self.SurfaceRunoff * self.QMMConv
        self.InitialStorage = self.FirstZoneDepth + self.UStoreDepth
        self.CellStorage = self.FirstZoneDepth + self.UStoreDepth

        # Determine actual water depth
        self.zi = pcr.max(
            0.0,
            self.FirstZoneThickness - self.FirstZoneDepth / (self.thetaS - self.thetaR),
        )
        # TOPOG_cqf type soil stuff
        self.f = (self.thetaS - self.thetaR) / self.M

    def dynamic(self):
        """
    Stuf that is done for each timestep
    
    
    Below a list of variables that can be save to disk as maps or as 
    timeseries (see ini file for syntax):
        
    :var self.SurfaceRunoff: Surface runoff in the kinematic wave [m^3/s]
    :var self.ActEvap: Actual EvapoTranspiration [mm]
    :var self.WaterLevel: Water level in the kinematic wave [m] (above the bottom)
    :var self.ActInfilt: Actual infiltration into the unsaturated zone [mm]
    :var self.CanopyStorage: actual canopystorage (only for subdaily timesteps) [mm]
    :var self.FirstZoneDepth: Amount of water in the saturated store [mm]
    :var self.UStoreDepth: Amount of water in the unsaturated store [mm]
    :var self.TSoil: Top soil temperature [oC]
    :var self.FirstZoneDepth: amount of available water in the saturated part of the soil [mm]
    :var self.UStoreDepth: amount of available water in the unsaturated zone [mm]
    :var self.Transfer: downward flux from unsaturated to saturated zone [mm]
    :var self.CapFlux: capilary flux from saturated to unsaturated zone [mm]
    
    
    Static variables:
        
    :var self.Altitude: The altitude of each cell [m]
    :var self.Bw: Width of the river [m]
    :var self.River: booolean map indicating the presence of a river [-]
    :var self.DLC: length of the river within a cell [m]
    :var self.ToCubic: Mutiplier to convert mm to m^3/s for fluxes
    """

        self.logger.debug(
            "Step: "
            + str(int(self.thestep + self._d_firstTimeStep))
            + "/"
            + str(int(self._d_nrTimeSteps))
        )
        self.thestep = self.thestep + 1

        self.Precipitation = pcr.cover(self.wf_readmap(self.P_mapstack, 0.0), 0)
        self.HP = pcr.cover(self.wf_readmap(self.HP_mapstack, 0.0), 0)
        # self.PotenEvap=pcr.cover(self.wf_readmap(self.PET_mapstack,0.0),0)
        self.Radiation = pcr.cover(self.wf_readmap(self.RAD_mapstack, 0.0), 0)
        # Inflow=pcr.cover(self.readmap(self.Inflow),0)
        self.Inflow = pcrut.readmapSave(self.Inflow_mapstack, 0.0)
        self.Seepage = pcrut.readmapSave(self.Seepage_mapstack, 0.0)
        # Inflow=pcr.spatial(pcr.scalar(0.0))
        self.Temperature = pcr.cover(self.wf_readmap(self.TEMP_mapstack, 0.0), 0)
        self.RH = pcr.cover(self.wf_readmap(self.RH_mapstack, 0.0), 0)
        self.WindSpeed = pcr.cover(self.wf_readmap(self.WindSpeed_mapstack, 0.0), 0)
        self.WaterCatch = pcr.cover(self.wf_readmap(self.WaterCatch_mapstack, 0.0), 0)

        for k, v in multdynapars.items():
            estr = k + "=" + k + "*" + str(v)
            self.logger.debug("Dynamic Parameter multiplication: " + estr)
            exec(estr)

        # PotEvap = self.PotenEvap #

        self.CanopyStorage = self.CanopyStorage + self.WaterCatch

        ShortWave = self.Radiation  # Change this to measure later!!!!
        ##########################################################################
        # Calculate Penman Montieth evaporation ##################################
        ##########################################################################

        # Estimate Soil Heat Flux from radiation and leaf area index
        G = self.Radiation * (1 - self.Albedo) * 1.0 / (self.LeafAreaIndex + 12)

        # Determine Esat and Delta according to Calder 1990
        Esat = 6.1078 * pcr.exp(17.2694 * self.Temperature / (self.Temperature + 237.3))
        Delta = Esat * 17.2694 * 237.3 / sqr(self.Temperature + 237.3)

        # Determine Eact using relative humidity
        Eact = self.RH * Esat / 100

        # Determine specific heat of air
        Lambda = 4185.5 * (751.78 - (0.5655 * (self.Temperature + 273.15)))

        # Now determine Gamma
        p = 900.0  # pressure in mb
        cp = 1005.0  # J/(kgK)
        Gamma = (cp * p) / (0.622 * Lambda)  #  mbC
        # density of dry air in kg/m^3
        rho = 1.201 * (290 * (p - 0.378 * Eact) / (1000 * (self.Temperature + 273.15)))
        # rho = hml_rho(p,Eact,Temperature);

        # At present set A to Radiation * (1- Albedo) and split according to
        # wetted part of the Canopy
        WetPart = pcr.min(1.0, self.CanopyStorage / self.Cmax)
        Atrans = (self.Radiation - G) * (1 - self.Albedo) * (1 - WetPart)
        A = (self.Radiation - G) * (1 - self.Albedo)
        Acanopy = (self.Radiation - G) * (1 - self.Albedo)
        Aint = Acanopy - Atrans

        # Potential in mm temp var, needed for check
        AtransMM = Atrans / (Delta + Gamma) / Lambda
        AintMM = Aint / (Delta + Gamma) / Lambda

        # Determine Ra as a function of Windspeed and canopy parameters
        # Calculates ra (aerodynamic resistance) according to Arnouds function
        # CQ Specific function!
        z = self.WindSpeedHeigth
        Zom = self.VegetationHeigth * 0.123
        Zoh = 0.25 * Zom
        d = 0.66 * self.VegetationHeigth

        Ra = (
            4.72
            * pcr.ln((z - d) / Zom)
            * pcr.ln((z - d) / Zoh)
            / (1 + 0.54 * self.WindSpeed)
        )

        # Now the actual formula, this is for Interception, rs is zero
        VPD = Esat - Eact
        n = Delta + Gamma
        # for interception Rs = 0;
        Rs = self.ZeroMap
        tmp = Rs / Ra
        nn = Delta + (Gamma * (1 + tmp))

        t = (Delta * A) + (rho * cp * VPD / Ra)
        # t = (Delta * Aint) + (rho * cp * VPD / Ra);
        EA = t / n

        PotEvap = EA / Lambda * self.timestepsecs  # now in mm

        # Now the actual formula, this is for Transpiration
        # Determine Rs seperate for Pasture and Forest (Hard coded, should be paramiterized in files later)
        # Reference equations
        InVPD = pcr.ifthenelse(VPD < 0.01, 0.01, VPD)
        InWave = pcr.ifthenelse(ShortWave < 5.0, 5.0, ShortWave)
        # CQ Specific function!
        PasRs = pcr.exp(1.05 * pcr.ln(InVPD) - 0.651 * pcr.ln(InWave) + 5.89)
        # CQ Specific function!
        ForRs = pcr.exp(0.867 * pcr.ln(InVPD) - 0.000831 * InWave + 2.81)

        Rs = pcr.max(
            0.5,
            pcr.min(1000, pcr.ifthenelse(pcr.scalar(self.LandUse) > 1.0, ForRs, PasRs)),
        )

        # No transpiration at nigth, this is of no use to the trees.
        Rs = pcr.ifthenelse(InWave < 10.0, 500, Rs)

        ##########################################################################
        # Interception according to a modified Rutter model with hourly timesteps#
        ##########################################################################

        p = self.CanopyGapFraction
        pt = 0.1 * p

        # Amount of P that falls on the canopy
        Pfrac = (1 - p - pt) * self.Precipitation

        # S cannot be larger than Cmax, no gravity drainage bolow that
        DD = pcr.ifthenelse(
            self.CanopyStorage > self.Cmax, self.Cmax - self.CanopyStorage, 0.0
        )
        self.CanopyStorage = self.CanopyStorage - DD

        # Add the precipitation that falls on the canopy to the store
        self.CanopyStorage = self.CanopyStorage + Pfrac

        # Now do the Evap, make sure the store does not get negative
        dC = -1 * pcr.min(self.CanopyStorage, PotEvap)
        self.CanopyStorage = self.CanopyStorage + dC

        LeftOver = PotEvap + dC
        # Amount of evap not used

        # Now drain the canopy storage again if needed...
        D = pcr.ifthenelse(
            self.CanopyStorage > self.Cmax, self.CanopyStorage - self.Cmax, 0.0
        )
        self.CanopyStorage = self.CanopyStorage - D

        # Calculate throughfall
        ThroughFall = DD + D + p * self.Precipitation
        StemFlow = self.Precipitation * pt

        # Calculate interception, this is NET Interception
        NetInterception = self.Precipitation - ThroughFall - StemFlow
        Interception = -dC

        # Determine Evnergy left over for transpiration
        Atrans = pcr.ifthenelse(
            self.CanopyStorage > 0.001,
            pcr.ifthenelse(PotEvap > 0, LeftOver / PotEvap * A, A),
            A,
        )

        t = (Delta * Atrans) + (rho * cp * VPD / Ra)
        # t = (Delta * A) + (rho * cp * VPD/ Ra);
        tmp = Rs / Ra
        n = Delta + (Gamma * (1 + tmp))

        EA = t / n

        PotTrans = EA / Lambda * self.timestepsecs  # now in mm

        RestPotEvap = PotTrans

        # TODOL bring timeseries export also to the framework
        # sample timeseries
        # Do runoff always
        # self.runTss.sample(Runoff)
        # self.levTss.sample(self.WaterLevel)

        ##########################################################################
        # Start with the soil calculations  ######################################
        ##########################################################################

        self.ExfiltWater = self.ZeroMap
        FreeWaterDepth = self.ZeroMap

        ##########################################################################
        # Determine infiltration into Unsaturated store...########################
        ##########################################################################
        # Add precipitation surplus  FreeWater storage...
        FreeWaterDepth = ThroughFall + StemFlow
        UStoreCapacity = self.FirstZoneCapacity - self.FirstZoneDepth - self.UStoreDepth

        # Runoff onto water boddies and river network
        self.RunoffOpenWater = self.RiverFrac * self.WaterFrac * FreeWaterDepth
        # self.RunoffOpenWater = self.ZeroMap
        FreeWaterDepth = FreeWaterDepth - self.RunoffOpenWater

        if self.RunoffGenSigmaFunction:
            self.AbsoluteGW = self.DemMax - (self.zi * self.GWScale)
            self.SubCellFrac = sCurve(self.AbsoluteGW, c=self.CC, a=self.Altitude + 1.0)
            self.SubCellRunoff = self.SubCellFrac * FreeWaterDepth
            self.SubCellGWRunoff = pcr.min(
                self.SubCellFrac * self.FirstZoneDepth,
                self.SubCellFrac
                * self.Slope
                * self.FirstZoneKsatVer
                * pcr.exp(-self.f * self.zi),
            )
            self.FirstZoneDepth = self.FirstZoneDepth - self.SubCellGWRunoff
            FreeWaterDepth = FreeWaterDepth - self.SubCellRunoff
        else:
            self.AbsoluteGW = self.DemMax - (self.zi * self.GWScale)
            self.SubCellFrac = pcr.spatial(pcr.scalar(0.0))
            self.SubCellGWRunoff = pcr.spatial(pcr.scalar(0.0))
            self.SubCellRunoff = pcr.spatial(pcr.scalar(0.0))

        # ----->>
        # First determine if the soil infiltration capacity can deal with the
        # amount of water
        # split between infiltration in undisturbed soil and compacted areas (paths)

        SoilInf = FreeWaterDepth * (1 - self.PathFrac)
        PathInf = FreeWaterDepth * self.PathFrac
        if self.modelSnow:
            soilInfRedu = pcr.ifthenelse(self.TSoil < 0.0, self.cf_soil, 1.0)
        else:
            soilInfRedu = 1.0
        MaxInfiltSoil = pcr.min(self.InfiltCapSoil * soilInfRedu, SoilInf)

        self.SoilInfiltExceeded = self.SoilInfiltExceeded + pcr.scalar(
            self.InfiltCapSoil * soilInfRedu < SoilInf
        )
        InfiltSoil = pcr.min(MaxInfiltSoil, UStoreCapacity)
        self.UStoreDepth = self.UStoreDepth + InfiltSoil
        UStoreCapacity = UStoreCapacity - InfiltSoil
        FreeWaterDepth = FreeWaterDepth - InfiltSoil
        # <-------
        MaxInfiltPath = pcr.min(self.InfiltCapPath * soilInfRedu, PathInf)
        # self.PathInfiltExceeded=self.PathInfiltExceeded + pcr.ifthenelse(self.InfiltCapPath < FreeWaterDepth, pcr.scalar(1), pcr.scalar(0))
        self.PathInfiltExceeded = self.PathInfiltExceeded + pcr.scalar(
            self.InfiltCapPath * soilInfRedu < PathInf
        )
        InfiltPath = pcr.min(MaxInfiltPath, UStoreCapacity)
        self.UStoreDepth = self.UStoreDepth + InfiltPath
        UStoreCapacity = UStoreCapacity - InfiltPath
        FreeWaterDepth = FreeWaterDepth - InfiltPath

        self.ActInfilt = InfiltPath + InfiltSoil

        self.InfiltExcess = pcr.ifthenelse(UStoreCapacity > 0.0, FreeWaterDepth, 0.0)
        self.CumInfiltExcess = self.CumInfiltExcess + self.InfiltExcess

        self.ActEvap, self.FirstZoneDepth, self.UStoreDepth, self.ActEvapUStore = actEvap_SBM(
            self.RootingDepth,
            self.zi,
            self.UStoreDepth,
            self.FirstZoneDepth,
            PotTrans,
            self.rootdistpar,
        )
        # self.ActEvap = self.ZeroMap
        # self.ActEvapUStore = self.ZeroMap
        ##########################################################################
        # Transfer of water from unsaturated to saturated store...################
        ##########################################################################
        self.zi = pcr.max(
            0.0,
            self.FirstZoneThickness - self.FirstZoneDepth / (self.thetaS - self.thetaR),
        )  # Determine actual water depth
        Ksat = self.FirstZoneKsatVer * pcr.exp(-self.f * self.zi)
        self.DeepKsat = self.FirstZoneKsatVer * pcr.exp(
            -self.f * self.FirstZoneThickness
        )

        # Determine saturation deficit. NB, as noted by Vertessy and Elsenbeer 1997
        # this deficit does NOT take into account the water in the unsaturated zone
        SaturationDeficit = self.FirstZoneCapacity - self.FirstZoneDepth

        # now the actual tranfer to the saturated store..
        self.Transfer = pcr.min(
            self.UStoreDepth,
            pcr.ifthenelse(
                SaturationDeficit <= 0.00001,
                0.0,
                Ksat * self.UStoreDepth / (SaturationDeficit + 1),
            ),
        )
        # Determine Ksat at base
        # DeepTransfer = pcr.min(self.UStoreDepth,ifthenelse (SaturationDeficit <= 0.00001, 0.0, DeepKsat * self.UStoreDepth/(SaturationDeficit+1)))

        # Now add leakage
        # Limit to MaxLeakage/day. Leakage percentage gets bigger if the
        # storm is bigger (macropores start kicking in...
        ActLeakage = pcr.max(
            0,
            pcr.min(self.MaxLeakage, self.Transfer * pcr.exp(0.01 * self.Transfer) / e),
        )
        self.Transfer = self.Transfer - ActLeakage
        # Now add leakage. to deeper groundwater
        # ActLeakage = pcr.cover(pcr.max(0,pcr.min(self.MaxLeakage* timestepsecs/basetimestep,ActLeakage)),0)

        # Now look if there is Seeapage

        # ActLeakage = pcr.ifthenelse(self.Seepage > 0.0, -1.0 * Seepage, ActLeakage)
        self.FirstZoneDepth = self.FirstZoneDepth + self.Transfer - ActLeakage
        self.UStoreDepth = self.UStoreDepth - self.Transfer

        # Determine % saturated
        # Sat = pcr.ifthenelse(self.FirstZoneDepth >= (self.FirstZoneCapacity*0.999), pcr.scalar(1.0), pcr.scalar(0.0))
        self.Sat = pcr.max(
            self.SubCellFrac,
            pcr.scalar(self.FirstZoneDepth >= (self.FirstZoneCapacity * 0.999)),
        )
        # PercSat = pcr.areaaverage(pcr.scalar(Sat),self.TopoId) * 100

        ##########################################################################
        # Horizontal (downstream) transport of water #############################
        ##########################################################################

        if self.waterdem:
            waterDem = self.Altitude - (self.zi * 0.001)
            waterLdd = pcr.lddcreate(waterDem, 1e35, 1e35, 1e35, 1e35)
            # waterLdd = pcr.lddcreate(waterDem,1,1,1,1)
            waterSlope = pcr.max(
                0.00001, pcr.slope(waterDem) * pcr.celllength() / self.reallength
            )

        self.zi = pcr.max(
            0.0,
            self.FirstZoneThickness - self.FirstZoneDepth / (self.thetaS - self.thetaR),
        )  # Determine actual water depth

        if self.waterdem:
            MaxHor = pcr.max(
                0.0,
                pcr.min(
                    self.FirstZoneKsatVer
                    * waterSlope
                    * pcr.exp(-SaturationDeficit / self.M),
                    self.FirstZoneDepth,
                ),
            )
            self.FirstZoneFlux = accucapacityflux(waterLdd, self.FirstZoneDepth, MaxHor)
            self.FirstZoneDepth = accucapacitystate(
                waterLdd, self.FirstZoneDepth, MaxHor
            )
        else:
            #
            # MaxHor = pcr.max(0,pcr.min(self.FirstZoneKsatVer * self.Slope * pcr.exp(-SaturationDeficit/self.M),self.FirstZoneDepth*(self.thetaS-self.thetaR))) * timestepsecs/basetimestep
            MaxHor = pcr.max(
                0.0,
                pcr.min(
                    self.FirstZoneKsatVer
                    * self.Slope
                    * pcr.exp(-SaturationDeficit / self.M),
                    self.FirstZoneDepth,
                ),
            )
            self.FirstZoneFlux = accucapacityflux(
                self.TopoLdd, self.FirstZoneDepth, MaxHor
            )
            self.FirstZoneDepth = accucapacitystate(
                self.TopoLdd, self.FirstZoneDepth, MaxHor
            )

        ##########################################################################
        # Determine returnflow from first zone          ##########################
        ##########################################################################
        self.ExfiltWaterFrac = sCurve(
            self.FirstZoneDepth, a=self.FirstZoneCapacity, c=5.0
        )
        self.ExfiltWater = self.ExfiltWaterFrac * (
            self.FirstZoneDepth - self.FirstZoneCapacity
        )
        # self.ExfiltWater=ifthenelse (self.FirstZoneDepth - self.FirstZoneCapacity > 0 , self.FirstZoneDepth - self.FirstZoneCapacity , 0.0)
        self.FirstZoneDepth = self.FirstZoneDepth - self.ExfiltWater

        # Re-determine UStoreCapacity
        UStoreCapacity = self.FirstZoneCapacity - self.FirstZoneDepth - self.UStoreDepth
        # Determine capilary rise
        self.zi = pcr.max(
            0.0,
            self.FirstZoneThickness - self.FirstZoneDepth / (self.thetaS - self.thetaR),
        )  # Determine actual water depth
        Ksat = self.FirstZoneKsatVer * pcr.exp(-self.f * self.zi)

        MaxCapFlux = pcr.max(
            0.0, pcr.min(Ksat, self.ActEvapUStore, UStoreCapacity, self.FirstZoneDepth)
        )
        # No capilary flux is roots are in water, max flux if very near to water, lower flux if distance is large
        CapFluxScale = pcr.ifthenelse(
            self.zi > self.RootingDepth,
            self.CapScale / (self.CapScale + self.zi - self.RootingDepth),
            0.0,
        )
        self.CapFlux = MaxCapFlux * CapFluxScale

        self.UStoreDepth = self.UStoreDepth + self.CapFlux
        self.FirstZoneDepth = self.FirstZoneDepth - self.CapFlux

        # org SurfaceWater = self.SurfaceRunoff * self.DCL * self.QMMConv # SurfaceWater (mm) from SurfaceRunoff (m3/s)
        SurfaceWater = (
            self.SurfaceRunoff * self.QMMConv
        )  # SurfaceWater (mm) from SurfaceRunoff (m3/s)
        self.CumSurfaceWater = self.CumSurfaceWater + SurfaceWater

        # Estimate water that may re-infiltrate
        if self.reInfilt:
            Reinfilt = pcr.max(
                0, pcr.min(SurfaceWater, pcr.min(self.InfiltCapSoil, UStoreCapacity))
            )
            self.CumReinfilt = self.CumReinfilt + Reinfilt
            self.UStoreDepth = self.UStoreDepth + Reinfilt
        else:
            Reinfilt = self.ZeroMap

        self.InwaterMM = pcr.max(
            0.0,
            self.ExfiltWater
            + FreeWaterDepth
            + self.SubCellRunoff
            + self.SubCellGWRunoff
            + self.RunoffOpenWater
            - Reinfilt,
        )
        self.Inwater = self.InwaterMM * self.ToCubic  # m3/s

        self.ExfiltWaterCubic = self.ExfiltWater * self.ToCubic
        self.SubCellGWRunoffCubic = self.SubCellGWRunoff * self.ToCubic
        self.SubCellRunoffCubic = self.SubCellRunoff * self.ToCubic
        self.InfiltExcessCubic = self.InfiltExcess * self.ToCubic
        self.FreeWaterDepthCubic = FreeWaterDepth * self.ToCubic
        self.ReinfiltCubic = -1.0 * Reinfilt * self.ToCubic
        self.Inwater = self.Inwater + self.Inflow  # Add abstractions/inflows in m^3/sec

        ##########################################################################
        # Runoff calculation via Kinematic wave ##################################
        ##########################################################################
        # per distance along stream
        q = self.Inwater / self.DCL
        # discharge (m3/s)
        self.SurfaceRunoff = pcr.kinematic(
            self.TopoLdd,
            self.SurfaceRunoff,
            q,
            self.Alpha,
            self.Beta,
            self.Tslice,
            self.timestepsecs,
            self.DCL,
        )  # m3/s
        self.SurfaceRunoffMM = (
            self.SurfaceRunoff * self.QMMConv
        )  # SurfaceRunoffMM (mm) from SurfaceRunoff (m3/s)
        self.updateRunOff()
        self.InflowKinWaveCell = pcr.upstream(self.TopoLdd, self.SurfaceRunoff)
        self.MassBalKinWave = (
            (self.KinWaveVolume - self.OldKinWaveVolume) / self.timestepsecs
            + self.InflowKinWaveCell
            + self.Inwater
            - self.SurfaceRunoff
        )

        Runoff = self.SurfaceRunoff

        # Updating
        # --------
        # Assume a tss file with as many columns as outpulocs. Start updating for each non-missing value and start with the
        # first column (nr 1). Assumes that outputloc and columns match!

        if self.updating:
            QM = pcr.timeinputscalar(updateFile, self.UpdateMap) * self.QMMConv

            # Now update the state. Just add to the Ustore
            # self.UStoreDepth =  result
            # No determine multiplication ratio for each gauge influence area.
            # For missing gauges 1.0 is assumed (no change).
            # UpDiff = pcr.areamaximum(QM,  self.UpdateMap) - pcr.areamaximum(self.SurfaceRunoffMM, self.UpdateMap)
            UpRatio = pcr.areamaximum(QM, self.UpdateMap) / pcr.areamaximum(
                self.SurfaceRunoffMM, self.UpdateMap
            )

            UpRatio = pcr.cover(pcr.areaaverage(UpRatio, self.TopoId), 1.0)
            # Now split between Soil and Kyn  wave
            UpRatioKyn = pcr.min(
                MaxUpdMult, pcr.max(MinUpdMult, (UpRatio - 1.0) * UpFrac + 1.0)
            )
            UpRatioSoil = pcr.min(
                MaxUpdMult, pcr.max(MinUpdMult, (UpRatio - 1.0) * (1.0 - UpFrac) + 1.0)
            )

            # update/nudge self.UStoreDepth for the whole upstream area,
            # not sure how much this helps or worsens things
            if UpdSoil:
                toadd = pcr.min(
                    (self.UStoreDepth * UpRatioSoil) - self.UStoreDepth,
                    StorageDeficit * 0.95,
                )
                self.UStoreDepth = self.UStoreDepth + toadd

            # Update the kinematic wave reservoir up to a maximum upstream distance
            # TODO:  add (much smaller) downstream updating also?
            MM = (1.0 - UpRatioKyn) / self.UpdMaxDist
            UpRatioKyn = MM * self.DistToUpdPt + UpRatioKyn

            self.SurfaceRunoff = self.SurfaceRunoff * UpRatioKyn
            self.SurfaceRunoffMM = (
                self.SurfaceRunoff * self.QMMConv
            )  # SurfaceRunoffMM (mm) from SurfaceRunoff (m3/s)
            self.updateRunOff()

            Runoff = self.SurfaceRunoff

        ##########################################################################
        # water balance ###########################################
        ##########################################################################

        # Single cell based water budget
        CellStorage = self.UStoreDepth + self.FirstZoneDepth
        DeltaStorage = CellStorage - self.InitialStorage
        OutFlow = self.FirstZoneFlux
        CellInFlow = pcr.upstream(self.TopoLdd, pcr.scalar(self.FirstZoneFlux))
        # CellWatBal = ActInfilt - self.ActEvap - self.ExfiltWater - ActLeakage + Reinfilt + IF - OutFlow + (OldCellStorage - CellStorage)
        # SumCellWatBal = SumCellWatBal + CellWatBal;

        self.CumOutFlow = self.CumOutFlow + OutFlow
        self.CumActInfilt = self.CumActInfilt + self.ActInfilt
        self.CumCellInFlow = self.CumCellInFlow + CellInFlow
        self.CumPrec = self.CumPrec + self.Precipitation
        self.CumEvap = self.CumEvap + self.ActEvap
        self.CumPotenEvap = self.CumPotenEvap + PotTrans
        self.CumInt = self.CumInt + Interception
        self.CumLeakage = self.CumLeakage + ActLeakage
        self.CumInwaterMM = self.CumInwaterMM + self.InwaterMM
        self.CumExfiltWater = self.CumExfiltWater + self.ExfiltWater
        # Water budget
        # self.watbal = self.CumPrec- self.CumEvap - self.CumInt - self.CumInwaterMM - DeltaStorage  - self.CumOutFlow + self.CumIF
        # self.watbal = self.CumActInfilt  - self.CumEvap - self.CumExfiltWater - DeltaStorage - self.CumOutFlow + self.CumIF
        self.watbal = (
            self.CumPrec
            + self.CumCellInFlow
            - self.CumOutFlow
            - self.CumEvap
            - self.CumLeakage
            - self.CumInwaterMM
            - self.CumInt
            - DeltaStorage
            + self.CumReinfilt
        )


def main():

    """
    Perform command line execution of the model.
    """
    caseName = "default_cqf"
    global multpars
    runId = "run_default"
    configfile = "wflow_cqf.ini"
    _lastTimeStep = 0
    _firstTimeStep = 1
    runinfoFile = "runinfo.xml"
    timestepsecs = 86400
    wflow_cloneMap = "wflow_subcatch.map"
    NoOverWrite = 1

    ## Main model starts here
    ########################################################################
    try:
        opts, args = getopt.getopt(
            sys.argv[1:], "XF:L:hC:Ii:v:S:T:WNR:u:s:EP:p:Xx:U:fOc:")
    except getopt.error as msg:
        pcrut.usage(msg)

    for o, a in opts:
        if o == "-P":
            exec("multpars =" + a, globals(), globals())
        if o == "-p":
            exec("multdynapars =" + a)
            exec("multdynapars =" + a, globals(), globals())
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
        if o == "-h":
            usage()
        if o == "-f":
            NoOverWrite = 0

    if _lastTimeStep < _firstTimeStep:
        usage()

    myModel = WflowModel(wflow_cloneMap, caseName, runId, configfile)
    dynModelFw = wf_DynamicFramework(
        myModel, _lastTimeStep, firstTimestep=_firstTimeStep
    )
    dynModelFw.createRunId(NoOverWrite=NoOverWrite)

    for o, a in opts:
        if o == "-X":
            configset(myModel.config, "model", "OverWriteInit", "1", overwrite=True)
        if o == "-I":
            configset(myModel.config, "model", "reinit", "1", overwrite=True)
        if o == "-i":
            configset(myModel.config, "model", "intbl", a, overwrite=True)
        if o == "-s":
            configset(myModel.config, "model", "timestepsecs", a, overwrite=True)
        if o == "-x":
            configset(myModel.config, "model", "sCatch", a, overwrite=True)
        if o == "-c":
            configset(myModel.config, "model", "configfile", a, overwrite=True)
        if o == "-M":
            configset(myModel.config, "model", "MassWasting", "1", overwrite=True)
        if o == "-N":
            configset(myModel.config, "model", "nolateral", "1", overwrite=True)
        if o == "-Q":
            configset(myModel.config, "model", "ExternalQbase", "1", overwrite=True)
        if o == "-U":
            configset(myModel.config, "model", "updateFile", a, overwrite=True)
            configset(myModel.config, "model", "updating", "1", overwrite=True)
        if o == "-u":
            print(a)
            exec("updateCols =" + a)
        if o == "-E":
            configset(myModel.config, "model", "reInfilt", "1", overwrite=True)
        if o == "-R":
            runId = a
        if o == "-W":
            configset(myModel.config, "model", "waterdem", "1", overwrite=True)

    dynModelFw._runInitial()
    dynModelFw._runResume()
    dynModelFw._runDynamic(_firstTimeStep, _lastTimeStep)
    dynModelFw._runSuspend()

    fp = open(caseName + "/" + runId + "/runinfo/configofrun.ini", "wb")
    myModel.config.write(fp)


if __name__ == "__main__":
    main()
