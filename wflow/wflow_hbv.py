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

# TODO: split off routing

"""
Run the wflow_hbv hydrological model..

usage:
wflow_hbv::

      [-h][-v level][-L logfile][-C casename][-R runId]
      [-c configfile][-T timesteps][-s seconds][-W][-E][-N][-U discharge]
      [-P parameter multiplication][-X][-l loglevel]

    -f: Force overwrite of existing results

    -T: Set end time of the run: yyyy-mm-dd hh:mm:ss

    -S: Set start time of the run: yyyy-mm-dd hh:mm:ss

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
        number of gauges in the wflow_gauges.map file.

    -u: list of gauges/columns to use in update. Format:
        -u [1 , 4 ,13]
        The above example uses column 1, 4 and 13

    -P: set parameter change string (e.g: -P "self.FC = self.FC * 1.6") for non-dynamic variables

    -p: set parameter change string (e.g: -P "self.Precipitation = self.Precipitation * 1.11") for
        dynamic variables

    -l: loglevel (most be one of DEBUG, WARNING, ERROR)

    -X overwrites the initial values at the end of each timestep


"""

import os.path

import pcraster.framework
import pcraster as pcr
import numpy as np

from wflow.wf_DynamicFramework import *
from wflow.wflow_adapt import *
from wflow.wflow_funcs import *

wflow = "wflow_hbv"


#: columns used in updating
updateCols = []  #: columns used in updating
""" Column used in updating """


def usage(*args):
    """
    Print usage information

    -  *args: command line arguments given
    """
    sys.stdout = sys.stderr
    for msg in args:
        print(msg)
    print(__doc__)
    sys.exit(0)


class WflowModel(pcraster.framework.DynamicModel):

    """
  The user defined model class.

  """

    def __init__(self, cloneMap, Dir, RunDir, configfile):
        pcraster.framework.DynamicModel.__init__(self)
        self.caseName = os.path.abspath(Dir)
        self.clonemappath = os.path.join(os.path.abspath(Dir), "staticmaps", cloneMap)
        pcr.setclone(self.clonemappath)
        self.runId = RunDir
        self.Dir = os.path.abspath(Dir)
        self.configfile = configfile
        self.SaveDir = os.path.join(self.Dir, self.runId)

    def updateRunOff(self):
        """
      Updates the kinematic wave reservoir
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

     :var self.SurfaceRunoff: Surface runoff in the kin-wave resrvoir [m^3/s]
     :var self.WaterLevel: Water level in the kin-wave resrvoir [m]
     :var self.DrySnow: Snow pack [mm]
     :var self.FreeWater:  Available free water [mm]
     :var self.UpperZoneStorage: Water in the upper zone [mm]
     :var self.LowerZoneStorage: Water in the lower zone [mm]
     :var self.SoilMoisture: Soil moisture [mm]
     :var self.InterceptionStorage: Amount of water on the Canopy [mm]

      """
        states = [
            "FreeWater",
            "SoilMoisture",
            "UpperZoneStorage",
            "LowerZoneStorage",
            "InterceptionStorage",
            "SurfaceRunoff",
            "WaterLevel",
            "DrySnow",
        ]

        if hasattr(self, "ReserVoirSimpleLocs"):
            states.append("ReservoirVolume")

        if hasattr(self, "LakeLocs"):
            states.append("LakeWaterLevel")
            
        if hasattr(self, "GlacierFrac"):
            states.append("GlacierStore")

        return states

    # The following are made to better connect to deltashell/openmi
    def supplyCurrentTime(self):
        """
      gets the current time in seconds after the start of the run

      Ouput:
          - time in seconds since the start of the model run
      """
        return self.currentTimeStep() * int(
            configget(self.config, "run", "timestepsecs", "86400")
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

        # Meteo and other forcing

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
        self.Seepage_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Seepage", "/inmaps/SE"
        )  # timeseries for rainfall "/inmaps/SE" # in/outflow locations (abstractions
        # Meteo and other forcing
        modelparameters.append(
            self.ParamType(
                name="Precipitation",
                stack=self.P_mapstack,
                type="timeseries",
                default=0.0,
                verbose=True,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="PotEvaporation",
                stack=self.PET_mapstack,
                type="timeseries",
                default=0.0,
                verbose=True,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="Temperature",
                stack=self.TEMP_mapstack,
                type="timeseries",
                default=10.0,
                verbose=True,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="Inflow",
                stack=self.Inflow_mapstack,
                type="timeseries",
                default=0.0,
                verbose=False,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="Seepage",
                stack=self.Seepage_mapstack,
                type="timeseries",
                default=0.0,
                verbose=False,
                lookupmaps=[],
            )
        )

        return modelparameters

    def suspend(self):
        """
      Suspends the model to disk. All variables needed to restart the model
      are saved to disk as pcraster maps. Use resume() to re-read them
    """

        self.logger.info("Saving initial conditions...")
        self.wf_suspend(os.path.join(self.SaveDir, "outstate"))

        if self.OverWriteInit:
            self.logger.info("Saving initial conditions over start conditions...")
            self.wf_suspend(os.path.join(self.SaveDir, "instate"))


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
    :var SUZ.tbl: Level over which K0 is used (100.0)
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

        pcr.setglobaloption("unittrue")

        self.thestep = pcr.scalar(0)
        self.basetimestep = 86400
        
        self.mv = -999
        self.count = 0

        #: files to be used in case of timesries (scalar) input to the model

        #: name of the tss file with precipitation data ("../intss/P.tss")
        self.precipTss = "../intss/P.tss"
        self.evapTss = (
            "../intss/PET.tss"
        )  #: name of the tss file with potential evap data ("../intss/PET.tss")
        self.tempTss = (
            "../intss/T.tss"
        )  #: name of the tss file with temperature  data ("../intss/T.tss")
        self.inflowTss = (
            "../intss/Inflow.tss"
        )  #: NOT TESTED name of the tss file with inflow data ("../intss/Inflow.tss")
        self.SeepageTss = (
            "../intss/Seepage.tss"
        )  #: NOT TESTED name of the tss file with seepage data ("../intss/Seepage.tss")"

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
        
        self.kinwaveIters = int(configget(self.config, "model", "kinwaveIters", "0"))         
        self.kinwaveTstep = int(configget(self.config, "model", "kinwaveTstep", "0"))     
        if self.kinwaveIters == 1:
            self.logger.info(
                "Using sub timestep for kinematic wave (iterate)"
            )
            if self.kinwaveTstep > 0:
                self.logger.info(
                    "Using a fixed timestep (seconds) for kinematic wave flow: " + str(self.kinwaveTstep)
                )

        self.sCatch = int(configget(self.config, "model", "sCatch", "0"))
        self.intbl = configget(self.config, "model", "intbl", "intbl")
        self.P_style = int(configget(self.config, "model", "P_style", "1"))
        self.PET_style = int(configget(self.config, "model", "PET_style", "1"))
        self.TEMP_style = int(configget(self.config, "model", "TEMP_style", "1"))

        self.modelSnow = int(configget(self.config, "model", "ModelSnow", "1"))
        sizeinmetres = int(configget(self.config, "layout", "sizeinmetres", "0"))
        alf = float(configget(self.config, "model", "Alpha", "60"))
        Qmax = float(configget(self.config, "model", "AnnualDischarge", "300"))
        self.UpdMaxDist = float(configget(self.config, "model", "UpdMaxDist", "100"))
        self.MaxUpdMult = float(configget(self.config, "model", "MaxUpdMult", "1.3"))
        self.MinUpdMult = float(configget(self.config, "model", "MinUpdMult", "0.7"))
        self.UpFrac = float(configget(self.config, "model", "UpFrac", "0.8"))
        self.ExternalQbase = int(configget(self.config, "model", "ExternalQbase", "0"))
        self.SetKquickFlow = int(configget(self.config, "model", "SetKquickFlow", "0"))
        self.MassWasting = int(configget(self.config, "model", "MassWasting", "0"))
        self.SubCatchFlowOnly = int(
            configget(self.config, "model", "SubCatchFlowOnly", "0")
        )
        self.NRiverMethod = int(configget(self.config, "model", "nrivermethod", "1"))

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
        wflow_riverwidth = configget(
            self.config, "model", "wflow_riverwidth", "staticmaps/wflow_riverwidth.map"
        )
        wflow_streamorder = configget(
            self.config, "model", "wflow_streamorder", "staticmaps/wflow_streamorder.map"
        )

        # 2: Input base maps ########################################################
        subcatch = pcr.ordinal(
            self.wf_readmap(os.path.join(self.Dir, wflow_subcatch), 0.0, fail=True)
        )  # Determines the area of calculations (all cells > 0)
        subcatch = pcr.ifthen(subcatch > 0, subcatch)
        if self.sCatch > 0:
            subcatch = pcr.ifthen(subcatch == sCatch, subcatch)

        self.Altitude = self.wf_readmap(
            os.path.join(self.Dir, wflow_dem), 0.0, fail=True
        ) * pcr.scalar(
            pcr.defined(subcatch)
        )  #: The digital elevation map (DEM)
        self.TopoLdd = self.wf_readmap(
            os.path.join(self.Dir, wflow_ldd), 0.0, fail=True
        )  #: The local drinage definition map (ldd)
        self.TopoId = pcr.ordinal(
            self.wf_readmap(os.path.join(self.Dir, wflow_subcatch), 0.0, fail=True)
        )  #: Map define the area over which the calculations are done (mask)
        self.River = pcr.cover(
            pcr.boolean(
                self.wf_readmap(os.path.join(self.Dir, wflow_river), 0.0, fail=True)
            ),
            0,
        )  #: river network map. Fro those cell that belong to a river a specific width is used in the kinematic wave caulations
        self.RiverLength = self.wf_readmap(
            os.path.join(self.Dir, wflow_riverlength), 0.0
        )
        # Factor to multiply riverlength with (defaults to 1.0)
        self.RiverLengthFac = self.wf_readmap(
            os.path.join(self.Dir, wflow_riverlength_fact), 1.0
        )

        # read landuse and soilmap and make sure there are no missing points related to the
        # subcatchment map. Currently sets the lu and soil type  type to 1
        self.LandUse = self.wf_readmap(
            os.path.join(self.Dir, wflow_landuse), 0.0, fail=True
        )  #: Map with lan-use/cover classes
        self.LandUse = pcr.cover(self.LandUse, pcr.nominal(pcr.ordinal(subcatch) > 0))
        self.Soil = self.wf_readmap(
            os.path.join(self.Dir, wflow_soil), 0.0, fail=True
        )  #: Map with soil classes
        self.Soil = pcr.cover(self.Soil, pcr.nominal(pcr.ordinal(subcatch) > 0))
        self.OutputLoc = self.wf_readmap(
            os.path.join(self.Dir, wflow_gauges), 0.0, fail=True
        )  #: Map with locations of output gauge(s)
        self.InflowLoc = pcr.nominal(
            self.wf_readmap(os.path.join(self.Dir, wflow_inflow), 0.0)
        )  #: Map with location of abstractions/inflows.
        self.SeepageLoc = self.wf_readmap(
            os.path.join(self.Dir, wflow_inflow), 0.0
        )  #: Seapage from external model (if configured)
        RiverWidth = self.wf_readmap(os.path.join(self.Dir, wflow_riverwidth), 0.0)

        # Temperature correction per cell to add
        self.TempCor = self.wf_readmap(
            os.path.join(
                self.Dir,
                configget(
                    self.config,
                    "model",
                    "TemperatureCorrectionMap",
                    "staticmap/swflow_tempcor.map",
                ),
            ),
            0.0,
        )

        if self.scalarInput:
            self.gaugesMap = self.wf_readmap(
                os.path.join(self.Dir, wflow_mgauges), 0.0, fail=True
            )  #: Map with locations of rainfall/evap/temp gauge(s). Only needed if the input to the model is not in maps
        self.OutputId = self.wf_readmap(
            os.path.join(self.Dir, wflow_subcatch), 0.0, fail=True
        )  # location of subcatchment

        self.ZeroMap = 0.0 * pcr.scalar(
            pcr.defined(self.Altitude)
        )  # map with only zero's

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
        self.Seepage_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Seepage", "/inmaps/SE"
        )  # timeseries for rainfall "/inmaps/SE" # in/outflow locations (abstractions)
        # For in memory override:
        self.P = self.ZeroMap
        self.PET = self.ZeroMap
        self.TEMP = self.ZeroMap
        # Set static initial values here #########################################

        self.Latitude = pcr.ycoordinate(pcr.boolean(self.Altitude))
        self.Longitude = pcr.xcoordinate(pcr.boolean(self.Altitude))

        self.logger.info("Linking parameters to landuse, catchment and soil...")

        self.Beta = pcr.scalar(0.6)  # For sheetflow
        # self.M=pcr.lookupscalar(self.Dir + "/" + modelEnv['intbl'] + "/M.tbl" ,self.LandUse,subcatch,self.Soil) # Decay parameter in Topog_sbm
        self.N = pcr.lookupscalar(
            self.Dir + "/" + self.intbl + "/N.tbl", self.LandUse, subcatch, self.Soil
        )  # Manning overland flow
        """ *Parameter:* Manning's N for all non-river cells """
        if self.NRiverMethod == 1:
            self.NRiver = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/N_River.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                0.036,
            )  # Manning river
        if self.NRiverMethod == 2:
            self.NRiver = self.readtblFlexDefault(
                self.Dir + "/" + self.intbl + "/N_River.tbl", 0.036, wflow_streamorder
            ) #Read from streamorder instead of landuse, subcatch and soil
        """ Manning's N for all cells that are marked as a river """

        self.wf_updateparameters()

        self.ReserVoirLocs = self.ZeroMap

        if hasattr(self, "ReserVoirSimpleLocs"):
            # Check if we have simple and or complex reservoirs
            tt_simple = pcr.pcr2numpy(self.ReserVoirSimpleLocs, 0.0)
            self.nrresSimple = tt_simple.max()
            self.ReserVoirLocs = self.ReserVoirLocs + pcr.cover(
                pcr.scalar(self.ReserVoirSimpleLocs), 0.0
            )
        else:
            self.nrresSimple = 0

        if hasattr(self, "LakeLocs"):
            #add parameter for lake threshold estimation
            self.estimatelakethresh = int(configget(self.config, "model", "estimatelakethresh", "0"))
            self.LakeAreasMap = pcr.nominal(self.LakeAreasMap)
            self.LakeLocs = pcr.nominal(self.LakeLocs)
            tt_lake = pcr.pcr2numpy(self.LakeLocs, 0.0)
            #self.nrlake = tt_lake.max()
            self.nrlake = np.size(np.where(tt_lake > 0.0)[0])
            self.ReserVoirLocs = self.ReserVoirLocs + pcr.cover(
                pcr.scalar(self.LakeLocs), 0.0
            )
            lake_area = pcr.cover(pcr.scalar(self.LakeAreasMap), 0.0)
            self.filter_P_PET = pcr.ifthenelse(
                lake_area > 0, lake_area * 0.0, self.filter_P_PET
            )

            # read files
            self.sh = {}
            lake_ids = pcr.ifthen(self.LakeStorFunc == 2, self.LakeLocs)
            np_lake_ids = pcr.pcr2numpy(lake_ids, 0)
            np_lake_ids_u = np.unique(np_lake_ids[np.nonzero(np_lake_ids)])
            if np.size(np_lake_ids_u) > 0:
                for item in np.nditer(np_lake_ids_u):
                    self.sh[int(item)] = np.loadtxt(
                        self.Dir
                        + "/"
                        + self.intbl
                        + "/Lake_SH_"
                        + str(item)
                        + ".tbl"
                    )
            self.hq = {}
            lake_ids = pcr.ifthen(self.LakeOutflowFunc == 1, self.LakeLocs)
            np_lake_ids = pcr.pcr2numpy(lake_ids, 0)
            np_lake_ids_u = np.unique(np_lake_ids[np.nonzero(np_lake_ids)])
            if np.size(np_lake_ids_u) > 0:
                for item in np.nditer(np_lake_ids_u):
                    self.hq[int(item)] = np.loadtxt(
                        self.Dir
                        + "/"
                        + self.intbl
                        + "/Lake_HQ_"
                        + str(item)
                        + ".tbl",
                        skiprows=3,
                    )
                    
            #Ini for the Modified Puls Approach (Burek et al., 2013, LISFLOOD)
            #Check which lakes uses the puls approach (LakeOutflowFunc = 3)
            #And if the corresponding Lake_e=2 and LakeStorFunc=2
            
            #Update Lake_b in ini if ResThreshold different from zero
            
            np_lakeoutflowfunc_old = pcr.pcr2numpy(self.LakeOutflowFunc, 0)
            self.LakeOutflowFunc = pcr.ifthenelse(
                    pcr.pcrand(self.LakeOutflowFunc == 3, self.LakeStorFunc == 1),
                    2,
                    self.LakeOutflowFunc
                    )
            self.LakeOutflowFunc = pcr.ifthenelse(
                    pcr.pcrand(self.LakeOutflowFunc == 3, self.Lake_e == 2.0),
                    2,
                    self.LakeOutflowFunc
                    )
            np_lakeoutflowfunc = pcr.pcr2numpy(self.LakeOutflowFunc, 0)
            if np_lakeoutflowfunc_old.sum() != np_lakeoutflowfunc.sum():
                self.logger.warning("Lake outflow modelling using the modified puls approach selected "+ 
                                    "but found contradictory arguments for LakeStorFunc/Lake_e: "+
                                    "using the general iteration method instead")
            

        else:
            self.nrlake = 0

        if (self.nrresSimple + self.nrlake) > 0:
            self.ReserVoirLocs = pcr.ordinal(self.ReserVoirLocs)
            self.logger.info(
                "A total of "
                + str(self.nrresSimple)
                + " simple reservoirs and "
                + str(self.nrlake)
                + " lakes found."
            )
            self.ReserVoirDownstreamLocs = pcr.downstream(
                self.TopoLdd, self.ReserVoirLocs
            )
            self.TopoLddOrg = self.TopoLdd
            self.TopoLdd = pcr.lddrepair(
                pcr.cover(
                    pcr.ifthen(pcr.boolean(self.ReserVoirLocs), pcr.ldd(5)),
                    self.TopoLdd,
                )
            )

        # HBV Soil params
        self.FC = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/FC.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            260.0,
        )
        self.BetaSeepage = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/BetaSeepage.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            1.8,
        )  # exponent in soil runoff generation equation
        self.LP = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/LP.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.53000,
        )  # fraction of Fieldcapacity below which actual evaporation=potential evaporation (LP)
        self.K4 = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/K4.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.02307,
        )  # Recession constant baseflow   #K4=0.07; BASEFLOW:LINEARRESERVOIR
        if self.SetKquickFlow:
            self.KQuickFlow = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/KQuickFlow.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                0.09880,
            )  # recession rate at flow HQ     #KHQ=0.2; OUTFLOWUPPERZONE_NONLINEARRESERVOIR
            self.SUZ = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/SUZ.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                100.0,
            )  # Level over wich K0 is used
            self.K0 = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/K0.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                0.3,
            )  # K0
        else:
            self.KHQ = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/KHQ.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                0.09880,
            )  # recession rate at flow HQ     #KHQ=0.2; OUTFLOWUPPERZONE_NONLINEARRESERVOIR
            self.HQ = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/HQ.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                3.27000,
            )  # high flow rate HQ for which recession rate of upper reservoir is known   #HQ=3.76;
            self.AlphaNL = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/AlphaNL.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                1.1,
            )  # measure of non-linearity of upper reservoir  #Alpha=1.6;

        self.PERC = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/PERC.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.4000,
        )  # percolation from Upper to Lowerzone (mm/day)
        self.CFR = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/CFR.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.05000,
        )  # refreezing efficiency constant in refreezing of freewater in snow
        # self.FoCfmax=self.readtblDefault(self.Dir + "/" + modelEnv['intbl'] + "/FoCfmax.tbl",self.LandUse,subcatch,self.Soil, 0.6000)  # correcton factor for snow melt/refreezing in forested and non-forested areas
        self.Pcorr = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/Pcorr.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            1.0,
        )  # correction factor for precipitation
        self.RFCF = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/RFCF.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            1.0,
        )  # correction factor for rainfall
        self.SFCF = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/SFCF.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            1.0,
        )  # correction factor for snowfall
        self.Cflux = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/Cflux.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            2.0,
        )  # maximum capillary rise from runoff response routine to soil moisture routine
        self.ICF = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/ICF.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            2.0,
        )  # maximum interception storage (in forested AND non-forested areas)
        self.CEVPF = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/CEVPF.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            1.0,
        )  # correction factor for potential evaporation (1.15 in in forested areas )
        self.EPF = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/EPF.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.0,
        )  # exponent of correction factor for evaporation on days with precipitation
        self.ECORR = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/ECORR.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            1.0,
        )  # evap correction
        # Soil Moisture  parameters
        self.ECALT = self.ZeroMap + 0.00000  # evaporation lapse per 100m
        # self.Ecorr=self.ZeroMap+1            # correction factor for evaporation

        # HBV Snow parameters
        # critical temperature for snowmelt and refreezing:  TTI= 1.000
        self.TTI = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/TTI.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            1.0,
        )
        # TT = -1.41934 # defines interval in which precipitation falls as rainfall and snowfall
        self.TT = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/TT.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            -1.41934,
        )
        # Cfmax = 3.75653 # meltconstant in temperature-index
        self.Cfmax = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/Cfmax.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            3.75653,
        )
        # WHC= 0.10000        # fraction of Snowvolume that can store water
        self.WHC = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/WHC.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.1,
        )

        # Determine real slope and cell length
        self.xl, self.yl, self.reallength = pcrut.detRealCellLength(
            self.ZeroMap, sizeinmetres
        )
        self.Slope = pcr.slope(self.Altitude)
        self.Slope = pcr.ifthen(
            pcr.boolean(self.TopoId),
            pcr.max(0.001, self.Slope * pcr.celllength() / self.reallength),
        )
        Terrain_angle = pcr.scalar(pcr.atan(self.Slope))
        temp = (
            pcr.catchmenttotal(pcr.cover(1.0), self.TopoLdd)
            * self.reallength
            * 0.001
            * 0.001
            * self.reallength
        )
        self.QMMConvUp = pcr.cover(self.timestepsecs * 0.001) / temp

        # Multiply parameters with a factor (for calibration etc) -P option in command line

        self.wf_multparameters()
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
        # Use supplied riverwidth if possible, else calulate
        RiverWidth = pcr.ifthenelse(RiverWidth <= 0.0, W, RiverWidth)
        #Use W instead of RiverWidth for reservoirs and lake cells
        if self.nrresSimple > 0:
            self.RiverWidth = pcr.ifthenelse(
                    pcr.cover(pcr.scalar(self.ReservoirSimpleAreas), 0.0) > 0.0,
                    W,
                    self.RiverWidth
                    )
        if self.nrlake > 0:
            self.RiverWidth = pcr.ifthenelse(
                    pcr.cover(pcr.scalar(self.LakeAreasMap), 0.0) > 0.0,
                    W,
                    self.RiverWidth
                    )

        self.SnowWater = self.ZeroMap

        # Which columns/gauges to use/ignore in kinematic wave updating
        self.UpdateMap = self.ZeroMap

        if self.updating:
            _tmp = pcr.pcr2numpy(self.OutputLoc, 0.0)
            gaugear = _tmp
            touse = numpy.zeros(gaugear.shape, dtype="int")

            for thecol in updateCols:
                idx = (gaugear == thecol).nonzero()
                touse[idx] = thecol

            self.UpdateMap = pcr.numpy2pcr(pcr.Nominal, touse, 0.0)
            # Calculate distance to updating points (upstream) annd use to scale the correction
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
            # self.DistToUpdPt = ldddist(self.TopoLdd,pcr.boolean(pcr.cover(self.OutputId,0.0)),1)
            # * self.reallength/celllength()

        # Initializing of variables
        self.logger.info("Initializing of model variables..")
        self.TopoLdd = pcr.lddmask(self.TopoLdd, pcr.boolean(self.TopoId))
        catchmentcells = pcr.maptotal(pcr.scalar(self.TopoId))

        # Limit lateral flow per subcatchment (make pits at all subcatch boundaries)
        # This is very handy for Ribasim etc...
        if self.SubCatchFlowOnly > 0:
            self.logger.info("Creating subcatchment-only drainage network (ldd)")
            ds = pcr.downstream(self.TopoLdd, self.TopoId)
            usid = pcr.ifthenelse(ds != self.TopoId, self.TopoId, 0)
            self.TopoLdd = pcr.lddrepair(
                pcr.ifthenelse(pcr.boolean(usid), pcr.ldd(5), self.TopoLdd)
            )

        # Used to seperate output per LandUse/management classes
        # OutZones = self.LandUse
        # pcr.report(self.reallength,"rl.map")
        # pcr.report(catchmentcells,"kk.map")
        self.QMMConv = self.timestepsecs / (
            self.reallength * self.reallength * 0.001
        )  # m3/s --> mm
        self.ToCubic = (
            self.reallength * self.reallength * 0.001
        ) / self.timestepsecs  # m3/s
        self.sumprecip = self.ZeroMap  #: accumulated rainfall for water balance
        self.sumevap = self.ZeroMap  #: accumulated evaporation for water balance
        self.sumrunoff = (
            self.ZeroMap
        )  #: accumulated runoff for water balance (weigthted for upstream area)
        self.sumlevel = self.ZeroMap  #: accumulated level for water balance
        self.sumpotevap = self.ZeroMap  # accumulated runoff for water balance
        self.sumsoilevap = self.ZeroMap
        self.sumtemp = self.ZeroMap  # accumulated runoff for water balance
        self.ForecQ_qmec = (
            self.ZeroMap
        )  # Extra inflow to kinematic wave reservoir for forcing in m^/sec
        self.KinWaveVolume = self.ZeroMap
        self.OldKinWaveVolume = self.ZeroMap
        self.Qvolume = self.ZeroMap
        self.Q = self.ZeroMap
        self.suminflow = self.ZeroMap
        # cntd
        self.FieldCapacity = self.FC  #: total water holding capacity of the soil
        self.Treshold = (
            self.LP * self.FieldCapacity
        )  # Threshold soilwaterstorage above which AE=PE
        # CatSurface=pcr.maptotal(pcr.scalar(pcr.ifthen(pcr.scalar(self.TopoId)>scalar(0.0),pcr.scalar(1.0))))                   # catchment surface (in  km2)

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

        # term for Alpha
        self.AlpTerm = pow((self.N / (pcr.sqrt(self.Slope))), self.Beta)
        # power for Alpha
        self.AlpPow = (2.0 / 3.0) * self.Beta
        # initial approximation for Alpha
        
        #Estimate LakeThreshold depending on outlet characteristics
        if (self.nrlake > 0 and self.estimatelakethresh == 1):
            #initial waterLevel
            level_map_path = self.Dir + "/instate/WaterLevelR.map"
            if os.path.exists(level_map_path):
                level_map = pcr.readmap(level_map_path)
            else:
                level_map = self.ZeroMap
            alphaR = self.AlpTermR * (pcr.downstream(self.TopoLdd, (self.Bw+2*level_map))) ** self.AlpPow
            outletRivLevel = alphaR * self.LakeAvgOut ** self.Beta
            #Lake Threshold = Lake Level - 130% River Level
            #130% = Adjustment linked to possible uncertainties in River Level estimation
            self.LakeThreshold = pcr.ifthenelse(
                    self.LakeThreshold > 0.0,
                    self.LakeThreshold,
                    pcr.max(self.LakeAvgLevel - 1.3 * outletRivLevel, 0.0)
                    )
            
            #Reupdate rating curve coefficient
            self.Lake_b = pcr.ifthenelse(
                    self.LakeOutflowFunc == 3,
                    self.LakeAvgOut / (self.LakeAvgLevel - self.LakeThreshold) ** 2,
                    self.Lake_b
                    )

        # calculate catchmentsize
        self.upsize = pcr.catchmenttotal(self.xl * self.yl, self.TopoLdd)
        self.csize = pcr.areamaximum(self.upsize, self.TopoId)
        
        
        ##### Set variables and framework for kinematic wave with iterations #####
        # convert pcr objects to numpy for kinemativ wave surface water
        np_zeros = pcr.pcr2numpy(self.ZeroMap, self.mv).ravel()
        np_2d_zeros = pcr.pcr2numpy(self.ZeroMap, self.mv)
        self.shape = np_2d_zeros.shape
        
        static_dtype = np.dtype(
                [('River', np.float64),
                 ('Beta', np.float64),
                 ('DCL', np.float64),
                 ('Bw', np.float64),
                 ('AlpPow', np.float64),
                 ('AlpTerm', np.float64)
                 ])
                    
        self.static = np.zeros(np_zeros.size, dtype=static_dtype)   
        self.static['River'] = pcr.pcr2numpy(self.River, self.mv).ravel()
        self.static['Beta'] = pcr.pcr2numpy(self.Beta, self.mv).ravel()
        self.static['DCL'] = pcr.pcr2numpy(self.DCL, self.mv).ravel()
        self.static['Bw'] = pcr.pcr2numpy(self.Bw, self.mv).ravel()
        self.static['AlpPow'] = pcr.pcr2numpy(self.AlpPow, self.mv).ravel()
        self.static['AlpTerm'] = pcr.pcr2numpy(self.AlpTerm, self.mv).ravel()
        
        dyn_dtype = np.dtype(
                [('SurfaceRunoff', np.float64),
                 ('Alpha', np.float64)
                 ])        
        
        self.dyn = np.zeros(np_zeros.size, dtype=dyn_dtype)
        
        # determine flow network and upstream nodes
        self.np_ldd = pcr.pcr2numpy(self.TopoLdd, self.mv)
        # initialize us-ds network for all cells
        self.nodes, self.nodes_up = set_dd(self.np_ldd)

        self.logger.info("End of initial section.")

    def default_summarymaps(self):
        """
      Returns a list of default summary-maps at the end of a run.
      This is model specific. You can also add them to the [summary]section of the ini file but stuff
      you think is crucial to the model should be listed here

       Example:

      """
        lst = [
            "self.Cfmax",
            "self.csize",
            "self.upsize",
            "self.TTI",
            "self.TT",
            "self.WHC",
            "self.Slope",
            "self.N",
            "self.xl",
            "self.yl",
            "self.reallength",
            "self.DCL",
            "self.Bw",
        ]

        return lst

    def resume(self):
        """ read initial state maps (they are output of a previous call to suspend()) """

        if self.reinit == 1:
            self.logger.info("Setting initial conditions to default")
            self.FreeWater = pcr.cover(0.0)  #: Water on surface (state variable [mm])
            self.SoilMoisture = self.FC  #: Soil moisture (state variable [mm])
            self.UpperZoneStorage = (
                0.2 * self.FC
            )  #: Storage in Upper Zone (state variable [mm])
            self.LowerZoneStorage = 1.0 / (
                3.0 * self.K4
            )  #: Storage in Uppe Zone (state variable [mm])
            self.InterceptionStorage = pcr.cover(
                0.0
            )  #: Interception Storage (state variable [mm])
            self.SurfaceRunoff = pcr.cover(
                0.0
            )  #: Discharge in kinimatic wave (state variable [m^3/s])
            self.WaterLevel = pcr.cover(
                0.0
            )  #: Water level in kinimatic wave (state variable [m])
            self.DrySnow = pcr.cover(0.0)  #: Snow amount (state variable [mm])
            if hasattr(self, "ReserVoirSimpleLocs"):
                self.ReservoirVolume = self.ResMaxVolume * self.ResTargetFullFrac
            if hasattr(self, "LakeLocs"):
                self.LakeWaterLevel = self.LakeAvgLevel
            if hasattr(self, "GlacierFrac"):
                self.GlacierStore = self.wf_readmap(
                    os.path.join(self.Dir, "staticmaps", "wflow_glacierstore.map"),
                    55.0 * 1000,
                )
        else:
            self.logger.info("Setting initial conditions from state files")
            self.wf_resume(os.path.join(self.Dir, "instate"))

        P = self.Bw + (2.0 * self.WaterLevel)
        self.Alpha = self.AlpTerm * pow(P, self.AlpPow)

        self.OldSurfaceRunoff = self.SurfaceRunoff

        self.SurfaceRunoffMM = self.SurfaceRunoff * self.QMMConv
        # Determine initial kinematic wave volume
        self.KinWaveVolume = self.WaterLevel * self.Bw * self.DCL
        self.OldKinWaveVolume = self.KinWaveVolume
        self.initstorage = (
            self.FreeWater
            + self.DrySnow
            + self.SoilMoisture
            + self.UpperZoneStorage
            + self.LowerZoneStorage
            + self.InterceptionStorage
        )

        if not self.SetKquickFlow:
            self.KQuickFlow = (self.KHQ ** (1.0 + self.AlphaNL)) * (
                self.HQ ** -self.AlphaNL
            )  # recession rate of the upper reservoir, KHQ*UHQ=HQ=kquickflow*(UHQ**alpha)

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
    :var self.InUpperZone: water inflow into Upper zone [mm]
    :var self.HBVSeepage: recharge to Upper zone [mm]
    :var self.DirectRunoff: direct runoff to Upper Zone [mm]
    :var self.BaseFlow: Specific runoff (baseflow part) per cell [mm]
    :var self.Percolation: actual percolation to the lower zone [mm]
    :var self.SoilMoisture: actual soil moisture [mm]
    :var se lf.QuickFlow: specific runoff (quickflow part) [mm]
    :var self.RealQuickFlow: specific runoff (quickflow), If K upper zone is precalculated [mm]
    :var self.CapFlux: capilary rise [mm]
    :var self.SurfaceRunoffMM: SurfaceRunoff in mm
    :var self.KinWaveVolume: Volume in the kinematic wave reservoir
    :var self.SurfaceWaterSupply: the negative Inflow (water demand) that could be met from the surfacewater [m^3/s]


    *Static variables*

    :var self.Altitude: The altitude of each cell [m]
    :var self.Bw: Width of the river [m]
    :var self.River: booolean map indicating the presence of a river [-]
    :var self.DLC: length of the river within a cell [m]
    :var self.ToCubic: Mutiplier to convert mm to m^3/s for fluxes
    """

        self.wf_updateparameters()  # read forcing an dynamic parameters
        self.Precipitation = pcr.max(0.0, self.Precipitation) * self.Pcorr

        # self.Precipitation=pcr.cover(self.wf_readmap(self.P_mapstack,0.0),0.0) * self.Pcorr
        # self.PotEvaporation=pcr.cover(self.wf_readmap(self.PET_mapstack,0.0),0.0)
        # self.Inflow=pcr.cover(self.wf_readmap(self.Inflow_mapstack,0.0,verbose=False),0.0)
        # These ar ALWAYS 0 at present!!!
        # self.Inflow=pcrut.readmapSave(self.Inflow_mapstack,0.0)
        if self.ExternalQbase:
            self.Seepage = pcr.cover(self.wf_readmap(self.Seepage_mapstack, 0.0), 0.0)
        else:
            self.Seepage = pcr.cover(0.0)
        self.Temperature = pcr.cover(self.wf_readmap(self.TEMP_mapstack, 10.0), 10.0)
        self.Temperature = self.Temperature + self.TempCor

        # Multiply input parameters with a factor (for calibration etc) -p option in command line (no also in ini)

        self.wf_multparameters()

        RainFrac = pcr.ifthenelse(
            1.0 * self.TTI == 0.0,
            pcr.ifthenelse(
                self.Temperature <= self.TT, pcr.scalar(0.0), pcr.scalar(1.0)
            ),
            pcr.min(
                (self.Temperature - (self.TT - self.TTI / 2.0)) / self.TTI,
                pcr.scalar(1.0),
            ),
        )
        RainFrac = pcr.max(
            RainFrac, pcr.scalar(0.0)
        )  # fraction of precipitation which falls as rain
        SnowFrac = 1.0 - RainFrac  # fraction of self.Precipitation which falls as snow

        self.Precipitation = (
            self.SFCF * SnowFrac * self.Precipitation
            + self.RFCF * RainFrac * self.Precipitation
        )  # different correction for rainfall and snowfall

        # Water onto the canopy
        Interception = pcr.min(
            self.Precipitation, self.ICF - self.InterceptionStorage
        )  #: Interception in mm/timestep
        self.InterceptionStorage = (
            self.InterceptionStorage + Interception
        )  #: Current interception storage
        self.Precipitation = self.Precipitation - Interception

        self.PotEvaporation = (
            pcr.exp(-self.EPF * self.Precipitation) * self.ECORR * self.PotEvaporation
        )  # correction for potential evaporation on wet days
        self.PotEvaporation = self.CEVPF * self.PotEvaporation  # Correct per landuse

        self.IntEvap = pcr.min(
            self.InterceptionStorage, self.PotEvaporation
        )  #: Evaporation from interception storage
        self.InterceptionStorage = self.InterceptionStorage - self.IntEvap

        # I nthe origal HBV code
        RestEvap = pcr.max(0.0, self.PotEvaporation - self.IntEvap)

        if hasattr(self, "ReserVoirComplexLocs"):
            self.ReserVoirPotEvap = self.PotEvaporation
            self.ReserVoirPrecip = self.Precipitation

            self.PotEvaporation = self.filter_P_PET * self.PotEvaporation
            self.Precipitation = self.filter_P_PET * self.Precipitation

        SnowFall = SnowFrac * self.Precipitation  #: snowfall depth
        RainFall = RainFrac * self.Precipitation  #: rainfall depth
        PotSnowMelt = pcr.ifthenelse(
            self.Temperature > self.TT,
            self.Cfmax * (self.Temperature - self.TT),
            pcr.scalar(0.0),
        )  # Potential snow melt, based on temperature
        PotRefreezing = pcr.ifthenelse(
            self.Temperature < self.TT,
            self.Cfmax * self.CFR * (self.TT - self.Temperature),
            0.0,
        )  # Potential refreezing, based on temperature

        Refreezing = pcr.ifthenelse(
            self.Temperature < self.TT, pcr.min(PotRefreezing, self.FreeWater), 0.0
        )  # actual refreezing
        self.SnowMelt = pcr.min(PotSnowMelt, self.DrySnow)  # actual snow melt
        self.DrySnow = (
            self.DrySnow + SnowFall + Refreezing - self.SnowMelt
        )  # dry snow content
        self.FreeWater = self.FreeWater - Refreezing  # free water content in snow
        MaxFreeWater = self.DrySnow * self.WHC
        self.FreeWater = self.FreeWater + self.SnowMelt + RainFall
        InSoil = pcr.max(
            self.FreeWater - MaxFreeWater, 0.0
        )  # abundant water in snow pack which goes into soil
        self.FreeWater = self.FreeWater - InSoil
        RainAndSnowmelt = RainFall + self.SnowMelt

        self.SnowCover = pcr.ifthenelse(self.DrySnow > 0, pcr.scalar(1), pcr.scalar(0))
        self.NrCell = pcr.areatotal(self.SnowCover, self.TopoId)
        

        # first part of precipitation is intercepted
        # Interception=pcr.min(InSoil,self.ICF-self.InterceptionStorage)#: Interception in mm/timestep
        # self.InterceptionStorage=self.InterceptionStorage+Interception #: Current interception storage
        # NetInSoil=InSoil-Interception
        NetInSoil = InSoil

        self.SoilMoisture = self.SoilMoisture + NetInSoil
        DirectRunoff = pcr.max(
            self.SoilMoisture - self.FieldCapacity, 0.0
        )  # if soil is filled to capacity: abundant water runs of directly
        self.SoilMoisture = self.SoilMoisture - DirectRunoff
        NetInSoil = NetInSoil - DirectRunoff  # net water which infiltrates into soil

        MaxSnowPack = 10000.0
        if self.MassWasting:
            # Masswasting of snow
            # 5.67 = tan 80 graden
            SnowFluxFrac = pcr.min(0.5, self.Slope / 5.67) * pcr.min(
                1.0, self.DrySnow / MaxSnowPack
            )
            MaxFlux = SnowFluxFrac * self.DrySnow
            self.DrySnow = accucapacitystate(self.TopoLdd, self.DrySnow, MaxFlux)
            self.FreeWater = accucapacitystate(
                self.TopoLdd, self.FreeWater, SnowFluxFrac * self.FreeWater
            )
        else:
            SnowFluxFrac = self.ZeroMap
            MaxFlux = self.ZeroMap
            
        if hasattr(self, "GlacierFrac"):
            """
            Run Glacier module and add the snowpack on-top of it.
            Estimate the fraction of snow turned into ice (HBV-light).
            Estimate glacier melt.
            glacierHBV function in wflow_lib.py
            """

            self.DrySnow, self.Snow2Glacier, self.GlacierStore, self.GlacierMelt = glacierHBV(
                self.GlacierFrac,
                self.GlacierStore,
                self.DrySnow,
                self.Temperature,
                self.G_TT,
                self.G_Cfmax,
                self.G_SIfrac,
                self.timestepsecs,
                self.basetimestep
            )
            # Convert to mm per grid cell and add to snowmelt
            self.GlacierMelt = self.GlacierMelt * self.GlacierFrac
            self.FreeWater = (
                self.FreeWater + self.GlacierMelt
            )

        # IntEvap=pcr.min(self.InterceptionStorage,self.PotEvaporation)  #: Evaporation from interception storage
        # self.InterceptionStorage=self.InterceptionStorage-IntEvap

        # I nthe origal HBV code
        # RestEvap = pcr.max(0.0,self.PotEvaporation-IntEvap)

        self.SoilEvap = pcr.ifthenelse(
            self.SoilMoisture > self.Treshold,
            pcr.min(self.SoilMoisture, RestEvap),
            pcr.min(
                self.SoilMoisture,
                pcr.min(
                    RestEvap, self.PotEvaporation * (self.SoilMoisture / self.Treshold)
                ),
            ),
        )
        #: soil evapotranspiration
        self.SoilMoisture = (
            self.SoilMoisture - self.SoilEvap
        )  # evaporation from soil moisture storage

        self.ActEvap = (
            self.IntEvap + self.SoilEvap
        )  #: Sum of evaporation components (IntEvap+SoilEvap)
        self.HBVSeepage = (
            (pcr.min(self.SoilMoisture / self.FieldCapacity, 1)) ** self.BetaSeepage
        ) * NetInSoil  # runoff water from soil
        self.SoilMoisture = self.SoilMoisture - self.HBVSeepage

        Backtosoil = pcr.min(
            self.FieldCapacity - self.SoilMoisture, DirectRunoff
        )  # correction for extremely wet periods: soil is filled to capacity
        self.DirectRunoff = DirectRunoff - Backtosoil
        self.SoilMoisture = self.SoilMoisture + Backtosoil
        self.InUpperZone = (
            self.DirectRunoff + self.HBVSeepage
        )  # total water available for runoff

        # Steps is always 1 at the moment
        # calculations for Upper zone
        self.UpperZoneStorage = (
            self.UpperZoneStorage + self.InUpperZone
        )  # incoming water from soil
        self.Percolation = pcr.min(
            self.PERC, self.UpperZoneStorage - self.InUpperZone / 2
        )  # Percolation
        self.UpperZoneStorage = self.UpperZoneStorage - self.Percolation
        self.CapFlux = self.Cflux * (
            ((self.FieldCapacity - self.SoilMoisture) / self.FieldCapacity)
        )  #: Capillary flux flowing back to soil
        self.CapFlux = pcr.min(self.UpperZoneStorage, self.CapFlux)
        self.CapFlux = pcr.min(self.FieldCapacity - self.SoilMoisture, self.CapFlux)
        self.UpperZoneStorage = self.UpperZoneStorage - self.CapFlux
        self.SoilMoisture = self.SoilMoisture + self.CapFlux

        if not self.SetKquickFlow:
            self.QuickFlow = pcr.min(
                pcr.ifthenelse(
                    self.Percolation < self.PERC,
                    0,
                    self.KQuickFlow
                    * (
                        (
                            self.UpperZoneStorage
                            - pcr.min(self.InUpperZone / 2, self.UpperZoneStorage)
                        )
                        ** (1.0 + self.AlphaNL)
                    ),
                ),
                self.UpperZoneStorage,
            )
            self.UpperZoneStorage = pcr.max(
                pcr.ifthenelse(
                    self.Percolation < self.PERC,
                    self.UpperZoneStorage,
                    self.UpperZoneStorage - self.QuickFlow,
                ),
                0,
            )
            # QuickFlow_temp = pcr.max(0,self.KQuickFlow*(self.UpperZoneStorage**(1.0+self.AlphaNL)))
            # self.QuickFlow = pcr.min(QuickFlow_temp,self.UpperZoneStorage)
            self.RealQuickFlow = self.ZeroMap
        else:
            self.QuickFlow = self.KQuickFlow * self.UpperZoneStorage
            self.RealQuickFlow = pcr.max(
                0, self.K0 * (self.UpperZoneStorage - self.SUZ)
            )
            self.UpperZoneStorage = (
                self.UpperZoneStorage - self.QuickFlow - self.RealQuickFlow
            )
        """Quickflow volume in mm/timestep"""
        # self.UpperZoneStorage=self.UpperZoneStorage-self.QuickFlow-self.RealQuickFlow

        # calculations for Lower zone
        self.LowerZoneStorage = self.LowerZoneStorage + self.Percolation
        self.BaseFlow = pcr.min(
            self.LowerZoneStorage, self.K4 * self.LowerZoneStorage
        )  #: Baseflow in mm/timestep
        self.LowerZoneStorage = self.LowerZoneStorage - self.BaseFlow
        # Direct runoff generation
        if self.ExternalQbase:
            DirectRunoffStorage = self.QuickFlow + self.Seepage + self.RealQuickFlow
        else:
            DirectRunoffStorage = self.QuickFlow + self.BaseFlow + self.RealQuickFlow

        self.InSoil = InSoil
        self.RainAndSnowmelt = RainAndSnowmelt
        self.NetInSoil = NetInSoil
        self.InwaterMM = pcr.max(0.0, DirectRunoffStorage)
        self.Inwater = self.InwaterMM * self.ToCubic
        
        self.Inflow = pcr.cover(self.Inflow, self.ZeroMap)
        # only run the reservoir module if needed

        if self.nrresSimple > 0:
            self.ReservoirVolume, self.Outflow, self.ResPercFull, self.DemandRelease = simplereservoir(
                self.ReservoirVolume,
                self.SurfaceRunoff,
                self.ResMaxVolume,
                self.ResTargetFullFrac,
                self.ResMaxRelease,
                self.ResDemand,
                self.ResTargetMinFrac,
                self.ReserVoirSimpleLocs,
                timestepsecs=self.timestepsecs,
            )
            self.OutflowDwn = pcr.upstream(
                self.TopoLddOrg, pcr.cover(self.Outflow, pcr.scalar(0.0))
            )
            self.Inflow = self.OutflowDwn + self.Inflow
        # else:
        #    self.Inflow= pcr.cover(self.Inflow,self.ZeroMap)

        if self.nrlake > 0:
            self.LakeWaterLevel, self.LakeOutflow, self.LakePrecip, self.LakeEvap, self.LakeVolume = naturalLake(
                self.LakeWaterLevel,
                self.LakeLocs,
                self.LinkedLakeLocs,
                self.LakeArea,
                self.LakeThreshold,
                self.LakeStorFunc,
                self.LakeOutflowFunc,
                self.sh,
                self.hq,
                self.Lake_b,
                self.Lake_e,
                self.RiverRunoff + self.LandRunoff + self.SubsurfaceFlow/1000/1000/1000/self.timestepsecs,
                self.ReserVoirPrecip,
                self.ReserVoirPotEvap,
                self.LakeAreasMap,
                self.wf_supplyJulianDOY(),
                timestepsecs=self.timestepsecs,
            )
            self.OutflowDwn = pcr.upstream(
                self.TopoLddOrg, pcr.cover(self.LakeOutflow, pcr.scalar(0.0))
            )
            self.Inflow = self.OutflowDwn + self.Inflow

        self.QuickFlowCubic = (self.QuickFlow + self.RealQuickFlow) * self.ToCubic
        self.BaseFlowCubic = self.BaseFlow * self.ToCubic

        self.SurfaceWaterSupply = pcr.ifthenelse(
            self.Inflow < 0.0,
            pcr.max(-1.0 * self.Inwater, self.SurfaceRunoff),
            self.ZeroMap,
        )
        self.Inwater = self.Inwater + pcr.ifthenelse(
            self.SurfaceWaterSupply > 0, -1.0 * self.SurfaceWaterSupply, self.Inflow
        )

        ##########################################################################
        # Runoff calculation via Kinematic wave ##################################
        ##########################################################################
        # per distance along stream
        q = self.Inwater / self.DCL + self.ForecQ_qmec / self.DCL
        self.OldSurfaceRunoff = self.SurfaceRunoff
        
        #Kinematic wave runoff iterations
        it_kinR=1
        if self.kinwaveIters == 1:
            if self.kinwaveTstep == 0:
                it_kinR = estimate_iterations_kin_wave(self.RiverRunoff, self.Beta, self.AlphaR, self.timestepsecs, self.DCL, self.mv)
            else:
                it_kinR = int(np.ceil(self.timestepsecs/self.kinwaveTstep))
            
        #Convert from pcr to numpy for the kinematic wave function
        q_np =  pcr.pcr2numpy(q,self.mv).ravel()        
        SurfaceRunoff = pcr.pcr2numpy(self.SurfaceRunoff,self.mv).ravel()
        self.dyn['Alpha'] = pcr.pcr2numpy(self.Alpha, self.mv).ravel()
        
        #Run the kinematic wave
        acc_flow = kin_wave(
                self.nodes,
                self.nodes_up,
                SurfaceRunoff,
                q_np,
                self.dyn['Alpha'],
                self.static['Beta'],
                self.static['DCL'],
                self.static['River'],
                self.static['Bw'],
                self.static['AlpTerm'],
                self.static['AlpPow'],
                self.timestepsecs,
                it_kinR) # m3/s

#        self.SurfaceRunoff = pcr.kinematic(
#            self.TopoLdd,
#            self.SurfaceRunoff,
#            q,
#            self.Alpha,
#            self.Beta,
#            self.Tslice,
#            self.timestepsecs,
#            self.DCL,
#        )  # m3/s
#        
        Qsurface = acc_flow/self.timestepsecs
        self.SurfaceRunoff = pcr.numpy2pcr(pcr.Scalar, np.copy(Qsurface).reshape(self.shape),self.mv)
                
        self.SurfaceRunoffMM = (
            self.SurfaceRunoff * self.QMMConv
        )  # SurfaceRunoffMM (mm) from SurfaceRunoff (m3/s)

        self.updateRunOff()
        InflowKinWaveCell = pcr.upstream(self.TopoLdd, self.SurfaceRunoff)
        self.MassBalKinWave = (
            (self.KinWaveVolume - self.OldKinWaveVolume) / self.timestepsecs
            + InflowKinWaveCell
            + self.Inwater
            - self.SurfaceRunoff
        )
        Runoff = self.SurfaceRunoff

        # Updating
        # --------
        # Assume a tss file with as many columns as outpulocs. Start updating for each non-missing value and start with the
        # first column (nr 1). Assumes that outputloc and columns match!

        if self.updating:
            QM = pcr.timeinputscalar(self.updateFile, self.UpdateMap) * self.QMMConv

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
            self.UpRatioKyn = pcr.min(
                self.MaxUpdMult,
                pcr.max(self.MinUpdMult, (UpRatio - 1.0) * self.UpFrac + 1.0),
            )
            UpRatioSoil = pcr.min(
                self.MaxUpdMult,
                pcr.max(self.MinUpdMult, (UpRatio - 1.0) * (1.0 - self.UpFrac) + 1.0),
            )

            # update/nudge self.UStoreDepth for the whole upstream area,
            # not sure how much this helps or worsens things
            UpdSoil = True
            if UpdSoil:
                toadd = (self.UpperZoneStorage * UpRatioSoil) - self.UpperZoneStorage
                self.UpperZoneStorage = self.UpperZoneStorage + toadd

            # Update the kinematic wave reservoir up to a maximum upstream distance
            # TODO:  add (much smaller) downstream updating also?
            MM = (1.0 - self.UpRatioKyn) / self.UpdMaxDist
            self.UpRatioKyn = MM * self.DistToUpdPt + self.UpRatioKyn

            self.SurfaceRunoff = self.SurfaceRunoff * self.UpRatioKyn
            self.SurfaceRunoffMM = (
                self.SurfaceRunoff * self.QMMConv
            )  # SurfaceRunoffMM (mm) from SurfaceRunoff (m3/s)
            self.updateRunOff()

            Runoff = self.SurfaceRunoff

        self.QCatchmentMM = self.SurfaceRunoff * self.QMMConvUp
        # self.RunoffCoeff = self.QCatchmentMM/catchmenttotal(self.Precipitation, self.TopoLdd)/catchmenttotal(pcr.cover(1.0), self.TopoLdd)

        self.sumprecip = (
            self.sumprecip + self.Precipitation
        )  # accumulated rainfall for water balance
        self.sumevap = (
            self.sumevap + self.ActEvap
        )  # accumulated evaporation for water balance
        self.sumsoilevap = self.sumsoilevap + self.SoilEvap
        self.sumpotevap = self.sumpotevap + self.PotEvaporation
        self.sumtemp = self.sumtemp + self.Temperature
        self.sumrunoff = (
            self.sumrunoff + self.InwaterMM
        )  # accumulated Cell runoff for water balance
        self.sumlevel = self.sumlevel + self.WaterLevel
        self.suminflow = self.suminflow + self.Inflow
        self.storage = (
            self.FreeWater
            + self.DrySnow
            + self.SoilMoisture
            + self.UpperZoneStorage
            + self.LowerZoneStorage
        )
        # + self.InterceptionStorage
        self.watbal = (
            (self.initstorage - self.storage)
            + self.sumprecip
            - self.sumsoilevap
            - self.sumrunoff
        )


# The main function is used to run the program from the command line


def main(argv=None):
    """
    Perform command line execution of the model.
    """
    global multpars
    global updateCols
    caseName = "default_hbv"
    runId = "run_default"
    configfile = "wflow_hbv.ini"
    LogFileName = "wflow.log"
    _lastTimeStep = 0
    _firstTimeStep = 0
    runinfoFile = "runinfo.xml"
    timestepsecs = 86400
    wflow_cloneMap = "wflow_subcatch.map"
    _NoOverWrite = 1
    loglevel = logging.DEBUG

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return

    ## Main model starts here
    ########################################################################
    try:
        opts, args = getopt.getopt(argv, "c:QXS:hC:Ii:T:R:u:s:P:p:Xx:U:fl:L:")
    except getopt.error as msg:
        pcrut.usage(msg)

    for o, a in opts:
        if o == "-C":
            caseName = a
        if o == "-R":
            runId = a
        if o == "-L":
            LogFileName = a
        if o == "-l":
            exec("loglevel = logging." + a)
        if o == "-c":
            configfile = a
        if o == "-s":
            timestepsecs = int(a)
        if o == "-h":
            usage()
        if o == "-f":
            _NoOverWrite = 0


    starttime = dt.datetime(1990, 1, 1)

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
        myModel, _lastTimeStep, firstTimestep=_firstTimeStep, datetimestart=starttime
    )
    dynModelFw.createRunId(
        NoOverWrite=_NoOverWrite,
        logfname=LogFileName,
        level=loglevel,
        model="wflow_hbv",
        doSetupFramework=False,
    )

    for o, a in opts:
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
        if o == "-X":
            configset(myModel.config, "model", "OverWriteInit", "1", overwrite=True)
        if o == "-I":
            configset(myModel.config, "run", "reinit", "1", overwrite=True)
        if o == "-i":
            configset(myModel.config, "model", "intbl", a, overwrite=True)
        if o == "-s":
            configset(myModel.config, "run", "timestepsecs", a, overwrite=True)
        if o == "-x":
            configset(myModel.config, "model", "sCatch", a, overwrite=True)
        if o == "-c":
            configset(myModel.config, "model", "configfile", a, overwrite=True)
        if o == "-M":
            configset(myModel.config, "model", "MassWasting", "0", overwrite=True)
        if o == "-Q":
            configset(myModel.config, "model", "ExternalQbase", "1", overwrite=True)
        if o == "-U":
            configset(myModel.config, "model", "updateFile", a, overwrite=True)
            configset(myModel.config, "model", "updating", "1", overwrite=True)
        if o == "-u":
            exec("zz =" + a)
            updateCols = zz
        if o == "-T":
            configset(myModel.config, "run", "endtime", a, overwrite=True)
        if o == "-S":
            configset(myModel.config, "run", "starttime", a, overwrite=True)

    dynModelFw.setupFramework()
    dynModelFw.logger.info("Command line: " + str(argv))
    dynModelFw._runInitial()

    dynModelFw._runResume()
    # dynModelFw._runDynamic(0,0)
    dynModelFw._runDynamic(_firstTimeStep, _lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()

    os.chdir("../../")


if __name__ == "__main__":
    main()
