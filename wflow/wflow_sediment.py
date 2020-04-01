#!/usr/bin/python

# Wflow is Free software, see below:
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
#
#@author: H. Boisgontier / Deltares 2018-2020
#
"""
Definition of the wflow_sediment model
The model is based on EUROSEM and ANSWERS for soil loss.
Inland routing uses Govers equation for general sediment transport and Yalin equation for particle class differentiation.
Transport and erosion in rivers and reservoirs are based on SWAT.
Ref: Beasley DB and Huggins LF. 1981. ANSWERS Users Manual. EPA-905/9-82-001:USEPA. Chicago, Il.
     Morgan et al. 1998. The European Soil Erosion Model (EUROSEM): documentation and user guide. Silsoe College, Cranfield University.
     Hessel R and Jetten V. 2007. Suitability of transport equations in modelling soil erosion for a small Loess Plateau catchment. Engineering Geology.
     Neitsch et al. 2011. Soil and Water Assessment Tool - Theoretical Documentation. Texas A&M University System.

---------------------------------------

Usage:
wflow_sediment  -C case -R Runid -c inifile -s seconds -T last timestep -S Firststimestep

    -C: set the name  of the case (directory) to run
    
    -R: set the name runId within the current case
    
    -c name of the config file (in the case directory)
    
    -s: Set the model timesteps in seconds
    
    -T: Set end time of the run: yyyy-mm-dd hh:mm:ss

    -S: Set start time of the run: yyyy-mm-dd hh:mm:ss

"""

import numpy as np
import os, sys
import os.path
import shutil, glob
import getopt

from wflow.wf_DynamicFramework import *
from wflow.wflow_adapt import *
from wflow.wflow_funcs import *

import pcraster.framework
import pcraster as pcr

wflow = "wflow_sediment: "


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args:
        print(msg)
    print(__doc__)
    sys.exit(0)
    
def estimate_sedriv_tranport_iter(Q, deltaT, h, x, W, mv = -999):
    """
    Estimate the number of iterations needed for sediment tranport in the river.
    Q: River runoff [m3/s]
    deltaT: model timestepsecs [s]
    h: river water level [m]
    x: river length [m]
    W: river width [m]
    """
    
    minTstep = pcr.ifthenelse(Q>0, h*x*W/Q, deltaT)
    minTstep = pcr.pcr2numpy(minTstep, mv)
    
    it = np.ceil(deltaT / np.amin(minTstep))
    #Maximum number of iterations
    it = max(it, 1000)
    
    return it


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


    def parameters(self):
        """
        Define all model parameters here that the framework should handle for the model
        See wf_updateparameters and the parameters section of the ini file
        If you use this make sure to all wf_updateparameters at the start of the dynamic section
        and at the start/end of the initial section
        
        *Meteo*
        :var Precipitation: Gross precipitation per timestep [mm]
        
        *Forcing from wflow_sbm*
        :var Interception: Rainfall intercepted by vegetation canopy [mm]
        :var RiverRunoff: Surface runoff in the kinematic wave in the river [m3/s]
        :var LandRunoff: Surface runoff in the kinematic wave inland [m3/s]
        :var WaterLevelR: Water level in the kinematic wave in the river [m]
        :var WaterLevelL: Water level in the kinematic wave inland [m]
        
      """

        modelparameters = []

        # Input time series from ini file
        self.P_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Precipitation", "/inmaps/P"
        )  # precipitation
        self.Int_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Interception", "/inmaps/int"
        )  # rainfall interception
        self.RR_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "RiverRunoff", "/inmaps/runR"
        )  # river runoff
        self.LR_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "LandRunoff", "/inmaps/runL"
        )  # land runoff
        self.WLR_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "WaterLevelR", "/inmaps/levKinR"
        )  # river water level
        self.WLL_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "WaterLevelL", "/inmaps/levKinL"
        )  # land water level

        # Meteo and other forcing from wflow_sbm
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
                name="Interception",
                stack=self.Int_mapstack,
                type="timeseries",
                default=0.0,
                verbose=False,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="RiverRunoff",
                stack=self.RR_mapstack,
                type="timeseries",
                default=0.0,
                verbose=True,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="LandRunoff",
                stack=self.LR_mapstack,
                type="timeseries",
                default=0.0,
                verbose=True,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="WaterLevelR",
                stack=self.WLR_mapstack,
                type="timeseries",
                default=0.0,
                verbose=False,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="WaterLevelL",
                stack=self.WLL_mapstack,
                type="timeseries",
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
      
      :var SedLoad: Total sediment concentration/load in the river [ton/cell/timestep]
      :var OutSedLoad: Total sediment load transported out of river cells [ton/cell/timestep]
      :var RivStoreSed: Deposited sediment storage in the river [ton/cell/timestep]
      """

        if self.RunRiverModel == 0:
            states = []
        else:
            states = [
                "SedLoad",
                "ClayLoad",
                "SiltLoad",
                "SandLoad",
                "SaggLoad",
                "LaggLoad",
                "GravLoad",
                "OutSedLoad",
                "OutClayLoad",
                "OutSiltLoad",
                "OutSandLoad",
                "OutSaggLoad",
                "OutLaggLoad",
                "OutGravLoad",
                "RivStoreSed",
                "RivStoreClay",
                "RivStoreSilt",
                "RivStoreSand",
                "RivStoreSagg",
                "RivStoreLagg",
                "RivStoreGrav",
            ]
            
            if hasattr(self, "ReserVoirSimpleLocs"):
                states.append("ReservoirSedStore")

            if hasattr(self, "LakeLocs"):
                states.append("LakeSedStore")

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

        return self.currentTimeStep() * self.timestepsecs

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
            self.wf_suspend(os.path.join(self.Dir, "instate"))

    def initial(self):

        """
    *Required*
    
    Initial part of the model, executed only once. It reads all static model
    information (parameters) and sets-up the variables used in modelling.
    
    This function is required. The contents is free. However, in order to
    easily connect to other models it is advised to adhere to the directory
    structure used in the other models.
    
    *Soil*
    :var PercentClay: mass fraction of clay content in the soil [%]
    :var PercentSilt: mass fraction of silt content in the soil [%]
    :var PercentOC: soil organic carbon content [%]
    :var BulkDensity: soil bulk density [kg/m3]
    :var PathFrac: Fraction of compacted area per grid cell [-]
    
    *Vegetation*
    :var CanopyHeight: height of the vegetation [m]
    
    *Erosion*
    :var ErosSpl: Exponent reducing the impact of splash erosion with water level [-]
    :var ErosOv: Coefficient for ANSWERS overland flow erosion [-]
    :var UsleC: USLE crop management factor [-]
    
    *River*
    :var D50River: Median particle diameter of the river bed and bank [mm]
    :var CovRiver: Factor representing bank erodibility reduction due to its vegetation [-]
    
    """
        #: pcraster option to calculate with units or cells. Not really an issue
        #: in this model but always good to keep in mind.
        pcr.setglobaloption("unittrue")

        self.timestepsecs = int(
            configget(self.config, "model", "timestepsecs", "86400")
        )
        self.reinit = int(configget(self.config, "run", "reinit", "0"))
        self.OverWriteInit = int(configget(self.config, "model", "OverWriteInit", "0"))
        sizeinmetres = int(configget(self.config, "layout", "sizeinmetres", "0"))
        self.basetimestep = 86400
        self.mv = -999

        # Reads all parameter from disk
        self.wf_updateparameters()

        # Set and get defaults from ConfigFile here ###################################

        self.intbl = configget(self.config, "model", "intbl", "intbl")
        self.RunRiverModel = int(configget(self.config, "model", "runrivermodel", "1"))
        self.UsleKMethod = int(configget(self.config, "model", "uslekmethod", "2"))
        self.RainErodMethod = int(
            configget(self.config, "model", "rainerodmethod", "1")
        )
        self.LandTransportMethod = int(
            configget(self.config, "model", "landtransportmethod", "1")
        )
        self.RivTransportMethod = int(
            configget(self.config, "model", "rivtransportmethod", "1")
        )
        if self.RunRiverModel == 1:
            self.LandTransportMethod = 1
            
        self.transportIters = int(configget(self.config, "model", "transportIters", "0"))
        self.transportRiverTstep = int(configget(self.config, "model", "transportRiverTstep", "0"))
        self.transportLandTstep = int(configget(self.config, "model", "transportLandTstep", "0"))        
        if self.transportIters == 1:
            self.logger.info(
                "Using sub timestep for sediment transport (iterate)"
            )
            if self.transportRiverTstep > 0:
                self.logger.info(
                    "Using a fixed timestep (seconds) for sediment transport in the river: " + str(self.transportRiverTstep)
                )
            if self.transportLandTstep > 0:
                self.logger.info(
                    "Using a fixed timestep (seconds) for sediment transport in overland flow: " + str(self.transportLandTstep)
                ) 
        
        self.dmClay = float(configget(self.config, "model", "dmClay", "2.0"))
        self.dmSilt = float(configget(self.config, "model", "dmSilt", "10.0"))
        self.dmSand = float(configget(self.config, "model", "dmSand", "200.0"))
        self.dmSagg = float(configget(self.config, "model", "dmSagg", "30.0"))
        self.dmLagg = float(configget(self.config, "model", "dmLagg", "500.0"))
        self.dmGrav = float(configget(self.config, "model", "dmGrav", "2000.0"))
        self.rhoSed = float(configget(self.config, "model", "rhoSed", "2650.0"))
        alf = float(configget(self.config, "model", "Alpha", "60"))
        Qmax = float(configget(self.config, "model", "AnnualDischarge", "300"))
        
        # static maps to use (normally default)
        wflow_dem = configget(
            self.config, "model", "wflow_dem", "staticmaps/wflow_dem.map"
        )
        wflow_landuse = configget(
            self.config, "model", "wflow_landuse", "staticmaps/wflow_landuse.map"
        )
        wflow_soil = configget(
            self.config, "model", "wflow_soil", "staticmaps/wflow_soil.map"
        )
        wflow_subcatch = configget(
            self.config, "model", "wflow_subcatch", "staticmaps/wflow_subcatch.map"
        )
        wflow_ldd = configget(
            self.config, "model", "wflow_ldd", "staticmaps/wflow_ldd.map"
        )
        wflow_river = configget(
            self.config, "model", "wflow_river", "staticmaps/wflow_river.map"
        )
        wflow_riverwidth = configget(
            self.config, "model", "wflow_riverwidth", "staticmaps/wflow_riverwidth.map"
        )
        wflow_riverlength = configget(
            self.config, "model", "wflow_riverlength", "staticmaps/wflow_riverlength.map"
        )
        wflow_riverlength_fact = configget(
            self.config,
            "model",
            "wflow_riverlength_fact",
            "staticmaps/wflow_riverlength_fact.map",
        )
#        wflow_riverslope = configget(
#            self.config, "model", "wflow_riverslope", "staticmaps/RiverSlope.map"
#        )
        wflow_streamorder = configget(
            self.config,
            "model",
            "wflow_streamorder",
            "staticmaps/wflow_streamorder.map",
        )

        # Soil
        wflow_clay = configget(
            self.config, "model", "wflow_clay", "staticmaps/percent_clay.map"
        )
        wflow_silt = configget(
            self.config, "model", "wflow_silt", "staticmaps/percent_silt.map"
        )
        if self.UsleKMethod == 3:
            wflow_oc = configget(
                self.config, "model", "wflow_oc", "staticmaps/percent_oc.map"
            )


        # 2: Input base maps ########################################################
        # Subcatchment map
        subcatch = pcr.ordinal(
            self.wf_readmap(os.path.join(self.Dir, wflow_subcatch), 0.0, fail=True)
        )  # Determines the area of calculations (all cells > 0)
        subcatch = pcr.ifthen(subcatch > 0, subcatch)
        
        self.Altitude = self.wf_readmap(
            os.path.join(self.Dir, wflow_dem), 0.0, fail=True
        )
#        self.LandUse = self.wf_readmap(
#            os.path.join(self.Dir, wflow_landuse), 0.0, fail=True
#        )
#        self.Soil = self.wf_readmap(os.path.join(self.Dir, wflow_soil), 0.0, fail=False)
        self.LandUse = pcr.ordinal(
            self.wf_readmap(os.path.join(self.Dir, wflow_landuse), 0.0, fail=True)
        )
        self.LandUse = pcr.cover(self.LandUse, pcr.ordinal(subcatch > 0))
        self.Soil = pcr.ordinal(
            self.wf_readmap(os.path.join(self.Dir, wflow_soil), 0.0, fail=True)
        )
        self.Soil = pcr.cover(self.Soil, pcr.ordinal(subcatch > 0))
        self.TopoId = self.wf_readmap(
            os.path.join(self.Dir, wflow_subcatch), 0.0, fail=True
        )
        self.TopoLdd = self.wf_readmap(
            os.path.join(self.Dir, wflow_ldd), 0.0, fail=True
        )
        self.River = self.wf_readmap(
            os.path.join(self.Dir, wflow_river), 0.0, fail=True
        )
        self.RiverWidth = self.wf_readmap(
            os.path.join(self.Dir, wflow_riverwidth), 0.0, fail=True
        )
        self.RiverLength = self.wf_readmap(
            os.path.join(self.Dir, wflow_riverlength), 0.0, fail=True
        )
        # Factor to multiply riverlength with (defaults to 1.0)
        self.RiverLengthFac = self.wf_readmap(
            os.path.join(self.Dir, wflow_riverlength_fact), 1.0
        )
#        self.RiverSlope = self.wf_readmap(
#            os.path.join(self.Dir, wflow_riverslope), 0.0, fail=True
#        )
        self.streamorder = self.wf_readmap(
            os.path.join(self.Dir, wflow_streamorder), 0.0, fail=True
        )
        
        self.PercentClay = self.wf_readmap(
            os.path.join(self.Dir, wflow_clay), 0.1, fail=True
        )
        self.PercentSilt = self.wf_readmap(
            os.path.join(self.Dir, wflow_silt), 0.1, fail=True
        )
        if self.UsleKMethod == 3:
            self.PercentOC = self.wf_readmap(
                os.path.join(self.Dir, wflow_oc), 0.1, fail=False
            )

        
        # Set static initial values here #########################################
        
        #Detachability of the soil (k) [g/J]
        self.ErosK = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/ErosK.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.6,
        )
        
        # USLE C factor map based on land use (from Gericke 2015, Soil loss estimation and empirical relationships for sediment delivery ratios of European river catchments)
        self.UsleC = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/USLE_C.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                0.01,
            )
        
        # Soil impervious area
        self.PathFrac = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/PathFrac.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.01,
        )

        # Soil model parameters
        self.ErosOv = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/eros_ov.tbl",
            self.LandUse,
            subcatch,
            self.Soil,
            0.90,
        )
        self.NRiver = self.readtblFlexDefault(
                self.Dir + "/" + self.intbl + "/N_River.tbl", 0.036, wflow_streamorder
            )

        if self.RunRiverModel == 1:
            # River model parameters
            self.D50River = self.readtblFlexDefault(
                self.Dir + "/" + self.intbl + "/D50_River.tbl", 0.050, wflow_streamorder
            )
            self.D50River = pcr.cover(self.D50River, pcr.scalar(0.0))
            self.CovRiver = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/cov_River.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                1.0,
            )
            
            # River particle size distribution (estimated with SWAT method)
            d50riv = self.Dir + "/" + self.runId + "/outsum/D50River.map"
            pcr.report(self.D50River, d50riv)
            self.FracClayRiv = self.readtblFlexDefault(
                self.Dir + "/" + self.intbl + "/ClayF_River.tbl", 0.15, d50riv
            )
            self.FracSiltRiv = self.readtblFlexDefault(
                self.Dir + "/" + self.intbl + "/SiltF_River.tbl", 0.65, d50riv
            )
            self.FracSandRiv = self.readtblFlexDefault(
                self.Dir + "/" + self.intbl + "/SandF_River.tbl", 0.15, d50riv
            )
            self.FracGravRiv = self.readtblFlexDefault(
                self.Dir + "/" + self.intbl + "/GravelF_River.tbl", 0.05, d50riv
            )

        """ Determine global variables """
        # Map with zeros
        self.ZeroMap = 0.0 * pcr.scalar(subcatch)  # map with only zero's

        # Determine real slope, cell length and cell area
        self.xl, self.yl, self.reallength = pcrut.detRealCellLength(
            self.ZeroMap, sizeinmetres
        )
        
        self.Slope = pcr.slope(self.Altitude)
        self.Slope = pcr.max(0.00001, self.Slope * pcr.celllength() / self.reallength)

        self.cellareaKm = (self.reallength / 1000.0) ** 2
        self.UpArea = pcr.accuflux(self.TopoLdd, self.cellareaKm)

        # Sine of the slope
        self.sinSlope = pcr.sin(pcr.atan(self.Slope))

        """ Determine variables for the soil loss model """

        # Canopy gap fraction based on LAI or input table
        if hasattr(self, "LAI"):
            if not hasattr(self, "Kext"):
                self.logger.error(
                    "Kext (canopy extinction coefficient) not defined! Needed becausee LAI is defined."
                )
                self.logger.error("Please add it to the modelparameters section. e.g.:")
                self.logger.error(
                    "Kext=inmaps/clim/LCtoExtinctionCoefficient.tbl,tbl,0.5,1,inmaps/clim/LC.map"
                )
            self.CanopyGapFraction = pcr.exp(-self.Kext * self.LAI)
        else:
            self.CanopyGapFraction = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/CanopyGapFraction.tbl",
                self.LandUse,
                self.TopoId,
                self.Soil,
                0.1,
            )

        # Determine sand content
        self.PercentSand = 100 - (self.PercentClay + self.PercentSilt)

        # Compute USLE K factor
        if self.UsleKMethod == 1:
            self.UsleK = self.wf_readmap(
                os.path.join(self.Dir, "staticmaps/USLE_K.map"), 0.1, fail=True
            )
        if self.UsleKMethod == 2:
            # Calculate USLE K factor (from Renard et al. 1997, with the geometric mean particle diameter Dg)
            self.Dg = pcr.exp(
                0.01
                * (
                    self.PercentClay * pcr.ln(0.001)
                    + self.PercentSilt * pcr.ln(0.025)
                    + self.PercentSand * pcr.ln(0.999)
                )
            )  # [mm]
            self.UsleK = 0.0034 + 0.0405 * pcr.exp(
                -1 / 2 * ((pcr.log10(self.Dg) + 1.659) / 0.7101) ** 2
            )
            # Remove possible outliers
            self.UsleK = pcr.max(0.0, self.UsleK)

        if self.UsleKMethod == 3:
            # Calculate USLE K factor (from Williams and Renard 1983, EPIC: a new method for assessing erosion's effect on soil productivity)
            self.SN = 1 - (self.PercentSand / 100)
            self.UsleK = (
                (
                    0.2
                    + 0.3
                    * pcr.exp(-0.0256 * self.PercentSand * (1 - self.PercentSilt / 100))
                )
                * (self.PercentSilt / (pcr.max(0.01, self.PercentClay + self.PercentSilt)))
                ** 0.3
                * (
                    1
                    - (0.25 * self.PercentOC)
                    / (self.PercentOC + pcr.exp(3.72 - 2.95 * self.PercentOC))
                )
                * (1 - (0.7 * self.SN) / (self.SN + pcr.exp(-5.51 + 22.9 * self.SN)))
            )
            # Remove possible outliers
            self.UsleK = pcr.max(0.0, self.UsleK)

        # Parameters for either EUROSEM or ANSWERS rainfall erosion
        if self.RainErodMethod == 1:  # EUROSEM
            # Canopy height [m]
            self.CanopyHeight = self.wf_readmap(
                os.path.join(self.Dir, "staticmaps/canopy_height.map"), 1.0, fail=True
            )
            # Coefficient
            self.ErosSpl = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/eros_spl.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                2.0,
            )
        if self.RainErodMethod == 2:  # ANSWERS
            # Coefficient
            self.ErosSpl = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/eros_spl.tbl",
                self.LandUse,
                subcatch,
                self.Soil,
                0.108,
            )

        """ Variables for inland transport model """
        # If the river model is run, use Yalin's equation with particle differentiation.
        # If just the soil loss model is run, use Govers equation without particle differentiation

        if self.LandTransportMethod == 2 or self.LandTransportMethod == 3:
            # Calculation of D50 and fraction of fine and very fine sand (fvfs) from Fooladmand et al, 2006
            self.PercentSand999 = (
                self.PercentSand * (999 - 25) / (1000 - 25)
            )  #%sand with a mean radius of 999um instead of 1000
            self.vd50 = pcr.ln(
                (1 / ((self.PercentClay + self.PercentSilt) / 100) - 1)
                / (1 / (self.PercentClay / 100) - 1)
            )
            self.wd50 = pcr.ln(
                (
                    1
                    / (
                        (self.PercentClay + self.PercentSilt + self.PercentSand999)
                        / 100
                    )
                    - 1
                )
                / (1 / (self.PercentClay / 100) - 1)
            )
            ad50 = 1 / (-3.727699)  # 1 / ln((25-1)/(999-1))
            bd50 = ad50 * 3.17805  # ad50 / ln((25-1)/1)
            self.cd50 = ad50 * pcr.ln(self.vd50 / self.wd50)
            self.ud50 = (-self.vd50) ** (1 - bd50) / (-self.wd50) ** (-bd50)

            self.D50 = 1 + (
                -1 / self.ud50 * pcr.ln(1 / (1 / (self.PercentClay / 100) - 1))
            ) ** (
                1 / self.cd50
            )  # [um]
            self.D50 = self.D50 / 1000  # [mm]

            #    #Fraction of fine/very fine sand (d=125um) and coarse sand
            #    self.percent_fvfs = 100 * 1 / (1+(1/(self.percent_clay/100)-1)*exp(-self.ud50*(125-1)**self.cd50)) - self.percent_silt - self.percent_clay #[%]
            #    self.percent_coars = self.percent_sand - self.percent_fvfs

            # Calculate Govers transport capacity coefficients
            self.cGovers = ((self.D50 * 1000 + 5) / 0.32) ** (-0.6)
            self.nGovers = ((self.D50 * 1000 + 5) / 300) ** 0.25

        elif self.LandTransportMethod == 1:
            # Determine sediment size distribution, estimated from primary particle size distribution (Foster et al., 1980)
            self.FracClay = 0.20 * self.PercentClay / 100
            self.FracSilt = 0.13 * self.PercentSilt / 100
            self.FracSand = self.PercentSand / 100 * (1 - self.PercentClay / 100) ** 2.4
            self.FracSagg = pcr.ifthenelse(
                self.PercentClay < 25,
                2.0 * self.PercentClay / 100,
                pcr.ifthenelse(
                    self.PercentClay > 50,
                    0.57,
                    0.28 * (self.PercentClay / 100 - 0.25) + 0.5,
                ),
            )
            self.FracLagg = (
                1.0 - self.FracClay - self.FracSilt - self.FracSand - self.FracSagg
            )

        """ Variables for the river transport model """

        if self.RunRiverModel == 1:

            # Parameters of Bagnold transport formula
            if self.RivTransportMethod == 2:
                self.cBagnold = self.readtblDefault(
                    self.Dir + "/" + self.intbl + "/c_Bagnold.tbl",
                    self.LandUse,
                    subcatch,
                    self.Soil,
                    0.0015,
                )
                self.expBagnold = self.readtblDefault(
                    self.Dir + "/" + self.intbl + "/exp_Bagnold.tbl",
                    self.LandUse,
                    subcatch,
                    self.Soil,
                    1.4,
                )

            # Parameters of Kodatie transport formula
            if self.RivTransportMethod == 3:
                self.aK = pcr.ifthenelse(
                    self.D50River <= 0.05,
                    281.4,
                    pcr.ifthenelse(
                        self.D50River <= 0.25,
                        2829.6,
                        pcr.ifthenelse(self.D50River <= 2, 2123.4, pcr.scalar(431884.8)),
                    ),
                )
                self.bK = pcr.ifthenelse(
                    self.D50River <= 0.05,
                    2.622,
                    pcr.ifthenelse(
                        self.D50River <= 0.25,
                        3.646,
                        pcr.ifthenelse(self.D50River <= 2, 3.300, pcr.scalar(1.0)),
                    ),
                )
                self.cK = pcr.ifthenelse(
                    self.D50River <= 0.05,
                    0.182,
                    pcr.ifthenelse(
                        self.D50River <= 0.25,
                        0.406,
                        pcr.ifthenelse(self.D50River <= 2, 0.468, pcr.scalar(1.0)),
                    ),
                )
                self.dK = pcr.ifthenelse(
                    self.D50River <= 0.05,
                    0.0,
                    pcr.ifthenelse(
                        self.D50River <= 0.25,
                        0.412,
                        pcr.ifthenelse(self.D50River <= 2, 0.613, pcr.scalar(2.0)),
                    ),
                )

            # Critical bed and bank shear stresses and erodibilities [N/m2] [m3/N.s]
            # Bank from Julian & Torres 2006 + Hanson & Simon 2001
            self.TCrBank = (
                0.1
                + 0.1779 * (100 * self.FracClayRiv + 100 * self.FracSiltRiv)
                + 0.0028 * (100 * self.FracClayRiv + 100 * self.FracSiltRiv) ** 2
                - 2.34
                * 10 ** (-5)
                * (100 * self.FracClayRiv + 100 * self.FracSiltRiv) ** 3
            ) * self.CovRiver
            self.KdBank = 0.2 * self.TCrBank ** (-0.5) * 10 ** (-6)
            # Bed from Shields diagram
            E = (
                (2.65 - 1)
                * 9.81
                * (self.D50River * 10 ** (-3)) ** 3
                / (1 * 10 ** (-12))
            ) ** (0.33)
            self.TCrBed = (
                (2.65 - 1)
                * 9.81
                * self.D50River
                * (
                    0.13 * E ** (-0.392) * pcr.exp(-0.015 * E ** 2)
                    + 0.045 * (1 - pcr.exp(-0.068 * E))
                )
            )
            self.KdBed = 0.2 * self.TCrBed ** (-0.5) * 10 ** (-6)

        """ Variables for lakes and reservoirs model """

        if hasattr(self, "ReserVoirSimpleLocs") or hasattr(
            self, "LakeLocs"
        ):
            self.ReserVoirLocs = self.ZeroMap
            self.filter_Eros_TC = self.ZeroMap + 1.0

        if hasattr(self, "ReserVoirSimpleLocs"):
            # Check if we have simple and or complex reservoirs
            self.ReserVoirSimpleLocs = pcr.nominal(self.ReserVoirSimpleLocs)
            self.ReservoirSimpleAreas = pcr.nominal(self.ReservoirSimpleAreas)
            tt_simple = pcr.pcr2numpy(self.ReserVoirSimpleLocs, 0.0)
            self.nrresSimple = np.size(np.where(tt_simple > 0.0)[0])
            self.ReserVoirLocs = self.ReserVoirLocs + pcr.cover(
                pcr.scalar(self.ReserVoirSimpleLocs), 0.0
            )
            
            res_area = pcr.cover(pcr.scalar(self.ReservoirSimpleAreas), 0.0)
            self.filter_Eros_TC = pcr.ifthenelse(
                pcr.boolean(pcr.cover(res_area, pcr.scalar(0.0))),
                res_area * 0.0,
                self.filter_Eros_TC,
            )
        else:
            self.nrresSimple = 0

        if hasattr(self, "LakeLocs"):
            self.LakeAreasMap = pcr.nominal(self.LakeAreasMap)
            self.LakeLocs = pcr.nominal(self.LakeLocs)
            tt_lake = pcr.pcr2numpy(self.LakeLocs, 0.0)
            #self.nrlake = tt_lake.max()
            self.nrlake = np.size(np.where(tt_lake > 0.0)[0])
            self.ReserVoirLocs = self.ReserVoirLocs + pcr.cover(
                pcr.scalar(self.LakeLocs), 0.0
            )
            lake_area = pcr.cover(pcr.scalar(self.LakeAreasMap), 0.0)
            self.filter_Eros_TC = pcr.ifthenelse(
                lake_area > 0, lake_area * 0.0, self.filter_Eros_TC
            )
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

            tt_filter = pcr.pcr2numpy(self.filter_Eros_TC, 1.0)
            self.filterResArea = tt_filter.min()
            
        
        #Compute length and width and correct for reservoirs/lakes
        # Factor on river length (self.RiverLengthFac) only used in combination with
        # calculated (by wflow_sbm) slope    
        # Set DCL to riverlength if that is longer that the basic length calculated from grid
        drainlength = detdrainlength(self.TopoLdd, self.xl, self.yl)
        # Multiply with Factor (taken from upscaling operation, defaults to 1.0 if no map is supplied)
        self.DCL = drainlength * pcr.max(1.0, self.RiverLengthFac)
        # Correct slope for extra length of the river in a gridcel
        riverslopecor = drainlength / self.DCL
        self.riverSlope = self.Slope * riverslopecor
        
        # If river slope available as map, also provide river length 
        self.riverSlope = pcr.max(
                pcr.scalar(0.00001), 
                self.wf_readmap(os.path.join(self.Dir, "staticmaps/RiverSlope.map"), self.riverSlope)
                )
        if os.path.isfile(os.path.join(self.Dir, wflow_riverlength)):
            self.DCL = self.RiverLength # m


        # Determine river width from DEM, upstream area and yearly average discharge
        # Scale yearly average Q at outlet with upstream are to get Q over whole catchment
        # Alf ranges from 5 to > 60. 5 for hardrock. large values for sediments
        # "Noah J. Finnegan et al 2005 Controls on the channel width of rivers:
        # Implications for modeling fluvial incision of bedrock"
        if (self.nrresSimple + self.nrlake) > 0:
            upstr = pcr.catchmenttotal(1, self.TopoLddOrg)
        else:
            upstr = pcr.catchmenttotal(1, self.TopoLdd)
        Qscale = upstr / pcr.mapmaximum(upstr) * Qmax
        W = (
            (alf * (alf + 2.0) ** (0.6666666667)) ** (0.375)
            * Qscale ** (0.375)
            * (pcr.max(0.0001, pcr.windowaverage(self.riverSlope, pcr.celllength() * 4.0)))
            ** (-0.1875)
            * self.NRiver ** (0.375)
        )
        # Use supplied riverwidth if possible, else calulate
        self.RiverWidth = pcr.ifthenelse(self.RiverWidth <= 0.0, W, self.RiverWidth)
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
        
        # water depth (m)
        # set width for kinematic wave to cell width for all cells
        self.Bw = detdrainwidth(self.TopoLdd, self.xl, self.yl)
        # However, in the main river we have real flow so set the width to the
        # width of the river

        self.Bw = pcr.ifthenelse(self.River, self.RiverWidth, self.Bw)

#        # Add rivers to the WaterFrac, but check with waterfrac map and correct
#        self.RiverFrac = pcr.min(
#            1.0,
#            pcr.ifthenelse(
#                self.River, (self.RiverWidth * self.DCL) / (self.xl * self.yl), 0
#            ),
#        )
#        
#        self.WaterFrac = pcr.max(self.WaterFrac - self.RiverFrac, 0)

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

        if self.reinit == 1:
            self.logger.warn("Setting initial states to default")
            if self.RunRiverModel == 1:
                self.SedLoad = self.ZeroMap
                self.ClayLoad = self.ZeroMap
                self.SiltLoad = self.ZeroMap
                self.SandLoad = self.ZeroMap
                self.SaggLoad = self.ZeroMap
                self.LaggLoad = self.ZeroMap
                self.GravLoad = self.ZeroMap

                self.OutSedLoad = self.ZeroMap
                self.OutClayLoad = self.ZeroMap
                self.OutSiltLoad = self.ZeroMap
                self.OutSandLoad = self.ZeroMap
                self.OutSaggLoad = self.ZeroMap
                self.OutLaggLoad = self.ZeroMap
                self.OutGravLoad = self.ZeroMap

                self.RivStoreSed = self.ZeroMap
                self.RivStoreClay = self.ZeroMap
                self.RivStoreSilt = self.ZeroMap
                self.RivStoreSand = self.ZeroMap
                self.RivStoreSagg = self.ZeroMap
                self.RivStoreLagg = self.ZeroMap
                self.RivStoreGrav = self.ZeroMap
                
                if hasattr(self, "ReserVoirSimpleLocs"):
                    self.ReservoirSedStore = self.wf_readmap(
                            os.path.join(self.Dir, "staticmaps", "ReservoirSedStore.map"),
                            0.0                        
                    )
                if hasattr(self, "LakeLocs"):
                    self.LakeSedStore = self.ZeroMap
        else:
            self.logger.info("Setting initial conditions from state files")
            self.wf_resume(os.path.join(self.Dir, "instate"))

    def default_summarymaps(self):
        """
      *Optional*

      Return a default list of variables to report as summary maps in the outsum dir.
      The ini file has more options, including average and sum
      """
        lst = [
            "self.PercentSand",
            "self.UsleK",
            "self.UsleC",
            "self.ErodK",
#            "self.D50River",
            "self.FracClayRiv",
            "self.FracSiltRiv",
            "self.FracSandRiv",
            "self.FracGravRiv",
            "self.filter_Eros_TC",
            "self.riverSlope",
            "self.RiverWidth",
            "self.RiverDem",
            "self.UpArea",
            "self.cellareaKm",
        ]

        return lst

    def dynamic(self):
        """
      *Required*
      
      This is where all the time dependent functions are executed. Time dependent
      output should also be saved here.
      
      Below a list of the main variables that can be saved to disk as maps or as
      timeseries.
      
      *Soil loss model*
      :var self.SedSpl: Sediment eroded by rainfall [ton]
      :var self.SedOv: Sediment eroded by overland flow [ton]
      :var self.SoilLoss: Total eroded sediment [ton]
      
      *Inland sediment routing*
      :var self.InLandSed: Inland eroded sediment in surface runoff reaching the river [ton]
      
      *River model*
      :var self.SedLoad: Sediment river load [ton]
      :var self.InSedLoad: Sediment load at the beginning of the time step [ton]
          (previous river load+outload from upstream river cells+load from land erosion)
      :var self.RivErodSed: Eroded sediment from resuspension, bed and bank erosion [ton]
      :var self.DepSedLoad: Sediment deposited in the river channel [ton]
      :var self.OutSedLoad: Sediment transported downstream [ton]   
      
      :var self.SedConc: Sediment concentration in the river [mg/L]
      :var self.SSConc: Suspended sediment concentration [mg/L]
      :var self.BedConc: Bed sediment concentration [mg/L]
      
      """

        self.wf_updateparameters()
        self.Precipitation = pcr.max(0.0, self.Precipitation)

        ##########################################################################
        # Soil loss model #############################
        ##########################################################################

        """ Splash / Rainfall erosion """
        # From EUROSEM
        if self.RainErodMethod == 1:
            # calculate rainfall intensity [mm/h]
            self.rintnsty = self.Precipitation / (self.timestepsecs / 3600)

            # Kinectic energy of direct throughfall [J/m^2/mm]
            # self.KeDirect = max(11.87 + 8.73 * log10(max(0.0001,self.rintnsty)), 0.0)  #basis used in USLE
            self.KeDirect = pcr.max(
                8.95 + 8.44 * pcr.log10(pcr.max(0.0001, self.rintnsty)), 0.0
            )  # variant, most used in distributed models

            # Kinetic energy of leaf drainage [J/m^2/mm]
            pheff = 0.5 * self.CanopyHeight  # [m]
            self.KeLeaf = pcr.max((15.8 * pheff ** 0.5) - 5.87, 0.0)

            # Depths of rainfall (total, leaf drainage, direct) [mm]
            # rdepth_tot = max(self.Precipitation/self.timestepsecs, 0.0)
            rDepthTot = pcr.max(self.Precipitation, 0.0)
            rDepthLeaf = pcr.max(rDepthTot * 0.1 * self.CanopyGapFraction, 0.0)  # stemflow
            rDepthDirect = pcr.max(
                rDepthTot - rDepthLeaf - self.Interception, 0.0
            )  # throughfall

            # Total kinetic energy by rainfall [J/m^2]
            self.KeTotal = (
                rDepthDirect * self.KeDirect + rDepthLeaf * self.KeLeaf
            ) * 0.001

            # Rainfall/Splash erosion
            self.SedSpl = (
                self.ErosK * self.KeTotal * pcr.exp(-self.ErosSpl * self.WaterLevelL)
            )  # [g/m^2]
            self.Sedspl = (
                self.cellareaKm * self.SedSpl
            )  # * self.timestepsecs # [ton/cell/timestep]

        if self.RainErodMethod == 2:
            # calculate rainfall intensity [mm/min]
            self.rintnsty = self.Precipitation / (self.timestepsecs / 60)

            # Splash erosion
            self.SedSpl = (
                0.108
                * self.UsleC
                * self.UsleK
                * self.cellareaKm
                * 10 ** 6
                * self.rintnsty ** 2
            )  # [kg/min]
            self.SedSpl = (
                self.timestepsecs / 60.0 * 10 ** (-3) * self.SedSpl
            )  # [ton/timestep]

        # Remove the impervious areas
        self.SedSpl = self.SedSpl * (1.0 - self.PathFrac)
        # Remove nodata values
        self.SedSpl = pcr.cover(self.SedSpl, self.ZeroMap)

        """ Overland flow erosion from ANSWERS"""
        # Only calculate overland flow erosion outside of river cells
        self.LandRunoffRate = self.LandRunoff * 60 / self.reallength  # [m2/min]

        # Overland flow erosion
        # For a wide range of slope, it is better to use the sine of slope rather than tangeant
        self.SedOv = (
            self.ErosOv
            * self.UsleC
            * self.UsleK
            * self.cellareaKm
            * 10 ** 6
            * self.sinSlope
            * self.LandRunoffRate
        )  # [kg/min]
        self.SedOv = (
            self.timestepsecs / 60.0 * 10 ** (-3) * self.SedOv
        )  # [ton/timestep]

        # Remove the impervious areas
        self.SedOv = self.SedOv * (1.0 - self.PathFrac)
        # Remove nodata values
        self.SedOv = pcr.cover(self.SedOv, self.ZeroMap)

        """ Total soil detachment """
        self.SoilLoss = self.SedSpl + self.SedOv  # [ton/cell/timestep]

        # Remove land erosion for reservoir cells
        if (self.nrresSimple + self.nrlake) > 0 and self.filterResArea == 0:
            self.SedSpl = self.filter_Eros_TC * self.SedSpl
            self.SedOv = self.filter_Eros_TC * self.SedOv
            self.SoilLoss = self.filter_Eros_TC * self.SoilLoss

        ##########################################################################
        # Inland sediment routing model              #############################
        ##########################################################################

        # If the river model is run, use Yalin's equation with particle differentiation.
        # If just the soil loss model is run, use Govers equation without particle differentiation

        """ Transport of overland flow with no particle differenciation using Govers equation"""

        if self.LandTransportMethod == 2 or self.LandTransportMethod == 3:
            if self.LandTransportMethod == 2:
                # Unit stream power
                self.velocityL = pcr.cover(
                    pcr.ifthenelse(
                        self.WaterLevelL > 0,
                        self.LandRunoff / (self.reallength * self.WaterLevelL),
                        0.0,
                    ),
                    0.0,
                )  # [m/s]
                self.omega = 10 * self.sinSlope * 100 * self.velocityL  # [cm/s]
                # self.omega = self.sinSlope * 100 * self.velocity #[cm/s]

                # Transport capacity from Govers, 1990
                self.TCf = pcr.ifthenelse(
                    self.omega > 0.4,
                    self.cGovers * (self.omega - 0.4) ** self.nGovers * 2650,
                    0.0,
                )  # [kg/m3]
                # self.TC = self.TCf / (1 - self.TCf/2650) #[kg/m3]
                # self.TC = max(self.TC, 2650)
                self.TC = (
                    self.TCf * self.LandRunoff * self.timestepsecs * 10 ** (-3)
                )  # [ton/cell/timestep]
                # Remove nodata values
                self.TC = pcr.cover(self.TC, self.ZeroMap)
                # Assume that eroded soil on lake cells all reach the river cells of the reservoir
                if (
                    self.nrresSimple + self.nrlake
                ) > 0 and self.filterResArea == 0:
                    self.TC = pcr.ifthenelse(
                        pcr.pcrand(self.filter_Eros_TC == 0, self.River == 0),
                        10 ** 9,
                        self.TC,
                    )

            elif self.LandTransportMethod == 3:
                # Transport capacity from Yalin
                self.delta = pcr.max(
                    self.WaterLevelL
                    * self.sinSlope
                    / (self.D50 * 10 ** (-3) * (self.rhoSed / 1000 - 1))
                    / 0.06
                    - 1,
                    0.0,
                )
                self.TC = (
                    self.reallength
                    / self.LandRunoff
                    * (self.rhoSed - 1000)
                    * self.D50
                    * 10 ** (-3)
                    * (9.81 * self.WaterLevelL * self.sinSlope)
                    * 0.635
                    * self.delta
                    * (
                        1
                        - pcr.ln(1+self.delta*2.45/(self.rhoSed / 1000)**0.4 *0.06** 0.5)
                        / self.delta
                        * 2.45
                        / (self.rhoSed / 1000) ** 0.4
                        * 0.06 ** 0.5
                    )
                )  # [kg/m3]
                self.TC = pcr.cover(
                        self.TC * self.LandRunoff * self.timestepsecs * 10 ** (-3), 
                        self.ZeroMap
                        )  # [ton/cell/timestep]
                
                # Assume that eroded soil on reservoir/lake cells all reach the river cells of the reservoir/lake
                if (self.nrresSimple + self.nrlake)>0 and self.filterResArea == 0:
                    self.TC = pcr.ifthenelse(
                        pcr.pcrand(self.filter_Eros_TC == 0, self.River == 0),
                        10 ** 9,
                        self.TC,
                    )

            # To get total sediment input from land into the river systems, river cells transport all sediment to the output (huge TC)
            self.TCRiv = pcr.cover(pcr.ifthenelse(self.River == 1, 10 ** 9, self.TC), self.TC)
            # Transported sediment over the land
            self.SedFlux = pcr.accucapacityflux(
                self.TopoLdd, self.SoilLoss, self.TCRiv
            )  # [ton/cell/tinestep]
            # Deposited sediment over the land
            self.SedDep = pcr.accucapacitystate(self.TopoLdd, self.SoilLoss, self.TCRiv)

            # Sediment amount reaching each river cell '''
            self.SedDep2 = pcr.accucapacitystate(self.TopoLdd, self.SoilLoss, self.TC)
            # Remove inland deposition
            self.OvSed = pcr.cover(pcr.ifthenelse(self.River == 1, self.SedDep2, 0.0), 0.0)


            """ Transport of overland flow with particle differenciation using Yalin equation"""
        else:
            # Determine the eroded amount of clay/silt/sand/aggregates on each cell [ton/cell/timestep]
            self.LandErodClay = self.SoilLoss * self.FracClay
            self.LandErodSilt = self.SoilLoss * self.FracSilt
            self.LandErodSand = self.SoilLoss * self.FracSand
            self.LandErodSagg = self.SoilLoss * self.FracSagg
            self.LandErodLagg = self.SoilLoss * self.FracLagg

            
            # Delta parameter of Yalin for each particle class
            deltaCoeff = (self.WaterLevelL* self.sinSlope / (10 ** (-6) 
                            * (self.rhoSed / 1000 - 1)) / 0.06)
            self.DClay = pcr.max((1/self.dmClay * deltaCoeff -1), 0.0)
            self.DSilt = pcr.max((1/self.dmSilt * deltaCoeff -1), 0.0)
            self.DSand = pcr.max((1/self.dmSand * deltaCoeff -1), 0.0)
            self.DSagg = pcr.max((1/self.dmSagg * deltaCoeff -1), 0.0)
            self.DLagg = pcr.max((1/self.dmLagg * deltaCoeff -1), 0.0)

            # Total transportability
            self.Dtot = self.DClay + self.DSilt + self.DSand + self.DSagg + self.DLagg

            # Yalin Transport capacity of overland flow for each particle class
            TCa = (self.reallength / self.LandRunoff * (self.rhoSed - 1000) * 10 ** (-6) 
                    * (9.81 * self.WaterLevelL * self.sinSlope))
            TCb = 2.45 / (self.rhoSed / 1000) ** 0.4 * 0.06 ** 0.5
            self.TCClay = (
                TCa * self.dmClay
                * self.DClay / self.Dtot
                * 0.635 * self.DClay
                * (1 - pcr.ln(1 + self.DClay * TCb) / self.DClay * TCb)
            )  # [kg/m3]
            self.TCClay = pcr.cover(
                    self.TCClay * self.LandRunoff * self.timestepsecs * 10 ** (-3), 
                    self.ZeroMap
                    ) # [ton/cell/timestep]

            self.TCSilt = (
                TCa * self.dmSilt
                * self.DSilt / self.Dtot
                * 0.635 * self.DClay
                * (1 - pcr.ln(1 + self.DSilt * TCb) / self.DSilt * TCb)
            )  # [kg/m3]
            self.TCSilt = pcr.cover(
                    self.TCSilt * self.LandRunoff * self.timestepsecs * 10 ** (-3), 
                    self.ZeroMap
                    ) # [ton/cell/timestep]

            self.TCSand = (
                TCa * self.dmSand
                * self.DSand / self.Dtot
                * 0.635 * self.DSand
                * (1 - pcr.ln(1 + self.DSand * TCb) / self.DSand * TCb)
            )  # [kg/m3]
            self.TCSand = pcr.cover(
                    self.TCSand * self.LandRunoff * self.timestepsecs * 10 ** (-3), 
                    self.ZeroMap) # [ton/cell/timestep]

            self.TCSagg = (
                TCa * self.dmSagg
                * self.DSagg / self.Dtot
                * 0.635 * self.DSagg
                * (1 - pcr.ln(1 + self.DSagg * TCb) / self.DSagg * TCb)
            )  # [kg/m3]
            self.TCSagg = pcr.cover(
                    self.TCSagg * self.LandRunoff * self.timestepsecs * 10 ** (-3), 
                    self.ZeroMap
                    ) # [ton/cell/timestep]

            self.TCLagg = (
                TCa * self.dmLagg
                * self.DLagg / self.Dtot
                * 0.635 * self.DLagg
                * (1 - pcr.ln(1 + self.DLagg * TCb) / self.DLagg * TCb)
            )  # [kg/m3]
            self.TCLagg = pcr.cover(
                    self.TCLagg * self.LandRunoff * self.timestepsecs * 10 ** (-3), 
                    self.ZeroMap
                    ) # [ton/cell/timestep]

            # Assume that eroded soil on lake cells all reach the river cells of the reservoir
            if (self.nrresSimple + self.nrlake) > 0:
                self.TCClay = pcr.cover(
                    pcr.ifthenelse(self.filter_Eros_TC == 0, 10 ** 9, self.TCClay),
                    self.TCClay,
                )
                self.TCClay = pcr.cover(
                    pcr.ifthenelse(self.River == 1, 0.0, self.TCClay), self.TCClay
                )
                self.TCSilt = pcr.cover(
                    pcr.ifthenelse(self.filter_Eros_TC == 0, 10 ** 9, self.TCSilt),
                    self.TCSilt,
                )
                self.TCSilt = pcr.cover(
                    pcr.ifthenelse(self.River == 1, 0.0, self.TCSilt), self.TCSilt
                )
                self.TCSand = pcr.cover(
                    pcr.ifthenelse(self.filter_Eros_TC == 0, 10 ** 9, self.TCSand),
                    self.TCSand,
                )
                self.TCSand = pcr.cover(
                    pcr.ifthenelse(self.River == 1, 0.0, self.TCSand), self.TCSand
                )
                self.TCSagg = pcr.cover(
                    pcr.ifthenelse(self.filter_Eros_TC == 0, 10 ** 9, self.TCSagg),
                    self.TCSagg,
                )
                self.TCSagg = pcr.cover(
                    pcr.ifthenelse(self.River == 1, 0.0, self.TCSagg), self.TCSagg
                )
                self.TCLagg = pcr.cover(
                    pcr.ifthenelse(self.filter_Eros_TC == 0, 10 ** 9, self.TCLagg),
                    self.TCLagg,
                )
                self.TCLagg = pcr.cover(
                    pcr.ifthenelse(self.River == 1, 0.0, self.TCLagg), self.TCLagg
                )

            # Eroded sediment in overland flow reaching the river system per particle class [ton/cell/timestep]
            self.InLandClay = pcr.cover(
                pcr.ifthenelse(
                    self.River == 1,
                    pcr.accucapacitystate(self.TopoLdd, self.LandErodClay, self.TCClay),
                    0.0,
                ),
                0.0,
            )
            self.InLandSilt = pcr.cover(
                pcr.ifthenelse(
                    self.River == 1,
                    pcr.accucapacitystate(self.TopoLdd, self.LandErodSilt, self.TCSilt),
                    0.0,
                ),
                0.0,
            )
            self.InLandSand = pcr.cover(
                pcr.ifthenelse(
                    self.River == 1,
                    pcr.accucapacitystate(self.TopoLdd, self.LandErodSand, self.TCSand),
                    0.0,
                ),
                0.0,
            )
            self.InLandSagg = pcr.cover(
                pcr.ifthenelse(
                    self.River == 1,
                    pcr.accucapacitystate(self.TopoLdd, self.LandErodSagg, self.TCSagg),
                    0.0,
                ),
                0.0,
            )
            self.InLandLagg = pcr.cover(
                pcr.ifthenelse(
                    self.River == 1,
                    pcr.accucapacitystate(self.TopoLdd, self.LandErodLagg, self.TCLagg),
                    0.0,
                ),
                0.0,
            )
            self.InLandSed = (
                self.InLandClay
                + self.InLandSilt
                + self.InLandSand
                + self.InLandSagg
                + self.InLandLagg
            )  # total


        ##########################################################################
        # River transport and processes              #############################
        ##########################################################################

        if self.RunRiverModel == 1:
            # River sediment loads are separated into different particle class.
            # Clay, silt and sand can both come from land, resuspension or river channel erosion.
            # Small and large aggregates only come from land erosion or resuspension.
            # Gravel only comes from resuspension or river channel erosion.
            
            ### Equations independant of the iterations ###
            """River erosion """
            # Hydraulic radius of the river [m] (rectangular channel)
            self.HydRad = (
                self.WaterLevelR
                * self.RiverWidth
                / (self.RiverWidth + 2 * self.WaterLevelR)
            )
            
            # Transport capacity
            # Engelund and Hansen transport formula
            if self.RivTransportMethod == 1:
                # Shields parameter
                self.ThetaShields = (
                    1000
                    * 9.81
                    * self.HydRad
                    * self.riverSlope
                    / ((self.rhoSed - 1000) * 9.81 * self.D50River / 1000)
                )
                self.vmean = pcr.ifthenelse(
                    self.WaterLevelR > 0,
                    self.RiverRunoff / (self.RiverWidth * self.WaterLevelR),
                    pcr.scalar(0.0),
                )
                self.vshear = (9.81 * self.HydRad * self.riverSlope) ** (0.5)
                #          self.Cw = 0.05 * (2650/(2650-1000)) * self.vmean * self.riverSlope \
                #              / ((2650-1000)/1000*9.81*self.D50River/1000)**(0.5) * self.ThetaShields**(0.5) # sediment concentration by weight
                self.Cw = pcr.min(
                    1.0,
                    pcr.ifthenelse(
                        pcr.pcrand(self.HydRad > 0, self.River == 1),
                        self.rhoSed/1000
                        * 0.05
                        * self.vmean
                        * self.vshear ** 3
                        / ((self.rhoSed/1000 - 1) ** 2 * 9.81 ** 2 *self.D50River*self.HydRad),
                        pcr.scalar(0.0),
                    ),
                )  # concentration by weight
                self.MaxSedLoad = pcr.max(
                    self.Cw / (self.Cw + (1 - self.Cw) * self.rhoSed/1000) * self.rhoSed/1000,
                    self.ZeroMap
                )  # [tons/m3]
            
            # Simplified Bagnold transport formula
            if self.RivTransportMethod == 2:
                self.MaxSedLoad = (
                    self.cBagnold
                    * (self.RiverRunoff / (self.WaterLevelR * self.RiverWidth))
                    ** self.expBagnold
                )  # [ton/m3]

            # Kodatie transport formula
            if self.RivTransportMethod == 3:
                self.vmean = pcr.ifthenelse(
                    self.WaterLevelR > 0,
                    self.RiverRunoff / (self.RiverWidth * self.WaterLevelR),
                    pcr.scalar(0.0),
                )
                self.MaxSedLoad = (
                    self.aK
                    * self.vmean ** (self.bK)
                    * self.WaterLevelR ** (self.cK)
                    * self.riverSlope ** (self.dK)
                ) * (
                    self.RiverWidth
                )  # [tons]
                
            # Yang transport formula
            if self.RivTransportMethod == 4:
                self.wsRiv = 411 * self.D50River ** 2 / 3600
                self.vshear = (9.81 * self.HydRad * self.riverSlope) ** (0.5)
                self.var1 = self.vshear * self.D50River / 1000 / (1.16 * 10 ** (-6))
                self.var2 = self.wsRiv * self.D50River / 1000 / (1.16 * 10 ** (-6))
                self.vcr = pcr.ifthenelse(
                    self.var1 >= 70,
                    2.05 * self.wsRiv,
                    self.wsRiv * (2.5 / (pcr.log10(self.var1) - 0.06) + 0.66),
                )

                # Sand equation
                self.logCppm = (
                    5.435
                    - 0.286 * pcr.log10(self.var2)
                    - 0.457 * pcr.log10(self.vshear / self.wsRiv)
                    + (
                        1.799
                        - 0.409 * pcr.log10(self.var2)
                        - 0.314 * pcr.log10(self.vshear / self.wsRiv)
                    )
                    * pcr.log10(
                        (
                            self.RiverRunoff / (self.RiverWidth * self.WaterLevelR)
                            - self.vcr
                        )
                        * self.riverSlope
                        / self.wsRiv
                    )
                )
                # Gravel equation
                self.logCppm = pcr.ifthenelse(
                    self.D50River < 2.0,
                    self.logCppm,
                    6.681
                    - 0.633 * pcr.log10(self.var2)
                    - 4.816 * pcr.log10(self.vshear / self.wsRiv)
                    + (
                        2.784
                        - 0.305 * pcr.log10(self.var2)
                        - 0.282 * pcr.log10(self.vshear / self.wsRiv)
                    )
                    * pcr.log10(
                        (
                            self.RiverRunoff / (self.RiverWidth * self.WaterLevelR)
                            - self.vcr
                        )
                        * self.riverSlope
                        / self.wsRiv
                    ),
                )
                self.Cw = 10 ** self.logCppm * 10**(-6)  # sediment concentration by weight
                self.MaxSedLoad = pcr.max(
                    self.Cw / (self.Cw + (1 - self.Cw) * self.rhoSed/1000) * self.rhoSed/1000,
                    self.ZeroMap
                )  # [tons/m3]

            # Molinas & Wu transport formula
            if self.RivTransportMethod == 5:
                self.wsRiv = 411 * self.D50River ** 2 / 3600
                self.psi = (
                    (self.rhoSed/1000 - 1)
                    * 9.81
                    * self.WaterLevelR
                    * self.wsRiv
                    * (pcr.log10(self.WaterLevelR / self.D50River)) ** 2
                ) ** (0.5)
                self.Cw = (
                    1430
                    * (0.86 + (self.psi) ** (0.5))
                    * (self.psi) ** (1.5)
                    / (0.016 + self.psi)
                    * 10 ** (-6)
                )  # weight
                self.MaxSedLoad = pcr.max(
                    self.Cw / (self.Cw + (1 - self.Cw) * self.rhoSed/1000) * self.rhoSed/1000,
                    self.ZeroMap
                )  # [tons/m3]
                
            
            # Repartition of the effective shear stress between the bank and the bed from Knight et al. 1984
            self.SFBank = pcr.ifthenelse(
                self.WaterLevelR > 0,
                pcr.exp(-3.230 * pcr.log10(self.RiverWidth / self.WaterLevelR + 3) + 6.146),
                pcr.scalar(0.0),
            )  # [%]
            # Effective shear stress on river bed and banks [N/m2]
            self.TEffBank = pcr.ifthenelse(
                self.WaterLevelR > 0,
                1000
                * 9.81
                * self.HydRad
                * self.riverSlope
                * self.SFBank
                / 100
                * (1 + self.RiverWidth / (2 * self.WaterLevelR)),
                pcr.scalar(0.0),
            )
            self.TEffBed = (
                1000
                * 9.81
                * self.HydRad
                * self.riverSlope
                * (1 - self.SFBank / 100)
                * (1 + 2 * self.WaterLevelR / self.RiverWidth)
            )
            
            """ Depositio/settling """
            # Fractions of deposited particles in river cells from the Einstein formula [-]
            # Particle fall velocity [m/s] from Stokes
            self.x = pcr.ifthenelse(
                self.RiverRunoff > 0,
                1.055 * self.DCL / (self.RiverRunoff / self.RiverWidth),
                pcr.scalar(0.0),
            )
            self.x = pcr.cover(self.x, self.ZeroMap)
            
            self.xClay = pcr.min(1.0, 1 - 1 / pcr.exp(self.x * (411 * (self.dmClay/1000) ** 2 / 3600)))
            self.xSilt = pcr.min(1.0, 1 - 1 / pcr.exp(self.x * (411 * (self.dmSilt/1000) ** 2 / 3600)))
            self.xSand = pcr.min(1.0, 1 - 1 / pcr.exp(self.x * (411 * (self.dmSand/1000) ** 2 / 3600)))
            self.xSagg = pcr.min(1.0, 1 - 1 / pcr.exp(self.x * (411 * (self.dmSagg/1000) ** 2 / 3600)))
            self.xLagg = pcr.min(1.0, 1 - 1 / pcr.exp(self.x * (411 * (self.dmLagg/1000) ** 2 / 3600)))
            self.xGrav = pcr.min(1.0, 1 - 1 / pcr.exp(self.x * (411 * (self.dmGrav/1000) ** 2 / 3600)))
            
            
            ### Time iterations of the river processes ###
            #Compute the number of iterations if necessary
            it_sedR = 1
            if self.transportIters == 1:
                if self.transportRiverTstep == 0:
                    it_sedR = estimate_sedriv_tranport_iter(self.RiverRunoff, 
                                                            self.timestepsecs, 
                                                            self.WaterLevelR, 
                                                            self.DCL, 
                                                            self.RiverWidth, 
                                                            self.mv
                                                            )
                else:
                    it_sedR = int(np.ceil(self.timestepsecs/self.transportRiverTstep))

            #Start iterations of the river transport
            for v in range(0,it_sedR):
                """ Initial concentration and input from land and upstream river cells """
                # Save the loads from the previous time step that remained in the cell
                self.OldSedLoad = self.SedLoad
                self.OldClayLoad = self.ClayLoad
                self.OldSiltLoad = self.SiltLoad
                self.OldSandLoad = self.SandLoad
                self.OldSaggLoad = self.SaggLoad
                self.OldLaggLoad = self.LaggLoad
                self.OldGravLoad = self.GravLoad
    
                # Input concentration including the old load, the incoming load from upstream river cells and the incoming load from land erosion
                self.InSedLoad = (
                    self.OldSedLoad
                    + pcr.upstream(self.TopoLdd, self.OutSedLoad)
                    + self.InLandSed / it_sedR
                )
                self.InClayLoad = (
                    self.OldClayLoad
                    + pcr.upstream(self.TopoLdd, self.OutClayLoad)
                    + self.InLandClay / it_sedR
                )
                self.InSiltLoad = (
                    self.OldSiltLoad
                    + pcr.upstream(self.TopoLdd, self.OutSiltLoad)
                    + self.InLandSilt / it_sedR
                )
                self.InSandLoad = (
                    self.OldSandLoad
                    + pcr.upstream(self.TopoLdd, self.OutSandLoad)
                    + self.InLandSand / it_sedR
                )
                self.InSaggLoad = (
                    self.OldSaggLoad
                    + pcr.upstream(self.TopoLdd, self.OutSaggLoad)
                    + self.InLandSagg / it_sedR
                )
                self.InLaggLoad = (
                    self.OldLaggLoad
                    + pcr.upstream(self.TopoLdd, self.OutLaggLoad)
                    + self.InLandLagg / it_sedR
                )
                self.InGravLoad = self.OldGravLoad + pcr.upstream(
                    self.TopoLdd, self.OutGravLoad
                )
    
                """ River erosion """    
                # Transport capacity in [tons]
                if self.RivTransportMethod != 3:
                    self.MaxSedLoad = self.MaxSedLoad * (
                        self.WaterLevelR * self.RiverWidth * self.DCL
                        + self.RiverRunoff * self.timestepsecs / it_sedR
                    )  # [ton]
    
                self.MaxSedLoad = pcr.cover(self.MaxSedLoad, self.ZeroMap)
    
                # Potential erosion rates of the bed and bank [t/cell/timestep] 
                #(assuming only one bank is eroding)
                self.ERBank = pcr.max(
                    0.0,
                    self.KdBank
                    * (self.TEffBank - self.TCrBank)
                    * (self.DCL * self.WaterLevelR)
                    * 1.4
                    * self.timestepsecs / it_sedR,
                )  # 1.4 is bank default bulk density
                self.ERBed = pcr.max(
                    0.0,
                    self.KdBed
                    * (self.TEffBed - self.TCrBed)
                    * (self.DCL * self.RiverWidth)
                    * 1.5
                    * self.timestepsecs / it_sedR,
                )  # 1.5 is bed default bulk density
                # Relative potential erosion rates of the bed and bank [-]
                self.RTEBank = pcr.ifthenelse(
                    self.ERBank + self.ERBed > 0.0,
                    self.ERBank / (self.ERBank + self.ERBed),
                    0.0,
                )
                self.RTEBed = 1.0 - self.RTEBank
    
                # Excess transport capacity [ton/cell/timestep]
                # Erosion only if the load is below the transport capacity of the flow
                self.SedEx = pcr.max(self.MaxSedLoad - self.InSedLoad, 0.0)
                # Bed and bank are eroded after the previously deposited material
                self.EffSedEx = pcr.max(self.SedEx - self.RivStoreSed, 0.0)
    
                # Bank erosion [ton/cell/timestep]
                self.BankSedLoad = pcr.ifthenelse(
                    self.EffSedEx == 0,
                    0.0,
                    pcr.ifthenelse(
                        self.EffSedEx * self.RTEBank <= self.ERBank,
                        self.EffSedEx * self.RTEBank,
                        self.ERBank,
                    ),
                )
                self.BankClay = self.FracClayRiv * self.BankSedLoad
                self.BankSilt = self.FracSiltRiv * self.BankSedLoad
                self.BankSand = self.FracSandRiv * self.BankSedLoad
                self.BankGrav = self.FracGravRiv * self.BankSedLoad
    
                # Bed erosion [ton/cell/timestep]
                self.BedSedLoad = pcr.ifthenelse(
                    self.EffSedEx == 0,
                    0.0,
                    pcr.ifthenelse(
                        self.EffSedEx * self.RTEBed <= self.ERBed,
                        self.EffSedEx * self.RTEBed,
                        self.ERBed,
                    ),
                )
                self.BedClay = self.FracClayRiv * self.BedSedLoad
                self.BedSilt = self.FracSiltRiv * self.BedSedLoad
                self.BedSand = self.FracSandRiv * self.BedSedLoad
                self.BedGrav = self.FracGravRiv * self.BedSedLoad
    
                # Erosion/degradation of the previously deposited sediment (from clay to gravel) [ton/cell/timestep]
                self.DegStoreClay = pcr.ifthenelse(
                    self.RivStoreClay >= self.SedEx, self.SedEx, self.RivStoreClay
                )
                self.RivStoreClay = self.RivStoreClay - self.DegStoreClay  # update store
                self.SedEx = pcr.max(
                    self.SedEx - self.DegStoreClay, 0.0
                )  # update amount of sediment that need to be degraded
    
                self.DegStoreSilt = pcr.ifthenelse(
                    self.RivStoreSilt >= self.SedEx, self.SedEx, self.RivStoreSilt
                )
                self.RivStoreSilt = self.RivStoreSilt - self.DegStoreSilt
                self.SedEx = pcr.max(self.SedEx - self.DegStoreSilt, 0.0)
    
                self.DegStoreSagg = pcr.ifthenelse(
                    self.RivStoreSagg >= self.SedEx, self.SedEx, self.RivStoreSagg
                )
                self.RivStoreSagg = self.RivStoreSagg - self.DegStoreSagg
                self.SedEx = pcr.max(self.SedEx - self.DegStoreSagg, 0.0)
    
                self.DegStoreSand = pcr.ifthenelse(
                    self.RivStoreSand >= self.SedEx, self.SedEx, self.RivStoreSand
                )
                self.RivStoreSand = self.RivStoreSand - self.DegStoreSand
                self.SedEx = pcr.max(self.SedEx - self.DegStoreSand, 0.0)
    
                self.DegStoreLagg = pcr.ifthenelse(
                    self.RivStoreLagg >= self.SedEx, self.SedEx, self.RivStoreLagg
                )
                self.RivStoreLagg = self.RivStoreLagg - self.DegStoreLagg
                self.SedEx = pcr.max(self.SedEx - self.DegStoreLagg, 0.0)
    
                self.DegStoreGrav = pcr.ifthenelse(
                    self.RivStoreGrav >= self.SedEx, self.SedEx, self.RivStoreGrav
                )
                self.RivStoreGrav = self.RivStoreGrav - self.DegStoreGrav
                self.SedEx = pcr.max(self.SedEx - self.DegStoreGrav, 0.0)
    
                self.DegStoreSed = (
                    self.DegStoreClay
                    + self.DegStoreSilt
                    + self.DegStoreSand
                    + self.DegStoreSagg
                    + self.DegStoreLagg
                    + self.DegStoreGrav
                )
                self.RivStoreSed = self.RivStoreSed - self.DegStoreSed
    
                # Sum all erosion sources per particle class
                self.RivErodSed = self.BankSedLoad + self.BedSedLoad + self.DegStoreSed
                self.RivErodClay = self.BankClay + self.BedClay + self.DegStoreClay
                self.RivErodSilt = self.BankSilt + self.BedSilt + self.DegStoreSilt
                self.RivErodSand = self.BankSand + self.BedSand + self.DegStoreSand
                self.RivErodSagg = self.DegStoreSagg
                self.RivErodLagg = self.DegStoreLagg
                self.RivErodGrav = self.BankGrav + self.BedGrav + self.DegStoreGrav
    
                # Assume that there is no erosion/resuspension in reservoir
                if (self.nrresSimple + self.nrlake) > 0 and self.filterResArea == 0:
                    self.RivErodSed = self.filter_Eros_TC * self.RivErodSed
                    self.RivErodClay = self.filter_Eros_TC * self.RivErodClay
                    self.RivErodSilt = self.filter_Eros_TC * self.RivErodSilt
                    self.RivErodSand = self.filter_Eros_TC * self.RivErodSand
                    self.RivErodSagg = self.filter_Eros_TC * self.RivErodSagg
                    self.RivErodLagg = self.filter_Eros_TC * self.RivErodLagg
                    self.RivErodGrav = self.filter_Eros_TC * self.RivErodGrav
    
                """ Deposition/settling """
                # If transport capacity is exceeded, the excess is deposited.
                # Else classic deposition with Einstein formula.
                self.TCEx = self.MaxSedLoad - self.InSedLoad
                self.FracDepEx = pcr.ifthenelse(
                    self.TCEx < 0, 1 - self.MaxSedLoad / self.InSedLoad, pcr.scalar(0.0)
                )
                
                # Fractions of deposited particles in river cells from the Einstein formula [-]
                # Particle fall velocity [m/s] from Stokes
                self.FracDepClay = pcr.ifthenelse(self.TCEx >= 0, self.xClay, self.FracDepEx)
                self.FracDepSilt = pcr.ifthenelse(self.TCEx >= 0, self.xSilt, self.FracDepEx)
                self.FracDepSand = pcr.ifthenelse(self.TCEx >= 0, self.xSand, self.FracDepEx)
                self.FracDepSagg = pcr.ifthenelse(self.TCEx >= 0, self.xSagg, self.FracDepEx)
                self.FracDepLagg = pcr.ifthenelse(self.TCEx >= 0, self.xLagg, self.FracDepEx)
                self.FracDepGrav = pcr.ifthenelse(self.TCEx >= 0, self.xGrav, self.FracDepEx)
    
                # Sediment deposited in the channel [ton/cell/timestep]
                self.DepClay = self.FracDepClay * (self.InClayLoad + self.RivErodClay)
                self.DepSilt = self.FracDepSilt * (self.InSiltLoad + self.RivErodSilt)
                self.DepSand = self.FracDepSand * (self.InSandLoad + self.RivErodSand)
                self.DepSagg = self.FracDepSagg * (self.InSaggLoad + self.RivErodSagg)
                self.DepLagg = self.FracDepLagg * (self.InLaggLoad + self.RivErodLagg)
                self.DepGrav = self.FracDepGrav * (self.InGravLoad + self.RivErodGrav)
                self.DepSedLoad = (
                    self.DepClay
                    + self.DepSilt
                    + self.DepSand
                    + self.DepSagg
                    + self.DepLagg
                    + self.DepGrav
                )
    
                # Assume that deposition happens only on reservoir output cells where the reservoir mass balance takes place
                if (self.nrresSimple + self.nrlake) > 0:
                    self.DepClay = self.filter_Eros_TC * self.DepClay
                    self.DepSilt = self.filter_Eros_TC * self.DepSilt
                    self.DepSand = self.filter_Eros_TC * self.DepSand
                    self.DepSagg = self.filter_Eros_TC * self.DepSagg
                    self.DepLagg = self.filter_Eros_TC * self.DepLagg
                    self.DepGrav = self.filter_Eros_TC * self.DepGrav
                    self.DepSedLoad = self.filter_Eros_TC * self.DepSedLoad
    
                # Deposition in lakes from Camp 1945
                if self.nrlake > 0:
                    self.VcRes = pcr.ifthenelse(
                        self.LakeArea > 0, self.RiverRunoff / self.LakeArea, self.ZeroMap
                    )
                    DCRes = 411 / 3600 / self.VcRes
                    self.DepClay = pcr.ifthenelse(
                        pcr.cover(self.LakeArea, pcr.scalar(0.0)) > 0,
                        self.InClayLoad * pcr.min(1.0, (DCRes * (self.dmClay/1000) ** 2)),
                        self.DepClay
                        )
                    self.DepSilt = pcr.ifthenelse(
                        pcr.cover(self.LakeArea, pcr.scalar(0.0)) > 0,
                        self.InSiltLoad * pcr.min(1.0, (DCRes * (self.dmSilt/1000) ** 2)),
                        self.DepSilt
                        )
                    self.DepSand = pcr.ifthenelse(
                        pcr.cover(self.LakeArea, pcr.scalar(0.0)) > 0,
                        self.InSandLoad * pcr.min(1.0, (DCRes * (self.dmSand/1000) ** 2)),
                        self.DepSand
                        )
                    self.DepSagg = pcr.ifthenelse(
                        pcr.cover(self.LakeArea, pcr.scalar(0.0)) > 0,
                        self.InSaggLoad * pcr.min(1.0, (DCRes * (self.dmSagg/1000) ** 2)),
                        self.DepSagg
                        )
                    self.DepLagg = pcr.ifthenelse(
                        pcr.cover(self.LakeArea, pcr.scalar(0.0)) > 0,
                        self.InLaggLoad * pcr.min(1.0, (DCRes * (self.dmLagg/1000) ** 2)),
                        self.DepLagg
                        )
                    self.DepGrav = pcr.ifthenelse(
                        pcr.cover(self.LakeArea, pcr.scalar(0.0)) > 0,
                        self.InGravLoad * pcr.min(1.0, (DCRes * (self.dmGrav/1000) ** 2)),
                        self.DepGrav
                        )
                    self.DepSedLoad = (self.DepClay
                            + self.DepSilt
                            + self.DepSand
                            + self.DepSagg
                            + self.DepLagg
                            + self.DepGrav)
    
                    self.LakeSedStore = pcr.ifthen(self.LakeArea > 0, self.DepSedLoad)
                    
                #Deposition in reservoir, simple retention coeff for now
                if self.nrresSimple > 0:
                    
                    self.DepClay = pcr.ifthenelse(
                            pcr.cover(self.ResSimpleArea, pcr.scalar(0.0)) > 0,
                            0.8*self.InClayLoad,
                            self.DepClay
                            )
                    self.DepSilt = pcr.ifthenelse(
                            pcr.cover(self.ResSimpleArea, pcr.scalar(0.0)) > 0,
                            0.8*self.InSiltLoad,
                            self.DepSilt
                            )
                    self.DepSand = pcr.ifthenelse(
                            pcr.cover(self.ResSimpleArea, pcr.scalar(0.0)) > 0,
                            0.9*self.InSandLoad,
                            self.DepSand
                            )
                    self.DepSagg = pcr.ifthenelse(
                            pcr.cover(self.ResSimpleArea, pcr.scalar(0.0)) > 0,
                            0.9*self.InSaggLoad,
                            self.DepSagg
                            )
                    self.DepLagg = pcr.ifthenelse(
                            pcr.cover(self.ResSimpleArea, pcr.scalar(0.0)) > 0,
                            self.InLaggLoad,
                            self.DepLagg
                            )
                    self.DepGrav = pcr.ifthenelse(
                            pcr.cover(self.ResSimpleArea, pcr.scalar(0.0)) > 0,
                            self.InGravLoad,
                            self.DepGrav
                            )
                    self.DepSedLoad = (self.DepClay
                            + self.DepSilt
                            + self.DepSand
                            + self.DepSagg
                            + self.DepLagg
                            + self.DepGrav)
                    
                    self.ReservoirSedStore = pcr.ifthen(self.ResSimpleArea > 0, 
                                                        self.DepSedLoad)
    
                # Update the river deposited sediment storage
                self.RivStoreSed = self.RivStoreSed + self.DepSedLoad
                self.RivStoreClay = self.RivStoreClay + self.DepClay
                self.RivStoreSilt = self.RivStoreSilt + self.DepSilt
                self.RivStoreSand = self.RivStoreSand + self.DepSand
                self.RivStoreSagg = self.RivStoreSagg + self.DepSagg
                self.RivStoreLagg = self.RivStoreLagg + self.DepLagg
                self.RivStoreGrav = self.RivStoreGrav + self.DepGrav
    
                """ Output loads """
                # Sediment transported out of the cell during the timestep [ton/cell/timestep]
                # 0 in case all sediment are deposited in the cell
                # Reduce the fraction so that there is still some sediment staying in the river cell
                self.OutFracWat = pcr.min(
                    self.RiverRunoff
                    * self.timestepsecs / it_sedR
                    / (
                        self.WaterLevelR * self.RiverWidth * self.DCL
                    ),
                    1.0,
                )
    
                self.OutSedLoad = pcr.cover(
                        pcr.max(
                            0.0,
                            (self.InSedLoad + self.RivErodSed - self.DepSedLoad)
                            * self.OutFracWat,
                        ),
                        self.ZeroMap)
                self.OutClayLoad = pcr.cover(
                        pcr.max(
                            0.0,
                            (self.InClayLoad + self.RivErodClay - self.DepClay)
                            * self.OutFracWat,
                        ),
                        self.ZeroMap)
                self.OutSiltLoad = pcr.cover(
                        pcr.max(
                            0.0,
                            (self.InSiltLoad + self.RivErodSilt - self.DepSilt)
                            * self.OutFracWat,
                        ),
                        self.ZeroMap)
                self.OutSandLoad = pcr.cover(
                        pcr.max(
                            0.0,
                            (self.InSandLoad + self.RivErodSand - self.DepSand)
                            * self.OutFracWat,
                        ),
                        self.ZeroMap)
                self.OutSaggLoad = pcr.cover(
                        pcr.max(
                            0.0,
                            (self.InSaggLoad + self.RivErodSagg - self.DepSagg)
                            * self.OutFracWat,
                        ),
                        self.ZeroMap)
                self.OutLaggLoad = pcr.cover(
                        pcr.max(
                            0.0,
                            (self.InLaggLoad + self.RivErodLagg - self.DepLagg)
                            * self.OutFracWat,
                        ),
                        self.ZeroMap)
                self.OutGravLoad = pcr.cover(
                        pcr.max(
                            0.0,
                            (self.InGravLoad + self.RivErodGrav - self.DepGrav)
                            * self.OutFracWat,
                        ),
                        self.ZeroMap)
    
                """ Mass balance, sediment concentration in each river cell """
                # Sediment load [ton/cell/timestep]
                # 0 in case all sediment are deposited in the cell
                self.SedLoad = pcr.max(
                    0.0,
                    self.InSedLoad + self.RivErodSed - self.DepSedLoad - self.OutSedLoad,
                )
                self.ClayLoad = pcr.max(
                    0.0,
                    self.InClayLoad + self.RivErodClay - self.DepClay - self.OutClayLoad,
                )
                self.SiltLoad = pcr.max(
                    0.0,
                    self.InSiltLoad + self.RivErodSilt - self.DepSilt - self.OutSiltLoad,
                )
                self.SandLoad = pcr.max(
                    0.0,
                    self.InSandLoad + self.RivErodSand - self.DepSand - self.OutSandLoad,
                )
                self.SaggLoad = pcr.max(
                    0.0,
                    self.InSaggLoad + self.RivErodSagg - self.DepSagg - self.OutSaggLoad,
                )
                self.LaggLoad = pcr.max(
                    0.0,
                    self.InLaggLoad + self.RivErodLagg - self.DepLagg - self.OutLaggLoad,
                )
                self.GravLoad = pcr.max(
                    0.0,
                    self.InGravLoad + self.RivErodGrav - self.DepGrav - self.OutGravLoad,
                )

            #End of the transport iterations
            
            # Sediment concentration [mg/L]
            # The suspended load is composed of clay and silt.
            # The bed load is composed of sand, small and large aggregates from overland flow and gravel coming from river bed erosion.
            # Conversion from load [ton] to concentration for rivers [mg/L]
            # self.ToConc = 10**(6) / (self.RiverWidth * self.WaterLevelR * self.DCL)
            self.ToConc = pcr.cover((10 ** (6) / (self.RiverRunoff * self.timestepsecs)), self.ZeroMap)
            self.SedConc = self.OutSedLoad * self.ToConc
            self.SSConc = (self.OutClayLoad + self.OutSiltLoad) * self.ToConc
            self.BedConc = (
                self.OutSandLoad
                + self.OutSaggLoad
                + self.OutLaggLoad
                + self.OutGravLoad
            ) * self.ToConc
            self.MaxSedConc = self.MaxSedLoad * self.ToConc


        # Area specific runoff (annual runoff / upstream area) [L/km2/s]
        # values between 5 and 25 for the Rhine
        self.SpecRun = (
            self.RiverRunoff * 1000 / pcr.accuflux(self.TopoLdd, self.cellareaKm)
        )


# The main function is used to run the program from the command line


def main(argv=None):
    """
    *Optional but needed it you want to run the model from the command line*
    
    Perform command line execution of the model. This example uses the getopt
    module to parse the command line options.
    
    The user can set the caseName, the runDir, the timestep and the configfile.
    """
    global multpars
    caseName = "default_sediment"
    runId = "run_default"
    configfile = "wflow_sediment.ini"
    _lastTimeStep = 0
    _firstTimeStep = 0
    LogFileName = 'wflow.log'
    
    timestepsecs = 86400
    wflow_cloneMap = "wflow_subcatch.map"

    # This allows us to use the model both on the command line and to call
    # the model using main function from another python script.

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return

    ########################################################################
    ## Process command-line options                                        #
    ########################################################################
    try:
        opts, args = getopt.getopt(argv, "C:S:T:c:s:R:L:hP:p:XIi:")
    except getopt.error as msg:
        pcrut.usage(msg)

    for o, a in opts:
        if o == "-C":
            caseName = a
        if o == "-R":
            runId = a
        if o == "-c":
            configfile = a
        if o == "-L":
            LogFileName = a
        if o == "-s":
            timestepsecs = int(a)
        if o == "-T":
            _lastTimeStep = int(a)
        if o == "-S":
            _firstTimeStep = int(a)
        if o == "-h":
            usage()
    
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
            NoOverWrite=False, 
            level=logging.DEBUG, 
            logfname=LogFileName, 
            model = "wflow_sediment",
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
        if o == "-T":
            configset(myModel.config, "run", "endtime", a, overwrite=True)
        if o == "-S":
            configset(myModel.config, "run", "starttime", a, overwrite=True)
    
    dynModelFw.setupFramework()
    dynModelFw._runInitial()
    dynModelFw._runResume()
    # dynModelFw._runDynamic(0,0)
    dynModelFw._runDynamic(_firstTimeStep, _lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()


if __name__ == "__main__":
    main()
