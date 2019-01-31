#!/usr/bin/python

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

import numpy
import os
import os.path
import shutil, glob
import getopt

from wflow.wf_DynamicFramework import *
from wflow.wflow_adapt import *

# import scipy


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
        self.SaveDir = os.path.join(self.Dir, self.runId)

    def smoothriv(self, river, subcatch, streamorder, param, wlgt):
        """
      Average a parameter (eg. slope, width...) over river cells only and by streamorder

      :param river: River map 
      :param subcatch: Subcatchment map 
      :param streamorder: Streamorder map
      :param param: river parameter to average
      :param wglt: window length in nb of cells

      :return: 
      """
        # Cover river map with zeros
        # Rivercov = ifthenelse(river == 1, 1.0, scalar(0.0))
        # Rivercov = ifthen(subcatch > 0, cover(boolean(river), 0))
        Rivercov = cover(boolean(river), 0)
        Rivercov = scalar(Rivercov)

        order3 = ifthenelse(streamorder == 3, param, scalar(0.0))
        riv3 = ifthenelse(streamorder == 3, Rivercov, scalar(0.0))
        order4 = ifthenelse(streamorder == 4, param, scalar(0.0))
        riv4 = ifthenelse(streamorder == 4, Rivercov, scalar(0.0))
        order5 = ifthenelse(streamorder == 5, param, scalar(0.0))
        riv5 = ifthenelse(streamorder == 5, Rivercov, scalar(0.0))
        order6 = ifthenelse(streamorder == 6, param, scalar(0.0))
        riv6 = ifthenelse(streamorder == 6, Rivercov, scalar(0.0))
        order7 = ifthenelse(streamorder == 7, param, scalar(0.0))
        riv7 = ifthenelse(streamorder == 7, Rivercov, scalar(0.0))
        order8 = ifthenelse(streamorder == 8, param, scalar(0.0))
        riv8 = ifthenelse(streamorder == 8, Rivercov, scalar(0.0))
        order9 = ifthenelse(streamorder == 9, param, scalar(0.0))
        riv9 = ifthenelse(streamorder == 9, Rivercov, scalar(0.0))

        unitmap = riv9 * 0.0 + 1.0

        # self.W = ifthenelse(self.streamorder == 1, windowaverage(order1, celllength()*2) * windowtotal(riv1, celllength()*2), self.W)
        param = ifthenelse(
            streamorder == 9,
            windowaverage(order9, celllength() * wlgt)
            * windowtotal(unitmap, celllength() * wlgt)
            / windowtotal(riv9, celllength() * wlgt),
            param,
        )
        param = ifthenelse(
            streamorder == 8,
            windowaverage(order8, celllength() * wlgt)
            * windowtotal(unitmap, celllength() * wlgt)
            / windowtotal(riv8, celllength() * wlgt),
            param,
        )
        param = ifthenelse(
            streamorder == 7,
            windowaverage(order7, celllength() * wlgt)
            * windowtotal(unitmap, celllength() * wlgt)
            / windowtotal(riv7, celllength() * wlgt),
            param,
        )
        param = ifthenelse(
            streamorder == 6,
            windowaverage(order6, celllength() * wlgt)
            * windowtotal(unitmap, celllength() * wlgt)
            / windowtotal(riv6, celllength() * wlgt),
            param,
        )
        param = ifthenelse(
            streamorder == 5,
            windowaverage(order5, celllength() * wlgt)
            * windowtotal(unitmap, celllength() * wlgt)
            / windowtotal(riv5, celllength() * wlgt),
            param,
        )
        param = ifthenelse(
            streamorder == 4,
            windowaverage(order4, celllength() * wlgt)
            * windowtotal(unitmap, celllength() * wlgt)
            / windowtotal(riv4, celllength() * wlgt),
            param,
        )
        param = ifthenelse(
            streamorder == 3,
            windowaverage(order3, celllength() * wlgt)
            * windowtotal(unitmap, celllength() * wlgt)
            / windowtotal(riv3, celllength() * wlgt),
            param,
        )

        return param

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
        :var SurfaceRunoff: Surface runoff in the kinematic wave [m3/s]
        :var WaterLevel: Water level in the kinematic wave [m]
        
      """

        modelparameters = []

        # Input time series from ini file
        self.P_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Precipitation", "/inmaps/P"
        )  # precipitation
        self.Int_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Interception", "/inmaps/int"
        )  # rainfall interception
        self.SR_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "SurfaceRunoff", "/inmaps/run"
        )  # surface runoff
        self.WL_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "WaterLevel", "/inmaps/levKin"
        )  # water level

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
                name="SurfaceRunoff",
                stack=self.SR_mapstack,
                type="timeseries",
                default=0.0,
                verbose=True,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="WaterLevel",
                stack=self.WL_mapstack,
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
        self.wf_suspend(self.SaveDir + "/outstate/")

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
        setglobaloption("unittrue")

        self.timestepsecs = int(
            configget(self.config, "model", "timestepsecs", "86400")
        )
        sizeinmetres = int(configget(self.config, "layout", "sizeinmetres", "0"))
        self.basetimestep = 86400

        # Reads all parameter from disk
        self.wf_updateparameters()

        """  Read static model parameters/maps from ini file """

        self.intbl = configget(self.config, "model", "intbl", "intbl")
        self.RunRiverModel = int(configget(self.config, "model", "runrivermodel", "1"))
        self.slopecorr = int(configget(self.config, "model", "slopecorr", "0"))
        self.UsleKMethod = int(configget(self.config, "model", "uslekmethod", "2"))
        self.UsleCMethod = int(configget(self.config, "model", "uslecmethod", "2"))
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

        # Static maps to use
        wflow_dem = configget(
            self.config, "model", "wflow_dem", "staticmaps/wflow_dem.map"
        )
        self.Altitude = self.wf_readmap(
            os.path.join(self.Dir, wflow_dem), 0.0, fail=True
        )
        wflow_landuse = configget(
            self.config, "model", "wflow_landuse", "staticmaps/wflow_landuse.map"
        )
        self.LandUse = self.wf_readmap(
            os.path.join(self.Dir, wflow_landuse), 0.0, fail=True
        )
        wflow_soil = configget(
            self.config, "model", "wflow_soil", "staticmaps/wflow_soil.map"
        )
        self.Soil = self.wf_readmap(os.path.join(self.Dir, wflow_soil), 0.0, fail=False)
        wflow_subcatch = configget(
            self.config, "model", "wflow_subcatch", "staticmaps/wflow_subcatch.map"
        )
        self.TopoId = self.wf_readmap(
            os.path.join(self.Dir, wflow_subcatch), 0.0, fail=True
        )
        wflow_Hype = configget(
            self.config, "model", "wflow_Hype", "staticmaps/SUBID-HYPE-Rhine.map"
        )
        self.HypeId = self.wf_readmap(
            os.path.join(self.Dir, wflow_Hype), 0.0, fail=False
        )
        wflow_ldd = configget(
            self.config, "model", "wflow_ldd", "staticmaps/wflow_ldd.map"
        )
        self.TopoLdd = self.wf_readmap(
            os.path.join(self.Dir, wflow_ldd), 0.0, fail=True
        )
        wflow_river = configget(
            self.config, "model", "wflow_river", "staticmaps/wflow_river.map"
        )
        self.River = self.wf_readmap(
            os.path.join(self.Dir, wflow_river), 0.0, fail=True
        )
        wflow_riverwidth = configget(
            self.config, "model", "wflow_riverwidth", "staticmaps/RiverWidth.map"
        )
        self.RiverWidth = self.wf_readmap(
            os.path.join(self.Dir, wflow_riverwidth), 0.0, fail=True
        )
        self.RiverWidth = ifthenelse(self.River == 1, self.RiverWidth, scalar(0.0))
        wflow_dcl = configget(self.config, "model", "wflow_dcl", "staticmaps/DCL.map")
        self.DCL = self.wf_readmap(
            os.path.join(self.Dir, wflow_dcl), 0.0, fail=True
        )  # Drain/River length
        wflow_streamorder = configget(
            self.config,
            "model",
            "wflow_streamorder",
            "staticmaps/wflow_streamorder.map",
        )
        self.streamorder = self.wf_readmap(
            os.path.join(self.Dir, wflow_streamorder), 0.0, fail=True
        )

        # Soil
        wflow_clay = configget(
            self.config, "model", "wflow_clay", "staticmaps/percent_clay.map"
        )
        self.PercentClay = self.wf_readmap(
            os.path.join(self.Dir, wflow_clay), 0.1, fail=True
        )
        wflow_silt = configget(
            self.config, "model", "wflow_silt", "staticmaps/percent_silt.map"
        )
        self.PercentSilt = self.wf_readmap(
            os.path.join(self.Dir, wflow_silt), 0.1, fail=True
        )
        wflow_oc = configget(
            self.config, "model", "wflow_oc", "staticmaps/percent_oc.map"
        )
        self.PercentOC = self.wf_readmap(
            os.path.join(self.Dir, wflow_oc), 0.1, fail=False
        )
        #    wflow_bulk = configget(self.config, "model", "wflow_bulk", "staticmaps/bulk_density.map")
        #    self.BulkDensity = self.wf_readmap(os.path.join(self.Dir,wflow_bulk), 0.1, fail=False)

        """ Read tbl parameters and link them to landuse, soil, subcatch...) """

        # Soil impervious area
        self.PathFrac = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/PathFrac.tbl",
            self.LandUse,
            self.TopoId,
            self.Soil,
            0.01,
        )

        # Soil model parameters
        self.ErosOv = self.readtblDefault(
            self.Dir + "/" + self.intbl + "/eros_ov.tbl",
            self.LandUse,
            self.TopoId,
            self.Soil,
            0.90,
        )

        if self.RunRiverModel == 1:
            # River model parameters
            self.D50River = self.readtblFlexDefault(
                self.Dir + "/" + self.intbl + "/D50_River.tbl", 0.050, wflow_streamorder
            )
            self.D50River = ifthenelse(self.River == 1, self.D50River, scalar(0.0))
            self.CovRiver = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/cov_River.tbl",
                self.LandUse,
                self.TopoId,
                self.Soil,
                1.0,
            )

        """ Determine global variables """
        # Map with zeros
        self.ZeroMap = 0.0 * self.Altitude

        # Subcatchment map
        self.subcatch = ordinal(self.TopoId)
        self.subcatch = ifthen(self.subcatch > 0, self.subcatch)

        # HYPE subcatchment map
        self.HYPEcatch = ordinal(self.HypeId)
        self.HYPEcatch = ifthen(self.HYPEcatch > 0, self.HYPEcatch)

        # Determine real slope, cell length and cell area
        self.xl, self.yl, self.reallength = pcrut.detRealCellLength(
            self.ZeroMap, sizeinmetres
        )
        self.Slope = slope(self.Altitude)
        self.Slope = max(0.00001, self.Slope * celllength() / self.reallength)

        self.cellareaKm = (self.reallength / 1000.0) ** 2
        self.UpArea = accuflux(self.TopoLdd, self.cellareaKm)

        # Sine of the slope
        self.sinSlope = sin(atan(self.Slope))

        # If necessary reduce dem and slope for river cells
        if self.slopecorr == 1:
            self.AltitudeMin = self.wf_readmap(
                os.path.join(self.Dir, "staticmaps/wflow_demmin.map"), 0.0, fail=True
            )
            # In 90m SRTM correct sea cells (-32768) with wflow_dem
            self.AltitudeMin = ifthenelse(
                self.AltitudeMin == -32768, scalar(0.0), self.AltitudeMin
            )
            # Take this minimum value only for river cells (rest is original DEM)
            self.RiverDem = ifthenelse(self.River == 1, self.AltitudeMin, scalar(0.0))
            self.RiverDem = ifthen(self.subcatch > 0, cover(self.RiverDem, scalar(0.0)))
            # Compute corrected slope for river cells
            self.Rivercov = ifthenelse(self.River == 1, 1.0, scalar(0.0))
            self.Rivercov = ifthen(self.subcatch > 0, cover(self.Rivercov, scalar(0.0)))
            self.nrupcell = upstream(self.TopoLdd, self.Rivercov)
            self.RiverSlope = (
                upstream(self.TopoLdd, self.RiverDem) / self.nrupcell
                - downstream(self.TopoLdd, self.RiverDem)
            ) / (2 * self.reallength)
            self.RiverSlope = max(0.0003, self.RiverSlope)

            # Cover the non river cells with original slope/dem
            self.RiverDem = cover(self.RiverDem, self.Altitude)
            self.RiverSlope = cover(self.RiverSlope, self.Slope)

        else:
            self.RiverDem = self.Altitude
            self.RiverSlope = self.Slope

        # Correct slope with drain length
        drainlength = detdrainlength(self.TopoLdd, self.xl, self.yl)
        riverslopecor = drainlength / self.DCL
        self.RiverSlope = self.RiverSlope * riverslopecor

        # Smooth river slope
        self.RiverSlope = self.smoothriv(
            self.River, subcatch, self.streamorder, self.RiverSlope, 6
        )

        #    #Calculate RiverWidth
        #    upstr = catchmenttotal(1, self.TopoLdd)
        #    Qmax = 2200
        #    Qscale = upstr / mapmaximum(upstr) * Qmax
        #    alf = 120
        #    self.N = self.readtblDefault(self.Dir + "/" + self.intbl + "/N.tbl",self.LandUse,self.TopoId,self.Soil,0.072)
        #    self.NRiver = self.readtblFlexDefault(self.Dir + "/" + self.intbl + "/N_River.tbl", 0.036, wflow_streamorder)
        #    self.N = ifthen(subcatch > 0, cover(ifthenelse(self.River, self.NRiver, self.N), self.N))
        #    self.W = (alf * (alf + 2.0) ** (0.6666666667)) ** (0.375)* Qscale ** (0.375) \
        #            * (max(0.0001, windowaverage(self.Slope, celllength() * 4.0))) ** (-0.1875)* self.N ** (0.375)
        #    self.WRiv = ifthen(self.River,((alf * (alf + 2.0) ** (0.6666666667)) ** (0.375)* Qscale ** (0.375)* \
        #                                   self.RiverSlope ** (-0.1875)* self.N ** (0.375)))
        #    #Smooth river width
        #    self.WRiv = self.smoothriv(self.River, subcatch, self.streamorder, self.WRiv, 4)
        #    self.W = ifthen(subcatch > 0, cover(ifthenelse(self.River, self.WRiv, self.W), self.W))
        #    self.RiverWidth = self.W

        """ Determine variables for the soil loss model """

        # Canopy gap fraction based on LAI or input table
        if hasattr(self, "LAI"):
            if not hasattr(self, "Kext"):
                logging.error(
                    "Kext (canopy extinction coefficient) not defined! Needed becausee LAI is defined."
                )
                logging.error("Please add it to the modelparameters section. e.g.:")
                logging.error(
                    "Kext=inmaps/clim/LCtoExtinctionCoefficient.tbl,tbl,0.5,1,inmaps/clim/LC.map"
                )
            self.CanopyGapFraction = exp(-self.Kext * self.LAI)
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

        # Calculate detachability of the soil (k) [g/J]
        self.ErosK = ifthenelse(
            pcrand(
                self.PercentClay >= 40.0,
                pcrand(self.PercentSand >= 20.0, self.PercentSand <= 45.0),
            ),
            2.0,
            ifthenelse(
                pcrand(
                    self.PercentClay >= 27.0,
                    pcrand(self.PercentSand >= 20.0, self.PercentSand <= 45.0),
                ),
                1.7,
                ifthenelse(
                    pcrand(self.PercentSilt <= 40.0, self.PercentSand <= 20.0),
                    2.0,
                    ifthenelse(
                        pcrand(self.PercentSilt > 40.0, self.PercentClay >= 40.0),
                        1.6,
                        ifthenelse(
                            pcrand(self.PercentClay >= 35.0, self.PercentSand >= 45.0),
                            1.9,
                            ifthenelse(
                                pcrand(
                                    self.PercentClay >= 27.0, self.PercentSand < 20.0
                                ),
                                1.6,
                                ifthenelse(
                                    pcrand(
                                        self.PercentClay <= 10.0,
                                        self.PercentSilt >= 80.0,
                                    ),
                                    1.2,
                                    ifthenelse(
                                        self.PercentSilt >= 50,
                                        1.5,
                                        ifthenelse(
                                            pcrand(
                                                self.PercentClay >= 7.0,
                                                pcrand(
                                                    self.PercentSand <= 52.0,
                                                    self.PercentSilt >= 28.0,
                                                ),
                                            ),
                                            2.0,
                                            ifthenelse(
                                                self.PercentClay >= 20.0,
                                                2.1,
                                                ifthenelse(
                                                    self.PercentClay
                                                    >= self.PercentSand - 70.0,
                                                    2.6,
                                                    ifthenelse(
                                                        self.PercentClay
                                                        >= (2.0 * self.PercentSand)
                                                        - 170.0,
                                                        3,
                                                        scalar(1.9),
                                                    ),
                                                ),
                                            ),
                                        ),
                                    ),
                                ),
                            ),
                        ),
                    ),
                ),
            ),
        )

        # Compute USLE K factor
        if self.UsleKMethod == 1:
            self.UsleC = self.wf_readmap(
                os.path.join(self.Dir, "staticmaps/USLE_K.map"), 0.1, fail=True
            )
        if self.UsleKMethod == 2:
            # Calculate USLE K factor (from Renard et al. 1997, with the geometric mean particle diameter Dg)
            self.Dg = exp(
                0.01
                * (
                    self.PercentClay * ln(0.001)
                    + self.PercentSilt * ln(0.025)
                    + self.PercentSand * ln(0.999)
                )
            )  # [mm]
            self.UsleK = 0.0034 + 0.0405 * exp(
                -1 / 2 * ((log10(self.Dg) + 1.659) / 0.7101) ** 2
            )
            # Remove possible outliers
            self.UsleK = max(0.0, self.UsleK)

        if self.UsleKMethod == 3:
            # Calculate USLE K factor (from Williams and Renard 1983, EPIC: a new method for assessing erosion's effect on soil productivity)
            self.SN = (1 - self.PercentSand) / 100
            self.UsleK = (
                (
                    0.2
                    + 0.3
                    * exp(-0.0256 * self.PercentSand * (1 - self.PercentSilt / 100))
                )
                * (self.PercentSilt / (max(0.01, self.PercentClay + self.PercentSilt)))
                ** 0.3
                * (
                    1
                    - (0.25 * self.PercentOC)
                    / (self.PercentOC + exp(3.72 - 2.95 * self.PercentOC))
                )
                * (1 - (0.7 * self.SN) / (self.SN + exp(-5.51 + 22.9 * self.SN)))
            )
            # Remove possible outliers
            self.UsleK = max(0.0, self.UsleK)

        # Compute USLE C factor
        if self.UsleCMethod == 1:
            self.UsleC = self.wf_readmap(
                os.path.join(self.Dir, "staticmaps/USLE_C.map"), 0.1, fail=True
            )
        if self.UsleCMethod == 2:
            # USLE C factor map based on land use (from Gericke 2015, Soil loss estimation and empirical relationships for sediment delivery ratios of European river catchments)
            self.UsleC = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/USLE_C.tbl",
                self.LandUse,
                self.TopoId,
                self.Soil,
                1.0,
            )

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
                self.TopoId,
                self.Soil,
                2.0,
            )
        if self.RainErodMethod == 2:  # ANSWERS
            # Coefficient
            self.ErosSpl = self.readtblDefault(
                self.Dir + "/" + self.intbl + "/eros_spl.tbl",
                self.LandUse,
                self.TopoId,
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
            self.vd50 = ln(
                (1 / ((self.PercentClay + self.PercentSilt) / 100) - 1)
                / (1 / (self.PercentClay / 100) - 1)
            )
            self.wd50 = ln(
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
            self.cd50 = ad50 * ln(self.vd50 / self.wd50)
            self.ud50 = (-self.vd50) ** (1 - bd50) / (-self.wd50) ** (-bd50)

            self.D50 = 1 + (
                -1 / self.ud50 * ln(1 / (1 / (self.PercentClay / 100) - 1))
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
            self.FracSagg = ifthenelse(
                self.PercentClay < 25,
                2.0 * self.PercentClay / 100,
                ifthenelse(
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
            # River particle size distribution (estimated with SWAT method) !!! D50 can be calibrated
            self.FracClayRiv = ifthenelse(
                self.D50River <= 0.005,
                0.65,
                ifthenelse(self.D50River > 2.0, 0.005, scalar(0.15)),
            )
            self.FracSiltRiv = ifthenelse(
                pcrand(self.D50River > 0.005, self.D50River <= 0.050),
                0.65,
                scalar(0.15),
            )
            self.FracSandRiv = ifthenelse(
                pcrand(self.D50River > 0.050, self.D50River <= 2.0), 0.65, scalar(0.15)
            )
            self.FracGravRiv = ifthenelse(self.D50River > 2.0, 0.65, scalar(0.05))

            # Parameters of Bagnold transport formula
            if self.RivTransportMethod == 2:
                self.cBagnold = self.readtblDefault(
                    self.Dir + "/" + self.intbl + "/c_Bagnold.tbl",
                    self.LandUse,
                    self.TopoId,
                    self.Soil,
                    0.0015,
                )
                self.expBagnold = self.readtblDefault(
                    self.Dir + "/" + self.intbl + "/exp_Bagnold.tbl",
                    self.LandUse,
                    self.TopoId,
                    self.Soil,
                    1.4,
                )

            # Parameters of Kodatie transport formula
            if self.RivTransportMethod == 3:
                self.aK = ifthenelse(
                    self.D50River <= 0.05,
                    281.4,
                    ifthenelse(
                        self.D50River <= 0.25,
                        2829.6,
                        ifthenelse(self.D50River <= 2, 2123.4, scalar(431884.8)),
                    ),
                )
                self.bK = ifthenelse(
                    self.D50River <= 0.05,
                    2.622,
                    ifthenelse(
                        self.D50River <= 0.25,
                        3.646,
                        ifthenelse(self.D50River <= 2, 3.300, scalar(1.0)),
                    ),
                )
                self.cK = ifthenelse(
                    self.D50River <= 0.05,
                    0.182,
                    ifthenelse(
                        self.D50River <= 0.25,
                        0.406,
                        ifthenelse(self.D50River <= 2, 0.468, scalar(1.0)),
                    ),
                )
                self.dK = ifthenelse(
                    self.D50River <= 0.05,
                    0.0,
                    ifthenelse(
                        self.D50River <= 0.25,
                        0.412,
                        ifthenelse(self.D50River <= 2, 0.613, scalar(2.0)),
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
                    0.13 * E ** (-0.392) * exp(-0.015 * E ** 2)
                    + 0.045 * (1 - exp(-0.068 * E))
                )
            )
            self.KdBed = 0.2 * self.TCrBed ** (-0.5) * 10 ** (-6)

        """ Variables for reservoirs model """

        if hasattr(self, "ReserVoirSimpleLocs") or hasattr(
            self, "ReserVoirComplexLocs"
        ):
            self.ReserVoirLocs = self.ZeroMap
            self.filter_Eros_TC = self.ZeroMap + 1.0

        if hasattr(self, "ReserVoirSimpleLocs"):
            # Check if we have simple and or complex reservoirs
            tt_simple = pcr2numpy(self.ReserVoirSimpleLocs, 0.0)
            self.nrresSimple = tt_simple.max()
            self.ReserVoirLocs = self.ReserVoirLocs + cover(
                scalar(self.ReserVoirSimpleLocs)
            )
            areamap = self.reallength * self.reallength
            res_area = areatotal(spatial(areamap), self.ReservoirSimpleAreas)

            resarea_pnt = ifthen(boolean(self.ReserVoirSimpleLocs), res_area)
            self.ResSimpleArea = ifthenelse(
                cover(self.ResSimpleArea, scalar(0.0)) > 0,
                self.ResSimpleArea,
                cover(resarea_pnt, scalar(0.0)),
            )
            self.filter_Eros_TC = ifthenelse(
                boolean(cover(res_area, scalar(0.0))),
                res_area * 0.0,
                self.filter_Eros_TC,
            )
        else:
            self.nrresSimple = 0

        if hasattr(self, "ReserVoirComplexLocs"):
            tt_complex = pcr2numpy(self.ReserVoirComplexLocs, 0.0)
            self.nrresComplex = tt_complex.max()
            self.ReserVoirLocs = self.ReserVoirLocs + cover(
                scalar(self.ReserVoirComplexLocs)
            )
            res_area = cover(scalar(self.ReservoirComplexAreas), 0.0)
            self.filter_Eros_TC = ifthenelse(
                res_area > 0, res_area * 0.0, self.filter_Eros_TC
            )
        else:
            self.nrresComplex = 0

        if (self.nrresSimple + self.nrresComplex) > 0:
            self.ReserVoirLocs = ordinal(self.ReserVoirLocs)
            self.logger.info(
                "A total of "
                + str(self.nrresSimple)
                + " simple reservoirs and "
                + str(self.nrresComplex)
                + " complex reservoirs found."
            )
            self.ReserVoirDownstreamLocs = downstream(self.TopoLdd, self.ReserVoirLocs)
            self.TopoLddOrg = lddrepair(
                cover(ifthen(boolean(self.ReserVoirLocs), ldd(5)), self.TopoLdd)
            )

            tt_filter = pcr2numpy(self.filter_Eros_TC, 1.0)
            self.filterResArea = tt_filter.min()

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
            "self.D50River",
            "self.FracClayRiv",
            "self.filter_Eros_TC",
            "self.RiverSlope",
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
      :var self.HYPESedCatch: Land eroded sediment reaching the river summed per HYPE subcatchment [ton]
      
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
        self.Precipitation = max(0.0, self.Precipitation)

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
            self.KeDirect = max(
                8.95 + 8.44 * log10(max(0.0001, self.rintnsty)), 0.0
            )  # variant, most used in distributed models

            # Kinetic energy of leaf drainage [J/m^2/mm]
            pheff = 0.5 * self.CanopyHeight  # [m]
            self.KeLeaf = max((15.8 * pheff ** 0.5) - 5.87, 0.0)

            # Depths of rainfall (total, leaf drainage, direct) [mm]
            # rdepth_tot = max(self.Precipitation/self.timestepsecs, 0.0)
            rDepthTot = max(self.Precipitation, 0.0)
            rDepthLeaf = max(rDepthTot * 0.1 * self.CanopyGapFraction, 0.0)  # stemflow
            rDepthDirect = max(
                rDepthTot - rDepthLeaf - self.Interception, 0.0
            )  # throughfall

            # Total kinetic energy by rainfall [J/m^2]
            self.KeTotal = (
                rDepthDirect * self.KeDirect + rDepthLeaf * self.KeLeaf
            ) * 0.001

            # Rainfall/Splash erosion
            self.SedSpl = (
                self.ErosK * self.KeTotal * exp(-self.ErosSpl * self.WaterLevel)
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
        self.SedSpl = ifthen(self.subcatch > 0, cover(self.SedSpl, 0.0))

        """ Overland flow erosion from ANSWERS"""
        # Only calculate overland flow erosion outside of river cells
        self.OvRun = cover(
            ifthenelse(self.River == 1, 0.0, self.SurfaceRunoff), self.SurfaceRunoff
        )  # [m3/s]
        self.OvRunRate = self.OvRun * 60 / self.reallength  # [m2/min]

        # Overland flow erosion
        # For a wide range of slope, it is better to use the sine of slope rather than tangeant
        self.SedOv = (
            self.ErosOv
            * self.UsleC
            * self.UsleK
            * self.cellareaKm
            * 10 ** 6
            * self.sinSlope
            * self.OvRunRate
        )  # [kg/min]
        self.SedOv = (
            self.timestepsecs / 60.0 * 10 ** (-3) * self.SedOv
        )  # [ton/timestep]

        # Remove the impervious areas
        self.SedOv = self.SedOv * (1.0 - self.PathFrac)
        # Remove nodata values
        self.SedOv = ifthen(self.subcatch > 0, cover(self.SedOv, 0.0))

        """ Total soil detachment """
        self.SoilLoss = self.SedSpl + self.SedOv  # [ton/cell/timestep]

        # Remove land erosion for reservoir cells
        if (self.nrresSimple + self.nrresComplex) > 0 and self.filterResArea == 0:
            self.SedSpl = self.filter_Eros_TC * self.SedSpl
            self.SedOv = self.filter_Eros_TC * self.SedOv
            self.SoilLoss = self.filter_Eros_TC * self.SoilLoss

        ##########################################################################
        # Inland sediment routing model #############################
        ##########################################################################

        # If the river model is run, use Yalin's equation with particle differentiation.
        # If just the soil loss model is run, use Govers equation without particle differentiation

        """ Transport of overland flow with no particle differenciation using Govers equation"""

        if self.LandTransportMethod == 2 or self.LandTransportMethod == 3:
            if self.LandTransportMethod == 2:
                # Unit stream power
                self.velocity = cover(
                    ifthenelse(
                        self.WaterLevel > 0,
                        self.OvRun / (self.reallength * self.WaterLevel),
                        0.0,
                    ),
                    0.0,
                )  # [m/s]
                self.omega = 10 * self.sinSlope * 100 * self.velocity  # [cm/s]
                # self.omega = self.sinSlope * 100 * self.velocity #[cm/s]

                # Transport capacity from Govers, 1990
                self.TCf = ifthenelse(
                    self.omega > 0.4,
                    self.cGovers * (self.omega - 0.4) ** self.nGovers * 2650,
                    0.0,
                )  # [kg/m3]
                # self.TC = self.TCf / (1 - self.TCf/2650) #[kg/m3]
                # self.TC = max(self.TC, 2650)
                self.TC = (
                    self.TCf * self.OvRun * self.timestepsecs * 10 ** (-3)
                )  # [ton/cell/timestep]
                # Remove nodata values
                self.TC = ifthen(self.subcatch > 0, cover(self.TC, 0.0))
                # Assume that eroded soil on lake cells all reach the river cells of the reservoir
                if (
                    self.nrresSimple + self.nrresComplex
                ) > 0 and self.filterResArea == 0:
                    self.TC = ifthenelse(
                        pcrand(self.filter_Eros_TC == 0, self.River == 0),
                        10 ** 9,
                        self.TC,
                    )

            elif self.LandTransportMethod == 3:
                # Transport capacity from Yalin
                self.OvLevel = cover(
                    ifthenelse(self.River == 1, 0.0, self.WaterLevel), self.WaterLevel
                )  # [m]
                self.delta = max(
                    self.OvLevel
                    * self.sinSlope
                    / (self.D50 * 10 ** (-3) * (2650 / 1000 - 1))
                    / 0.06
                    - 1,
                    0.0,
                )
                self.TC = (
                    self.reallength
                    / self.OvRun
                    * (2650 - 1000)
                    * self.D50
                    * 10 ** (-3)
                    * (9.81 * self.OvLevel * self.sinSlope)
                    * 0.635
                    * self.delta
                    * (
                        1
                        - ln(1 + self.delta * 2.45 / (2650 / 1000) ** 0.4 * 0.06 ** 0.5)
                        / self.delta
                        * 2.45
                        / (2650 / 1000) ** 0.4
                        * 0.06 ** 0.5
                    )
                )  # [kg/m3]
                self.TC = ifthen(
                    self.subcatch > 0,
                    cover(self.TC * self.OvRun * self.timestepsecs * 10 ** (-3), 0.0),
                )  # [ton/cell/timestep]
                # Assume that eroded soil on lake cells all reach the river cells of the reservoir
                if (
                    self.nrresSimple + self.nrresComplex
                ) > 0 and self.filterResArea == 0:
                    self.TC = ifthenelse(
                        pcrand(self.filter_Eros_TC == 0, self.River == 0),
                        10 ** 9,
                        self.TC,
                    )

            # To get total sediment input from land into the river systems, river cells transport all sediment to the output (huge TC)
            self.TCRiv = cover(ifthenelse(self.River == 1, 10 ** 9, self.TC), self.TC)
            # Transported sediment over the land
            self.SedFlux = accucapacityflux(
                self.TopoLdd, self.SoilLoss, self.TCRiv
            )  # [ton/cell/tinestep]
            # Deposited sediment over the land
            self.SedDep = accucapacitystate(self.TopoLdd, self.SoilLoss, self.TCRiv)

            # Sediment amount reaching each river cell '''
            self.SedDep2 = accucapacitystate(self.TopoLdd, self.SoilLoss, self.TC)
            # Remove inland deposition
            self.OvSed = cover(ifthenelse(self.River == 1, self.SedDep2, 0.0), 0.0)

            # Sum the results per subcatchment (input for D-WAQ) [kg/ha/timestep]
            self.HYPEOvSedCatch = areatotal(
                self.OvSed, self.HYPEcatch
            )  # [ton/cell/timestep]
            self.HYPEOvSedCatch = (
                self.HYPEOvSedCatch * 1000 / (self.cellareaKm * 100)
            )  # [kg/ha/timestep]

            """ Transport of overland flow with particle differenciation using Yalin equation"""
        else:
            # Determine the eroded amount of clay/silt/sand/aggregates on each cell [ton/cell/timestep]
            self.LandErodClay = self.SoilLoss * self.FracClay
            self.LandErodSilt = self.SoilLoss * self.FracSilt
            self.LandErodSand = self.SoilLoss * self.FracSand
            self.LandErodSagg = self.SoilLoss * self.FracSagg
            self.LandErodLagg = self.SoilLoss * self.FracLagg

            # Water level of overland flow
            self.OvLevel = cover(
                ifthenelse(self.River == 1, 0.0, self.WaterLevel), self.WaterLevel
            )  # [m]

            # Delta parameter of Yalin for each particle class
            self.DClay = max(
                self.OvLevel
                * self.sinSlope
                / (2 * 10 ** (-6) * (2650 / 1000 - 1))
                / 0.06
                - 1,
                0.0,
            )
            self.DSilt = max(
                self.OvLevel
                * self.sinSlope
                / (10 * 10 ** (-6) * (2650 / 1000 - 1))
                / 0.06
                - 1,
                0.0,
            )
            self.DSand = max(
                self.OvLevel
                * self.sinSlope
                / (200 * 10 ** (-6) * (2650 / 1000 - 1))
                / 0.06
                - 1,
                0.0,
            )
            self.DSagg = max(
                self.OvLevel
                * self.sinSlope
                / (30 * 10 ** (-6) * (2650 / 1000 - 1))
                / 0.06
                - 1,
                0.0,
            )
            self.DLagg = max(
                self.OvLevel
                * self.sinSlope
                / (500 * 10 ** (-6) * (2650 / 1000 - 1))
                / 0.06
                - 1,
                0.0,
            )

            # Total transportability
            self.Dtot = self.DClay + self.DSilt + self.DSand + self.DSagg + self.DLagg

            # Yalin Transport capacity of overland flow for each particle class
            self.TCClay = (
                self.reallength
                / self.OvRun
                * (2650 - 1000)
                * 2
                * 10 ** (-6)
                * (9.81 * self.OvLevel * self.sinSlope)
                * self.DClay
                / self.Dtot
                * 0.635
                * self.DClay
                * (
                    1
                    - ln(1 + self.DClay * 2.45 / (2650 / 1000) ** 0.4 * 0.06 ** 0.5)
                    / self.DClay
                    * 2.45
                    / (2650 / 1000) ** 0.4
                    * 0.06 ** 0.5
                )
            )  # [kg/m3]
            self.TCClay = ifthen(
                self.subcatch > 0,
                cover(self.TCClay * self.OvRun * self.timestepsecs * 10 ** (-3), 0.0),
            )  # [ton/cell/timestep]

            self.TCSilt = (
                self.reallength
                / self.OvRun
                * (2650 - 1000)
                * 10
                * 10 ** (-6)
                * (9.81 * self.OvLevel * self.sinSlope)
                * self.DSilt
                / self.Dtot
                * 0.635
                * self.DSilt
                * (
                    1
                    - ln(1 + self.DSilt * 2.45 / (2650 / 1000) ** 0.4 * 0.06 ** 0.5)
                    / self.DSilt
                    * 2.45
                    / (2650 / 1000) ** 0.4
                    * 0.06 ** 0.5
                )
            )  # [kg/m3]
            self.TCSilt = ifthen(
                self.subcatch > 0,
                cover(self.TCSilt * self.OvRun * self.timestepsecs * 10 ** (-3), 0.0),
            )  # [ton/cell/timestep]

            self.TCSand = (
                self.reallength
                / self.OvRun
                * (2650 - 1000)
                * 200
                * 10 ** (-6)
                * (9.81 * self.OvLevel * self.sinSlope)
                * self.DSand
                / self.Dtot
                * 0.635
                * self.DSand
                * (
                    1
                    - ln(1 + self.DSand * 2.45 / (2650 / 1000) ** 0.4 * 0.06 ** 0.5)
                    / self.DSand
                    * 2.45
                    / (2650 / 1000) ** 0.4
                    * 0.06 ** 0.5
                )
            )  # [kg/m3]
            self.TCSand = ifthen(
                self.subcatch > 0,
                cover(self.TCSand * self.OvRun * self.timestepsecs * 10 ** (-3), 0.0),
            )  # [ton/cell/timestep]

            self.TCSagg = (
                self.reallength
                / self.OvRun
                * (2650 - 1000)
                * 30
                * 10 ** (-6)
                * (9.81 * self.OvLevel * self.sinSlope)
                * self.DSagg
                / self.Dtot
                * 0.635
                * self.DSagg
                * (
                    1
                    - ln(1 + self.DSagg * 2.45 / (2650 / 1000) ** 0.4 * 0.06 ** 0.5)
                    / self.DSagg
                    * 2.45
                    / (2650 / 1000) ** 0.4
                    * 0.06 ** 0.5
                )
            )  # [kg/m3]
            self.TCSagg = ifthen(
                self.subcatch > 0,
                cover(self.TCSagg * self.OvRun * self.timestepsecs * 10 ** (-3), 0.0),
            )  # [ton/cell/timestep]

            self.TCLagg = (
                self.reallength
                / self.OvRun
                * (2650 - 1000)
                * 500
                * 10 ** (-6)
                * (9.81 * self.OvLevel * self.sinSlope)
                * self.DLagg
                / self.Dtot
                * 0.635
                * self.DLagg
                * (
                    1
                    - ln(1 + self.DLagg * 2.45 / (2650 / 1000) ** 0.4 * 0.06 ** 0.5)
                    / self.DLagg
                    * 2.45
                    / (2650 / 1000) ** 0.4
                    * 0.06 ** 0.5
                )
            )  # [kg/m3]
            self.TCLagg = ifthen(
                self.subcatch > 0,
                cover(self.TCLagg * self.OvRun * self.timestepsecs * 10 ** (-3), 0.0),
            )  # [ton/cell/timestep]

            # Assume that eroded soil on lake cells all reach the river cells of the reservoir
            if (self.nrresSimple + self.nrresComplex) > 0:
                self.TCClay = cover(
                    ifthenelse(self.filter_Eros_TC == 0, 10 ** 9, self.TCClay),
                    self.TCClay,
                )
                self.TCClay = cover(
                    ifthenelse(self.River == 1, 0.0, self.TCClay), self.TCClay
                )
                self.TCSilt = cover(
                    ifthenelse(self.filter_Eros_TC == 0, 10 ** 9, self.TCSilt),
                    self.TCSilt,
                )
                self.TCSilt = cover(
                    ifthenelse(self.River == 1, 0.0, self.TCSilt), self.TCSilt
                )
                self.TCSand = cover(
                    ifthenelse(self.filter_Eros_TC == 0, 10 ** 9, self.TCSand),
                    self.TCSand,
                )
                self.TCSand = cover(
                    ifthenelse(self.River == 1, 0.0, self.TCSand), self.TCSand
                )
                self.TCSagg = cover(
                    ifthenelse(self.filter_Eros_TC == 0, 10 ** 9, self.TCSagg),
                    self.TCSagg,
                )
                self.TCSagg = cover(
                    ifthenelse(self.River == 1, 0.0, self.TCSagg), self.TCSagg
                )
                self.TCLagg = cover(
                    ifthenelse(self.filter_Eros_TC == 0, 10 ** 9, self.TCLagg),
                    self.TCLagg,
                )
                self.TCLagg = cover(
                    ifthenelse(self.River == 1, 0.0, self.TCLagg), self.TCLagg
                )

            # Eroded sediment in overland flow reaching the river system per particle class [ton/cell/timestep]
            self.InLandClay = cover(
                ifthenelse(
                    self.River == 1,
                    accucapacitystate(self.TopoLdd, self.LandErodClay, self.TCClay),
                    0.0,
                ),
                0.0,
            )
            self.InLandSilt = cover(
                ifthenelse(
                    self.River == 1,
                    accucapacitystate(self.TopoLdd, self.LandErodSilt, self.TCSilt),
                    0.0,
                ),
                0.0,
            )
            self.InLandSand = cover(
                ifthenelse(
                    self.River == 1,
                    accucapacitystate(self.TopoLdd, self.LandErodSand, self.TCSand),
                    0.0,
                ),
                0.0,
            )
            self.InLandSagg = cover(
                ifthenelse(
                    self.River == 1,
                    accucapacitystate(self.TopoLdd, self.LandErodSagg, self.TCSagg),
                    0.0,
                ),
                0.0,
            )
            self.InLandLagg = cover(
                ifthenelse(
                    self.River == 1,
                    accucapacitystate(self.TopoLdd, self.LandErodLagg, self.TCLagg),
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

            # Sediment input reaching the river system per particle class and per HYPE subcatchment [kg/ha/timestep]
            self.ClayCatch = (
                areatotal(self.InLandClay, self.HYPEcatch)
                * 1000
                / (self.cellareaKm * 100)
            )
            self.SiltCatch = (
                areatotal(self.InLandSilt, self.HYPEcatch)
                * 1000
                / (self.cellareaKm * 100)
            )
            self.SandCatch = (
                areatotal(self.InLandSand, self.HYPEcatch)
                * 1000
                / (self.cellareaKm * 100)
            )
            self.SaggCatch = (
                areatotal(self.InLandSagg, self.HYPEcatch)
                * 1000
                / (self.cellareaKm * 100)
            )
            self.LaggCatch = (
                areatotal(self.InLandLagg, self.HYPEcatch)
                * 1000
                / (self.cellareaKm * 100)
            )
            self.HYPESedCatch = (
                self.ClayCatch
                + self.SiltCatch
                + self.SandCatch
                + self.SaggCatch
                + self.LaggCatch
            )  # total

        ##########################################################################
        # River transport and processes #############################
        ##########################################################################

        if self.RunRiverModel == 1:
            # River sediment loads are separated into different particle class.
            # Clay, silt and sand can both come from land, resuspension or river channel erosion.
            # Small and large aggregates only come from land erosion or resuspension.
            # Gravel only comes from resuspension or river channel erosion.

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
                + upstream(self.TopoLdd, self.OutSedLoad)
                + self.InLandSed
            )
            self.InClayLoad = (
                self.OldClayLoad
                + upstream(self.TopoLdd, self.OutClayLoad)
                + self.InLandClay
            )
            self.InSiltLoad = (
                self.OldSiltLoad
                + upstream(self.TopoLdd, self.OutSiltLoad)
                + self.InLandSilt
            )
            self.InSandLoad = (
                self.OldSandLoad
                + upstream(self.TopoLdd, self.OutSandLoad)
                + self.InLandSand
            )
            self.InSaggLoad = (
                self.OldSaggLoad
                + upstream(self.TopoLdd, self.OutSaggLoad)
                + self.InLandSagg
            )
            self.InLaggLoad = (
                self.OldLaggLoad
                + upstream(self.TopoLdd, self.OutLaggLoad)
                + self.InLandLagg
            )
            self.InGravLoad = self.OldGravLoad + upstream(
                self.TopoLdd, self.OutGravLoad
            )

            """ River erosion """
            # Hydraulic radius of the river [m] (rectangular channel)
            self.HydRad = (
                self.WaterLevel
                * self.RiverWidth
                / (self.RiverWidth + 2 * self.WaterLevel)
            )

            # Transport capacity
            # Engelund and Hansen transport formula
            if self.RivTransportMethod == 1:
                # Shields parameter
                self.ThetaShields = (
                    1000
                    * 9.81
                    * self.HydRad
                    * self.RiverSlope
                    / ((2650 - 1000) * 9.81 * self.D50River / 1000)
                )
                self.vmean = ifthenelse(
                    pcrand(self.WaterLevel > 0, self.River == 1),
                    self.SurfaceRunoff / (self.RiverWidth * self.WaterLevel),
                    scalar(0.0),
                )
                self.vshear = (9.81 * self.HydRad * self.RiverSlope) ** (0.5)
                #          self.Cw = 0.05 * (2650/(2650-1000)) * self.vmean * self.RiverSlope \
                #              / ((2650-1000)/1000*9.81*self.D50River/1000)**(0.5) * self.ThetaShields**(0.5) # sediment concentration by weight
                self.Cw = min(
                    1.0,
                    ifthenelse(
                        pcrand(self.HydRad > 0, self.River == 1),
                        2.65
                        * 0.05
                        * self.vmean
                        * self.vshear ** 3
                        / ((2.65 - 1) ** 2 * 9.81 ** 2 * self.D50River * self.HydRad),
                        scalar(0.0),
                    ),
                )  # concentration by weight
                self.MaxSedLoad = max(
                    0.0, self.Cw / (self.Cw + (1 - self.Cw) * 2.65) * 2.65
                )  # [tons/m3]
                self.MaxSedLoad = self.MaxSedLoad * (
                    self.WaterLevel * self.RiverWidth * self.DCL
                    + self.SurfaceRunoff * self.timestepsecs
                )  # [ton]

            # Simplified Bagnold transport formula
            if self.RivTransportMethod == 2:
                self.MaxSedLoad = (
                    self.cBagnold
                    * (self.SurfaceRunoff / (self.WaterLevel * self.RiverWidth))
                    ** self.expBagnold
                )  # [ton/m3]
                self.MaxSedLoad = self.MaxSedLoad * (
                    self.WaterLevel * self.RiverWidth * self.DCL
                    + self.SurfaceRunoff * self.timestepsecs
                )  # [ton]

            # Kodatie transport formula
            if self.RivTransportMethod == 3:
                self.vmean = ifthenelse(
                    pcrand(self.WaterLevel > 0, self.River == 1),
                    self.SurfaceRunoff / (self.RiverWidth * self.WaterLevel),
                    scalar(0.0),
                )
                self.MaxSedLoad = (
                    self.aK
                    * self.vmean ** (self.bK)
                    * self.WaterLevel ** (self.cK)
                    * self.RiverSlope ** (self.dK)
                ) * (
                    self.RiverWidth
                )  # [tons]

            # Yang transport formula
            if self.RivTransportMethod == 4:
                self.wsRiv = 411 * self.D50River ** 2 / 3600
                self.vshear = (9.81 * self.HydRad * self.RiverSlope) ** (0.5)
                self.var1 = self.vshear * self.D50River / 1000 / (1.16 * 10 ** (-6))
                self.var2 = self.wsRiv * self.D50River / 1000 / (1.16 * 10 ** (-6))
                self.vcr = ifthenelse(
                    self.var1 >= 70,
                    2.05 * self.wsRiv,
                    self.wsRiv * (2.5 / (log10(self.var1) - 0.06) + 0.66),
                )

                # Sand equation
                self.logCppm = (
                    5.435
                    - 0.286 * log10(self.var2)
                    - 0.457 * log10(self.vshear / self.wsRiv)
                    + (
                        1.799
                        - 0.409 * log10(self.var2)
                        - 0.314 * log10(self.vshear / self.wsRiv)
                    )
                    * log10(
                        (
                            self.SurfaceRunoff / (self.RiverWidth * self.WaterLevel)
                            - self.vcr
                        )
                        * self.RiverSlope
                        / self.wsRiv
                    )
                )
                # Gravel equation
                self.logCppm = ifthenelse(
                    self.D50River < 2.0,
                    self.logCppm,
                    6.681
                    - 0.633 * log10(self.var2)
                    - 4.816 * log10(self.vshear / self.wsRiv)
                    + (
                        2.784
                        - 0.305 * log10(self.var2)
                        - 0.282 * log10(self.vshear / self.wsRiv)
                    )
                    * log10(
                        (
                            self.SurfaceRunoff / (self.RiverWidth * self.WaterLevel)
                            - self.vcr
                        )
                        * self.RiverSlope
                        / self.wsRiv
                    ),
                )
                self.Cw = 10 ** self.logCppm * 10 ** (
                    -6
                )  # sediment concentration by weight
                self.MaxSedLoad = max(
                    0.0, self.Cw / (self.Cw + (1 - self.Cw) * 2.65) * 2.65
                )  # [tons/m3]
                self.MaxSedLoad = self.MaxSedLoad * (
                    self.WaterLevel * self.RiverWidth * self.DCL
                    + self.SurfaceRunoff * self.timestepsecs
                )  # [ton]

            # Molinas & Wu transport formula
            if self.RivTransportMethod == 5:
                self.wsRiv = 411 * self.D50River ** 2 / 3600
                self.psi = (
                    (2.65 - 1)
                    * 9.81
                    * self.WaterLevel
                    * self.wsRiv
                    * (log10(self.WaterLevel / self.D50River)) ** 2
                ) ** (0.5)
                self.Cw = (
                    1430
                    * (0.86 + (self.psi) ** (0.5))
                    * (self.psi) ** (1.5)
                    / (0.016 + self.psi)
                    * 10 ** (-6)
                )  # weight
                self.MaxSedLoad = max(
                    0.0, self.Cw / (self.Cw + (1 - self.Cw) * 2.65) * 2.65
                )  # [tons/m3]
                self.MaxSedLoad = self.MaxSedLoad * (
                    self.WaterLevel * self.RiverWidth * self.DCL
                    + self.SurfaceRunoff * self.timestepsecs
                )  # [ton]

            self.MaxSedLoad = ifthen(self.subcatch > 0, cover(self.MaxSedLoad, 0.0))

            # Repartition of the effective shear stress between the bank and the bed from Knight et al. 1984
            self.SFBank = ifthenelse(
                pcrand(self.WaterLevel > 0, self.River == 1),
                exp(-3.230 * log10(self.RiverWidth / self.WaterLevel + 3) + 6.146),
                scalar(0.0),
            )  # [%]
            # Effective shear stress on river bed and banks [N/m2]
            self.TEffBank = ifthenelse(
                pcrand(self.WaterLevel > 0, self.River == 1),
                1000
                * 9.81
                * self.HydRad
                * self.RiverSlope
                * self.SFBank
                / 100
                * (1 + self.RiverWidth / (2 * self.WaterLevel)),
                scalar(0.0),
            )
            self.TEffBed = (
                1000
                * 9.81
                * self.HydRad
                * self.RiverSlope
                * (1 - self.SFBank / 100)
                * (1 + 2 * self.WaterLevel / self.RiverWidth)
            )
            # Potential erosion rates of the bed and bank [t/cell/timestep] (assuming only one bank is eroding)
            self.ERBank = max(
                0.0,
                self.KdBank
                * (self.TEffBank - self.TCrBank)
                * (self.DCL * self.WaterLevel)
                * 1.4
                * self.timestepsecs,
            )  # 1.4 is bank default bulk density
            self.ERBed = max(
                0.0,
                self.KdBed
                * (self.TEffBed - self.TCrBed)
                * (self.DCL * self.RiverWidth)
                * 1.5
                * self.timestepsecs,
            )  # 1.5 is bed default bulk density
            # Relative potential erosion rates of the bed and bank [-]
            self.RTEBank = ifthenelse(
                self.ERBank + self.ERBed > 0.0,
                self.ERBank / (self.ERBank + self.ERBed),
                0.0,
            )
            self.RTEBed = 1.0 - self.RTEBank

            # Excess transport capacity [ton/cell/timestep]
            # Erosion only if the load is below the transport capacity of the flow
            self.SedEx = max(self.MaxSedLoad - self.InSedLoad, 0.0)
            # Bed and bank are eroded after the previously deposited material
            self.EffSedEx = max(self.SedEx - self.RivStoreSed, 0.0)

            # Bank erosion [ton/cell/timestep]
            self.BankSedLoad = ifthenelse(
                self.EffSedEx == 0,
                0.0,
                ifthenelse(
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
            self.BedSedLoad = ifthenelse(
                self.EffSedEx == 0,
                0.0,
                ifthenelse(
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
            self.DegStoreClay = ifthenelse(
                self.RivStoreClay >= self.SedEx, self.SedEx, self.RivStoreClay
            )
            self.RivStoreClay = self.RivStoreClay - self.DegStoreClay  # update store
            self.SedEx = max(
                self.SedEx - self.DegStoreClay, 0.0
            )  # update amount of sediment that need to be degraded

            self.DegStoreSilt = ifthenelse(
                self.RivStoreSilt >= self.SedEx, self.SedEx, self.RivStoreSilt
            )
            self.RivStoreSilt = self.RivStoreSilt - self.DegStoreSilt
            self.SedEx = max(self.SedEx - self.DegStoreSilt, 0.0)

            self.DegStoreSagg = ifthenelse(
                self.RivStoreSagg >= self.SedEx, self.SedEx, self.RivStoreSagg
            )
            self.RivStoreSagg = self.RivStoreSagg - self.DegStoreSagg
            self.SedEx = max(self.SedEx - self.DegStoreSagg, 0.0)

            self.DegStoreSand = ifthenelse(
                self.RivStoreSand >= self.SedEx, self.SedEx, self.RivStoreSand
            )
            self.RivStoreSand = self.RivStoreSand - self.DegStoreSand
            self.SedEx = max(self.SedEx - self.DegStoreSand, 0.0)

            self.DegStoreLagg = ifthenelse(
                self.RivStoreLagg >= self.SedEx, self.SedEx, self.RivStoreLagg
            )
            self.RivStoreLagg = self.RivStoreLagg - self.DegStoreLagg
            self.SedEx = max(self.SedEx - self.DegStoreLagg, 0.0)

            self.DegStoreGrav = ifthenelse(
                self.RivStoreGrav >= self.SedEx, self.SedEx, self.RivStoreGrav
            )
            self.RivStoreGrav = self.RivStoreGrav - self.DegStoreGrav
            self.SedEx = max(self.SedEx - self.DegStoreGrav, 0.0)

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
            if (self.nrresSimple + self.nrresComplex) > 0 and self.filterResArea == 0:
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
            self.FracDepEx = ifthenelse(
                self.TCEx < 0, 1 - self.MaxSedLoad / self.InSedLoad, scalar(0.0)
            )

            # Fractions of deposited particles in river cells from the Einstein formula [-]
            # Particle fall velocity [m/s] from Stokes
            self.x = ifthenelse(
                pcrand(self.SurfaceRunoff > 0, self.River == 1),
                1.055 * self.DCL / (self.SurfaceRunoff / self.RiverWidth),
                scalar(0.0),
            )
            self.x = ifthen(self.subcatch > 0, cover(self.x, 0.0))

            self.FracDepClay = ifthenelse(
                self.TCEx >= 0,
                min(1.0, 1 - 1 / exp(self.x * (411 * 0.002 ** 2 / 3600))),
                self.FracDepEx,
            )
            self.FracDepSilt = ifthenelse(
                self.TCEx >= 0,
                min(1.0, 1 - 1 / exp(self.x * (411 * 0.010 ** 2 / 3600))),
                self.FracDepEx,
            )
            self.FracDepSand = ifthenelse(
                self.TCEx >= 0,
                min(1.0, 1 - 1 / exp(self.x * (411 * 0.200 ** 2 / 3600))),
                self.FracDepEx,
            )
            self.FracDepSagg = ifthenelse(
                self.TCEx >= 0,
                min(1.0, 1 - 1 / exp(self.x * (411 * 0.030 ** 2 / 3600))),
                self.FracDepEx,
            )
            self.FracDepLagg = ifthenelse(
                self.TCEx >= 0,
                min(1.0, 1 - 1 / exp(self.x * (411 * 0.500 ** 2 / 3600))),
                self.FracDepEx,
            )
            self.FracDepGrav = ifthenelse(
                self.TCEx >= 0,
                min(1.0, 1 - 1 / exp(self.x * (411 * 2.000 ** 2 / 3600))),
                self.FracDepEx,
            )

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
            if (self.nrresSimple + self.nrresComplex) > 0:
                self.DepClay = self.filter_Eros_TC * self.DepClay
                self.DepSilt = self.filter_Eros_TC * self.DepSilt
                self.DepSand = self.filter_Eros_TC * self.DepSand
                self.DepSagg = self.filter_Eros_TC * self.DepSagg
                self.DepLagg = self.filter_Eros_TC * self.DepLagg
                self.DepGrav = self.filter_Eros_TC * self.DepGrav
                self.DepSedLoad = self.filter_Eros_TC * self.DepSedLoad

            # Deposition in reservoirs from Camp 1945
            if (self.nrresSimple + self.nrresComplex) > 0 and self.filterResArea == 0:
                self.VcRes = ifthenelse(
                    self.ReserVoirLocs > 0, self.SurfaceRunoff / self.ResArea, 0.0
                )
                self.DepClay = cover(
                    ifthen(
                        self.ReserVoirLocs > 0,
                        self.InClayLoad
                        * min(1.0, (411 / 3600 * 0.002 ** 2) / self.VcRes),
                    ),
                    self.DepClay,
                )
                self.DepSilt = cover(
                    ifthen(
                        self.ReserVoirLocs > 0,
                        self.InSiltLoad
                        * min(1.0, (411 / 3600 * 0.010 ** 2) / self.VcRes),
                    ),
                    self.DepSilt,
                )
                self.DepSand = cover(
                    ifthen(
                        self.ReserVoirLocs > 0,
                        self.InSandLoad
                        * min(1.0, (411 / 3600 * 0.200 ** 2) / self.VcRes),
                    ),
                    self.DepSand,
                )
                self.DepSagg = cover(
                    ifthen(
                        self.ReserVoirLocs > 0,
                        self.InSaggLoad
                        * min(1.0, (411 / 3600 * 0.030 ** 2) / self.VcRes),
                    ),
                    self.DepSagg,
                )
                self.DepLagg = cover(
                    ifthen(
                        self.ReserVoirLocs > 0,
                        self.InLaggLoad
                        * min(1.0, (411 / 3600 * 0.500 ** 2) / self.VcRes),
                    ),
                    self.DepLagg,
                )
                self.DepGrav = cover(
                    ifthen(
                        self.ReserVoirLocs > 0,
                        self.InGravLoad
                        * min(1.0, (411 / 3600 * 2.000 ** 2) / self.VcRes),
                    ),
                    self.DepGrav,
                )
                self.DepSedLoad = cover(
                    ifthen(
                        self.ReserVoirLocs > 0,
                        self.DepClay
                        + self.DepSilt
                        + self.DepSand
                        + self.DepSagg
                        + self.DepLagg
                        + self.DepGrav,
                    ),
                    self.DepSedLoad,
                )

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
            self.OutFracWat = min(
                self.SurfaceRunoff
                * self.timestepsecs
                / (
                    self.WaterLevel * self.RiverWidth * self.DCL
                    + self.SurfaceRunoff * self.timestepsecs
                ),
                0.999999,
            )

            self.OutSedLoad = ifthen(
                self.subcatch > 0,
                cover(
                    max(
                        0.0,
                        (self.InSedLoad + self.RivErodSed - self.DepSedLoad)
                        * self.OutFracWat,
                    ),
                    0.0,
                ),
            )
            self.OutClayLoad = ifthen(
                self.subcatch > 0,
                cover(
                    max(
                        0.0,
                        (self.InClayLoad + self.RivErodClay - self.DepClay)
                        * self.OutFracWat,
                    ),
                    0.0,
                ),
            )
            self.OutSiltLoad = ifthen(
                self.subcatch > 0,
                cover(
                    max(
                        0.0,
                        (self.InSiltLoad + self.RivErodSilt - self.DepSilt)
                        * self.OutFracWat,
                    ),
                    0.0,
                ),
            )
            self.OutSandLoad = ifthen(
                self.subcatch > 0,
                cover(
                    max(
                        0.0,
                        (self.InSandLoad + self.RivErodSand - self.DepSand)
                        * self.OutFracWat,
                    ),
                    0.0,
                ),
            )
            self.OutSaggLoad = ifthen(
                self.subcatch > 0,
                cover(
                    max(
                        0.0,
                        (self.InSaggLoad + self.RivErodSagg - self.DepSagg)
                        * self.OutFracWat,
                    ),
                    0.0,
                ),
            )
            self.OutLaggLoad = ifthen(
                self.subcatch > 0,
                cover(
                    max(
                        0.0,
                        (self.InLaggLoad + self.RivErodLagg - self.DepLagg)
                        * self.OutFracWat,
                    ),
                    0.0,
                ),
            )
            self.OutGravLoad = ifthen(
                self.subcatch > 0,
                cover(
                    max(
                        0.0,
                        (self.InGravLoad + self.RivErodGrav - self.DepGrav)
                        * self.OutFracWat,
                    ),
                    0.0,
                ),
            )

            """ Mass balance, sediment concentration in each river cell """
            # Sediment load [ton/cell/timestep]
            # 0 in case all sediment are deposited in the cell
            self.SedLoad = max(
                0.0,
                self.InSedLoad + self.RivErodSed - self.DepSedLoad - self.OutSedLoad,
            )
            self.ClayLoad = max(
                0.0,
                self.InClayLoad + self.RivErodClay - self.DepClay - self.OutClayLoad,
            )
            self.SiltLoad = max(
                0.0,
                self.InSiltLoad + self.RivErodSilt - self.DepSilt - self.OutSiltLoad,
            )
            self.SandLoad = max(
                0.0,
                self.InSandLoad + self.RivErodSand - self.DepSand - self.OutSandLoad,
            )
            self.SaggLoad = max(
                0.0,
                self.InSaggLoad + self.RivErodSagg - self.DepSagg - self.OutSaggLoad,
            )
            self.LaggLoad = max(
                0.0,
                self.InLaggLoad + self.RivErodLagg - self.DepLagg - self.OutLaggLoad,
            )
            self.GravLoad = max(
                0.0,
                self.InGravLoad + self.RivErodGrav - self.DepGrav - self.OutGravLoad,
            )

            # Sediment concentration [mg/L]
            # The suspended load is composed of clay and silt.
            # The bed load is composed of sand, small and large aggregates from overland flow and gravel coming from river bed erosion.
            # Conversion from load [ton] to concentration for rivers [mg/L]
            # self.ToConc = 10**(6) / (self.RiverWidth * self.WaterLevel * self.DCL)
            self.ToConc = 10 ** (6) / (self.SurfaceRunoff * self.timestepsecs)
            self.SedConc = self.OutSedLoad * self.ToConc
            self.SSConc = (self.OutClayLoad + self.OutSiltLoad) * self.ToConc
            self.BedConc = (
                self.OutSandLoad
                + self.OutSaggLoad
                + self.OutLaggLoad
                + self.OutGravLoad
            ) * self.ToConc
            self.MaxSedConc = self.MaxSedLoad * self.ToConc

        """ Check / summary maps """
        # Checking processes contribution to sediment budget
        self.PrecipitationCatch = areatotal(self.Precipitation, self.HYPEcatch)
        # self.InLandSedCatch = areatotal(self.InLandSed, self.HYPEcatch)
        # self.MaxSedLoadCatch = areatotal(self.MaxSedLoad, self.HYPEcatch)
        # self.SedExCatch = areatotal(self.SedEx, self.HYPEcatch)
        # self.BankSedLoadCatch = areatotal(self.BankSedLoad, self.HYPEcatch)
        # self.BedSedLoadCatch = areatotal(self.BedSedLoad, self.HYPEcatch)
        # self.DegStoreSedCatch = areatotal(self.DegStoreSed, self.HYPEcatch)
        # self.DepSedLoadCatch = areatotal(self.DepSedLoad, self.HYPEcatch)
        # self.OutSedLoadCatch = areatotal(self.OutSedLoad, self.HYPEcatch)
        # self.SedLoadCatch = areatotal(self.SedLoad, self.HYPEcatch)

        # Area specific runoff
        # annual runoff / upstream (area)
        # L/km2/s
        # values between 5 and 25 for the Rhine
        self.SpecRun = (
            self.SurfaceRunoff * 1000 / accuflux(self.TopoLdd, self.cellareaKm)
        )  # [L/km2/s]


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
    timestepsecs = 86400
    wflow_cloneMap = "wflow_subcatch.map"

    # This allows us to use the model both on the command line and to call
    # the model using main function from another python script.

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
    # dynModelFw._runDynamic(0,0)
    dynModelFw._runDynamic(_firstTimeStep, _lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()


if __name__ == "__main__":
    main()
