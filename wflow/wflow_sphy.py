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


import sys
import os
import os.path
import getopt
import datetime as dt

from wflow.wf_DynamicFramework import *
from wflow.wflow_funcs import *
from wflow.wflow_adapt import *
import pcraster as pcr

# from wflow.wflow_lib import reporting

# import scipy
# import pcrut


wflow = "wflow_sphy"


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


class WflowModel(DynamicModel):

    """
  The user defined model class.

  """

    def __init__(self, cloneMap, Dir, RunDir, configfile):
        DynamicModel.__init__(self)
        self.caseName = os.path.abspath(Dir)
        self.clonemappath = os.path.join(os.path.abspath(Dir), "staticmaps", cloneMap)
        setclone(self.clonemappath)
        self.runId = RunDir
        self.Dir = os.path.abspath(Dir)
        self.configfile = configfile
        self.SaveDir = os.path.join(self.Dir, self.runId)

    def stateVariables(self):
        """
      returns a list of state variables that are essential to the model.
      This list is essential for the resume and suspend functions to work.

      This function is specific for each model and **must** be present.

     :var self.RootWater: Rootwater content [mm]
     :var self.SubWater: Water content in subsoil [mm]
     :var self.CapRise: Capillary rise [mm]
     :var self.RootDrain: Drainage from rootzone [mm]
     :var self.Scanopy: Water stored in canopy [mm]
     :var self.GwRecharge: Groundwater recharge [mm]
     :var self.BaseR: Baseflow [mm]
     :var self.Gw: Groundwater storage [mm]
     :var self.H_gw: Groundwater table height in [m] below surface
     :var self.SubDrain: Drainage from the subsoil [mm]
     :var self.SnowStore: Snow storage [mm]
     :var self.SnowWatStore: Water stored in snowpack [mm]
     :var self.GlacFrac: Glacier fraction map [-]
     :var self.Q/Rain/Snow/Glac/BaseRAold: Initial routed total / rainfall / 
             snow / glacier / baseflow runoff [m3/s]
     :var self.StorRes: Total storage in lakes and reservoirs [m3]
     :var self.Rain/Snow/Glac/BaseRAstor: Storage of rain / snow / glacier /
             baseflow waters in lakes and reservoirs [m3]

      """
        states = ["RootWater", "SubWater", "CapRise", "RootDrain"]

        if self.DynVegFLAG == 1:
            states.append("Scanopy")
        if self.GroundFLAG == 1:
            states.extend(["GwRecharge", "BaseR", "Gw", "H_gw"])
        else:
            states.append("SubDrain")
        if self.SnowFLAG == 1:
            states.extend(["SnowStore", "SnowWatStore"])
        if self.GlacFLAG == 1:
            states.append("GlacFrac")
        if self.RoutFLAG == 1 or self.LakeFLAG == 1 or self.ResFLAG == 1:
            states.append("QRAold")
            if self.RainRA_FLAG == 1:
                states.append("RainRAold")
            if self.SnowRA_FLAG == 1:
                states.append("SnowRAold")
            if self.GlacRA_FLAG == 1:
                states.append("GlacRAold")
            if self.BaseRA_FLAG == 1:
                states.append("BaseRAold")
        if self.LakeFLAG == 1 or self.ResFLAG == 1:
            states.append("StorRES")
            if self.RainRA_FLAG == 1:
                states.append("RainRAstor")
            if self.SnowRA_FLAG == 1:
                states.append("SnowRAstor")
            if self.GlacRA_FLAG == 1:
                states.append("GlacRAstor")
            if self.BaseRA_FLAG == 1:
                states.append("BaseRAstor")

        return states

    # The following are made to better connect to deltashell/openmi
    def supplyCurrentTime(self):  # - this may not be required
        """
      gets the current time in seconds after the start of the run

      Ouput:
          - time in seconds since the start of the model run
      """
        return self.currentTimeStep() * int(
            configget(self.config, "model", "timestepsecs", "86400")
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

        self.Prec_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Prec", "/inmaps/prec"
        )  # timeseries for rainfall
        self.Tair_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Tair", "/inmaps/tavg"
        )
        self.Tmax_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Tmax", "/inmaps/tmax"
        )
        self.Tmin_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Tmin", "/inmaps/tmin"
        )

        # Meteo and other forcing
        modelparameters.append(
            self.ParamType(
                name="Prec",
                stack=self.Prec_mapstack,
                type="timeseries",
                default=0.0,
                verbose=True,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="Tair",
                stack=self.Tair_mapstack,
                type="timeseries",
                default=0.0,
                verbose=True,
                lookupmaps=[],
            )
        )

        modelparameters.append(
            self.ParamType(
                name="Tmax",
                stack=self.Tmax_mapstack,
                type="timeseries",
                default=0.0,
                verbose=True,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="Tmin",
                stack=self.Tmin_mapstack,
                type="timeseries",
                default=0.0,
                verbose=True,
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

        if self.fewsrun:
            self.logger.info("Saving initial conditions for FEWS...")
            self.wf_suspend(os.path.join(self.Dir, "outstate"))

    def initial(self):

        global statistics
        global multpars
        global updateCols

        setglobaloption("unittrue")

        self.thestep = scalar(0)

        self.logger.info("running for " + str(self.nrTimeSteps()) + " timesteps")

        # Set and get defaults from ConfigFile here ###################################
        self.reinit = int(configget(self.config, "run", "reinit", "0"))
        self.fewsrun = int(configget(self.config, "run", "fewsrun", "0"))
        self.OverWriteInit = int(configget(self.config, "model", "OverWriteInit", "0"))

        # Print model info
        print(
            "The Spatial Processes in HYdrology (SPHY) model is "
            "developed and owned by FutureWater, Wageningen, The Netherlands"
        )
        print("Version 2.1")
        print(" ")

        # Read the modules to be used
        self.GlacFLAG = int(configget(self.config, "model", "GlacFLAG", "0"))
        self.SnowFLAG = int(configget(self.config, "model", "SnowFLAG", "0"))
        self.RoutFLAG = int(configget(self.config, "model", "RoutFLAG", "0"))
        self.ResFLAG = int(configget(self.config, "model", "ResFLAG", "0"))
        self.LakeFLAG = int(configget(self.config, "model", "LakeFLAG", "0"))
        self.DynVegFLAG = int(configget(self.config, "model", "DynVegFLAG", "0"))
        self.GroundFLAG = int(configget(self.config, "model", "GroundFLAG", "0"))
        # Optional modules
        if self.GroundFLAG == 0:
            self.SeepStatFLAG = int(configget(self.config, "model", "SeepStatic", "0"))
        if self.DynVegFLAG == 0:
            self.KcStatFLAG = int(configget(self.config, "model", "KCstatic", "0"))
        self.ETREF_FLAG = int(configget(self.config, "model", "ETREF_FLAG", "0"))

        if self.RoutFLAG == 1 or self.LakeFLAG == 1 or self.ResFLAG == 1:
            self.SnowRA_FLAG = int(configget(self.config, "model", "SnowRA_FLAG", "0"))
            self.RainRA_FLAG = int(configget(self.config, "model", "RainRA_FLAG", "0"))
            self.GlacRA_FLAG = int(configget(self.config, "model", "GlacRA_FLAG", "0"))
            self.BaseRA_FLAG = int(configget(self.config, "model", "BaseRA_FLAG", "0"))

        # import the required modules
        import datetime, calendar
        from wflow.sphy import reporting as reporting
        from wflow.sphy import timecalc as timecalc
        from wflow.sphy import ET as ET
        from wflow.sphy import rootzone as rootzone
        from wflow.sphy import subzone as subzone

        # from wflow.wflow_lib import *
        from math import pi

        # -standard python modules
        self.datetime = datetime
        self.calendar = calendar
        self.pi = pi
        # -FW defined modules
        self.reporting = reporting
        self.timecalc = timecalc
        self.ET = ET
        self.rootzone = rootzone
        self.subzone = subzone
        del datetime, calendar, pi, reporting, timecalc, ET, rootzone, subzone

        # -import additional modules if required
        if self.GlacFLAG == 1:
            self.SnowFLAG = 1
            self.GroundFLAG = 1
            import wflow.sphy.glacier as glacier  # glacier melting processes

            self.glacier = glacier
            del glacier
        if self.SnowFLAG == 1:
            import wflow.sphy.snow as snow  # snow melt processes

            self.snow = snow
            del snow
        if self.RoutFLAG == 1:
            import wflow.sphy.routing as routing  # simple routing scheme

            self.routing = routing
            del routing
        if self.LakeFLAG == 1:
            import wflow.sphy.lakes as lakes  # import lake module

            self.lakes = lakes
            del lakes
        if self.ResFLAG == 1:
            import wflow.sphy.reservoirs as reservoirs  # import reservoir module

            self.reservoirs = reservoirs
            del reservoirs
        if self.LakeFLAG == 1 or self.ResFLAG == 1:
            import wflow.sphy.advanced_routing as advanced_routing  # overwrite the simple routing scheme

            self.routing = advanced_routing
            del advanced_routing
            self.RoutFLAG = 0
        if self.DynVegFLAG == 1:
            import wflow.sphy.dynamic_veg as dynamic_veg  # dynamic crop growth using ndvi or kc time-series

            self.dynamic_veg = dynamic_veg
            del dynamic_veg
        if self.GroundFLAG == 1:
            import wflow.sphy.groundwater as groundwater  # groundwater storage as third storage layer. This is used instead of a fixed bottomflux

            self.groundwater = groundwater
            del groundwater

        # -set the global options
        setglobaloption("radians")
        # -set the 2000 julian date number
        self.julian_date_2000 = 2451545
        # -set the option to calculate the fluxes in mm for the upstream area
        self.mm_rep_FLAG = int(configget(self.config, "model", "mm_rep_FLAG", "1"))
        # -convert flow from m3/s to mm
        self.ToMM = 1000 * 3600 * 24 / cellarea()

        # static maps to use (normally default)
        wflow_dem = configget(self.config, "model", "dem", "staticmaps/dem.map")
        wflow_slope = configget(self.config, "model", "slope", "staticmaps/slope.map")
        wflow_locs = configget(
            self.config, "model", "locations", "staticmaps/outlets.map"
        )
        wflow_landuse = configget(
            self.config, "model", "landuse", "staticmaps/landuse.map"
        )
        wflow_ldd = configget(self.config, "model", "flowdir", "staticmaps/ldd.map")

        wflow_rootF = configget(
            self.config, "model", "RootFieldMap", "staticmaps/root_field.map"
        )
        wflow_rootS = configget(
            self.config, "model", "RootSatMap", "staticmaps/root_sat.map"
        )
        wflow_rootD = configget(
            self.config, "model", "RootDryMap", "staticmaps/root_dry.map"
        )
        wflow_rootW = configget(
            self.config, "model", "RootWiltMap", "staticmaps/root_wilt.map"
        )
        wflow_rootK = configget(
            self.config, "model", "RootKsatMap", "staticmaps/root_ksat.map"
        )
        wflow_subF = configget(
            self.config, "model", "SubFieldMap", "staticmaps/sub_field.map"
        )
        wflow_subS = configget(
            self.config, "model", "SubSatMap", "staticmaps/sub_sat.map"
        )
        wflow_subK = configget(
            self.config, "model", "SubKsat", "staticmaps/sub_ksat.map"
        )

        # 2: Input base maps ########################################################
        self.DEM = self.wf_readmap(os.path.join(self.Dir, wflow_dem), 0.0, fail=True)
        self.Slope = self.wf_readmap(
            os.path.join(self.Dir, wflow_slope), 0.0, fail=True
        )
        self.Locations = self.wf_readmap(
            os.path.join(self.Dir, wflow_locs), 0.0, fail=True
        )
        self.LandUse = self.wf_readmap(
            os.path.join(self.Dir, wflow_landuse), 1.0, fail=True
        )
        self.FlowDir = self.wf_readmap(
            os.path.join(self.Dir, wflow_ldd), 0.0, fail=False
        )

        self.RootFieldMap = self.wf_readmap(
            os.path.join(self.Dir, wflow_rootF), 0.35, fail=False
        )
        self.RootSatMap = self.wf_readmap(
            os.path.join(self.Dir, wflow_rootS), 0.45, fail=False
        )
        self.RootDryMap = self.wf_readmap(
            os.path.join(self.Dir, wflow_rootD), 0.1, fail=False
        )
        self.RootWiltMap = self.wf_readmap(
            os.path.join(self.Dir, wflow_rootW), 0.2, fail=False
        )
        self.RootKsat = self.wf_readmap(
            os.path.join(self.Dir, wflow_rootK), 20.0, fail=False
        )
        self.SubSatMap = self.wf_readmap(
            os.path.join(self.Dir, wflow_subS), 0.4, fail=False
        )
        self.SubFieldMap = self.wf_readmap(
            os.path.join(self.Dir, wflow_subF), 0.35, fail=False
        )
        self.SubKsat = self.wf_readmap(
            os.path.join(self.Dir, wflow_subK), 10, fail=False
        )

        # Set static initial values here #########################################
        self.ZeroMap = 0.0 * scalar(self.DEM)  # map with only zero's
        self.Latitude = ycoordinate(boolean(self.ZeroMap))
        self.Longitude = xcoordinate(boolean(self.ZeroMap))

        # Read parameters NEW Method
        self.logger.info("Linking parameters to landuse, catchment and soil...")
        self.wf_updateparameters()

        self.wf_multparameters()

        # Set static initial variables
        self.RootDrainVel = self.RootKsat * self.Slope

        if self.GroundFLAG == 0:
            self.SubDrainVel = self.SubKsat * self.Slope

        # -calculate soil properties
        self.RootField = self.RootFieldMap * self.RootDepthFlat
        self.RootSat = self.RootSatMap * self.RootDepthFlat
        self.RootDry = self.RootDryMap * self.RootDepthFlat
        self.RootWilt = self.RootWiltMap * self.RootDepthFlat
        self.SubSat = self.SubSatMap * self.SubDepthFlat
        self.SubField = self.SubFieldMap * self.SubDepthFlat
        self.RootTT = (self.RootSat - self.RootField) / self.RootKsat
        self.SubTT = (self.SubSat - self.SubField) / self.SubKsat

        # soil max and soil min for scaling of gwl if groundwater module is not used
        if self.GroundFLAG == 0:
            self.SoilMax = self.RootSat + self.SubSat
            self.SoilMin = self.RootDry + self.SubField

        if self.ETREF_FLAG == 0:
            from wflow.sphy import hargreaves

            self.Hargreaves = hargreaves
            del hargreaves

        setglobaloption("matrixtable")
        # -read lake maps and parameters if lake module is used
        if self.LakeFLAG == 1:
            # nominal map with lake IDs
            self.LakeID = cover(self.LakeID, 0)
            # lookup table with function for each lake (exp, 1-order poly, 2-order poly, 3-order poly)
            LakeFunc_Tab = os.path.join(self.Dir, "lake_function.tbl")
            # lookup table with Qh-coeficients for each lake
            LakeQH_Tab = os.path.join(self.Dir, "lake_QH.tbl")
            # lookup table with Sh-coeficients for each lake
            LakeSH_Tab = os.path.join(self.Dir, "lake_SH.tbl")
            # lookup table with hS-coeficients for each lake
            LakeHS_Tab = os.path.join(self.Dir, "lake_HS.tbl")
            # create lake coefficient maps
            self.LakeQH_Func = lookupnominal(LakeFunc_Tab, 1, self.LakeID)
            self.LakeSH_Func = lookupnominal(LakeFunc_Tab, 2, self.LakeID)
            self.LakeHS_Func = lookupnominal(LakeFunc_Tab, 3, self.LakeID)
            # Read QH coefficients
            self.LakeQH_exp_a = lookupscalar(LakeQH_Tab, 1, self.LakeID)
            self.LakeQH_exp_b = lookupscalar(LakeQH_Tab, 2, self.LakeID)
            self.LakeQH_pol_b = lookupscalar(LakeQH_Tab, 3, self.LakeID)
            self.LakeQH_pol_a1 = lookupscalar(LakeQH_Tab, 4, self.LakeID)
            self.LakeQH_pol_a2 = lookupscalar(LakeQH_Tab, 5, self.LakeID)
            self.LakeQH_pol_a3 = lookupscalar(LakeQH_Tab, 6, self.LakeID)
            # Read SH coefficients
            self.LakeSH_exp_a = lookupscalar(LakeSH_Tab, 1, self.LakeID)
            self.LakeSH_exp_b = lookupscalar(LakeSH_Tab, 2, self.LakeID)
            self.LakeSH_pol_b = lookupscalar(LakeSH_Tab, 3, self.LakeID)
            self.LakeSH_pol_a1 = lookupscalar(LakeSH_Tab, 4, self.LakeID)
            self.LakeSH_pol_a2 = lookupscalar(LakeSH_Tab, 5, self.LakeID)
            self.LakeSH_pol_a3 = lookupscalar(LakeSH_Tab, 6, self.LakeID)
            # Read HS coefficients
            self.LakeHS_exp_a = lookupscalar(LakeHS_Tab, 1, self.LakeID)
            self.LakeHS_exp_b = lookupscalar(LakeHS_Tab, 2, self.LakeID)
            self.LakeHS_pol_b = lookupscalar(LakeHS_Tab, 3, self.LakeID)
            self.LakeHS_pol_a1 = lookupscalar(LakeHS_Tab, 4, self.LakeID)
            self.LakeHS_pol_a2 = lookupscalar(LakeHS_Tab, 5, self.LakeID)
            self.LakeHS_pol_a3 = lookupscalar(LakeHS_Tab, 6, self.LakeID)
            # -Qfrac for lake or servoir cells should be zero, else 1
            self.QFRAC = ifthenelse(self.LakeID != 0, scalar(0), 1)

        # -read reservoir maps and parameters if reservoir module is used
        if self.ResFLAG == 1:
            # nominal map with reservoir IDs
            self.ResID = cover(self.ResID, 0)
            # lookup table with operational scheme to use (simple or advanced)
            ResFunc_Tab = os.path.join(self.Dir, "res_id")
            # Reservoir function
            self.ResFunc = cover(lookupscalar(ResFunc_Tab, 1, self.ResID), 0)
            try:
                # lookup table with coefficients for simple reservoirs
                ResSimple_Tab = os.path.join(self.Dir, "reservoir_simple")
                # Read coefficients for simple reservoirs
                self.ResKr = lookupscalar(ResSimple_Tab, 1, self.ResID)
                self.ResSmax = (
                    lookupscalar(ResSimple_Tab, 2, self.ResID) * 10 ** 6
                )  # convert to m3
                self.ResSimple = True
            except:
                self.ResSimple = False
            try:
                # lookup table with coefficients for advanced reservoirs
                ResAdvanced_Tab = os.path.join(self.Dir, "reservoir_advanced")
                # Read coefficients for advanced reservoirs
                self.ResEVOL = (
                    lookupscalar(ResAdvanced_Tab, 1, self.ResID) * 10 ** 6
                )  # convert to m3
                self.ResPVOL = (
                    lookupscalar(ResAdvanced_Tab, 2, self.ResID) * 10 ** 6
                )  # convert to m3
                self.ResMaxFl = (
                    lookupscalar(ResAdvanced_Tab, 3, self.ResID) * 10 ** 6
                )  # convert to m3/d
                self.ResDemFl = (
                    lookupscalar(ResAdvanced_Tab, 4, self.ResID) * 10 ** 6
                )  # convert to m3/d
                self.ResFlStart = lookupscalar(ResAdvanced_Tab, 5, self.ResID)
                self.ResFlEnd = lookupscalar(ResAdvanced_Tab, 6, self.ResID)
                self.ResAdvanced = True
            except:
                self.ResAdvanced = False
            # -Qfrac for reservoir cells should be zero, else 1
            if self.LakeFLAG == 1:
                self.QFRAC = ifthenelse(self.ResID != 0, scalar(0), self.QFRAC)
            else:
                self.QFRAC = ifthenelse(self.ResID != 0, scalar(0), 1)

    def default_summarymaps(self):  ##-maybe not needed. check later
        """
      Returns a list of default summary-maps at the end of a run.
      This is model specific. You can also add them to the [summary]section of the ini file but stuff
      you think is crucial to the model should be listed here

       Example:

      """
        # lst = ['self.Cfmax','self.csize','self.upsize','self.TTI','self.TT','self.WHC',
        #       'self.Slope','self.N','self.xl','self.yl','self.reallength','self.DCL','self.Bw',]
        if self.GlacFLAG == 1:
            lst = ["self.GlacFrac"]
        else:
            lst = []

        return lst

    def resume(self):
        """ read initial state maps (they are output of a previous call to suspend()) """

        if self.reinit == 1:  # -to be defined for sphy model state variables!!!
            self.logger.info("Setting initial conditions to default (zero!)")
            # -initial soil properties
            # -initial rootwater content
            self.RootWater = self.RootField
            self.SubWater = self.SubField
            # -initial capillary rise
            self.CapRise = self.ZeroMap + 3.0
            # -initial drainage from rootzone
            self.RootDrain = self.ZeroMap + 3.0
            if self.DynVegFLAG == 1:
                # -initial canopy storage
                self.Scanopy = self.ZeroMap
                # -initial ndvi if first map is not provided
                # self.ndviOld = scalar((self.NDVImax + self.NDVImin) / 2)
            # -initial groundwater properties
            if self.GroundFLAG == 1:
                self.GwRecharge = self.ZeroMap + 2.0
                self.BaseR = self.ZeroMap + 1.0
                self.Gw = self.ZeroMap + 1500
                self.H_gw = self.ZeroMap + 3.0
            else:
                self.SubDrain = self.ZeroMap + 3.0
            # -initial snow properties
            if self.SnowFLAG == 1:
                self.SnowStore = self.ZeroMap
                self.SnowWatStore = self.ZeroMap
            # -initial glacier properties
            if self.GlacFLAG == 1:
                self.GlacFrac = self.ZeroMap
            # -initial routed total runoff
            if self.RoutFLAG == 1 or self.LakeFLAG == 1 or self.ResFLAG == 1:
                self.QRAold = self.ZeroMap
                self.RainRAold = self.ZeroMap
                self.SnowRAold = self.ZeroMap
                self.GlacRAold = self.ZeroMap
                self.BaseRAold = self.ZeroMap
            # -initial storage in lakes and reservoirs
            if self.LakeFLAG == 1 or self.ResFLAG == 1:
                # -Read initial storages from table/reservoir file
                self.StorRES = self.ZeroMap
                self.RainRAstor = self.ZeroMap
                self.SnowRAstor = self.ZeroMap
                self.GlacRAstor = self.ZeroMap
                self.BaseRAstor = self.ZeroMap

                # -Read initial storages from table/reservoir file
                if self.LakeFLAG == 1:
                    LakeStor_Tab = os.path.join(self.Dir, "lake_id.tbl")
                    if os.path.exists(LakeStor_Tab):
                        self.StorRES = (
                            self.StorRES
                            + cover(lookupscalar(LakeStor_Tab, 1, self.LakeID), 0)
                            * 10 ** 6
                        )  # convert to m3
                        self.RainRAstor = (
                            self.RainRAstor
                            + cover(lookupscalar(LakeStor_Tab, 2, self.LakeID), 0)
                            * 10 ** 6
                        )
                        self.SnowRAstor = (
                            self.SnowRAstor
                            + cover(lookupscalar(LakeStor_Tab, 3, self.LakeID), 0)
                            * 10 ** 6
                        )
                        self.GlacRAstor = (
                            self.GlacRAstor
                            + cover(lookupscalar(LakeStor_Tab, 4, self.LakeID), 0)
                            * 10 ** 6
                        )
                        self.BaseRAstor = (
                            self.BaseRAstor
                            + cover(lookupscalar(LakeStor_Tab, 5, self.LakeID), 0)
                            * 10 ** 6
                        )
                    else:
                        self.logger.debug(
                            "Initial default state lake_id.tbl not found returning 0.0"
                        )
                if self.ResFLAG == 1:
                    ResStor_Tab = os.path.join(self.Dir, "reservoir_id.tbl")
                    if os.path.exists(ResStor_Tab):
                        self.StorRES = (
                            self.StorRES
                            + cover(lookupscalar(ResStor_Tab, 2, self.ResID), 0)
                            * 10 ** 6
                        )
                        self.RainRAstor = (
                            self.RainRAstor
                            + cover(lookupscalar(ResStor_Tab, 3, self.ResID), 0)
                            * 10 ** 6
                        )
                        self.SnowRAstor = (
                            self.SnowRAstor
                            + cover(lookupscalar(ResStor_Tab, 4, self.ResID), 0)
                            * 10 ** 6
                        )
                        self.GlacRAstor = (
                            self.GlacRAstor
                            + cover(lookupscalar(ResStor_Tab, 5, self.ResID), 0)
                            * 10 ** 6
                        )
                        self.BaseRAstor = (
                            self.BaseRAstor
                            + cover(lookupscalar(ResStor_Tab, 6, self.ResID), 0)
                            * 10 ** 6
                        )
                    else:
                        self.logger.debug(
                            "Initial default state res_id.tbl not found returning 0.0"
                        )

        else:
            self.wf_resume(os.path.join(self.Dir, "instate"))

        # -initial water storage in rootzone + subsoil
        self.SoilWater = self.RootWater + self.SubWater

        if self.SnowFLAG == 1:
            self.TotalSnowStore = self.SnowStore + self.SnowWatStore

    def dynamic(self):

        self.wf_updateparameters()  # read forcing an dynamic parameters

        # Snow and glacier fraction settings
        if self.GlacFLAG == 0:
            self.GlacFrac = scalar(0)
        if self.SnowFLAG == 0:
            self.SnowStore = scalar(0)
        SnowFrac = ifthenelse(self.SnowStore > 0, scalar(1 - self.GlacFrac), 0)
        RainFrac = ifthenelse(self.SnowStore == 0, scalar(1 - self.GlacFrac), 0)

        # -Read the precipitation time-series
        self.Precip = self.Prec
        self.PrecipF = self.Prec * (1 - self.GlacFrac)

        # -Temperature and reference evapotranspiration
        Temp = self.Tair
        if self.ETREF_FLAG == 0:
            TempMax = self.Tmax
            TempMin = self.Tmin
            self.ETref = self.Hargreaves.Hargreaves(
                pcr, self.Hargreaves.extrarad(self, pcr), Temp, TempMax, TempMin
            )

        # -Interception and effective precipitation
        # -Update canopy storage
        if self.DynVegFLAG == 1:
            # -fill missing ndvi values with ndvi base
            self.NDVI = ifthenelse(defined(self.NDVI) == 1, self.NDVI, self.NDVIbase)
            # -calculate the vegetation parameters
            vegoutput = self.dynamic_veg.Veg_function(
                pcr,
                self.NDVI,
                self.FPARmax,
                self.FPARmin,
                self.LAImax,
                self.NDVImin,
                self.NDVImax,
                self.KCmin,
                self.KCmax,
            )
            # -Kc
            self.Kc = vegoutput[0]
            # -Update canopy storage
            self.Scanopy = self.Scanopy + self.Precip
            # -interception and effective precipitation
            intercep = self.dynamic_veg.Inter_function(
                pcr, self.Scanopy, vegoutput[1], self.ETref
            )
            # -interception
            self.Int = intercep[0]
            # -interception corrected for fraction
            self.IntF = self.Int * (1 - self.GlacFrac)
            # -effective precipitation
            self.Precip = intercep[1]
            # -effective precipitation corrected for fraction
            self.PrecipEF = self.Precip * (1 - self.GlacFrac)

            # -canopy storage
            self.Scanopy = intercep[2]

        # Snow and rain
        if self.SnowFLAG == 1:
            # -Snow and rain differentiation
            self.Snow = ifthenelse(Temp >= self.Tcrit, 0, self.Precip)
            self.SnowF = self.Snow * (1 - self.GlacFrac)
            self.Rain = ifthenelse(Temp < self.Tcrit, 0, self.Precip)

            # -Snow melt
            PotSnowMelt = self.snow.PotSnowMelt(pcr, Temp, self.DDFS)
            self.ActSnowMelt = self.snow.ActSnowMelt(pcr, self.SnowStore, PotSnowMelt)
            self.ActSnowMeltF = self.ActSnowMelt * SnowFrac
            # -Update snow store
            self.SnowStore = self.snow.SnowStoreUpdate(
                pcr,
                self.SnowStore,
                self.Snow,
                self.ActSnowMelt,
                Temp,
                self.SnowWatStore,
            )
            # -Caclulate the maximum amount of water that can be stored in snowwatstore
            MaxSnowWatStore = self.snow.MaxSnowWatStorage(self.SnowSC, self.SnowStore)
            OldSnowWatStore = self.SnowWatStore
            # -Calculate the actual amount of water stored in snowwatstore
            self.SnowWatStore = self.snow.SnowWatStorage(
                pcr,
                Temp,
                MaxSnowWatStore,
                self.SnowWatStore,
                self.ActSnowMelt,
                self.Rain,
            )
            # -Changes in total water storage in snow (SnowStore and SnowWatStore)
            OldTotalSnowStore = self.TotalSnowStore
            self.TotalSnowStore = self.snow.TotSnowStorage(
                self.SnowStore, self.SnowWatStore, SnowFrac, RainFrac
            )
            # -Snow runoff
            self.SnowR = self.snow.SnowR(
                pcr,
                self.SnowWatStore,
                MaxSnowWatStore,
                self.ActSnowMelt,
                self.Rain,
                OldSnowWatStore,
                SnowFrac,
            )
        else:
            self.Rain = self.Precip
            self.SnowR = 0
            OldTotalSnowStore = 0
            self.TotalSnowStore = 0
        self.RainF = self.Rain * (1 - self.GlacFrac)

        # -Glacier calculations
        if self.GlacFLAG == 1:
            # -Glacier melt from clean ice glaciers
            GlacCIMelt = self.glacier.GlacCDMelt(pcr, Temp, self.DDFG, self.GlacFracCI)
            # -Glacier melt from debris covered glaciers
            GlacDCMelt = self.glacier.GlacCDMelt(pcr, Temp, self.DDFDG, self.GlacFracDB)
            # -Total melt from glaciers
            GlacMelt = self.glacier.GMelt(GlacCIMelt, GlacDCMelt)
            self.GlacMelt = GlacMelt
            self.GlacMeltF = self.GlacMelt * self.GlacFrac
            # -Glacier runoff
            self.GlacR = self.glacier.GlacR(self.GlacF, self.GlacMelt, self.GlacFrac)
            # -Glacier percolation to groundwater
            self.GlacPerc = self.glacier.GPerc(self.GlacF, self.GlacMelt, self.GlacFrac)
        else:
            self.GlacR = 0
            self.GlacMelt = 0
            self.GlacPerc = 0

        # -Potential evapotranspiration (THIS SHOULD STILL BE IMPROVED WITH DYNAMIC VEGETATION MODULE)
        self.ETpot = self.ET.ETpot(self.ETref, self.Kc)
        self.ETpotF = self.ETpot * RainFrac

        # -Rootzone calculations
        self.RootWater = (
            self.RootWater + ifthenelse(RainFrac > 0, self.Rain, 0) + self.CapRise
        )
        # -Rootzone runoff
        RootRunoff = self.rootzone.RootRunoff(
            pcr, RainFrac, self.RootWater, self.RootSat
        )
        self.RootWater = self.RootWater - RootRunoff
        # -Actual evapotranspiration
        etreddry = max(
            min((self.RootWater - self.RootDry) / (self.RootWilt - self.RootDry), 1), 0
        )
        self.ETact = self.ET.ETact(
            pcr, self.ETpot, self.RootWater, self.RootSat, etreddry, RainFrac
        )
        # -Actual evapotranspiration, corrected for rain fraction
        self.ActETact = self.ETact * RainFrac
        # -Update rootwater content
        self.RootWater = max(self.RootWater - self.ETact, 0)
        # -Rootwater drainage
        self.RootDrain = self.rootzone.RootDrainage(
            pcr,
            self.RootWater,
            self.RootDrain,
            self.RootField,
            self.RootSat,
            self.RootDrainVel,
            self.RootTT,
        )
        # -Update rootwater content
        self.RootWater = self.RootWater - self.RootDrain
        # -Rootwater percolation
        self.rootperc = self.rootzone.RootPercolation(
            pcr, self.RootWater, self.SubWater, self.RootField, self.RootTT, self.SubSat
        )
        self.rootpercF = self.rootperc * (1 - self.GlacFrac)
        # -Update rootwater content
        self.RootWater = self.RootWater - self.rootperc

        # -Sub soil calculations
        self.SubWater = self.SubWater + self.rootperc
        if self.GroundFLAG == 0:
            self.SubWater = min(max(self.SubWater - self.SeePage, 0), self.SubSat)
        # -Capillary rise
        self.CapRise = self.subzone.CapilRise(
            pcr,
            self.SubField,
            self.SubWater,
            self.CapRiseMax,
            self.RootWater,
            self.RootSat,
            self.RootField,
        )
        self.CapRiseF = self.CapRise * (1 - self.GlacFrac)

        # -Update sub soil water content
        self.SubWater = self.SubWater - self.CapRise
        if (
            self.GroundFLAG == 1
        ):  # sub percolation will be calculated instead of subdrainage
            subperc = self.subzone.SubPercolation(
                pcr, self.SubWater, self.SubField, self.SubTT, self.Gw, self.GwSat
            )
            self.ActSubPerc = subperc * (1 - self.GlacFrac)
            # -Update sub soil water content
            self.SubWater = self.SubWater - subperc
        else:  # sub drainage will be calculated instead of sub percolation
            self.SubDrain = self.subzone.SubDrainage(
                pcr,
                self.SubWater,
                self.SubField,
                self.SubSat,
                self.SubDrainVel,
                self.SubDrain,
                self.SubTT,
            )
            # -Update sub soil water content
            self.SubWater = self.SubWater - self.SubDrain

        # -Changes in soil water storage
        OldSoilWater = self.SoilWater
        self.SoilWater = (self.RootWater + self.SubWater) * (1 - self.GlacFrac)

        # -Rootzone runoff
        self.RootR = RootRunoff * RainFrac
        # -Rootzone drainage
        self.RootD = self.RootDrain * (1 - self.GlacFrac)
        # -Rain runoff
        self.RainR = self.RootR + self.RootD

        # -Groundwater calculations
        if self.GroundFLAG == 1:
            GwOld = self.Gw
            # -Groundwater recharge
            self.GwRecharge = self.groundwater.GroundWaterRecharge(
                pcr, self.deltaGw, self.GwRecharge, self.ActSubPerc, self.GlacPerc
            )
            # -Update groundwater storage
            self.Gw = self.Gw + self.GwRecharge
            # -Baseflow
            self.BaseR = self.groundwater.BaseFlow(
                pcr, self.Gw, self.BaseR, self.GwRecharge, self.BaseThresh, self.alphaGw
            )
            # -Update groundwater storage
            self.Gw = self.Gw - self.BaseR
            # -Calculate groundwater level
            self.H_gw = self.groundwater.HLevel(
                pcr, self.H_gw, self.alphaGw, self.GwRecharge, self.YieldGw
            )
            self.GWL = (
                (self.SubDepthFlat + self.RootDepthFlat + self.GwDepth) / 1000
                - self.H_gw
            ) * (-1)
        else:
            # -Use drainage from subsoil as baseflow
            self.BaseR = self.SubDrain
            # -Groundwater level as scaled between min and max measured gwl
            SoilAct = self.RootWater + self.SubWater
            SoilRel = (SoilAct - self.SoilMin) / (
                self.SoilMax - self.SoilMin
            )  # scale between 0 (dry) and 1 (wet)
            self.GWL = self.GWL_base - (SoilRel - 0.5) * self.GWL_base

        self.TotRF = self.BaseR + self.RainR + self.SnowR + self.GlacR

        # -Water balance
        if self.GroundFLAG == 1:
            self.waterbalance = (
                self.Precip * (1 - self.GlacFrac)
                + self.GlacMelt * self.GlacFrac
                - self.ActETact
                - self.GlacR
                - self.SnowR
                - self.RainR
                - self.BaseR
                - (self.SoilWater - OldSoilWater)
                - (self.TotalSnowStore - OldTotalSnowStore)
                - (self.Gw - GwOld)
            )
        elif self.GroundFLAG == 0:
            self.waterbalance = (
                self.Precip
                - self.ActETact
                - self.SeePage
                - self.SnowR
                - self.RainR
                - self.BaseR
                - (self.SoilWater - OldSoilWater)
                - (self.TotalSnowStore - OldTotalSnowStore)
            )

        # -Routing for lake and/or reservoir modules
        if self.LakeFLAG == 1 or self.ResFLAG == 1:
            # -Update storage in lakes/reservoirs (m3) with specific runoff
            self.StorRES = self.StorRES + ifthenelse(
                self.QFRAC == 0,
                0.001
                * cellarea()
                * (self.BaseR + self.RainR + self.GlacR + self.SnowR),
                0,
            )
            OldStorage = self.StorRES
            # -Calculate lake/reservoir outflow volumes
            if self.LakeFLAG == 1 and self.ResFLAG == 1:
                tempvar = self.lakes.UpdateLakeHStore(self, pcr, pcrm)
                LakeLevel = tempvar[0]
                self.StorRES = tempvar[1]
                LakeQ = self.lakes.QLake(self, pcr, LakeLevel)
                ResQ = self.reservoirs.QRes(self, pcr)
                self.Qout = ifthenelse(
                    self.ResID != 0, ResQ, ifthenelse(self.LakeID != 0, LakeQ, 0)
                )
            elif self.LakeFLAG == 1:
                tempvar = self.lakes.UpdateLakeHStore(self, pcr, pcrm)
                LakeLevel = tempvar[0]
                self.StorRES = tempvar[1]
                self.Qout = self.lakes.QLake(self, pcr, LakeLevel)
            else:
                self.Qout = self.reservoirs.QRes(self, pcr)

            # -Calculate volume available for routing (=outflow lakes/reservoir + cell specific runoff)
            RunoffVolume = upstream(self.FlowDir, self.Qout) + ifthenelse(
                self.QFRAC == 0,
                0,
                0.001
                * cellarea()
                * (self.BaseR + self.RainR + self.GlacR + self.SnowR),
            )
            # -Routing of total flow
            tempvar = self.routing.ROUT(
                self, pcr, RunoffVolume, self.QRAold, self.Qout, self.StorRES
            )
            self.StorRES = tempvar[0]
            self.Q = tempvar[1]
            self.Qin = tempvar[2]
            self.QRAold = self.Q

            # -Routing of individual contributers
            # -Snow routing
            if self.SnowRA_FLAG == 1 and self.SnowFLAG == 1:
                self.SnowRAstor = self.SnowRAstor + ifthenelse(
                    self.QFRAC == 0, self.SnowR * 0.001 * cellarea(), 0
                )
                cQfrac = cover(self.SnowRAstor / OldStorage, 0)
                self.cQout = cQfrac * self.Qout
                cRunoffVolume = upstream(self.FlowDir, self.cQout) + ifthenelse(
                    self.QFRAC == 0, 0, 0.001 * cellarea() * self.SnowR
                )
                tempvar = self.routing.ROUT(
                    self,
                    pcr,
                    cRunoffVolume,
                    self.SnowRAold,
                    self.cQout,
                    self.SnowRAstor,
                )
                self.SnowRAstor = tempvar[0]
                self.SnowRA = tempvar[1]
                self.cQin = tempvar[2]
                self.SnowRAold = self.SnowRA

            # -Rain routing
            if self.RainRA_FLAG == 1:
                self.RainRAstor = self.RainRAstor + ifthenelse(
                    self.QFRAC == 0, self.RainR * 0.001 * cellarea(), 0
                )
                cQfrac = cover(self.RainRAstor / OldStorage, 0)
                self.cQout = cQfrac * self.Qout
                cRunoffVolume = upstream(self.FlowDir, cQout) + ifthenelse(
                    self.QFRAC == 0, 0, 0.001 * cellarea() * self.RainR
                )
                tempvar = self.routing.ROUT(
                    self,
                    pcr,
                    cRunoffVolume,
                    self.RainRAold,
                    self.cQout,
                    self.RainRAstor,
                )
                self.RainRAstor = tempvar[0]
                self.RainRA = tempvar[1]
                self.cQin = tempvar[2]
                self.RainRAold = self.RainRA

            # -Glacier routing
            if self.GlacRA_FLAG == 1 and self.GlacFLAG == 1:
                self.GlacRAstor = self.GlacRAstor + ifthenelse(
                    self.QFRAC == 0, self.GlacR * 0.001 * cellarea(), 0
                )
                cQfrac = cover(self.GlacRAstor / OldStorage, 0)
                self.cQout = cQfrac * self.Qout
                cRunoffVolume = upstream(self.FlowDir, self.cQout) + ifthenelse(
                    self.QFRAC == 0, 0, 0.001 * cellarea() * self.GlacR
                )
                tempvar = self.routing.ROUT(
                    self,
                    pcr,
                    cRunoffVolume,
                    self.GlacRAold,
                    self.cQout,
                    self.GlacRAstor,
                )
                self.GlacRAstor = tempvar[0]
                self.GlacRA = tempvar[1]
                self.cQin = tempvar[2]
                self.GlacRAold = self.GlacRA

            # -Baseflow routing
            if self.BaseRA_FLAG == 1:
                self.BaseRAstor = self.BaseRAstor + ifthenelse(
                    self.QFRAC == 0, self.BaseR * 0.001 * cellarea(), 0
                )
                cQfrac = cover(self.BaseRAstor / OldStorage, 0)
                self.cQout = cQfrac * self.Qout
                cRunoffVolume = upstream(self.FlowDir, self.cQout) + ifthenelse(
                    self.QFRAC == 0, 0, 0.001 * cellarea() * self.BaseR
                )
                tempvar = self.routing.ROUT(
                    self,
                    pcr,
                    cRunoffVolume,
                    self.BaseRAold,
                    self.cQout,
                    self.BaseRAstor,
                )
                self.BaseRAstor = tempvar[0]
                self.BaseRA = tempvar[1]
                self.cQin = tempvar[2]
                self.BaseRAold = self.BaseRA

        # -Normal routing module
        elif self.RoutFLAG == 1:
            # -Rout total runoff
            self.Q = self.routing.ROUT(
                pcr,
                self.BaseR + self.RainR + self.GlacR + self.SnowR,
                self.QRAold,
                self.FlowDir,
                self.kx,
            )
            self.QRAold = self.Q
            # -Snow routing
            if self.SnowRA_FLAG == 1 and self.SnowFLAG == 1:
                self.SnowRA = self.routing.ROUT(
                    pcr, self.SnowR, self.SnowRAold, self.FlowDir, self.kx
                )
                self.SnowRAold = self.SnowRA

            # -Rain routing
            if self.RainRA_FLAG == 1:
                self.RainRA = self.routing.ROUT(
                    pcr, self.RainR, self.RainRAold, self.FlowDir, self.kx
                )
                self.RainRAold = self.RainRA
            # -Glacier routing
            if self.GlacRA_FLAG == 1 and self.GlacFLAG == 1:
                self.GlacRA = self.routing.ROUT(
                    pcr, self.GlacR, self.GlacRAold, self.FlowDir, self.kx
                )
                self.GlacRAold = self.GlacRA
            # -Baseflow routing
            if self.BaseRA_FLAG == 1:
                self.BaseRA = self.routing.ROUT(
                    pcr, self.BaseR, self.BaseRAold, self.FlowDir, self.kx
                )
                self.BaseRAold = self.BaseRA


# The main function is used to run the program from the command line


def main(argv=None):
    """
    Perform command line execution of the model.
    """
    global multpars
    global updateCols
    caseName = "wflow_ganga_sphy"
    runId = "run_default"
    configfile = "wflow_sphy.ini"
    LogFileName = "wflow.log"
    _lastTimeStep = 0
    _firstTimeStep = 0
    fewsrun = False
    runinfoFile = "runinfo.xml"
    timestepsecs = 86400
    wflow_cloneMap = "wflow_subcatch.map"
    NoOverWrite = 1
    loglevel = logging.DEBUG

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return

    ## Main model starts here
    ########################################################################
    try:
        opts, args = getopt.getopt(argv, "c:QXS:F:hC:Ii:T:R:u:s:P:p:Xx:U:fl:L:",['version'])
    except getopt.error as msg:
        pcrut.usage(msg)

    for o, a in opts:
        if o == "-F":
            runinfoFile = a
            fewsrun = True
        if o == "--version":
            import wflow
            print("wflow version: ", wflow.__version__)
            sys.exit(0)
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
            NoOverWrite = 0

    if fewsrun:
        ts = getTimeStepsfromRuninfo(runinfoFile, timestepsecs)
        starttime = getStartTimefromRuninfo(runinfoFile)
        if ts:
            _lastTimeStep = ts  # * 86400/timestepsecs
            _firstTimeStep = 1
        else:
            print("Failed to get timesteps from runinfo file: " + runinfoFile)
            sys.exit(2)
    else:
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
        NoOverWrite=NoOverWrite,
        logfname=LogFileName,
        level=loglevel,
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
