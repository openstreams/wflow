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


import numpy
import sys
import os
import os.path
import shutil, glob
import getopt
import datetime as dt

from wflow.wf_DynamicFramework import *
from wflow.wf_DynamicFramework import *
from wflow.wflow_adapt import *
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
        print msg
    print __doc__
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

    def updateRunOff(self):  # - this may not be required
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
            "RootWater",
            "SubWater",
            "CapRise",
            "RootDrain",
            "SubDrain",
            "GwRecharge",
            "Gw",
            "H_gw",
            "SnowWatStore",
        ]  ## -> complete list with required states

        # if hasattr(self,'ReserVoirSimpleLocs'):   ## -> add states that may be required if certain modules are used
        # states.append('ReservoirVolume')

        # if hasattr(self,'ReserVoirComplexLocs'):
        # states.append('ReservoirWaterLevel')

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
            self.config, "inputmapstacks", "Prec", "/inmaps/Prec"
        )  # timeseries for rainfall
        self.Tair_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Tair", "/inmaps/Tair"
        )  # timeseries for rainfall "/inmaps/TEMP"          # global radiation
        self.Tmax_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Tmax", "/inmaps/Tmax"
        )  # timeseries for rainfall "/inmaps/TEMP"          # global radiation
        self.Tmin_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Tmin", "/inmaps/Tmin"
        )  # timeseries for rainfall "/inmaps/TEMP"          # global radiation

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

        #: files to be used in case of timesries (scalar) input to the model

        # #: name of the tss file with precipitation data ("../intss/P.tss")
        # self.precipTss = "../intss/P.tss"
        # self.evapTss="../intss/PET.tss" #: name of the tss file with potential evap data ("../intss/PET.tss")
        # self.tempTss="../intss/T.tss" #: name of the tss file with temperature  data ("../intss/T.tss")
        # self.inflowTss="../intss/Inflow.tss" #: NOT TESTED name of the tss file with inflow data ("../intss/Inflow.tss")
        # self.SeepageTss="../intss/Seepage.tss" #: NOT TESTED name of the tss file with seepage data ("../intss/Seepage.tss")"

        self.logger.info("running for " + str(self.nrTimeSteps()) + " timesteps")

        # Set and get defaults from ConfigFile here ###################################
        # self.scalarInput = int(configget(self.config,"model","ScalarInput","0"))
        # self.Tslice = int(configget(self.config,"model","Tslice","1"))
        # self.interpolMethod = configget(self.config,"model","InterpolationMethod","inv")
        self.reinit = int(configget(self.config, "run", "reinit", "0"))
        self.fewsrun = int(configget(self.config, "run", "fewsrun", "0"))
        self.OverWriteInit = int(configget(self.config, "model", "OverWriteInit", "0"))
        # self.updating = int(configget(self.config,"model","updating","0"))
        # self.updateFile = configget(self.config,"model","updateFile","no_set")

        # self.sCatch = int(configget(self.config,"model","sCatch","0"))
        # self.intbl = configget(self.config,"model","intbl","intbl")
        # self.P_style = int(configget(self.config,"model","P_style","1"))
        # self.PET_style = int(configget(self.config,"model","PET_style","1"))
        # self.TEMP_style = int(configget(self.config,"model","TEMP_style","1"))

        # self.modelSnow = int(configget(self.config,"model","ModelSnow","1"))
        # sizeinmetres = int(configget(self.config,"layout","sizeinmetres","0"))
        # alf = float(configget(self.config,"model","Alpha","60"))
        # Qmax = float(configget(self.config,"model","AnnualDischarge","300"))
        # self.UpdMaxDist =float(configget(self.config,"model","UpdMaxDist","100"))
        # self.MaxUpdMult =float(configget(self.config,"model","MaxUpdMult","1.3"))
        # self.MinUpdMult =float(configget(self.config,"model","MinUpdMult","0.7"))
        # self.UpFrac =float(configget(self.config,"model","UpFrac","0.8"))
        # self.ExternalQbase=int(configget(self.config,'model','ExternalQbase','0'))
        # self.SetKquickFlow=int(configget(self.config,'model','SetKquickFlow','0'))
        # self.MassWasting = int(configget(self.config,"model","MassWasting","0"))
        # self.SubCatchFlowOnly = int(configget(self.config, 'model', 'SubCatchFlowOnly', '0'))

        # Print model info
        print "The Spatial Processes in HYdrology (SPHY) model is " "developed and owned by FutureWater, Wageningen, The Netherlands"
        print "Version 2.1"
        print " "

        # Read the modules to be used
        self.GlacFLAG = int(configget(self.config, "MODULES", "GlacFLAG", "0"))
        self.SnowFLAG = int(configget(self.config, "MODULES", "SnowFLAG", "0"))
        self.RoutFLAG = int(configget(self.config, "MODULES", "RoutFLAG", "0"))
        self.ResFLAG = int(configget(self.config, "MODULES", "ResFLAG", "0"))
        self.LakeFLAG = int(configget(self.config, "MODULES", "LakeFLAG", "0"))
        self.DynVegFLAG = int(configget(self.config, "MODULES", "DynVegFLAG", "0"))
        self.GroundFLAG = int(configget(self.config, "MODULES", "GroundFLAG", "0"))

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
            import lakes  # import lake module

            self.lakes = lakes
            del lakes
        if self.ResFLAG == 1:
            import reservoirs  # import reservoir module

            self.reservoirs = reservoirs
            del reservoirs
        if self.LakeFLAG == 1 or self.ResFLAG == 1:
            import advanced_routing  # overwrite the simple routing scheme

            self.routing = advanced_routing
            del advanced_routing
            self.RoutFLAG = 0
        if self.DynVegFLAG == 1:
            import dynamic_veg  # dynamic crop growth using ndvi or kc time-series

            self.dynamic_veg = dynamic_veg
            del dynamic_veg
        if self.GroundFLAG == 1:
            import wflow.sphy.groundwater as groundwater  # groundwater storage as third storage layer. This is used instead of a fixed bottomflux

            self.groundwater = groundwater
            del groundwater

        # -read the input and output directories from the configuration file
        # self.inpath = config.get('DIRS', 'inputdir')
        # self.inpathforcingT = config.get('DIRS','inputforcingdirT')
        # self.inpathforcingP = config.get('DIRS','inputforcingdirP')
        # self.outpath = config.get('DIRS', 'outputdir')

        # self.starttime = configget(self.config,"run","starttime","0")
        # ds = dt.datetime.strptime(self.starttime, '%Y-%m-%d %H:%M:%S %Z')
        # self.endtime = configget(self.config,"run","endtime","0")
        # de = dt.datetime.strptime(self.endtime, '%Y-%m-%d %H:%M:%S %Z')

        # #-set the timing criteria
        # sy = ds.year
        # sm = ds.month
        # sd = ds.day
        # ey = de.year
        # em = de.month
        # ed = de.day
        # self.startdate = self.datetime.datetime(sy,sm,sd)
        # self.enddate = self.datetime.datetime(ey,em,ed)

        # #-get start date of first forcing file in forcing directory
        # syF = config.getint('TIMING', 'startyear_F')
        # smF = config.getint('TIMING', 'startmonth_F')
        # sdF = config.getint('TIMING', 'startday_F')
        # self.startdateF = self.datetime.datetime(syF, smF, sdF)

        # -set the global options
        setglobaloption("radians")
        # -set the 2000 julian date number
        self.julian_date_2000 = 2451545
        # -set the option to calculate the fluxes in mm for the upstream area
        self.mm_rep_FLAG = int(configget(self.config, "REPORTING", "mm_rep_FLAG", "1"))

        # #-setting clone map
        # clonemap = self.inpath + config.get('GENERAL','mask')  ##->check
        # setclone(clonemap)
        # self.clone = readmap(clonemap)

        # -read general maps
        self.DEM = readmap(
            os.path.join(
                self.Dir,
                "staticmaps",
                configget(self.config, "GENERAL", "dem", "dem.map"),
            )
        )  # -> This has to be implemented for all readmap functions
        self.Slope = readmap(
            os.path.join(
                self.Dir,
                "staticmaps",
                configget(self.config, "GENERAL", "Slope", "slope.map"),
            )
        )
        self.Locations = readmap(
            os.path.join(
                self.Dir,
                "staticmaps",
                configget(self.config, "GENERAL", "locations", "outlets.map"),
            )
        )

        # -read soil maps
        # self.Soil = readmap(self.inpath + config.get('SOIL','Soil'))
        self.RootFieldMap = readmap(
            os.path.join(
                self.Dir,
                "staticmaps",
                configget(self.config, "SOIL", "RootFieldMap", "root_field.map"),
            )
        )
        self.RootSatMap = readmap(
            os.path.join(
                self.Dir,
                "staticmaps",
                configget(self.config, "SOIL", "RootSatMap", "root_sat.map"),
            )
        )
        self.RootDryMap = readmap(
            os.path.join(
                self.Dir,
                "staticmaps",
                configget(self.config, "SOIL", "RootDryMap", "root_dry.map"),
            )
        )
        self.RootWiltMap = readmap(
            os.path.join(
                self.Dir,
                "staticmaps",
                configget(self.config, "SOIL", "RootWiltMap", "root_wilt.map"),
            )
        )
        self.RootKsat = readmap(
            os.path.join(
                self.Dir,
                "staticmaps",
                configget(self.config, "SOIL", "RootKsat", "root_ksat.map"),
            )
        )
        self.SubSatMap = readmap(
            os.path.join(
                self.Dir,
                "staticmaps",
                configget(self.config, "SOIL", "SubSatMap", "sub_sat.map"),
            )
        )
        self.SubFieldMap = readmap(
            os.path.join(
                self.Dir,
                "staticmaps",
                configget(self.config, "SOIL", "SubFieldMap", "sub_field.map"),
            )
        )
        self.SubKsat = readmap(
            os.path.join(
                self.Dir,
                "staticmaps",
                configget(self.config, "SOIL", "SubKsat", "sub_ksat.map"),
            )
        )
        self.RootDrainVel = self.RootKsat * self.Slope

        # -Read and set the soil parameters
        pars = ["CapRiseMax", "RootDepthFlat", "SubDepthFlat"]
        for i in pars:
            try:
                # setattr(self, i, readmap(self.inpath + config.get('SOILPARS',i)))
                setattr(
                    self,
                    i,
                    readmap(
                        os.path.join(
                            self.Dir,
                            "staticmaps",
                            configget(self.config, "SOILPARS", i, i),
                        )
                    ),
                )
            except:
                # setattr(self, i, config.getfloat('SOILPARS',i))
                setattr(self, i, float(configget(self.config, "SOILPARS", i, i)))
        if (
            self.GroundFLAG == 0
        ):  # if groundwater module is not used, read seepage and gwl_base
            self.SeepStatFLAG = config.getint("SOILPARS", "SeepStatic")
            if self.SeepStatFLAG == 0:  # set the seepage map series
                self.Seepmaps = self.inpath + config.get("SOILPARS", "SeePage")
            else:  # -set a static map or value for seepage
                try:
                    self.SeePage = readmap(
                        self.inpath + config.get("SOILPARS", "SeePage")
                    )
                except:
                    self.SeePage = config.getfloat("SOILPARS", "SeePage")
            try:
                self.GWL_base = readmap(
                    self.inpath + config.get("SOILPARS", "GWL_base")
                )
            except:
                self.GWL_base = config.getfloat("SOILPARS", "GWL_base")

            self.SubDrainVel = self.SubKsat * self.Slope
        else:  # if groundwater module is used, then read the groundwater soil parameters
            pars = ["GwDepth", "GwSat", "deltaGw", "BaseThresh", "alphaGw", "YieldGw"]
            for i in pars:
                try:
                    setattr(
                        self, i, readmap(self.inpath + config.get("GROUNDW_PARS", i))
                    )
                except:
                    setattr(
                        self, i, float(configget(self.config, "GROUNDW_PARS", i, i))
                    )

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

        # -read the crop coefficient table if the dynamic vegetation module is not used
        if self.DynVegFLAG == 0:
            self.KcStatFLAG = int(
                configget(self.config, "LANDUSE", "KCstatic", "kc.tbl")
            )
            if self.KcStatFLAG == 1:
                # -read land use map and kc table
                self.LandUse = readmap(
                    os.path.join(
                        self.Dir,
                        "staticmaps",
                        configget(self.config, "LANDUSE", "LandUse", "landuse.map"),
                    )
                )
                self.kc_table = os.path.join(
                    self.Dir,
                    "staticmaps",
                    configget(self.config, "LANDUSE", "CropFac", "kc.tbl"),
                )
                self.Kc = lookupscalar(self.kc_table, self.LandUse)
            else:
                # -set the kc map series
                self.Kcmaps = self.inpath + config.get("LANDUSE", "KC")
        # -Use the dynamic vegetation module
        else:
            # -set the ndvi map series to be read
            self.ndvi = self.inpath + config.get("DYNVEG", "NDVI")
            # -read the vegetation parameters
            pars = [
                "NDVImax",
                "NDVImin",
                "NDVIbase",
                "KCmax",
                "KCmin",
                "LAImax",
                "FPARmax",
                "FPARmin",
            ]
            for i in pars:
                try:
                    setattr(self, i, readmap(self.inpath + config.get("DYNVEG", i)))
                except:
                    setattr(self, i, config.getfloat("DYNVEG", i))

        # -read and set glacier maps and parameters if glacier module is used
        if self.GlacFLAG == 1:
            # self.GlacFracCI = readmap(self.inpath + config.get('GLACIER','GlacFracCI'))
            # self.GlacFracDB = readmap(self.inpath + config.get('GLACIER','GlacFracDB'))
            self.GlacFracCI = readmap(
                os.path.join(
                    self.Dir,
                    "staticmaps",
                    configget(
                        self.config, "GLACIER", "GlacFracCI", "glacier_clean.map"
                    ),
                )
            )
            self.GlacFracDB = readmap(
                os.path.join(
                    self.Dir,
                    "staticmaps",
                    configget(
                        self.config, "GLACIER", "GlacFracDB", "glacier_debris.map"
                    ),
                )
            )
            pars = ["DDFG", "DDFDG", "GlacF"]
            for i in pars:
                try:
                    setattr(self, i, readmap(self.inpath + config.get("GLACIER", i)))
                except:
                    # setattr(self, i, config.getfloat('GLACIER',i))
                    setattr(self, i, float(configget(self.config, "GLACIER", i, i)))

        # -read and set snow maps and parameters if snow modules are used
        if self.SnowFLAG == 1:
            pars = ["Tcrit", "SnowSC", "DDFS"]
            for i in pars:
                try:
                    setattr(self, i, readmap(self.inpath + config.get("SNOW", i)))
                except:
                    #    			setattr(self, i, float(configget(self.config,'SNOW',i,i)))
                    setattr(self, i, float(configget(self.config, "SNOW", i, i)))

        # -read and set climate forcing and the calculation of etref

        self.Prec_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Prec", "/inmaps/Prec"
        )  # timeseries for rainfall
        self.Tair_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Tair", "/inmaps/Tair"
        )  # timeseries for rainfall "/inmaps/TEMP"          # global radiation
        self.Tmax_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Tmax", "/inmaps/Tmax"
        )  # timeseries for rainfall "/inmaps/TEMP"          # global radiation
        self.Tmin_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "Tmin", "/inmaps/Tmin"
        )  # timeseries for rainfall "/inmaps/TEMP"          # global radiation

        # self.Prec = self.inpathforcingP + config.get('CLIMATE','Prec')
        # self.Tair = self.inpathforcingT + config.get('CLIMATE','Tair')
        # self.ETREF_FLAG = config.getint('ETREF','ETREF_FLAG')   ##-> for now should be zero.
        self.ETREF_FLAG = int(
            configget(self.config, "ETREF", "ETREF_FLAG", 0)
        )  ##-> for now should be zero.
        # -determine the use of a given etref time-series or calculate etref using Hargreaves
        if self.ETREF_FLAG == 1:
            self.ETref = self.inpath + config.get("ETREF", "ETref")
        else:
            # self.Lat = readmap(self.inpath + config.get('ETREF','Lat'))
            self.Lat = readmap(
                os.path.join(
                    self.Dir,
                    "staticmaps",
                    configget(self.config, "ETREF", "Lat", "latitude.map"),
                )
            )
            self.Tmax_mapstack = self.Dir + configget(
                self.config, "inputmapstacks", "Tmax", "/inmaps/Tmax"
            )  # timeseries for rainfall "/inmaps/TEMP"          # global radiation
            self.Tmin_mapstack = self.Dir + configget(
                self.config, "inputmapstacks", "Tmin", "/inmaps/Tmin"
            )  # timeseries for rainfall "/inmaps/TEMP"          # global radiation
            # self.Gsc = config.getfloat('ETREF', 'Gsc')
            self.Gsc = float(configget(self.config, "ETREF", "Gsc", 0.0820))
            from wflow.sphy import hargreaves

            self.Hargreaves = hargreaves
            del hargreaves

        # -read and set routing maps and parameters
        if self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1:
            # self.FlowDir = readmap(self.inpath + config.get('ROUTING','flowdir'))
            self.FlowDir = readmap(
                os.path.join(
                    self.Dir,
                    "staticmaps",
                    configget(self.config, "ROUTING", "flowdir", "ldd.map"),
                )
            )
            try:
                self.kx = readmap(self.inpath + config.get("ROUTING", "kx"))
            except:
                # self.kx = config.getfloat('ROUTING','kx')
                self.kx = float(configget(self.config, "ROUTING", "kx", 1))

        setglobaloption("matrixtable")
        # -read lake maps and parameters if lake module is used
        if self.LakeFLAG == 1:
            # nominal map with lake IDs
            self.LakeID = cover(readmap(self.inpath + config.get("LAKE", "LakeId")), 0)
            # lookup table with function for each lake (exp, 1-order poly, 2-order poly, 3-order poly)
            LakeFunc_Tab = self.inpath + config.get("LAKE", "LakeFunc")
            # lookup table with Qh-coeficients for each lake
            LakeQH_Tab = self.inpath + config.get("LAKE", "LakeQH")
            # lookup table with Sh-coeficients for each lake
            LakeSH_Tab = self.inpath + config.get("LAKE", "LakeSH")
            # lookup table with hS-coeficients for each lake
            LakeHS_Tab = self.inpath + config.get("LAKE", "LakeHS")
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
            # -read water level maps and parameters if available
            try:
                self.UpdateLakeLevel = readmap(
                    self.inpath + config.get("LAKE", "updatelakelevel")
                )
                self.LLevel = self.inpath + config.get("LAKE", "LakeFile")
                print "measured lake levels will be used to update lake storage"
            except:
                pass

        # -read reservior maps and parameters if reservoir module is used
        if self.ResFLAG == 1:
            # nominal map with reservoir IDs
            self.ResID = cover(
                readmap(self.inpath + config.get("RESERVOIR", "ResId")), 0
            )
            # lookup table with operational scheme to use (simple or advanced)
            ResFunc_Tab = self.inpath + config.get("RESERVOIR", "ResFuncStor")
            # Reservoir function
            self.ResFunc = cover(lookupscalar(ResFunc_Tab, 1, self.ResID), 0)
            try:
                # lookup table with coefficients for simple reservoirs
                ResSimple_Tab = self.inpath + config.get("RESERVOIR", "ResSimple")
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
                ResAdvanced_Tab = self.inpath + config.get("RESERVOIR", "ResAdv")
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

        # -below is original initial part from sphy

        # -get the correct forcing file number, depending on the start date of your simulation
        # -and the start date of the first forcing file in your forcing directory.
        # self.counter = (self.startdate - self.startdateF).days
        # #-initial date
        # self.curdate = self.startdate
        # self.curdate = configget(self.config,"run","starttime","0")
        # print self.curdate
        # -initial soil properties
        # -initial rootwater content

        self.RootWater = self.RootField
        self.SubWater = self.SubField

        # if not config.get('SOIL_INIT','RootWater'):
        # self.RootWater = self.RootField
        # else:
        # try:
        # self.RootWater = config.getfloat('SOIL_INIT','RootWater')
        # except:
        # self.RootWater = readmap(self.inpath + config.get('SOIL_INIT','RootWater'))
        # #-initial water content in subsoil
        # if not config.get('SOIL_INIT','SubWater'):
        # self.SubWater = self.SubField
        # else:
        # try:
        # self.SubWater = config.getfloat('SOIL_INIT','SubWater')
        # except:
        # self.SubWater = readmap(self.inpath + config.get('SOIL_INIT','SubWater'))
        # -initial water storage in rootzone + subsoil
        self.SoilWater = self.RootWater + self.SubWater
        # -initial capillary rise
        self.CapRise = configget(self.config, "SOIL_INIT", "CapRise", 3)
        # try:
        # self.CapRise = config.getfloat('SOIL_INIT','CapRise')
        # except:
        # self.CapRise = readmap(self.inpath + config.get('SOIL_INIT','CapRise'))
        # -initial drainage from rootzone
        self.RootDrain = configget(self.config, "SOIL_INIT", "RootDrain", 3)
        # try:
        # self.RootDrain = config.getfloat('SOIL_INIT','RootDrain')
        # except:
        # self.RootDrain = readmap(self.inpath + config.get('SOIL_INIT','RootDrain'))

        if self.DynVegFLAG == 1:
            # -initial canopy storage
            self.Scanopy = 0
            # -initial ndvi if first map is not provided
            self.ndviOld = scalar((self.NDVImax + self.NDVImin) / 2)
        elif self.KcStatFLAG == 0:
            # -set initial kc value to one, if kc map is not available for first timestep
            self.KcOld = scalar(1)

        # -initial groundwater properties
        if self.GroundFLAG == 1:
            self.GwRecharge = float(
                configget(self.config, "GROUNDW_INIT", "GwRecharge", 0)
            )
            self.BaseR = float(configget(self.config, "GROUNDW_INIT", "BaseR", 1))
            self.Gw = float(configget(self.config, "GROUNDW_INIT", "Gw", 1))
            self.H_gw = float(configget(self.config, "GROUNDW_INIT", "H_gw", 1))
            # #-initial groundwater recharge
            # try:
            # self.GwRecharge = config.getfloat('GROUNDW_INIT','GwRecharge')
            # except:
            # self.GwRecharge = readmap(self.inpath + config.get('GROUNDW_INIT','GwRecharge'))
            # #-initial baseflow
            # try:
            # self.BaseR = config.getfloat('GROUNDW_INIT','BaseR')
            # except:
            # self.BaseR = readmap(self.inpath + config.get('GROUNDW_INIT','BaseR'))
            # #-initial groundwater storage
            # try:
            # self.Gw = config.getfloat('GROUNDW_INIT','Gw')
            # except:
            # self.Gw = readmap(self.inpath + config.get('GROUNDW_INIT','Gw'))
            # #-initial groundwater level
            # try:
            # self.H_gw = config.getfloat('GROUNDW_INIT','H_gw')
            # except:
            # self.H_gw = readmap(self.inpath + config.get('GROUNDW_INIT','H_gw'))
            # self.H_gw = max((self.RootDepthFlat + self.SubDepthFlat + self.GwDepth)/1000 - self.H_gw, 0)
        # else:
        # #-initial drainage from subsoil
        # try:
        # self.SubDrain = config.getfloat('SOIL_INIT','SubDrain')
        # except:
        # self.SubDrain = readmap(self.inpath + config.get('SOIL_INIT','SubDrain'))
        # #-initial seepage value if seepage map series is used
        # if self.SeepStatFLAG == 0:
        # self.SeepOld = scalar(0)

        # -initial snow properties
        if self.SnowFLAG == 1:
            try:
                # self.SnowStore = config.getfloat('SNOW_INIT','SnowIni')
                self.SnowStore = float(
                    configget(self.config, "SNOW_INIT", "SnowIni", 0)
                )
            except:
                self.SnowStore = readmap(
                    self.inpath + config.get("SNOW_INIT", "SnowIni")
                )
                # -initial water stored in snowpack
            try:
                self.SnowWatStore = float(
                    configget(self.config, "SNOW_INIT", "SnowWatStore", 0)
                )
                # self.SnowWatStore = config.getfloat('SNOW_INIT','SnowWatStore')
            except:
                self.SnowWatStore = readmap(
                    self.inpath + config.get("SNOW_INIT", "SnowWatStore")
                )
            self.TotalSnowStore = self.SnowStore + self.SnowWatStore

        # -initial glacier properties
        if self.GlacFLAG == 1:
            # try:
            # self.GlacFrac = config.getfloat('GLACIER_INIT','GlacFrac')
            # except:
            # self.GlacFrac = readmap(self.inpath + config.get('GLACIER_INIT','GlacFrac'))
            self.GlacFrac = readmap(
                os.path.join(
                    self.Dir,
                    "staticmaps",
                    configget(
                        self.config, "GLACIER_INIT", "GlacFrac", "glacierfraction.map"
                    ),
                )
            )
            print self.GlacFrac
        # -initial routed total runoff and of individual components
        if self.RoutFLAG == 1 or self.LakeFLAG == 1 or self.ResFLAG == 1:
            # -initial routed total runoff
            try:
                self.QRAold = config.getfloat("ROUT_INIT", "QRA_init")
            except:
                try:
                    self.QRAold = readmap(
                        self.inpath + config.get("ROUT_INIT", "QRA_init")
                    )
                except:
                    self.QRAold = 0
                    # -initial routed runoff	for the individual components
            pars = ["RainRA", "SnowRA", "GlacRA", "BaseRA"]
            self.RainRAold = 0
            self.SnowRAold = 0
            self.GlacRAold = 0
            self.BaseRAold = 0
            self.RainRA_FLAG = True
            self.SnowRA_FLAG = True
            self.GlacRA_FLAG = True
            self.BaseRA_FLAG = True
            # for i in pars:
            # try:
            # setattr(self, i + 'old', readmap(self.inpath + config.get('ROUT_INIT', i + '_init')))
            # setattr(self, i + '_FLAG', True)
            # except:
            # try:
            # #setattr(self, i + 'old', config.getfloat('ROUT_INIT', i + '_init'))
            # setattr(self, i + '_FLAG', True)
            # print RainRA_init
            # except:
            # setattr(self, i + '_FLAG', False)

        # -initial storage in lakes and reservoirs
        if self.LakeFLAG == 1 or self.ResFLAG == 1:
            # -Read initial storages from table/reservoir file
            if self.LakeFLAG == 1:
                LakeStor_Tab = self.inpath + config.get("LAKE", "LakeStor")
                self.StorRES = (
                    cover(lookupscalar(LakeStor_Tab, 1, self.LakeID), 0) * 10 ** 6
                )  # convert to m3
                # -Qfrac for lake cells should be zero, else 1
                self.QFRAC = ifthenelse(self.LakeID != 0, scalar(0), 1)
            if self.ResFLAG == 1:
                ResStor_Tab = self.inpath + config.get("RESERVOIR", "ResFuncStor")
                ResStor = (
                    cover(lookupscalar(ResStor_Tab, 2, self.ResID), 0) * 10 ** 6
                )  # convert to m3
                try:
                    self.StorRES = self.StorRES + ResStor
                    # -Qfrac for reservoir cells should be zero, else 1
                    self.QFRAC = ifthenelse(self.ResID != 0, scalar(0), self.QFRAC)
                except:
                    self.StorRES = ResStor
                    # -Qfrac for reservoir cells should be zero, else 1
                    self.QFRAC = ifthenelse(self.ResID != 0, scalar(0), 1)

                    # -initial storage in lakes/reservoirs of individual flow components
            pars = ["RainRA", "SnowRA", "GlacRA", "BaseRA"]
            for i in pars:
                column = pars.index(
                    i
                )  # identify column to be read from lake or reservoir table
                try:  # -try to sum the storages read from the lake and reservoir tables if both thse modules are used
                    setattr(
                        self,
                        i + "stor",
                        (
                            cover(
                                lookupscalar(LakeStor_Tab, column + 2, self.LakeID), 0
                            )
                            + cover(
                                lookupscalar(ResStor_Tab, column + 3, self.ResID), 0
                            )
                        )
                        * 10 ** 6,
                    )
                    if eval("self." + i + "_FLAG"):
                        setattr(self, i + "_FLAG", True)
                    else:
                        setattr(self, i + "_FLAG", False)
                except:
                    try:  # -try to read the storages from the lake table
                        setattr(
                            self,
                            i + "stor",
                            cover(
                                lookupscalar(LakeStor_Tab, column + 2, self.LakeID), 0
                            )
                            * 10 ** 6,
                        )
                        if eval("self." + i + "_FLAG"):
                            setattr(self, i + "_FLAG", True)
                        else:
                            setattr(self, i + "_FLAG", False)
                    except:  # -try to read the storages from the reservoir table
                        try:
                            setattr(
                                self,
                                i + "stor",
                                cover(
                                    lookupscalar(ResStor_Tab, column + 3, self.ResID), 0
                                )
                                * 10 ** 6,
                            )
                            if eval("self." + i + "_FLAG"):
                                setattr(self, i + "_FLAG", True)
                            else:
                                setattr(self, i + "_FLAG", False)
                        except:
                            setattr(self, i + "_FLAG", False)

        # -Initial values for reporting and setting of time-series
        # -set time-series reporting for mm flux from upstream area for prec and eta
        if self.mm_rep_FLAG == 1 and (
            self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1
        ):
            self.PrecSubBasinTSS = TimeoutputTimeseries(
                "PrecSubBasinTSS", self, self.Locations, noHeader=False
            )
            self.ETaSubBasinTSS = TimeoutputTimeseries(
                "ETaSubBasinTSS", self, self.Locations, noHeader=False
            )
        if self.GlacFLAG == 1:
            pars = [
                "wbal",
                "GWL",
                "TotPrec",
                "TotPrecF",
                "TotPrecEF",
                "TotIntF",
                "TotRain",
                "TotRainF",
                "TotETpot",
                "TotETpotF",
                "TotETact",
                "TotETactF",
                "TotSnow",
                "TotSnowF",
                "TotSnowMelt",
                "TotSnowMeltF",
                "TotGlacMelt",
                "TotGlacMeltF",
                "TotRootRF",
                "TotRootDF",
                "TotRootPF",
                "TotSubPF",
                "TotCapRF",
                "TotGlacPercF",
                "TotGwRechargeF",
                "TotRainRF",
                "TotBaseRF",
                "TotSnowRF",
                "TotGlacRF",
                "TotRF",
                "RainRAtot",
                "SnowRAtot",
                "GlacRAtot",
                "BaseRAtot",
                "QallRAtot",
            ]
            # -set time-series reporting for mm fluxes from upstream area if glacier and routing/reservoir modules are used
            if self.mm_rep_FLAG == 1 and (
                self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1
            ):
                self.GMeltSubBasinTSS = TimeoutputTimeseries(
                    "GMeltSubBasinTSS", self, self.Locations, noHeader=False
                )
                self.QSNOWSubBasinTSS = TimeoutputTimeseries(
                    "QSNOWSubBasinTSS", self, self.Locations, noHeader=False
                )
                self.QRAINSubBasinTSS = TimeoutputTimeseries(
                    "QRAINSubBasinTSS", self, self.Locations, noHeader=False
                )
                self.QGLACSubBasinTSS = TimeoutputTimeseries(
                    "QGLACSubBasinTSS", self, self.Locations, noHeader=False
                )
                self.QBASFSubBasinTSS = TimeoutputTimeseries(
                    "QBASFSubBasinTSS", self, self.Locations, noHeader=False
                )
                self.QTOTSubBasinTSS = TimeoutputTimeseries(
                    "QTOTSubBasinTSS", self, self.Locations, noHeader=False
                )
        elif self.SnowFLAG == 1:
            if self.GroundFLAG == 1:
                pars = [
                    "wbal",
                    "GWL",
                    "TotPrec",
                    "TotPrecF",
                    "TotPrecEF",
                    "TotIntF",
                    "TotRain",
                    "TotRainF",
                    "TotETpot",
                    "TotETpotF",
                    "TotETact",
                    "TotETactF",
                    "TotSnow",
                    "TotSnowF",
                    "TotSnowMelt",
                    "TotSnowMeltF",
                    "TotRootRF",
                    "TotRootDF",
                    "TotRootPF",
                    "TotSubPF",
                    "TotCapRF",
                    "TotGwRechargeF",
                    "TotRainRF",
                    "TotBaseRF",
                    "TotSnowRF",
                    "TotRF",
                    "RainRAtot",
                    "SnowRAtot",
                    "BaseRAtot",
                    "QallRAtot",
                ]
                # -set time-series reporting for mm fluxes from upstream area if snow, groundwater and routing/reservoir modules are used
                if self.mm_rep_FLAG == 1 and (
                    self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1
                ):
                    self.QSNOWSubBasinTSS = TimeoutputTimeseries(
                        "QSNOWSubBasinTSS", self, self.Locations, noHeader=False
                    )
                    self.QRAINSubBasinTSS = TimeoutputTimeseries(
                        "QRAINSubBasinTSS", self, self.Locations, noHeader=False
                    )
                    self.QBASFSubBasinTSS = TimeoutputTimeseries(
                        "QBASFSubBasinTSS", self, self.Locations, noHeader=False
                    )
                    self.QTOTSubBasinTSS = TimeoutputTimeseries(
                        "QTOTSubBasinTSS", self, self.Locations, noHeader=False
                    )
            else:
                pars = [
                    "wbal",
                    "GWL",
                    "TotPrec",
                    "TotPrecF",
                    "TotPrecEF",
                    "TotIntF",
                    "TotRain",
                    "TotRainF",
                    "TotETpot",
                    "TotETpotF",
                    "TotETact",
                    "TotETactF",
                    "TotSnow",
                    "TotSnowF",
                    "TotSnowMelt",
                    "TotSnowMeltF",
                    "TotRootRF",
                    "TotRootDF",
                    "TotRootPF",
                    "TotSubDF",
                    "TotCapRF",
                    "TotSeepF",
                    "TotRainRF",
                    "TotSnowRF",
                    "TotRF",
                    "RainRAtot",
                    "SnowRAtot",
                    "BaseRAtot",
                    "QallRAtot",
                ]
                # -set time-series reporting for mm fluxes from upstream area if snow and routing/reservoir modules are used
                if self.mm_rep_FLAG == 1 and (
                    self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1
                ):
                    self.SeepSubBasinTSS = TimeoutputTimeseries(
                        "SeepSubBasinTSS", self, self.Locations, noHeader=False
                    )
                    self.QSNOWSubBasinTSS = TimeoutputTimeseries(
                        "QSNOWSubBasinTSS", self, self.Locations, noHeader=False
                    )
                    self.QRAINSubBasinTSS = TimeoutputTimeseries(
                        "QRAINSubBasinTSS", self, self.Locations, noHeader=False
                    )
                    self.QBASFSubBasinTSS = TimeoutputTimeseries(
                        "QBASFSubBasinTSS", self, self.Locations, noHeader=False
                    )
                    self.QTOTSubBasinTSS = TimeoutputTimeseries(
                        "QTOTSubBasinTSS", self, self.Locations, noHeader=False
                    )
        else:
            if self.GroundFLAG == 1:
                pars = [
                    "wbal",
                    "GWL",
                    "TotPrec",
                    "TotPrecF",
                    "TotPrecEF",
                    "TotIntF",
                    "TotRain",
                    "TotRainF",
                    "TotETpot",
                    "TotETpotF",
                    "TotETact",
                    "TotETactF",
                    "TotRootRF",
                    "TotRootDF",
                    "TotRootPF",
                    "TotSubPF",
                    "TotCapRF",
                    "TotGwRechargeF",
                    "TotRainRF",
                    "TotBaseRF",
                    "TotRF",
                    "RainRAtot",
                    "BaseRAtot",
                    "QallRAtot",
                ]
                # -set time-series reporting for mm fluxes from upstream area if groundwater and routing/reservoir modules are used
                if self.mm_rep_FLAG == 1 and (
                    self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1
                ):
                    self.QRAINSubBasinTSS = TimeoutputTimeseries(
                        "QRAINSubBasinTSS", self, self.Locations, noHeader=False
                    )
                    self.QBASFSubBasinTSS = TimeoutputTimeseries(
                        "QBASFSubBasinTSS", self, self.Locations, noHeader=False
                    )
                    self.QTOTSubBasinTSS = TimeoutputTimeseries(
                        "QTOTSubBasinTSS", self, self.Locations, noHeader=False
                    )
            else:
                pars = [
                    "wbal",
                    "GWL",
                    "TotPrec",
                    "TotPrecF",
                    "TotPrecEF",
                    "TotIntF",
                    "TotRain",
                    "TotRainF",
                    "TotETpot",
                    "TotETpotF",
                    "TotETact",
                    "TotETactF",
                    "TotRootRF",
                    "TotRootDF",
                    "TotRootPF",
                    "TotSubDF",
                    "TotCapRF",
                    "TotSeepF",
                    "TotRainRF",
                    "TotRF",
                    "RainRAtot",
                    "BaseRAtot",
                    "QallRAtot",
                ]
                # -set time-series reporting for mm fluxes from upstream area if routing/reservoir modules are used
                if self.mm_rep_FLAG == 1 and (
                    self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1
                ):
                    self.SeepSubBasinTSS = TimeoutputTimeseries(
                        "SeepSubBasinTSS", self, self.Locations, noHeader=False
                    )
                    self.QRAINSubBasinTSS = TimeoutputTimeseries(
                        "QRAINSubBasinTSS", self, self.Locations, noHeader=False
                    )
                    self.QBASFSubBasinTSS = TimeoutputTimeseries(
                        "QBASFSubBasinTSS", self, self.Locations, noHeader=False
                    )
                    self.QTOTSubBasinTSS = TimeoutputTimeseries(
                        "QTOTSubBasinTSS", self, self.Locations, noHeader=False
                    )
        # -remove routing output from reported list of parameters if these modules are not used
        if self.RoutFLAG == 0 and self.ResFLAG == 0 and self.LakeFLAG == 0:
            rpars = ["RainRAtot", "SnowRAtot", "GlacRAtot", "BaseRAtot", "QallRAtot"]
            for i in rpars:
                try:
                    j = pars.index(i)
                    del pars[j]
                except:
                    pass
        # -set reporting options and read initial values
        for i in pars:
            mapoutops = configget(
                self.config, "REPORTING", i + "_mapoutput", i + "_mapoutput"
            )
            # mapoutops = config.get('REPORTING', i+'_mapoutput')
            TSoutops = configget(
                self.config, "REPORTING", i + "_TSoutput", i + "_TSoutput"
            )
            # TSoutops = config.get('REPORTING', i+'_TSoutput')
            if mapoutops == "NONE" and TSoutops == "NONE":
                print i + " will NOT be reported"
            else:
                print i + " will be reported"
                fname = configget(self.config, "REPORTING", i + "_fname", i + "_fname")
                # fname = config.get('REPORTING', i+'_fname')
                setattr(self, i + "_fname", fname)
                # 				try:
                # 					setattr(self, i, readmap(self.inpath + config.get('INITTOT', i)))
                # 				except:
                # 					try:
                # 						setattr(self, i, config.getfloat('INITTOT', i))
                # 					except:
                # 						setattr(self, i, 0.)
                setattr(
                    self, i, 0.
                )  # use this instead of the commented part above, because it is more logical to always zero as initial condition for reporting
                if mapoutops != "NONE":
                    mapoutops = mapoutops.split(",")
                    for j in mapoutops:
                        if j == "D":
                            setattr(self, i + "_Day", eval("self." + i))
                            setattr(self, i + "_Day_map", 1)
                        elif j == "M":
                            setattr(self, i + "_Month", eval("self." + i))
                            setattr(self, i + "_Month_map", 1)
                        elif j == "Y":
                            setattr(self, i + "_Year", eval("self." + i))
                            setattr(self, i + "_Year_map", 1)
                        else:
                            setattr(self, i + "_Final", eval("self." + i))
                            setattr(self, i + "_Final_map", 1)
                if TSoutops != "NONE":
                    TSoutops = TSoutops.split(",")
                    for j in TSoutops:
                        if j == "D":
                            setattr(self, i + "_Day", eval("self." + i))
                            setattr(
                                self,
                                i + "_DayTS",
                                eval(
                                    'TimeoutputTimeseries("'
                                    + fname
                                    + "DTS"
                                    + '", self, self.Locations, noHeader=False)'
                                ),
                            )
                        elif j == "M":
                            setattr(self, i + "_Month", eval("self." + i))
                            setattr(
                                self,
                                i + "_MonthTS",
                                eval(
                                    'TimeoutputTimeseries("'
                                    + fname
                                    + "MTS"
                                    + '", self, self.Locations, noHeader=False)'
                                ),
                            )
                        elif j == "Y":
                            setattr(self, i + "_Year", eval("self." + i))
                            setattr(
                                self,
                                i + "_YearTS",
                                eval(
                                    'TimeoutputTimeseries("'
                                    + fname
                                    + "YTS"
                                    + '", self, self.Locations, noHeader=False)'
                                ),
                            )

        # -set reporting of water balances for lakes
        if self.LakeFLAG == 1 and config.getint("REPORTING", "Lake_wbal") == 1:
            self.LakeInTSS = pcrm.TimeoutputTimeseries(
                "LakeInTSS", self, self.LakeID, noHeader=True
            )
            self.LakeOutTSS = pcrm.TimeoutputTimeseries(
                "LakeOutTSS", self, self.LakeID, noHeader=True
            )
            self.LakeStorTSS = pcrm.TimeoutputTimeseries(
                "LakeStorTSS", self, self.LakeID, noHeader=True
            )
            if (
                self.RainRA_FLAG == 1
            ):  # -set reporting of water balances for individual components
                self.LakeRainInTSS = pcrm.TimeoutputTimeseries(
                    "LakeRainInTSS", self, self.LakeID, noHeader=True
                )
                self.LakeRainOutTSS = pcrm.TimeoutputTimeseries(
                    "LakeRainOutTSS", self, self.LakeID, noHeader=True
                )
                self.LakeRainStorTSS = pcrm.TimeoutputTimeseries(
                    "LakeRainStorTSS", self, self.LakeID, noHeader=True
                )
            if self.SnowRA_FLAG == 1:
                self.LakeSnowInTSS = pcrm.TimeoutputTimeseries(
                    "LakeSnowInTSS", self, self.LakeID, noHeader=True
                )
                self.LakeSnowOutTSS = pcrm.TimeoutputTimeseries(
                    "LakeSnowOutTSS", self, self.LakeID, noHeader=True
                )
                self.LakeSnowStorTSS = pcrm.TimeoutputTimeseries(
                    "LakeSnowStorTSS", self, self.LakeID, noHeader=True
                )
            if self.GlacRA_FLAG == 1:
                self.LakeGlacInTSS = pcrm.TimeoutputTimeseries(
                    "LakeGlacInTSS", self, self.LakeID, noHeader=True
                )
                self.LakeGlacOutTSS = pcrm.TimeoutputTimeseries(
                    "LakeGlacOutTSS", self, self.LakeID, noHeader=True
                )
                self.LakeGlacStorTSS = pcrm.TimeoutputTimeseries(
                    "LakeGlacStorTSS", self, self.LakeID, noHeader=True
                )
            if self.BaseRA_FLAG == 1:
                self.LakeBaseInTSS = pcrm.TimeoutputTimeseries(
                    "LakeBaseInTSS", self, self.LakeID, noHeader=True
                )
                self.LakeBaseOutTSS = pcrm.TimeoutputTimeseries(
                    "LakeBaseOutTSS", self, self.LakeID, noHeader=True
                )
                self.LakeBaseStorTSS = pcrm.TimeoutputTimeseries(
                    "LakeBaseStorTSS", self, self.LakeID, noHeader=True
                )
        # -set reporting of water balances for reservoirs
        if self.ResFLAG == 1 and config.getint("REPORTING", "Res_wbal") == 1:
            self.ResInTSS = pcrm.TimeoutputTimeseries(
                "ResInTSS", self, self.ResID, noHeader=True
            )
            self.ResOutTSS = pcrm.TimeoutputTimeseries(
                "ResOutTSS", self, self.ResID, noHeader=True
            )
            self.ResStorTSS = pcrm.TimeoutputTimeseries(
                "ResStorTSS", self, self.ResID, noHeader=True
            )
            if (
                self.RainRA_FLAG == 1
            ):  # -set reporting of water balances for individual components
                self.ResRainInTSS = pcrm.TimeoutputTimeseries(
                    "ResRainInTSS", self, self.ResID, noHeader=True
                )
                self.ResRainOutTSS = pcrm.TimeoutputTimeseries(
                    "ResRainOutTSS", self, self.ResID, noHeader=True
                )
                self.ResRainStorTSS = pcrm.TimeoutputTimeseries(
                    "ResRainStorTSS", self, self.ResID, noHeader=True
                )
            if self.SnowRA_FLAG == 1:
                self.ResSnowInTSS = pcrm.TimeoutputTimeseries(
                    "ResSnowInTSS", self, self.ResID, noHeader=True
                )
                self.ResSnowOutTSS = pcrm.TimeoutputTimeseries(
                    "ResSnowOutTSS", self, self.ResID, noHeader=True
                )
                self.ResSnowStorTSS = pcrm.TimeoutputTimeseries(
                    "ResSnowStorTSS", self, self.ResID, noHeader=True
                )
            if self.GlacRA_FLAG == 1:
                self.ResGlacInTSS = pcrm.TimeoutputTimeseries(
                    "ResGlacInTSS", self, self.ResID, noHeader=True
                )
                self.ResGlacOutTSS = pcrm.TimeoutputTimeseries(
                    "ResGlacOutTSS", self, self.ResID, noHeader=True
                )
                self.ResGlacStorTSS = pcrm.TimeoutputTimeseries(
                    "ResGlacStorTSS", self, self.ResID, noHeader=True
                )
            if self.BaseRA_FLAG == 1:
                self.ResBaseInTSS = pcrm.TimeoutputTimeseries(
                    "ResBaseInTSS", self, self.ResID, noHeader=True
                )
                self.ResBaseOutTSS = pcrm.TimeoutputTimeseries(
                    "ResBaseOutTSS", self, self.ResID, noHeader=True
                )
                self.ResBaseStorTSS = pcrm.TimeoutputTimeseries(
                    "ResBaseStorTSS", self, self.ResID, noHeader=True
                )

        # if self.scalarInput:
        # self.gaugesMap=self.wf_readmap(os.path.join(self.Dir , wflow_mgauges),0.0,fail=True) #: Map with locations of rainfall/evap/temp gauge(s). Only needed if the input to the model is not in maps
        # self.OutputId=self.wf_readmap(os.path.join(self.Dir , wflow_subcatch),0.0,fail=True)       # location of subcatchment

        self.ZeroMap = 0.0 * scalar(defined(self.DEM))  # map with only zero's

        # For in memory override:
        # self.Prec, self.Tair, self.Tmax, self.Tmin = self.ZeroMap

        # Set static initial values here #########################################
        self.Latitude = ycoordinate(boolean(self.ZeroMap))
        self.Longitude = xcoordinate(boolean(self.ZeroMap))

        # self.logger.info("Linking parameters to landuse, catchment and soil...")

        # self.Beta = scalar(0.6) # For sheetflow
        # #self.M=lookupscalar(self.Dir + "/" + modelEnv['intbl'] + "/M.tbl" ,self.LandUse,subcatch,self.Soil) # Decay parameter in Topog_sbm
        # self.N=lookupscalar(self.Dir + "/" + self.intbl + "/N.tbl",self.LandUse,subcatch,self.Soil)  # Manning overland flow
        # """ *Parameter:* Manning's N for all non-river cells """
        # self.NRiver=lookupscalar(self.Dir + "/" + self.intbl + "/N_River.tbl",self.LandUse,subcatch,self.Soil)  # Manning river
        # """ Manning's N for all cells that are marked as a river """

        self.wf_updateparameters()

        # Multiply parameters with a factor (for calibration etc) -P option in command line

        self.wf_multparameters()

    def default_summarymaps(self):  ##-maybe not needed. check later
        """
      Returns a list of default summary-maps at the end of a run.
      This is model specific. You can also add them to the [summary]section of the ini file but stuff
      you think is crucial to the model should be listed here

       Example:

      """
        # lst = ['self.Cfmax','self.csize','self.upsize','self.TTI','self.TT','self.WHC',
        #       'self.Slope','self.N','self.xl','self.yl','self.reallength','self.DCL','self.Bw',]
        lst = ["self.GlacFrac"]

        return lst

    def resume(self):
        """ read initial state maps (they are output of a previous call to suspend()) """

        if self.reinit == 1:  # -to be defined for sphy model state variables!!!
            self.logger.info("Setting initial conditions to default (zero!)")
            self.RootWater = RootFieldMap
            self.SubWater = SubFieldMap
            self.CapRise = 3
            self.RootDrain = 3
            self.SubDrain = 3
            self.GwRecharge = 2
            self.BaseR = 1
            self.Gw = 1500
            self.H_gw = 3
            self.SnowIni = cover(0.0)
            self.SnowWatStore = cover(0.0)
            self.QRA_init = cover(0.0)
            self.RainRA_init = cover(0.0)
            self.BaseRA_init = cover(0.0)
            self.SnowRA_init = cover(0.0)
            self.GlacRA_init = cover(0.0)
            # self.FreeWater =  cover(0.0) #: Water on surface (state variable [mm])
            # self.SoilMoisture =  self.FC #: Soil moisture (state variable [mm])
            # self.UpperZoneStorage = 0.2 * self.FC #: Storage in Upper Zone (state variable [mm])
            # self.LowerZoneStorage = 1.0/(3.0 * self.K4) #: Storage in Uppe Zone (state variable [mm])
            # self.InterceptionStorage = cover(0.0) #: Interception Storage (state variable [mm])
            # self.SurfaceRunoff = cover(0.0) #: Discharge in kinimatic wave (state variable [m^3/s])
            # self.WaterLevel = cover(0.0) #: Water level in kinimatic wave (state variable [m])
            # self.DrySnow=cover(0.0) #: Snow amount (state variable [mm])
            # if hasattr(self, 'ReserVoirSimpleLocs'):
            #    self.ReservoirVolume = self.ResMaxVolume * self.ResTargetFullFrac
            # if hasattr(self, 'ReserVoirComplexLocs'):
            #    self.ReservoirWaterLevel = cover(0.0)
        else:
            self.wf_resume(os.path.join(self.Dir, "instate"))

    def dynamic(self):

        self.wf_updateparameters()  # read forcing an dynamic parameters
        # self.counter+=1

        # print str(self.curdate.day)+'-'+str(self.curdate.month)+'-'+str(self.curdate.year)+'  t = '+str(self.counter)

        # Snow and glacier fraction settings
        if self.GlacFLAG == 0:
            self.GlacFrac = 0
        if self.SnowFLAG == 0:
            self.SnowStore = scalar(0)
        SnowFrac = ifthenelse(self.SnowStore > 0, scalar(1 - self.GlacFrac), 0)
        RainFrac = ifthenelse(self.SnowStore == 0, scalar(1 - self.GlacFrac), 0)

        # -Read the precipitation time-series
        Precip = self.Prec
        # -Report Precip
        self.reporting.reporting(self, pcr, "TotPrec", Precip)
        self.reporting.reporting(self, pcr, "TotPrecF", Precip * (1 - self.GlacFrac))

        # -Temperature and reference evapotranspiration
        Temp = self.Tair
        if self.ETREF_FLAG == 0:
            TempMax = self.Tmax
            TempMin = self.Tmin
            ETref = self.Hargreaves.Hargreaves(
                pcr, self.Hargreaves.extrarad(self, pcr), Temp, TempMax, TempMin
            )
        else:
            ETref = readmap(generateNameT(self.ETref, self.counter))

            # -Interception and effective precipitation
            # -Update canopy storage
        if self.DynVegFLAG == 1:
            # -try to read the ndvi map series. If not available, then use ndvi old
            try:
                ndvi = readmap(pcrm.generateNameT(self.ndvi, self.counter))
                self.ndviOld = ndvi
            except:
                ndvi = self.ndviOld
                # -fill missing ndvi values with ndvi base
            ndvi = ifthenelse(defined(ndvi) == 1, ndvi, self.NDVIbase)
            # -calculate the vegetation parameters
            vegoutput = self.dynamic_veg.Veg_function(
                pcr,
                ndvi,
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
            self.Scanopy = self.Scanopy + Precip
            # -interception and effective precipitation
            intercep = self.dynamic_veg.Inter_function(
                pcr, self.Scanopy, vegoutput[1], ETref
            )
            # -interception
            Int = intercep[0]
            # -report interception corrected for fraction
            self.reporting.reporting(self, pcr, "TotIntF", Int * (1 - self.GlacFrac))
            # -effective precipitation
            Precip = intercep[1]
            # -Report effective precipitation corrected for fraction
            self.reporting.reporting(
                self, pcr, "TotPrecEF", Precip * (1 - self.GlacFrac)
            )
            # -canopy storage
            self.Scanopy = intercep[2]
        elif self.KcStatFLAG == 0:
            # -Try to read the KC map series
            try:
                self.Kc = readmap(pcrm.generateNameT(self.Kcmaps, self.counter))
                self.KcOld = self.Kc
            except:
                self.Kc = self.KcOld
                # -report mm effective precipitation for sub-basin averages
        if self.mm_rep_FLAG == 1 and (
            self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1
        ):
            self.PrecSubBasinTSS.sample(
                catchmenttotal(Precip * (1 - self.GlacFrac), self.FlowDir)
                / catchmenttotal(1, self.FlowDir)
            )

            # Snow and rain
        if self.SnowFLAG == 1:
            # -Snow and rain differentiation
            Snow = ifthenelse(Temp >= self.Tcrit, 0, Precip)
            Rain = ifthenelse(Temp < self.Tcrit, 0, Precip)
            # -Report Snow
            self.reporting.reporting(self, pcr, "TotSnow", Snow)
            self.reporting.reporting(self, pcr, "TotSnowF", Snow * (1 - self.GlacFrac))
            # -Snow melt
            PotSnowMelt = self.snow.PotSnowMelt(pcr, Temp, self.DDFS)
            ActSnowMelt = self.snow.ActSnowMelt(pcr, self.SnowStore, PotSnowMelt)
            # -Report snow melt
            self.reporting.reporting(self, pcr, "TotSnowMelt", ActSnowMelt)
            self.reporting.reporting(self, pcr, "TotSnowMeltF", ActSnowMelt * SnowFrac)
            # -Update snow store
            self.SnowStore = self.snow.SnowStoreUpdate(
                pcr, self.SnowStore, Snow, ActSnowMelt, Temp, self.SnowWatStore
            )
            # -Caclulate the maximum amount of water that can be stored in snowwatstore
            MaxSnowWatStore = self.snow.MaxSnowWatStorage(self.SnowSC, self.SnowStore)
            OldSnowWatStore = self.SnowWatStore
            # -Calculate the actual amount of water stored in snowwatstore
            self.SnowWatStore = self.snow.SnowWatStorage(
                pcr, Temp, MaxSnowWatStore, self.SnowWatStore, ActSnowMelt, Rain
            )
            # -Changes in total water storage in snow (SnowStore and SnowWatStore)
            OldTotalSnowStore = self.TotalSnowStore
            self.TotalSnowStore = self.snow.TotSnowStorage(
                self.SnowStore, self.SnowWatStore, SnowFrac, RainFrac
            )
            # -Snow runoff
            SnowR = self.snow.SnowR(
                pcr,
                self.SnowWatStore,
                MaxSnowWatStore,
                ActSnowMelt,
                Rain,
                OldSnowWatStore,
                SnowFrac,
            )
            # -Report Snow runoff
            self.reporting.reporting(self, pcr, "TotSnowRF", SnowR)
        else:
            Rain = Precip
            SnowR = 0
            OldTotalSnowStore = 0
            self.TotalSnowStore = 0
            # -Report Rain
        self.reporting.reporting(self, pcr, "TotRain", Rain)
        self.reporting.reporting(self, pcr, "TotRainF", Rain * (1 - self.GlacFrac))

        # -Glacier calculations
        if self.GlacFLAG == 1:
            # -Glacier melt from clean ice glaciers
            GlacCIMelt = self.glacier.GlacCDMelt(pcr, Temp, self.DDFG, self.GlacFracCI)
            # -Glacier melt from debris covered glaciers
            GlacDCMelt = self.glacier.GlacCDMelt(pcr, Temp, self.DDFDG, self.GlacFracDB)
            # -Total melt from glaciers
            GlacMelt = self.glacier.GMelt(GlacCIMelt, GlacDCMelt)
            self.GlacMelt = GlacMelt
            # -Report glacier melt
            self.reporting.reporting(self, pcr, "TotGlacMelt", GlacMelt)
            self.reporting.reporting(
                self, pcr, "TotGlacMeltF", GlacMelt * self.GlacFrac
            )
            if self.mm_rep_FLAG == 1 and (
                self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1
            ):
                self.GMeltSubBasinTSS.sample(
                    catchmenttotal(GlacMelt * self.GlacFrac, self.FlowDir)
                    / catchmenttotal(1, self.FlowDir)
                )
                # -Glacier runoff
            GlacR = self.glacier.GlacR(self.GlacF, GlacMelt, self.GlacFrac)
            # -Report glacier runoff
            self.reporting.reporting(self, pcr, "TotGlacRF", GlacR)
            # -Glacier percolation to groundwater
            GlacPerc = self.glacier.GPerc(self.GlacF, GlacMelt, self.GlacFrac)
            # -Report glacier percolation to groundwater
            self.reporting.reporting(self, pcr, "TotGlacPercF", GlacPerc)
        else:
            GlacR = 0
            GlacMelt = 0
            GlacPerc = 0

            # -Potential evapotranspiration (THIS SHOULD STILL BE IMPROVED WITH DYNAMIC VEGETATION MODULE)
        ETpot = self.ET.ETpot(ETref, self.Kc)
        # -Report ETpot
        self.reporting.reporting(self, pcr, "TotETpot", ETpot)
        self.reporting.reporting(self, pcr, "TotETpotF", ETpot * RainFrac)

        # -Rootzone calculations
        self.RootWater = (
            self.RootWater + ifthenelse(RainFrac > 0, Rain, 0) + self.CapRise
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
        ETact = self.ET.ETact(
            pcr, ETpot, self.RootWater, self.RootSat, etreddry, RainFrac
        )
        # -Report the actual evapotranspiration
        self.reporting.reporting(self, pcr, "TotETact", ETact)
        # -Actual evapotranspiration, corrected for rain fraction
        ActETact = ETact * RainFrac
        # -Report the actual evapotranspiration, corrected for rain fraction
        self.reporting.reporting(self, pcr, "TotETactF", ActETact)
        if self.mm_rep_FLAG == 1 and (
            self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1
        ):
            self.ETaSubBasinTSS.sample(
                catchmenttotal(ActETact, self.FlowDir) / catchmenttotal(1, self.FlowDir)
            )
            # -Update rootwater content
        self.RootWater = max(self.RootWater - ETact, 0)
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
        rootperc = self.rootzone.RootPercolation(
            pcr, self.RootWater, self.SubWater, self.RootField, self.RootTT, self.SubSat
        )
        # -Report rootzone percolation, corrected for fraction
        self.reporting.reporting(self, pcr, "TotRootPF", rootperc * (1 - self.GlacFrac))
        # -Update rootwater content
        self.RootWater = self.RootWater - rootperc

        # -Sub soil calculations
        self.SubWater = self.SubWater + rootperc
        if self.GroundFLAG == 0:
            if self.SeepStatFLAG == 0:
                try:
                    self.SeePage = readmap(
                        pcrm.generateNameT(self.Seepmaps, self.counter)
                    )
                    self.SeepOld = self.SeePage
                except:
                    self.SeePage = self.SeepOld
                    # -Report seepage
            self.reporting.reporting(self, pcr, "TotSeepF", scalar(self.SeePage))
            self.SubWater = min(max(self.SubWater - self.SeePage, 0), self.SubSat)
            if self.mm_rep_FLAG == 1 and (
                self.RoutFLAG == 1 or self.ResFLAG == 1 or self.LakeFLAG == 1
            ):
                self.SeepSubBasinTSS.sample(
                    catchmenttotal(self.SeePage, self.FlowDir)
                    / catchmenttotal(1, self.FlowDir)
                )
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
        # -Report capillary rise, corrected for fraction
        self.reporting.reporting(
            self, pcr, "TotCapRF", self.CapRise * (1 - self.GlacFrac)
        )
        # -Update sub soil water content
        self.SubWater = self.SubWater - self.CapRise
        if (
            self.GroundFLAG == 1
        ):  # sub percolation will be calculated instead of subdrainage
            subperc = self.subzone.SubPercolation(
                pcr, self.SubWater, self.SubField, self.SubTT, self.Gw, self.GwSat
            )
            ActSubPerc = subperc * (1 - self.GlacFrac)
            # -Report the subzone percolation, corrected for the fraction
            self.reporting.reporting(self, pcr, "TotSubPF", ActSubPerc)
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
            # -Report drainage from subzone
            self.reporting.reporting(self, pcr, "TotSubDF", self.SubDrain)
            # -Update sub soil water content
            self.SubWater = self.SubWater - self.SubDrain

            # -Changes in soil water storage
        OldSoilWater = self.SoilWater
        self.SoilWater = (self.RootWater + self.SubWater) * (1 - self.GlacFrac)

        # -Rootzone runoff
        RootR = RootRunoff * RainFrac
        # -Report rootzone runoff, corrected for fraction
        self.reporting.reporting(self, pcr, "TotRootRF", RootR)
        # -Rootzone drainage
        RootD = self.RootDrain * (1 - self.GlacFrac)
        # -Report rootzone drainage, corrected for fraction
        self.reporting.reporting(self, pcr, "TotRootDF", RootD)
        # -Rain runoff
        RainR = RootR + RootD
        # -Report rain runoff
        self.reporting.reporting(self, pcr, "TotRainRF", RainR)

        # -Groundwater calculations
        if self.GroundFLAG == 1:
            GwOld = self.Gw
            # -Groundwater recharge
            self.GwRecharge = self.groundwater.GroundWaterRecharge(
                pcr, self.deltaGw, self.GwRecharge, ActSubPerc, GlacPerc
            )
            # -Report groundwater recharge
            self.reporting.reporting(self, pcr, "TotGwRechargeF", self.GwRecharge)
            # -Update groundwater storage
            self.Gw = self.Gw + self.GwRecharge
            # -Baseflow
            self.BaseR = self.groundwater.BaseFlow(
                pcr, self.Gw, self.BaseR, self.GwRecharge, self.BaseThresh, self.alphaGw
            )
            # -Report Baseflow
            self.reporting.reporting(self, pcr, "TotBaseRF", self.BaseR)
            # -Update groundwater storage
            self.Gw = self.Gw - self.BaseR
            # -Calculate groundwater level
            self.H_gw = self.groundwater.HLevel(
                pcr, self.H_gw, self.alphaGw, self.GwRecharge, self.YieldGw
            )
            # -Report groundwater
            self.reporting.reporting(
                self,
                pcr,
                "GWL",
                (
                    (self.SubDepthFlat + self.RootDepthFlat + self.GwDepth) / 1000
                    - self.H_gw
                )
                * -1,
            )

        else:
            # -Use drainage from subsoil as baseflow
            self.BaseR = self.SubDrain
            # -Groundwater level as scaled between min and max measured gwl
            SoilAct = self.RootWater + self.SubWater
            SoilRel = (SoilAct - self.SoilMin) / (
                self.SoilMax - self.SoilMin
            )  # scale between 0 (dry) and 1 (wet)
            GWL = self.GWL_base - (SoilRel - 0.5) * self.GWL_base
            # -Report groundwater
            self.reporting.reporting(self, pcr, "GWL", GWL)

            # -Report Total runoff
        self.reporting.reporting(self, pcr, "TotRF", self.BaseR + RainR + SnowR + GlacR)

        # -Water balance
        if self.GroundFLAG == 1:
            waterbalance = (
                Precip * (1 - self.GlacFrac)
                + GlacMelt * self.GlacFrac
                - ActETact
                - GlacR
                - SnowR
                - RainR
                - self.BaseR
                - (self.SoilWater - OldSoilWater)
                - (self.TotalSnowStore - OldTotalSnowStore)
                - (self.Gw - GwOld)
            )
        elif self.GroundFLAG == 0:
            waterbalance = (
                Precip
                - ActETact
                - self.SeePage
                - SnowR
                - RainR
                - self.BaseR
                - (self.SoilWater - OldSoilWater)
                - (self.TotalSnowStore - OldTotalSnowStore)
            )
        self.reporting.reporting(self, pcr, "wbal", waterbalance)

        # -Routing for lake and/or reservoir modules
        if self.LakeFLAG == 1 or self.ResFLAG == 1:
            # -Update storage in lakes/reservoirs (m3) with specific runoff
            self.StorRES = self.StorRES + ifthenelse(
                self.QFRAC == 0,
                0.001 * cellarea() * (self.BaseR + RainR + GlacR + SnowR),
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
                Qout = ifthenelse(
                    self.ResID != 0, ResQ, ifthenelse(self.LakeID != 0, LakeQ, 0)
                )
            elif self.LakeFLAG == 1:
                tempvar = self.lakes.UpdateLakeHStore(self, pcr, pcrm)
                LakeLevel = tempvar[0]
                self.StorRES = tempvar[1]
                Qout = self.lakes.QLake(self, pcr, LakeLevel)
            else:
                Qout = self.reservoirs.QRes(self, pcr)

                # -Calculate volume available for routing (=outflow lakes/reservoir + cell specific runoff)
            RunoffVolume = upstream(self.FlowDir, Qout) + ifthenelse(
                self.QFRAC == 0,
                0,
                0.001 * cellarea() * (self.BaseR + RainR + GlacR + SnowR),
            )
            # -Routing of total flow
            tempvar = self.routing.ROUT(
                self, pcr, RunoffVolume, self.QRAold, Qout, self.StorRES
            )
            self.StorRES = tempvar[0]
            Q = tempvar[1]
            Qin = tempvar[2]
            self.QRAold = Q
            self.reporting.reporting(self, pcr, "QallRAtot", Q)
            # -report flux in mm
            if self.mm_rep_FLAG == 1:
                self.QTOTSubBasinTSS.sample(
                    ((Q * 3600 * 24) / catchmenttotal(cellarea(), self.FlowDir)) * 1000
                )
                # -report lake and reservoir waterbalance
            if self.LakeFLAG == 1 and config.getint("REPORTING", "Lake_wbal") == 1:
                self.LakeInTSS.sample(Qin)
                self.LakeOutTSS.sample(Qout)
                self.LakeStorTSS.sample(self.StorRES)
            if self.ResFLAG == 1 and config.getint("REPORTING", "Res_wbal") == 1:
                self.ResInTSS.sample(Qin)
                self.ResOutTSS.sample(Qout)
                self.ResStorTSS.sample(self.StorRES)

                # -Routing of individual contributers
                # -Snow routing
            if self.SnowRA_FLAG == 1 and self.SnowFLAG == 1:
                self.SnowRAstor = self.SnowRAstor + ifthenelse(
                    self.QFRAC == 0, SnowR * 0.001 * cellarea(), 0
                )
                cQfrac = cover(self.SnowRAstor / OldStorage, 0)
                cQout = cQfrac * Qout
                cRunoffVolume = upstream(self.FlowDir, cQout) + ifthenelse(
                    self.QFRAC == 0, 0, 0.001 * cellarea() * SnowR
                )
                tempvar = self.routing.ROUT(
                    self, pcr, cRunoffVolume, self.SnowRAold, cQout, self.SnowRAstor
                )
                self.SnowRAstor = tempvar[0]
                SnowRA = tempvar[1]
                cQin = tempvar[2]
                self.SnowRAold = SnowRA
                self.reporting.reporting(self, pcr, "SnowRAtot", SnowRA)
                if self.mm_rep_FLAG == 1:
                    self.QSNOWSubBasinTSS.sample(
                        (
                            (SnowRA * 3600 * 24)
                            / catchmenttotal(cellarea(), self.FlowDir)
                        )
                        * 1000
                    )
                    # -report lake and reservoir waterbalance
                if self.LakeFLAG == 1 and config.getint("REPORTING", "Lake_wbal") == 1:
                    self.LakeSnowInTSS.sample(cQin)
                    self.LakeSnowOutTSS.sample(cQout)
                    self.LakeSnowStorTSS.sample(self.SnowRAstor)
                if self.ResFLAG == 1 and config.getint("REPORTING", "Res_wbal") == 1:
                    self.ResSnowInTSS.sample(cQin)
                    self.ResSnowOutTSS.sample(cQout)
                    self.ResSnowStorTSS.sample(self.SnowRAstor)
                    # -Rain routing
            if self.RainRA_FLAG == 1:
                self.RainRAstor = self.RainRAstor + ifthenelse(
                    self.QFRAC == 0, RainR * 0.001 * cellarea(), 0
                )
                cQfrac = cover(self.RainRAstor / OldStorage, 0)
                cQout = cQfrac * Qout
                cRunoffVolume = upstream(self.FlowDir, cQout) + ifthenelse(
                    self.QFRAC == 0, 0, 0.001 * cellarea() * RainR
                )
                tempvar = self.routing.ROUT(
                    self, pcr, cRunoffVolume, self.RainRAold, cQout, self.RainRAstor
                )
                self.RainRAstor = tempvar[0]
                RainRA = tempvar[1]
                cQin = tempvar[2]
                self.RainRAold = RainRA
                self.reporting.reporting(self, pcr, "RainRAtot", RainRA)
                if self.mm_rep_FLAG == 1:
                    self.QRAINSubBasinTSS.sample(
                        (
                            (RainRA * 3600 * 24)
                            / catchmenttotal(cellarea(), self.FlowDir)
                        )
                        * 1000
                    )
                    # -report lake and reservoir waterbalance
                if self.LakeFLAG == 1 and config.getint("REPORTING", "Lake_wbal") == 1:
                    self.LakeRainInTSS.sample(cQin)
                    self.LakeRainOutTSS.sample(cQout)
                    self.LakeRainStorTSS.sample(self.RainRAstor)
                if self.ResFLAG == 1 and config.getint("REPORTING", "Res_wbal") == 1:
                    self.ResRainInTSS.sample(cQin)
                    self.ResRainOutTSS.sample(cQout)
                    self.ResRainStorTSS.sample(self.RainRAstor)
                    # -Glacier routing
            if self.GlacRA_FLAG == 1 and self.GlacFLAG == 1:
                self.GlacRAstor = self.GlacRAstor + ifthenelse(
                    self.QFRAC == 0, GlacR * 0.001 * cellarea(), 0
                )
                cQfrac = cover(self.GlacRAstor / OldStorage, 0)
                cQout = cQfrac * Qout
                cRunoffVolume = upstream(self.FlowDir, cQout) + ifthenelse(
                    self.QFRAC == 0, 0, 0.001 * cellarea() * GlacR
                )
                tempvar = self.routing.ROUT(
                    self, pcr, cRunoffVolume, self.GlacRAold, cQout, self.GlacRAstor
                )
                self.GlacRAstor = tempvar[0]
                GlacRA = tempvar[1]
                cQin = tempvar[2]
                self.GlacRAold = GlacRA
                self.reporting.reporting(self, pcr, "GlacRAtot", GlacRA)
                if self.mm_rep_FLAG == 1:
                    self.QGLACSubBasinTSS.sample(
                        (
                            (GlacRA * 3600 * 24)
                            / catchmenttotal(cellarea(), self.FlowDir)
                        )
                        * 1000
                    )
                    # -report lake and reservoir waterbalance
                if self.LakeFLAG == 1 and config.getint("REPORTING", "Lake_wbal") == 1:
                    self.LakeGlacInTSS.sample(cQin)
                    self.LakeGlacOutTSS.sample(cQout)
                    self.LakeGlacStorTSS.sample(self.GlacRAstor)
                if self.ResFLAG == 1 and config.getint("REPORTING", "Res_wbal") == 1:
                    self.ResGlacInTSS.sample(cQin)
                    self.ResGlacOutTSS.sample(cQout)
                    self.ResGlacStorTSS.sample(self.GlacRAstor)
                    # -Baseflow routing
            if self.BaseRA_FLAG == 1:
                self.BaseRAstor = self.BaseRAstor + ifthenelse(
                    self.QFRAC == 0, self.BaseR * 0.001 * cellarea(), 0
                )
                cQfrac = cover(self.BaseRAstor / OldStorage, 0)
                cQout = cQfrac * Qout
                cRunoffVolume = upstream(self.FlowDir, cQout) + ifthenelse(
                    self.QFRAC == 0, 0, 0.001 * cellarea() * self.BaseR
                )
                tempvar = self.routing.ROUT(
                    self, pcr, cRunoffVolume, self.BaseRAold, cQout, self.BaseRAstor
                )
                self.BaseRAstor = tempvar[0]
                BaseRA = tempvar[1]
                cQin = tempvar[2]
                self.BaseRAold = BaseRA
                self.reporting.reporting(self, pcr, "BaseRAtot", BaseRA)
                if self.mm_rep_FLAG == 1:
                    self.QBASFSubBasinTSS.sample(
                        (
                            (BaseRA * 3600 * 24)
                            / catchmenttotal(cellarea(), self.FlowDir)
                        )
                        * 1000
                    )
                    # -report lake and reservoir waterbalance
                if self.LakeFLAG == 1 and config.getint("REPORTING", "Lake_wbal") == 1:
                    self.LakeBaseInTSS.sample(cQin)
                    self.LakeBaseOutTSS.sample(cQout)
                    self.LakeBaseStorTSS.sample(self.BaseRAstor)
                if self.ResFLAG == 1 and config.getint("REPORTING", "Res_wbal") == 1:
                    self.ResBaseInTSS.sample(cQin)
                    self.ResBaseOutTSS.sample(cQout)
                    self.ResBaseStorTSS.sample(self.BaseRAstor)

                    # -Normal routing module
        elif self.RoutFLAG == 1:
            # -Rout total runoff
            Q = self.routing.ROUT(
                pcr,
                self.BaseR + RainR + GlacR + SnowR,
                self.QRAold,
                self.FlowDir,
                self.kx,
            )
            self.QRAold = Q
            self.reporting.reporting(self, pcr, "QallRAtot", Q)
            if self.mm_rep_FLAG == 1:
                self.QTOTSubBasinTSS.sample(
                    ((Q * 3600 * 24) / catchmenttotal(cellarea(), self.FlowDir)) * 1000
                )
                # -Snow routing
            if self.SnowRA_FLAG == 1 and self.SnowFLAG == 1:
                SnowRA = self.routing.ROUT(
                    pcr, SnowR, self.SnowRAold, self.FlowDir, self.kx
                )
                self.SnowRAold = SnowRA
                self.reporting.reporting(self, pcr, "SnowRAtot", SnowRA)
                if self.mm_rep_FLAG == 1:
                    self.QSNOWSubBasinTSS.sample(
                        (
                            (SnowRA * 3600 * 24)
                            / catchmenttotal(cellarea(), self.FlowDir)
                        )
                        * 1000
                    )
                    # -Rain routing
            if self.RainRA_FLAG == 1:
                RainRA = self.routing.ROUT(
                    pcr, RainR, self.RainRAold, self.FlowDir, self.kx
                )
                self.RainRAold = RainRA
                self.reporting.reporting(self, pcr, "RainRAtot", RainRA)
                if self.mm_rep_FLAG == 1:
                    self.QRAINSubBasinTSS.sample(
                        (
                            (RainRA * 3600 * 24)
                            / catchmenttotal(cellarea(), self.FlowDir)
                        )
                        * 1000
                    )
                    # -Glacier routing
            if self.GlacRA_FLAG == 1 and self.GlacFLAG == 1:
                GlacRA = self.routing.ROUT(
                    pcr, GlacR, self.GlacRAold, self.FlowDir, self.kx
                )
                self.GlacRAold = GlacRA
                self.reporting.reporting(self, pcr, "GlacRAtot", GlacRA)
                if self.mm_rep_FLAG == 1:
                    self.QGLACSubBasinTSS.sample(
                        (
                            (GlacRA * 3600 * 24)
                            / catchmenttotal(cellarea(), self.FlowDir)
                        )
                        * 1000
                    )
                    # -Baseflow routing
            if self.BaseRA_FLAG == 1:
                BaseRA = self.routing.ROUT(
                    pcr, self.BaseR, self.BaseRAold, self.FlowDir, self.kx
                )
                self.BaseRAold = BaseRA
                self.reporting.reporting(self, pcr, "BaseRAtot", BaseRA)
                if self.mm_rep_FLAG == 1:
                    self.QBASFSubBasinTSS.sample(
                        (
                            (BaseRA * 3600 * 24)
                            / catchmenttotal(cellarea(), self.FlowDir)
                        )
                        * 1000
                    )


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
        opts, args = getopt.getopt(argv, "c:QXS:F:hC:Ii:T:R:u:s:P:p:Xx:U:fl:L:")
    except getopt.error, msg:
        pcrut.usage(msg)

    for o, a in opts:
        if o == "-F":
            runinfoFile = a
            fewsrun = True

        if o == "-C":
            caseName = a
        if o == "-R":
            runId = a
        if o == "-L":
            LogFileName = a
        if o == "-l":
            exec "loglevel = logging." + a
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
            print "Failed to get timesteps from runinfo file: " + runinfoFile
            sys.exit(2)
    else:
        starttime = dt.datetime(1990, 01, 01)

    if _lastTimeStep < _firstTimeStep:
        print "The starttimestep (" + str(
            _firstTimeStep
        ) + ") is smaller than the last timestep (" + str(_lastTimeStep) + ")"
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
            exec "zz =" + a
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
