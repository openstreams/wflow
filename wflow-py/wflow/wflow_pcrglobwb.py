#!/usr/bin/python
# Wflow is Free software, see below:
#
# Copyright (c) J. Schellekens/Deltares 2005-2014
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
Run the wflow_pcrglobwb hydrological model..

usage

::

    wflow_pcrglobwb [-h][-v level][-F runinfofile][-L logfile][-C casename][-R runId]
          [-c configfile][-T last_step][-S first_step][-s seconds][-W][-E][-N][-U discharge]
          [-P parameter multiplication][-X][-f][-I][-i tbl_dir][-x subcatchId][-u updatecols]
          [-p inputparameter multiplication][-l loglevel]


    -X: save state at the end of the run over the initial conditions at the start

    -f: Force overwrite of existing results

    -T: Set end time of the run: yyyy-mm-dd hh:mm:ss

    -S: Set start time of the run: yyyy-mm-dd hh:mm:ss

    -s: Set the model timesteps in seconds

    -I: re-initialize the initial model conditions with default

    -i: Set input table directory (default is intbl)

    -x: Apply multipliers (-P/-p ) for subcatchment only (e.g. -x 1)

    -C: set the name  of the case (directory) to run

    -R: set the name runId within the current case

    -L: set the logfile

    -E: Switch on reinfiltration of overland flow

    -c: name of wflow the configuration file (default: Casename/wflow_sbm.ini).

    -h: print usage information

    -W: If set, this flag indicates that an ldd is created for the water level
        for each timestep. If not the water is assumed to flow according to the
        DEM. Wflow will run a lot slower with this option. Most of the time
        (shallow soil, steep topography) you do not need this option. Also, if you
        need it you migth actually need another model.

    -U: The argument to this option should be a .tss file with measured discharge in
        [m^3/s] which the progam will use to update the internal state to match
        the measured flow. The number of columns in this file should match the
        number of gauges in the wflow\_gauges.map file.

    -u: list of gauges/columns to use in update. Format:
        -u [1 , 4 ,13]
        The above example uses column 1, 4 and 13

    -P: set parameter change string (e.g: -P "self.FC = self.FC * 1.6") for non-dynamic variables

    -p: set parameter change string (e.g: -P "self.Precipitation = self.Precipitation * 1.11") for
        dynamic variables


    -l: loglevel (most be one of DEBUG, WARNING, ERROR)



"""


# import pcrut
import sys
import os
import os.path
import getopt


from wflow.wf_DynamicFramework import *
from wflow.wflow_funcs import *
from wflow.wflow_adapt import *

from wflow.pcrglobwb import landSurface
from wflow.pcrglobwb import groundwater
from wflow.pcrglobwb import routing

wflow = "wflow_pcrglobwb: "


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


class Struct(object):
    def __init__(self, **entries):
        self.__dict__.update(entries)


def getLandSurfaceStates(landSurface):
    if landSurface.numberOfSoilLayers == 2:
        for coverType in landSurface.coverTypes:
            setattr(
                landSurface,
                "interceptStor_" + str(coverType),
                landSurface.landCoverObj[coverType].interceptStor,
            )
            setattr(
                landSurface,
                "snowCoverSWE_" + str(coverType),
                landSurface.landCoverObj[coverType].snowCoverSWE,
            )
            setattr(
                landSurface,
                "snowFreeWater_" + str(coverType),
                landSurface.landCoverObj[coverType].snowFreeWater,
            )
            setattr(
                landSurface,
                "topWaterLayer_" + str(coverType),
                landSurface.landCoverObj[coverType].topWaterLayer,
            )
            setattr(
                landSurface,
                "storUpp_" + str(coverType),
                landSurface.landCoverObj[coverType].storUpp,
            )
            setattr(
                landSurface,
                "storLow_" + str(coverType),
                landSurface.landCoverObj[coverType].storLow,
            )
            setattr(
                landSurface,
                "interflow_" + str(coverType),
                landSurface.landCoverObj[coverType].interflow,
            )

    if landSurface.numberOfSoilLayers == 3:
        for coverType in landSurface.coverTypes:
            setattr(
                landSurface,
                "interceptStor_" + str(coverType),
                landSurface.landCoverObj[coverType].interceptStor,
            )
            setattr(
                landSurface,
                "snowCoverSWE_" + str(coverType),
                landSurface.landCoverObj[coverType].snowCoverSWE,
            )
            setattr(
                landSurface,
                "snowFreeWater_" + str(coverType),
                landSurface.landCoverObj[coverType].snowFreeWater,
            )
            setattr(
                landSurface,
                "topWaterLayer_" + str(coverType),
                landSurface.landCoverObj[coverType].topWaterLayer,
            )
            setattr(
                landSurface,
                "storUpp000005_" + str(coverType),
                landSurface.landCoverObj[coverType].storUpp000005,
            )
            setattr(
                landSurface,
                "storUpp005030_" + str(coverType),
                landSurface.landCoverObj[coverType].storUpp005030,
            )
            setattr(
                landSurface,
                "storLow030150_" + str(coverType),
                landSurface.landCoverObj[coverType].storLow030150,
            )
            setattr(
                landSurface,
                "interflow_" + str(coverType),
                landSurface.landCoverObj[coverType].interflow,
            )


def setLandSurfaceStates(landSurface):
    if landSurface.numberOfSoilLayers == 2:
        for coverType in landSurface.coverTypes:
            landSurface.landCoverObj[coverType].interceptStor = (
                landSurface.interceptStor + "_" + str(coverType)
            )
            landSurface.landCoverObj[coverType].snowCoverSWE = (
                landSurface.snowCoverSWE + "_" + str(coverType)
            )
            landSurface.landCoverObj[coverType].snowFreeWater = (
                landSurface.snowFreeWater + "_" + str(coverType)
            )
            landSurface.landCoverObj[coverType].topWaterLayer = (
                landSurface.topWaterLayer + "_" + str(coverType)
            )
            landSurface.landCoverObj[coverType].storUpp = (
                landSurface.storUpp + "_" + str(coverType)
            )
            landSurface.landCoverObj[coverType].storLow = (
                landSurface.storLow + "_" + str(coverType)
            )
            landSurface.landCoverObj[coverType].interflow = (
                landSurface.interflow + "_" + str(coverType)
            )

    if landSurface.numberOfSoilLayers == 3:
        for coverType in landSurface.coverTypes:
            landSurface.landCoverObj[coverType].interceptStor = getattr(
                landSurface, "interceptStor_" + str(coverType)
            )
            landSurface.landCoverObj[coverType].snowCoverSWE = getattr(
                landSurface, "snowCoverSWE_" + str(coverType)
            )
            landSurface.landCoverObj[coverType].snowFreeWater = getattr(
                landSurface, "snowFreeWater_" + str(coverType)
            )
            landSurface.landCoverObj[coverType].topWaterLayer = getattr(
                landSurface, "topWaterLayer_" + str(coverType)
            )
            landSurface.landCoverObj[coverType].storUpp000005 = getattr(
                landSurface, "storUpp000005_" + str(coverType)
            )
            landSurface.landCoverObj[coverType].storUpp005030 = getattr(
                landSurface, "storUpp005030_" + str(coverType)
            )
            landSurface.landCoverObj[coverType].storLow030150 = getattr(
                landSurface, "storLow030150_" + str(coverType)
            )
            landSurface.landCoverObj[coverType].interflow = getattr(
                landSurface, "interflow_" + str(coverType)
            )


class WflowModel(DynamicModel):

    """
  The user defined model class.

  """

    def __init__(self, cloneMap, Dir, RunDir, configfile, staticmaps):
        DynamicModel.__init__(self)

        self.caseName = os.path.abspath(Dir)
        self.runId = RunDir
        self.Dir = os.path.abspath(Dir)
        self.staticmaps = os.path.join(self.Dir, staticmaps)
        self.clonemappath = os.path.join(os.path.abspath(Dir), staticmaps, cloneMap)
        setclone(self.clonemappath)
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

        states = [
            "landSurface.interceptStor_forest",
            "landSurface.interceptStor_grassland",
            "landSurface.snowCoverSWE_forest",
            "landSurface.snowCoverSWE_grassland",
            "landSurface.snowFreeWater_forest",
            "landSurface.snowFreeWater_grassland",
            "landSurface.topWaterLayer_forest",
            "landSurface.topWaterLayer_grassland",
            "landSurface.interflow_forest",
            "landSurface.interflow_grassland",
            "groundwater.storGroundwater",
            "groundwater.storGroundwaterFossil",
            "groundwater.avgAbstraction",
            "groundwater.avgAllocation",
            "groundwater.avgAllocationShort",
            "groundwater.avgNonFossilAllocation",
            "groundwater.avgNonFossilAllocationShort",
            "groundwater.relativeGroundwaterHead",
            "groundwater.baseflow",
            "routing.timestepsToAvgDischarge",
            "routing.channelStorage",
            "routing.readAvlChannelStorage",
            "routing.avgDischarge",
            "routing.m2tDischarge",
            "routing.avgBaseflow",
            "routing.riverbedExchange",
            "routing.avgDischargeShort",
            "routing.subDischarge",
            "routing.waterBodyStorage",
            "routing.avgInflow",
            "routing.avgOutflow",
        ]

        if (
            configget(self.config, "landSurfaceOptions", "includeIrrigation", "False")
            == "True"
        ):
            states += [
                "landSurface.interceptStor_irrPaddy",
                "landSurface.interceptStor_irrNonPaddy",
                "landSurface.snowCoverSWE_irrPaddy",
                "landSurface.snowCoverSWE_irrNonPaddy",
                "landSurface.snowFreeWater_irrPaddy",
                "landSurface.snowFreeWater_irrNonPaddy",
                "landSurface.topWaterLayer_irrPaddy",
                "landSurface.topWaterLayer_irrNonPaddy",
                "landSurface.interflow_irrPaddy",
                "landSurface.interflow_irrNonPaddy",
            ]

        if self.landSurface.numberOfSoilLayers == 2:
            states += [
                "landSurface.storUpp_forest",
                "landSurface.storUpp_grassland",
                "landSurface.storLow_forest",
                "landSurface.storLow_grassland",
            ]
            if (
                configget(
                    self.config, "landSurfaceOptions", "includeIrrigation", "False"
                )
                == "True"
            ):
                states += [
                    "landSurface.storUpp_irrPaddy",
                    "landSurface.storUpp_irrNonPaddy",
                    "landSurface.storLow_irrPaddy",
                    "landSurface.storLow_irrNonPaddy",
                ]

        if self.landSurface.numberOfSoilLayers == 3:
            states += [
                "landSurface.storUpp000005_forest",
                "landSurface.storUpp000005_grassland",
                "landSurface.storUpp005030_forest",
                "landSurface.storUpp005030_grassland",
                "landSurface.storLow030150_forest",
                "landSurface.storLow030150_grassland",
            ]
            if (
                configget(
                    self.config, "landSurfaceOptions", "includeIrrigation", "False"
                )
                == "True"
            ):
                states += [
                    "landSurface.storUpp000005_irrPaddy",
                    "landSurface.storUpp000005_irrNonPaddy",
                    "landSurface.storUpp005030_irrPaddy",
                    "landSurface.storUpp005030_irrNonPaddy",
                    "landSurface.storLow030150_irrPaddy",
                    "landSurface.storLow030150_irrNonPaddy",
                ]

        return states

    # The following are made to better connect to deltashell/openmi
    def supplyCurrentTime(self):
        """
      gets the current time in seconds after the start of the run

      Ouput:
          - time in seconds since the start of the model run
      """
        return self.currentTimeStep() * int(
            configget(self.config, "model", "timestepsecs", "86400")
        )

    def supplyTimeInfo(self):

        timeInfo = {}

        timeInfo["timeStepPCR"] = self.currentTimeStep()
        timeInfo["day"] = self.wf_supplyCurrentDateTime() + timedelta(days=1)
        timeInfo["fulldate"] = "%04i-%02i-%02i" % (
            timeInfo["day"].year,
            timeInfo["day"].month,
            timeInfo["day"].day,
        )
        timeInfo["month"] = timeInfo["day"].month
        timeInfo["year"] = timeInfo["day"].year
        timeInfo["yesterday"] = timeInfo["day"] - timedelta(days=1)
        timeInfo["doy"] = timeInfo["day"].timetuple().tm_yday
        timeInfo["isLastDayOfYear"] = (
            timeInfo["day"] + timedelta(days=1)
        ).timetuple().tm_yday == 1
        timeInfo["endMonth"] = (timeInfo["day"] + timedelta(days=1)).day == 1
        timeInfo["monthIdx"] = self.monthIdx
        timeInfo["endYear"] = (
            timeInfo["day"] + timedelta(days=1)
        ).timetuple().tm_yday == 1
        timeInfo["annuaIdx"] = self.annuaIdx

        # to fix as 'real' object?
        self.timeInfo = Struct(**timeInfo)

        return self.timeInfo

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
        # meteooptions in PCRGLOBWB is replaced with the WFlow mapstacks

        # Meteo and other forcing

        self.P_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "precipitation", "/inmaps/P"
        )  # timeseries for rainfall
        self.PET_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "referencePotET", "/inmaps/PET"
        )  # timeseries for rainfall"/inmaps/PET"          # potential evapotranspiration
        self.TEMP_mapstack = self.Dir + configget(
            self.config, "inputmapstacks", "temperature", "/inmaps/TEMP"
        )  # timeseries for rainfall "/inmaps/TEMP"          # global radiation
        # Meteo and other forcing
        modelparameters.append(
            self.ParamType(
                name="precipitation",
                stack=self.P_mapstack,
                type="timeseries",
                default=0.0,
                verbose=True,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="referencePotET",
                stack=self.PET_mapstack,
                type="timeseries",
                default=0.0,
                verbose=True,
                lookupmaps=[],
            )
        )
        modelparameters.append(
            self.ParamType(
                name="temperature",
                stack=self.TEMP_mapstack,
                type="timeseries",
                default=10.0,
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
            self.wf_suspend(self.SaveDir + "/instate/")

    def initial(self):

        # from wflow.pcrglobwb import landSurface
        # from wflow.pcrglobwb import groundwater
        # from wflow.pcrglobwb import routing

        self.OverWriteInit = int(configget(self.config, "model", "OverWriteInit", "0"))

        initialState = None

        landmask = configget(
            self.config, "globalOptions", "landmask", "wflow_landmask.map"
        )
        lddMap = configget(self.config, "routingOptions", "lddMap", "wflow_ldd.map")

        wflow_landmask = self.wf_readmap(
            os.path.join(self.staticmaps, landmask), 0.0, fail=True
        )
        wflow_ldd = ldd(
            self.wf_readmap(os.path.join(self.staticmaps, lddMap), 0.0, fail=True)
        )

        self.monthIdx = 0
        self.annuaIdx = 0

        startTime = self.wf_supplyStartDateTime()

        self.landSurface = landSurface.LandSurface(
            self.config,
            wflow_landmask,
            self.Dir,
            self.staticmaps,
            self.clonemappath,
            startTime,
            initialState,
        )
        self.groundwater = groundwater.Groundwater(
            self.config,
            wflow_landmask,
            initialState,
            self.Dir,
            self.staticmaps,
            self.clonemappath,
        )
        self.routing = routing.Routing(
            self.config,
            initialState,
            wflow_ldd,
            self.Dir,
            self.staticmaps,
            self.clonemappath,
        )

        self.wf_updateparameters()

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
            self.logger.info("Setting initial conditions to default (zero!)")

        else:
            self.wf_resume(os.path.join(self.Dir, "instate"))
            setLandSurfaceStates(self.landSurface)

    def dynamic(self):

        self.wf_updateparameters()

        self.currTimeStep = self.supplyTimeInfo()

        if self.currTimeStep.isLastDayOfYear:
            self.annuaIdx = self.annuaIdx + 1

        if self.currTimeStep.endMonth:
            self.monthIdx = self.monthIdx + 1

        meteo = {}
        meteo["precipitation"] = self.precipitation
        meteo["temperature"] = self.temperature
        meteo["referencePotET"] = self.referencePotET

        # to FIX as 'real' object?
        self.meteo = Struct(**meteo)

        self.landSurface.update(
            self.meteo, self.groundwater, self.routing, self.currTimeStep, self.logger
        )

        self.groundwater.update(self.landSurface, self.routing, self.currTimeStep)

        self.routing.update(
            self.landSurface, self.groundwater, self.currTimeStep, self.meteo
        )

        getLandSurfaceStates(self.landSurface)


def main(argv=None):

    """
    Perform command line execution of the model.
    """
    global multpars
    global updateCols
    caseName = "default_pcrglobwb"
    runId = "run_default"
    configfile = "wflow_pcrglobwb.ini"
    staticmaps = "staticmaps"
    LogFileName = "wflow.log"
    _lastTimeStep = 0
    _firstTimeStep = 0

    runinfoFile = "runinfo.xml"
    timestepsecs = 86400
    wflow_cloneMap = "wflow_clone.map"
    _NoOverWrite = 1
    loglevel = logging.DEBUG

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return

    ########################################################################
    ## Process command-line options                                        #
    ########################################################################
    try:
        opts, args = getopt.getopt(argv, "XL:hC:Ii:v:S:T:WR:u:s:EP:p:Xx:U:fOc:l:d:")
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
        if o == "-h":
            usage()
        if o == "-f":
            _NoOverWrite = 0
        if o == "-l":
            exec("loglevel = logging." + a)
        if o == "-d":
            staticmaps = a

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

    myModel = WflowModel(wflow_cloneMap, caseName, runId, configfile, staticmaps)
    dynModelFw = wf_DynamicFramework(
        myModel, _lastTimeStep, firstTimestep=_firstTimeStep, datetimestart=starttime
    )
    dynModelFw.createRunId(
        NoOverWrite=_NoOverWrite,
        level=loglevel,
        logfname=LogFileName,
        model="wflow_pcrglobwb",
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
            zz = []
            exec("zz =" + a)
            updateCols = zz
        if o == "-E":
            configset(myModel.config, "model", "reInfilt", "1", overwrite=True)
        if o == "-R":
            runId = a
        if o == "-W":
            configset(myModel.config, "model", "waterdem", "1", overwrite=True)
        if o == "-T":
            configset(myModel.config, "run", "endtime", a, overwrite=True)
        if o == "-S":
            configset(myModel.config, "run", "starttime", a, overwrite=True)

    dynModelFw.setupFramework()
    dynModelFw._runInitial()
    # dynModelFw._runResume()
    # dynModelFw._runDynamic(0, 0)
    dynModelFw._runDynamic(_firstTimeStep, _lastTimeStep)
    dynModelFw._runSuspend()
    dynModelFw._wf_shutdown()


if __name__ == "__main__":
    main()
