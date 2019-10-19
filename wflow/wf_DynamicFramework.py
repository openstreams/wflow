"""
wf_DynamicFramework
-------------------

This is a replacement for the standard pcraster/python DynamicFramwork class.\
It provides extra functionality to simplify linking the models build in the framework
with other models. The provided functionality allows external programs to control
and interrogate the model.

"""

# TODO: Remove command-line options from models such as -F that is now in the ini
# TODO: Fix timestep not forewarding in BMI runs (for reading writing maps)

import calendar
import configparser
import csv
import datetime
import glob
import logging
import shutil
import traceback
from collections import namedtuple
from functools import reduce

import numpy as np
from wflow.wf_netcdfio import *
import pcraster as pcr
import pcraster.framework

from wflow import pcrut
from wflow import wflow_adapt
from wflow.wflow_lib import *
from wflow import __version__
import time  # last to prevent clobbering by *


def log_uncaught_exceptions(ex_cls, ex, tb):
    global logging
    logging.error("".join(traceback.format_tb(tb)))
    logging.error("{0}: {1}".format(ex_cls, ex))


sys.excepthook = log_uncaught_exceptions

logging.getLogger("foo").addHandler(logging.NullHandler())


class runDateTimeInfo:
    """
    class to maintain and retrieve date/time info of the model run.  IN order to support
    difefrent views on date/time the class supports both a step (each input time is timestep) and
    an interval base method (each model timestep is the interval between two input timesteps)

    """

    def __init__(
        self,
        datetimestart=dt.datetime(1990, 1, 1),
        datetimeend=dt.datetime(1990, 1, 5),
        timestepsecs=86400,
        mode="steps",
    ):
        self.runStartTime = datetimestart
        self.runEndTime = datetimeend
        self.timeStepSecs = timestepsecs
        self.currentTimeStep = 0
        self.lastTimeStep = 0
        self.startadjusted = 0
        self.startendadjusted = 0
        self.currentmode = mode
        self.callstopupdate = 0

        if mode == "steps":
            self.runStateTime = self.runStartTime - datetime.timedelta(
                seconds=self.timeStepSecs
            )
        else:
            self.runStateTime = self.runStartTime

        self.setByBMI = False
        self.currentDateTime = self.runStateTime
        self.outPutStartTime = self.runStateTime + datetime.timedelta(
            seconds=self.timeStepSecs
        )
        self.runTimeSteps = (
            calendar.timegm(self.runEndTime.utctimetuple())
            - calendar.timegm(self.runStateTime.utctimetuple())
        ) / self.timeStepSecs
        self.currentMonth = self.currentDateTime.month
        self.currentYday = self.currentDateTime.timetuple().tm_yday
        self.currentHour = self.currentDateTime.hour
        self.nextDateTime = self.currentDateTime + datetime.timedelta(
            seconds=self.timeStepSecs
        )
        self.lastTimeStep = self.runTimeSteps + self.currentTimeStep

    def __str__(self):
        a = self.__dict__

        return str(a)

    def update(
        self,
        timestepsecs=None,
        datetimestart=None,
        datetimeend=None,
        currentTimeStep=None,
        currentDatetime=None,
        runTimeSteps=None,
        mode="steps",
        incrementStep=False,
        setByBMI=False,
    ):
        """
        Updates the content of the framework date/time object. Use only one input parameter per call. or runTimeSteps and datatimestart at the same time
        use the mode option to switch between steps and intervals ('steps' or 'intervals')

        :param timestepsecs:
        :param datetimestart: data time start of the input data
        :param datetimeend:
        :param currentTimeStep:
        :param currentDatetime:
        :return:
        """
        self.currentmode = mode
        self.callstopupdate = self.callstopupdate + 1

        if setByBMI:
            self.setByBMI = True
        if timestepsecs and not runTimeSteps:
            self.timeStepSecs = timestepsecs
            self.runTimeSteps = (
                calendar.timegm(self.runEndTime.utctimetuple())
                - calendar.timegm(self.runStateTime.utctimetuple())
            ) / self.timeStepSecs

            if self.currentmode == "steps":
                self.runStateTime = self.runStartTime - datetime.timedelta(
                    seconds=self.timeStepSecs
                )

            self.outPutStartTime = self.runStateTime + datetime.timedelta(
                seconds=self.timeStepSecs
            )
        elif timestepsecs and runTimeSteps:
            self.timeStepSecs = timestepsecs
            self.runTimeSteps = runTimeSteps

        if datetimestart:
            self.currentTimeStep = 1

            # if self.startadjusted
            if self.currentmode == "steps":
                self.runStartTime = datetimestart
                self.startadjusted = 0
                self.runStateTime = self.runStartTime - datetime.timedelta(
                    seconds=self.timeStepSecs
                )
            else:
                # self.runStartTime = datetimestart + datetime.timedelta(seconds=self.timeStepSecs)
                self.runStartTime = (
                    datetimestart
                )  # + datetime.timedelta(seconds=self.timeStepSecs)
                self.startadjusted = 1
                self.runStateTime = (
                    self.runStartTime
                )  # - datetime.timedelta(seconds=self.timeStepSecs)

            self.currentDateTime = self.runStateTime
            self.outPutStartTime = self.currentDateTime + datetime.timedelta(
                seconds=self.timeStepSecs
            )
            self.runTimeSteps = (
                calendar.timegm(self.runEndTime.utctimetuple())
                - calendar.timegm(self.runStateTime.utctimetuple())
            ) / self.timeStepSecs

            if self.runTimeSteps < 1:  # End time before start time
                self.runTimeSteps = 1
                self.runEndTime = self.runStateTime + datetime.timedelta(
                    seconds=self.timeStepSecs * self.runTimeSteps
                )

        if datetimestart and runTimeSteps:

            self.currentTimeStep = 1
            self.currentDateTime = self.runStartTime
            if self.currentmode == "steps":
                self.runStartTime = datetimestart
                self.startadjusted = 0
                self.runStateTime = self.runStartTime - datetime.timedelta(
                    seconds=self.timeStepSecs
                )
            else:
                self.runStartTime = (
                    datetimestart
                )  # + datetime.timedelta(seconds=self.timeStepSecs)
                self.startadjusted = 1
                self.runStateTime = self.runStartTime

            self.outPutStartTime = self.runStateTime + datetime.timedelta(
                seconds=self.timeStepSecs
            )
            self.currentDateTime = self.runStartTime
            self.runEndTime = self.runStateTime + datetime.timedelta(
                seconds=self.timeStepSecs * runTimeSteps
            )

        if datetimeend:
            self.runEndTime = datetimeend
            self.runTimeSteps = (
                calendar.timegm(self.runEndTime.utctimetuple())
                - calendar.timegm(self.runStateTime.utctimetuple())
            ) / self.timeStepSecs
            if self.runTimeSteps < 1:  # End time before start time
                self.runTimeSteps = 1
                self.runStartTime = self.runEndTime - datetime.timedelta(
                    seconds=self.timeStepSecs * self.runTimeSteps
                )

        if currentTimeStep and currentTimeStep != self.currentTimeStep:
            self.currentTimeStep = currentTimeStep
            self.currentDateTime = self.runStateTime + datetime.timedelta(
                seconds=self.timeStepSecs * (self.currentTimeStep - 1)
            )

        if incrementStep:
            self.currentTimeStep = self.currentTimeStep + 1
            self.currentDateTime = self.currentDateTime + datetime.timedelta(
                seconds=self.timeStepSecs
            )

        if currentDatetime:
            self.currentDateTime = currentDatetime
            self.currentTimeStep = (
                calendar.timegm(self.currentDateTime.utctimetuple())
                - calendar.timegm(self.runStateTime.utctimetuple())
            ) / self.timeStepSecs + 1

        self.nextDateTime = self.currentDateTime + datetime.timedelta(
            seconds=self.timeStepSecs
        )
        self.lastTimeStep = self.runTimeSteps
        self.currentMonth = self.currentDateTime.month
        self.currentYday = self.currentDateTime.timetuple().tm_yday
        self.currentHour = self.currentDateTime.hour


class wf_exchnageVariables:
    """
    List of exchange variables
    The style determined how they are used
    - 1: read from file like normal
    - 2: set by the api in mem (for consistency this is style 0 in the ini file)
    """

    def __init__(self):
        self.vars = []

    def varexists(self, name):

        exists = 0
        for item in self.vars:
            if item[0] == name:
                exists = 1

        return exists

    def addvar(self, name, role, unit):

        if not self.varexists(name):
            if unit == "0":
                unit = "mm/timestep"
            elif unit == "1":
                unit = "m^3/sec"
            elif unit == "2":
                unit = "ma"
            elif unit == "3":
                unit = "degree Celcius"
            elif unit == "4":
                unit = "mm"
            elif unit == "5":
                unit = "-"

            tvar = [name, role, unit]
            self.vars.append(tvar)

    def getvars(self):

        return self.vars

    def getvarStyle(self, name):
        """
        returns 2 if this is a input variable to be set from api otherwise 1
        ( in the ini 0 is for in memory variables)
        A bit confusing!!!
        """
        for xx in self.vars:
            if xx.__contains__(name):
                if xx[1] == 0:
                    return 2
                else:
                    return 1
        return 1


class wf_online_stats:
    def __init__(self):
        """

        :param invarname:
        :param mode:
        :param points:
        :param filename:
        """
        self.count = {}
        self.rangecount = {}
        self.result = {}
        self.mode = {}
        self.points = {}
        self.filename = {}
        self.statvarname = {}

    def addstat(self, name, mode="mean", points=30, filename=None):
        """

        :param name:
        :param mode:
        :param points:
        :param filename:
        :return:
        """
        self.statvarname[name] = name + "_" + mode + "_" + str(points)
        self.mode[name] = mode
        self.points[name] = points
        self.count[name] = 0
        self.rangecount[name] = 0
        self.filename[name] = filename

    def getstat(self, data, name):
        """

        :param data:
        :param name:
        :return:
        """
        if self.count[name] == 0:
            self.result[name] = data
        else:
            if self.mode[name] == "mean":
                self.result[name] = (
                    self.result[name] * (self.points[name] - 1) / self.points[name]
                    + data / self.points[name]
                )

        self.count[name] = self.count[name] + 1

        return pcr.scalar(self.result[name])


class wf_sumavg:
    def __init__(self, varname, mode="sum", filename=None):
        """
        Class to hold variable in the usermodel that must be averaged summed etc.
        """
        if filename == None:
            filename = varname
        self.mode = mode
        self.varname = varname
        self.filename = filename
        self.data = []
        self.count = 0
        self.result = []
        self.availtypes = ["sum", "avg", "min", "max"]

    def add_one(self, data):
        """
        Ad a map (timmestep)
        """
        if self.count == 0:
            self.data = data
        else:
            if self.mode == "sum" or self.mode == "avg":
                self.data = self.data + data
            if self.mode == "max":
                self.data = pcr.max(self.data, data)
            if self.mode == "min":
                self.data = pcr.min(self.data, data)
        self.count = self.count + 1

    def finalise(self):
        """
        Perform final calculations if needed (average, etc) and assign to the
        result variable
        """
        if hasattr(self.data, "isSpatial"):
            if self.mode == "sum" or self.mode == "min" or self.mode == "max":
                self.result = self.data
            if self.mode == "avg":
                self.result = self.data / self.count


class wf_OutputTimeSeriesArea:
    def __init__(self, area, oformat="csv", areafunction="average", tformat="steps"):
        """
        Replacement timeseries output function for the pcraster framework

        area - an area-map to average from
        oformat  - format of the output file (csv, txt, tss, only csv and tss at the moment)
        tformat - steps of datetime (format of the timsteps/stamp)

        Step 1: make average of variable using the areaverage function
        Step 2: Sample the values from the areas (remember the index so we can do it faster lateron)
        step 3: store them in order
        """

        self.steps = 0
        self.timeformat = tformat
        self.area = area
        self.areanp = pcr.pcr2numpy(area, 0).copy()
        self.oformat = oformat
        self.areafunction = areafunction
        """ average, total, minimum, maximum, majority"""

        self.flatarea, self.idx = np.unique(self.areanp, return_index=True)
        # print self.flatarea
        # self.flatarea = self.flatarea[np.isfinite(self.flatarea)]
        # self.idx = self.idx[np.isfinite(self.flatarea)]
        self.fnamelist = []
        self.writer = []
        self.ofile = []

    def closeall(self):
        """
        Close all open filepointers
        """

        for fp in self.ofile:
            fp.close()

        self.fnamelist = []
        self.writer = []
        self.ofile = []

    def writestep(self, variable, fname, timestep=None, dtobj=None):
        """
        write a single timestep

        variable - pcraster map to save to tss
        fname - name of the timeseries file
        """
        # Add new file if not already present
        if fname not in self.fnamelist:
            bufsize = 1  # Implies line buffered
            self.fnamelist.append(fname)

            self.ofile.append(open(fname, "w", bufsize, newline="\n"))
            if self.oformat == "csv":  # Always the case
                self.writer.append(csv.writer(self.ofile[-1]))
                self.ofile[-1].write("# Timestep,")
                self.writer[-1].writerow(self.flatarea)
            elif self.oformat == "tss":  # test
                self.writer.append(csv.writer(self.ofile[-1], delimiter=" "))
                self.ofile[-1].write("timeseries scalar\n")
                self.ofile[-1].write(str(len(self.flatarea) + 1) + "\n")
                self.ofile[-1].write("timestep\n")
                for idd in self.flatarea:
                    self.ofile[-1].write(str(idd) + "\n")
            else:
                print("Not implemented yet")

        self.steps = self.steps + 1

        tmpvar = pcr.spatial(pcr.scalar(variable))
        if self.areafunction == "average":
            self.resmap = pcr.areaaverage(tmpvar, pcr.nominal(self.area))
        elif self.areafunction == "total":
            self.resmap = pcr.areatotal(tmpvar, pcr.nominal(self.area))
        elif self.areafunction == "maximum":
            self.resmap = pcr.areamaximum(tmpvar, pcr.nominal(self.area))
        elif self.areafunction == "minimum":
            self.resmap = pcr.areaminimum(tmpvar, pcr.nominal(self.area))
        elif self.areafunction == "majority":
            self.resmap = pcr.areamajority(tmpvar, pcr.nominal(self.area))
        else:
            self.resmap = pcr.areaaverage(tmpvar, pcr.nominal(self.area))

        self.remap_np = pcr.pcr2numpy(self.resmap, 0)
        self.flatres = self.remap_np.flatten()[self.idx]

        thiswriter = self.fnamelist.index(fname)
        if dtobj and self.timeformat == "datetime":
            self.writer[thiswriter].writerow([str(dtobj)] + self.flatres.tolist())
        elif timestep:
            self.writer[thiswriter].writerow([timestep] + self.flatres.tolist())
        else:
            self.writer[thiswriter].writerow([self.steps] + self.flatres.tolist())
            # self.flatres = np.insert(self.flatres,0,self.steps)


class wf_DynamicFramework(pcraster.framework.frameworkBase.FrameworkBase):
    ## \brief Constructor
    #
    # \param userModel class containing the user model
    # \param lastTimeStep last timestep to run
    # \param firstTimestep sets the starting timestep of the model (optional,
    #        default is 1)
    #
    def __init__(
        self,
        userModel,
        lastTimeStep=0,
        firstTimestep=1,
        datetimestart=dt.datetime(1990, 1, 1),
        timestepsecs=86400,
    ):
        pcraster.framework.frameworkBase.FrameworkBase.__init__(self)

        self.ParamType = namedtuple(
            "ParamType", "name stack type default verbose lookupmaps"
        )
        self.modelparameters = []  # list of model parameters
        self.modelparameters_changes_once = {}
        self.modelparameters_changes_timestep = {}
        self.exchnageitems = wf_exchnageVariables()
        self.setQuiet(True)
        self.reinit = 0
        self._d_model = userModel
        self._testRequirements()

        dte = datetimestart + datetime.timedelta(
            seconds=(lastTimeStep - firstTimestep) * timestepsecs
        )

        self.DT = runDateTimeInfo(
            timestepsecs=timestepsecs,
            datetimestart=datetimestart,
            datetimeend=dte,
            mode="steps",
        )

        if lastTimeStep != 0:
            if firstTimestep == 0:
                firstTimestep = 1
            self.DT.update(runTimeSteps=(lastTimeStep - firstTimestep))
            self.DT.update(currentTimeStep=firstTimestep - 1)

        self.setviaAPI = {}
        # Flag for each variable. If 1 it is set by the API before this timestep. Reset is done at the end of each timestep

        if firstTimestep > lastTimeStep:
            msg = (
                "Cannot run dynamic framework: Start timestep smaller than end timestep"
            )
            raise pcraster.framework.frameworkBase.FrameworkError(msg)

        # fttb
        self._addMethodToClass(self._readmapNew)
        self._addMethodToClass(self._reportNew)
        self._addMethodToClass(self.wf_suspend)
        self._addMethodToClass(self.wf_resume)
        self._addMethodToClass(self.wf_readmap)
        self._addMethodToClass(self.wf_multparameters)
        self._addMethodToClass(self.wf_readmapClimatology)
        self._addMethodToClass(self.readtblDefault)
        self._addMethodToClass(self.readtblLayersDefault)
        self._addMethodToClass(self.readtblFlexDefault)
        self._addMethodToClass(self.wf_supplyVariableNamesAndRoles)
        self._addMethodToClass(self.wf_updateparameters)
        self._addMethodToClass(self.wf_savesummarymaps)
        self._addMethodToClass(self.wf_supplyStartTimeDOY)
        self._addMethodToClass(self.wf_supplyStartTime)
        self._addMethodToClass(self.wf_supplyJulianDOY)
        self._addMethodToClass(self.wf_supplyStartDateTime)
        self._addMethodToClass(self.wf_supplyCurrentDateTime)
        self._addMethodToClass(self.wf_supplyEndTime)
        self._addAttributeToClass("ParamType", self.ParamType)
        self._addAttributeToClass("timestepsecs", self.DT.timeStepSecs)
        self._addAttributeToClass("__version__", __version__)
        # self._addAttributeToClass("__release__", __release__)
        # self._addAttributeToClass("__build__", __build__)
        self.skipfirsttimestep = 0

        if firstTimestep == 0:
            firstTimestep = 1

        # self._userModel()._setNrTimeSteps(lastTimeStep - firstTimestep + 1)
        # self._userModel()._setNrTimeSteps(self.DT.runTimeSteps)
        # self._d_firstTimestep = 1
        # self._userModel()._setFirstTimeStep(1)
        # self._d_lastTimestep = self.DT.runTimeSteps
        # self.APIDebug = 0
        # self._userModel().currentdatetime = self.DT.currentDateTime
        # self._userModel()._setCurrentTimeStep(firstTimestep)

        self._update_time_from_DT()

        self.TheClone = (
            pcr.scalar(pcr.xcoordinate((pcr.spatial(pcr.boolean(1.0))))) * 0.0
        )

    def _update_time_from_DT(self):
        """

        :return:
        """

        self._userModel()._setNrTimeSteps(int(self.DT.runTimeSteps))
        self._d_firstTimestep = 1
        self._userModel()._setFirstTimeStep(1)
        self._d_lastTimestep = self.DT.runTimeSteps
        self.APIDebug = 0
        self._userModel().currentdatetime = self.DT.currentDateTime

        if self.DT.currentTimeStep == 0:
            self._userModel()._setCurrentTimeStep(1)
        else:
            self._userModel()._setCurrentTimeStep(int(self.DT.currentTimeStep))
        self._userModel().timestepsecs = self.DT.timeStepSecs

    def wf_multparameters(self):
        """

        :return:
        """
        if self._userModel()._inDynamic():
            for cmdd in self.modelparameters_changes_timestep:
                var = cmdd.replace("self._userModel().", "").strip()
                if not hasattr(self._userModel(), var):
                    self.logger.error(
                        "Variable change (apply_timestep) could not be applied to "
                        + str(var)
                    )
                else:
                    setattr(
                        self._userModel(),
                        var,
                        getattr(self._userModel(), var)
                        * float(
                            self.modelparameters_changes_timestep[cmdd].split("*")[1]
                        ),  # self.modelparameters_changes_timestep[cmdd],
                    )
                    self.logger.warning(
                        "Variable change (apply_timestep) applied to "
                        + str(var)
                        + " with factor"
                        + self.modelparameters_changes_timestep[cmdd].split("*")[1]
                    )

        if self._userModel()._inInitial():
            #            import pdb; pdb.set_trace()
            for cmdd in self.modelparameters_changes_once:
                var = cmdd.replace("self._userModel().", "").strip()
                
                # for List objects in topoflex
                if '[' in var:
                    listnr = var.split('[')[-1].split(']')[0]
                    varname = var.split('[')[0]
                    mapmult = getattr(self._userModel(), varname)[int(listnr)] * float(self.modelparameters_changes_once[cmdd].split('*')[1])
                    getattr(self._userModel(), varname)[int(listnr)] = mapmult                    

                    self.logger.warning(
                        "Variable change (apply_once) applied to "
                        + str(var) + " with factor" + self.modelparameters_changes_once[cmdd].split('*')[1]
                    )
                
                # for List objects in sbm
                elif var.split('_')[-1].isdigit():
                    mapmult = getattr(self._userModel(), var.split('_')[0])[int(var.split('_')[-1])] * float(self.modelparameters_changes_once[cmdd].split("*")[1])
                    getattr(self._userModel(),var.split('_')[0])[int(var.split('_')[-1])] = mapmult

                    self.logger.warning(
                    "Variable change (apply_once) applied to "
                    + str(var.split('_')[0]) + str([int(var.split('_')[-1])]) + " with factor" + self.modelparameters_changes_once[cmdd].split('*')[1]
                    )

                elif not hasattr(self._userModel(), var):
                    self.logger.error(
                        "Variable change (apply_once) could not be applied to "
                        + str(var)
                    )
                else:
                    setattr(
                        self._userModel(),
                        var,
                        getattr(self._userModel(), var)
                        * float(self.modelparameters_changes_once[cmdd].split("*")[1]),
                    )
                    self.logger.warning(
                        "Variable change (apply_once) applied to "
                        + str(var)
                        + " with factor"
                        + self.modelparameters_changes_once[cmdd].split("*")[1]
                    )

    def wf_updateparameters(self):
        """
        Update the model Parameters (can be used in static and dynamic part of the model)

        It does this by looking at the parameters listed in [parameters] section in the
        ini file and those defined in the parameters() function in the actual model
        (defined by the model developer).

        :return nothing:
        """

        for par in self.modelparameters:

            if self._userModel()._inInitial():
                if par.type == "tbl" or par.type == "tblsparse":
                    if not hasattr(self._userModel(), par.name):
                        self._userModel().logger.info(
                            "Initial: Adding " + par.name + " to model."
                        )
                    tblname = os.path.join(self._userModel().Dir, par.stack)
                    theparmap = self.readtblFlexDefault(
                        tblname, par.default, *par.lookupmaps
                    )
                    setattr(self._userModel(), par.name, theparmap)

                if par.type == "statictbl":
                    if not hasattr(self._userModel(), par.name):
                        self._userModel().logger.info(
                            "Adding " + par.name + " to model."
                        )
                    tblname = os.path.join(self._userModel().Dir, par.stack)
                    theparmap = self.readtblDefault(
                        tblname,
                        self._userModel().LandUse,
                        self._userModel().TopoId,
                        self._userModel().Soil,
                        par.default,
                    )
                    setattr(self._userModel(), par.name, theparmap)

                if par.type == "staticmap":
                    if not hasattr(self._userModel(), par.name):
                        self._userModel().logger.info(
                            "Adding " + par.name + " to model."
                        )
                    fname = os.path.join(self._userModel().Dir, par.stack)
                    fileName, fileExtension = os.path.splitext(fname)
                    if fileExtension == ".map":
                        theparmap = self.wf_readmap(
                            fname, par.default, fail=int(par.verbose)
                        )
                    else:
                        self._userModel().logger.error(
                            fname + " Does not have a .map extension"
                        )

                    setattr(self._userModel(), par.name, theparmap)

            if self._userModel()._inDynamic() or self._userModel()._inInitial():
                if par.type == "timeseries":
                    if not hasattr(self._userModel(), par.name):
                        self._userModel().logger.info(
                            "Adding " + par.name + " to model."
                        )

                    theparmap = self.wf_readmap(
                        os.path.join(self._userModel().caseName, par.stack),
                        par.default,
                        verbose=int(par.verbose),
                    )
                    theparmap = pcr.cover(theparmap, par.default)
                    setattr(self._userModel(), par.name, theparmap)

                if par.type == "monthlyclim":
                    if not hasattr(self._userModel(), par.name):
                        self._userModel().logger.info(
                            "Adding " + par.name + " to model."
                        )
                    theparmap = self.wf_readmapClimatology(
                        os.path.join(self._userModel().caseName, par.stack),
                        kind=1,
                        default=par.default,
                        verbose=int(par.verbose),
                    )
                    theparmap = pcr.cover(theparmap, par.default)
                    setattr(self._userModel(), par.name, theparmap)

                if par.type == "tblmonthlyclim":
                    if not hasattr(self._userModel(), par.name):
                        self._userModel().logger.info(
                            "Initial: Adding " + par.name + " to model."
                        )
                    month = self.DT.currentDateTime.month
                    ptex = os.path.splitext(par.stack)
                    newName = ptex[0] + "_" + str(month) + ptex[1]
                    tblname = os.path.join(self._userModel().Dir, newName)
                    theparmap = self.readtblFlexDefault(
                        tblname, par.default, *par.lookupmaps
                    )
                    setattr(self._userModel(), par.name, theparmap)

                if par.type == "hourlyclim":
                    if not hasattr(self._userModel(), par.name):
                        self._userModel().logger.info(
                            "Adding " + par.name + " to model."
                        )
                    print("hourlyclim has " + par.name + par.stack)
                    print("not been implemented yet")

                if par.type == "dailyclim":
                    if not hasattr(self._userModel(), par.name):
                        self._userModel().logger.info(
                            par.name + " is not defined yet, adding anyway."
                        )
                    theparmap = self.wf_readmapClimatology(
                        os.path.join(self._userModel().caseName, par.stack),
                        kind=2,
                        default=par.default,
                        verbose=int(par.verbose),
                    )
                    setattr(self._userModel(), par.name, theparmap)

            if self._userModel()._inDynamic():
                if par.type == "tss":
                    if not hasattr(self._userModel(), par.name):
                        self._userModel().logger.info(
                            par.name + " is not defined yet, adding anyway."
                        )
                    theparmap = self.wf_timeinputscalar(
                        os.path.join(self._userModel().caseName, par.stack),
                        os.path.join(self._userModel().caseName, par.lookupmaps[0]),
                        par.default,
                    )
                    setattr(self._userModel(), par.name, theparmap)

                if par.type == "tblts":
                    if not hasattr(self._userModel(), par.name):
                        self._userModel().logger.info(
                            "Adding " + par.name + " to model."
                        )
                    tblname = os.path.join(
                        self._userModel().Dir,
                        par.stack + "_" + str(self._userModel().currentStep),
                    )
                    theparmap = self.readtblFlexDefault(
                        tblname, par.default, *par.lookupmaps
                    )
                    setattr(self._userModel(), par.name, theparmap)

                if par.type == "tblsparse":
                    if not hasattr(self._userModel(), par.name):
                        self._userModel().logger.info(
                            "Adding " + par.name + " to model."
                        )

                    tblname = os.path.join(
                        self._userModel().Dir,
                        par.stack + "_" + str(self._userModel().currentStep),
                    )
                    # Only added a new table if available
                    if os.path.exists(tblname):
                        theparmap = self.readtblFlexDefault(
                            tblname, par.default, *par.lookupmaps
                        )
                        setattr(self._userModel(), par.name, theparmap)

        self.setviaAPI = {}

    def wf_supplyStartTimeDOY(self):
        DOY = self.DT.runStartTime.utctimetuple().tm_yday

        return DOY

    def wf_supplyJulianDOY(self):

        JDOY = self.DT.currentYday - (
            calendar.isleap(self.DT.currentDateTime.timetuple().tm_year)
            and self.DT.currentYday > 60
        )

        return JDOY

    def wf_timeinputscalar(self, tssfile, areamap, default):
        """

        :param tssfile:
        :param areamap:
        :return: tss converted to a map
        """
        return pcr.cover(pcr.timeinputscalar(tssfile, pcr.nominal(areamap)), default)

    def _wf_shutdown(self):
        """
        Makes sure the logging closed
        """
        if hasattr(self, "NcOutput"):
            self.NcOutput.finish()

        fp = open(
            os.path.join(
                self._userModel().caseName, self._userModel().runId, "configofrun.ini"
            ),
            "w",
        )
        self._userModel().config.write(fp)

        fp.close()

        for key, value in self.oscv.items():
            value.closeall()

    def loggingSetUp(
        self, caseName, runId, logfname, model, modelversion, level=pcrut.logging.INFO
    ):
        """
        Sets up the logging system assuming we are in the runId directory
        """

        # Set logging
        logfile = os.path.join(caseName, runId, logfname)
        logger = pcrut.setlogger(logfile, model, thelevel=level)
        logger.info(
            model + " " + modelversion + " Case: " + caseName + " Runid: " + runId
        )

        return logger

    def readtblDefault(self, pathtotbl, landuse, subcatch, soil, default):
        """
        First check if a prepared  maps of the same name is present
        in the staticmaps directory. next try to
        read a tbl file to match a landuse, catchment and soil map. Returns
        the default value if the tbl file is not found.

        Finally check of a tbl file exists with a .mult postfix (e.g. Cmax.tbl.mult) and apply the
        multiplication to the loaded data.

        Input:
            -  pathtotbl: full path to table file
            -  landuse: landuse map
            -  subcatch: subcatchment map
            -  soil: soil map
            -  default: default value


        Output:
            - map constructed from tbl file or map with default value

        .. todo::

            Add checking for missing values
        """

        mapname = os.path.join(
            os.path.dirname(pathtotbl),
            "../staticmaps",
            os.path.splitext(os.path.basename(pathtotbl))[0] + ".map",
        )
        if os.path.exists(mapname):
            self.logger.info("reading map parameter file: " + mapname)
            rest = pcr.cover(pcr.readmap(mapname), default)
        else:
            if os.path.isfile(pathtotbl):
                rest = pcr.lookupscalar(pathtotbl, landuse, subcatch, soil)
                self.logger.info("Creating map from table: " + pathtotbl)
            else:
                self.logger.warning(
                    "tbl file not found ("
                    + pathtotbl
                    + ") returning default value: "
                    + str(default)
                )
                rest = pcr.spatial(pcr.cover(pcr.scalar(default)))

            cmask = self._userModel().TopoId

            cmask = pcr.ifthen(cmask > 0, cmask)
            totalzeromap = pcr.pcr2numpy(
                pcr.maptotal(pcr.scalar(pcr.defined(cmask))), 0
            )
            resttotal = pcr.pcr2numpy(pcr.maptotal(pcr.scalar(pcr.defined(rest))), 0)

            if resttotal[0, 0] < totalzeromap[0, 0]:
                self.logger.error(
                    "Not all catchment cells have a value for ["
                    + pathtotbl
                    + "] : "
                    + str(resttotal[0, 0])
                    + "!="
                    + str(totalzeromap[0, 0])
                )

        # Apply multiplication table if present
        multname = os.path.dirname(pathtotbl) + ".mult"
        if os.path.exists(multname):
            multfac = pcr.lookupscalar(multname, landuse, subcatch, soil)
            rest = rest * multfac
            self.logger.info("Applying multiplication from table: " + multname)

        return rest

    def readtblLayersDefault(self, pathtotbl, landuse, subcatch, soil, n, default):
        """
        First check if a prepared  maps of the same name is present
        in the staticmaps directory. next try to
        read a tbl file to match a landuse, catchment and soil map. Returns
        the default value if the tbl file is not found.

        Finally check of a tbl file exists with a .mult postfix (e.g. Cmax.tbl.mult) and apply the
        multiplication to the loaded data.

        Input:
            -  pathtotbl: full path to table file
            -  landuse: landuse map
            -  subcatch: subcatchment map
            -  soil: soil map
            -  default: default value


        Output:
            - map constructed from tbl file or map with default value

        .. todo::

            Add checking for missing values
        """

        mapname = os.path.join(
            os.path.dirname(pathtotbl),
            "../staticmaps",
            os.path.splitext(os.path.basename(pathtotbl))[0] + "_" + str(n) + ".map",
        )

        if os.path.exists(mapname):
            self.logger.info("reading map parameter file: " + mapname)
            rest = pcr.cover(pcr.readmap(mapname), default)
        else:
            if os.path.isfile(pathtotbl):
                rest = pcr.lookupscalar(
                    pathtotbl, landuse, subcatch, soil, pcr.cover(0.0) + n
                )
                self.logger.info("Creating map from table: " + pathtotbl)
            else:
                self.logger.warning(
                    "tbl file not found ("
                    + pathtotbl
                    + ") returning default value: "
                    + str(default)
                )
                rest = pcr.spatial(pcr.cover(pcr.scalar(default)))

            cmask = self._userModel().TopoId

            cmask = pcr.ifthen(cmask > 0, cmask)
            totalzeromap = pcr.pcr2numpy(
                pcr.maptotal(pcr.scalar(pcr.defined(cmask))), 0
            )
            resttotal = pcr.pcr2numpy(pcr.maptotal(pcr.scalar(pcr.defined(rest))), 0)

            if resttotal[0, 0] < totalzeromap[0, 0]:
                self.logger.warning(
                    "Not all catchment cells have a value for ["
                    + pathtotbl
                    + "] : "
                    + str(resttotal[0, 0])
                    + "!="
                    + str(totalzeromap[0, 0])
                )

        # Apply multiplication table if present
        multname = os.path.dirname(pathtotbl) + ".mult"
        if os.path.exists(multname):
            multfac = pcr.lookupscalar(multname, landuse, subcatch, soil)
            rest = rest * multfac
            self.logger.info("Applying multiplication from table: " + multname)

        return rest

    def readtblFlexDefault(self, pathtotbl, default, *args):
        """
        First check if a prepared  maps of the same name is present
        in the staticmaps directory. next try to
        read a tbl file to match  a number of maps. Returns
        the default value if the tbl file is not found.

        Finally check of a tbl file exists with a .mult postfix (e.g. Cmax.tbl.mult) and apply the
        multiplication to the loaded data.

        Input:
            -  pathtotbl: full path to table file
            -  default: default value
            - *args: maps for the lookup table directly fed to lookupscalar

        Output:
            - map constructed from tbl file or map with default value

        .. todo::

            Add checking for missing values
        """

        mapname = (
            os.path.dirname(pathtotbl)
            + "/../staticmaps/"
            + os.path.splitext(os.path.basename(pathtotbl))[0]
            + ".map"
        )
        if os.path.exists(mapname):
            self.logger.info("Reading map parameter file: " + mapname)
            rest = pcr.cover(pcr.readmap(mapname), default)
        else:
            if os.path.isfile(pathtotbl):
                newargs = []
                args = list(args)
                for mapje in args:
                    if len(os.path.splitext(mapje)[1]) > 1:  # We have an extension...
                        newargs.append(os.path.join(self._userModel().caseName, mapje))
                        # we specify a full map
                    else:
                        # Assume we have monthly climatology as no extension is present
                        theparmap = self.wf_readmapClimatology(
                            os.path.join(self._userModel().caseName, mapje),
                            kind=1,
                            default=default,
                            verbose=True,
                        )
                        theparmap = pcr.cover(theparmap, default)
                        newargs.append(theparmap)

                for lmap in newargs:
                    if not os.path.exists(lmap):
                        rest = pcr.spatial(pcr.scalar(default))
                        self.logger.debug(
                            "map file not found ("
                            + lmap
                            + ") returning default value: "
                            + str(default)
                        )
                    else:
                        rest = pcr.lookupscalar(pathtotbl, *newargs)
            else:
                self.logger.debug(
                    "tbl file not found ("
                    + pathtotbl
                    + ") returning default value: "
                    + str(default)
                )
                rest = pcr.spatial(pcr.scalar(default))

                # cmask = self._userModel().TopoId

                # cmask = pcr.ifthen(cmask > 0,cmask)
                # totalzeromap = pcr.pcr2numpy(pcr.maptotal(pcr.scalar(pcr.defined(cmask))),0)
                # resttotal = pcr.pcr2numpy(pcr.maptotal(pcr.scalar(pcr.defined(rest))),0)

                # if resttotal[0,0] < totalzeromap[0,0]:
                #    self.logger.warning("Not all catchment cells have a value for [" + pathtotbl + "] : " + str(resttotal[0,0]) + "!=" + str(totalzeromap[0,0]))

        # Apply multiplication table if present
        multname = os.path.dirname(pathtotbl) + ".mult"
        if os.path.exists(multname):
            multfac = pcr.lookupscalar(multname, *args)
            rest = rest * multfac
            self.logger.info("Applying multiplication from table: " + multname)

        return rest

    def createRunId(
        self,
        intbl="intbl",
        logfname="wflow.log",
        NoOverWrite=True,
        model="model",
        modelVersion="no version",
        level=pcrut.logging.DEBUG,
        doSetupFramework=True,
    ):
        """
        Create runId dir and copy table files to it
        Also changes the working dir to the case/runid directory
        """

        caseName = self._userModel().caseName
        runId = self._userModel().runId
        if modelVersion == "no version":
            modelVersion = __version__

        configfile = self._userModel().configfile
        if not os.path.isdir(caseName + "/" + runId):
            os.makedirs(caseName + "/" + runId + "/outmaps/")
            os.makedirs(caseName + "/" + runId + "/outstate/")
            os.makedirs(caseName + "/" + runId + "/outsum/")
            os.makedirs(caseName + "/" + runId + "/intbl/")
            os.makedirs(caseName + "/" + runId + "/runinfo/")
        else:
            if os.path.exists(caseName + "/" + runId + "/run.tss"):
                if NoOverWrite:
                    print(
                        "ERROR: refusing to overwrite an existing run: "
                        + caseName
                        + "/"
                        + runId
                        + "/run.tss"
                    )
                    sys.exit(1)

        for file in glob.glob(caseName + "/" + intbl + "/*.tbl"):
            shutil.copy(file, caseName + "/" + runId + "/" + intbl)
        try:
            shutil.copy(
                caseName + "/" + configfile, caseName + "/" + runId + "/runinfo"
            )
        except:
            print("Cannot find config file: " + caseName + "/" + configfile)

        self._userModel().logger = self.loggingSetUp(
            caseName, runId, logfname, model, modelVersion, level=level
        )

        self.logger = self._userModel().logger

        self.logger.info(
            "Initialise framework version: " + __version__
        )  # + "(" + __release__ + ")")

        global logging
        logging = self.logger

        self._userModel().config = self.iniFileSetUp(caseName, runId, configfile)
        modelnamefromobject = self._userModel().__module__
        self.modelname = configget(
            self._userModel().config, "model", "modeltype", "not set"
        )

        if self.modelname == "not set":
            self.logger.warning(
                "Ini file does not contain model name, assuming " + modelnamefromobject
            )
            self.modelname = modelnamefromobject

        if modelnamefromobject != self.modelname:
            self.logger.warning(
                "Ini file made for "
                + self.modelname
                + " but found "
                + modelnamefromobject
                + " in code."
            )

        self.runlengthdetermination = configget(
            self._userModel().config, "run", "runlengthdetermination", "steps"
        )
        self.DT.update(
            timestepsecs=int(
                configget(self._userModel().config, "run", "timestepsecs", "86400")
            ),
            mode=self.runlengthdetermination,
            runTimeSteps=self.DT.runTimeSteps,
        )
        self._update_time_from_DT()

        if doSetupFramework:
            self.setupFramework()

    def _initAPIVars(self):
        """
        Sets vars in the API that are forcing variables to the model
        """
        apivars = self.wf_supplyVariableNamesAndRoles()

        for var in apivars:
            if not hasattr(self._userModel(), var[0]):
                setattr(self._userModel(), var[0], self.TheClone)

    def setuptimeInfo(self):
        """

        :return:
        """
        from dateutil import parser

        st = configget(self._userModel().config, "run", "starttime", "None")

        # self.skipfirsttimestep =  int(configget(self._userModel().config, 'run', 'skipfirst', "0"))

        # Assume that we have set this via BMI
        if self.DT.setByBMI:
            self.logger.info(
                "Not reading time from ini file, assuming it is set by BMI or otherwise (calls = "
                + str(self.DT.callstopupdate)
                + ")"
            )
        else:
            if st == "None":  # try from the runinfo file
                rinfo_str = configget(
                    self._userModel().config, "run", "runinfo", "None"
                )
                rinfo = os.path.join(self._userModel().Dir, rinfo_str)
                self.DT.update(
                    timestepsecs=int(
                        configget(
                            self._userModel().config, "run", "timestepsecs", "86400"
                        )
                    ),
                    mode=self.runlengthdetermination,
                    runTimeSteps=self.DT.runTimeSteps,
                )
                self._update_time_from_DT()
                if rinfo_str != "None":
                    self.DT.update(
                        datetimestart=wflow_adapt.getStartTimefromRuninfo(rinfo),
                        mode=self.runlengthdetermination,
                    )
                    self.DT.update(
                        datetimeend=wflow_adapt.getEndTimefromRuninfo(rinfo),
                        mode=self.runlengthdetermination,
                    )
                    self._update_time_from_DT()
                    # add one step to start time if it is the same s the state time
                    # if self.skipfirsttimestep:
                    #    self.logger.debug("Skipping first timestep...")
                    #    self.DT.skiptime()

                    self._userModel().currentdatetime = self.DT.currentDateTime

                    self.DT.update(
                        timestepsecs=int(
                            configget(
                                self._userModel().config, "run", "timestepsecs", "86400"
                            )
                        ),
                        mode=self.runlengthdetermination,
                    )
                    self.DT.update(
                        currentTimeStep=self.DT.currentTimeStep,
                        mode=self.runlengthdetermination,
                    )
                    self._update_time_from_DT()
                else:
                    self.DT.update(
                        datetimestart=parser.parse("1990-01-01 00:00:00 GMT"),
                        mode=self.runlengthdetermination,
                    )
                    self.logger.info(
                        "Not enough information in the [run] section. Need start and end time or a runinfo.xml file.... Reverting to default date/time"
                    )
            else:
                self.DT.update(
                    datetimestart=parser.parse(st), mode=self.runlengthdetermination
                )
                self.DT.update(
                    currentTimeStep=self.DT.currentTimeStep,
                    mode=self.runlengthdetermination,
                )
                # if self.skipfirsttimestep:
                #    self.logger.debug("Skipping first timestep...")
                #    self.DT.skiptime()

                self._userModel().currentdatetime = self.DT.currentDateTime
                ed = configget(self._userModel().config, "run", "endtime", "None")
                if ed != "None":
                    self.DT.update(
                        datetimeend=parser.parse(ed), mode=self.runlengthdetermination
                    )
                    self.DT.update(
                        timestepsecs=int(
                            configget(
                                self._userModel().config, "run", "timestepsecs", "86400"
                            )
                        ),
                        mode=self.runlengthdetermination,
                    )
                    self.DT.update(
                        currentTimeStep=self.DT.currentTimeStep,
                        mode=self.runlengthdetermination,
                    )
                else:
                    self.logger.error(
                        "No end time given with start time: [run] endtime = " + ed
                    )
                    sys.exit(1)

                self._update_time_from_DT()

    def setupFramework(self):
        """
        Second step, after setting the log file and reading the ini file get data from config, setup
         IO etc

        :return:
        """

        self._initAPIVars()
        self.framework_setup = True
        caseName = self._userModel().caseName
        runId = self._userModel().runId
        self.outputFormat = int(
            configget(self._userModel().config, "framework", "outputformat", "1")
        )
        self.APIDebug = int(
            configget(
                self._userModel().config, "framework", "debug", str(self.APIDebug)
            )
        )
        self.ncfile = configget(
            self._userModel().config, "framework", "netcdfinput", "None"
        )
        self.ncinfilestates = configget(
            self._userModel().config, "framework", "netcdfstatesinput", "None"
        )
        self.ncoutfile = configget(
            self._userModel().config, "framework", "netcdfoutput", "None"
        )
        self.ncoutfilestatic = configget(
            self._userModel().config, "framework", "netcdfstaticoutput", "None"
        )
        self.ncoutfilestates = configget(
            self._userModel().config, "framework", "netcdfstatesoutput", "None"
        )
        self.ncfilestatic = configget(
            self._userModel().config, "framework", "netcdfstaticinput", "None"
        )
        self.EPSG = configget(
            self._userModel().config, "framework", "EPSG", "EPSG:4326"
        )
        self.ncfileformat = configget(
            self._userModel().config, "framework", "netcdf_format", "NETCDF4"
        )
        self.ncfilecompression = configget(
            self._userModel().config, "framework", "netcdf_zlib", "True"
        )
        self.ncfiledigits = configget(
            self._userModel().config,
            "framework",
            "netcdf_least_significant_digit",
            "None",
        )

        if self.ncfiledigits == "None":
            self.ncfiledigits = None
        else:
            self.ncfiledigits = int(self.ncfiledigits)

        if self.ncfilecompression == "True":
            self.ncfilecompression = True
        else:
            self.ncfilecompression = False

        # Set the re-init hint for the local model
        self.reinit = int(
            configget(self._userModel().config, "run", "reinit", str(self.reinit))
        )
        self._userModel().reinit = self.reinit
        # Now finally set the start end time. First check if set in ini otherwise check if the ini defines
        # a runinfo file

        self.setuptimeInfo()

        # Setup all the netCDF files that may be used for input/output
        if self.ncfile != "None":
            varlst = []
            if hasattr(self._userModel(), "parameters"):
                for ms in self._userModel().parameters():
                    if ms.type == "timeseries":
                        varlst.append(os.path.basename(ms.stack))

            mstacks = configsection(self._userModel().config, "inputmapstacks")
            for ms in mstacks:
                varlst.append(
                    os.path.basename(
                        configget(
                            self._userModel().config, "inputmapstacks", ms, "None"
                        )
                    )
                )

            self.logger.debug(
                "Found following input variables to get from netcdf file: "
                + str(varlst)
            )
            self.NcInput = netcdfinput(
                os.path.join(caseName, self.ncfile), self.logger, varlst
            )

        # Meta info for netcdf files
        meta = {}
        meta["caseName"] = caseName
        meta["runId"] = runId
        meta["wflow_version"] = __version__
        # meta['wflow_release'] = __release__
        # meta['wflow_build'] = __build__
        meta["wflow_ini"] = self._userModel().configfile
        if hasattr(sys, "frozen"):
            meta["wflow_exe"] = "True"
        else:
            meta["wflow_exe"] = "False"

        try:
            metafrom_config = dict(self._userModel().config.items("netcdfmetadata"))
        except:
            metafrom_config = {}

        meta.update(metafrom_config)

        if self.ncinfilestates != "None":
            smaps = self._userModel().stateVariables()
            maps = [s + ".map" for s in smaps]
            self.logger.debug(
                "Found following input states to get from netcdf file: " + str(maps)
            )
            self.NcInputStates = netcdfinputstates(
                os.path.join(caseName, self.ncinfilestates), self.logger, maps
            )

        if self.ncfilestatic != "None":
            self.NcInputStatic = netcdfinputstatic(
                os.path.join(caseName, self.ncfilestatic), self.logger
            )

        if self.ncoutfile != "None":  # Ncoutput
            buffer = int(
                configget(
                    self._userModel().config, "framework", "netcdfwritebuffer", "50"
                )
            )

            self.NcOutput = netcdfoutput(
                os.path.join(caseName, runId, self.ncoutfile),
                self.logger,
                self.DT.outPutStartTime,
                self.DT.runTimeSteps,
                maxbuf=buffer,
                metadata=meta,
                EPSG=self.EPSG,
                timestepsecs=self.DT.timeStepSecs,
                Format=self.ncfileformat,
                zlib=self.ncfilecompression,
                least_significant_digit=self.ncfiledigits,
            )

        if self.ncoutfilestatic != "None":  # Ncoutput
            self.NcOutputStatic = netcdfoutputstatic(
                os.path.join(caseName, runId, self.ncoutfilestatic),
                self.logger,
                self.DT.runEndTime,
                1,
                timestepsecs=self.DT.timeStepSecs,
                maxbuf=1,
                metadata=meta,
                EPSG=self.EPSG,
                Format=self.ncfileformat,
                zlib=self.ncfilecompression,
                least_significant_digit=self.ncfiledigits,
            )

        if self.ncoutfilestates != "None":  # Ncoutput
            self.NcOutputState = netcdfoutputstatic(
                os.path.join(caseName, runId, self.ncoutfilestates),
                self.logger,
                self.DT.runEndTime,
                1,
                timestepsecs=self.DT.timeStepSecs,
                maxbuf=1,
                metadata=meta,
                EPSG=self.EPSG,
                Format=self.ncfileformat,
                zlib=self.ncfilecompression,
                least_significant_digit=self.ncfiledigits,
            )

        # Add the on-line statistics
        self.onlinestat = wf_online_stats()

        rollingvars = configsection(self._userModel().config, "rollingmean")
        for thisvar in rollingvars:
            try:
                thisvarnoself = thisvar.split("self.")[1]
            except:
                logging.error("Entry in ini invalid: " + thisvar)
                raise ValueError
            pts = int(self._userModel().config.get("rollingmean", thisvar))
            self.onlinestat.addstat(thisvarnoself, points=pts)

        # and set the var names
        for key in self.onlinestat.statvarname:
            setattr(
                self._userModel(), self.onlinestat.statvarname[key], self.TheClone * 0.0
            )

        # Fill the summary (stat) list from the ini file
        self.statslst = []
        _type = wf_sumavg(None)
        for sttype in _type.availtypes:
            _maps = configsection(self._userModel().config, "summary_" + sttype)
            for thismap in _maps:
                thismapname = os.path.join(
                    caseName,
                    runId,
                    "outsum",
                    self._userModel().config.get("summary_" + sttype, thismap),
                )
                try:
                    thismap = thismap.split("self.")[1]
                except:
                    logging.error("Entry in ini invalid: " + thismap)
                    raise ValueError

                self.statslst.append(
                    wf_sumavg(thismap, mode=sttype, filename=thismapname)
                )

        # Get model parameters from model object
        if hasattr(self._userModel(), "parameters"):
            self.modelparameters = self._userModel().parameters()
        else:
            self.modelparameters = []
        # Read extra model parameters from ini file
        modpars = configsection(self._userModel().config, "modelparameters")
        for par in modpars:
            aline = self._userModel().config.get("modelparameters", par)
            vals = aline.split(",")
            if len(vals) >= 4:
                # check if par already present
                present = par in [xxx[0] for xxx in self.modelparameters]
                if present:
                    pos = [xxx[0] for xxx in self.modelparameters].index(par)
                    # Check if the existing definition is static, in that case append, otherwise overwrite
                    if "static" in self.modelparameters[pos].type:
                        self._userModel().logger.debug(
                            "Creating extra parameter specification for par: "
                            + par
                            + " ("
                            + str(vals)
                            + ")"
                        )
                        self.modelparameters.append(
                            self.ParamType(
                                name=par,
                                stack=vals[0],
                                type=vals[1],
                                default=float(vals[2]),
                            ),
                            verbose=vals[3],
                            lookupmaps=vals[4:],
                        )
                    else:
                        self._userModel().logger.debug(
                            "Updating existing parameter specification for par: "
                            + par
                            + " ("
                            + str(vals)
                            + ")"
                        )
                        self.modelparameters[pos] = self.ParamType(
                            name=par,
                            stack=vals[0],
                            type=vals[1],
                            default=float(vals[2]),
                            verbose=vals[3],
                            lookupmaps=vals[4:],
                        )
                else:
                    self._userModel().logger.debug(
                        "Creating parameter specification for par: "
                        + par
                        + " ("
                        + str(vals)
                        + ")"
                    )
                    self.modelparameters.append(
                        self.ParamType(
                            name=par,
                            stack=vals[0],
                            type=vals[1],
                            default=float(vals[2]),
                            verbose=vals[3],
                            lookupmaps=vals[4:],
                        )
                    )
            else:
                logging.error("Parameter line in ini not valid: " + aline)
                raise ValueError

        varchanges = configsection(self._userModel().config, "variable_change_once")
        for chvar in varchanges:
            a = chvar.replace("self", "self._userModel()")
            self.modelparameters_changes_once[a] = (
                self._userModel()
                .config.get("variable_change_once", chvar)
                .replace("self", "self._userModel()")
            )

        varchanges = configsection(self._userModel().config, "variable_change_timestep")
        for chvar in varchanges:
            a = chvar.replace("self", "self._userModel()")
            self.modelparameters_changes_timestep[a] = (
                self._userModel()
                .config.get("variable_change_timestep", chvar)
                .replace("self", "self._userModel()")
            )

        # Now gather all the csv/tss/txt etc timeseries output objects
        # Print .ini defined outputmaps per timestep

        checktss = configsection(self._userModel().config, "outputtss")
        if len(checktss) > 0:
            self.logger.warning(
                "Found a outputtss section. This is NOT used anymore in this version. Please use outputtss_0 .. n"
            )

        self.oscv = {}
        self.samplenamecsv = {}
        self.varnamecsv = {}
        for tsformat in ["csv", "tss"]:
            secnr = 0
            toprint = [None]

            while len(toprint) > 0:
                thissection = "output" + tsformat + "_" + str(secnr)
                toprint = configsection(self._userModel().config, thissection)
                secnr = secnr + 1
                samplemapname = os.path.join(
                    caseName,
                    configget(
                        self._userModel().config, thissection, "samplemap", "None"
                    ),
                )
                areafunction = configget(
                    self._userModel().config, thissection, "function", "average"
                )
                timeformat = configget(
                    self._userModel().config, thissection, "timeformat", "steps"
                )
                if "None" not in samplemapname:
                    try:
                        self.samplemap = self.wf_readmap(samplemapname, 0.0, fail=True)
                        idd = tsformat + ":" + samplemapname + ":" + areafunction
                        self.oscv[idd] = wf_OutputTimeSeriesArea(
                            self.samplemap,
                            oformat=tsformat,
                            areafunction=areafunction,
                            tformat=timeformat,
                        )
                        self.logger.info(
                            "Adding "
                            + tsformat
                            + " output at "
                            + samplemapname
                            + " function: "
                            + areafunction
                        )
                    except:
                        self.logger.warning(
                            "Could not read sample id-map for timeseries: "
                            + samplemapname
                        )
                        self.logger.warning(sys.exc_info())
                    for a in toprint:
                        if (
                            "samplemap" not in a
                            and "function" not in a
                            and "timeformat" not in a
                        ):
                            b = a.replace("self", "self._userModel()")
                            fn = os.path.join(
                                caseName,
                                runId,
                                self._userModel().config.get(thissection, a),
                            )
                            self.samplenamecsv[fn] = idd
                            self.varnamecsv[fn] = b

    def wf_suspend(self, directory):
        """
        Suspend the state variables to disk as .map files
        Also saves the summary maps

        """

        self._incrementIndentLevel()
        self._traceIn("suspend")

        allvars = self._userModel().stateVariables()

        for var in allvars:
            try:
                fname = os.path.join(directory, var).replace("\\", "/") + ".map"
                # savevar = getattr(self._userModel(), var)
                savevar = reduce(getattr, var.split("."), self._userModel())

                try:  # Check if we have a list of maps
                    b = len(savevar)
                    a = 0
                    for z in savevar:
                        fname = (
                            os.path.join(directory, var + "_" + str(a)).replace(
                                "\\", "/"
                            )
                            + ".map"
                        )

                        self.reportState(
                            pcr.cover(z), fname, style=1, gzipit=False, longname=fname
                        )
                        a = a + 1
                except:
                    # thevar = getattr(self._userModel(), var)
                    thevar = reduce(getattr, var.split("."), self._userModel())
                    self.reportState(
                        thevar, fname, style=1, gzipit=False, longname=fname
                    )
            except:
                self.logger.warning("Problem saving state variable: " + var)
                # self.logger.warning(execstr)
                self.logger.warning(sys.exc_info())

        # Save the summary maps
        self.wf_savesummarymaps()
        self._traceOut("suspend")
        self._decrementIndentLevel()

    def wf_saveTimeSeries(self):
        """
        Print .ini defined output csv/tss timeseries per timestep
        """

        for a in self.samplenamecsv:
            found = 1
            try:
                if "+" in self.varnamecsv[a]:
                    a_ = self.varnamecsv[a].split("+")
                    tmpvar = pcr.cover(0.0)
                    for i in np.arange(0, len(a_)):
                        tmpvar = tmpvar + reduce(
                            getattr,
                            a_[i].strip().replace("self._userModel().", "").split("."),
                            self._userModel(),
                        )

                else:
                    # this is added for flextopo -- because list of variables for different classes
                    if '[' in self.varnamecsv[a].replace("self._userModel().", ""):
                        listnr = self.varnamecsv[a].replace("self._userModel().", "").split('[')[-1].split(']')[0]
                        varname = self.varnamecsv[a].replace("self._userModel().", "").split('[')[0]
                        tmpvar = getattr(self._userModel(),varname)[int(listnr)]
                        
                    else:
                        tmpvar = reduce(
                            getattr,
                            self.varnamecsv[a].replace("self._userModel().", "").split("."),
                            self._userModel(),
                        )
                        
            except:
                found = 0
                self.logger.fatal(
                    "Cannot find: " + self.varnamecsv[a] + " variable not in model."
                )
                sys.exit(1)

            self.oscv[self.samplenamecsv[a]].writestep(
                tmpvar,
                a,
                timestep=self.DT.currentTimeStep - 1,
                dtobj=self.DT.currentDateTime,
            )

    def wf_savesummarymaps(self):
        """
        Saves the maps defined in the summary section to disk
        [summary] # Single values or end values
        [summary_sum] # accumulative maps over the model run
        [summary_avg] # average of maps over the model run
        [summary_max] # max of maps over the model run
        [summary_min] # min of maps over the model run
        """

        self._userModel().logger.info("Saving summary maps to disk...")
        toprint = configsection(self._userModel().config, "summary")
        for a in toprint:
            b = a.replace("self.", "")
            try:
                #below if statement for topoflex lists code 
                if '[' in str(b):
                    listnr = str(b).split('[')[-1].split(']')[0]
                    varname = str(b).split('[')[0]
                    pcrmap = getattr(self._userModel(),varname)[int(listnr)]
                else:
                    pcrmap = getattr(self._userModel(), b)

                self.reportStatic(
                    pcrmap,
                    os.path.join(
                        self._userModel().Dir,
                        self._userModel().runId,
                        "outsum",
                        self._userModel().config.get("summary", a),
                    ),
                    style=1,
                )

            except:
                self._userModel().logger.warning(
                    "Could not find or save the configured summary map:" + a
                )

        # Check of the usermodel has a list of summary maps defined and save those
        if hasattr(self._userModel(), "default_summarymaps"):
            for a in self._userModel().default_summarymaps():
                b = a.replace("self.", "")
                if hasattr(self._userModel(), b):
                    pcrmap = getattr(self._userModel(), b)
                    # pcr.report( pcrmap , os.path.join(self._userModel().Dir, self._userModel().runId, "outsum", b + ".map" ))
                    self.reportStatic(
                        pcrmap,
                        os.path.join(
                            self._userModel().Dir,
                            self._userModel().runId,
                            "outsum",
                            b + ".map",
                        ),
                        style=1,
                    )

        # These are the ones in the _sum _average etc sections
        for a in range(0, len(self.statslst)):
            self.statslst[a].finalise()
            if hasattr(self.statslst[a].result, "isSpatial"):
                data = self.statslst[a].result
                fname = self.statslst[a].filename
                if hasattr(data, "isSpatial"):
                    # report (data,fname)
                    self.reportStatic(data, fname, style=1)

    def wf_savedynMaps(self):
        """
        Save the maps defined in the ini file for the dynamic section

        .. todo::

            Save maps to be used in memory at startup and do not call the ini file each time

        """
        # Print .ini defined outputmaps per timestep
        toprint = configsection(self._userModel().config, "outputmaps")

        self.logger.info("saving maps")

        for a in toprint:
            report = False
            # possible to add variables
            if "+" in a:
                a_ = a.split("+")
                thevar = pcr.cover(0.0)
                for i in np.arange(0, len(a_)):
                    # check for nested objects
                    if len(a_[i].replace("self.", "").split(".")) > 1:
                        if hasattr(
                            (
                                self._userModel(),
                                a_[i].replace("self.", "").split(".")[0],
                            )
                            and reduce(
                                getattr,
                                a_[i].replace("self.", "").split("."),
                                self._userModel(),
                            )
                            is not None
                        ):

                            thevar = thevar + reduce(
                                getattr,
                                a_[i].replace("self.", "").split("."),
                                self._userModel(),
                            )
                            report = True

                    elif hasattr(self._userModel(), a_[i].strip().replace("self.", "")):
                        thevar = thevar + getattr(
                            self._userModel(), a_[i].strip().replace("self.", "")
                        )
                        report = True

                    else:
                        report = False
                        break

            else:
                # check for nested objects
                if len(a.replace("self.", "").split(".")) > 1:
                    if (
                        hasattr(self._userModel(), a.replace("self.", "").split(".")[0])
                        and reduce(
                            getattr,
                            a.replace("self.", "").split("."),
                            self._userModel(),
                        )
                        is not None
                    ):
                        thevar = reduce(
                            getattr,
                            a.replace("self.", "").split("."),
                            self._userModel(),
                        )
                        report = True
                
                #add lines below for topoflex
                elif '[' in str(a.replace("self.", "")):
                    listnr = str(a.replace("self.", "")).split('[')[-1].split(']')[0]
                    varname = str(a.replace("self.", "")).split('[')[0]
                    thevar = getattr(self._userModel(),varname)[int(listnr)]
                    report = True

                elif hasattr(self._userModel(), a.replace("self.", "")):
                    thevar = getattr(self._userModel(), a.replace("self.", ""))
                    report = True

            if report == True:
                if type(thevar) is list:
                    a = self._userModel().config.get("outputmaps", a)
                    for i in np.arange(0, len(thevar)):
                        thename = a + "_" + str(i) + "_"
                        self._reportNew(
                            thevar[i],
                            os.path.join(
                                self._userModel().Dir,
                                self._userModel().runId,
                                "outmaps",
                                thename,
                            ),
                            longname=thename,
                        )
                else:
                    self._reportNew(
                        thevar,
                        os.path.join(
                            self._userModel().Dir,
                            self._userModel().runId,
                            "outmaps",
                            self._userModel().config.get("outputmaps", a),
                        ),
                        longname=a,
                    )
            else:
                self.logger.warning("outputmap " + a + " not found in usermodel")

    def wf_resume(self, directory):
        """
        Resumes the state variables from disk as .map files (or arrays of maps files using
        a _? postfix)

        """
        self._incrementIndentLevel()
        self._traceIn("resume")
        allvars = self._userModel().stateVariables()

        for var in allvars:
            # First try to read a stack of state files
            stop = 0
            nr = 0

            while stop == 0:
                name = os.path.join(directory, var + "_" + str(nr) + ".map").replace(
                    "\\", "/"
                )

                if os.path.exists(name):
                    tvar = self.wf_readmap(
                        name,
                        0.0,
                        ncfilesource=self.ncinfilestates,
                        fail=True,
                        silent=True,
                    )
                    if nr == 0:
                        setattr(self._userModel(), var, [])
                    getattr(self._userModel(), var).append(tvar)
                    nr = nr + 1
                else:
                    if nr > 0:
                        self.logger.info(
                            "state variable "
                            + str(var)
                            + " contains "
                            + str(nr)
                            + " state files (stack)"
                        )
                    stop = 1

            if nr == 0:
                try:
                    mpath = os.path.join(directory, var + ".map").replace("\\", "/")
                    tvar = self.wf_readmap(mpath, 0.0, ncfilesource=self.ncinfilestates)

                    # check for nested objects
                    if "." in var:
                        attrs = var.split(".")
                        c = getattr(self._userModel(), attrs[0])
                        setattr(c, attrs[1], tvar)
                    else:
                        setattr(self._userModel(), var, tvar)
                except:
                    self.logger.error(
                        "problem while reading state variable from disk: "
                        + mpath
                        + " Suggest to use the -I option to restart"
                    )
                    sys.exit(1)

        self._traceOut("resume")
        self._decrementIndentLevel()

    def wf_QuickSuspend(self):
        """
        Save the state variable of the current timestep in memory
        it uses the wf_supplyVariableNamesAndRoles() function to find them.
        The variables are inserted into the model object
        This function is normally called as part of the run. Normally there is
        no need to call it directly.
        """
        allvars = self._userModel().stateVariables()

        for var in allvars:
            try:
                setattr(
                    self._userModel(),
                    var + "_laststep",
                    reduce(getattr, var.split("."), self._userModel()),
                )
            except:
                self.logger.warning("Problem saving state variable: " + var)

    def wf_QuickResume(self):
        """
        Resumes the state variable of the previous timestep in memory
        it uses the wf_supplyVariableNamesAndRoles() function to find them.
        The variables are inserted into the model object

        """
        allvars = self._userModel().stateVariables()

        for var in allvars:
            setattr(
                self._userModel(), var, getattr(self._userModel(), var + "_laststep")
            )

        ts = self._userModel().currentTimeStep()
        self._userModel()._setCurrentTimeStep(ts)

        self.DT.update(currentTimeStep=ts)
        self._userModel().currentdatetime = self.DT.currentDateTime
        self.logger.debug(
            "Going one timestep back, redoing: "
            + str(ts)
            + " "
            + str(self.DT.currentDateTime)
        )

    def iniFileSetUp(self, caseName, runId, configfile):
        """
        Reads .ini file and returns a config object.

        Input:
            - caseName - dir with case
            - runId - run dir within case
            - configfile - name of the configfile (.ini type)

        Output:
            - python config object

        """

        config = configparser.ConfigParser()
        config.optionxform = str
        if os.path.exists(os.path.join(caseName, configfile)):
            config.read(os.path.join(caseName, configfile))
        else:
            self.logger.error(
                "Cannot open ini file: " + os.path.join(caseName, configfile)
            )
            sys.exit(1)

        return config

    def wf_setValuesAsNumpy(self, mapname, values):
        """
        set a map with values from a numpy array. Current settings for
        dimensions are assumed. if the name of the maps contains the string "LDD" or "ldd"
        the maps is assumed to be an LDD maps and an lddrepair call is made,
        assume -999 as missing value
        also checks if map contains int at end, then map is part of List object:
        map_0 is part of self.map[0]

        Input:
            - mapname - string with name of map
            - values - numpy array

        :returns: 1 if the map was present, 0 if a new map was created
        """

        arpcr = pcr.numpy2pcr(pcr.Scalar, np.flipud(values).copy(), -999)

        self.setviaAPI[mapname] = 1

        if hasattr(self._userModel(), mapname):

            if "LDD" in mapname.upper():
                setattr(self._userModel(), mapname, pcr.lddrepair(pcr.ldd(arpcr)))
            else:
                if mapname.split('_')[-1].isdigit():
                    getattr(self._userModel(), mapname.split('_')[0])[int(mapname.split('_')[-1])] = arpcr
                else:
                    setattr(self._userModel(), mapname, arpcr)
            return 1
        else:
            self.logger.debug(
                mapname + " is not defined in the usermodel: setting anyway"
            )
            if mapname.split('_')[-1].isdigit():
                getattr(self._userModel(), mapname.split('_')[0])[int(mapname.split('_')[-1])] = arpcr
            else:
                setattr(self._userModel(), mapname, arpcr)
            return 0

    def wf_setValuesAsPcrMap(self, mapname, pcrmap):
        """
        set a map with values from a pcrmap. Current settings for
        dimensions are assumed.
        also checks if map contains int at end, then map is part of List object:
        map_0 is part of self.map[0]

        Input:
            - mapname - string with name of map
            - pcrmap - pcraster map

        :returns: 1 if the map was present, 0 if a new map was created
        """

        arpcr = pcrmap
        self.setviaAPI[mapname] = 1

        if hasattr(self._userModel(), mapname):
            if mapname.split('_')[-1].isdigit():
                getattr(self._userModel(), mapname.split('_')[0])[int(mapname.split('_')[-1])] = arpcr
            else:
                setattr(self._userModel(), mapname, arpcr)
            return 1
        else:
            self.logger.debug(
                mapname + " is not defined in the usermodel: setting anyway"
            )
            if mapname.split('_')[-1].isdigit():
                getattr(self._userModel(), mapname.split('_')[0])[int(mapname.split('_')[-1])] = arpcr
            else:
                setattr(self._userModel(), mapname, arpcr)
            return 0

    def wf_setValues(self, mapname, values):
        """
        set a map with values from a python list or a single scalar
        value. In case a single value is specified the value will be distributed
        uniformly over the map. Current settings for
        dimensions are assumed.
        also checks if map contains int at end, then map is part of List object:
        map_0 is part of self.map[0]

        Input:
            - mapname - string with name of map
            - values - single list of value of length rows * cols or a single
               scalar

        :returns: 1 if the map was present, 0 if a new map was created
        """
        self.setviaAPI[mapname] = 1
        if isinstance(values, list):
            ar = np.array(values)
            ar.reshape(getrows(), getcols())
            arpcr = pcr.numpy2pcr(pcr.Scalar, ar.reshape(getrows(), getcols()), -999)
        else:
            self.logger.debug("Setting single value: " + str(values))
            arpcr = pcr.cover(pcr.scalar(values))

        if hasattr(self._userModel(), mapname):
            if mapname.split('_')[-1].isdigit():
                getattr(self._userModel(), mapname.split('_')[0])[int(mapname.split('_')[-1])] = arpcr
            else:
                setattr(self._userModel(), mapname, arpcr)
            return 1
        else:
            self.logger.debug(
                mapname + " is not defined in the usermodel: setting anyway"
            )
            if mapname.split('_')[-1].isdigit():
                getattr(self._userModel(), mapname.split('_')[0])[int(mapname.split('_')[-1])] = arpcr
            else:
                setattr(self._userModel(), mapname, arpcr)
            return 0

    def wf_setValueRowCol(self, mapname, value, row, col):
        """
        set single value in a map on row, col (0 based). All other values in the
        map remain the same. Numbering starts at the upper left corner.
        also checks if map contains int at end, then map is part of List object:
        map_0 is part of self.map[0]

        Input:
            - mapname - string with name of map
            - row - row to set the value in
            - col - column to set the value in
            - values - single scalar

        :returns: 1 if the map was present, 0 if nothing was done
        """
        self.setviaAPI[mapname] = 1
        if hasattr(self._userModel(), mapname):
            if mapname.split('_')[-1].isdigit():
                ar = pcr.pcr2numpy(getattr(self._userModel(), mapname.split('_')[0])[int(mapname.split('_')[-1])])
            else:
                ar = pcr.pcr2numpy(getattr(self._userModel(), mapname), -999)
            ar[row, col] = value
            arpcr = pcr.numpy2pcr(pcr.Scalar, ar.copy(), -999)
            if mapname.split('_')[-1].isdigit():
                getattr(self._userModel(), mapname.split('_')[0])[int(mapname.split('_')[-1])] = arpcr
            else:
                setattr(self._userModel(), mapname, arpcr)
            return 1
        else:
            self.logger.debug(
                mapname + " is not defined in the usermodel. Doing nothing"
            )
            return 0

    def wf_setValue(self, mapname, value, xcor, ycor):
        """
        set single value in a map on xcor, ycor (0 based). All other values in the
        map remain the same.
        Also checks if map contains int at end, then map is part of List object:
        map_0 is part of self.map[0]

        Input:
            - mapname - string with name of map
            - xcor - xcor to set the value in
            - ycor - ycor to set the value in
            - value - single scalar

        :returns: 1 if the map was present, 0 if nothing was done
        """
        self.setviaAPI[mapname] = 1
        if hasattr(self._userModel(), mapname):
            if mapname.split('_')[-1].isdigit():
                pcrmap = getattr(self._userModel(), mapname.split('_')[0])[int(mapname.split('_')[-1])]
            else:
                pcrmap = getattr(self._userModel(), mapname)
            ar = pcr.pcr2numpy(pcr.scalar(pcrmap), -999)
            row, col = getRowColPoint(pcrmap, xcor, ycor)
            ar[row, col] = value
            save("tt.np", ar)
            pcrmap = pcr.numpy2pcr(pcr.Scalar, ar.copy(), -999)
            pcr.report(pcrmap, "zz.map")
            if mapname.split('_')[-1].isdigit():
                getattr(self._userModel(), mapname.split('_')[0])[int(mapname.split('_')[-1])] = pcrmap
            else:
                setattr(self._userModel(), mapname, pcrmap)
            return 1
        else:
            self.logger.debug(mapname + " is not defined in the usermodel")
            return 0

    def wf_setValueLdd(self, mapname, value, xcor, ycor):
        """
        set single value in an ldd on xcor, ycor (0 based). All other values in the
        map remain the same. Calls lddrepair to ensure the ldd is sound

        Input:
            - mapname of tipy ldd - string with name of map
            - xcor - xcor to set the value in
            - ycor - ycor to set the value in
            - values - single scalar (see pcraster ldddescription for legal values)
                        e.g. use 5 for setting a pit

        :returns: 1 if the map was present, 0 if nothing was done
        """
        self.setviaAPI[mapname] = 1
        if hasattr(self._userModel(), mapname):
            pcrmap = getattr(self._userModel(), mapname)
            ar = pcr.pcr2numpy(pcrmap, -999)
            row, col = getRowColPoint(pcrmap, xcor, ycor)
            ar[row, col] = value
            arpcr = pcr.numpy2pcr(pcr.Scalar, ar.copy(), -999)
            setattr(self._userModel(), mapname, pcr.lddrepair(pcr.ldd(arpcr)))
            return 1
        else:
            self.logger.debug(
                mapname + " is not defined in the usermodel, doing nothing"
            )
            return 0

    def wf_multParameterValues(self, mapname, value):
        """
        multiply a parameter map with a single scalar
        value. Current settings for dimensions are assumed.

        This method must be called *after* the runinitial() method

        Input:
            - mapname - string with name of map
            - value - single scalar

        :returns: 1 if the map was present, 0 if nothing was done
        """

        arpcr = pcr.cover(value)
        self.setviaAPI[mapname] = 1
        if hasattr(self._userModel(), mapname):
            setattr(
                self._userModel(), mapname, arpcr * getattr(self._userModel(), mapname)
            )
            return 1
        else:
            self.logger.debug(
                mapname + " is not defined in the usermodel, doing nothing"
            )
            return 0

    def wf_multParameterValuesArea(self, mapname, value, areacode, areamapname):
        """
        multiply a parameter map with a single scalar
        value for area with id area only. Current settings for dimensions are assumed.

        This method must be called *after* the runinitial() method

        Input:
            - mapname - string with name of map
            - value - single scalar
            - areacode - id of the area in the areamap
            - areamapname - name of the areamap

        :returns: 1 if the map was present, 0 if nothing was done
        """

        arpcr = pcr.cover(value)
        self.setviaAPI[mapname] = 1
        if hasattr(self._userModel(), mapname):
            setattr(
                self._userModel(),
                mapname,
                pcr.ifthenelse(
                    getattr(self._userModel(), areamapname) == str(areacode),
                    arpcr * getattr(self._userModel(), areamapname),
                    getattr(self._userModel(), areamapname),
                ),
            )

            return 1
        else:
            self.logger.debug(
                mapname + " is not defined in the usermodel, doing nothing"
            )
            return 0

    def wf_setParameterValues(self, mapname, values):
        """
        set a parameter map with values from a python list or a single scalar
        value. In case a single value is specified the value will be distributed
        uniformly over the map. Current settings for dimensions are assumed.

        This method must be called _after_ the runinitial() method

        Input:
            - mapname - string with name of map
            - values - single list of value of length rows * cols or a single
               scalar

        :returns: 1 if the map was present, 0 if nothing was done
        """
        self.setviaAPI[mapname] = 1
        if isinstance(values, list):
            ar = np.array(values)

            ar.reshape(getrows(), getcols())
            arpcr = pcr.numpy2pcr(pcr.Scalar, ar.reshape(getrows(), getcols()), -999)
        else:
            arpcr = pcr.cover(values)

        if hasattr(self._userModel(), mapname):
            setattr(self._userModel(), mapname, arpcr)
            return 1
        else:
            self.logger.debug(
                mapname + " is not defined in the usermodel, doing nothing"
            )
            return 0

    def wf_supplyParameterAsList(self, mapname):
        """
        Returns a python list for the specified parameter map and the current
        timestep. If the maps is not dynamic the current status of the map is
        returned.

        Input:
            - mapname (string)

        Output:
            - list
        """
        if hasattr(self._userModel(), mapname):
            retval = pcr.pcr2numpy(getattr(self._userModel(), mapname), -999)
            return retval.flatten().tolist()
        else:
            self.logger.debug(
                mapname + " is not defined in the usermodel, returning empty list"
            )
            return []

    def wf_supplyMapAsList(self, mapname):
        """
        Returns a python list for the specified map and the current
        timestep. If the maps is not dynamic the current status of the map is
        returned which may be undefined for maps that are filled with data
        at the end of a run


        Input:
            - mapname (string)

        Output:
            - list
        """

        if hasattr(self._userModel(), mapname):
            retval = pcr.pcr2numpy(getattr(self._userModel(), mapname), -999)
            if self.APIDebug:
                self.logger.debug("wf_supplyMapAsList returning: " + mapname)

            return retval.flatten().tolist()
        else:
            self.logger.warning(
                mapname + " is not defined in the usermodel, returning empty list"
            )
            return []

    def wf_supplyMapAsNumpy(self, mapname):
        """
        Returns a numpy array (matrix) for the specified map and the current
        timestep. If the maps is not dynamic the current status of the map is
        returns which may be undefined for maps that are filled with data
        at the end of a run.
        Also checks if map contains int at end, then map is part of List object:
        map_0 is part of self.map[0]
        Missing value is -999

        Input:
            - mapname (string)

        Output:
            - numpy array
        """
        
        if hasattr(self._userModel(), mapname):
            if mapname.split('_')[-1].isdigit():
                pcrmap = getattr(self._userModel(), mapname.split('_')[0])[int(mapname.split('_')[-1])]
            else:
                pcrmap = getattr(self._userModel(), mapname)
            if isinstance(pcrmap, pcraster._pcraster.Field):
                tt = pcr.pcr2numpy(pcrmap, -999.0)
                retval = np.flipud(tt).copy()
            else:
                if type(pcrmap) == np.ndarray:
                    retval = pcrmap
                else:
                    retval = np.array(pcrmap)

            if self.APIDebug:
                self.logger.debug("wf_supplyMapAsNumpy returning: " + mapname)
        else:
            self.logger.warning(
                mapname + " is not defined in the usermodel, returning empty array"
            )
            return []

        return retval

    def wf_supplyMapXAsNumpy(self):
        """

        :return x-coordinates of the current clone map:

        Missing value is -999
        """

        x = pcr.xcoordinate((pcr.spatial(pcr.boolean(1.0))))
        retval = pcr.pcr_as_numpy(x).copy()

        return np.flipud(retval).copy()

    def wf_supplyMapYAsNumpy(self):
        """

        :return y-coordinates of the current clone map:

        Missing value is -999
        """

        y = pcr.ycoordinate((pcr.spatial(pcr.boolean(1.0))))
        retval = pcr.pcr_as_numpy(y).copy()

        return np.flipud(retval).copy()

    def wf_supplyMapZAsNumpy(self):
        """

        :return z-coordinates of the current clone map:

        Assumes an Altitude map is present, otherwise return empty numpy
        Missing value is -999
        """

        if hasattr(self._userModel(), "Altitude"):
            retval = getattr(self._userModel(), "Altitude")

            return np.flipud(pcr.pcr2numpy(retval, -999)).copy()
        else:
            self.logger.warning(
                "Altitude is not defined in the usermodel, returning empty list"
            )
            return []

    def wf_supplyMapOrigin(self):
        """

        :return: lower left corner of the map as X, Y

        """
        a = pcr.boolean(1)

        Y = self.wf_supplyMapYAsNumpy()
        X = self.wf_supplyMapXAsNumpy()

        return np.array([X.flatten.min(), Y.flatten.min()])

    def wf_supplyMapAsPcrMap(self, mapname):
        """
        Returns a pcrmap for the specified map and the current
        timestep. If the maps is not dynamic the current staus of the map is
        returns which may be undefined for maps that are filled with data
        at the end of a run
        Missing value is -999

        Input:
            - mapname (string)

        Output:
            - numpy array
        """
        if hasattr(self._userModel(), mapname):
            retval = getattr(self._userModel(), mapname)
            if self.APIDebug:
                self.logger.debug("wf_supplyMapAsNumpy returning: " + mapname)
        else:
            self.logger.warning(
                mapname + " is not defined in the usermodel, returning empty list"
            )
            return []

        return retval

    def wf_supplyGridDim(self):
        """
        return the dimension of the current model grid as list::

         [ Xul, Yul, xsize, ysize, rows, cols, Xlr, Ylr]
        """

        return getgridparams()

    def wf_supplyVariableNamesAndRoles(self):
        """
        Returns a list of variables
        List of list with the following structure::

            [[ name, role, unit]
            [ name, role, unit]
            ...
            ]
            role: 0 = input (to the model)
                  1 = is output (from the model)
                  2 = input/output (state information)
                  3 = model parameter
            unit: 0 = mm/timestep
                  1 = m^3/sec
                  2 = m
                  3 = degree Celcius
                  4 = mm
                  5 = -

        The first time this function is called the exchangeitems object is filled
        with data from the ini file.
        """

        res = self.exchnageitems.getvars()

        # Fill object with data from ini file
        # TODO: clean up!!
        if len(res) == 0:
            API = configsection(self._userModel().config, "API")
            for a in API:
                tt = []
                line = self._userModel().config.get("API", a)
                tt.append(a)
                tt.append(int(line.split(",")[0]))
                tt.append((line.split(",")[1]))
                res.append(tt)
                self.exchnageitems.addvar(tt[0], tt[1], tt[2])

        if hasattr(self._userModel(), "supplyVariableNamesAndRoles"):
            if self.APIDebug:
                res = self._userModel().supplyVariableNamesAndRoles()
                self.logger.debug(
                    "wf_supplyVariableNamesAndRoles from usermodel: " + str(res)
                )
            return res
        else:
            if self.APIDebug:
                self.logger.debug(
                    "wf_supplyVariableNamesAndRoles from framework: " + str(res)
                )
            return res

    def wf_supplyVariableNames(self):
        """
        returns:
            - the a list of variable names
        """
        varlist = self.wf_supplyVariableNamesAndRoles()
        ret = list(range(len(varlist)))
        for ss in range(len(varlist)):
            ret[ss] = varlist[ss][0]

        if self.APIDebug:
            self.logger.debug("wf_supplyVariableNames from framework: " + str(ret))

        return ret

    def wf_supplyVariableRoles(self):
        """
        returns:
            - the a list of variable roles
        """
        varlist = self.wf_supplyVariableNamesAndRoles()
        ret = list(range(len(varlist)))
        for ss in range(len(varlist)):
            ret[ss] = varlist[ss][1]

        if self.APIDebug:
            self.logger.debug("wf_supplyVariableRoles from framework: " + str(ret))

        return ret

    def wf_supplyVariableCount(self):
        """
        returns:
            - the number of exchangable variables
        """
        varlist = self.wf_supplyVariableNamesAndRoles()

        if self.APIDebug:
            self.logger.debug(
                "wf_supplyVariableCount from framework: " + str(len(varlist))
            )

        return len(varlist)

    def wf_supplyVariableUnits(self):
        """
        returns:
            - the a list of variable units
        """
        varlist = self.wf_supplyVariableNamesAndRoles()
        ret = list(range(len(varlist)))
        for ss in range(len(varlist)):
            ret[ss] = varlist[ss][2]

        if self.APIDebug:
            self.logger.debug("wf_supplyVariableUnits from framework: " + str(ret))

        return ret

    def wf_supplyEndTime(self):
        """
        gets the end time of the model run
        :return: current time as seconds since epoch
        """
        seconds_since_epoch = calendar.timegm(self.DT.runEndTime.utctimetuple())

        return seconds_since_epoch

    def wf_supplyStartTime(self):
        """
        gets the start time of the model run
        :return: current time as seconds since epoch
        """
        seconds_since_epoch = calendar.timegm(self.DT.runStartTime.utctimetuple())

        return seconds_since_epoch

    def wf_supplyCurrentTime(self):
        """
        gets the current time in seconds after the start of the run
        Assumed daily timesteps if not defined in the user model

        Output:
           - current model time (since start of the run)

        """
        dtt = self.DT.currentDateTime
        seconds_since_epoch = calendar.timegm(dtt.utctimetuple())

        return seconds_since_epoch

    def wf_supplyCurrentDateTime(self):
        """
        gets the current time in seconds after the start of the run
        Assumed daily timesteps if not defined in the user model

        Output:
           - current model time (since start of the run)

        """
        dtt = self.DT.currentDateTime

        return dtt

    def wf_supplyStartDateTime(self):
        """
        gets the current time in seconds after the start of the run
        Assumed daily timesteps if not defined in the user model

        Output:
           - current model time (since start of the run)

        """
        dtt = self.DT.runStartTime

        return dtt

    def wf_supplyEpoch(self):
        """
        Supplies the time epoch as a CF string
        Output:
           - current model time (since start of the run)
        """
        epoch = time.gmtime(0)

        epochstr = "seconds since %04d-%02d-%02d %02d:%02d:%02d.0 00:00" % (
            epoch.tm_year,
            epoch.tm_mon,
            epoch.tm_mday,
            epoch.tm_hour,
            epoch.tm_min,
            epoch.tm_sec,
        )
        return epochstr

    def wf_supplyRowCol(self, mapname, xcor, ycor):
        """
        returns a tuple (Row,Col) for the given X and y coordinate

        Input:
            - mapname
            - xcor
            - ycor

        Output:
            - tuple with row, col
        """
        pt = getRowColPoint(mapname, xcor, ycor)
        if self.APIDebug:
            self.logger.debug("wf_supplyRowCol from framework: " + str(pt))

        return pt

    def wf_supplyScalar(self, mapname, xcor, ycor):
        """
        returns a single value for the x and y coordinates from the
        map given uses getValAtPoint(in_map,xcor,ycor) from terrain_lib.py

        Input:
            - mapname
            - xcor
            - ycor

        Output:
            - value at location xcor, ycor

        """
        pcmap = getattr(self._userModel(), mapname)
        pt = getValAtPoint(pcmap, xcor, ycor)

        if self.APIDebug:
            self.logger.debug("wf_supplyScalar from framework: " + str(pt))

        return pt

    def wf_supplyScalarRowCol(self, mapname, row, col):
        """
        returns a single value for row and col from the
        map given (zero based).

        Input:
            - mapname
            - xcor
            - ycor

        Output:
            - value at location row, col
        """
        pcmap = getattr(self._userModel(), mapname)
        ret = pcr.cellvalue(pcmap, row + 1, col + 1)

        return ret[0]

    def _userModel(self):
        """ Returns the class provided by the user """
        return self._d_model

    def _runDynamic(self, firststep, laststep):
        """
        Runs the dynamic model from firststep to laststep

        Input:

            :ivar firststep: first timestep of the model run
            :ivar laststep: last timestep of the model run

        """
        self._userModel()._setInDynamic(True)

        #
        if firststep == 0:
            step = self._d_firstTimestep
        else:
            step = firststep

        if laststep == 0:
            laststep = self._d_lastTimestep

        self._userModel()._setNrTimeSteps(int(laststep))

        self.DT.update(
            currentTimeStep=self.DT.currentTimeStep, mode=self.runlengthdetermination
        )

        while step <= self._userModel().nrTimeSteps():
            self._incrementIndentLevel()
            self._atStartOfTimeStep(step)
            self._userModel()._setCurrentTimeStep(step)

            if hasattr(self._userModel(), "dynamic"):
                self._incrementIndentLevel()
                self._traceIn("dynamic")
                self._userModel().dynamic()
                self._traceOut("dynamic")
                self._decrementIndentLevel()
                # Save state variables in memory
                self.wf_QuickSuspend()

            # Make the summary variables
            for a in range(0, len(self.statslst)):
                data = getattr(self._userModel(), self.statslst[a].varname)
                self.statslst[a].add_one(data)

            # Online statistics (rolling mean for now)
            for key in self.onlinestat.statvarname:
                stvar = self.onlinestat.getstat(getattr(self._userModel(), key), key)
                # stvar = self.onlinestat.getstat(pcr.cover(self.DT.currentTimeStep * 1.0), key)
                setattr(self._userModel(), self.onlinestat.statvarname[key], stvar)

            # Increment one timesteps
            self.DT.update(incrementStep=True, mode=self.runlengthdetermination)
            self._userModel().currentdatetime = self.DT.currentDateTime

            self.wf_savedynMaps()
            self.wf_saveTimeSeries()

            self.logger.debug(
                "timestep: "
                + str(self._userModel().currentTimeStep())
                + "/"
                + str(self.DT.lastTimeStep)
                + " ("
                + str(self.DT.currentDateTime)
                + ")"
            )

            self._timeStepFinished()
            self._decrementIndentLevel()
            step += 1
            self.setviaAPI = {}

        self._userModel()._setInDynamic(False)

    ## \brief Re-implemented from ShellScript.
    #
    # Runs a dynamic user model.
    def run(self):
        """ Runs the dynamic model for all timesteps """

        self._atStartOfScript()
        if hasattr(self._userModel(), "resume"):
            if self._userModel().firstTimeStep() == 1:
                self._runInitial()
            else:
                self._runResume()
        else:
            self._runInitial()

        self._runDynamic()

        # only execute this section while running filter frameworks
        if hasattr(self._userModel(), "suspend") and hasattr(
            self._userModel(), "filterPeriod"
        ):
            self._runSuspend()

        return 0

    def reportStatic(self, variable, name, style=1, gzipit=False, longname=None):
        """

        :param variable:
        :param name:
        :param style:
        :param gzipit:
        :param longname:
        :return:
        """
        if longname == None:
            longname = name
        path = name

        if self.outputFormat == 1:

            if not hasattr(self, "NcOutputStatic"):
                pcr.report(variable, path)
                if gzipit:
                    Gzip(path, storePath=True)
            else:
                self.NcOutputStatic.savetimestep(1, variable, var=name, name=longname)

        elif self.outputFormat == 2:
            np.savez(path, pcr.pcr2numpy(variable, -999))
        elif self.outputFormat == 3:
            np.savez(path, pcr.pcr2numpy(variable, -999))
        elif self.outputFormat == 4:
            np.savetxt(path, pcr.pcr2numpy(variable, -999), fmt="%0.6g")

    def reportState(self, variable, name, style=1, gzipit=False, longname=None):
        """

        :param variable:
        :param name:
        :param style:
        :param gzipit:
        :param longname:
        :return:
        """
        if longname == None:
            longname = name
        path = name

        if self.outputFormat == 1:
            if not hasattr(self, "NcOutputState"):
                pcr.report(variable, path)
                if gzipit:
                    Gzip(path, storePath=True)
            else:
                self.NcOutputState.savetimestep(1, variable, var=name, name=longname)

        elif self.outputFormat == 2:
            np.savez(path, pcr.pcr2numpy(variable, -999))
        elif self.outputFormat == 3:
            np.savez(path, pcr.pcr2numpy(variable, -999))
        elif self.outputFormat == 4:
            np.savetxt(path, pcr.pcr2numpy(variable, -999), fmt="%0.6g")

    def _reportNew(self, variable, name, style=1, gzipit=False, longname=None):
        """
        outputformat: (set in the [framework] section of the init file).
            1: pcraster
            2: numpy (compressed)
            3: matlab
            4: numpy text files (large and slow)

            ..
            # Example:

            [framework]
            outputformat = 4
        """

        if longname == None:
            longname = name
        head, tail = os.path.split(name)

        # if re.search("\.", tail):
        #  msg = "File extension given in '" + name + "' not allowed, provide filename without extension"
        #  raise FrameworkError(msg)

        directoryPrefix = ""
        nameSuffix = ".map"
        newName = ""

        if hasattr(self._userModel(), "_inStochastic"):
            if self._userModel()._inStochastic():
                if self._userModel()._inPremc():
                    newName = name + nameSuffix
                elif self._userModel()._inPostmc():
                    newName = name + nameSuffix
                else:
                    directoryPrefix = str(self._userModel().currentSampleNumber())

        if self._userModel()._inInitial():
            newName = name + nameSuffix

        if hasattr(self._userModel(), "_inDynamic"):
            if self._userModel()._inDynamic() or self._inUpdateWeight():
                newName = pcraster.framework.generateNameT(
                    name, self._userModel().currentTimeStep()
                )

        if newName == "":  # For files from suspend
            newName = name

        path = os.path.join(directoryPrefix, newName)

        if self.outputFormat == 1:
            if not hasattr(self, "NcOutput"):
                pcr.report(variable, path)
                if gzipit:
                    Gzip(path, storePath=True)
            else:
                self.NcOutput.savetimestep(
                    self._userModel().currentTimeStep(),
                    variable,
                    var=name,
                    name=longname,
                )

        elif self.outputFormat == 2:
            np.savez(path, pcr.pcr2numpy(variable, -999))
        elif self.outputFormat == 3:
            np.savez(path, pcr.pcr2numpy(variable, -999))
        elif self.outputFormat == 4:
            np.savetxt(path, pcr.pcr2numpy(variable, -999), fmt="%0.6g")

    def wf_readmapClimatology(self, name, kind=1, default=0.0, verbose=1):
        """
        Read a climatology map. The current date/time is converted to:
        1: a month and the file for the current month is returned
        2: days of year and the file for the current day is returned
        3: hour of day and the file for the current hours is returned

        :param name: name if the mapstack
        :param kind: type of the climatology
        :return: a map
        """

        # Assume the variable is via the API (replaces the
        if os.path.basename(name) in self.setviaAPI:
            self.setviaAPI.pop(os.path.basename(name))
            self.logger.debug(
                os.path.basename(name)
                + " set via API, not reading from file, using memory copy"
            )
            return getattr(self._userModel(), os.path.basename(name))

        # TODO: Add support for netcdf files
        directoryPrefix = ""
        if kind == 1:
            month = self.DT.currentDateTime.month
            newName = pcraster.framework.generateNameT(name, month)
            path = os.path.join(directoryPrefix, newName)
            if os.path.isfile(path):
                mapje = pcr.readmap(path)
                return mapje
            else:
                if verbose:
                    self.logger.warning(
                        "Climatology data ("
                        + path
                        + ") for timestep not present, returning "
                        + str(default)
                    )

                return pcr.scalar(default)
        elif kind == 2:
            yday = self.DT.currentDateTime.timetuple().tm_yday
            newName = pcraster.framework.generateNameT(name, yday)
            path = os.path.join(directoryPrefix, newName)
            if os.path.isfile(path):
                mapje = pcr.readmap(path)
                return mapje
            else:
                if verbose:
                    self.logger.warning(
                        "Climatology data ("
                        + path
                        + ") for timestep not present, returning "
                        + str(default)
                    )

                return pcr.scalar(default)

        else:
            self.logger.error(
                "This Kind of climatology not implemented yet: " + str(kind)
            )

    def wf_readmap(
        self,
        name,
        default,
        verbose=True,
        fail=False,
        ncfilesource="not set",
        silent=False,
    ):
        """
          Adjusted version of readmapNew. the style variable is used to indicated
          how the data is read::

              1 - default: reads pcrmaps
              2 - memory: assumes the map is made available (in memory) using
              the in-memory interface


        """
        directoryPrefix = ""
        nameSuffix = ".map"
        newName = ""

        varname = os.path.basename(name)

        # find if this is an exchnageitem
        thevars = self.exchnageitems.getvars()
        if len(thevars) == 0:
            self.wf_supplyVariableNamesAndRoles()

        style = self.exchnageitems.getvarStyle(varname)

        # set this for initial (before the model is actually running)
        if os.path.splitext(name)[1] == ".map":
            newName = name
        else:
            newName = name + nameSuffix

        # Assume the variable is via the API (replaces the
        if os.path.basename(name) in self.setviaAPI:
            self.setviaAPI.pop(os.path.basename(name))
            self.logger.debug(
                os.path.basename(name)
                + " set via API, not reading from file, using memory copy"
            )
            return getattr(self._userModel(), os.path.basename(name))

        if hasattr(self._userModel(), "_inStochastic"):
            if self._userModel()._inStochastic():
                if self._userModel()._inPremc() or self._userModel()._inPostmc():
                    newName = name + nameSuffix
                else:
                    directoryPrefix = str(self._userModel().currentSampleNumber())

        if hasattr(self._userModel(), "_inInitial"):
            if self._userModel()._inInitial():
                if os.path.splitext(name)[1] == ".map":
                    newName = name
                else:
                    newName = name + nameSuffix

        if self._inResume():
            if os.path.splitext(name)[1] == ".map":
                newName = name
            else:
                newName = name + nameSuffix

        if hasattr(self._userModel(), "_inDynamic"):
            if self._userModel()._inDynamic() or self._inUpdateWeight():
                timestep = self._userModel().currentTimeStep()
                # print timestep
                if "None" not in self.ncfile:
                    newName = name
                else:
                    newName = pcraster.framework.generateNameT(name, timestep)

        if style == 1:  # Normal reading of mapstack from DISK per via or via netcdf
            path = os.path.join(directoryPrefix, newName)
            assert path is not ""

            if self._userModel()._inDynamic():
                if "None" not in self.ncfile:
                    retval, succ = self.NcInput.gettimestep(
                        self._userModel().currentTimeStep(),
                        self.logger,
                        tsdatetime=self.DT.nextDateTime,
                        var=varname,
                        shifttime=self.DT.startadjusted,
                    )
                    if succ:
                        return retval
                    else:
                        return pcr.cover(pcr.scalar(default))

                if os.path.isfile(path):
                    mapje = pcr.readmap(path)
                    return mapje
                else:
                    if verbose:
                        self.logger.debug(
                            "Input data ("
                            + os.path.abspath(path)
                            + ") for timestep not present, returning "
                            + str(default)
                        )
                    if fail:
                        self.logger.error(
                            "Required map: "
                            + os.path.abspath(path)
                            + " not found, exiting.."
                        )
                        raise ValueError("Input map not found")
                    return pcr.cover(pcr.scalar(default))

            elif self._userModel()._inInitial():
                if "None" not in self.ncfilestatic:
                    retval, succ = self.NcInputStatic.gettimestep(
                        1, self.logger, var=varname
                    )
                    if succ:
                        return retval
                    else:
                        if fail:
                            if not silent:
                                self.logger.error(
                                    "Required map: "
                                    + os.path.abspath(path)
                                    + " not found in "
                                    + self.ncfilestatic
                                    + "  exiting.."
                                )
                            raise ValueError(
                                "Input static variable not found in netcdf"
                            )
                        else:
                            return self.TheClone + default

                if os.path.isfile(path):
                    mapje = pcr.readmap(path)
                    return mapje
                else:
                    if verbose:
                        self.logger.debug(
                            "Static input data ("
                            + os.path.abspath(path)
                            + ")  not present, returning "
                            + str(default)
                        )
                    if fail:
                        if not silent:
                            self.logger.error(
                                "Required map: "
                                + os.path.abspath(path)
                                + " not found, exiting.."
                            )
                        raise ValueError("Input static variable not found")
                    return self.TheClone + default

            elif self._inResume():
                if ncfilesource == self.ncinfilestates and ncfilesource not in "None":
                    retval, succ = self.NcInputStates.gettimestep(
                        1, self.logger, var=varname, tsdatetime=self.DT.runStateTime
                    )
                    if succ:
                        return retval
                    else:
                        if fail:
                            if not silent:
                                self.logger.error(
                                    "Required map: "
                                    + os.path.abspath(path)
                                    + " not found, exiting.."
                                )
                            raise ValueError("Input state variable not found")

                        return self.TheClone + default

                if os.path.isfile(path):
                    mapje = pcr.readmap(path)
                    return mapje
                else:
                    if verbose:
                        self.logger.debug(
                            "State input data ("
                            + os.path.abspath(path)
                            + ")  not present, returning "
                            + str(default)
                        )
                    if fail:
                        if not silent:
                            self.logger.error(
                                "Required map: "
                                + os.path.abspath(path)
                                + " not found, exiting.."
                            )
                        raise ValueError("Input state variable not found")
                    return pcr.cover(pcr.scalar(default))
            else:  # Assuming we are in pre-or post loop within the framwork
                if "None" not in self.ncfilestatic:
                    retval, succ = self.NcInputStatic.gettimestep(
                        1, self.logger, var=varname
                    )
                    if succ:
                        return retval
                    else:
                        if fail:
                            if not silent:
                                self.logger.error(
                                    "Required map: "
                                    + os.path.abspath(path)
                                    + " not found in "
                                    + self.ncfilestatic
                                    + "  exiting.."
                                )
                            raise ValueError("Input variable not found in netcdf")
                        else:
                            return self.TheClone + default

                if os.path.isfile(path):
                    mapje = pcr.readmap(path)
                    return mapje
                else:
                    if verbose:
                        self.logger.debug(
                            "Static input data ("
                            + os.path.abspath(path)
                            + ")  not present, returning "
                            + str(default)
                        )
                    if fail:
                        if not silent:
                            self.logger.error(
                                "Required map: "
                                + os.path.abspath(path)
                                + " not found, exiting.."
                            )
                        raise ValueError("Input variable not found")
                    return self.TheClone + default

        elif style == 2:  # Assuming they are set in memory by the API
            #
            # first get basename (last bit of path)
            name = os.path.basename(name)
            if hasattr(self._userModel(), name):
                return pcr.cover(getattr(self._userModel(), name), pcr.scalar(default))
            else:
                self.logger.warning(
                    "Variable: " + name + " not set by API, returning default"
                )
                setattr(self._userModel(), name, pcr.cover(pcr.scalar(default)))
                return getattr(self._userModel(), name)
        else:
            self.logger.warning(
                "Unknown style ("
                + str(style)
                + ") for variable: "
                + name
                + ", returning default"
            )
            return self.TheClone + default

    ## \brief testing the requirements for the dynamic framework
    #
    # To use the dynamic framework the user must implement the following methods
    # in this class:
    # - either "run" or "initial" and "dynamic"
    def _testRequirements(self):
        if hasattr(self._userModel(), "_userModel"):
            msg = "The _userModel method is deprecated and obsolete"
            self.showWarning(msg)

        if not hasattr(self._userModel(), "dynamic") and not hasattr(
            self._userModel(), "run"
        ):
            msg = "Cannot run dynamic framework: Implement dynamic method"
            raise frameworkBase.FrameworkError(msg)

        if not hasattr(self._userModel(), "initial"):
            if self._debug():
                self.showWarning("No initial section defined.")

        if not hasattr(self._userModel(), "stateVariables"):
            if self._debug():
                self.showWarning("No stateVariables defined in usermodel.")

    def setQuiet(self, quiet):
        self._d_quiet = quiet
        self._d_quietProgressDots = quiet
