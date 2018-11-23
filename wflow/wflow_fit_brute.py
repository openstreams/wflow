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
Fit a wflow\_ hydrological model using scipy.leastsq.

usage

::
    
    wflow_fit -M ModelName [-h][-F runinfofile][-C casename]
          [-c configfile][-T last_step][-S first_step][-s seconds]
          

    -M: model to fit (e.g. wflow_sbm, wflow_hbv, wflow_cqf)          
                      
    -T: Set last timestep
    
    -S: Set the start timestep (default = 1)
    
    -C: set the name  of the case (directory) to run
    
    -R: set the name runId within the current case
    
    -U: save the map after each step ti the input (staticmaps) dir so 
        that next steps (colums) use the previous results
    
    -c: name of wflow the configuration file (default: Casename/wflow_sbm.ini). 
    
    -h: print usage information
    
   
    
For this program to work you must add a [fit] section to the
ini file of the program to fit (e.g. the wflow\_hbv program)



$Author: schelle $
$Id: wflow_fit_brute.py 785 2013-09-16 12:50:31Z schelle $
$Rev: 785 $

"""


import csv
import getopt
import os.path
import sys

import numpy as np
import pylab
import wflow.pcrut as pcrut
import wflow.stats as stats


# TODO: do not read results from file
# TODO: allow framework to be silent (no debug lines)
# TODO: filter on the points to use (high/low/regression etc)
# See http://www.netlib.org/minpack/lmdif.f


def configget(config, section, var, default):
    """
    gets parameter from config file and returns a default value
    if the parameter is not found
    """
    try:
        ret = config.get(section, var)
    except:
        print("returning default (" + default + ") for " + section + ":" + var)
        ret = default

    return ret


class wfmodel_fit_API:
    """
    Class that initializes and runs a wflow model
    """

    def __init__(
        self,
        startTime,
        stopTime,
        casename,
        runId="_fitrun",
        modeltofit="wflow_sbm",
        config="wflow_sbm.ini",
        clonemap="wflow_subcatch.map",
    ):

        # try:
        # self.WF = __import__(modeltofit, globals(), locals(), [], -1)
        # except:
        self.NS = []
        self.BIAS = []
        self.CORR = []
        self.MABSE = []
        self.range = {}

        try:
            mod = __import__("wflow." + modeltofit, globals(), locals(), [], -1)
            self.WF = mod.wflow_sbm
        except:
            mod = __import__(modeltofit, globals(), locals(), [], -1)
            self.WF = mod

        self.results = []
        self.runId = runId
        self.caseName = casename
        self.stopTime = stopTime
        self.startTime = startTime
        configfile = config
        self.pars = []
        self.calibpars = []
        wflow_cloneMap = clonemap

        self.myModel = self.WF.WflowModel(
            wflow_cloneMap, self.caseName, self.runId, configfile
        )
        # initialise the framework
        self.dynModelFw = self.WF.wf_DynamicFramework(
            self.myModel, self.stopTime, self.startTime
        )
        # Load model config from files and check directory structure
        self.dynModelFw.createRunId(NoOverWrite=False, level=30)
        # self.dynModelFw.logger.setLevel(20)
        self.log = self.dynModelFw._userModel().logger
        self.conf = self.dynModelFw._userModel().config
        self.log.log(45, "Initialising fit module...")

        # Gets all para_0 to n parameters to be fitted
        #!TODO: add area code to parameter, only do that area
        # TODO: add columns in the measued file to fit to (shoudl mach column in the simulated file)
        item = 0
        while configget(self.conf, "fit", "parameter_" + str(item), "__") != "__":
            parstr = configget(self.conf, "fit", "parameter_" + str(item), "M").split(
                ":"
            )
            if len(parstr) == 2:
                par = parstr[0]
                exec("ar = np.array(" + parstr[1] + ")")
                self.range[par] = ar
            else:
                par = parstr[0]
                self.range = None

            self.calibpars.append(par)
            self.pars.append(1.0)
            item = item + 1
        self.qmeasname = configget(self.conf, "fit", "Q", "calib.tss")
        self.epsfcn = float(configget(self.conf, "fit", "epsfcn", "0.00001"))
        self.ftol = float(configget(self.conf, "fit", "ftol", "0.0001"))
        self.xtol = float(configget(self.conf, "fit", "xtol", "0.0001"))
        self.gtol = float(configget(self.conf, "fit", "gtol", "0.0001"))
        self.factor = float(configget(self.conf, "fit", "factor", "100.0"))

        exec("self.ColSimS  = " + configget(self.conf, "fit", "ColSim", "[1]"))
        exec("self.ColMeasS = " + configget(self.conf, "fit", "ColMeas", "[1]"))
        self.WarmUpSteps = int(configget(self.conf, "fit", "WarmUpSteps", "1"))
        self.AreaMapName = configget(self.conf, "fit", "areamap", "wflow_catchment.map")
        self.AreaMap = self.WF.readmap(os.path.join(self.caseName, self.AreaMapName))
        exec("self.AreaCodeS = " + configget(self.conf, "fit", "areacode", "[1]"))

        # Shift columns as the maps are one bases and the cols 0 based
        i = 0
        for a in self.ColSimS:
            self.ColSimS[i] = self.ColSimS[i] - 1
            self.ColMeasS[i] = self.ColMeasS[i] - 1
            i = i + 1

        self.ColSim = self.ColSimS[0]
        self.ColMeas = self.ColMeasS[0]
        self.AreaCode = self.AreaCodeS[0]

    def multVarWithPar(self, pars):
        """
        Multiply a parameter in the model with the fit parameters.
        Use a map to limit the area to adjust
        """
        i = 0
        for j in pars:
            self.log.info(
                "Areacode: "
                + str(self.AreaCode)
                + " Multiplying parameter: "
                + self.calibpars[i]
                + " with: "
                + str(j)
            )
            # self.dynModelFw.wf_multParameterValues(self.calibpars[i],j)
            themappcr = self.dynModelFw.wf_supplyMapAsPcrMap(self.calibpars[i])
            zz = self.WF.ifthenelse(
                self.AreaMap == int(self.AreaCode),
                self.WF.boolean(1),
                self.WF.boolean(0),
            )
            # self.WF.report(zz,self.calibpars[i] + "_area.map")
            themappcr = self.WF.ifthenelse(
                self.AreaMap == int(self.AreaCode), themappcr * j, themappcr
            )
            # self.WF.report(themappcr,self.calibpars[i] + str(j) + ".map")
            self.dynModelFw.wf_setValuesAsPcrMap(self.calibpars[i], themappcr)
            i = i + 1

    def saveinitpars(self):
        self.dynModelFw._runInitial()  # Runs initial part
        i = 0
        for j in self.pars:
            self.log.info("Saving parameter (initial values): " + self.calibpars[i])
            strr_org = (
                "self.WF.report(self.dynModelFw._userModel()."
                + self.calibpars[i]
                + ',"'
                + self.caseName
                + "/"
                + self.runId
                + "/"
                + self.calibpars[i]
                + '_org.map")'
            )
            exec(strr_org)
            i = i + 1

    def run(self, pars):
        """
        Run the model for the number of timesteps.
        """

        # Run the initial part of the model (reads parameters and sets initial values)
        self.dynModelFw._runInitial()  # Runs initial part
        # self.dynModelFw.wf_multParameterValues('M',pars[0])
        self.multVarWithPar(pars)

        self.dynModelFw._runResume()  # gets the state variables

        for ts in range(self.startTime, self.stopTime):
            self.dynModelFw._runDynamic(ts, ts)  # runs for all timesteps

        # save output state
        self.dynModelFw._runSuspend()
        self.dynModelFw._wf_shutdown()
        results, head = pcrut.readtss(
            os.path.join(self.caseName, self.runId, "run.tss")
        )
        return results[self.WarmUpSteps :, self.ColSim].astype(np.float64)

    def savemaps(self, pars, savetoinput=False):
        """
        Ssave the adjusted (and original) parameter maps
        
        """

        # To get the original values of the parameters
        self.dynModelFw._runInitial()
        # !!!!!!!!!! Not sure if the last version of the par is the best fit!!
        i = 0
        for j in pars:
            self.log.log(45, "Saving parameter: " + self.calibpars[i])
            exec("newmap = self.dynModelFw._userModel()." + self.calibpars[i])
            newmap = self.WF.ifthenelse(
                self.AreaMap == self.AreaCode, newmap * j, newmap
            )
            strr_new = (
                "self.WF.report(newmap,"
                + '"'
                + self.caseName
                + "/"
                + self.runId
                + "/"
                + self.calibpars[i]
                + "_"
                + str(self.ColSim)
                + "_"
                + str(self.ColMeas)
                + "_"
                + str(self.AreaCode)
                + '.map")'
            )
            if savetoinput:
                self.log.log(45, "Saving adjusted map to input!!")
                str_save = (
                    "self.WF.report(newmap,"
                    + '"'
                    + self.caseName
                    + "/staticmaps/"
                    + self.calibpars[i]
                    + '.map")'
                )
                exec(str_save)

            exec(strr_new)
            i = i + 1

    def shutdown(self):
        """
        Shutdown the model 
        
        """

        self.dynModelFw._wf_shutdown()


def bruteforce(mimo, qmeas, caseName, runId):
    """
    """
    nowpars = mimo.pars
    print(nowpars)
    i = 0
    results = []

    execstr = "mimo.range[mimo.calibpars[" + str(i) + "]]"
    print("=========")
    for i in range(0, len(mimo.calibpars)):
        if i > 0:
            execstr = execstr + ", mimo.range[mimo.calibpars[" + str(i) + "]]"
        print(execstr)
        for parmult in mimo.range[mimo.calibpars[i]]:
            nowpars[i] = parmult

    exec("combi = itertools.product(" + execstr + ")")

    while 1:
        try:
            thisparset = next(combi)
            print("Parameters:" + str(thisparset))
            q = mimo.run(thisparset)

            res = q - qmeas
            # only resturn non-nan values
            resnonan = res[~np.isnan(res)]
            mimo.log.log(45, "Parameters now: " + str(thisparset))
            pylab.plot(q)
            mimo.NS.append(
                stats.get_nash_sutcliffe(
                    qmeas[~np.isnan(res)], q[~np.isnan(res)], NoData=np.nan
                )
            )
            mimo.BIAS.append(
                stats.get_bias(qmeas[~np.isnan(res)], q[~np.isnan(res)], NoData=np.nan)
            )
            mimo.CORR.append(
                stats.get_correlation(
                    qmeas[~np.isnan(res)], q[~np.isnan(res)], NoData=np.nan
                )
            )
            mimo.MABSE.append(
                stats.get_mean_absolute_error(
                    qmeas[~np.isnan(res)], q[~np.isnan(res)], NoData=np.nan
                )
            )
            mimo.log.log(45, "NS: " + str(mimo.NS[-1]))
            mimo.log.log(45, "BIAS: " + str(mimo.BIAS[-1]))
            mimo.log.log(45, "CORR: " + str(mimo.CORR[-1]))
            mimo.log.log(45, "MABSE: " + str(mimo.MABSE[-1]))
            results.append(
                [
                    thisparset,
                    stats.get_nash_sutcliffe(
                        qmeas[~np.isnan(res)], q[~np.isnan(res)], NoData=np.nan
                    ),
                    stats.get_bias(
                        qmeas[~np.isnan(res)], q[~np.isnan(res)], NoData=np.nan
                    ),
                ]
            )

            pylab.savefig(os.path.join(caseName, runId, str(mimo.ColSim) + "fit.png"))
            # zz = errfuncFIT(thisparset,qmeas,mimo,caseName,runId)
        except:
            print("Unexpected error:", sys.exc_info()[0])
            break

    return results


def errfuncFIT(pars, qmeas, mimo, caseName, runId):
    q = mimo.run(pars)
    res = q - qmeas
    # only resturn non-nan values
    resnonan = res[~np.isnan(res)]

    mimo.log.log(45, "Parameters now: " + str(pars))
    pylab.plot(q)

    mimo.NS.append(
        stats.get_nash_sutcliffe(
            qmeas[~np.isnan(res)], q[~np.isnan(res)], NoData=np.nan
        )
    )
    mimo.BIAS.append(
        stats.get_bias(qmeas[~np.isnan(res)], q[~np.isnan(res)], NoData=np.nan)
    )
    mimo.CORR.append(
        stats.get_correlation(qmeas[~np.isnan(res)], q[~np.isnan(res)], NoData=np.nan)
    )
    mimo.MABSE.append(
        stats.get_mean_absolute_error(
            qmeas[~np.isnan(res)], q[~np.isnan(res)], NoData=np.nan
        )
    )
    mimo.log.log(45, "NS: " + str(mimo.NS[-1]))
    mimo.log.log(45, "BIAS: " + str(mimo.BIAS[-1]))
    mimo.log.log(45, "CORR: " + str(mimo.CORR[-1]))
    mimo.log.log(45, "MABSE: " + str(mimo.MABSE[-1]))

    pylab.savefig(os.path.join(caseName, runId, str(mimo.ColSim) + "fit.png"))
    return resnonan


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args:
        print(msg)
    print(__doc__)
    sys.exit(0)


def printresults(pp, a, b, c, d, calibpars, fname, model):

    ff = open(fname, "w")

    i = 0
    print("Optimised parameter multiplication values:", file=ff)
    if np.iterable(pp):
        for par in pp:
            print("Parameter " + calibpars[i] + " = " + str(par), file=ff)
            i = i + 1
    else:
        print("Parameter " + calibpars[0] + " = " + str(pp), file=ff)

    print("Estimate of the jacobian around the solution: " + str(a), file=ff)
    for dtc in b:
        print(dtc + " = " + str(b[dtc]), file=ff)

    if d in [1, 2, 3, 4]:
        print("A solution was found (" + str(d) + ")", file=ff)
        print(c, file=ff)
    else:
        print("No solution was found (" + str(d) + ")", file=ff)
        print(c, file=ff)

    print("NS: " + str(model.NS), file=ff)
    print("BIAS: " + str(model.BIAS), file=ff)
    print("CORR: " + str(model.CORR), file=ff)
    print("MABSE: " + str(model.MABSE), file=ff)
    ff.close()


def main(argv=None):

    caseName = "not_set"
    _lastTimeStep = 10
    _firstTimeStep = 1
    fitname = "wflow_fit.res"
    runId = "_fitrun"
    # theModel = 'wflow_cqf'
    theModel = "wflow_sbm"
    configfile = None
    saveResults = False
    fitmethod = "fmin"

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return

    opts, args = getopt.getopt(argv, "C:S:T:c:s:R:hM:U")

    for o, a in opts:
        if o == "-C":
            caseName = a
        if o == "-c":
            configfile = a
        if o == "-s":
            timestepsecs = int(a)
        if o == "-T":
            _lastTimeStep = int(a)
        if o == "-S":
            _firstTimeStep = int(a)
        if o == "-R":
            runId = a
        if o == "-M":
            theModel = a
        if o == "-U":
            saveResults = True
        if o == "-h":
            usage()

    if configfile == None:
        configfile = theModel + ".ini"

    mimo = wfmodel_fit_API(
        _firstTimeStep,
        _lastTimeStep,
        caseName,
        runId,
        modeltofit=theModel,
        config=configfile,
    )
    pars = mimo.pars
    diag = mimo.pars

    catchment = 0
    mimo.saveinitpars()
    for catch in mimo.ColSimS:
        fitname = str(catch) + "_wflow_fit.res"
        mimo.NS = []
        mimo.BIAS = []
        mimo.CORR = []
        mimo.MABSE = []
        # set the catchment
        # print "........> Catchment: " _ str(catchment)
        mimo.ColSim = mimo.ColSimS[catchment]
        mimo.ColMeas = mimo.ColMeasS[catchment]
        mimo.AreaCode = mimo.AreaCodeS[catchment]
        # print mimo.AreaCode
        # print mimo.ColMeas
        # print mimo.ColSim

        qmeas, header = pcrut.readtss(os.path.join(caseName, mimo.qmeasname))
        qmeas = qmeas.astype(np.float64)
        qmeas = qmeas[
            _firstTimeStep - 1 + mimo.WarmUpSteps : _lastTimeStep - 1, mimo.ColMeas
        ]
        lstr = (
            "Currently fitting... Sim: "
            + str(mimo.ColSim)
            + " Meas: "
            + str(mimo.ColMeas)
            + " Area: "
            + str(mimo.AreaCode)
        )
        mimo.log.log(45, lstr)
        pylab.plot(qmeas, "+")
        pylab.title(
            "Sim: "
            + str(mimo.ColSim)
            + " Meas: "
            + str(mimo.ColMeas)
            + " Area: "
            + str(mimo.AreaCode)
        )

        # bruteforce(mimo)
        res = bruteforce(mimo, qmeas, caseName, runId)
        # print pp
        # pylab.plot(mimo.run(pp),"r",linewidth=2.0)
        # printresults(pp,a,b,c,d,mimo.calibpars,os.path.join(caseName,runId,fitname),mimo)
        catchment = catchment + 1
        pylab.clf()
        # mimo.results.append([catchment,pp,a,b,c,d,mimo.NS,mimo.BIAS,mimo.CORR,mimo.MABSE])
        # mimo.savemaps(pp,saveResults)

    mimo.shutdown()

    f = open(os.path.join(caseName, runId, "wflow_fit.csv"), "wb")
    writer = csv.writer(f)
    writer.writerows(res)
    f.close()
    # print pp


if __name__ == "__main__":
    main()
