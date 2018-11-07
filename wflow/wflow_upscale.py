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
wflow_upscale -- resample a working wflow model to a lower resolution version

Usage::

    -C CaseName
    -N NewCaseName
    -r resample factor
    -I skip input mapstacks if specified
    -f force overwrite an existing model
    -M maxcpu
       maximum number of cpu's/cores to use (default = 4)

The script uses the pcraster resample program to reduce the maps. The original
river network is used to force the river network in the reduced version of the
model. Nevertheless it may be needed to manually adjust the locations of
the gauges in the gauges.col file.

A more sophisticated method for resampling is implemented in the wflow_prepare
scripts. In general, if you need a high and a low resolution model it is best to
run the wflow_prepare scripts twice to create the different models.
"""


from wflow.wflow_lib import *


import os, sys, shlex, time
import os.path
import glob
import getopt
import subprocess


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args:
        print(msg)
    print(__doc__)
    sys.exit(0)


def runCommands(commands, maxCpu):
    """
    Runs a list of processes dividing
    over maxCpu number of cores.
    """

    def removeFinishedProcesses(processes):
        """ given a list of (commandString, process),
            remove those that have completed and return the result
        """
        newProcs = []
        for pollCmd, pollProc in processes:
            retCode = pollProc.poll()
            if retCode == None:
                # still running
                newProcs.append((pollCmd, pollProc))
            elif retCode != 0:
                # failed
                raise Exception("Command %s failed" % pollCmd)
            else:
                print("Command %s completed successfully" % pollCmd)
        return newProcs

    processes = []
    for command in commands:
        command = command.replace(
            "\\", "/"
        )  # otherwise shlex.split removes all path separators
        proc = subprocess.Popen(shlex.split(command))
        procTuple = (command, proc)
        processes.append(procTuple)
        while len(processes) >= maxCpu:
            time.sleep(0.2)
            processes = removeFinishedProcesses(processes)

    # wait for all processes
    while len(processes) > 0:
        time.sleep(0.5)
        processes = removeFinishedProcesses(processes)
    print("All processes in que (" + str(len(commands)) + ") completed.")


def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:], "fhC:N:Ir:M:")
    except getopt.error as msg:
        usage(msg)

    factor = 1
    Verbose = 1
    inmaps = True
    force = False
    caseName = "rhineNew"
    caseNameNew = "rhineNew_resamp"
    maxcpu = 4

    for o, a in opts:
        if o == "-C":
            caseName = a
        if o == "-N":
            caseNameNew = a
        if o == "-r":
            factor = int(a)
        if o == "-I":
            inmaps = False
        if o == "-h":
            usage()
        if o == "-f":
            force = True
        if o == "-M":
            maxcpu = int(a)

    dirs = ["/intbl/", "/inmaps/", "/staticmaps/", "/intss/", "/instate/", "/outstate/"]
    ext_to_copy = ["*.tss", "*.tbl", "*.col", "*.xml"]
    if os.path.isdir(caseNameNew) and not force:
        print("Refusing to write into an existing directory:" + caseNameNew)
        sys.exit()

    if not os.path.isdir(caseNameNew):
        for ddir in dirs:
            os.makedirs(caseNameNew + ddir)
        for inifile in glob.glob(caseName + "/*.ini"):
            shutil.copy(inifile, inifile.replace(caseName, caseNameNew))

    for ddir in dirs:
        allcmd = []
        for mfile in glob.glob(caseName + ddir + "/*.map"):
            if "_ldd.map" not in mfile:
                mstr = (
                    "resample -r "
                    + str(factor)
                    + " "
                    + mfile
                    + " "
                    + mfile.replace(caseName, caseNameNew)
                )
                # print mstr
                allcmd.append(mstr)
            # os.system(mstr)
        runCommands(allcmd, maxcpu)

        if inmaps:
            allcmd = []
            for mfile in glob.glob(caseName + ddir + "/*.[0-9][0-9][0-9]"):
                mstr = (
                    "resample -r "
                    + str(factor)
                    + " "
                    + mfile
                    + " "
                    + mfile.replace(caseName, caseNameNew)
                )
                if not os.path.exists(mfile.replace(caseName, caseNameNew)):
                    # print mstr
                    allcmd.append(mstr)
                    # os.system(mstr)
                else:
                    print("skipping " + mfile.replace(caseName, caseNameNew))
            runCommands(allcmd, maxcpu)

        for ext in ext_to_copy:
            for mfile in glob.glob(caseName + ddir + ext):
                shutil.copy(mfile, mfile.replace(caseName, caseNameNew))

    # Because the ldd cannot be resampled this way we have to recreate
    # in including the subcatchments that are derived from it

    print("recreating static maps ...")
    # Create new ldd using old river network
    dem = readmap(caseNameNew + "/staticmaps/wflow_dem.map")
    # orig low res river
    riverburn = readmap(caseNameNew + "/staticmaps/wflow_river.map")
    # save it
    report(riverburn, caseNameNew + "/staticmaps/wflow_riverburnin.map")
    demburn = cover(ifthen(boolean(riverburn), dem - 600), dem)
    print("Creating ldd...")
    ldd = lddcreate_save(
        caseNameNew + "/staticmaps/wflow_ldd.map", demburn, True, 10.0e35
    )
    ## Find catchment (overall)
    outlet = find_outlet(ldd)
    sub = subcatch(ldd, outlet)
    report(sub, caseNameNew + "/staticmaps/wflow_catchment.map")
    report(outlet, caseNameNew + "/staticmaps/wflow_outlet.map")
    # os.system("col2map --clone " + caseNameNew + "/staticmaps/wflow_subcatch.map " + caseNameNew + "/staticmaps/gauges.col " + caseNameNew + "/staticmaps/wflow_gauges.map")
    gmap = readmap(caseNameNew + "/staticmaps/wflow_gauges.map")
    scatch = subcatch(ldd, gmap)
    report(scatch, caseNameNew + "/staticmaps/wflow_subcatch.map")


if __name__ == "__main__":
    main()
