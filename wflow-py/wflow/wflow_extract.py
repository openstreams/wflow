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
wflow_extract -- cut out a section of working wflow model and optionally resample to a higher resolution

Usage::

    -C CaseName
    -N NewCaseName
    -c clone map
    -I skip input mapstacks if specified
    -f force overwrite an existing model

The script uses the pcraster resample program to reduce the maps.

"""


from wflow.wflow_lib import *
import wflow.pcrut as pcrut
import sys
import os
import os.path
import glob
import getopt


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args:
        print msg
    print __doc__
    sys.exit(0)


def main(argv=None):

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return

    opts, args = getopt.getopt(argv, "fhC:N:Ir:c:")

    factor = 1
    Verbose = 1
    inmaps = True
    force = False
    caseName = "rhineNew"
    caseNameNew = "rhineNew_resampaa"
    cloneMap = "clone.map"

    for o, a in opts:
        if o == "-C":
            caseName = a
        if o == "-c":
            cloneMap = a
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

    dirs = ["/intbl/", "/inmaps/", "/staticmaps/", "/intss/", "/instate/", "/outstate/"]
    if os.path.isdir(caseNameNew) and not force:
        print "Refusing to write into an existing directory:" + caseNameNew
        sys.exit()

    if not os.path.isdir(caseNameNew):
        for ddir in dirs:
            os.makedirs(caseNameNew + ddir)
        for inifile in glob.glob(caseName + "/*.ini"):
            shutil.copy(inifile, inifile.replace(caseName, caseNameNew))

    for ddir in dirs:
        for mfile in glob.glob(caseName + ddir + "/*.map"):
            mstr = (
                "resample --clone "
                + cloneMap
                + " "
                + mfile
                + " "
                + mfile.replace(caseName, caseNameNew)
            )
            print mstr
            os.system(mstr)
        if inmaps:
            for mfile in glob.glob(caseName + ddir + "/*.[0-9][0-9][0-9]"):
                mstr = (
                    "resample --clone "
                    + cloneMap
                    + " "
                    + mfile
                    + " "
                    + mfile.replace(caseName, caseNameNew)
                )
                if not os.path.exists(mfile.replace(caseName, caseNameNew)):
                    print mstr
                    os.system(mstr)
                else:
                    print "skipping " + mfile.replace(caseName, caseNameNew)
        for mfile in glob.glob(caseName + ddir + "*.tbl"):
            shutil.copy(mfile, mfile.replace(caseName, caseNameNew))
        for mfile in glob.glob(caseName + ddir + "*.col"):
            shutil.copy(mfile, mfile.replace(caseName, caseNameNew))
        for mfile in glob.glob(caseName + ddir + "*.tss"):
            shutil.copy(mfile, mfile.replace(caseName, caseNameNew))

    print "recreating static maps ..."
    # Create new ldd using old river network
    dem = readmap(caseNameNew + "/staticmaps/wflow_dem.map")
    # orig low res river
    riverburn = readmap(caseNameNew + "/staticmaps/wflow_river.map")
    # save it
    report(riverburn, caseNameNew + "/staticmaps/wflow_riverburnin.map")
    demburn = cover(ifthen(boolean(riverburn), dem - 600), dem)
    print "Creating ldd..."
    ldd = lddcreate_save(
        caseNameNew + "/staticmaps/wflow_ldd.map", demburn, True, 10.0E35
    )
    #
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
