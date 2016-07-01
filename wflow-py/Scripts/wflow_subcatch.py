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
wflow_subcatch -- extract a subcatchment from a large model

Usage::

    -C CaseName
    -N NewCaseName
    -s id of subcatchment
    -I skip input mapstacks if specified
    -f force overwrite an existing model
    -M maxcpu
       maximum number of cpu's/cores to use (default = 4)


"""


from wflow.wflow_lib import *
import wflow.pcrut as pcrut

    
import os, sys, shlex, time
import os.path
import glob
import getopt
import subprocess



def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
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
            if retCode==None:
                # still running
                newProcs.append((pollCmd, pollProc))
            elif retCode!=0:
                # failed
                raise Exception("Command %s failed" % pollCmd)
            else:
                print "Command %s completed successfully" % pollCmd
        return newProcs

    processes = []
    for command in commands:
        command = command.replace('\\','/') # otherwise shlex.split removes all path separators
        proc =  subprocess.Popen(shlex.split(command))
        procTuple = (command, proc)
        processes.append(procTuple)
        while len(processes) >= maxCpu:
            time.sleep(.2)
            processes = removeFinishedProcesses(processes)

    # wait for all processes
    while len(processes)>0:
        time.sleep(0.5)
        processes = removeFinishedProcesses(processes)
    print "All processes in que (" + str(len(commands)) + ") completed."


def main():
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'fhC:N:I:s:M:')
    except getopt.error, msg:
        usage(msg)

    factor = 1
    Verbose=1
    inmaps = True
    force = False
    caseName = "thecase"
    caseNameNew = "thecase_resamp"
    maxcpu = 4

    for o, a in opts:
        if o == '-C': caseName = a
        if o == '-N': caseNameNew = a
        if o == '-s': subcatch = int(a)
        if o == '-I': inmaps = False
        if o == '-h': usage()
        if o == '-f': force = True
        if o == '-M': maxcpu = int(a)

    dirs = ['/intbl/',  '/staticmaps/', '/intss/', '/instate/', '/outstate/','/inmaps/' ,'/inmaps/clim/', '/intbl/clim/']
    ext_to_copy = ['*.tss','*.tbl','*.col','*.xml']
    if os.path.isdir(caseNameNew) and not force:
        print "Refusing to write into an existing directory:" + caseNameNew
        exit()

    #ddir = []
    dirs = []
    for (path, thedirs, files) in os.walk(caseName):
        print path
        dirs.append(path)

    if not os.path.isdir(caseNameNew):
        for ddir in dirs:
            os.makedirs(ddir.replace(caseName,caseNameNew))
        for inifile in glob.glob(caseName + "/*.ini"):
            shutil.copy(inifile, inifile.replace(caseName,caseNameNew))


    # read subcatchment map
    x, y, subcatchmap, FillVal = readMap(os.path.join(caseName,'staticmaps','wflow_subcatch.map'), 'PCRaster')
    for ddir in dirs:
        print ddir
        allcmd = []
        for mfile in glob.glob(ddir + '/*.map'):
            x, y, data, FillVal = readMap(mfile,'PCRaster')
            xn, yn, datan = cutMapById(data,subcatchmap,subcatch,x,y,FillVal)
            ofile = mfile.replace(caseName,caseNameNew)
            if data.dtype == np.int32 or  data.dtype == np.uint8:
                writeMap(ofile,'PCRaster',xn,yn,datan.astype(np.int32),FillVal)
            else:
                writeMap(ofile, 'PCRaster', xn, yn, datan, FillVal)

            # Assume ldd and repair
            if data.dtype == np.uint8:
                myldd = ldd(readmap(ofile))
                myldd = lddrepair(myldd)
                report(myldd,ofile)

        for mfile in glob.glob(ddir + '/*.[0-9][0-9][0-9]'):
            x, y, data, FillVal = readMap(mfile,'PCRaster')
            xn, yn, datan = cutMapById(data,subcatchmap,subcatch,x,y,FillVal)
            ofile = mfile.replace(caseName,caseNameNew)
            if data.dtype == np.int32 or  data.dtype == np.uint8:
                writeMap(ofile,'PCRaster',xn,yn,datan.astype(np.int32),FillVal)
            else:
                writeMap(ofile, 'PCRaster', xn, yn, datan, FillVal)


        for ext in ext_to_copy:
            for mfile in glob.glob(os.path.join(ddir, ext)):
                shutil.copy(mfile, mfile.replace(caseName,caseNameNew))

        # Copy ini files
        for mfile in glob.glob(os.path.join(caseName,'*.ini')):
            shutil.copy(mfile, mfile.replace(caseName, caseNameNew))





if __name__ == "__main__":
    main()
