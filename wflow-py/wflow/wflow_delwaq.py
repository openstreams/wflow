#!/usr/bin/python
"""

Introduction
------------

Simple export library for pcraster/python delwaq link. The module can be
used to export an ldd to a delwaq pointer file and fill the input arrays.
The library includes a command-line interface that allows you to setup
a delwaq model and feed it with forcing data.

.. warning::
    
    This is an experimental version. A complete redesign is needed as this 
    version is unstable and very poorly structured!



the wflow run should have saved at least the folowing mapstacks::

    - self.OldKinWaveVolume=vol
    - self.WaterLevel=lev
    - self.SurfaceRunoff=run
    - self.Inwater=inw  (or the different components that make up this flux)

The script always sets-up at least two Substances, Initial and Check. Initial
is present everywhere at startup and the concentration is zero
in all inputs. Check is not present at startup and set to 1 in all inputs.

The script takes an areamap that can be used to tag water as it enters the
model. This can be a landuse map, a subcatchment map etc. Furthermore water
can also be tagged based on the flux into the model.

The naming of the sustances is as follows: "Area" areamap_class inflow_number 
    
Command line options::
    
    -C: caseDir - set the wflow case directory to use
    -R: runId - set the wflow runId to use
    -T: Set last timestep
    -O: set starttime ('%Y-%m-%d% H:%M:%S')
    -a: Also write dynamic area data if this option is set
    -j: if this option is set the static data is not generated (default is on) 
    -A: sets the areamap used to specify the fraction sources. This can be
        a subcatcment map, a soil type map, a land use maps etc. Default is:
        staticmaps/wflow_subcatch.map (relative to the caseDir directory)
    -D: delwaqdir - set the basedir to create the delwaq schematisation in 
    -S: sourcemap - name of the wflow output map to use as source. 
        it should be a variable that flows into the kinematic wave routine
        inw is normally used as it contain all water per cell that flows into
        the kinematic wave function.
        Use multiple -S options to include multiple maps
    -s: Set the model timesteps in seconds (default 86400)
    -F: if set the model is expected to be run by FEWS. It will determine
        the timesteps from the runinfo.xml file and save the output initial
        conditions to an alternate location. The runinfo.xml file should be located
        in the inmaps directory of the case. 
    -c: Name of the wflow configuration file

  
.. todo::

    add support for a coarser delwaq network based on supplied map.      
    
    
.. todo::
    
    Test option to seperate construction of network from filling of the input
    arrays
    
.. todo::
    
    Ad support to not only follow the kinematic wave reservoir but also
    the flow trough the soil reservoirs. Basically make three layers:
        
        #. kinematic wave reservoir (surface water)
        #. unsaturated store (only vertical flow)
        #. saturated store (horizontal and vertical flow)
    

$Author: schelle $
$Id: wflow_delwaq.py 813 2013-10-07 09:18:36Z schelle $
$Rev: 813 $        
"""
import  wflow.wflow_adapt as wflow_adapt
from  wflow.wf_DynamicFramework import *


 
from datetime import *
import os
import os.path
import shutil, glob
import getopt


logger = ""
volumeMapStack="vol"
runoffMapStack="run"
waterlevelMapStack="lev"

def dw_WriteNrSegments(fname,nr):
    """ 
    Writes the number of segments to B3 file 
    
    B3\_nrofseg.inc
    """
    exfile = open(fname,'w')
    print >>exfile,";Written by dw_WriteNrSegments"
    print >>exfile,str(nr) + " ; nr of segments"
    exfile.close()


def dw_WriteNrExChnages(fname,nr):
    """ 
    Writes the number of exchnages to file (number of rows in the pointer file)
    
    B4\_nrofexch.inc
    """
    exfile = open(fname,'w')
    print >>exfile,";Written by dw_WriteNrExChnages"
    print >>exfile,str(nr) + " 0 0 ; x, y, z direction"
    exfile.close()


def dw_WriteBoundData(fname,areas):
    """ 
    writes B5\_bounddata.inc
    """  
    
    areas = sorted(areas,reverse=True)
    exfile = open(fname,'w')
    print >>exfile,";Written by dw_WriteBoundData"
    for i in areas:
        print >>exfile, "ITEM \'Area_%s\'" %  (i)
        print >>exfile, "CONCENTRATION  \'Area_%s\' \'Check\' \'Initial\'" %  (i)
        print >>exfile, "DATA"
        print >>exfile, "1.0  1.0  0.0"
        print >>exfile, ""
 
    exfile.close()

    
def dw_WriteInitials(fname,inmaps):
    """
    B8_initials.inc
    """
    
    maps = ['Initial','Check']
    exfile = open(fname,'w')
    print >>exfile,"INITIALS"
    for rr in inmaps:
        print >>exfile,"'" + rr + "'", 
    for rr in maps:
        print >>exfile,"'" + rr + "'",        
    print >>exfile
    print >>exfile,"DEFAULTS"
    for rr in inmaps:
        print >>exfile,str(0.0) + " ",
    for rr in maps:
        print >>exfile,str(1.0) + " ",
    print >>exfile    
    exfile.close()
    
    
def dw_WriteBoundlist(fname,pointer,areas,of,inflowtypes):
    """ 
    Writes the boundary list file
    B5\_boundlist.inc
    Numbering is abs(exchnage id)
    
    Input: 
        - fname, pointer
        
    .. todo::
        
        - add labeling of different inflows ( the information is already present)
    """
    totareas = areas
    exfile = open(fname,'w')
    print >>exfile,";Written by dw_WriteBoundlist"
    print >>exfile,";'NodeID' 'Number' 'Type'"
    nr_inflowtypes = len(inflowtypes)
    
    #for i in range(nr_inflowtypes-1):
    #    totareas = vstack((totareas,areas))
    totareas = areas

    arid = 0
    for i in range(len(pointer)):
        if pointer[i,1] < 0:
            print >>exfile,"'BD_" + str(absolute(pointer[i,1])) +  "'  '" + str(absolute(pointer[i,1])) + "'" + " 'Outflow'"
        elif   pointer[i,0] < 0:
            #ar = int(absolute(totareas[arid]))
            ar = totareas[arid]
            print >>exfile,"'BD_" +str(absolute(pointer[i,0])) + "' " + "'" + str(absolute(pointer[i,0])) + "'" + " 'Area_" + str(ar) + "'"
            arid = arid + 1
                
    exfile.close()    
 

def dw_WritePointer(fname,pointer):
    """ 
    WRites the pointer file
    B4\_pointer.inc
    """
    exfile = open(fname,'w')
    print >>exfile,";Written by dw_WritePointer"
    print >>exfile,";nr of pointers is: ", str(pointer.shape[0])
    savetxt(exfile,pointer,fmt='%10.0f')
    exfile.close()   
            

def dw_WriteSegmentOrExchangeData(ttime,fname,datablock,boundids,WriteAscii=True):
    """
    Writes a timestep to a segment/exchange data file (appends to an existing
    file or creates a new one). 
        
    Input:
        - time - time for this timestep  
        - fname - File path of the segment/exchange data file</param>
        - datablock - array with data
        - boundids to write more than 1 block
        - WriteAscii - set to 1 to alse make an ascii dump
        
    """
    # First convert the array to a 32 bit float
    totareas = datablock
    for i in range(boundids-1):
        totareas = vstack((totareas,datablock))
    
    artow= array(totareas,dtype=float32).copy()
    timear= array(ttime,dtype=int32)
    if os.path.isfile(fname): # append to existing file
        fp = open(fname,'ab')
        tstr = timear.tostring() + artow.tostring()
        fp.write(tstr) 
        if WriteAscii:
            fpa = open(fname+".asc",'a')
            timear.tofile(fpa,format="%d\t", sep=":")
            artow.tofile(fpa,format="%10.8f", sep="\t")
            fpa.write('\n')
    else:
        fp = open(fname,'wb')
        tstr = timear.tostring() + artow.tostring()
        fp.write(tstr)
        if WriteAscii:
            fpa = open(fname+".asc",'w')
            timear.tofile(fpa,format="%d\t", sep=":")
            artow.tofile(fpa,format="%10.8f", sep="\t") 
            fpa.write('\n')
        
    fp.close()
    if WriteAscii:
        fpa.close()


#TODO: Add exta column with boundary labels (of the inflows)

def dw_mkDelwaqPointers(ldd,amap,difboun,layers):
    """
    An ldd is used to determine the from-to relations for delwaq using
    the PCraster up/downstreams commands.
    *amap* is used to link boundaries to the segments for delwaq (negative 
    numbers). These are area based boundaries. Diffboun is a 
    python dictionary with inflows for each
    cell.
    
    Input:
        - ldd
        - map to determine the active points)
        - difboun : number of inflow boundaries per cell
        - layers [nr of soil layers (only vertical flow)]. 
    
    .. note:: Only one layer at present (layers must be 1)
        
    Output:
        - pointer, fromto, outflows, boundaries, segment
        - matrix with 4 colums: from to, zero, zero.
        - catchid

    .. note::  use savetxt("pointer.inc",pointer,fmt='%10.0f') to save this
        for use with delwaq
       
       
    .. note:: The pointers list first contains the "internal" fluxes in
        the kinematic wave reservoir, next the fluxes (1-n) into the 
        kinematic wave reservoir.
            
        
    .. todo:: 
        Add exta column with boundary labels (of the inflows)    
        
    """
    # Firts make sure there is at least on outflow in the model
    ptid = uniqueid(boolean(amap))
    flowto = downstream(ldd,ptid)
    # Fix if downsteam is no pit.In that case flowto is missing, set it so itself
    hasflowto = defined(flowto)
    flowto = ifthenelse(defined(ptid) != hasflowto, ptid, flowto)

    # find all upstream cells (these must be set negative)
    upbound = upstream(ldd,1.0)
    upbound = ifthen(amap > 0, upbound)
    # Find the lower boundaries (and pits). These flow to themselves
    
    
    # make into flatted numpy arrays
    np_ptid = pcr2numpy(ptid,NaN).flatten()
    np_flowto = pcr2numpy(flowto,NaN).flatten()
    np_catchid = pcr2numpy(scalar(amap),NaN).flatten()
    np_upbound = pcr2numpy(upbound,NaN).flatten()

    # remove all non-active cells
    np_catchid = np_catchid[np_catchid > 0.0]
    np_upbound = np_upbound[isfinite(np_upbound)]
    np_flowto = np_flowto[isfinite(np_flowto)]
    np_ptid = np_ptid[isfinite(np_ptid)]
    np_flowto= np_flowto.reshape(len(np_flowto),1)
    np_ptid= np_ptid.reshape(len(np_ptid),1)
    np_catchid= np_catchid.reshape(len(np_catchid),1) * -1
    # Now make catchid a list
    np_catchid = np_catchid.flatten()
    np_catchid = array(int_(np_catchid),dtype='|S').tolist()
    # find all downstream segments (flowto == ptid)
    # now set the flowto points (outflows, usually just one) also to negative
    lowerck = absolute(np_ptid) == absolute(np_flowto)
    # mak epointer matrix and add to zero zolumns
    orgpointer = hstack((np_ptid,np_flowto,zeros((len(np_flowto),1)),zeros((len(np_flowto),1))))
    pointer = orgpointer.copy()
    extraboun = []
    # Add the inflow boundaries here.
    cells = pointer[:,0]
    cells = cells.reshape(len(cells),1)
    bounid = cells.copy()
    zzerocol = zeros((len(np_flowto),1))
    inflowId = bounid.copy()
    
    
    # outflow to pointer
    # point -> - point
    lopt = np_ptid[lowerck]
    lopt = lopt.reshape(len(lopt),1)
    zerocol = zeros((len(lopt),1))
    lowerids = arange(1,len(lopt) + 1) * -1
    #of = hstack((lopt,lopt * -1.0,zerocol,zerocol))
    lowerids = lowerids.reshape(len(lowerids),1)
    of = hstack((lopt,lowerids,zerocol,zerocol))
    
    #pointer = vstack((pointer,of))
    
    # Now remove double pointer to itself and replace by lower boundary
    lowerck = pointer[:,0] == pointer[:,1]   
    pointer[lowerck,:] = of 
    start = absolute(lowerids.min()) + 1
    bouns = 1
    for idd in range(1,difboun + 1):
        inflowId[:] = idd
        bounid = arange(start,(len(cells)+start)).reshape((len(cells),1)) * -1.0
        if bouns == 1:
            extraboun = hstack((bounid,cells,zzerocol,zzerocol))        
        else:
            extraboun = vstack((extraboun,hstack((bounid,cells,zzerocol,zzerocol))))
        bouns = bouns +1
        start = start + len(cells)
    
    res = []
    for idd in range(1,difboun + 1):
        ct = list(np_catchid)
        print "ct: "
        print unique(ct)
        for i in range(0,len(np_catchid)):
            ct[i] = np_catchid[i] + "_" + str(idd)
        res = res + ct
    print unique(res)
    np_catchid = res
    #pointer = vstack((pointer,extraboun))
    # now catchment id's
    #zerocol = zeros((len(np_catchid),1))
    #extraboun= hstack((np_catchid,cells,zerocol,zerocol))
    #print np_catchid
    
    
    if len(extraboun) > 0:
        pointer = vstack((pointer,extraboun)) 
      
    return ptid, flowto, pointer, orgpointer[:,0], of[:,0:2], extraboun[:,0:1].flatten(), np_ptid.flatten(), np_catchid



def dw_pcrToDataBlock(pcrmap):
    """
    Converts a pcrmap to a numpy array.that is flattend and from which
    missing values are removed. Used for generating delwaq data
    """    
    ttar = pcr2numpy(pcrmap,NaN).flatten()
    ttar = ttar[isfinite(ttar)]
    
    return ttar
            

 
def readTS(name, ts):
    """
    Read a pcraster map for a timestep without using the dynamic framework
    
    """
    mname = os.path.basename(name)
    #  now generate timestep
    tsje = "%0.11d" % ts
    ff = mname + tsje[len(mname):]
    ff = ff[:8] + "." + ff[8:]
    name = os.path.dirname(name) + "/" + ff
    mapje = readmap(name)
    
    return mapje

def dw_CreateDwRun(thedir):
    """"
    create the dir to save delwaq info in
    """
    if not os.path.isdir(thedir):
        os.makedirs(thedir + "/fixed/")
        os.makedirs(thedir + "/includes_deltashell/")
        os.makedirs(thedir + "/includes_flow/")
        os.makedirs(thedir + "/debug/")
    if os.path.exists(thedir + "/includes_flow/area.dat"):
        os.remove(thedir + "/includes_flow/area.dat")
    if os.path.exists(thedir + "/includes_flow/flow.dat"):
        os.remove(thedir + "/includes_flow/flow.dat")
    if os.path.exists(thedir + "/includes_flow/volume.dat"):
        os.remove(thedir + "/includes_flow/volume.dat")
    if os.path.exists(thedir + "/includes_flow/length.dat"):
        os.remove(thedir + "/includes_flow/length.dat")
    if os.path.exists(thedir + "/includes_flow/surface.dat"):
        os.remove(thedir + "/includes_flow/surface.dat")        
    if os.path.exists(thedir + "/includes_flow/area.dat.asc"):
        os.remove(thedir + "/includes_flow/area.dat.asc")
    if os.path.exists(thedir + "/includes_flow/flow.dat.asc"):
        os.remove(thedir + "/includes_flow/flow.dat.asc")
    if os.path.exists(thedir + "/includes_flow/volume.dat.asc"):
        os.remove(thedir + "/includes_flow/volume.dat.asc")
    if os.path.exists(thedir + "/includes_flow/length.dat.asc"):
        os.remove(thedir + "/includes_flow/length.dat.asc")
    if os.path.exists(thedir + "/includes_flow/surface.dat.asc"):
        os.remove(thedir + "/includes_flow/surface.dat.asc")          


def dw_Write_Times(dwdir,T0,timeSteps,timeStepSec):
    """
    Writes B1_T0.inc, B2_outputtimers.inc, B2_sysclock.inc and /B2_simtimers.inc
    Assumes daily timesteps for now!
    """
    # B1_T0.inc
    exfile = open(dwdir + "/B1_T0.inc",'w')
    print >>exfile, "\'T0: " + T0.strftime("%Y.%m.%d %H:%M:%S") + "  (scu=       1s)\'"
    exfile.close()

    # B2_outputtimers.inc
    timeRange  = timedelta(seconds=timeStepSec * timeSteps)
    
    days = int(timeStepSec / 86400) 
    hours = int(timeStepSec / 3600) 
    minutes = int(timeStepSec / 60) 
    seconds = int(timeStepSec - minutes*60)
    minutes -= hours*60
    hours -= days*24   
    timestepstring = "  %03d%02d%02d%02d" % (days, hours, minutes, seconds)
    
    exfile = open(dwdir + "/B2_outputtimers.inc",'w')
    etime = T0 + timeRange
    print >>exfile, "  " + T0.strftime("%Y/%m/%d-%H:%M:%S") + "  " + etime.strftime("%Y/%m/%d-%H:%M:%S") + timestepstring
    print >>exfile, "  " + T0.strftime("%Y/%m/%d-%H:%M:%S")  + "  " + etime.strftime("%Y/%m/%d-%H:%M:%S") + timestepstring
    print >>exfile, "  " + T0.strftime("%Y/%m/%d-%H:%M:%S")  + "  " + etime.strftime("%Y/%m/%d-%H:%M:%S") + timestepstring 
    exfile.close()
    
    #B2_simtimers.inc
    exfile = open(dwdir + "/B2_simtimers.inc",'w')
    print >>exfile, "  " + T0.strftime("%Y/%m/%d-%H:%M:%S")
    print >>exfile, "  " + etime.strftime("%Y/%m/%d-%H:%M:%S")
    print >>exfile, "  0 ; timestep constant"
    print >>exfile, "; dddhhmmss format for timestep"
    print >>exfile, timestepstring + " ; timestep"
    exfile.close()
    
    #B2_sysclock.inc
    exfile = open(dwdir + "/B2_sysclock.inc",'w')
    print >>exfile,"%7d \'DDHHMMSS\' \'DDHHMMSS\'  ; system clock" % timeStepSec
    exfile.close()


def dw_Write_Substances(fname,areas):
    """
    Writes the B1_sublist.inc file
    input:
        
        it writes substances for the areas and an initial and mass balance 
        check substance
        
    """

    exfile = open(fname,'w')
    areas = sorted(areas,reverse=True)
    print >>exfile,"; number of active and inactive substances"
    print >>exfile,"%d         0" % (len(areas) + 2)
    print >>exfile,"; active substances"
    print >>exfile, "1             \'Initial\' ; "
    print >>exfile, "2             'Check' ; "
    j = 2
    for i in areas:
        j = j + 1
        print >>exfile, "%d            \'Area_%s\'" %  (j,i)
    print >>exfile,"; passive substances"
    
        
    exfile.close()
    
    
def dw_Write_B2_outlocs(fname,gauges,segs):
    """
    Write an output loc file based on the wflow_gauges
    map.
    """
    segs = ifthenelse(gauges > 0,segs,NaN)
    gauges = ifthenelse(gauges > 0,scalar(gauges),NaN)
    np_gauges = pcr2numpy(gauges,NaN).flatten()
    np_segs = pcr2numpy(segs,NaN).flatten()
        
    np_gauges = np_gauges[isfinite(np_gauges)]
    np_segs = np_segs[isfinite(np_segs)]
    
    if len(np_segs) != len(np_gauges):
        logger.error("Gauges and segments do not match!")

    pts = size(np_segs)
    exfile = open(fname,'w')
    print >>exfile,"%d ; nr of locations" % pts
    print >>exfile,"; \'outlocname\' numberofsegments segment list"
    i = 0
    for loc in np_gauges:
        print >>exfile," \'%d\' 1 %d" % (loc, np_segs[i])
        i = i + 1
    exfile.close()


   

def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)
                
pointer = ""




def main():
    caseId = "Ahr_DW/"
    caseId = "default_hbv"
    runId = "run_default"
    dwdir = "dw_rhine"
    areamap = "staticmaps/wflow_subcatch.map"
    timeSteps = 1
    timestepsecs = 86400
    configfile = "wflow_sbm.ini"
    sourcesMap = []
    fewsrun = False
    WriteAscii=False
    Write_Dynamic= False
    Write_Structure = True
    T0 = datetime.strptime("2000-01-01 00:00:00",'%Y-%m-%d %H:%M:%S')

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'adD:C:R:S:hT:F:s:O:A:jc:')
    except getopt.error, msg:
        pcrut.usage(msg)
    
    for o, a in opts:
        if o == '-F': 
            runinfoFile = a
            fewsrun = True
        if o == '-C':caseId = a
        if o == '-R': runId = a
        if o == '-D': dwdir = a
        if o == '-d': Write_Dynamic = True
        if o == '-f': Write_Structure = False
        if o == '-s': timestepsecs = int(a)
        if o == '-S': sourcesMap.append(a)
        if o == '-h': usage()
        if o == '-T': timeSteps = int(a)
        if o == '-A': areamap = a.strip()
        if o == '-c': configfile = a.strip()
        if o == '-O': T0 = datetime.strptime(a,'%Y-%m-%d %H:%M:%S')
        
    

   
    global pointer
    dw_CreateDwRun(dwdir)


    config = ConfigParser.SafeConfigParser()
    config.optionxform = str
    config.read(caseId + "/" + configfile)
    
    timestepsecs = int(configget(config,"model","timestepsecs",str(timestepsecs)))
    
    if fewsrun: 
        timeSteps =  wflow_adapt.getTimeStepsfromRuninfo(runinfoFile,timestepsecs)
        T0 =  wflow_adapt.getStartTimefromRuninfo(runinfoFile)
        print timeSteps
    
    #: we need one delwaq calculation timesteps less than hydrology
    # timeSteps = timeSteps # need one more hydrological timestep as dw timestep
    firstTimeStep = 0
        
        
    logger = pcrut.setlogger(dwdir + "/debug/wflow_delwaq.log","wflow_delwaq") 
    #caseid = "default_hbv"
    logger.info("T0 of run: " + str(T0))
    boundids = len(sourcesMap)  # extra number of exchnages for all bounds

    #Number of exchnages is elements minus number of outflows!!  
    
    # Get subcatchment data
    logger.info("Reading basemaps")
    
    setclone(caseId + configget(config,"model","wflow_subcatch","//staticmaps/wflow_subcatch.map"))
    print caseId + "/staticmaps/wflow_subcatch.map"
    amap = scalar(readmap(caseId +  "/" + areamap))
    #amap = scalar(readmap(caseId + "/staticmaps/wflow_subcatch.map"))
    modelmap = readmap(caseId + configget(config,"model","wflow_subcatch","/staticmaps/wflow_subcatch.map"))
        # get ldd
    ldd = readmap(caseId + configget(config,"model","wflow_ldd","/staticmaps/wflow_ldd.map"))
    gauges = readmap(caseId + configget(config,"model","wflow_gauges","/staticmaps/wflow_gauges.map"))

    # Some models yield a reallength.map, others a rl.map.
    rl_map_file = caseId + "/" + runId + "/outsum/rl.map"
    if not os.path.exists(rl_map_file):
        rl_map_file = caseId + "/" + runId + "/outsum/reallength.map"
    cellsize = float(pcr2numpy(readmap(rl_map_file),NaN)[0,0])
    logger.info("Cellsize model: " + str(cellsize))
    
    # Limit areas map to modelmap (subcatchments)
    amap = ifthen(modelmap > 0, amap)
    ldd = ifthen(amap > 0, ldd)
    report(amap,dwdir +"/debug/area.map")
    report(ldd,dwdir +"/debug/ldd.map")
    report(modelmap,dwdir +"/debug/modelmap.map")
    
    thecells = pcr2numpy(modelmap,NaN).flatten()
    nrcells = len(thecells)
    nractcells = len(thecells[isfinite(thecells)])
    
    logger.info("Total number gridcells (including inactive): " + str(nrcells))        
    logger.info("Total number of used gridcells: " + str(nractcells))
    
        # find all upstream cells (these must be set negative)
    upbound = upstream(ldd,1.0)
    upbound = ifthen(upbound == 0, upbound)
    upar=pcr2numpy(scalar(upbound),NaN).flatten()
    logger.info("Number of upstream cells (without upstream connection): " + str(len(upar[isfinite(upar)])))
    report(upbound,dwdir +"/debug/upbound.map")
    
   

    if Write_Structure:    
        # get pointer an boundaries from ldd, subcatch and defined boundaries (P only now)
        ptid, flowto, pointer, fromto, of , bounds, segments, areas = dw_mkDelwaqPointers(ldd,amap,boundids,1)

        save(dwdir +"/debug/pointer.npy",pointer)
        save(dwdir +"/debug/fromto.npy",fromto)
        save(dwdir +"/debug/of.npy",of)
        save(dwdir +"/debug/bounds.npy",bounds)
        save(dwdir +"/debug/segments.npy",segments)
        save(dwdir +"/debug/areas.npy",areas)
    
        # Write id maps to debug area
        report(ptid,dwdir + "/debug/ptid.map")
        report(flowto,dwdir + "/debug/flowto.map")
        logger.info("Unique areas: " + str(unique(areas)))
        #logger.info("Number of area inflows: " + str(len(areas) * boundids))
        logger.info("Number of segments: " + str(len(segments.flatten())))
        logger.info("Number of internal flows: " + str(len(fromto.flatten())))
        logger.info("outflow  ids: " + str(of))
        logger.info("source maps: " + str(sourcesMap))    
        
        NOOF = of.shape[0]
        NOSQ = segments.shape[0]
        NOQ = pointer.shape[0]
        
        dw_WriteNrSegments(dwdir + "/includes_deltashell/B3_nrofseg.inc",NOSQ)
        # Write pointer file
        #TODO: add sources maps here (no only one source supported)
        dw_WritePointer(dwdir + "/includes_deltashell/B4_pointer.inc",pointer)
        # Write the number of exchanges
        dw_WriteNrExChnages(dwdir + "/includes_deltashell/B4_nrofexch.inc",NOQ)
        dw_WriteBoundlist(dwdir + "/includes_deltashell/B5_boundlist.inc",pointer,areas,of,sourcesMap)
        dw_WriteBoundData(dwdir + "/includes_deltashell/B5_bounddata.inc",unique(areas))
    
        dw_WriteInitials(dwdir + "/includes_deltashell/B8_initials.inc",sourcesMap)
        dw_Write_Substances(dwdir + "/includes_deltashell/B1_sublist.inc",unique(areas))    
        dw_Write_B2_outlocs(dwdir + "/includes_deltashell/B2_outlocs.inc",gauges,ptid)
        

    

    internalflowwidth = readmap(caseId + "/" + runId + "/outsum/Bw.map")
    internalflowlength = readmap(caseId + "/" + runId + "/outsum/DCL.map")
    surface_map = internalflowwidth * internalflowlength
    surface_block = dw_pcrToDataBlock(surface_map)
    logger.info("Writing surface.dat. Nr of points: " + str(size(surface_block)))
    dw_WriteSegmentOrExchangeData(0,dwdir + '/includes_flow/surface.dat',surface_block,1,WriteAscii)
            
    
    # create dummy length file
    length_block = zeros(pointer.shape[0] * 2) + 0.5
    # write  length file
    logger.info("Writing length.dat. Nr of points: " + str(size(length_block)))
    dw_WriteSegmentOrExchangeData(0,dwdir + '/includes_flow/length.dat',length_block,1,WriteAscii)
    
    # mask to filter out inactive segments
    zero_map = 0.0*scalar(ptid)

    
    ts = 1
    
    if Write_Dynamic:       
        dw_Write_Times(dwdir + "/includes_deltashell/",T0,timeSteps-1,timestepsecs)
        
        for i in range(firstTimeStep,timeSteps * timestepsecs,timestepsecs):
            volume_map = readTS(caseId + "/" + runId + "/outmaps/" + volumeMapStack,ts)
            volume_block = dw_pcrToDataBlock(volume_map)
            
            # volume for each timestep and number of segments
           
            logger.info("Writing volumes.dat. Nr of points: " + str(size(volume_block)))
            dw_WriteSegmentOrExchangeData(i,dwdir + '/includes_flow/volume.dat',volume_block,1,WriteAscii)

            
            # Now write arreas  
            
            
        
            # Now write the flows (exchnages)
            # First read the flows in the kinematic wave reservoir (internal exchnages)
            flow = readTS(caseId + "/" + runId + "/outmaps/" + runoffMapStack,ts)
            flow_block_Q = dw_pcrToDataBlock(flow)
            # now the inw
            flowblock = flow_block_Q
            
            
            wlevel = readTS(caseId + "/" + runId + "/outmaps/" + waterlevelMapStack,ts)
            areadyn = wlevel * internalflowwidth
            area_block_Q = dw_pcrToDataBlock(areadyn)
            area_block = area_block_Q
            
            # Now read the inflows in each segment (water that enters the kinamatic 
            # wave reservoir). Also write the areas
            for source in sourcesMap:
                logger.info("Step: " + str(ts) + " source: " + str(source))
                thesource = readTS(caseId + "/" + runId + "/outmaps/" + source ,ts)
                thesource = zero_map + thesource
                flow_block_IN = dw_pcrToDataBlock(thesource)
                flowblock = hstack((flowblock,flow_block_IN))
                area_block = hstack((area_block,surface_block))
                
            logger.info("Writing flow.dat. Nr of points: " + str(size(flowblock)))
            dw_WriteSegmentOrExchangeData(i,dwdir + '/includes_flow/flow.dat',flowblock,1,WriteAscii)
            logger.info("Writing area.dat. Nr of points: " + str(size(area_block)))
            dw_WriteSegmentOrExchangeData(i,dwdir + '/includes_flow/area.dat',area_block,1,WriteAscii)            
            
            ts = ts + 1
        
        """
        Write last volume block with current kinwavevol
        """  
        ts = ts -1    
        i = i + 1
        logger.info("Writing last step..")
        
        
        logger.info("Writing area.dat. Nr of points: " + str(size(area_block)))
        dw_WriteSegmentOrExchangeData(i,dwdir + '/includes_flow/area.dat',area_block,1,WriteAscii)
       
        #logger.info("Writing surface.dat. Nr of points: " + str(size(surface_block)))
        #dw_WriteSegmentOrExchangeData(i,dwdir + '/includes_flow/surface.dat',surface_block,1,WriteAscii)
        
        logger.info("Writing flow.dat. Nr of points: " + str(size(flowblock)))
        dw_WriteSegmentOrExchangeData(i,dwdir + '/includes_flow/flow.dat',flowblock,1,WriteAscii)
 
        volume_map = readTS(caseId + "/" + runId + "/outmaps/voln",ts)
        volume_block = dw_pcrToDataBlock(volume_map)
        logger.info("Writing volumes.dat. Nr of points: " + str(size(volume_block)))
        dw_WriteSegmentOrExchangeData(i,dwdir + '/includes_flow/volume.dat',volume_block,1,WriteAscii)
            
    
    
if __name__ == "__main__":
    main()
