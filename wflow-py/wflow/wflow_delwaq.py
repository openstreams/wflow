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
import time
import struct
import shutil
import __builtin__


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
 

def dw_WritePointer(fname,pointer,binary=False):
    """ 
    WRites the pointer file
    B4\_pointer.inc
    """
    if not binary:
        # Write ASCII file
        exfile = open(fname,'w')
        print >>exfile,";Written by dw_WritePointer"
        print >>exfile,";nr of pointers is: ", str(pointer.shape[0])
        savetxt(exfile,pointer,fmt='%10.0f')
        exfile.close()   
    else:
        # Write binary file
        f = open(fname, 'wb')
        for i in range(pointer.shape[0]):
            f.write(struct.pack('4i',*pointer[i,:]))
        f.close()


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
    # prepare hydfile directory
    comdir = os.sep.join([thedir,'com'])
    if os.path.isdir(comdir):
        shutil.rmtree(comdir)
    os.mkdir(comdir)



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


def dw_GetGridDimensions(ptid_map):
    """
    Returns number of cells in 1st and 2nd grid directions.

    input:
    - ptid_map : PCRaster map with unique id's
    """
    # find number of cells in m and n directions
    zero_map = scalar(ptid_map) * 0.0
    allx = dw_pcrToDataBlock(xcoordinate(boolean(cover(zero_map + 1,1))))
    i = 0
    diff = round(__builtin__.abs(allx[i] - allx[i+1]), 5)
    diff_next = diff
    while diff_next == diff:
        i += 1
        diff_next = __builtin__.abs(allx[i] - allx[i+1])
        diff_next = round(diff_next, 5)
    m = i+1
    n = allx.shape[0] / m
    m, n = n, m
    return m,n


def dw_WriteGridFiles(fname,ptid_map, pointer):
    """
    Writes Delwaq indices (*.lga) and coordinates (*.cco) files.

    input:
    - fname    : output file name (without file extension)
    - ptid_map : PCRaster map with unique id's
    """
    # horizontal dimensions
    m,n = dw_GetGridDimensions(ptid_map)
    # number of layers
    nolay = 1

    # prepare cco data
    zero_map = scalar(ptid_map) * 0.0
    setglobaloption('coorul')
    xxul = dw_pcrToDataBlock(xcoordinate(boolean(cover(zero_map + 1,1))))
    yyul = dw_pcrToDataBlock(ycoordinate(boolean(cover(zero_map + 1,1))))
    setglobaloption('coorlr')
    xxlr = dw_pcrToDataBlock(xcoordinate(boolean(cover(zero_map + 1,1))))
    yylr = dw_pcrToDataBlock(ycoordinate(boolean(cover(zero_map + 1,1))))
    setglobaloption('coorcentre')
    # reshape to 2d array
    xx = xxul.reshape((m,n))
    yy = yyul.reshape((m,n))
    # add extra row and columns for last corner coordinates
    xx = insert(xx,m,0,axis=0)
    xx = insert(xx,n,0,axis=1)
    yy = insert(yy,m,0,axis=0)
    yy = insert(yy,n,0,axis=1)
    # reshape lower right coordinates
    xxlr = xxlr.reshape((m,n))
    yylr = yylr.reshape((m,n))
    # use lower right coordinates were possible
    xx[1:,-1] = xxlr[:,-1]
    xx[-1,1:] = xxlr[-1,:]
    yy[1:,-1] = yylr[:,-1]
    yy[-1,1:] = yylr[-1,:]
    # fill top right and bottom left coordinates
    xx[0,-1]  = xx[1,-1]
    xx[-1,0]  = xx[-2,0]
    yy[0,-1]  = yy[0,-2]
    yy[-1,0]  = yy[-1,1]
    # add extra row and columns for dummy coordinates
    xx = insert(xx,m+1,0,axis=0)
    xx = insert(xx,n+1,0,axis=1)
    yy = insert(yy,m+1,0,axis=0)
    yy = insert(yy,n+1,0,axis=1)

    # prepare lga data
    lga = dw_pcrToDataBlock(scalar(cover(ptid_map,0)))
    # reshape lga to 2D array
    lga = lga.reshape((m,n))
    # add frame of zeros around lga matrix
    lga = insert(lga,m,0,axis=0)
    lga = insert(lga,n,0,axis=1)
    lga = insert(lga,0,0,axis=0)
    lga = insert(lga,0,0,axis=1)
    # find boundary nodes
    outflows = pointer[ pointer[:,1] < 0, 0:2]
    # allocate negative segment number to outflow boundary cells
    for segnr, boundnr in outflows:
        w = where(lga == segnr)
        i,j = w[0][0], w[1][0]
        # choose first zero neighbour and make it
        if lga[i,j-1] == 0:
            lga[i,j-1] = boundnr
        elif lga[i,j+1] == 0:
            lga[i,j+1] = boundnr
        elif lga[i-1,j] == 0:
            lga[i-1,j] = boundnr
        elif lga[i+1,j] == 0:
            lga[i+1,j] = boundnr

    # update grid dimensions
    m = m+2
    n = n+2
    # reshape to 1D arrays
    xx  = xx.reshape( (m*n) )
    yy  = yy.reshape( (m*n) )
    lga = lga.reshape( (m*n) )

    # write cco file
    print("cco.shape: {}",format(xx.shape))
    f = open(fname + '.cco','wb')
    f.write(array([m, n, ], dtype=numpy.int).tostring())
    f.write(array([xx[0], yy[0]], dtype=numpy.float32).tostring())
    f.write(array([0, 0, nolay, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=numpy.int).tostring())
    f.write(hstack((xx, yy)).tostring())
    f.close()

    # write lga file
    nosegh = dw_pcrToDataBlock(ptid_map).shape[0]
    print("lga.shape: {}",format(lga.shape))
    print("m,n : {},{}".format(m,n))
    noq = pointer.shape[0]
    ints = array([n, m, nosegh, nolay, noq, 0, 0])
    ints = hstack((ints, lga))
    ints = array(ints, dtype=numpy.int)
    f = open(fname + '.lga','wb')
    f.write(ints.tostring())
    f.close()


def dw_WriteSurfaceFile(fname,m,n,noseg,block):
    """
    Generates a Delwaq surface (*.srf) file.
    """
    f = open(fname, 'wb')
    f.write(struct.pack('i',m))
    f.write(struct.pack('i',n))
    f.write(struct.pack('i',noseg))
    f.write(struct.pack('i',noseg))
    f.write(struct.pack('i',noseg))
    f.write(struct.pack('i',0))
    f.write(struct.pack('%if'%len(block), *block))
    f.close()


def dw_WriteAttributesFile(fname, noseg):
    """
    Generates a Delwaq atributes (*.atr) file.

    input:
    - fname : file name to write to
    - noseg : number of delwaq segments
    """
    buff  = ""
    buff += "         ; DELWAQ_COMPLETE_ATTRIBUTES\n"
    buff += "    2    ; two blocks with input\n"
    buff += "    1    ; number of attributes, they are :\n"
    buff += "    1    ;  '1' is active '0' is no\n"
    buff += "    1    ; data follows in this fil\n"
    buff += "    1    ; all data is given without defaults\n"
    buff += ";    layer:            1\n"
    buff += "%i*1\n"%noseg
    buff += "    1    ; number of attributes, they are :\n"
    buff += "    2    ;  '1' has surface '3' has bottom\n"
    buff += "         ;  '0' has both    '2' has none\n"
    buff += "    1    ; data follows in this file\n"
    buff += "    1    ; all data is given without defaults\n"
    buff += ";    layer:            1\n"
    buff += "%i*0\n"%noseg
    buff += "    0    ; no time dependent attributes\n"
    f = open(fname, 'w')
    f.write(buff)
    f.close()


def dw_WriteHydFile(fname, d):
    """
    Generates a Delwaq *.hyd file.

    d is dict holding all the required data:
        - d['runid']  : current run id
        - d['tref']   : reference time of simulation as datetime
        - d['tstart'] : start time of simulation as datetime
        - d['tstop']  : stop time of simulation as datetime
        - d['tstep']  : timestep of simulation as timedelta
        - d['m']      : number of grid cells in 1st direction
        - d['n']      : number of grid cells in 2nd direction
    """
    def datetime2str(dt):
        return "{:04}{:02}{:02}{:02}{:02}{:02}".format(dt.year,dt.month,dt.day,dt.hour,dt.minute,dt.second)
    def timedelta2str(td):
        return "{:04}{:02}{:02}{}".format(0,0,td.days,time.strftime('%H%M%S',time.gmtime(td.seconds)))
    buff  = ""
    buff += "task      full-coupling\n"
    buff += "geometry  curvilinear-grid\n"
    buff += "horizontal-aggregation       automatic\n"
    buff += "minimum-vert-diffusion-used  no\n"
    buff += "vertical-diffusion           calculated\n"
    buff += "description\n"
    buff += "'%-60s'\n"%'Generated by Wflow'
    buff += "'%s'\n"%(' ' * 60)
    buff += "'%s'\n"%(' ' * 60)
    buff += "end-description\n"
    buff += "reference-time           '%s'\n"%(datetime2str(d['tref']))
    buff += "hydrodynamic-start-time  '%s'\n"%(datetime2str(d['tstart']))
    buff += "hydrodynamic-stop-time   '%s'\n"%(datetime2str(d['tstop']))
    buff += "hydrodynamic-timestep    '%s'\n"%(timedelta2str(d['tstep']))
    buff += "conversion-ref-time      '%s'\n"%(datetime2str(d['tref']))
    buff += "conversion-start-time    '%s'\n"%(datetime2str(d['tstart']))
    buff += "conversion-stop-time     '%s'\n"%(datetime2str(d['tstop']))
    buff += "conversion-timestep      '%s'\n"%(timedelta2str(d['tstep']))
    buff += "grid-cells-first-direction %7i\n"%d['m']
    buff += "grid-cells-second-direction %6i\n"%d['n']
    buff += "number-hydrodynamic-layers       1\n"
    buff += "number-water-quality-layers      1\n"
    buff += "hydrodynamic-file        none\n"
    buff += "aggregation-file         none\n"
    buff += "grid-indices-file        '%s.lga'\n"%d['runid']
    buff += "grid-coordinates-file    '%s.cco'\n"%d['runid']
    buff += "volumes-file             '%s.vol'\n"%d['runid']
    buff += "areas-file               '%s.are'\n"%d['runid']
    buff += "flows-file               '%s.flo'\n"%d['runid']
    buff += "pointers-file            '%s.poi'\n"%d['runid']
    buff += "lengths-file             '%s.len'\n"%d['runid']
    buff += "salinity-file            none\n"
    buff += "temperature-file         none\n"
    buff += "vert-diffusion-file      none\n"
    buff += "surfaces-file            '%s.srf'\n"%d['runid']
    buff += "depths-file              none\n"
    buff += "total-grid-file          none\n"
    buff += "discharges-file          none\n"
    buff += "chezy-coefficients-file  none\n"
    buff += "shear-stresses-file      none\n"
    buff += "walking-discharges-file  none\n"
    buff += "attributes-file          '%s.atr'\n"%d['runid']
    buff += "constant-dispersion\n"
    buff += "   first-direction    0.0000E+00\n"
    buff += "   second-direction   0.0000E+00\n"
    buff += "   third-direction    0.0000E+00\n"
    buff += "end-constant-dispersion\n"
    buff += "hydrodynamic-layers\n"
    buff += "          1.000\n"
    buff += "end-hydrodynamic-layers\n"
    buff += "water-quality-layers\n"
    buff += "        1.000\n"
    buff += "end-water-quality-layers\n"
    buff += "discharges\n"
    buff += "end-discharges\n"
    f = open(fname, 'w')
    f.write(buff)
    f.close()
   

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
    
    wflow_subcatch = caseId + "/" + configget(config,"model","wflow_subcatch","/staticmaps/wflow_subcatch.map")
    setclone(wflow_subcatch)
    amap = scalar(readmap(caseId +  "/" + areamap))
    modelmap = readmap(wflow_subcatch)
    ldd = readmap(caseId + "/" + configget(config,"model","wflow_ldd","/staticmaps/wflow_ldd.map"))
    gauges = readmap(caseId + "/" + configget(config,"model","wflow_gauges","/staticmaps/wflow_gauges.map"))

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
    
    # write static data for hyd-file set
    comroot = os.sep.join([dwdir,'com',runId])
    mmax, nmax = dw_GetGridDimensions(ptid)
    dw_WritePointer(comroot+'.poi',pointer,binary=True)
    dw_WriteSurfaceFile(comroot+'.srf',mmax,nmax,NOSQ,surface_block)
    dw_WriteSegmentOrExchangeData(0,comroot+'.len',length_block,1,WriteAscii)
    dw_WriteGridFiles(comroot, ptid, pointer)

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
            
            # write dynamic data for hyd-file set
            dw_WriteSegmentOrExchangeData(i,comroot+'.vol',volume_block,1,WriteAscii)
            dw_WriteSegmentOrExchangeData(i,comroot+'.flo',flowblock,1,WriteAscii)
            dw_WriteSegmentOrExchangeData(i,comroot+'.are',area_block,1,WriteAscii)
            
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

        # for hyd-file set
        dw_WriteSegmentOrExchangeData(i,comroot+'.are',area_block,1,WriteAscii)
        dw_WriteSegmentOrExchangeData(i,comroot+'.flo',flowblock,1,WriteAscii)
        dw_WriteSegmentOrExchangeData(i,comroot+'.vol',volume_block,1,WriteAscii)
        
        # Generate attribute file
        atr_file = comroot+'.atr'
        logger.info("Writing attribute file to '%s'"%atr_file)
        dw_WriteAttributesFile(atr_file, NOSQ)

        # Generate hyd-file 
        hyd_file = comroot+'.hyd'
        logger.info("Writing hyd-file to '%s'"%hyd_file)
        hydinfo = {}
        hydinfo['runid'] = runId
        hydinfo['tref'] = T0
        hydinfo['tstart'] = T0
        hydinfo['tstop'] = T0 + timedelta(seconds=(timeSteps-1) * timestepsecs )
        hydinfo['tstep'] = timedelta(seconds=timestepsecs)
        hydinfo['m'], hydinfo['n'] = mmax, nmax
        dw_WriteHydFile(hyd_file, hydinfo)

    
    
if __name__ == "__main__":
    main()
