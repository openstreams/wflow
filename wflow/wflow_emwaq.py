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
Coupling between wflow and D-Emission or D-WAQ.

usage

::

    wflow_emwaq [-C casename][-R runId][-D waqDir][-c configfile]
    [-n netCDF][-T last_step][-S first_step][-s seconds]

    -C: set the name of the wflow_sbm case (directory)

    -R: set the name runId with the wflow_sbm outputs
    
    -D: set the name of the basedir to create the EM/WAQ schematisation in 

    -c: name of wflow the configuration file (default: Casename/wflow_sbm.ini).
    
    -n: name of the wflow netCDF output file (expected in: Casename/runId/).
    If not present, mapstacks will be used.
    
    -i: write an ASCII copy of the binary files
    
    -u: write the structure files
    
    -y: write the dynamic files
    
    -f: write fraction files
    
    -e: emission or water quality coupling
    
    -F: FEWS adapter files
    
    -a: aggregation of results


"""

import numpy as np
import pcraster as pcr

# import pcrut
import sys
import os
import os.path
import shutil
import getopt
import builtins
import struct
import time as t

from wflow.wf_DynamicFramework import *
from wflow.wflow_funcs import *
from wflow.wflow_adapt import *
from wflow import wf_netcdfio
from wflow import pcrut

import pandas as pd
import math


wflow = "wflow_emwaq: "
logger = ""


def usage(*args):
    sys.stdout = sys.stderr
    """Way"""
    for msg in args:
        print(msg)
    print(__doc__)
    sys.exit(0)
    
def dw_CreateDwRun(thedir):
    """"
    Create the dir to save delwaq info in
    Delete flow files if the same runid is used
    Create extra folders if Fews option is enabled.
    """
    
    if not os.path.isdir(thedir):
        #os.makedirs(thedir + "/fixed/")
        os.makedirs(thedir + "/includes_deltashell/")
        os.makedirs(thedir + "/includes_fews/")
        os.makedirs(thedir + "/includes_flow/")
        os.makedirs(thedir + "/debug/")
    if os.path.exists(thedir + "/includes_flow/"):
        flowDir = thedir + "/includes_flow/"
        filelist = os.listdir(flowDir)
        for f in filelist:
            os.remove(os.path.join(flowDir, f))
    # prepare hydfile directory for FEWS
    comdir = os.sep.join([thedir, "com"])
    if os.path.isdir(comdir):
        shutil.rmtree(comdir)
    os.mkdir(comdir)

def dw_pcrToDataBlock(pcrmap, amap, agg = False):
    """
    Converts a pcrmap to a numpy array that is flattend and from which
    missing values are removed using a modelmap. Used for generating delwaq data.
    Can also aggregate data to the subcatchment level by reducing the numpy array.
    
    Input:
        - pcrmap: PCRaster map to turn into a numpy array.
                  If aggregation, map with the same value for all cells in subcatchment.
        - amap: PCRaster map with the model extent to delete NA values. Usually subcatchment map.
                If aggregation, must be the subcatchment map.
        - agg: If False, the numpy array represent values for all model cells.
               If True, returns one value per subcatchment in the numpy array.
    """
    ttar = pcr.pcr2numpy(pcrmap, np.NaN).flatten()
    model = pcr.pcr2numpy(amap, np.NaN).flatten()
    ttar = ttar[np.isfinite(model)]
    
    #If aggregation by subcathcment return the sorted aggregated numpy array of pcrmap
    if agg:
        model = model[np.isfinite(model)]
        model, idx = np.unique(model, return_index=True)
        ttar = ttar[idx]

    return ttar

########################################################################
## Functions for structure  files                                      #
########################################################################
    
def dw_mkDelwaqPointers(ldd, amap, gmap, comp, flux, Emission, Aggregation, bd_id):
    """
    An ldd is used to determine the from-to relations for delwaq using
    the PCraster up/downstreams commands.
    *amap* is used to link boundaries to the segments for delwaq (negative 
    numbers). These are area based boundaries. Diffboun is a 
    python dictionary with inflows for each
    cell.
    
    Input:
        - ldd
        - amap: subcatchment map (also used to determine the active points)
        - gmap: gauges map used for aggregation
        - comp : dataframe of compartments list and properties
        - flux : dataframe of fluxes list and properties
        - Emission : swith for EM or WAQ coupling
        - Aggregation : switch for aggregation by subcatchments
        - bd_id : type of boundary ID numbering (one per boundary or one per boundary and per segment)
    
        
    Output:
        - ptid: map with unique id for each wflow_sbm cell
        - pointer: indicates flow direction from cell id to downstrean cell id
        - pointer_lables: array with number determining the type of flow (lateral, vertical, boundary) for the pointer
        - flow_label: array with the names of the compartment/boundary the flow moves from and to
        - np_ptid.flatten(): array with all the cell ids
        - np_cpid.flatten(): array with the id of each segments (segment = one cell from one compartment)
        - comp_label.flatten(): array with compartments name corresponding to np_cpid
        - flow_var: wflow variable corresponding to the pointer flux
        
 
        
    """
    #For Emission, suppress some of the fluxes that are in the hydrology file but not the pointer file
    if Emission:
        pointerflux = flux.EmPointerHyd.isin(["P","PH","PT","PHT"])
        flux = flux[pointerflux]
    
    ########################################################################
    ## Transform input grids (amap, ldd)                                   #
    ########################################################################
    #Maps with ids for each cell and the flow direction derived from ldd
    ptid = pcr.uniqueid(pcr.boolean(amap))
    flowto = pcr.downstream(ldd, ptid)
    # Firts make sure there is at least one outflow in the model
    # Fix if downsteam is no pit.In that case flowto is missing, set it so itself
    hasflowto = pcr.defined(flowto)
    flowto = pcr.ifthenelse(pcr.defined(ptid) != hasflowto, ptid, flowto)
    
    # find all upstream cells (these must be set negative)
    upbound = pcr.upstream(ldd, 1.0)
    upbound = pcr.ifthen(amap > 0, upbound)
    # Find the lower boundaries (and pits). These flow to themselves
    
    # make into flatted numpy arrays
    np_ptid = pcr.pcr2numpy(ptid, np.NaN).flatten()
    np_flowto = pcr.pcr2numpy(flowto, np.NaN).flatten()
    # Array of subcatchment id
    np_catchid = pcr.pcr2numpy(pcr.scalar(amap), -999).flatten()
    np_upbound = pcr.pcr2numpy(upbound, np.NaN).flatten()
    # Array with gauge id
    np_gmap = pcr.pcr2numpy(pcr.scalar(gmap), np.NaN).flatten()

    # remove all non-active cells and reshape
    np_catchid = np_catchid[np_catchid > 0.0]
    np_upbound = np_upbound[np.isfinite(np_upbound)]
    np_flowto = np_flowto[np.isfinite(np_flowto)]
    np_gmap = np_gmap[np.isfinite(np_ptid)]
    np_ptid = np_ptid[np.isfinite(np_ptid)]
    np_flowto = np_flowto.reshape(len(np_flowto), 1)
    np_ptid = np_ptid.reshape(len(np_ptid), 1)
    np_catchid = np_catchid.reshape(len(np_catchid), 1)
    np_gmap = np_gmap.reshape(len(np_ptid), 1)
    
    if Aggregation:
        #Number of subcatchments
        NOSC = pcr.mapmaximum(amap)
        NOSC = int(np.unique(pcr.pcr2numpy(NOSC, np.NaN))[0])
        #Numpy array with cell (here subcatchment) id
        np_ptid = np.arange(1, NOSC+1).reshape(NOSC, 1)
        #Numpy array with subcatchment id
        np_catchid = np_ptid
        # Downstream subcatchment
        flowtosc = pcr.downstream(ldd, amap)
        flowtosc = pcr.ifthenelse(pcr.defined(amap) != pcr.defined(flowtosc), amap, flowtosc)
        np_flowtosc = pcr.pcr2numpy(pcr.scalar(flowtosc), -999).flatten()
        np_flowtosc = np_flowtosc[np_flowtosc > 0.0]
        np_flowtosc = np_flowtosc.reshape(len(np_flowtosc), 1)
        # Find the downstream subcatchment from the gauge
        np_flowto = []
        for i in np.arange(1,NOSC+1):
            np_flowto = np.append(np_flowto, np_flowtosc[np_gmap == i])  
        np_flowto = np_flowto.reshape(len(np_flowto), 1)                    
    
    ########################################################################
    ## Add compartments and lateral fluxes                                 #
    ########################################################################
    #Name of the compartments to which each cell belong
    comp_label = np.repeat(comp.ID[0],len(np_ptid)).reshape(len(np_ptid),1)
    
    #Add other compartments (cell id + catch id + downstream cell)
    np_cpid = np.arange(1,len(np_ptid)*len(comp)+1).reshape(len(np_ptid)*len(comp),1)
    np_catchidcp = np.tile(np_catchid[:,0], len(comp)).reshape(len(np_ptid)*len(comp),1)
    np_flowtocp = np_flowto
    
    if len(comp) > 1:
        for i in np.arange(1,len(comp)):
            np_flowtocp = np.vstack((np_flowtocp, np_flowto+i*len(np_ptid)))
            comp_label_i = np.repeat(comp.ID[i],len(np_ptid)).reshape(len(np_ptid),1)
            comp_label = np.vstack((comp_label, comp_label_i))
    
    #Set flowto to zero if a compartment has no lateral flux
    latflux = flux.From == flux.To
    compWithLat = flux.From[latflux].astype(str)
    compWithoutLat = ~comp["ID"].isin(compWithLat)
    compWithoutLat = comp.ID[compWithoutLat]

    for n in compWithoutLat:
        np_flowtocp[comp_label == n] = 0 
    
    #Initialize pointer with lateral fluxes and add two zero columns  
    NOLF = len(comp_label[np_flowtocp != 0]) 
    pointer = np.hstack(
            (np_cpid[np_flowtocp != 0].reshape(NOLF,1), 
            np_flowtocp[np_flowtocp != 0].reshape(NOLF,1),
            np.zeros((NOLF, 1)), np.zeros((NOLF, 1)))
            )

    # Pointer labels:
    #    negative: outflow boundary
    #    zero    : internal flow
    #    positive: inflow boundary
    pointer_labels = np.zeros((len(pointer[:,0])), dtype=np.int)
    #Flow labels: indicate from which compartment/boundary the flow moves from and to
    flow_label = (np.char.array(comp_label[np_flowtocp != 0]) + np.char.array(np.repeat("2", NOLF)) +
                  np.char.array(comp_label[np_flowtocp != 0]) + np.char.array(np.repeat("_", NOLF)) +
                  np.char.array(np_catchidcp[np_flowtocp != 0].astype(np.unicode)))
    flow_label = flow_label.reshape(len(flow_label), 1)
    
    #Flow var: indicate the corresponding mapstack name of the flux   
    compWithLat_var = flux.mapstack[latflux].astype(str)
    flow_var = []
    for i in np.arange(0, len(comp)):
        if np.isin(comp.ID[i], compWithLat.values):
            flow_var = np.append(flow_var, compWithLat_var[compWithLat.values == comp.ID[i]].values)
    
    ########################################################################
    ## Remove outflows from pointer                                        #
    ########################################################################
    
    # find all downstream segments (flowto == ptid)
    # now set the flowto points (outflows, usually just one) also to negative
    lowerck = np.absolute(pointer[:,0]) == np.absolute(pointer[:,1])
    
    # outflow to pointer
    # point -> - point
    lopt = pointer[lowerck,0]
    lopt = lopt.reshape(len(lopt), 1)
    zerocol = np.zeros((len(lopt), 1))
    lowerids = np.arange(1, len(lopt) + 1) * -1
    # of = hstack((lopt,lopt * -1.0,zerocol,zerocol))
    lowerids = lowerids.reshape(len(lowerids), 1)
    of = np.hstack((lopt, lowerids, zerocol, zerocol))
    
    # Now remove double pointer to itself and replace by lower boundary
    lowerck = pointer[:, 0] == pointer[:, 1]
    pointer[lowerck, :] = of
    pointer_labels[lowerck] = -1
    #Update flow labels
    lo_label = comp_label[np_flowtocp != 0].reshape(len(comp_label[np_flowtocp != 0]),1)
    lo_label = lo_label[lowerck,0].reshape(len(lopt),1)
    flow_label[lowerck]=(np.char.array(lo_label) + np.char.array(np.repeat("2Outflow", len(lopt)).reshape(len(lopt),1)))
    
    ########################################################################
    ## Add vertical fluxes between compartments                            #
    ########################################################################
    #Find fluxes between compartments in the flux dataframe
    compflux =  (flux.From.isin(comp.ID) & flux.To.isin(comp.ID) & (flux.From != flux.To))
    compflux = flux.loc[compflux,:]
    
    if not compflux.empty:
        #Initialize loop
        comps = 1
        vertcompflux = [] #pointer of the vertical fluxes between compartments
        vertcomp_label = [] #corresponding exchange labels
        zzerocol = np.zeros((len(np_ptid), 1), dtype=np.int)

        for i in np.arange(0,len(compflux)):
            #Get id of the compartment cells exchanging the flux
            comp_flowfrom = np_cpid[comp_label == compflux.From.iloc[i]]
            comp_flowfrom = comp_flowfrom.reshape(len(comp_flowfrom),1)
            comp_flowto = np_cpid[comp_label == compflux.To.iloc[i]]
            comp_flowto = comp_flowto.reshape(len(comp_flowto),1)
            vertcomp_label_i = (np.char.array(np.repeat((compflux.From.iloc[i]+"2"+compflux.To.iloc[i]+"_"), len(comp_flowfrom))) +
                                np.char.array(np_catchid.astype(np.unicode).flatten()))
            vertcomp_label_i = vertcomp_label_i.reshape(len(vertcomp_label_i),1)
            if comps == 1:
                vertcompflux = np.hstack((comp_flowfrom, comp_flowto, zzerocol, zzerocol))
                vertcomp_label = vertcomp_label_i
            else:
                vertcompflux = np.vstack((vertcompflux, np.hstack((comp_flowfrom, comp_flowto, zzerocol, zzerocol))))
                vertcomp_label = np.vstack((vertcomp_label, vertcomp_label_i))
            pointer_labels = np.hstack((pointer_labels, zzerocol[:, 0]))
            comps = comps + 1
    
        pointer = np.vstack((pointer, vertcompflux))
        flow_label = np.vstack((flow_label, vertcomp_label))
        flow_var = np.append(flow_var, compflux.mapstack.values)

    ########################################################################
    ## Add fluxes from/to boundaries                                       #
    ########################################################################
    #Find fluxes from/to boundaries in the flux dataframe
    boundflux =  (~flux.From.isin(comp.ID) | ~flux.To.isin(comp.ID))
    boundflux = flux.loc[boundflux,:]
    
    if not boundflux.empty:
        #Initialize loop
        if lowerids.size == 0:
            start = 1
        else:
            start = np.absolute(lowerids.min()) + 1
        bouns = 1
        extraboun=[]
        bound_label = []
        zzerocol = np.zeros((len(np_ptid), 1), dtype=np.int)

        #Different numbering depending if the flow is from or to the boundary
        for i in np.arange(0, len(boundflux)):
            if bd_id == "segments":
                bounid = np.arange(start, (len(np_ptid) + start)).reshape((len(np_ptid), 1)) * -1.0
            else:
                bounid = np.repeat(start, len(np_ptid)).reshape(len(np_ptid),1) * -1.0
            if np.bool(np.isin(boundflux.From.iloc[i],comp.ID)):
                compid = np_cpid[comp_label == boundflux.From.iloc[i]]
                compid = compid.reshape(len(compid),1)
                if bouns == 1:
                    extraboun = np.hstack((compid, bounid, zzerocol, zzerocol))
                else:
                    extraboun = np.vstack((extraboun, np.hstack((compid, bounid, zzerocol, zzerocol))))
            else:
                compid = np_cpid[comp_label == boundflux.To.iloc[i]]
                compid = compid.reshape(len(compid),1)
                if bouns == 1:
                    extraboun = np.hstack((bounid, compid, zzerocol, zzerocol))
                else:
                    extraboun = np.vstack((extraboun, np.hstack((bounid, compid, zzerocol, zzerocol))))
            if bd_id == "segments":
                bound_label_i = (np.char.array(np.repeat((boundflux.From.iloc[i]+"2"+boundflux.To.iloc[i]+"_"), len(np_ptid))) +
                             np.char.array(np_catchid.astype(np.unicode).flatten()))
            else:
                bound_label_i = np.char.array(np.repeat((boundflux.From.iloc[i]+"2"+boundflux.To.iloc[i]), len(np_ptid))) 
            bound_label_i = bound_label_i.reshape(len(bound_label_i),1)
            if bouns == 1: bound_label = bound_label_i        
            else: bound_label = np.vstack((bound_label, bound_label_i))                    
            pointer_labels = np.hstack((pointer_labels, zzerocol[:, 0] + bouns))
            bouns = bouns + 1
            if bd_id == "segments":
                start = start + len(np_ptid)
            else:
                start = start + 1
        
        pointer = np.vstack((pointer, extraboun))
        flow_label = np.vstack((flow_label, bound_label))
        flow_var = np.append(flow_var, boundflux.mapstack.values)
        
    
    return ptid, pointer, pointer_labels, flow_label, np_ptid.flatten(), np_cpid.flatten(), comp_label.flatten(), flow_var

    
def dw_WriteNrSegments(fname, nr):
    """ 
    Writes the number of segments to B3 file 
    
    B3\_nrofseg.inc
    """
    exfile = open(fname, "w")
    print(";Written by dw_WriteNrSegments", file=exfile)
    print(str(nr) + " ; nr of segments", file=exfile)
    exfile.close()

def dw_WriteAttributes(fname, comp, cells):
    """ 
    Writes the complete delwaq attributes to B3 file 
    
    B3\_attributes.inc
    """
    nb_cells = len(cells)
    
    exfile = open(fname, "w")
    print(";Written by dw_WriteAttributes", file=exfile)
    print("      ; DELWAQ_COMPLETE_ATTRIBUTES", file=exfile)
    print(" 1    ; one block with input", file=exfile)
    print(" 2    ; number of attributes, they are :", file=exfile)
    print("     1     2", file=exfile)
    print(" 1    ; file option in this file", file=exfile)
    print(" 1    ; option without defaults", file=exfile)
    for i in np.arange(0, len(comp)):
        print("     %d*%d%d ; %s" % (nb_cells, comp.At2[i], comp.At1[i], comp.ID[i]), file=exfile)
    print(" 0    ; no time dependent attributes", file=exfile)
    
    
    exfile.close()
    
def dw_WriteNrExChanges(fname, nr):
    """ 
    Writes the number of exchnages to file (number of rows in the pointer file)
    
    B4\_nrofexch.inc
    """
    exfile = open(fname, "w")
    print(";Written by dw_WriteNrExChnages", file=exfile)
    print(str(nr) + " 0 0 ; x, y, z direction", file=exfile)
    exfile.close()

def dw_WritePointer(fname, pointer, binary=False):
    """ 
    WRites the pointer file
    B4\_pointer.inc
    """
    if not binary:
        # Write ASCII file
        exfile = open(fname, "w")
        print(";Written by dw_WritePointer", file=exfile)
        print(";nr of pointers is: ", str(pointer.shape[0]), file=exfile)
        np.savetxt(exfile, pointer, fmt="%10.0f")
        exfile.close()
    else:
        # Write binary file
        f = open(fname, "wb")
        for i in range(pointer.shape[0]):
            f.write(struct.pack("4i", *np.int_(pointer[i, :])))
        f.close()

def dw_WriteBoundlist(fname, flow_labels, pointer_labels, bd_id):
    """ 
    Writes the boundary list file
    B5\_boundlist.inc
    Numbering is abs(exchnage id)
    
    Input: 
        - fname, pointer
        
    """
    exfile = open(fname, "w")
    print(";Written by dw_WriteBoundlist", file=exfile)
    print(";'NodeID' 'Number' 'Type'", file=exfile)
    
    boundloc = pointer_labels != 0
    boundflow_labels = np.char.array(flow_labels[boundloc].flatten())
    
    if bd_id == "segments":
        nrBoundflow = np.arange(1,len(boundflow_labels)+1)
        #Number = char.array(nrBoundflow.astype(numpy.unicode))
    
    else: 
        #Number of outflows
        Outflow = np.char.count(boundflow_labels, "Outflow")
        nrOutflow = len(Outflow[Outflow==1])
        #Number and name of other boundflows
        other_bound = boundflow_labels[Outflow==0]
        bd_labels, idx = np.unique(other_bound, return_index = True)
        bd_labels = other_bound[np.sort(idx)]
        nrBoundflow = len(bd_labels)
        
        #Total boundflow number and name
        if nrOutflow > 0:
            nrBoundflow = nrBoundflow + nrOutflow
            out_labels = boundflow_labels[Outflow==1]
            boundflow_labels = np.append(out_labels, bd_labels)
        else:
            boundflow_labels = bd_labels
        boundflow_labels = np.char.array(boundflow_labels)
            
        #Inputs for the list
        nrBoundflow = np.arange(1,nrBoundflow+1)
               
    Number = np.char.array(nrBoundflow.astype(np.unicode))
    BD = np.char.array(np.repeat("'BD_", len(Number)))
    quote = np.char.array(np.repeat("'", len(Number)))
    NodeID = BD + Number + quote
    Number = quote + Number + quote
    boundflow_labels = quote + boundflow_labels + quote

    boundlist = np.hstack((NodeID.reshape(len(NodeID),1), Number.reshape(len(Number),1), boundflow_labels.reshape(len(Number),1)))
    
    np.savetxt(exfile, boundlist, fmt="%.20s")
    
    exfile.close()
    
    return boundflow_labels


########################################################################
## Functions for monitoring files                                      #
########################################################################

def dw_mkDelwaqMonitoring(mon_points, mon_areas, gmap, amap, comp, cells, segments, seg_labels):
    """
    Calculates the number of points/monitoring areas and returns the stations list
    """
    np_gmap = pcr.pcr2numpy(gmap, np.NaN).flatten()    
    np_amap = pcr.pcr2numpy(amap, np.NaN).flatten()
    np_gmap = np_gmap[np.isfinite(np_amap)]
    np_amap = np_amap[np.isfinite(np_amap)]
    
    if mon_points == "segments":
        NOPT = segments.shape[0]
    elif mon_points == "gauges":
        NOPT = len(np_gmap[np.isfinite(np_gmap)])*len(comp)
    else: NOPT = 0

    if mon_areas == "subcatch":
        NOAR = len(np.unique(np_amap))
    elif mon_areas == "compartments":
        NOAR = len(comp)
    else: NOAR = 0
    
    if mon_points == "gauges":
        ID = "Gauge" 
        gauges_comp = np.tile(np_gmap, len(comp))
        cell_ID = np.int_(gauges_comp[np.isfinite(gauges_comp)])
        comp_list = seg_labels[np.isfinite(gauges_comp)]
        point_list = np.int_(segments[np.isfinite(gauges_comp)])    
    elif mon_points == "segments":
        ID = "Cell"
        cell_ID = np.int_(np.tile(cells, len(comp)))
        comp_list = seg_labels
        point_list = np.int_(segments)
    
    if NOPT > 0:
        outlocname = (np.char.array(np.repeat("'", NOPT)) + np.char.array(np.repeat(ID, NOPT)) + 
                      np.char.array(cell_ID.astype(np.unicode)) + np.char.array(np.repeat("_", NOPT)) +
                      np.char.array(comp_list) + np.char.array(np.repeat("'", NOPT)))
        onecol = np.repeat(1, NOPT).reshape(NOPT,1)
        balcol = np.repeat("NO_BALANCE", NOPT).reshape(NOPT,1)
        stations = np.hstack((outlocname.reshape(NOPT,1), balcol, onecol, point_list.reshape(NOPT,1)))
        stations_balance = np.hstack((outlocname.reshape(NOPT,1), onecol, point_list.reshape(NOPT,1)))
    else:
        stations = []
        stations_balance = []
    
    return NOPT, NOAR, stations, stations_balance
    

def dw_WriteNrMonAreas(fname, NOPT, NOAR, mon_points, mon_areas):
    """
    Write an output loc file with the number of monitoring locations (points+areas)
    map.
    """
    pts = NOPT + NOAR
    
    exfile = open(fname, "w")
    print(";Written by dw_NrMonAreas", file=exfile)
    print(";%d monitoring points (%s)" % (NOPT, mon_points), file=exfile)
    print(";%d monitoring areas (%s)" % (NOAR, mon_areas), file=exfile)
    print("%d ; nr of monitoring points/areas" % pts, file=exfile)
    exfile.close()

def dw_WriteStations(fname, stations):
    """
    Write an output loc file based on the monitoring points (stations)
    map.
    """
 
    exfile = open(fname, "w")
    print(";Written by dw_WriteStations", file=exfile)
    np.savetxt(exfile, stations, fmt="%.20s")
    exfile.close()
    
def dw_WriteMonAreas(fname, NOAR, mon_areas, amap, comp, segments, seg_labels, cells):
    """
    Write an output loc file based on the monitoring areas
    map.
    """
    np_amap = dw_pcrToDataBlock(amap, amap)
#    np_amap = pcr2numpy(amap, NaN).flatten()
#    np_amap = np_amap[isfinite(np_amap)]
    
    exfile = open(fname, "w")
    print(";Written by dw_WriteMonAreas", file=exfile)
    
    if mon_areas == "subcatch":
        area_ID = (np.char.array(np.repeat("'", NOAR)) + np.char.array(np.repeat("Subcatch", NOAR)) + 
            np.char.array(np.int_(np.unique(np_amap)).astype(np.unicode)) + np.char.array(np.repeat("'", NOAR)))
        for i in np.arange(1, NOAR+1):
            area_len = len(np_amap[np_amap == i])*len(comp)
            subcatch_comp = np.tile(np_amap, len(comp))
            area_list = np.int_(segments[subcatch_comp == i])
            #area_list = area_list.reshape(len(comp), len(np_amap[np_amap == i]))
            
            #Reshape area list as max characters /line in ASCII file is 1000
            #Max allowed number ID of cells has 20 characters -> 50 cells / row
            NTOT = len(area_list)
            #Number of complete rows
            NCOMP = int(len(area_list)/50)
            area_list_1 = area_list[0:NCOMP*50].reshape(NCOMP, 50)  
            area_list_2 = area_list[NCOMP*50:NTOT]
            area_list_2 = area_list_2.reshape(1, len(area_list_2))
            
            print(area_ID[i-1], area_len, file = exfile, sep = "        ")
            np.savetxt(exfile, area_list_1, fmt="%10.20s")
            np.savetxt(exfile, area_list_2, fmt="%10.20s")
    else:
        area_ID, idx = np.unique(seg_labels, return_index = True)
        area_ID = seg_labels[np.sort(idx)]
        area_name = (np.char.array(np.repeat("'", NOAR)) + np.char.array(area_ID) + 
                   np.char.array(np.repeat("'", NOAR)))
        area_len = len(cells)
        for i in range(0, NOAR):
            area_list = np.int_(segments[seg_labels == area_ID[i]])
            #area_list = area_list.reshape(1, len(area_list))
            
            #Reshape area list as max characters /line in ASCII file is 1000
            #Max allowed number ID of cells has 20 characters -> 50 cells / row
            NTOT = len(area_list)
            #Number of complete rows
            NCOMP = int(len(area_list)/50)
            area_list_1 = area_list[0:NCOMP*50].reshape(NCOMP, 50)  
            area_list_2 = area_list[NCOMP*50:NTOT]
            area_list_2 = area_list_2.reshape(1, len(area_list_2))
            
            print(area_name[i], area_len, file = exfile, sep = "        ")
            np.savetxt(exfile, area_list_1, fmt="%10.20s")
            np.savetxt(exfile, area_list_2, fmt="%10.20s")
    
    
    exfile.close()

########################################################################
## Functions for emission structure files                              #
########################################################################

def dw_WriteNrSegL(fname, NOSQ, NOCP):
    """
    Write the number of segments per layer to file
    nrofsegl.inc
    """
    exfile = open(fname, "w")
    nb_cells = int(NOSQ/NOCP)
    print("CONSTANTS nsca DATA        %d" % nb_cells, file=exfile)
    print("CONSTANTS nrec DATA        %d" % NOCP, file=exfile)
    print("PARAMETERS Surf ALL DATA        %d*1.0" % NOSQ, file=exfile)
    
    exfile.close()
    
def dw_WriteGeometry(fname, cells, amap, reallength, PathFrac, WaterFrac, Aggregation, WriteAscii, NOCP):
    """
    Write the geometry to file
    geometry.inc
    """
    cell_list = np.int_(cells)
    #Reshape cell list as max characters /line in ASCII file is 1000
    #Max allowed number ID of cells has 20 characters -> 50 cells / row
    NTOT = len(cell_list)
#    if NTOT >= 100:
#        #Number of complete rows
#        NCOMP = int(len(cell_list)/100)
#        cell_list_1 = cell_list[0:NCOMP*100].reshape(NCOMP, 100)  
#        cell_list_2 = cell_list[NCOMP*100:NTOT]
#        cell_list_2 = cell_list_2.reshape(1, len(cell_list_2))
#    else:
#        cell_list = cell_list.reshape(1, len(cell_list))
    
    #Transforms properties into flat numpy array
    np_amap = pcr.pcr2numpy(amap, np.NaN).flatten()
    reallength = pcr.pcr2numpy(reallength, np.NaN).flatten()
    PathFrac = pcr.pcr2numpy(PathFrac, np.NaN).flatten()
    WaterFrac = pcr.pcr2numpy(WaterFrac, np.NaN).flatten()
    #Reduce pathfrac and waterfrac to active cells
    reallength = reallength[np.isfinite(np_amap)]
    PathFrac = PathFrac[np.isfinite(np_amap)]
    WaterFrac = WaterFrac[np.isfinite(np_amap)]
    np_amap = np_amap[np.isfinite(np_amap)]
    
    realarea = reallength * reallength
    WaterFrac =  np.around(WaterFrac, decimals=3)
    fPaved = np.around(PathFrac * (np.repeat(1.0, len(PathFrac)) - WaterFrac), decimals=3)
    fUnPaved = np.around((np.repeat(1.0, len(PathFrac)) - PathFrac)*(np.repeat(1.0, len(PathFrac)) - WaterFrac), decimals=3)
    
    if Aggregation:
        realareasc = []
        WaterFracsc = []
        fPavedsc = []
        fUnPavedsc = []
        for i in range(1, len(cells)+1):
            sc_area = sum(realarea[np_amap == i])
            realareasc = np.append(realareasc, sc_area)
            WaterFracsc = np.append(WaterFracsc, sum(WaterFrac[np_amap == i]*realarea[np_amap == i])/sc_area)
            fPavedsc = np.append(fPavedsc, sum(fPaved[np_amap == i]*realarea[np_amap == i])/sc_area)
            fUnPavedsc = np.append(fUnPavedsc, sum(fUnPaved[np_amap == i]*realarea[np_amap == i])/sc_area)
        realarea = realareasc
        WaterFrac = WaterFracsc
        fPaved = fPavedsc
        fUnPaved = fUnPavedsc
    
    
    outlocname = (np.char.array(np.repeat("; ", len(cells))) + np.char.array(np.repeat("Cell", len(cells))) + 
                  np.char.array((np.int_(cells)).astype(np.unicode)))
    
    if WriteAscii:
        geometry_data = np.hstack((realarea.reshape(len(cells),1), fPaved.reshape(len(cells),1),
                            fUnPaved.reshape(len(cells),1), WaterFrac.reshape(len(cells),1),
                            outlocname.reshape(len(cells),1)))
        #Put zeros for the remaining compartments
        if NOCP > 1:
            extra_data = str(int(len(cells)*(NOCP-1)*4)) + "*0.0"
        exfile = open(fname + "geometry.inc", "w")        
        print("PARAMETERS TotArea fPaved fUnpaved fOpenWater ALL", file=exfile)
#        print("PARAMETERS TotArea fPaved fUnpaved fOpenWater SEGMENTS", file=exfile)
#        if NTOT >= 100:
#            savetxt(exfile, cell_list_1, fmt="%6.10s")
#            savetxt(exfile, cell_list_2, fmt="%6.10s")
#        else:
#            savetxt(exfile, cell_list, fmt="%6.10s")
        print("DATA", file=exfile)
        np.savetxt(exfile, geometry_data, fmt="%10.10s")
        if NOCP >1: print(extra_data, file=exfile)            
        exfile.close()
    
    # Write binary file
    #Flatten the geometry data and repeat them for each compartment
    geometry_data = np.hstack((realarea.reshape(len(cells),1), fPaved.reshape(len(cells),1),
                            fUnPaved.reshape(len(cells),1), WaterFrac.reshape(len(cells),1)))
    geometry_data = np.tile(geometry_data.flatten(), NOCP)
    artow = np.array(geometry_data, dtype=np.float32).copy()
    #Define dummy time
    timear = np.array(0, dtype=np.int32)
    #Open and write the data
    fp = open(fname + "geometry.bin", "wb")
    tstr = timear.tostring() + artow.tostring()
    fp.write(tstr)
    fp.close()
    
    #Write corresponding def file of the geometry
    fpa = open(fname + "parameters.inc", "w")        
    print("PARAMETERS TotArea fPaved fUnpaved fOpenWater", file=fpa)
    fpa.close()
    
    
    
def dw_WriteHydDef(fname, cells, flux):
    """
    Write the hydrology definition file to
    hydrology.inc
    """
    #Define the names of the fluxes saved in hydrology.bin
    #Add an f to the name for the fluxes defined as a fraction of the volume
    fracname = np.repeat("", len(flux))
    fracname[flux.EmFraction == 1] = "f"
    hydflux = flux.EmPointerHyd.isin(["H","PH","HT","PHT"])
    hydfluxname = flux.Name[hydflux].values.astype(str)
    fracname = fracname[hydflux]
    
    hydfluxname = (np.char.array(fracname) + np.char.array(hydfluxname))
    #Add TotalFlow
    hydfluxname = np.append(hydfluxname, "TotalFlow")
    hydfluxname = hydfluxname.reshape(1, len(hydfluxname))
    
#    cell_list = int_(cells)
#    #Reshape cell list as max characters /line in ASCII file is 1000
#    #Max allowed number ID of cells has 10 characters -> 100 cells / row
#    NTOT = len(cell_list)
#    if NTOT >= 100:
#        #Number of complete rows
#        NCOMP = int(len(cell_list)/100)
#        cell_list_1 = cell_list[0:NCOMP*100].reshape(NCOMP, 100)  
#        cell_list_2 = cell_list[NCOMP*100:NTOT]
#        cell_list_2 = cell_list_2.reshape(1, len(cell_list_2))
#    else:
#        cell_list = cell_list.reshape(1, len(cell_list))
    
    
    #Create and print file
    exfile = open(fname, "w")
    
    print("SEG_FUNCTIONS", file=exfile)
    np.savetxt(exfile, hydfluxname, fmt="%10.15s")
#    print("SEGMENTS", file=exfile)
#    if NTOT >= 100:
#        savetxt(exfile, cell_list_1, fmt="%6.10s")
#        savetxt(exfile, cell_list_2, fmt="%6.10s")
#    else:
#        savetxt(exfile, cell_list, fmt="%6.10s")
    
    exfile.close()
    
    

########################################################################
## Functions for fraction files                                        #
########################################################################

def dw_WriteBoundData(fname, boundflow_labels):
    """ 
    writes B5\_bounddata.inc
    """
    
    exfile = open(fname, "w")
    print(";Written by dw_WriteBoundData", file=exfile)
    for i in boundflow_labels:
        print("ITEM '%s'" % (i), file=exfile)
        print("CONCENTRATION  '%s' 'Check' 'Initial'" % (i), file=exfile)
        print("DATA", file=exfile)
        print("1.0  1.0  0.0", file=exfile)
        print("", file=exfile)

    exfile.close()


def dw_WriteInitials(fname, boundflow_labels):
    """
    B8_initials.inc
    """

    maps = ["Initial", "Check"]
    exfile = open(fname, "w")
    print("INITIALS", file=exfile)
    for rr in boundflow_labels:
        print("'" + rr + "'", end=" ", file=exfile)
    for rr in maps:
        print("'" + rr + "'", end=" ", file=exfile)
    print(file=exfile)
    print("DEFAULTS", file=exfile)
    for rr in boundflow_labels:
        print(str(0.0) + " ", end=" ", file=exfile)
    for rr in maps:
        print(str(1.0) + " ", end=" ", file=exfile)
    print(file=exfile)
    exfile.close()

def dw_Write_Substances(fname, boundflow_labels):
    """
    Writes the B1_sublist.inc file
    input:
        
        it writes substances for the areas and an initial and mass balance 
        check substance
        
    """
    
    exfile = open(fname, "w")
    print("; number of active and inactive substances", file=exfile)
    print("%d         0" % (len(boundflow_labels) + 2), file=exfile)
    print("; active substances", file=exfile)
    print("1             'Initial' ; ", file=exfile)
    print("2             'Check' ; ", file=exfile)
    j = 2
    for i in boundflow_labels:
        j = j + 1
        print("%d            '%s'" % (j, i), file=exfile)
    print("; passive substances", file=exfile)

    exfile.close()


########################################################################
## Functions for FEWS files                                            #
########################################################################

def dw_WriteFractionParameterIds(fname, boundflow_labels):
    """
    Writes the qualifier_fractions.csv file
    input:
        
        it writes qualifiers for the areas
        
    """    
    exfile = open(fname,'w')
    print("id;name;allow_missing;group;type;unit;value_resolution", file=exfile)
    print("Initial;Fraction of initial water;TRUE;WaterFraction;instantaneous;-;0.001", file=exfile) 
    for i in boundflow_labels:
        print("%s;Fraction of water from source %s;TRUE;WaterFraction;instantaneous;-;0.001" %  (i, i), file=exfile) 
        
    exfile.close()
    
def dw_Write_B2_outlocs(fname, gauges, segs):
    """
    Write an output loc file based on the wflow_gauges
    map.
    """
    segs = pcr.ifthenelse(gauges > 0, segs, np.NaN)
    gauges = pcr.ifthenelse(gauges > 0, pcr.scalar(gauges), np.NaN)
    np_gauges = pcr.pcr2numpy(gauges, np.NaN).flatten()
    np_segs = pcr.pcr2numpy(segs, np.NaN).flatten()
        
    np_gauges = np_gauges[np.isfinite(np_gauges)]
    np_segs = np_segs[np.isfinite(np_segs)]
    
    if len(np_segs) != len(np_gauges):
        logger.error("Gauges and segments do not match!")

    pts = np.size(np_segs)
    exfile = open(fname, "w")
    print("%d ; nr of locations" % pts, file=exfile)
    print("; 'outlocname' numberofsegments segment list", file=exfile)
    i = 0
    for loc in np_gauges:
        print(" '%d' 1 %d" % (loc, np_segs[i]), file=exfile)
        i = i + 1
    exfile.close()

def dw_GetGridDimensions(ptid_map):
    """
    Returns number of cells in 1st and 2nd grid directions.

    input:
    - ptid_map : PCRaster map with unique id's
    """
    # find number of cells in m and n directions
    zero_map = pcr.scalar(ptid_map) * 0.0
    allx = dw_pcrToDataBlock(pcr.xcoordinate(pcr.boolean(pcr.cover(zero_map + 1, 1))), zero_map)
    i = 0
    diff = round(builtins.abs(allx[i] - allx[i + 1]), 5)
    diff_next = diff
    while diff_next == diff:
        i += 1
        diff_next = builtins.abs(allx[i] - allx[i + 1])
        diff_next = round(diff_next, 5)
    m = i + 1
    n = allx.shape[0] / m
    m, n = n, m
    return m, n

def dw_WriteSurfaceFile(fname, block):
    """
    Generates a Delwaq surface (*.srf) file.
    """
    f = open(fname, "wb")
    f.write(struct.pack("i", 0))
    f.write(struct.pack("%if" % len(block), *block))
    f.close()
    
def dw_WriteWaqGeom(fname, ptid_map, ldd_map, FEWS):
    """
    Writes Delwaq netCDF geometry file (*_waqgeom.nc).

    input:
    - fname    : output file name (without file extension)
    - ptid_map : PCRaster map with unique id's
    """
    # Get coordinates

    zero_map = pcr.scalar(ptid_map) * 0.0
    pcr.setglobaloption("coorul")  # upper-left cell corners
    xxul = pcr.pcr2numpy(pcr.xcoordinate(pcr.boolean(pcr.cover(zero_map + 1, 1))), -1)
    yyul = pcr.pcr2numpy(pcr.ycoordinate(pcr.boolean(pcr.cover(zero_map + 1, 1))), -1)
    pcr.setglobaloption("coorlr")  # lower-right cell corners
    xxlr = pcr.pcr2numpy(pcr.xcoordinate(pcr.boolean(pcr.cover(zero_map + 1, 1))), -1)
    yylr = pcr.pcr2numpy(pcr.ycoordinate(pcr.boolean(pcr.cover(zero_map + 1, 1))), -1)

    # Convert pcr maps to numpy arrays

    np_ptid = pcr.pcr2numpy(ptid_map, -1)
    np_ldd = pcr.pcr2numpy(ldd_map, -1)
    np_ldd[np_ldd == 255] = 0

    # Number of segments in horizontal dimension

    nosegh = int(np.max(np_ptid))

    # Waqgeom dimensions

    n_net_node = 0
    n_net_link = 0
    n_net_link_pts = 2
    n_net_elem = nosegh
    n_net_elem_max_node = 4 # all elements are rectangles
    n_flow_link = nosegh - 1 # one per element, except for outlet
    n_flow_link_pts = 2

    # Prepare waqgeom data structures

    nodes_x = []
    nodes_y = []
    nodes_z = []
    net_links = []
    elem_nodes = np.zeros((n_net_elem, n_net_elem_max_node), dtype=np.int)
    flow_links = np.zeros((n_flow_link, n_flow_link_pts), dtype=np.int)
    flow_link_x = np.zeros((n_flow_link), dtype=np.float)
    flow_link_y = np.zeros((n_flow_link), dtype=np.float)
	
    # prepare all coordinates for grid csv
    nodes_x_all = []
    nodes_y_all = []
	
    # Keep track of nodes and links as dataset grows 

    i_node = 0  # index of last node
    i_flink = 0 # index of last flow link

    # PCR cell id's start at 1, we need it zero based

    np_ptid = np_ptid - 1

    # Wflow map dimensions

    m, n = np_ptid.shape

    # Helper function

    def add_node(i, j, corner):
        # Get coordinates
        if corner == UL:
            x = xxul[i, j]
            y = yyul[i, j]
        elif corner == LR:
            x = xxlr[i, j]
            y = yylr[i, j]
        elif corner == UR:
            x = xxlr[i, j]
            y = yyul[i, j]
        elif corner == LL:
            x = xxul[i, j]
            y = yylr[i, j]
        else:
            assert 0
        # Add node coordinates
        nodes_x.append(x)
        nodes_y.append(y)
        nodes_z.append(0)

    def add_all_nodes(i, j, corner):
        # Get coordinates
        if corner == UL:
            x = xxul[i, j]
            y = yyul[i, j]
        elif corner == LR:
            x = xxlr[i, j]
            y = yylr[i, j]
        elif corner == UR:
            x = xxlr[i, j]
            y = yyul[i, j]
        elif corner == LL:
            x = xxul[i, j]
            y = yylr[i, j]
        else:
            assert(0)
        # Add node coordinates
        nodes_x_all.append(x)
        nodes_y_all.append(y)
		
    # Cell corners

    UL, UR, LR, LL = 0, 1, 2, 3

    # Process all cells from upper-left to lower-right

    for i in range(m):
        for j in range(n):
            # Current element index
            i_elem = int(np_ptid[i, j])
            if i_elem < 0:
                # Skip inactive segment
                continue

            # Get index of neighbouring elements that could have been processed before

            if i == 0:
                i_elem_up_left  = -1
                i_elem_up       = -1
                i_elem_up_right = -1
            elif j == 0:
                i_elem_up_left  = -1
                i_elem_up = int(np_ptid[i - 1, j])
                i_elem_up_right = int(np_ptid[i - 1, j + 1])
            elif j == n - 1:
                i_elem_up_left = int(np_ptid[i - 1, j - 1])
                i_elem_up = int(np_ptid[i - 1, j])
                i_elem_up_right = -1
            else:
                i_elem_up_left = int(np_ptid[i - 1, j - 1])
                i_elem_up = int(np_ptid[i - 1, j])
                i_elem_up_right = int(np_ptid[i - 1, j + 1])
            
            if j == 0:
                i_elem_left = -1
            else:
                i_elem_left = int(np_ptid[i, j - 1])
            
            # Update nodes:
            # If left or upper neighbours are active, some nodes of current cell
            # have been added already.

            # UL node
            if i_elem_left < 0 and i_elem_up_left < 0 and i_elem_up < 0:
                add_node(i, j, UL)
                elem_nodes[i_elem, UL] = i_node
                i_node += 1
            elif i_elem_left >= 0:
                elem_nodes[i_elem, UL] = elem_nodes[i_elem_left, UR]
            elif i_elem_up_left >= 0:
                elem_nodes[i_elem, UL] = elem_nodes[i_elem_up_left, LR]
            elif i_elem_up >= 0:
                elem_nodes[i_elem, UL] = elem_nodes[i_elem_up, LL]

            # UR node
            if i_elem_up < 0 and i_elem_up_right < 0:
                add_node(i, j, UR)
                elem_nodes[i_elem, UR] = i_node
                i_node += 1
            elif i_elem_up >= 0:
                elem_nodes[i_elem, UR] = elem_nodes[i_elem_up, LR]
            elif i_elem_up_right >= 0:
                elem_nodes[i_elem, UR] = elem_nodes[i_elem_up_right, LL]
            if i_elem_up < 0:
                # add UL-UR link
                net_links.append((elem_nodes[i_elem, UL], elem_nodes[i_elem, UR]))

            # LL node
            if i_elem_left < 0:
                add_node(i, j, LL)
                elem_nodes[i_elem, LL] = i_node
                i_node += 1
                # add UL-LL link
                net_links.append((elem_nodes[i_elem, UL], elem_nodes[i_elem, LL]))
            else:
                elem_nodes[i_elem, LL] = elem_nodes[i_elem_left, LR]

            # LR node
            add_node(i, j, LR)
            add_all_nodes(i, j, LR)
            elem_nodes[i_elem, LR] = i_node
            i_node += 1
            # add LL-LR link
            net_links.append((elem_nodes[i_elem, LL], elem_nodes[i_elem, LR]))
            # add UR-LR link
            net_links.append((elem_nodes[i_elem, UR], elem_nodes[i_elem, LR]))

            # Update flow links based on local drain direction
            # TODO: diagonal flow links between cells that have only one node in common?
            
            direction = np_ldd[i, j]
            i_other = -1
            if direction == 1:
                i_other = np_ptid[i + 1, j - 1]  # to lower left
            elif direction == 2:
                i_other = np_ptid[i + 1, j]  # to lower
            elif direction == 3:
                i_other = np_ptid[i + 1, j + 1]  # to lower right
            elif direction == 4:
                i_other = np_ptid[i, j - 1]  # to left
            elif direction == 6:
                i_other = np_ptid[i, j + 1]  # to right
            elif direction == 7:
                i_other = np_ptid[i - 1, j - 1]  # to upper right
            elif direction == 8:
                i_other = np_ptid[i - 1, j]  # to upper
            elif direction == 9:
                i_other = np_ptid[i - 1, j + 1]  # to upper left
            if i_other >= 0:
                flow_links[i_flink, :] = i_elem, i_other
                i_flink += 1

	# Convert data to numpy arrays
    nodes_x_all = np.array(nodes_x_all)
    nodes_y_all = np.array(nodes_y_all)
    
    nodes_x = np.array(nodes_x)
    nodes_y = np.array(nodes_y)
    nodes_z = np.array(nodes_z)
    net_links = np.array(net_links)
    
    # Update dimensions

    n_net_node = nodes_x.shape[0]
    n_net_link = net_links.shape[0]

    # Create netCDF file in classic format

    f = netCDF4.Dataset(fname + "_waqgeom.nc", "w", format="NETCDF3_CLASSIC")

    # Create dimensions

    f.createDimension("dim", 1)
    f.createDimension("nNetNode", n_net_node)
    f.createDimension("nNetLink", n_net_link)
    f.createDimension("nNetLinkPts", n_net_link_pts)
    f.createDimension("nNetElem", n_net_elem)
    f.createDimension("nNetElemMaxNode", n_net_elem_max_node)
    f.createDimension("nFlowLink", n_flow_link)
    f.createDimension("nFlowLinkPts", n_flow_link_pts)

    # Create variables

    v_msh = f.createVariable("mesh", "i4", ("dim",))
    v_pcs = f.createVariable("projected_coordinate_system", "i4", ())
    v_nnx = f.createVariable("NetNode_x", "f8", ("nNetNode",))
    v_nny = f.createVariable("NetNode_y", "f8", ("nNetNode",))
    v_nnz = f.createVariable("NetNode_z", "f8", ("nNetNode",))
    v_nlk = f.createVariable("NetLink", "i4", ("nNetLink", "nNetLinkPts"))
    v_nen = f.createVariable(
        "NetElemNode", "i4", ("nNetElem", "nNetElemMaxNode"), fill_value=0
    )
    v_flk = f.createVariable("FlowLink", "i4", ("nFlowLink", "nFlowLinkPts"))
    v_flt = f.createVariable("FlowLinkType", "i4", ("nFlowLink",))
    v_flx = f.createVariable("FlowLink_xu", "f8", ("nFlowLink",))
    v_fly = f.createVariable("FlowLink_yu", "f8", ("nFlowLink",))

    # Variable attributes

    v_msh.long_name = "Delft3D FM aggregated mesh"
    v_msh.cf_role = "mesh_topology"
    v_msh.topology_dimension = "2 d"
    v_msh.node_coordinates = "NetNode_x NetNode_y"
    v_msh.face_node_connectivity = "NetElemNode"
    v_msh.edge_node_connectivity = "NetLink"
    v_msh.edge_face_connectivity = "FlowLink"
    v_msh.face_dimension = "nNetElem"
    v_msh.edge_dimension = "nNetLink"
    v_msh.node_dimension = "nNetNode"
    v_msh.face_coordinates = "Face_x Face_y"
    v_msh.edge_coordinates = "FlowLink_xu FlowLink_yu"

    # v_pcs.name = "Unknown projected"
    v_pcs.epsg = 4326
    v_pcs.grid_mapping_name = "Unknown projected"
    v_pcs.longitude_of_prime_meridian = 0.
    # v_pcs.semi_major_axis = 6378137.
    # v_pcs.semi_minor_axis = 6356752.314245
    v_pcs.inverse_flattening = 298.257223563
    v_pcs.epsg_code = "EPSG:4326"
    v_pcs.value = "value is equal to EPSG code"

    v_nnx.units = "degrees_east"
    v_nnx.standard_name = "longitude"
    v_nnx.long_name = "longitude"

    v_nny.units = "degrees_north"
    v_nny.standard_name = "latitude"
    v_nny.long_name = "latitude"

    v_nnz.units = "m"
    v_nnz.positive = "up"
    v_nnz.standard_name = "sea_floor_depth"
    v_nnz.long_name = "Bottom level at net nodes (flow element's corners)"
    v_nnz.coordinates = "NetNode_x NetNode_y"

    v_nlk.long_name = "link between two netnodes"
    v_nlk.start_index = 1

    v_nen.long_name = "Net element defined by nodes"
    v_nen.start_index = 1
    # v_nen._FillValue = 0

    v_flk.long_name = "link/interface between two flow elements"
    v_flk.start_index = 1

    v_flt.long_name = "type of flowlink"
    v_flt.valid_range = 1, 2
    v_flt.flag_values = 1, 2
    v_flt.flag_meanings = "link_between_1D_flow_elements link_between_2D_flow_elements"

    v_flx.units = "degrees_east"
    v_flx.standard_name = "longitude"
    v_flx.long_name = "x-Coordinate of velocity point on flow link."

    v_fly.units = "degrees_north"
    v_fly.standard_name = "latitude"
    v_fly.long_name = "y-Coordinate of velocity point on flow link."

    # Global attributes

    f.institution = "Deltares"
    f.references = "http://www.deltares.nl"
    time_string = t.strftime("%b %d %Y, %H:%M:%S")
    f.source = "Wflow, Deltares, %s." % time_string
    offset_s = -t.altzone
    offset_m = int((offset_s % 3600) / 60)
    offset_h = int((offset_s / 60 - offset_m) / 60)
    time_string = t.strftime("%Y-%m-%dT%H:%M:%S") + "+%02i%02i" % (
        offset_h,
        offset_m,
    )
    f.history = "Created on %s, wflow_delwaq.py" % time_string
    f.Conventions = "CF-1.6 UGRID-0.9"

    # Data

    v_nnx[:] = nodes_x
    v_nny[:] = nodes_y
    v_nnz[:] = nodes_z
    v_nlk[:, :] = net_links + 1  # uses 1-based indexes
    v_nen[:, :] = elem_nodes + 1  # uses 1-based indexes
    v_flk[:, :] = flow_links + 1  # uses 1-based indexes
    v_flt[:] = 2
    v_flx[:] = 0
    v_fly[:] = 0

    f.close()

    if FEWS:
        coordinates = pd.DataFrame({'X': nodes_x_all, 'Y': nodes_y_all})
        dw_WriteFractionsGrid(fname + "_fraction_grid.csv", coordinates)

def dw_WriteFractionsGrid(fname, coordinates):
    """
    Write csv file for FEWS grids with x, y coordinates
    """
    f = open(fname, 'w')	
    coordinates.to_csv(f, sep=';', index = False, header=True)
    f.close()
    
def dw_WriteBndFile(fname, ptid_map, pointer, pointer_labels, flow_labels):
    """
    Writes Delwaq *.bnd file.

    input:
    - fname          : output file name (without file extension)
    - ptid_map       : PCRaster map with unique id's
    - pointer        : delwaq pointers
    - pointer_labels : numpy array with pointer types
    - areas          : area id per inflow
    - source_ids     : list of source names

    A unique boundary is generated per source for all segments in a given area.
    A unique boundary is generated for each outflow.
    """
    buff = ""
    np_ptid = pcr.pcr2numpy(ptid_map, -1)
    #area_ids = unique(areas)
    
    boundloc = pointer_labels > 0
    boundflow_labels = np.unique(np.char.array(flow_labels[boundloc].flatten()))

    # Upper-left and lower-right Coordinates

    zero_map = pcr.scalar(ptid_map) * 0.0
    pcr.setglobaloption("coorul")  # upper-left cell corners
    xxul = pcr.pcr2numpy(pcr.xcoordinate(pcr.boolean(pcr.cover(zero_map + 1, 1))), -1)
    yyul = pcr.pcr2numpy(pcr.ycoordinate(pcr.boolean(pcr.cover(zero_map + 1, 1))), -1)
    pcr.setglobaloption("coorlr")  # lower-right cell corners
    xxlr = pcr.pcr2numpy(pcr.xcoordinate(pcr.boolean(pcr.cover(zero_map + 1, 1))), -1)
    yylr = pcr.pcr2numpy(pcr.ycoordinate(pcr.boolean(pcr.cover(zero_map + 1, 1))), -1)

    # Map dimensions

    m, n = np_ptid.shape

    # Build grid cell index lookup

    cell_indexes = {}
    for i in range(m):
        for j in range(n):
            if np_ptid[i, j] > 0:
                cell_indexes[np_ptid[i, j]] = (i, j)

    # Counter for number of boundaries
    
    n_boundaries = 0

    # Outflows

    for i_count, i_pointer in enumerate(np.where(pointer_labels < 0)[0]):
        segnum = pointer[i_pointer, 0]
        bndnum = pointer[i_pointer, 1]
        buff += "Outflow_%i\n" % (i_count + 1)
        buff += "1\n"
        n_boundaries += 1

        # Find a cell edge with no active neighbour

        i, j = cell_indexes[segnum]
        if i == 0 or np_ptid[i - 1, j] < 0:
            # first row or upper neighbour inactive: use upper edge
            point_a = xxul[i, j], yyul[i, j]
            point_b = xxlr[i, j], yyul[i, j]
        elif j == 0 or np_ptid[i, j - 1] < 0:
            # first column or left neighbour inactive: use left edge
            point_a = xxul[i, j], yylr[i, j]
            point_b = xxul[i, j], yyul[i, j]
        elif i == m - 1 or np_ptid[i + 1, j] < 0:
            # last row or lower neighbour inactive: use lower edge
            point_a = xxul[i, j], yylr[i, j]
            point_b = xxlr[i, j], yylr[i, j]
        elif j == n - 1 or np_ptid[i, j + 1]:
            # last column or right neighbour inactive: use right edge
            point_a = xxlr[i, j], yyul[i, j]
            point_b = xxlr[i, j], yylr[i, j]
        else:
            # no inactive neighbour: use upper left corner
            point_a = xxul[i, j], yyul[i, j]
            point_b = point_a

        buff += "%i %e %e %e %e\n" % (
            bndnum,
            point_a[0],
            point_a[1],
            point_b[0],
            point_b[1],
        )

    
    flow_labelsf = flow_labels.flatten()
    for bnd_id in boundflow_labels:
        bnd_pointer = pointer[flow_labelsf == bnd_id,:]
        buff += "%s\n" % (bnd_id)
        buff += "%i\n" % (len(bnd_pointer[:, 0]))
        n_boundaries += 1
        for i_pointer in np.arange(0, len(bnd_pointer[:, 0])):
            if bnd_pointer[i_pointer, 0] < 0.0:
                segnum = bnd_pointer[i_pointer, 1]
                bndnum = bnd_pointer[i_pointer, 0]
            else:
                segnum = bnd_pointer[i_pointer, 0]
                bndnum = bnd_pointer[i_pointer, 1]
            # Compute center coordinates of cell
            i, j = cell_indexes[segnum]
            x = (xxul[i, j] + xxlr[i, j]) * 0.5
            y = (yyul[i, j] + yylr[i, j]) * 0.5
            buff += "%i %e %e %e %e\n" % (bndnum, x, y, x, y)
                
            
    # Write file
    f = open(fname + ".bnd", "w")
    f.write("%i\n" % n_boundaries)
    f.write(buff)
    f.close()
    
def dw_WriteAttributesFile(fname, noseg):
    """
    Generates a Delwaq atributes (*.atr) file.

    input:
    - fname : file name to write to
    - noseg : number of delwaq segments
    """
    line_length = 100
    n_lines = noseg // line_length
    remaining_length = noseg % line_length

    buff = ""
    buff += "         ; DELWAQ_COMPLETE_ATTRIBUTES\n"
    buff += "    2    ; two blocks with input\n"
    buff += "    1    ; number of attributes, they are :\n"
    buff += "    1    ;  '1' is active '0' is no\n"
    buff += "    1    ; data follows in this fil\n"
    buff += "    1    ; all data is given without defaults\n"
    buff += ";    layer:            1\n"
    for iline in range(n_lines):
        buff += " ".join(["1" for _ in range(line_length)])
        buff += "\n"
    buff += " ".join(["1" for _ in range(remaining_length)])
    buff += "\n"

    buff += "    1    ; number of attributes, they are :\n"
    buff += "    2    ;  '1' has surface '3' has bottom\n"
    buff += "         ;  '0' has both    '2' has none\n"
    buff += "    1    ; data follows in this file\n"
    buff += "    1    ; all data is given without defaults\n"
    buff += ";    layer:            1\n"
    for iline in range(n_lines):
        buff += " ".join(["0" for _ in range(line_length)])
        buff += "\n"
    buff += " ".join(["0" for _ in range(remaining_length)])
    buff += "\n"

    buff += "    0    ; no time dependent attributes\n"
    f = open(fname, "w")
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
        return "{:04}{:02}{:02}{:02}{:02}{:02}".format(
            dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second
        )

    def timedelta2str(td):
        return "{:04}{:02}{:02}{}".format(
            0, 0, td.days, t.strftime("%H%M%S", t.gmtime(td.seconds))
        )

    buff = ""
    buff += "task      full-coupling\n"
    buff += "geometry  unstructured\n"
    buff += "horizontal-aggregation       no\n"
    buff += "minimum-vert-diffusion-used  no\n"
    buff += "vertical-diffusion           calculated\n"
    buff += "description\n"
    buff += "'%-60s'\n" % "Generated by Wflow"
    buff += "'%s'\n" % (" " * 60)
    buff += "'%s'\n" % (" " * 60)
    buff += "end-description\n"
    buff += "reference-time           '%s'\n" % (datetime2str(d["tref"]))
    buff += "hydrodynamic-start-time  '%s'\n" % (datetime2str(d["tstart"]))
    buff += "hydrodynamic-stop-time   '%s'\n" % (datetime2str(d["tstop"]))
    buff += "hydrodynamic-timestep    '%s'\n" % (timedelta2str(d["tstep"]))
    buff += "conversion-ref-time      '%s'\n" % (datetime2str(d["tref"]))
    buff += "conversion-start-time    '%s'\n" % (datetime2str(d["tstart"]))
    buff += "conversion-stop-time     '%s'\n" % (datetime2str(d["tstop"]))
    buff += "conversion-timestep      '%s'\n" % (timedelta2str(d["tstep"]))
    buff += "grid-cells-first-direction              %7i\n" % d["noseg"]
    buff += "grid-cells-second-direction             %7i\n" % 1
    buff += "number-hydrodynamic-layers              %7i\n" % 1
    buff += "number-horizontal-exchanges             %7i\n" % d["noqh"]
    buff += "number-vertical-exchanges               %7i\n" % d["noqv"]
    buff += "number-water-quality-segments-per-layer %7i\n" % d["nosegh"]
    buff += "number-water-quality-layers      1\n"
    buff += "hydrodynamic-file        none\n"
    buff += "aggregation-file         none\n"
    buff += "boundaries-file          '%s.bnd'\n" % d["runid"]
    buff += "waqgeom-file             '%s_waqgeom.nc'\n" % d["runid"]
    buff += "volumes-file             '%s.vol'\n" % d["runid"]
    buff += "areas-file               '%s.are'\n" % d["runid"]
    buff += "flows-file               '%s.flo'\n" % d["runid"]
    buff += "pointers-file            '%s.poi'\n" % d["runid"]
    buff += "lengths-file             '%s.len'\n" % d["runid"]
    buff += "salinity-file            none\n"
    buff += "temperature-file         none\n"
    buff += "vert-diffusion-file      none\n"
    buff += "horizontal-surfaces-file '%s.srf'\n" % d["runid"]
    buff += "depths-file              none\n"
    buff += "discharges-file          none\n"
    buff += "chezy-coefficients-file  none\n"
    buff += "shear-stresses-file      none\n"
    buff += "walking-discharges-file  none\n"
    buff += "attributes-file          '%s.atr'\n" % d["runid"]
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
    f = open(fname, "w")
    f.write(buff)
    f.close()

########################################################################
## Functions for dynamic files                                         #
########################################################################
    
def dw_WriteSegmentOrExchangeData(ttime, fname, datablock, boundids, WriteAscii=True):
    """
    Writes a timestep to a segment/exchange data file (appends to an existing
    file or creates a new one). 
        
    Input:
        - time - time for this timestep  
        - fname - File path of the segment/exchange data file</param>
        - datablock - array with data
        - boundids to write more than 1 block
        - WriteAscii - set to 1 to also make an ascii checkfile
        
    """
    #Supress potential NaN values to avoid error (replaced by -1.0)
    datablock[np.isnan(datablock)] = -1.0
    # Convert the array to a 32 bit float
    totareas = datablock
    for i in range(boundids - 1):
        totareas = np.vstack((totareas, datablock))

    artow = np.array(totareas, dtype=np.float32).copy()
    timear = np.array(ttime, dtype=np.int32)
            
    
    if os.path.isfile(fname):  # append to existing file
        fp = open(fname, "ab")
        tstr = timear.tostring() + artow.tostring()
        fp.write(tstr)
        if WriteAscii:
            fpa = open(fname + ".asc", "a")
            timear.tofile(fpa, format="%d\t", sep=":")
            artow.tofile(fpa, format="%10.8f", sep="\t")
            fpa.write("\n")
    else:
        fp = open(fname, "wb")
        tstr = timear.tostring() + artow.tostring()
        fp.write(tstr)
        if WriteAscii:
            fpa = open(fname + ".asc", "w")
            timear.tofile(fpa, format="%d\t", sep=":")
            artow.tofile(fpa, format="%10.8f", sep="\t")
            fpa.write("\n")

    fp.close()
    if WriteAscii:
        fpa.close()
        
def dw_Write_Times(dwdir, T0, timeSteps, timeStepSec, configfile):
    """
    Writes B1_timestamp.inc, B2_outputtimes.inc, (B2_sysclock.inc), B2_timers.inc, B10_simtime.inc.
    Assumes daily timesteps for now!
    """
    # B1_timestamp.inc
    exfile = open(dwdir + "/B1_timestamp.inc", "w")
    print(
        "'T0: " + T0.strftime("%Y.%m.%d %H:%M:%S") + "  (scu=       1s)'", file=exfile
    )
    exfile.close()

    # B2_outputtimes.inc
    timeRange = timedelta(seconds=timeStepSec * timeSteps)

    days = int(timeStepSec / 86400)
    hours = int(timeStepSec / 3600)
    minutes = int(timeStepSec / 60)
    seconds = int(timeStepSec - minutes * 60)
    minutes -= hours * 60
    hours -= days * 24
    timestepstring = "  %03d%02d%02d%02d" % (days, hours, minutes, seconds)

    exfile = open(dwdir + "/B2_outputtimes.inc", "w")
    etime = T0 + timeRange
    print(
        "  "
        + T0.strftime("%Y/%m/%d-%H:%M:%S")
        + "  "
        + etime.strftime("%Y/%m/%d-%H:%M:%S")
        + timestepstring
        + " ; mon/bal",
        file=exfile,
    )
    print(
        "  "
        + T0.strftime("%Y/%m/%d-%H:%M:%S")
        + "  "
        + etime.strftime("%Y/%m/%d-%H:%M:%S")
        + timestepstring
        + " ; map",
        file=exfile,
    )
    print(
        "  "
        + T0.strftime("%Y/%m/%d-%H:%M:%S")
        + "  "
        + etime.strftime("%Y/%m/%d-%H:%M:%S")
        + timestepstring
        + " ; his",
        file=exfile,
    )
    exfile.close()

    # B2_timers.inc
    exfile = open(dwdir + "/B2_timers.inc", "w")
    print("  " + T0.strftime("%Y/%m/%d-%H:%M:%S") + " ; start time", file=exfile)
    print("  " + etime.strftime("%Y/%m/%d-%H:%M:%S") + " ; stop time", file=exfile)
    print("  0 ; timestep constant", file=exfile)
    print("; dddhhmmss format for timestep", file=exfile)
    print(timestepstring + " ; timestep", file=exfile)
    exfile.close()

    # B2_sysclock.inc
    exfile = open(dwdir + "/B2_sysclock.inc", "w")
    print("%7d 'DDHHMMSS' 'DDHHMMSS'  ; system clock" % timeStepSec, file=exfile)
    exfile.close()
    
#    # B10_simtime.inc
#    exfile = open(dwdir + "/B10_simtime.inc", "w")
#    print("period 'sim'", file=exfile)
#    print("  suffix  'sim'", file=exfile)
#    print("  start-time  '" + T0.strftime("%Y/%m/%d-%H:%M:%S") + "'", file=exfile)
#    print("  stop-time  '" + etime.strftime("%Y/%m/%d-%H:%M:%S") + "'", file=exfile)
#    print("end-period", file=exfile)
#    exfile.close()
    
    #Statistics timers (block 10)
    secname = 'outputstat_'
    secnr = 0
    thissection = secname + str(secnr)
    toprint = configsection(configfile, thissection)
    
    while len(toprint) > 0:
        statname = configget(
            configfile, thissection, "name", "sim"
        )
        startperiod = configget(
            configfile, thissection, "start-period", T0.strftime("%Y/%m/%d-%H:%M:%S")
        )
        endperiod = configget(
            configfile, thissection, "end-period", etime.strftime("%Y/%m/%d-%H:%M:%S")
        )
        operation = configget(
            configfile, thissection, "output-operation", "STADSC"
        )
        
        
        exfile = open(dwdir + '/B10_' + statname + 'time.inc', "w")
        print("period '" + statname + "'", file=exfile)
        print("  suffix  '"+statname+"'", file=exfile)
        print("  start-time  '" + startperiod + "'", file=exfile)
        print("  stop-time  '" + endperiod + "'", file=exfile)
        print("end-period", file=exfile)
        print("", file=exfile)
                
        print("output-operation '" + operation + "'", file=exfile)
        
        for option in toprint:
            if (
                "name" not in option
                and "start-period" not in option
                and "end-period" not in option
                and "output-operation" not in option
            ):
                optionname = configget(configfile,thissection, option, "")
                lineprint = "  " + option +"  '" + optionname + "'"
                print(lineprint, file = exfile)

        print("  suffix  ''", file=exfile)
        print("end-output-operation", file=exfile)
        exfile.close()
        
        secnr = secnr + 1
        thissection = secname + str(secnr)
        toprint = configsection(configfile, thissection)
        

    
def dw_WriteCSVdata(timestepsecs, fname, varcsv, boundids, WriteAscii=False):
    """ 
    Writes a segment/exchange data file from aggregated csv tables outputs from wflow_sbm.
    Contrary to dw_WriteSegmentOrExchangeData, writes all the timesteps at the same time.
    
    Input:
        - timestepsecs: model timestep in seconds
        - fname: File path of the segment/exchange data file
        - varcsv: dataframe with exchange data
        - boundids: to write more than one block
        - WriteAscii: set to 1 to also make an ascii checkfile
    """    
    #Repeat the flow block for each compartment (boundids, needed for Emission)
    totareas = varcsv
    for i in range(boundids - 1):
        totareas = pd.concat([totareas, varcsv], axis = 1)
    # Write binary file
    f = open(fname, "wb")
    for i in range(totareas.shape[0]):
        tstep = np.array(i*timestepsecs, dtype = np.int32)
        vararr = np.array(totareas.loc[i].values, dtype=np.float32)
        f.write(tstep.tostring() + vararr.tostring())
    f.close()
    
    if WriteAscii:
        # Write ASCII file
        exfile = open(fname + ".asc", "w")
        for i in range(totareas.shape[0]):
            tstep = np.array(i*timestepsecs, dtype = np.int32)
            vararr = np.array(totareas.loc[i].values, dtype=np.float32)
            tstep.tofile(exfile, format="%d\t", sep=":")
            vararr.tofile(exfile, format="%10.8f", sep="\t")
            exfile.write("\n")
        exfile.close()
    
def _readTS(name, ts):
    """
    Read a pcraster map for a timestep without using the dynamic framework 
    """
    mname = os.path.basename(name)
    #  now generate timestep
    tsje = "%0.11d" % ts
    ff = mname + tsje[len(mname) :]
    ff = ff[:8] + "." + ff[8:]
    name = os.path.dirname(name) + "/" + ff
    mapje = pcr.readmap(name)

    return mapje

def read_timestep(nc, var, timestep, logger, caseId, runId):
    """
    Returns a map of the given variable at the given timestep.
    """
    if nc is not None:
        pcrmap, succes = nc.gettimestep(timestep, logger, var=var)
        assert succes
        return pcrmap
    else:
        return _readTS(caseId + "/" + runId + "/outmaps/" + var, timestep)
    

########################################################################
## Functions for emission data                                         #
########################################################################

def dw_WriteEmiData(caseId, outdir, ptid, cells, modelmap, emiData, comp, Aggregation, logger):
    """
    Convert PCRaster emission maps into include files.
    """
    
    for l in range(len(emiData)):
        #Read the pcr map and convert it to numpy
        emimap = pcr.readmap((caseId + "/" + emiData.Fileloc[l]))
        
        #Transform data for aggregation (default is average)
        if Aggregation:
            if emiData.AggType[l] == "total":
                emimap = pcr.areatotal(emimap, modelmap, Aggregation)
            else: 
                emimap = pcr.areaaverage(emimap, modelmap, Aggregation)
        np_emimap = dw_pcrToDataBlock(emimap, ptid, Aggregation)
        #Transforms emi data to zeros instead of NaN
        np_emimap[np.isnan(np_emimap)] = 0.0
        
        #Create emi data block for binary data
        emi_block = []
        for i in np.arange(0, len(comp)):
            if comp.ID[i] == emiData.To[l]:
                emi_block = np.append(emi_block,np_emimap)
            else:
                emi_block = np.append(emi_block, np.repeat(0, len(cells)))
        
        logger.info("Writing " + emiData.Name[l] + " data")
        dw_WriteSegmentOrExchangeData(
            0, outdir + emiData.Name[l] +".dat", emi_block, 1, WriteAscii=False
        )
        
        #For ASCII, maximum length per line is 1000, need to split emission data
        #Need to reshape the emi data to fit
        #Calculate max size of characters (+1 for space)
        MaxChar = np.amax([len(x) for x in np_emimap.astype(str)])+1
        #Number of values per line is floor(1000/MaxChar)
        NPRow = math.floor(1000/MaxChar)
        #Total number of values to write
        NTOT = len(np_emimap)
        #Reshape if necessary
        if NTOT > NPRow:
            #Number of complete rows
            NCRow = int(len(np_emimap)/NPRow)
            np_emimap_1 = np_emimap[0:NCRow*NPRow].reshape(NCRow, NPRow)  
            np_emimap_2 = np_emimap[NCRow*NPRow:NTOT]
            np_emimap_2 = np_emimap_2.reshape(1, len(np_emimap_2))
            
        #Write the ASCII file
        exfile = open((outdir + emiData.Name[l] +".inc"), "w")
        if NTOT > NPRow:
            np.savetxt(exfile, np_emimap_1.astype(str), fmt="%s")
            np.savetxt(exfile, np_emimap_2.astype(str), fmt="%s")
        else:
            np.savetxt(exfile, np.emimap.astype(str), fmt="%s")        
        #Put zeros for the remaining compartments
        if len(comp) > 1:                
            extra_data = str(int(len(cells)*(len(comp)-1)*4)) + "*0.0"
            print(extra_data, file=exfile)
        exfile.close()        
            


########################################################################
## Main function                                                       #
########################################################################

def main(argv=None):
    """
    Perform command line execution of the model.
    """
    
    from dateutil import parser
    
    ########################################################################
    ## Default settings                                                    #
    ########################################################################
    caseId = "default_sbm"
    global multpars
    runId = "run_default"
    dwdir = "default_waq"
    configfile = "wflow_emwaq.ini"
    _lastTimeStep = 0
    _firstTimeStep = 0
    LogFileName = "wflow.log"

    runinfoFile = "runinfo.xml"
    timeSteps = 1
    timestepsecs = 86400
    
    # Default model options
    WriteAscii = False
    WriteDynamic = True
    WriteStructure = True
    Fraction = False
    Emission = False
    CalcTotflow = False
    Fews = False
    Aggregation = False


    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return
    
    ########################################################################
    ## Process command-line for directories                                #
    ########################################################################                                           
    try:
        opts, args = getopt.getopt(argv, "C:R:D:c:n:s:T:S:iuyfeFa")
    except getopt.error as msg:
        pcrut.usage(msg)

    for o, a in opts:
        if o == "-C":
            caseId = a
        if o == "-R":
            runId = a
        if o == "-D":
            dwdir = a
        if o == "-c":
            configfile = a

    
    #Initialize waq output directories
    dw_CreateDwRun(dwdir)
    
    logger = pcrut.setlogger(dwdir + "/debug/wflow_emwaq.log", "wflow_emwaq")
        
    ########################################################################
    ## Process ini file                                                    #
    ########################################################################
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(caseId + "/" + configfile)
    
    #Model options
    write_ascii = int(configget(config, "model", "write_ascii", "0"))
    if write_ascii == 1: WriteAscii = True
    write_structure = int(configget(config, "model", "write_structure", "1"))
    if write_structure == 0: WriteStructure = False
    write_dynamic = int(configget(config, "model", "write_dynamic", "1"))
    if write_dynamic == 0: WriteDynamic = False
    fraction = int(configget(config, "model", "fraction", "0"))
    if fraction == 1: Fraction = True
    emission = int(configget(config, "model", "emission", "0"))
#    logger.info("Emission =" + str(emission))
    if emission == 1: Emission = True
    calc_totflow = int(configget(config, "model", "calc_totflow", "0"))
    if calc_totflow == 1: CalcTotflow = True
    fews = int(configget(config, "model", "fews", "0"))
    if fews == 1: Fews = True
    aggregation = int(configget(config, "model", "aggregation", "0"))
    if aggregation == 1: Aggregation = True
    
    #NetCDF or PCRaster or CSV
    inputType = configget(config, "model", "input_type", "netcdf")
    if inputType == "netcdf":
        netcdffile = configget(config, "model", "netcdfinput", "wflow2waq.nc")

    #Template ini file
    template_ini_path = configget(config, "model", "template_ini", "None")
    
    #Monitoring options
    mon_points = configget(config, "model", "mon_points", "gauges")
    mon_areas = configget(config, "model", "mon_areas", "subcatch")
    
    #Boundaries ID creation
    bd_id = configget(config, "model", "bd_id", "boundaries")
    
    #Layout option
    sizeinmetres = int(configget(config, "layout", "sizeinmetres", "0"))
    
    #Time options
    timestepsecs = int(configget(config, "run", "timestepsecs", str(timestepsecs)))

    st = configget(config, "run", "starttime", "None")
    runlengthdetermination = configget(config, "run", "runlengthdetermination", "steps")
    
    if st == "None":  # try from the runinfo file
        rinfo_str = configget(config, "run", "runinfo", "None")
        if rinfo_str != "None":
            T0 = wflow_adapt.getStartTimefromRuninfo(caseId + "/" + rinfo_str)
            datetimeend = wflow_adapt.getEndTimefromRuninfo(caseId + "/" + rinfo_str)
        else:
            logger.error(
                "Not enough information in the [run] section. Need start and end time or a runinfo.xml file...."
            )
            sys.exit(1)
    else:
        T0 = parser.parse(st)
        ed = configget(config, "run", "endtime", "None")
        if ed != "None":
            datetimeend = parser.parse(ed)
        else:
            logger.error("No end time given with start time: [run] endtime = " + ed)
            sys.exit(1)

    if runlengthdetermination == "steps":
        runStateTime = T0 #- datetime.timedelta(seconds=timestepsecs)
    else:
        runStateTime = T0 + datetime.timedelta(seconds=timestepsecs)

    timeSteps = (
        calendar.timegm(datetimeend.utctimetuple())
        - calendar.timegm(runStateTime.utctimetuple())
    ) / timestepsecs
    
    firstTimeStep = 0
    
    logger.info("T0 of wflow run: " + str(T0))
    
    #Reading basemaps
    logger.info("Reading basemaps")

    wflow_subcatch = (
        caseId + "/"
        + configget(config, "model", "wflow_subcatch", "/staticmaps/wflow_subcatch.map")
    )
    pcr.setclone(wflow_subcatch)
    amap = pcr.scalar(pcr.readmap(wflow_subcatch))
    modelmap = pcr.readmap(wflow_subcatch)
    ldd = pcr.readmap(
        caseId + "/"
        + configget(config, "model", "wflow_ldd", "/staticmaps/wflow_ldd.map")
    )
    gauges = pcr.readmap(
        caseId + "/"
        + configget(config, "model", "wflow_gauges", "/staticmaps/wflow_gauges.map")
    )
#    reallength = pcr.readmap(
#        caseId + "/" + runId + "/" 
#        + configget(config, "model", "wflow_reallength", "/outsum/reallength.map")
#    )
    xl, yl, reallength = pcrut.detRealCellLength(
            amap, sizeinmetres
        )
    cellsize = float(pcr.pcr2numpy(reallength, np.NaN)[0,0])
    logger.info("Cellsize model: " + str(cellsize))
    
    PathFrac = pcr.readmap(
        caseId + "/" + runId + "/" 
        + configget(config, "model", "wflow_pathfrac", "/outsum/PathFrac.map")
    )
    WaterFrac = pcr.readmap(
        caseId + "/" + runId + "/" 
        + configget(config, "model", "wflow_waterfrac", "/outsum/WaterFrac.map")
    )
    
    # Limit areas map to modelmap (subcatchments)
    amap = pcr.ifthen(modelmap > 0, amap)
    ldd = pcr.ifthen(amap > 0, ldd)
    gmap = pcr.ifthen(modelmap > 0, pcr.scalar(gauges))
    pcr.report(amap, dwdir + "/debug/area.map")
    pcr.report(ldd, dwdir + "/debug/ldd.map")
    pcr.report(modelmap, dwdir + "/debug/modelmap.map")
    pcr.report(gmap, dwdir + "/debug/gmap.map")

    thecells = pcr.pcr2numpy(modelmap, np.NaN).flatten()
    nrcells = len(thecells)
    nractcells = len(thecells[np.isfinite(thecells)])

    logger.info("Total number gridcells (including inactive): " + str(nrcells))
    logger.info("Total number of used gridcells: " + str(nractcells))
    
    # find all upstream cells (these must be set negative)
    upbound = pcr.upstream(ldd, 1.0)
    upbound = pcr.ifthen(upbound == 0, upbound)
    upar = pcr.pcr2numpy(pcr.scalar(upbound), np.NaN).flatten()
    logger.info(
        "Number of upstream cells (without upstream connection): "
        + str(len(upar[np.isfinite(upar)]))
    )
    pcr.report(upbound, dwdir + "/debug/upbound.map")
    
    
    ########################################################################
    ## Read compartments/fluxes/boundaries inputs                          #
    ########################################################################
    sepcsv = configget(config, "inputcsv", "sepcsv", ",")
    compFile = caseId + "/" + configget(config, "inputcsv", "compartments", "/csv/compartments.csv")
    comp = pd.read_csv(compFile, sep = sepcsv, header=0)
    boundFile = caseId + "/" + configget(config, "inputcsv", "boundaries", "/csv/boundaries.csv")
    bound = pd.read_csv(boundFile, sep = sepcsv, header=0)
    fluxFile = caseId + "/" + configget(config, "inputcsv", "fluxes", "/csv/fluxes.csv")
    flux = pd.read_csv(fluxFile, sep = sepcsv, header=0)
    emiFile = caseId + "/" + configget(config, "inputcsv", "emissions", "/csv/emissions.csv")
    if os.path.exists(emiFile):
        includeEmi = True
        emiData = pd.read_csv(emiFile, sep = sepcsv, header=0)
    else:
        includeEmi = False
    
    logger.info("Nb of compartments: " + str(len(comp)))
    logger.info("Nb of boundaries: " + str(len(bound)))
    logger.info("Nb of fluxes: " + str(len(flux)))
    
    
    ########################################################################
    ## Process command-line options                                        #
    ######################################################################## 
    for o, a in opts: 
#        if o == "-n":
#            netcdffile = a
#        if o == "-s":
#            timestepsecs = int(a)
#        if o == "-T":
#            _lastTimeStep = int(a)
#        if o == "-S":
#            _fisrtTimeStep = int(a)
        if o == "-i":
            WriteAscii = True
        if o == "-u":
            WriteStructure = False
        if o == "-y":
            WriteDynamic = False
        if o == "-f":
            Fraction = True
        if o == "-e":
            Emission = True
        if o == "-F":
            Fews = True
        if o == "-a":
            Aggregation = True

    #Settings warning or errors
    if Aggregation and inputType != "csv":
        logger.error("Aggregation of wflow_sbm results is not implemented with netcdf or PCRaster maps, use csv files instead.")
    if inputType == "csv" and not Aggregation:
        logger.warning("CSV input file detected but aggregation option is not set. Make sure CSV inputs are saved for all cells and not aggregated cells.")
    if Aggregation and (mon_points == "gauges" or mon_areas == "subcatch"):
        logger.error("Chosen monitoring option(s) are not available with aggregation option")
    
    
    ########################################################################
    ## Create all waq inputs from wflow_sbm                                #
    ########################################################################
    global pointer
    
    ######## Structure files ########
    #Get pointer and boundaries from ldd, subcatch and defined boundaries
    logger.info("Making pointer...")
    ptid, pointer, pointer_labels, flow_labels, cells, segments, seg_labels, flow_var = dw_mkDelwaqPointers(
        ldd, amap, gmap, comp, flux, Emission, Aggregation, bd_id
    )
    
    if WriteStructure:    
        # Write id maps to debug area
        pcr.report(ptid, dwdir + "/debug/ptid.map")
        logger.info("Number of segments: " + str(len(segments.flatten())))
        logger.info(
            "Number of internal flows: " + str(len(pointer_labels[pointer_labels == 0]))
        )
        
        #Number of segments and of exchanges
        NOSQ = segments.shape[0]
        NOQ = pointer.shape[0]
        
        #Write structure files related to pointer
        logger.info("Writing structure files...")
        dw_WriteNrSegments(dwdir + "/includes_deltashell/B3_nrofseg.inc", NOSQ)
        dw_WriteAttributes(dwdir + "/includes_deltashell/B3_attributes.inc", comp, cells)
        dw_WritePointer(dwdir + "/includes_deltashell/B4_pointer.inc", pointer)
        #Write the number of exchanges
        dw_WriteNrExChanges(dwdir + "/includes_deltashell/B4_nrofexch.inc", NOQ)
        boundflow_labels = dw_WriteBoundlist(
                dwdir + "/includes_deltashell/B5_boundlist.inc", flow_labels, pointer_labels, bd_id
        )
        
        #Convert emssions map data to inc
        if includeEmi:
            dw_WriteEmiData(caseId, dwdir + "/includes_flow/", ptid, cells, modelmap, emiData, comp, Aggregation, logger)
        
        #Fraction files
        if Fraction:
            dw_WriteBoundData(
                dwdir + "/includes_deltashell/B5_bounddata.inc", boundflow_labels
            )
            dw_WriteInitials(dwdir + "/includes_deltashell/B8_initials.inc", boundflow_labels)
            dw_Write_Substances(
                dwdir + "/includes_deltashell/B1_sublist.inc", boundflow_labels
            )
            if Fews:
                dw_WriteFractionParameterIds(dwdir + "/includes_fews/fraction_parameterids.csv", boundflow_labels)
        
        #Monitoring files
        logger.info("Writing monitoring files...")        
        NOPT, NOAR, stations, stations_balance = dw_mkDelwaqMonitoring(mon_points, mon_areas, gmap, amap,
                                                     comp, cells, segments, seg_labels)
        dw_WriteNrMonAreas(dwdir + "/includes_deltashell/B2_nrofmon.inc", NOPT, NOAR, mon_points, mon_areas)
        if NOPT != 0:
            dw_WriteStations(dwdir + "/includes_deltashell/B2_stations.inc", stations)
            dw_WriteStations(dwdir + "/includes_deltashell/B2_stations-balance.inc", stations_balance)
        if NOAR != 0:
            dw_WriteMonAreas(dwdir + "/includes_deltashell/B2_monareas.inc", NOAR, 
                             mon_areas, amap, comp, segments, seg_labels, cells)
        
        #Emission files
        if Emission:
            logger.info("Writing additional files for emission")
            dw_WriteNrSegL(dwdir + "/includes_deltashell/nrofsegl.inc", NOSQ, len(comp))
            dw_WriteGeometry(dwdir + "/includes_deltashell/", cells, amap, reallength, PathFrac, WaterFrac, Aggregation, WriteAscii, len(comp))
            dw_WriteHydDef(dwdir + "/includes_deltashell/hydrology.inc", cells, flux)
            
    
    ######## Static data ########
    #For D-WAQ
    #if not Emission:
    #Surface data [m^2]
    #For surface runoff, the surface is adjusted in river cells to the real river size.
    #For other compartments, surface should be cell size.
    if WriteDynamic:
        internalflowwidth = pcr.readmap(caseId + "/" + runId + "/outsum/Bw.map")
        internalflowlength = pcr.readmap(caseId + "/" + runId + "/outsum/DCL.map")
        surfaceWater_map = internalflowwidth * internalflowlength
        if Aggregation:
            surfaceWater_map = pcr.areatotal(surfaceWater_map, modelmap)
        surfaceWater_block = dw_pcrToDataBlock(surfaceWater_map, amap, Aggregation)
    
        realarea = reallength * reallength
        if Aggregation:
            realarea = pcr.areatotal(realarea, modelmap)
        realarea_block = dw_pcrToDataBlock(realarea, amap, Aggregation)
        
        surface_block = []
        for i in np.arange(0, len(comp)):
            if comp.At2[i] == 0:
                surface_block = np.append(surface_block, surfaceWater_block)
            else:
                surface_block = np.append(surface_block, realarea_block)
         
        if not Emission:           
            logger.info("Writing surface.dat. Nr of points: " + str(np.size(surface_block)))
            dw_WriteSegmentOrExchangeData(
                0, dwdir + "/includes_flow/surface.dat", surface_block, 1, WriteAscii
            )
                    
            #Manning
            #Manning for Surface Water, zero for other compartments
            Manning_block = []
            for i in np.arange(0, len(comp)):
                if comp.At2[i] == 0:
                    Manning_map = pcr.readmap(caseId + "/" + runId + "/outsum/N.map")
                    if Aggregation:
                        Manning_map = pcr.areaaverage(Manning_map, modelmap)
                    Manning_block = np.append(Manning_block, dw_pcrToDataBlock(Manning_map, amap, Aggregation))
                else:
                    Manning_block = np.append(Manning_block, np.repeat(0, len(cells)))
            logger.info("Writing Manning n coefficient")
            dw_WriteSegmentOrExchangeData(
                0, dwdir + "/includes_flow/manning.dat", Manning_block, 1, WriteAscii
            )
        
    #FEWS files
    if Fews:
        logger.info("Writing additional files for FEWS")
        np.save(dwdir + "/debug/pointer.npy", pointer)
        np.save(dwdir + "/debug/segments.npy", segments)
        dw_Write_B2_outlocs(dwdir + "/includes_deltashell/B2_outlocs.inc", gauges, ptid)
        
        # write static data for hyd-file set
        comroot = os.sep.join([dwdir, "com", runId])
        mmax, nmax = dw_GetGridDimensions(reallength)
        logger.info("mmax, nmax: " + str(mmax) +","+ str(nmax))
        dw_WritePointer(comroot + ".poi", pointer, binary=True)
        if not Emission and WriteDynamic:
            dw_WriteSurfaceFile(comroot + ".srf", surface_block)
            #dw_WriteSegmentOrExchangeData(0, comroot + ".len", length_block, 1, WriteAscii)
        logger.info("Writing waq geometry file")
        dw_WriteWaqGeom(comroot, ptid, ldd, Fews)
        logger.info("Writing fews boundary file")
        dw_WriteBndFile(comroot, ptid, pointer, pointer_labels, flow_labels)
       
    
    ######## Dynamic data ########
    if WriteDynamic:
        #Set up netcdf or mapstack of volumes and fluxes
        var_list = np.append(comp.mapstack.values, flux.mapstack.values)
        #Remove compartments with no volumes (mapstack = ZeroMap)
        var_list = var_list[var_list != "ZeroMap"]
        #Separate the variables that are sum of wflow fluxes
        var_sublist = var_list[np.char.find(var_list.astype(np.str), "+") != -1]
        if len(var_sublist) > 0:
            var_list = var_list[np.char.find(var_list.astype(np.str), "+") == -1]
            for v in range(len(var_sublist)):
                wflowVars = var_sublist[v].split("+")
                var_list = np.append(var_list, wflowVars)
        if inputType == "netcdf":
            pcr.setglobaloption('coorcentre')
            netcdffile = caseId + "/" + runId + "/" + netcdffile
            nc = wf_netcdfio.netcdfinput(
                netcdffile, logger, var_list
            )
        else:
            nc = None
        
        # mask to filter out inactive segments
        zeroMap = 0.0 * pcr.scalar(amap)
        minVolumeMap = zeroMap + 0.0001
        
        #Start dynamic section
        ts = 1
    
        #Time inputs
        logger.info("Writing timers")
        #dw_Write_Times(dwdir + "/includes_deltashell/", T0, timeSteps - 1, timestepsecs)
        dw_Write_Times(dwdir + "/includes_deltashell/", T0, timeSteps, timestepsecs, config)
        
        logger.info("Writing dynamic data...")
        #Start loop over time steps to read and write netcdf or PCRaster maps
        if (inputType == "netcdf" or inputType == "map"):
            for i in range(firstTimeStep, (int(timeSteps) + 1) * int(timestepsecs), timestepsecs):
                logger.info("Timestep :" + str(ts) + "/" + str(int(timeSteps)+1))
                
                #Volumes of compartments (for D-WAQ)
                if not Emission:
                    #Wflow writes volumes at the end of the timestep and DWAQ needs them at the beginning -> offset
                    #In order to avoid zero volumes, a basic minimum value of 0.0001 m3 is added to all volumes
                    volume_block = []
                    for c in range(0, len(comp)):
                    
                        #Try fisrt time step from states else use zeromap (should be same config as given in wflow_sbm)
                        if i == firstTimeStep:
                            #In the states the unsaturated store is divided in several layers
                            if comp.wflowVar[c] == "UStoreDepth" :
                                #Get the number of layers
                                UStoreLayerThickness = configget(config, "model", "UStoreLayerThickness", "0")
                                if UStoreLayerThickness != "0":
                                    maxLayers = len(UStoreLayerThickness.split(",")) + 1
                                else:
                                    maxLayers = 1
                                volume_map = zeroMap*0.0
                                #Read and sum state maps of each layer of the unsaturated store
                                for USi in np.arange(0, maxLayers):
                                    volume_map_path = caseId + "/instate/UStoreLayerDepth_" + str(USi) + ".map"
                                    if os.path.exists(volume_map_path):
                                        volume_map = volume_map + pcr.readmap(volume_map_path)
                            
                            #If KinWaveVolume is used instead of WaterLevel:
                            elif (comp.wflowVar[c] == "KinWaveVolumeR+KinWaveVolumeL" or comp.wflowVar[c] == "KinWaveVolumeL+KinWaveVolumeR"):
                                volume_map_path1 = caseId + "/instate/" + "WaterLevelR" + ".map"
                                volume_map_path2 = caseId + "/instate/" + "WaterLevelL" + ".map"
                                if os.path.exists(volume_map_path):
                                    volume_map1 = pcr.readmap(volume_map_path1)
                                    volume_map2 = pcr.readmap(volume_map_path2)
                                    volume_map = (volume_map1 + volume_map2) * surfaceWater_map
                                else:
                                    volume_map = zeroMap*0.0
                            #Other compartments    
                            else:
                                volume_map = zeroMap*0.0
                                wflowVars = comp.wflowVar[c].split("+")
                                for wv in range(len(wflowVars)):
                                    volume_map_path = caseId + "/instate/" + wflowVars[wv] + ".map"
                                    if os.path.exists(volume_map_path):
                                        volume_map = volume_map + pcr.readmap(volume_map_path)
                    
                        #Other timesteps
                        else:
                            volume_map = zeroMap*0.0
                            wflowVars = comp.mapstack[c].split("+")
                            for wv in range(len(wflowVars)):
                                volume_map = volume_map + read_timestep(nc, wflowVars[wv], (ts - 1), logger, caseId, runId)    
                    
                        #If needed, convert volume in mm or m to m3
                        if comp.Unit[c] == "mm" or comp.Unit[c] == "m":
                            #Adjusted area for the surface water
                            if comp.At2[c] == 0:
                                volume_map = volume_map * surfaceWater_map 
                            #TODO: for soil, recalculate real volume (with bulk volume and porosity)
                            #Now, volume = equivalent water volume
                            else:
                                volume_map = volume_map * realarea
                            if comp.Unit[c] == "mm":
                                volume_map = volume_map / 1000
                        #Add minimum volume and aggregates the compartments
                        volume_map = volume_map + minVolumeMap
                        volume_block = np.append(volume_block, dw_pcrToDataBlock(volume_map, amap))
                
                    #logger.info("Writing volumes.dat. Nr of points: " + str(size(volume_block)))
                    dw_WriteSegmentOrExchangeData(
                        float(i), dwdir + "/includes_flow/volume.dat", volume_block, 1, WriteAscii
                    )
                
                #Flow data, exchange area and velocity
                flow_block = []
                area_block = []
                conv_block = []
                velocity_block = []
                
                #In D-WAQ, flow are the fluxes defined as in the pointer file
                #In D-Emission, flow are the fluxes needed in the hydrology file and for total flow
                if Emission:
                    hydflux = ~flux.EmPointerHyd.isin(["P"])
                    flow_var = flux.mapstack[hydflux].values.astype(str)
                    #Initialize total flow map
                    totalflow_map = zeroMap*0.0
                
                for f in np.arange(0, len(flow_var)):
                    #Flow
                    #logger.info("Reading flow:" + str(flow_var[f]))
                    #For the last timestep, volume is defined, flow is assumed to be equal to the previous timestep
                    flow_map = zeroMap*0.0
                    wflowVars = flow_var[f].split("+")
                    if i == (int(timeSteps) + 1) * int(timestepsecs):
                        for v in range(len(wflowVars)):
                            flow_map = flow_map + read_timestep(nc, wflowVars[v], (ts - 1), logger, caseId, runId)
                    else:
                        for v in range(len(wflowVars)):
                            flow_map = flow_map + read_timestep(nc, wflowVars[v], ts, logger, caseId, runId)
                    #For emission add the read flow to total flow if needed
                    #Add the read flow to flow_block if part of the hydrology.bin
                    #For D-WAQ, all read flows go in the flow_block
                    if Emission:
#                        totflowlist= array(["RunoffOpenWater", "InfiltExcessSoil", "InfiltExcessPath", 
#                                         "ExfiltFromUstore","ExfiltFromSat", "ActEvapOpenWater"])
                        addtotflow_f = flux.EmPointerHyd[flux.mapstack == flow_var[f]].values
                        if (bool(np.isin(addtotflow_f[0], np.array(["T","PT","HT","PHT"])))):
                            #logger.info("Adding " + str(flow_var[f]) + " to total flow")
                            totalflow_map = totalflow_map + flow_map
                        addflowblock_f = flux.EmPointerHyd[flux.mapstack == flow_var[f]].values
                        if (bool(np.isin(addflowblock_f[0], np.array(["H","PH","HT","PHT"])))):
                            #logger.info("Adding " + str(flow_var[f]) + " to flow block")
                            #If needed convert flow unit to fraction by dividing by the from comp volume
                            #Contrary to volume.dat for D-WAQ, here the volume needed is the one at the
                            #end of the timestep.
                            convfrac_f = flux.EmFraction[flux.mapstack == flow_var[f]].values
                            unit_f = flux.Unit[flux.mapstack == flow_var[f]].values
                            if convfrac_f[0] == 1:
                                #logger.info("Converting " + str(flow_var[f]) + " to fraction")
                                #Find and read the corresponding volume
                                fromcomp = flux.From[flux.mapstack == flow_var[f]].values
                                volume_var = comp.mapstack[comp.ID == fromcomp[0]].values
                                volume_map = zeroMap*0.0
                                wflowVars = volume_var[0].split("+")
                                for v in range(len(wflowVars)):
                                    volume_map = volume_map + read_timestep(nc, wflowVars[v], ts, logger, caseId, runId)
                                #Convert to fraction and correct for zero volumes
                                flow_map = pcr.ifthenelse(volume_map != 0.0, flow_map / volume_map, pcr.scalar(0.0))
                            #Else convert from mm to m3/s
                            if (convfrac_f[0] == 0 and  unit_f[0] == "mm"):
                                flow_map = flow_map * reallength * reallength / (1000 * timestepsecs)
                            flow_block = np.append(flow_block, dw_pcrToDataBlock(flow_map, amap))
                    else:
                        flow_block = np.append(flow_block, dw_pcrToDataBlock(flow_map, amap))
                    
                    #Exchange area and velocity for D-WAQ
                    #Also create a conversion vector for flow                 
                    if not Emission:
                        #For D-WAQ flow must be converted from mm to m3/s
                        unit_f = flux.Unit[flux.mapstack == flow_var[f]].values
                        if unit_f[0] == "mm":
                            conv = np.repeat(1, len(cells))
                        else:
                            conv = np.repeat(0, len(cells))                
                        conv_block = np.append(conv_block, conv)
                    
    
                        #Exchange Area [m^2]
                        #For lateral flows (only surface runoff and saturated store), the exchange area is derived from water level
                        wflowVar_f = flux.wflowVar[flux.mapstack == flow_var[f]].values
                        if (wflowVar_f[0] == "RiverRunoff+LandRunoff" or wflowVar_f[0] == "LandRunoff+RiverRunoff") :
                            if i == firstTimeStep:
                                level_map_path1 = caseId + "/instate/WaterLevelR.map"
                                level_map_path2 = caseId + "/instate/WaterLevelL.map"
                                if os.path.exists(level_map_path1):
                                    level_map = pcr.readmap(level_map_path1)
                                else:
                                    level_map = zeroMap*0.0
                                if os.path.exists(level_map_path2):
                                    level_map = level_map + pcr.readmap(level_map_path2)
                            else:
                                level_mapstack = comp.mapstack[comp.wflowVar == "WaterLevelR+WaterLevelL"].values
                                if len(level_mapstack) == 0:
                                    level_mapstack = comp.mapstack[comp.wflowVar == "WaterLevelL+WaterLevelR"].values
                                level_map = zeroMap*0.0
                                wflowVars = level_mapstack[0].split("+")
                                for v in range(len(wflowVars)):
                                    level_map = level_map + read_timestep(nc, wflowVars[v], (ts - 1), logger, caseId, runId)
                            area_map = level_map * pcr.readmap(caseId + "/" + runId + "/outsum/Bw.map")
                            #Spot lateral flow to derive velocity from area and flow block
                            velocity = np.repeat(1, len(cells))
                        #TODO: define real exchange area for soil (bulk density, porosity)
                        elif wflowVar_f[0] == "SatWaterFlux" :
                            if i == firstTimeStep:
                                level_map_path = caseId + "/instate/SatWaterDepth.map"
                                if os.path.exists(level_map_path):
                                    level_map = pcr.readmap(level_map_path)
                                else:
                                    level_map = zeroMap*0.0
                            else:
                                level_mapstack = comp.mapstack[comp.wflowVar == "SatWaterDepth"].values
                                level_map = read_timestep(nc, level_mapstack[0], (ts - 1), logger, caseId, runId)
                            area_map = level_map * reallength
                            velocity = np.repeat(2, len(cells))
                        else:
                            area_map = reallength * reallength
                            velocity = np.repeat(0, len(cells))
                    
                        area_block = np.append(area_block, dw_pcrToDataBlock(area_map, amap))
                        velocity_block = np.append(velocity_block, velocity)
                    
                    
                        #Convert flow from mm to m3/s
                        conv_block = conv_block == 1
                        flow_block[conv_block] = (flow_block[conv_block] / 1000) * area_block[conv_block] / timestepsecs
                    
                        #Calculate velocity
                        #Correct nul areas to avoid dision by 0 errors
                        area_corr = area_block
                        flow_corr = flow_block
                        area_corr[area_corr == 0] = 0.1
                        flow_corr[area_corr == 0] = 0.0
                        velocity = []
                        for c in np.arange(0, len(comp)):
                            if (comp.wflowVar[c] == "LandRunoff+RiverRunoff" or comp.wflowVar[c] == "RiverRunoff+LandRunoff"):
                                calcv = velocity_block == 1
                                velocity = np.append(velocity, flow_corr[calcv] / area_corr[calcv])
                            elif comp.wflowVar[c] == "SatWaterFlux":
                                calcv = velocity_block == 2
                                velocity = np.append(velocity, flow_corr[calcv] / area_corr[calcv])
                            else:
                                velocity = np.append(velocity, np.repeat(0, len(cells)))
                
                #end of the loop reading from the flux list
    
                #Write dynamic data
                if Emission:
                    #For emission, convert and add totalflow map to flow block before writting
                    totalflow_map = totalflow_map * reallength * reallength / (1000 * timestepsecs) #mm to m3/s
                    flow_block = np.append(flow_block, dw_pcrToDataBlock(totalflow_map, amap))
                    #logger.info("Writing flow.dat. Nr of points: " + str(size(flow_block)))
                    dw_WriteSegmentOrExchangeData(
                        float(i), dwdir + "/includes_flow/hydrology.bin", flow_block, len(comp), WriteAscii
                    )   
                else:
                    #logger.info("Writing flow.dat. Nr of points: " + str(size(flow_block)))
                    dw_WriteSegmentOrExchangeData(
                        float(i), dwdir + "/includes_flow/flow.dat", flow_block, 1, WriteAscii
                    )
                    #logger.info("Writing area.dat. Nr of points: " + str(size(area_block)))
                    dw_WriteSegmentOrExchangeData(
                        float(i), dwdir + "/includes_flow/area.dat", area_block, 1, WriteAscii
                    )
                    #logger.info("Writing velocity.dat. Nr of points: " + str(size(velocity)))
                    dw_WriteSegmentOrExchangeData(
                        float(i), dwdir + "/includes_flow/velocity.dat", velocity, 1, WriteAscii
                    )
                    
                    if Fews:
                        # for hyd-file set
                        dw_WriteSegmentOrExchangeData(i, comroot + ".are", area_block, 1, WriteAscii)
                        dw_WriteSegmentOrExchangeData(i, comroot + ".flo", flow_block, 1, WriteAscii)
                        dw_WriteSegmentOrExchangeData(i, comroot + ".vol", volume_block, 1, WriteAscii)
                    
                                
                ts = ts + 1
            #End time loop
        #End dynamic data for netcdf or map input files
        
        # Dynamic data from csv files
        if inputType == "csv":
            # Volumes of the compartments
            if not Emission:
                for c in np.arange(0, len(comp)):
                    # Read volume csv table
                    wflowVars = comp.mapstack[c].split("+")
                    volume_path = caseId + "/" + runId + "/" + wflowVars[0] + ".csv"
                    volume_csv = pd.read_csv(volume_path, sep = ",", header=0)
                    #Suppress timestep and 0 column
                    if "0" in volume_csv: del volume_csv["0"]
                    if "# Timestep" in volume_csv: del volume_csv["# Timestep"]
                    
                    #Add potential sum of variables
                    if len(wflowVars) > 1:
                        for v in range(len(wflowVars)-1):
                           volume_path = caseId + "/" + runId + "/" + wflowVars[v+1] + ".csv"
                           volume_csv2 = pd.read_csv(volume_path, sep = ",", header=0)
                           if "0" in volume_csv2: del volume_csv2["0"]
                           if "# Timestep" in volume_csv2: del volume_csv2["# Timestep"]
                           volume_csv = volume_csv + volume_csv2
                           
                    #Convert unit
                    if comp.Unit[c] == "mm" or comp.Unit[c] == "m":
                        # Adjusted area for the surface water
                        if comp.At2[c] == 0:
                            volume_csv = volume_csv * surfaceWater_block
                        else:
                            volume_csv = volume_csv * realarea_block
                        if comp.Unit[c] == "mm":
                            volume_csv = volume_csv / 1000
                    #Add first row from states
                    #In the states the unsaturated store is divided in several layers
                    if comp.wflowVar[c] == "UStoreDepth" :
                        #Get the number of layers
                        UStoreLayerThickness = configget(config, "model", "UStoreLayerThickness", "0")
                        if UStoreLayerThickness != "0":
                            maxLayers = len(UStoreLayerThickness.split(",")) + 1
                        else:
                            maxLayers = 1
                        volume_map = zeroMap*0.0
                        #Read and sum state maps of each layer of the unsaturated store
                        for USi in np.arange(0, maxLayers):
                            volume_map_path = caseId + "/instate/UStoreLayerDepth_" + str(USi) + ".map"
                            if os.path.exists(volume_map_path):
                                volume_map = volume_map + pcr.readmap(volume_map_path)
                        # Convert from mm to m3
                        volume_map = volume_map * reallength * reallength / 1000
                    
                    #If KinWaveVolume is used instead of WaterLevel:
                    elif (comp.wflowVar[c] == "KinWaveVolumeR+KinWaveVolumeL" or comp.wflowVar[c] == "KinWaveVolumeL+KinWaveVolumeR"):
                        volume_map_path1 = caseId + "/instate/" + "WaterLevelR" + ".map"
                        volume_map_path2 = caseId + "/instate/" + "WaterLevelL" + ".map"
                        if os.path.exists(volume_map_path1):
                            volume_map = pcr.readmap(volume_map_path1)
                        else:
                            volume_map = zeroMap*0.0
                        if os.path.exists(volume_map_path2):
                            volume_map = volume_map + pcr.readmap(volume_map_path2)
                        #Convert from m to m3
                        volume_map = volume_map * internalflowwidth * internalflowlength
                    #Other compartments    
                    else:
                        volume_map = zeroMap*0.0
                        wflowVars = comp.wflowVar[c].split("+")
                        for v in len(wflowVars):
                            volume_map_path = caseId + "/instate/" + wflowVars[v] + ".map"
                            if os.path.exists(volume_map_path):
                                volume_map = volume_map + pcr.readmap(volume_map_path)
                        #Convert from mm to m3
                        volume_map = volume_map * reallength * reallength / 1000
                    #Aggregate the state map
                    volume_map = pcr.areatotal(volume_map, modelmap)
                    volume_state = dw_pcrToDataBlock(volume_map, amap, Aggregation)
                    #Add as first row of csv table
                    volume_csv.loc[-1] = volume_state
                    volume_csv.index = volume_csv.index+1
                    volume_csv.sort_index(inplace=True) 
                    #Add minimum volume
                    volume_csv = volume_csv + 0.0001
                    #Concatenate the different volume tables
                    if c == 0:
                        volume_block = volume_csv
                    else:
                        volume_block = pd.concat([volume_block,volume_csv], axis = 1)
                #Write the volume file
                dw_WriteCSVdata(timestepsecs, dwdir + "/includes_flow/volume.dat", volume_block, 1, WriteAscii)
                  
            #Flow, exchange area and velocity
            if Emission:
                #In the emssion model some fluxes are in the pointer but not hydrology file
                hydflux = ~flux.EmPointerHyd.isin(["P"])
                flow_var = flux.mapstack[hydflux].values.astype(str)
                
            for f in np.arange(0, len(flow_var)):
                #Read flow csv table
                wflowVars = flow_var[f].split("+")
                flow_path = caseId + "/" + runId + "/" + wflowVars[0] + ".csv"
                flow_csv = pd.read_csv(flow_path, sep = ",", header=0)
                #Suppress timestep and 0 column
                if "0" in flow_csv: del flow_csv["0"]
                if "# Timestep" in flow_csv: del flow_csv["# Timestep"]
                
                #Add potential sum of variables
                if len(wflowVars) > 1:
                    for v in range(len(wflowVars)-1):
                       flow_path = caseId + "/" + runId + "/" + wflowVars[v+1] + ".csv"
                       flow_csv2 = pd.read_csv(flow_path, sep = ",", header=0)
                       if "0" in flow_csv2: del flow_csv2["0"]
                       if "# Timestep" in flow_csv2: del flow_csv2["# Timestep"]
                       flow_csv = flow_csv + flow_csv2
                
                #Get exchange area
                wflowVar_f = flux.wflowVar[flux.mapstack == flow_var[f]].values
                if (wflowVar_f[0] == "RiverRunoff+LandRunoff" or wflowVar_f[0] == "LandRunoff+RiverRunoff") :
                    level_mapstack = comp.mapstack[(comp.wflowVar == "KinWaveVolumeR+KinWaveVolumeL" or 
                                                    comp.wflowVar == "KinWaveVolumeL+KinWaveVolumeR")].values
                    levVars = level_mapstack[0].split("+")
                    level_path = caseId + "/" + runId + "/" + levVars[0] + "g.csv"
                    level_csv = pd.read_csv(level_path, sep = ",", header=0)
                    if "0" in level_csv: del level_csv["0"]
                    if "# Timestep" in level_csv: del level_csv["# Timestep"]
                    if len(levVars) > 1:
                        for l in range(len(levVars)-1):
                           level_path = caseId + "/" + runId + "/" + levVars[l+1] + ".csv"
                           level_csv2 = pd.read_csv(level_path, sep = ",", header=0)
                           if "0" in level_csv2: del level_csv2["0"]
                           if "# Timestep" in level_csv2: del level_csv2["# Timestep"]
                           level_csv = level_csv + level_csv2
                       
                    DCL_sample = dw_pcrToDataBlock(internalflowlength, gmap, Aggregation)
                    area_csv = level_csv / DCL_sample
                elif wflowVar_f[0] == "SatWaterFlux":
                    level_mapstack = comp.mapstack[comp.wflowVar == "SatWaterDepth"].values
                    level_path = caseId + "/" + runId + "/" + level_mapstack[0] + "g.csv"
                    level_csv = pd.read_csv(level_path, sep = ",", header=0)
                    if "0" in level_csv: del level_csv["0"]
                    if "# Timestep" in level_csv: del level_csv["# Timestep"]
                    area_csv = level_csv / 1000 * realarea_block^(1/2)
                else:
                    area_csv = flow_csv * 0.0 + realarea_block
                        
                #Convert flow unit
                unit_f = flux.Unit[flux.mapstack == flow_var[f]].values
                if Emission:
                    #Some fluxes need to be converted as a fraction of the from comp volume
                    convfrac_f = flux.EmFraction[flux.mapstack == flow_var[f]].values[0]
                else: convfrac_f = 0
                if unit_f[0] == "mm" and convfrac_f == 0:
                    flow_csv = flow_csv * area_csv / (1000 * timestepsecs)
                elif convfrac_f == 1:
                    fromcomp = flux.From[flux.mapstack == flow_var[f]].values
                    volume_mapstack = comp.mapstack[comp.ID == fromcomp[0]].values
                    volVars = volume_mapstack[0].split("+")
                    volume_csv = flow_csv * 0.0
                    for vl in len(volVars):
                        volume_path = caseId + "/" + runId + "/" + volVars[vl] + ".csv"
                        volume_csv2 = pd.read_csv(volume_path, sep = ",", header=0)
                        if "0" in volume_csv2: del volume_csv2["0"]
                        if "# Timestep" in volume_csv2: del volume_csv2["# Timestep"]
                        volume_csv = volume_csv + volume_csv2
                    #For the corresponding emission fluxes, normally both flow and volume are in mm
                    #Handle possible zero volumes
                    flow_csv[volume_csv == 0] = 0.0
                    volume_csv[volume_csv == 0] = 0.1
                    flow_csv = flow_csv / volume_csv                    
                
                #Get velocity
                if not Emission:
                    velocity_csv = flow_csv / area_csv
                #Add final line (copy of the last one, there is one more final volume)
                last_t = flow_csv.shape[0]
                flow_csv.loc[last_t] = flow_csv.loc[(last_t-1)]
                area_csv.loc[last_t] = area_csv.loc[(last_t-1)]
                if not Emission:
                    velocity_csv.loc[last_t] = velocity_csv.loc[(last_t-1)]
                #For emission, add to total_flow
                if Emission :
                    addtotflow_f = flux.EmPointerHyd[flux.mapstack == flow_var[f]].values
                    if (bool(np.isin(addtotflow_f[0], np.array(["T","PT","HT","PHT"])))):                            
                        try:
                            totalflow_csv = totalflow_csv + flow_csv
                        except:
                            totalflow_csv = flow_csv    
                
                #Concatenate flow, area, velocity
                #For emission, some fluxes area needed for total_flow but are not saved in the hydrology file
                if Emission:
                    addflowblock_f = flux.EmPointerHyd[flux.mapstack == flow_var[f]].values
                    if (bool(np.isin(addflowblock_f[0], np.array(["H","PH","HT","PHT"])))):
                        addflow = True
                    else: addflow = False
                else: addflow = True
                if addflow:
                    if f == 0:
                        flow_block = flow_csv
                        area_block = area_csv
                        if not Emission: velocity_block = velocity_csv
                    else:
                        flow_block = pd.concat([flow_block, flow_csv], axis = 1)
                        area_block = pd.concat([area_block, area_csv], axis = 1)
                        if not Emission:
                            velocity_block = pd.concat([velocity_block, velocity_csv], axis = 1)
                
            #Save flow, area, velocity
            if Emission:
                #Concatenate total_flow to flow_block for emission
                flow_block = pd.concat([flow_block, totalflow_csv], axis = 1)
                #Write the hydrology file
                dw_WriteCSVdata(timestepsecs, dwdir + "/includes_flow/hydrology.bin", flow_block, len(comp), WriteAscii)
                
            else:
                #Write the flow, area and velocity files
                dw_WriteCSVdata(timestepsecs, dwdir + "/includes_flow/flow.dat", flow_block, 1, WriteAscii)
                dw_WriteCSVdata(timestepsecs, dwdir + "/includes_flow/area.dat", area_block, 1, WriteAscii)
                dw_WriteCSVdata(timestepsecs, dwdir + "/includes_flow/velocity.dat", velocity_block, 1, WriteAscii)
        
        if Fews:
            # Generate attribute file
            atr_file = comroot + ".atr"
            logger.info("Writing attribute file to '%s'" % atr_file)
            dw_WriteAttributesFile(atr_file, NOSQ)
    
            # Generate hyd-file 
            hyd_file = comroot + "_unstructured.hyd"
            logger.info("Writing hyd-file to '%s'" % hyd_file)
            hydinfo = {}
            hydinfo["runid"] = runId
            hydinfo["tref"] = T0
            hydinfo["tstart"] = T0
            hydinfo["tstop"] = T0 + timedelta(seconds=(timeSteps - 1) * timestepsecs)
            hydinfo["tstep"] = timedelta(seconds=timestepsecs)
            hydinfo["noseg"] = NOSQ
            hydinfo["nosegh"] = NOSQ
            hydinfo["noqh"] = pointer.shape[0]
            hydinfo["noqv"] = 0
            dw_WriteHydFile(hyd_file, hydinfo)
        
    #End write dynamic
    
    if template_ini_path != 'None':                                                                                                
        logger.info('Writing model ini file with template from ' + template_ini_path)
        if Emission: model_ini = 'EM'
        else: model_ini = 'WAQ'
        model_ini_path  = os.path.join(dwdir, model_ini + '.inp')
        if os.path.exists(template_ini_path):
            if os.path.exists(model_ini_path):
                os.remove(model_ini_path)
            shutil.copy(template_ini_path, model_ini_path)
        else:
            logger.error('Could not find specified template ini file!')        


if __name__ == "__main__":
    main()
