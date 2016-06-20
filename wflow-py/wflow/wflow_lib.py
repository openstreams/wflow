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
wflow_lib - terrain analysis and hydrological library
-----------------------------------------------------

The goal of this module is to make a series functions to upscale maps (DEM)
and to  maintain as much of the information in a detailled dem when upscaling
to a coarser DEM. These include:

    - river length (per cell)
    - river network location
    - elevation distribution
    - other terrain analysis

the wflow_prepare scripts use this library extensively.

$Author: schelle $
$Id: wflow_lib.py 808 2013-10-04 19:42:43Z schelle $
$Rev: 808 $
"""



import getopt
import os
import os.path
import sys


import osgeo.gdal as gdal
from osgeo.gdalconst import *
from pcraster import *
from pcraster.framework import *
import scipy
import numpy as np
import netCDF4 as nc4
import gzip, zipfile


def idtoid(sourceidmap, targetidmap,valuemap):
    """
    tranfer the values from valuemap at the point id's in sourceidmap to the areas in targetidmap.

    :param pointmap:
    :param areamap:
    :param valuemap:
    :return:
    """

    _area = pcr2numpy(targetidmap,0.0).copy()
    _pt = pcr2numpy(sourceidmap,0.0).copy()
    _val = pcr2numpy(valuemap,0.0).copy()

    for val in np.unique(_pt):
        if val > 0:  #
            _area[_area == val] = np.mean(_val[_pt == val])

    retmap = numpy2pcr(Scalar,_area,0.0)

    return retmap


def simplereservoir(storage, inflow, maxstorage, target_perc_full, maximum_Q, demand, minimum_full_perc, ReserVoirLocs, timestepsecs=86400):
    """

    :param storage: initial storage m^3
    :param inflow: inflow m^3/s
    :param maxstorage: maximum storage (above which water is spilled) m^3
    :param target_perc_full: target fraction full (of max storage) -
    :param maximum_Q: maximum Q to release m^3/s if below spillway
    :param demand: water demand (all combined) m^3/s
    :param minimum_full_perc: target minimum full fraction (of max storage) -
    :param ReserVoirLocs: map with reservoir locations
    :param timestepsecs: timestep of the model in seconds (default = 86400)
    :return: storage (m^3), outflow (m^3/s), PercentageFull (0-1), Release (m^3/sec)
    """

    inflow = ifthen(boolean(ReserVoirLocs), inflow)
    oldstorage = storage
    storage = storage + (inflow * timestepsecs)
    percfull = ((storage + oldstorage) * 0.5) / maxstorage
    # first determine minimum (environmental) flow using a simple sigmoid curve to scale for target level
    fac = sCurve(percfull, a=minimum_full_perc, c=30.0)
    demandRelease = fac * demand * timestepsecs
    storage = storage - demandRelease

    # Re-determine percfull
    percfull = ((storage + oldstorage) * 0.5) / maxstorage

    wantrel = max(0.0, storage - (maxstorage * target_perc_full))
    # Assume extra maximum Q if spilling
    overflowQ = (percfull - 1.0) * (storage - maxstorage)
    torelease = min(wantrel, overflowQ + maximum_Q * timestepsecs)
    storage = storage - torelease
    outflow = (torelease + demandRelease) / timestepsecs
    percfull = storage / maxstorage

    return storage, outflow, percfull, demandRelease/timestepsecs


Verbose=0

def lddcreate_save(lddname, dem, force, corevolume=1E35, catchmentprecipitation=1E35, corearea=1E35, outflowdepth=1E35):
    """
    Creates an ldd if a file does not exists or if the force flag is used

    input: 
        - lddname (name of the ldd to create)
        - dem (actual dem)
        - force (boolean to force recreation of the ldd)
        - outflowdepth (set to 10.0E35 normally but smaller if needed)
        
    Output:
        - the LDD
        
    """
    if os.path.exists(lddname) and not force:
        if Verbose:
            print("Returning existing ldd", lddname)
            return readmap(lddname)
    else:
        if Verbose:
            print("Creating ldd", lddname)
            LDD = lddcreate(dem, 10.0E35, outflowdepth, 10.0E35, 10.0E35)
            report(LDD, lddname)
            return LDD

        
def configget(config,section,var,default):
    """
    
    Gets a string from a config file (.ini) and returns a default value if
    the key is not found. If the key is not found it also sets the value 
    with the default in the config-file
    
    Input:
        - config - python ConfigParser object
        - section - section in the file
        - var - variable (key) to get
        - default - default string
        
    Returns:
        - string - either the value from the config file or the default value
        
        
    """
    Def = False
    try:
        ret = config.get(section,var)
    except:
        Def = True
        ret = default
        configset(config,section,var,default, overwrite=False)
    
    default = Def
    return ret       


def configset(config,section,var,value, overwrite=False):
    """   
    Sets a string in the in memory representation of the config object
    Deos NOT overwrite existing values if overwrite is set to False (default)
    
    Input:
        - config - python ConfigParser object
        - section - section in the file
        - var - variable (key) to set
        - value - the value to set
        - overwrite (optional, default is False)
   
    Returns:
        - nothing
        
    """
    
    if not config.has_section(section):
        config.add_section(section)
        config.set(section,var,value)
    else:     
        if not config.has_option(section,var):
            config.set(section,var,value)
        else:
            if overwrite:
                config.set(section,var,value)
                

def configsection(config,section):
    """
    gets the list of keys in a section
    
    Input:
        - config
        - section
        
    Output:
        - list of keys in the section
    """
    try:
        ret = config.options(section)
    except:
        ret = []
        
    return ret


def getrows():
    """
    returns the number of rows in the current map
    
    Input:
        - -
        
    Output:
        - nr of rows in the current clonemap as a scalar
    """
    a = pcr2numpy(celllength(),numpy.nan).shape[0]
    
    return a

def getcols():
    """
    returns the number of columns in the current map
    
    Input:
        - -
        
    Output:
        - nr of columns in the current clonemap as a scalar
    """
    a = pcr2numpy(celllength(),numpy.nan).shape[1]
    
    return a   

def getgridparams():
    """ return grid parameters in a python friendly way
    
    Output:    
        [ Xul, Yul, xsize, ysize, rows, cols]
    
        - xul - x upper left centre
        - yul - y upper left centre
        - xsize - size of a cell in x direction
        - ysize - size of a cell in y direction
        - cols - number of columns
        - rows - number of rows
        - xlr
        - ylr
    """
    
    # x and Y are the same for now
    xy = pcr2numpy(celllength(),numpy.nan)[0,0]
    xu = pcr2numpy(xcoordinate(1),numpy.nan)[0,0]
    yu = pcr2numpy(ycoordinate(1),numpy.nan)[0,0]
    ylr = pcr2numpy(ycoordinate(1),numpy.nan)[getrows()-1,getcols()-1]
    xlr = pcr2numpy(xcoordinate(1),numpy.nan)[getrows()-1,getcols()-1]
    
    return [xu, yu, xy, xy, getrows(), getcols(),xlr,ylr]
        

def snaptomap(points,mmap):
    """
    Snap the points in _points_ to nearest non missing
    values in _mmap_. Can be used to move gauge locations
    to the nearest rivers.
    
    Input: 
        - points - map with points to move
        - mmap - map with points to move to

    Return: 
        - map with shifted points
    """
    points = cover(points,0)
    # Create unique id map of mmap cells
    unq = nominal(cover(uniqueid(defined(mmap)),scalar(0.0)))
    # Now fill holes in mmap map with lues indicating the closes mmap cell.
    dist_cellid = scalar(spreadzone(unq,0,1))
    # Get map with values at location in points with closes mmap cell
    dist_cellid = ifthenelse(points > 0, dist_cellid, 0)
    # Spread this out
    dist_fill = spreadzone(nominal(dist_cellid),0,1)
    # Find the new (moved) locations
    npt = uniqueid(boolean(ifthen(dist_fill == unq, unq)))
    # Now recreate the original value in the points maps
    ptcover = spreadzone(cover(points,0),0,1)
    # Now get the org point value in the pt map
    nptorg = ifthen(npt > 0, ptcover)
    
    
    return nptorg

def riverlength(ldd,order):
    """
    Determines the length of a river using the ldd. 
    only determined for order and higher.

    Input: 
        - ldd, order (streamorder)
        
    Returns: 
        - totallength,lengthpercell, streamorder
    """    
    strorder=streamorder(ldd)
    strorder=ifthen(strorder >= ordinal(order),strorder)
    dist=max(celllength(),ifthen(boolean(strorder),downstreamdist(ldd)))
    
    return catchmenttotal(cover(dist,0),ldd), dist, strorder
    

def upscale_riverlength(ldd,order, factor):
    """
    Upscales the riverlength using 'factor'
    The resulting maps can be resampled (e.g. using resample.exe) by factor and should
    include the accurate length as determined with the original higher 
    resolution maps.  This function is **depricated**,  
    use are_riverlength instead as this version
    is very slow for large maps
    
    Input: 
        - ldd 
        - minimum streamorder to include
        
    Output: 
        - distance per factor cells 
    """
    
    strorder=streamorder(ldd)
    strorder=ifthen(strorder >= order,strorder)
    dist=cover(max(celllength(),ifthen(boolean(strorder),downstreamdist(ldd))),0)   
    totdist=max(ifthen(boolean(strorder),windowtotal(ifthen(boolean(strorder),dist),celllength() * factor)),dist)
    
    return totdist

def area_riverlength_factor(ldd, Area, Clength):
    """
    ceates correction factors for riverlength for 
    the largest streamorder in each area
    
    Input: 
        - ldd
        - Area
        - Clength (1d length of a cell (sqrt(Area))
        
    Output:
        - distance per area
    
    """
    strorder=streamorder(ldd)
    strordermax=areamaximum(strorder,Area)
    dist = downstreamdist(ldd)
    # count nr of strorder cells in each area
    nr = areatotal(ifthen(strorder == strordermax,dist),Area)
    #N = sqrt(areatotal(scalar(boolean(Area)),Area))
    N = Clength
    factor = nr/N
    
    
    return factor

def area_river_burnin(ldd, dem, order,Area):
  """
  Calculates the lowest values in as DEM for each erea in an area map for 
  river of order *order*
  
  Input: 
      - ldd
      - dem
      - order
      - Area map
      
  Output:
      - dem 
  """
  strorder = streamorder(ldd)
  strordermax=areamaximum(strorder,Area)
  maxordcell = ifthen(strordermax > order, strordermax)
  riverdem = areaminimum(dem,Area)
  
  return ifthen(boolean(maxordcell),riverdem)
  

def area_percentile(inmap,area,n,order,percentile):
  """
  calculates percentile of inmap per area
  n is the number of points in each area, 
  order, the sorter order of inmap per area (output of
  areaorder(inmap,area))
  n is the output of areatotal(spatial(scalar(1.0)),area)
  
  Input:
      - inmap
      - area map
      - n
      - order (riverorder)
      - percentile 
      
  Output:
      - percentile map
      
  """
  i = rounddown((n * percentile)/100.0 + 0.5) # index in order map
  perc = ifthen(i == order, inmap)
  
  return areaaverage(perc,area)

    
def find_outlet(ldd):
    """
    Tries to find the outlet of the largest catchment in the Ldd
    
    Input: 
        - Ldd
        
    Output: 
        - outlet map (single point in the map)
    """
    largest = mapmaximum(catchmenttotal(spatial(scalar(1.0)),ldd))
    outlet = ifthen(catchmenttotal(1.0,ldd) == largest,spatial(scalar(1.0)))
    
    return outlet
    

def subcatch(ldd,outlet):
    """
    Determines a subcatchment map using LDD and outlet(s). In the resulting
    subcatchment map the i's of the catchment are determiend by the id's of
    the outlets.
    
    Input:
        - ldd
        - Outlet - maps with points for each outlet. 
        
    Output:
        - map of subcatchments
    """
    subcatch=subcatchment(ldd, ordinal(outlet))
    
    return subcatch
    
def areastat(Var,Area):
    """
    Calculate several statistics of *Var* for each unique id in *Area*
    
    Input: 
        - Var
        - Area
        
    Output: 
        - Standard_Deviation,Average,Max,Min
    
    """
    Avg = areaaverage(Var,Area)
    Sq = (Var - Avg)**2
    N = areatotal(spatial(cellarea()),Area)/cellarea()
    Sd = (areatotal(Sq,Area)/N)**0.5
    Max = areamaximum(Var,Area)
    Min = areaminimum(Var,Area)
    
    return Sd,Avg,Max,Min
    


def checkerboard(mapin,fcc):
    """
    checkerboard create a checkerboard map with unique id's in a 
    fcc*fcc cells area. The resulting map can be used
    to derive statistics for (later) upscaling of maps (using the fcc factor)
    
    .. warning: use with unitcell to get most reliable results!
    
    Input: 
        - map (used to determine coordinates)
        - fcc (size of the areas in cells)
        
    Output: 
        - checkerboard type map
    """
    msker = defined(mapin)
    ymin = mapminimum(ycoordinate(msker))
    yc = (ycoordinate((msker))-ymin)/celllength()
    yc = rounddown(yc/fcc)
    #yc = yc/fcc
    xmin = mapminimum(xcoordinate((msker)))
    xc = (xcoordinate((msker)) - xmin)/celllength()
    xc = rounddown(xc/fcc)
    #xc = xc/fcc
    
    yc = yc * (mapmaximum(xc) + 1.0)
    
    xy = ordinal(xc + yc)
    
    return xy
   

def subcatch_order_a(ldd,oorder):
    """
    Determines subcatchments using the catchment order
    
    This version uses the last cell BELOW order to derive the
    catchments. In general you want the \_b version
    
    Input:
        - ldd
        - order - order to use
    
    Output:
        - map with catchment for the given streamorder
    """
    outl = find_outlet(ldd)
    large = subcatchment(ldd,boolean(outl))   
    stt = streamorder(ldd)
    sttd = downstream(ldd,stt)
    pts = ifthen((scalar(sttd) - scalar(stt)) > 0.0,sttd)
    dif = upstream(ldd,cover(ifthen(large,uniqueid(boolean(ifthen(stt == ordinal(oorder), pts)))),0))
    dif = cover(scalar(outl),dif) # Add catchment outlet
    dif = ordinal(uniqueid(boolean(dif)))
    sc = subcatchment(ldd,dif)
    
    return sc, dif, stt

def subcatch_order_b(ldd,oorder,Largest=False,sizelimit=0):
    """
    Determines subcatchments using the catchment order

    This version uses the bottommost cell of order
    If Largest is true the analysis is only done for the largest basin
    found in the ldd
    
    Input:
        - ldd
        - oorder - order to use
        - largest - toggle, default = False
        - sizelimit - smallest catchments to include, default is all (sizelimit=0) in number of cells
    
    Output:
        - map with catchment for the given streamorsder
    """
    outl = find_outlet(ldd)
    large = subcatchment(ldd,boolean(outl))   
    stt = streamorder(ldd)
    sttd = downstream(ldd,stt)
    pts = ifthen((scalar(sttd) - scalar(stt)) > 0.0,sttd)
    if Largest:
        dif = ifthen(large,uniqueid(boolean(ifthen(stt == ordinal(oorder), pts))))
    else:
        dif = uniqueid(cover(boolean(pit(ldd)),boolean(ifthen(stt == ordinal(oorder), pts))))
    
    scsize = catchmenttotal(1,ldd)
    dif = ordinal(uniqueid(boolean(ifthen(scsize >= sizelimit,dif))))
    sc = subcatchment(ldd,dif)
    
    return sc, dif, stt


def getRowColPoint(in_map,xcor,ycor):
    """
    returns the row and col in a map at the point given.
    Works but is rather slow.
    
    Input:
        - in_map - map to determine coordinates from
        - xcor - x coordinate
        - ycor - y coordinate
        
    Output:
        - row, column
    """
    x = pcr2numpy(xcoordinate(boolean(scalar(in_map) + 1.0)),numpy.nan)
    y = pcr2numpy(ycoordinate(boolean(scalar(in_map) + 1.0)),numpy.nan)
    XX = pcr2numpy(celllength(),0.0)
    tolerance = 0.5 # takes a single point

    diffx = x - xcor
    diffy = y - ycor
    col_ =  numpy.absolute(diffx) <= (XX[0,0] * tolerance)  # cellsize
    row_ =  numpy.absolute(diffy) <= (XX[0,0] * tolerance)# cellsize
    point = (col_ * row_)
    
    
    return point.argmax(0).max(), point.argmax(1).max()

def getValAtPoint(in_map,xcor,ycor):
    """
    returns the value in a map at the point given.
    works but is rather slow.
    
    Input:
        - in_map - map to determine coordinates from
        - xcor - x coordinate
        - ycor - y coordinate
        
    Output:
        - value
    """
    x = pcr2numpy(xcoordinate(defined(in_map)),numpy.nan)
    y = pcr2numpy(ycoordinate(defined(in_map)),numpy.nan)
    XX = pcr2numpy(celllength(),0.0)
    themap =pcr2numpy(in_map,numpy.nan)
    tolerance = 0.5 # takes a single point

    diffx = x - xcor
    diffy = y - ycor
    col_ =  numpy.absolute(diffx) <= (XX[0,0] * tolerance)  # cellsize
    row_ =  numpy.absolute(diffy) <= (XX[0,0] * tolerance)# cellsize
    point = (col_ * row_)
    pt = point.argmax()
    
    return themap.ravel()[pt]
   

def points_to_map(in_map,xcor,ycor,tolerance):
    """
    Returns a map with non zero values at the points defined
    in X, Y pairs. It's goal is to replace the pcraster col2map program. 
    
    tolerance should be 0.5 to select single points
    Performance is not very good and scales linear with the number of points
    
    
    Input:
        - in_map - map to determine coordinates from
        - xcor - x coordinate (array or single value)
        - ycor - y coordinate (array or single value)
        - tolerance - tolerance in cell units. 0.5 selects a single cell\
        10 would select a 10x10 block of cells
        
    Output:
        - Map with values burned in. 1 for first point, 2 for second and so on
    """
    point = in_map * 0.0
    
    x = pcr2numpy(xcoordinate(defined(in_map)),numpy.nan)
    y = pcr2numpy(ycoordinate(defined(in_map)),numpy.nan)
    XX = pcr2numpy(celllength(),0.0)
    
    # simple check to use both floats and numpy arrays
    try:
        c = xcor.ndim
    except:
        xcor = numpy.array([xcor])
        ycor = numpy.array([ycor])

    # Loop over points and "burn in" map
    for n in range(0,xcor.size):
        if Verbose:
            print(n)
        diffx = x - xcor[n]
        diffy = y - ycor[n]        
        col_ =  numpy.absolute(diffx) <= (XX[0,0] * tolerance)  # cellsize
        row_ =  numpy.absolute(diffy) <= (XX[0,0] * tolerance)# cellsize
        point =  point + numpy2pcr(Scalar,((col_ * row_) * (n+1)),numpy.nan)
    
    return ordinal(point)


def detdrainlength(ldd,xl,yl):
    """
    Determines the drainaige length (DCL) for a non square grid
    
    Input:
        - ldd - drainage network
        - xl - length of cells in x direction
        - yl - length of cells in y direction
        
    Output:
        - DCL
    """
    # take into account non-square cells    
    # if ldd is 8 or 2 use Ylength
    # if ldd is 4 or 6 use Xlength
    draindir = scalar(ldd)
    slantlength = sqrt(xl**2 + yl**2)
    drainlength = ifthenelse(draindir == 2,yl,
                             ifthenelse(draindir == 8,yl,
                                        ifthenelse(draindir == 4, xl,
                                                   ifthenelse(draindir == 6,xl,slantlength))))

                                                           
    return drainlength   

def detdrainwidth(ldd,xl,yl):
    """
    Determines width of drainage over DEM for a non square grid
    
    Input:
        - ldd - drainage network
        - xl - length of cells in x direction
        - yl - length of cells in y direction
        
    Output:
        - DCL
    """
    # take into account non-square cells    
    # if ldd is 8 or 2 use Xlength
    # if ldd is 4 or 6 use Ylength
    draindir = scalar(ldd)
    slantwidth = (xl + yl) * 0.5
    drainwidth = ifthenelse(draindir == 2,xl,
                             ifthenelse(draindir == 8,xl,
                                        ifthenelse(draindir == 4, yl,
                                                   ifthenelse(draindir == 6,yl,slantwidth))))
    return drainwidth


def classify(inmap,lower=[0,10,20,30],upper=[10,20,30,40],classes=[2,2,3,4]):
    """
    classify a scaler maps accroding to the boundaries given in classes.

    """

    result=ordinal(cover(-1))
    for l, u, c in zip(lower, upper,classes):
        result = cover(ifthen(inmap >= l,ifthen(inmap < u,ordinal(c))),result)

    return ifthen(result >=0,result)


def derive_HAND(dem, ldd, accuThreshold, rivers=None, basin=None):
    """
    Function derives Height-Above-Nearest-Drain.
    See http://www.sciencedirect.com/science/article/pii/S003442570800120X
    Input:
        dem -- pcraster object float32, elevation data
        ldd -- pcraster object direction, local drain directions
        accuThreshold -- upstream amount of cells as threshold for river
            delineation
        rivers=None -- you can provide a rivers layer here. Pixels that are 
                        identified as river should have a value > 0, other
                        pixels a value of zero.
        basin=None -- set a boolean pcraster map where areas with True are estimated using the nearest drain in ldd distance
                        and areas with False by means of the nearest friction distance. Friction distance estimated using the 
                        upstream area as weight (i.e. drains with a bigger upstream area have a lower friction)
                        the spreadzone operator is used in this case.
    Output:
        hand -- pcraster bject float32, height, normalised to nearest stream
        dist -- distance to nearest stream measured in cell lengths
            according to D8 directions
    """
    if rivers is None:
        stream = ifthenelse(accuflux(ldd, 1) >= accuThreshold,
                                boolean(1), boolean(0))
    else:
        stream = boolean(cover(rivers, 0))
    
    height_river = ifthenelse(stream, ordinal(dem*100), 0)
    if basin is None:
        up_elevation = scalar(subcatchment(ldd, height_river))
    else:
        drainage_surf = ifthen(rivers, accuflux(ldd, 1))
        weight = 1./scalar(spreadzone(cover(ordinal(drainage_surf), 0), 0, 0))
        up_elevation = ifthenelse(basin, scalar(subcatchment(ldd, height_river)), scalar(spreadzone(height_river, 0, weight)))
        # replace areas outside of basin by a spread zone calculation.
    hand = max(scalar(ordinal(dem*100))-up_elevation, 0)/100
    dist = ldddist(ldd, stream, 1)
        
    return hand, dist

    
def sCurve(X,a=0.0,b=1.0,c=1.0):
    """
    sCurve function:
        
    Input:
        - X input map
        - C determines the steepness or "stepwiseness" of the curve. 
          The higher C the sharper the function. A negative C reverses the function.
        - b determines the amplitude of the curve
        - a determines the centre level (default = 0)
    
    Output:
        - result
    """
    try:
        s = 1.0/(b + exp(-c * (X-a)))
    except:
        s = 1.0 / (b + np.exp(-c * (X - a)))
    return s

def sCurveSlope(X,a=0.0,b=1.0,c=1.0):
    """
    First derivative of the sCurve defined by a,b,c at point X
    
    Input:
        - X - value to calculate for
        - a 
        - b
        - c
    
    Output:
        - first derivative (slope) of the curve at point X
    """
    sc = sCurve(X,a=a,b=b,c=c)
    slope = sc * (1 - sc)
    return slope    


def Gzip(fileName, storePath=False, chunkSize=1024*1024):
        """
        Usage: Gzip(fileName, storePath=False, chunksize=1024*1024)
        Gzip the given file to the given storePath and then remove the file.
        A chunk size may be selected. Default is 1 megabyte
        Input:
            fileName:   file to be GZipped
            storePath:  destination folder. Default is False, meaning the file will be zipped to its own folder
            chunkSize:  size of chunks to write. If set too large, GZip will fail with memory problems
        """
        import gzip
        if not storePath:
            pathName = os.path.split(fileName)[0]
            fileName = os.path.split(fileName)[1]
            curdir   = os.path.curdir
            os.chdir(pathName)
        # open files for reading / writing
        r_file = open(fileName, 'rb')
        w_file = gzip.GzipFile(fileName + '.gz', 'wb', 9)
        dataChunk = r_file.read(chunkSize)
        while dataChunk:
            w_file.write(dataChunk)
            dataChunk = r_file.read(chunkSize)
        w_file.flush()
        w_file.close()
        r_file.close()
        os.unlink(fileName) #We don't need the file now
        if not storePath:
            os.chdir(curdir)



# These come from GLOFRIS_Utils

def Gzip(fileName, storePath=False, chunkSize=1024*1024):
        """
        Usage: Gzip(fileName, storePath=False, chunksize=1024*1024)
        Gzip the given file to the given storePath and then remove the file.
        A chunk size may be selected. Default is 1 megabyte
        Input:
            fileName:   file to be GZipped
            storePath:  destination folder. Default is False, meaning the file will be zipped to its own folder
            chunkSize:  size of chunks to write. If set too large, GZip will fail with memory problems
        """
        if not storePath:
            pathName = os.path.split(fileName)[0]
            fileName = os.path.split(fileName)[1]
            curdir   = os.path.curdir
            os.chdir(pathName)
        # open files for reading / writing
        r_file = open(fileName, 'rb')
        w_file = gzip.GzipFile(fileName + '.gz', 'wb', 9)
        dataChunk = r_file.read(chunkSize)
        while dataChunk:
            w_file.write(dataChunk)
            dataChunk = r_file.read(chunkSize)
        w_file.flush()
        w_file.close()
        r_file.close()
        os.unlink(fileName) #We don't need the file now
        if not storePath:
            os.chdir(curdir)

def zipFiles(fileList, fileTarget):
    """
    Usage: zipFiles(fileList, fileTarget)
    zip the given list of files to the given target file
    Input:
        fileList:   list of files to be zipped
        fileTarget: target zip-file
    """
    zout = zipfile.ZipFile(fileTarget, "w", compression=zipfile.ZIP_DEFLATED)
    for fname in fileList:
        zout.write(fname, arcname=os.path.split(fname)[1])
    zout.close()



def readMap(fileName, fileFormat):
    """ Read geographical file into memory
    """
    # Open file for binary-reading
    mapFormat = gdal.GetDriverByName(fileFormat)
    mapFormat.Register()
    ds = gdal.Open(fileName)
    if ds is None:
        print 'Could not open ' + fileName + '. Something went wrong!! Shutting down'
        sys.exit(1)
        # Retrieve geoTransform info
    geotrans = ds.GetGeoTransform()
    originX = geotrans[0]
    originY = geotrans[3]
    resX    = geotrans[1]
    resY    = geotrans[5]
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x = numpy.linspace(originX+resX/2,originX+resX/2+resX*(cols-1),cols)
    y = numpy.linspace(originY+resY/2,originY+resY/2+resY*(rows-1),rows)
    # Retrieve raster
    RasterBand = ds.GetRasterBand(1) # there's only 1 band, starting from 1
    data = RasterBand.ReadAsArray(0,0,cols,rows)
    FillVal = RasterBand.GetNoDataValue()
    RasterBand = None
    ds = None
    return x, y, data, FillVal

def writeMap(fileName, fileFormat, x, y, data, FillVal):
    """ Write geographical data into file"""

    verbose = False
    gdal.AllRegister()
    driver1 = gdal.GetDriverByName('GTiff')
    driver2 = gdal.GetDriverByName(fileFormat)

        # Processing
    if verbose:
        print 'Writing to temporary file ' + fileName + '.tif'
    # Create Output filename from (FEWS) product name and data and open for writing
    TempDataset = driver1.Create(fileName + '.tif',data.shape[1],data.shape[0],1,gdal.GDT_Float32)
    # Give georeferences
    xul = x[0]-(x[1]-x[0])/2
    yul = y[0]+(y[0]-y[1])/2
    TempDataset.SetGeoTransform( [ xul, x[1]-x[0], 0, yul, 0, y[1]-y[0] ] )
    # get rasterband entry
    TempBand = TempDataset.GetRasterBand(1)
    # fill rasterband with array
    TempBand.WriteArray(data,0,0)
    TempBand.FlushCache()
    TempBand.SetNoDataValue(FillVal)
    # Create data to write to correct format (supported by 'CreateCopy')
    if verbose:
        print 'Writing to ' + fileName + '.map'
    outDataset = driver2.CreateCopy(fileName, TempDataset, 0)
    TempDataset = None
    outDataset = None
    if verbose:
        print 'Removing temporary file ' + fileName + '.tif'
    os.remove(fileName + '.tif');

    if verbose:
        print 'Writing to ' + fileName + ' is done!'

