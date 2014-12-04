"""
Collection of utils for working with the
PCRaster python bindings.

"""

import os.path
import numpy
from math import *
import sys
import csv


from pcraster import *
from pcraster.framework import *
import logging
import logging.handlers


  
    

def lattometres(lat):
    """"
    Determines the length of one degree lat/long at a given latitude (in meter).
    Code taken from http:www.nga.mil/MSISiteContent/StaticFiles/Calculators/degree.html
    Input: map with lattitude values for each cell
    Returns: length of a cell lat, length of a cell long
    """
    #radlat = spatial(lat * ((2.0 * math.pi)/360.0))
    #radlat = lat * (2.0 * math.pi)/360.0
    radlat = spatial(lat) # pcraster cos/sin work in degrees!
    
    
    m1 = 111132.92        # latitude calculation term 1
    m2 = -559.82        # latitude calculation term 2
    m3 = 1.175            # latitude calculation term 3
    m4 = -0.0023        # latitude calculation term 4
    p1 = 111412.84        # longitude calculation term 1
    p2 = -93.5            # longitude calculation term 2
    p3 = 0.118            # longitude calculation term 3
    # # Calculate the length of a degree of latitude and longitude in meters
    
    latlen = m1 + (m2 * cos(2.0 * radlat)) + (m3 * cos(4.0 * radlat)) + (m4 * cos(6.0 * radlat))
    longlen = (p1 * cos(radlat)) + (p2 * cos(3.0 * radlat)) + (p3 * cos(5.0 * radlat))
        
    return latlen, longlen  
    
def detRealCellLength(ZeroMap,sizeinmetres):
    """
    Determine cellength. Always returns the length
    in meters.
    """
    
    if sizeinmetres:
            reallength = celllength()
            xl = celllength()
            yl = celllength()
    else:
        aa = ycoordinate(boolean(cover(ZeroMap + 1,1)))
        yl, xl = lattometres(aa)
           
        xl = xl * celllength()
        yl = yl * celllength()
        # Average length for surface area calculations. 
        
        reallength = (xl + yl) * 0.5
        
    return xl,yl,reallength


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)




def setlogger(logfilename,loggername, thelevel=logging.INFO):
    """
    Set-up the logging system and return a logger object. Exit if this fails
    """

    try:    
        #create logger
        logger = logging.getLogger(loggername)
        if not isinstance(thelevel, int):
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(thelevel)
        ch = logging.FileHandler(logfilename,mode='w')
        console = logging.StreamHandler()
        console.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
        #create formatter
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        #add formatter to ch
        ch.setFormatter(formatter)
        console.setFormatter(formatter)
        #add ch to logger
        logger.addHandler(ch)
        logger.addHandler(console)
        logger.debug("File logging to " + logfilename)
        return logger
    except IOError:
        print "ERROR: Failed to initialize logger with logfile: " + logfilename
        sys.exit(2)
    
def readmapSave(pathtomap, default):
    """
    Adpatation of readmap that returns a default map if the map cannot be found
    """
    if os.path.isfile(pathtomap):
        return readmap(pathtomap)
    else:
        return scalar(default)
        
def readtss(nname):
    """Reads a RCraster .tss file into a numpy array. 
    Error handling is minimal. The first column that
    contains the timestep is not returned.
    returns:
        matrix with data
        header
    """
    
    head = []
    if os.path.exists(nname):
        # determine the number of header lines to skip from
        # the second line in the tss file
        ifile = open(nname,"r")
        line = ifile.readline()
        line = ifile.readline()
        toskip = int(line) + 2
        line = ifile.readline()
        # read header
        for i in range(toskip - 3):
            line = ifile.readline()
            head.append(line.strip())
            
        ifile.close()
        
        
        mat = numpy.loadtxt(nname,skiprows=toskip)
        mis = mat == 1e+31
        mat[mis] = numpy.nan
        
        if len(mat.shape) > 1:
            return mat[:,1:], head
            #dumm = mat[:,1:].copy()
            #return numpy.vstack((dumm,mat[:,1:])), head
        else:
            return mat[1:], head
            #dumm = mat[1:].copy()
            #dumm[:] = 1E-31
            #return numpy.vstack((dumm,mat[1:])), head
    else:
        print nname + " does not exists."
    
    return

def interpolategauges(inputmap,method):
    """"
    Interpolate time series gauge data onto a grid using different methods
    inputmap: map with points data for a single timeste
    method: string indicating the method
        inv
        pol
        
    input: inputmap, method
    returns: interpolated map
    """   
     
    if method == "inv":
        result = inversedistance(1,inputmap,3,0,0)
    elif method == "pol":
        Unq = uniqueid(boolean(inputmap+1))
        result = spreadzone(ordinal(cover(Unq,0)),0,1)
        result = areaaverage(inputmap,result); 
    else:
        Unq = uniqueid(boolean(inputmap+1))
        result = spreadzone(ordinal(cover(Unq,0)),0,1)
        result = areaaverage(inputmap,result);
    
    return result 




_tableToMap_LastTbl = {}
_tableToMap_LastMap = {}


def tableToMapSparse (step, table, map):
    """Reads a pcraster.tbl file for step and assigns using the map in map.
    The behaviour of is a bit similar to the timeinputSparse
    command but in this case for both the tbl file and the map file.
   
    Input: step (int), table (string, path, without the .tbl extension), map
          (ordinal map, without the .map extension)

    How to create your maps/tables:
    if the maps for the timstep is not found the
    defaultmap (without the step appended) is returned or the last map
    that has been found.

    How to use this:
    - if you create the following maps and step ranges from 1 to 400 and
    the map name is "LAI":
    LAI.map
    LAI10.map
    LAI120.map
    LAI300.map

    - LAI.map will be use for step between 1 and 9 (default)
    - LAI10.map will be used between 10 and 119
    etc....
    The same holds for the tables.
    
    
    """
    global _tableToMap_LastTbl
    global _tableToMap_LastMap
    global debug
    
    # construct filenames
    fname_map = map + str(step) + ".map"
    fname_tbl = table + str(step) + ".tbl"

    
    if os.path.exists(fname_map):
	    print "found: " + fname_map
	    _tableToMap_LastMap[map] = step
  	    
    if os.path.exists(fname_tbl):
	    print "found: " + fname_tbl
	    _tableToMap_LastTbl[table] = step

    if _tableToMap_LastTbl.has_key(table) == False:
	    fname_tbl = table + ".tbl"
    else:
	    fname_tbl = table + str(_tableToMap_LastTbl[table]) + ".tbl"

    if _tableToMap_LastMap.has_key(map) == False:
	    fname_map = map + ".map"
    else:
	    fname_map = map + str(_tableToMap_LastMap[map]) + ".map"
        
    
    rmat = lookupscalar(str(fname_tbl),str(fname_map))

    return rmat



# def correctrad(Day,Hour,Lat,Lon,Slope,Aspect,Altitude):
#     """ Determines radiation over a DEM assuming clear sky"""
#
#     Sc  = 1367.0          # Solar constant (Gates, 1980) [W/m2]
#     Trans   = 0.6             # Transmissivity tau (Gates, 1980)
#     pi = 3.1416
#     AtmPcor = pow(((288.0-0.0065*Altitude)/288.0),5.256)
#     Lat = Lat * pi/180
#     ##########################################################################
#     # Calculate Solar Angle and correct radiation ############################
#     ##########################################################################
#     # Solar geometry
#     # ----------------------------
#     # SolDec  :declination sun per day  between +23 & -23 [deg]
#     # HourAng :hour angle [-] of sun during day
#     # SolAlt  :solar altitude [deg], height of sun above horizon
#     # SolDec  = -23.4*cos(360*(Day+10)/365);
#     # Now added a new function that should work on all latitudes!
#     theta    =(Day-1)*2 * pi/365  # day expressed in radians
#
#     SolDec =0.006918-0.399912 * cos(theta)+0.070257 * sin(theta) -   0.006758 * cos(2*theta)+0.000907 * sin(2*theta) -  0.002697 * cos(3*theta)+0.001480 * sin(3*theta)
#
#     HourAng = 180/pi * 15*(Hour-12.01)
#     SolAlt  = scalar(asin(scalar(sin(Lat)*sin(SolDec)+cos(Lat)*cos(SolDec)*cos(HourAng))))
#
#     # Solar azimuth
#     # ----------------------------
#     # SolAzi  :angle solar beams to N-S axes earth [deg]
#     SolAzi = scalar(acos((sin(SolDec)*cos(Lat)-cos(SolDec)* sin(Lat)*cos(HourAng))/cos(SolAlt)))
#     SolAzi = ifthenelse(Hour <= 12, SolAzi, 360 - SolAzi)
#
#     # Surface azimuth
#     # ----------------------------
#     # cosIncident :cosine of angle of incident; angle solar beams to angle surface
#     cosIncident = sin(SolAlt)*cos(Slope)+cos(SolAlt)*sin(Slope)*cos(SolAzi-Aspect)
#
#
#     # Critical angle sun
#     # ----------------------------
#     # HoriAng  :tan maximum angle over DEM in direction sun, 0 if neg
#     # CritSun  :tan of maximum angle in direction solar beams
#     # Shade    :cell in sun 1, in shade 0
#     # NOTE: for a changing DEM in time use following 3 statements and put a #
#     #       for the 4th CritSun statement
#     HoriAng   = cover(horizontan(Altitude,directional(SolAzi)),0)
#     #HoriAng   = horizontan(Altitude,directional(SolAzi))
#     HoriAng   = ifthenelse(HoriAng < 0, scalar(0), HoriAng)
#     CritSun   = ifthenelse(SolAlt > 90, scalar(0), scalar(atan(HoriAng)))
#     Shade   = SolAlt > CritSun
#
#     # Radiation outer atmosphere
#     # ----------------------------
#     OpCorr = Trans**((sqrt(1229+(614*sin(SolAlt))**2) -614*sin(SolAlt))*AtmPcor)    # correction for air masses [-]
#     Sout   = Sc*(1+0.034*cos(360*Day/365)) # radiation outer atmosphere [W/m2]
#     Snor   = Sout*OpCorr                   # rad on surface normal to the beam [W/m2]
#
#     # Radiation at DEM
#     # ----------------------------
#     # Sdir   :direct sunlight on a horizontal surface [W/m2] if no shade
#     # Sdiff  :diffuse light [W/m2] for shade and no shade
#     # Stot   :total incomming light Sdir+Sdiff [W/m2] at Hour
#     # Radiation :avg of Stot(Hour) and Stot(Hour-HourStep)
#     # NOTE: PradM only valid for HourStep & DayStep = 1
#     SdirCor   = ifthenelse(Snor*cosIncident*scalar(Shade)<0,0.0,Snor*cosIncident*scalar(Shade))
#     Sdir   = ifthenelse(Snor*cosIncident<0,0.0,Snor*cosIncident)
#     Sdiff  = ifthenelse(Sout*(0.271-0.294*OpCorr)*sin(SolAlt)<0, 0.0, Sout*(0.271-0.294*OpCorr)*sin(SolAlt))
#     AtmosDiffFrac = ifthenelse(Sdir > 0, Sdiff/Sdir, 1)
#
#     # Stot   = cover(Sdir+Sdiff,windowaverage(Sdir+Sdiff,3));     # Rad [W/m2]
#     Stot   = Sdir + Sdiff                                     	    # Rad [W/m2]
#     StotCor   = SdirCor + Sdiff                                   # Rad [W/m2]
#
#
#     return scalar(SolAlt)
#
#
# def GenRadMaps(SaveDir,SaveName,Lat,Lon,Slope,Aspect,Altitude,debug):
#     """ Generates daily radiation maps for a whole year.
#     It does so by running correctrad for a whole year with hourly
#     steps and averaging this per day."""
#
#
#     for Day in range(1,366):
#     	avgrad = 0.0 * Altitude
#     	for Hour in range(6,19):
#     		print "day: " + str(Day) + " Hour: " + str(Hour)
#     		crad = correctrad(Day,Hour,Lat,Lon,Slope,Aspect,Altitude)
#     		if debug:
#     			nr = "%0.3d" % Hour
#     			report(crad,SaveDir + "/" + str(Day) + "/" + SaveName + "00000." + nr)
#     		avgrad=avgrad + crad
#
#     	nr = "%0.3d" % Day
#     	report(avgrad/13.0,SaveDir + "/" + SaveName + "00000." + nr)
#



    
     

    