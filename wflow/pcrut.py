"""
Collection of utils for working with the
PCRaster python bindings.

"""

import logging
import logging.handlers
import os.path
import math
import sys

import numpy
import pcraster as pcr


def lattometres(lat):
    """"
    Determines the length of one degree lat/long at a given latitude (in meter).
    Code taken from http:www.nga.mil/MSISiteContent/StaticFiles/Calculators/degree.html
    Input: map with lattitude values for each cell
    Returns: length of a cell lat, length of a cell long
    """
    # radlat = pcr.spatial(lat * ((2.0 * math.pi)/360.0))
    # radlat = lat * (2.0 * math.pi)/360.0
    pcr.setglobaloption("degrees")
    radlat = pcr.spatial(lat)  # pcraster cos/sin work in degrees!

    m1 = 111132.92  # latitude calculation term 1
    m2 = -559.82  # latitude calculation term 2
    m3 = 1.175  # latitude calculation term 3
    m4 = -0.0023  # latitude calculation term 4
    p1 = 111412.84  # longitude calculation term 1
    p2 = -93.5  # longitude calculation term 2
    p3 = 0.118  # longitude calculation term 3
    # # Calculate the length of a degree of latitude and longitude in meters

    latlen = (
        m1
        + (m2 * pcr.cos(2.0 * radlat))
        + (m3 * pcr.cos(4.0 * radlat))
        + (m4 * pcr.cos(6.0 * radlat))
    )
    longlen = (
        (p1 * pcr.cos(radlat))
        + (p2 * pcr.cos(3.0 * radlat))
        + (p3 * pcr.cos(5.0 * radlat))
    )

    return latlen, longlen


def detRealCellLength(ZeroMap, sizeinmetres):
    """
    Determine cellength. Always returns the length
    in meters.
    """

    if sizeinmetres:
        reallength = pcr.celllength()
        xl = pcr.celllength()
        yl = pcr.celllength()
    else:
        aa = pcr.ycoordinate(pcr.boolean(pcr.cover(ZeroMap + 1, 1)))
        yl, xl = lattometres(aa)

        xl = xl * pcr.celllength()
        yl = yl * pcr.celllength()
        # Average length for surface area calculations.

        reallength = (xl + yl) * 0.5

    return xl, yl, reallength


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args:
        print(msg)
    print(__doc__)
    sys.exit(0)


def setlogger(logfilename, loggername, thelevel=logging.INFO):
    """
    Set-up the logging system and return a logger object. Exit if this fails
    """

    try:
        # create logger
        logger = logging.getLogger(loggername)
        if not isinstance(thelevel, int):
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(thelevel)
        ch = logging.FileHandler(logfilename, mode="w")
        console = logging.StreamHandler()
        console.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
        # create formatter
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(module)s - %(levelname)s - %(message)s"
        )
        # add formatter to ch
        ch.setFormatter(formatter)
        console.setFormatter(formatter)
        # add ch to logger
        logger.addHandler(ch)
        logger.addHandler(console)
        logger.debug("File logging to " + logfilename)
        return logger
    except IOError:
        print("ERROR: Failed to initialize logger with logfile: " + logfilename)
        sys.exit(2)


def readmapSave(pathtomap, default):
    """
    Adpatation of readmap that returns a default map if the map cannot be found
    """
    if os.path.isfile(pathtomap):
        return pcr.readmap(pathtomap)
    else:
        return pcr.scalar(default)


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
        ifile = open(nname, "r")
        line = ifile.readline()
        line = ifile.readline()
        toskip = int(line) + 2
        line = ifile.readline()
        # read header
        for i in range(toskip - 3):
            line = ifile.readline()
            head.append(line.strip())

        ifile.close()

        mat = numpy.loadtxt(nname, skiprows=toskip)
        mis = mat == 1e31
        mat[mis] = numpy.nan

        if len(mat.shape) > 1:
            return mat[:, 1:], head
            # dumm = mat[:,1:].copy()
            # return numpy.vstack((dumm,mat[:,1:])), head
        else:
            return mat[1:], head
            # dumm = mat[1:].copy()
            # dumm[:] = 1E-31
            # return numpy.vstack((dumm,mat[1:])), head
    else:
        print(nname + " does not exists.")

    return


def interpolategauges(inputmap, method):
    """"
    Interpolate time series gauge data onto a grid using different methods
    inputmap: map with points data for a single timestep
    method: string indicating the method
        inv
        pol
        
    input: inputmap, method
    returns: interpolated map
    """

    if method == "inv":
        result = pcr.inversedistance(1, inputmap, 3, 0, 0)
    elif method == "pol":
        Unq = pcr.uniqueid(pcr.boolean(inputmap + 1))
        result = pcr.spreadzone(pcr.ordinal(pcr.cover(Unq, 0)), 0, 1)
        result = pcr.areaaverage(inputmap, result)
    else:
        Unq = pcr.uniqueid(pcr.boolean(inputmap + 1))
        result = pcr.spreadzone(pcr.ordinal(pcr.cover(Unq, 0)), 0, 1)
        result = pcr.areaaverage(inputmap, result)

    return result


def tableToMapSparse(step, table, map):
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
    global debug

    if not hasattr(tableToMapSparse, "_tableToMap_LastTbl"):
        _tableToMap_LastTbl = {}
        _tableToMap_LastMap = {}

    # construct filenames
    fname_map = map + str(step) + ".map"
    fname_tbl = table + str(step) + ".tbl"

    if os.path.exists(fname_map):
        print("found: " + fname_map)
        tableToMapSparse._tableToMap_LastMap[map] = step

    if os.path.exists(fname_tbl):
        print("found: " + fname_tbl)
        tableToMapSparse._tableToMap_LastTbl[table] = step

    if (table in tableToMapSparse._tableToMap_LastTbl) == False:
        fname_tbl = table + ".tbl"
    else:
        fname_tbl = table + str(tableToMapSparse._tableToMap_LastTbl[table]) + ".tbl"

    if (map in tableToMapSparse._tableToMap_LastMap) == False:
        fname_map = map + ".map"
    else:
        fname_map = map + str(tableToMapSparse._tableToMap_LastMap[map]) + ".map"

    rmat = pcr.lookupscalar(str(fname_tbl), str(fname_map))

    return rmat
