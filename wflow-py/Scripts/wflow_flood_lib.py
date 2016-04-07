# -*- coding: utf-8 -*-
"""
Created on Wed Jul 08 16:12:29 2015

@author: winsemi

$Id: wflow_flood_lib.py $
$Date: 2016-04-07 12:05:38 +0200 (Thu, 7 Apr 2016) $
$Author: winsemi $
$Revision: $
$HeadURL:  $
$Keywords: $

"""

import sys
import os
import ConfigParser
import logging
import logging.handlers

import numpy as np

from osgeo import osr, gdal, gdalconst
import pcraster as pcr

def setlogger(logfilename, logReference, verbose=True):
    """
    Set-up the logging system. Exit if this fails
    """
    try:
        #create logger
        logger = logging.getLogger(logReference)
        logger.setLevel(logging.DEBUG)
        ch = logging.handlers.RotatingFileHandler(logfilename,maxBytes=10*1024*1024, backupCount=5)
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
        return logger, ch
    except IOError:
        print "ERROR: Failed to initialize logger with logfile: " + logfilename
        sys.exit(1)

def closeLogger(logger, ch):
    logger.removeHandler(ch)
    ch.flush()
    ch.close()
    return logger, ch 
    
def close_with_error(logger, ch, msg):
    logger.error(msg)
    logger, ch = closeLogger(logger, ch)
    del logger, ch
    sys.exit(1)

def open_conf(fn):
    config = ConfigParser.SafeConfigParser()
    config.optionxform = str

    if os.path.exists(fn):
        config.read(fn)
    else:
        print "Cannot open config file: " + fn
        sys.exit(1)
        
    return config

def configget(config, section, var, default, datatype='str'):
    """
    Gets a string from a config file (.ini) and returns a default value if
    the key is not found. If the key is not found it also sets the value
    with the default in the config-file

    Input:
        - config - python ConfigParser object
        - section - section in the file
        - var - variable (key) to get
        - default - default value
        - datatype='str' - can be set to 'boolean', 'int', 'float' or 'str'

    Returns:
        - value (str, boolean, float or int) - either the value from the config file or the default value
    """
    Def = False
    try:
        if datatype == 'int':
            ret = config.getint(section, var)
        elif datatype == 'float':
            ret = config.getfloat(section, var)
        elif datatype == 'boolean':
            ret = config.getboolean(section, var)
        else:
            ret = config.get(section, var)
    except:
        Def = True
        ret = default
        configset(config, section, var, str(default), overwrite=False)

    default = Def
    return ret

def get_gdal_extent(filename):
    ''' Return list of corner coordinates from a dataset'''
    ds = gdal.Open(filename, gdal.GA_ReadOnly)
    gt = ds.GetGeoTransform()
    # 'top left x', 'w-e pixel resolution', '0', 'top left y', '0', 'n-s pixel resolution (negative value)'
    nx, ny = ds.RasterXSize, ds.RasterYSize
    xmin = np.float64(gt[0])
    ymin = np.float64(gt[3]) +np.float64(ny) * np.float64(gt[5])
    xmax = np.float64(gt[0]) + np.float64(nx) * np.float64(gt[1])
    ymax = np.float64(gt[3])
    ds = None
    return xmin, ymin, xmax, ymax

def get_gdal_geotransform(filename):
    ''' Return geotransform of dataset'''
    ds = gdal.Open(filename, gdal.GA_ReadOnly)
    if ds is None:
        logging.warning('Could not open {:s} Shutting down').format(filename)
        sys.exit(1)
    # Retrieve geoTransform info
    gt = ds.GetGeoTransform()
    ds = None
    return gt

def get_gdal_axes(filename, logging=logging):
    geotrans = get_gdal_geotransform(filename)
    # Retrieve geoTransform info
    originX = geotrans[0]
    originY = geotrans[3]
    resX = geotrans[1]
    resY = geotrans[5]

    ds = gdal.Open(filename, gdal.GA_ReadOnly)
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x = np.linspace(originX+resX/2, originX+resX/2+resX*(cols-1), cols)
    y = np.linspace(originY+resY/2, originY+resY/2+resY*(rows-1), rows)
    ds = None
    return x, y

def get_gdal_fill(filename, logging=logging):
    ds = gdal.Open(filename, gdal.GA_ReadOnly)
    if ds is None:
        logging.warning('Could not open {:s} Shutting down').format(filename)
        sys.exit(1)
    # Retrieve geoTransform info
    geotrans = get_gdal_geotransform(filename)
    # Retrieve geoTransform info
    originX = geotrans[0]
    originY = geotrans[3]
    resX = geotrans[1]
    resY = geotrans[5]

    ds = None
    return fill_value

def get_gdal_projection(filename, logging=logging):
    ds = gdal.Open(filename, gdal.GA_ReadOnly)
    if ds is None:
        logging.warning('Could not open {:s} Shutting down').format(filename)
        sys.exit(1)
    WktString = ds.GetProjection()
    srs = osr.SpatialReference()
    srs.ImportFromWkt(WktString)
    ds = None
    return srs

def get_gdal_rasterband(filename, band=1, logging=logging):
    """

    :param filename: GDAL compatible raster file to read from
    :param band: band number (default=1)
    :param logging: logging object
    :return: gdal dataset object, gdal rasterband object
    """
    ds = gdal.Open(filename)
    if ds is None:
        logging.warning('Could not open {:s} Shutting down').format(filename)
        sys.exit(1)
    # Retrieve geoTransform info
    return ds, ds.GetRasterBand(band)   # there's only 1 band, starting from 1

def prepare_gdal(filename, x, y, format='GTiff', logging=logging,
                 gdal_type=gdal.GDT_Float32, zlib=True, srs=None):
    # prepare geotrans
    xul = x[0] - (x[1] - x[0]) / 2
    xres = x[1] - x[0]
    yul = y[0] + (y[0] - y[1]) / 2
    yres = y[1] - y[0]
    geotrans = [xul, xres, 0, yul, 0, yres]

    gdal.AllRegister()
    driver = gdal.GetDriverByName('GTiff')
    # Processing
    logging.info(str('Preparing file {:s}').format(filename))
    if zlib:
        ds = driver.Create(filename, len(x),
                                     len(y), 1, gdal_type,
                                     ['COMPRESS=DEFLATE'])
    else:
        ds = driver.Create(filename, len(x),
                                     len(y), 1, gdal_type)
    ds.SetGeoTransform(geotrans)
    if srs:
        ds.SetProjection(srs.ExportToWkt())
    # get rasterband entry
    logging.info('Prepared {:s}'.format(filename))

    return ds

def gdal_warp(src_filename, clone_filename, dst_filename, gdal_type=gdalconst.GDT_Float32,
              gdal_interp=gdalconst.GRA_Bilinear, format='GTiff', ds_in=None, override_src_proj=None):
    """
    Equivalent of the gdalwarp executable, commonly used on command line.
    The function prepares from a source file, a new file, that has the same
    extent and projection as a clone file.
    The clone file should contain the correct projection.
    The same projection will then be produced for the target file.
    If the clone does not have a projection, EPSG:4326 (i.e. WGS 1984 lat-lon)
    will be assumed.

    :param src_filename: string - file with data that will be warped
    :param clone_filename: string - containing clone file (with projection information)
    :param dst_filename: string - destination file (will have the same extent/projection as clone)
    :param gdal_type: - data type to use for output file (default=gdalconst.GDT_Float32)
    :param gdal_interp: - interpolation type used (default=gdalconst.GRA_Bilinear)
    :param format: - GDAL data format to return (default='GTiff')
    :return: No parameters returned, instead a file is prepared
    """
    if ds_in is None:
        src = gdal.Open(src_filename, gdalconst.GA_ReadOnly)
    else:
        src = ds_in
    src_proj = src.GetProjection()
    if override_src_proj is not None:
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(override_src_proj)
        src_proj = srs.ExportToWkt()
    src_nodata = src.GetRasterBand(1).GetNoDataValue()
    # replace nodata value temporarily for some other value
    src.GetRasterBand(1).SetNoDataValue(np.nan)
    # We want a section of source that matches this:
    clone_ds = gdal.Open(clone_filename, gdalconst.GA_ReadOnly)
    clone_proj = clone_ds.GetProjection()
    if not clone_proj:
        # assume a WGS 1984 projection
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        clone_proj = srs.ExportToWkt()
    clone_geotrans = clone_ds.GetGeoTransform()
    wide = clone_ds.RasterXSize
    high = clone_ds.RasterYSize
    # Output / destination
    dst_mem = gdal.GetDriverByName('MEM').Create('', wide, high, 1, gdal_type)
    dst_mem.SetGeoTransform(clone_geotrans)
    dst_mem.SetProjection(clone_proj)
    if not(src_nodata is None):
        dst_mem.GetRasterBand(1).SetNoDataValue(src_nodata)


    # Do the work, UUUUUUGGGGGHHHH: first make a nearest neighbour interpolation with the nodata values
    # as actual values and determine which indexes have nodata values. This is needed because there is a bug in
    # gdal.ReprojectImage, nodata values are not included and instead replaced by zeros! This is not ideal and if
    # a better solution comes up, it should be replaced.

    gdal.ReprojectImage(src, dst_mem, src_proj, clone_proj, gdalconst.GRA_NearestNeighbour)
    data = dst_mem.GetRasterBand(1).ReadAsArray(0, 0)
    idx = np.where(data==src_nodata)
    # now remove the dataset
    del data

    # now do the real transformation and replace the values that are covered by NaNs by the missing value
    if not(src_nodata is None):
        src.GetRasterBand(1).SetNoDataValue(src_nodata)

    gdal.ReprojectImage(src, dst_mem, src_proj, clone_proj, gdal_interp)
    data = dst_mem.GetRasterBand(1).ReadAsArray(0, 0)
    data[idx] = src_nodata
    dst_mem.GetRasterBand(1).WriteArray(data, 0, 0)

    if format=='MEM':
        return dst_mem
    else:
        # retrieve numpy array of interpolated values
        # write to final file in the chosen file format
        gdal.GetDriverByName(format).CreateCopy(dst_filename, dst_mem, 0)

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
        stream = pcr.ifthenelse(pcr.accuflux(ldd, 1) >= accuThreshold,
                                pcr.boolean(1), pcr.boolean(0))
    else:
        stream = pcr.boolean(pcr.cover(rivers, 0))

    height_river = pcr.ifthenelse(stream, pcr.ordinal(dem*100), 0)
    if basin is None:
        up_elevation = pcr.scalar(pcr.subcatchment(ldd, height_river))
    else:
        drainage_surf = pcr.ifthen(rivers, pcr.accuflux(ldd, 1))
        weight = 1./pcr.scalar(pcr.spreadzone(pcr.cover(pcr.ordinal(drainage_surf), 0), 0, 0))
        up_elevation = pcr.ifthenelse(basin, pcr.scalar(pcr.subcatchment(ldd, height_river)), pcr.scalar(pcr.spreadzone(height_river, 0, weight)))
        # replace areas outside of basin by a spread zone calculation.
    hand = pcr.max(pcr.scalar(pcr.ordinal(dem*100))-up_elevation, 0)/100
    dist = pcr.ldddist(ldd, stream, 1)

    return hand, dist

def subcatch_stream(ldd, threshold):
    """
    Derive catchments based upon strahler threshold
    Input:
        ldd -- pcraster object direction, local drain directions
        threshold -- integer, strahler threshold, subcatchments ge threshold are
                 derived
    output:
        stream_ge -- pcraster object, streams of strahler order ge threshold
        subcatch -- pcraster object, subcatchments of strahler order ge threshold

    """
    # derive stream order

    stream = pcr.streamorder(ldd)
    stream_ge = pcr.ifthen(stream >= threshold, stream)
    stream_up_sum = pcr.ordinal(pcr.upstream(ldd, pcr.cover(pcr.scalar(stream_ge), 0)))
    # detect any transfer of strahler order, to a higher strahler order.
    transition_strahler = pcr.ifthenelse(pcr.downstream(ldd, stream_ge) != stream_ge, pcr.boolean(1),
                                         pcr.ifthenelse(pcr.nominal(ldd) == 5, pcr.boolean(1), pcr.ifthenelse(pcr.downstream(ldd, pcr.scalar(stream_up_sum)) > pcr.scalar(stream_ge), pcr.boolean(1),
                                                                                           pcr.boolean(0))))

    # make unique ids (write to file)
    transition_unique = pcr.ordinal(pcr.uniqueid(transition_strahler))

    # derive upstream catchment areas (write to file)
    subcatch = pcr.nominal(pcr.subcatchment(ldd, transition_unique))
    return stream_ge, subcatch

def volume_spread(ldd, hand, subcatch, volume, volume_thres=0., area_multiplier=1., iterations=15):
    """
    Estimate 2D flooding from a 1D simulation per subcatchment reach
    Input:
        ldd -- pcraster object direction, local drain directions
        hand -- pcraster object float32, elevation data normalised to nearest drain
        subcatch -- pcraster object ordinal, subcatchments with IDs
        volume -- pcraster object float32, scalar flood volume (i.e. m3 volume outside the river bank within subcatchment)
        volume_thres=0. -- scalar threshold, at least this amount of m3 of volume should be present in a catchment
        area_multiplier=1. -- in case the maps are not in m2, set a multiplier other than 1. to convert
        iterations=15 -- number of iterations to use
    Output:
        inundation -- pcraster object float32, scalar inundation estimate
    """
    #initial values
    pcr.setglobaloption("unittrue")
    dem_min = pcr.areaminimum(hand, subcatch)  # minimum elevation in subcatchments
    # pcr.report(dem_min, 'dem_min.map')
    dem_norm = hand - dem_min
    # pcr.report(dem_norm, 'dem_norm.map')
    # surface of each subcatchment
    surface = pcr.areaarea(subcatch)*area_multiplier
    # pcr.report(surface, 'surface.map')

    error_abs = pcr.scalar(1e10)  # initial error (very high)
    volume_catch = pcr.areatotal(volume, subcatch)
    # pcr.report(volume_catch, 'volume_catch.map')

    depth_catch = volume_catch/surface
    # pcr.report(depth_catch, 'depth_catch.map')

    dem_max = pcr.ifthenelse(volume_catch > volume_thres, pcr.scalar(32.),
                             pcr.scalar(0))  # bizarre high inundation depth
    dem_min = pcr.scalar(0.)
    for n in range(iterations):
        print('Iteration: {:02d}'.format(n + 1))
        #####while np.logical_and(error_abs > error_thres, dem_min < dem_max):
        dem_av = (dem_min + dem_max)/2
        # pcr.report(dem_av, 'dem_av00.{:03d}'.format(n + 1))
        # compute value at dem_av
        average_depth_catch = pcr.areaaverage(pcr.max(dem_av - dem_norm, 0), subcatch)
        # pcr.report(average_depth_catch, 'depth_c0.{:03d}'.format(n + 1))
        error = pcr.cover((depth_catch-average_depth_catch)/depth_catch, depth_catch*0)
        # pcr.report(error, 'error000.{:03d}'.format(n + 1))
        dem_min = pcr.ifthenelse(error > 0, dem_av, dem_min)
        dem_max = pcr.ifthenelse(error <= 0, dem_av, dem_max)
    # error_abs = np.abs(error)  # TODO: not needed probably, remove
    inundation = pcr.max(dem_av - dem_norm, 0)
    return inundation

def gdal_writemap(file_name, file_format, x, y, data, fill_val, zlib=False,
                  gdal_type=gdal.GDT_Float32, resolution=None, srs=None, logging=logging):
    """ Write geographical file from numpy array
    Dependencies are osgeo.gdal and numpy
    Input:
        file_name: -- string: reference path to GDAL-compatible file
        file_format: -- string: file format according to GDAL acronym
        (see http://www.gdal.org/formats_list.html)
        x: -- 1D np-array: x-axis, or (if only one value), top-left x-coordinate
        y: -- 1D np-array: y-axis, or (if only one value), top-left y-coordinate
        data: -- 2D np-array: raster data
        fill_val: -- float: fill value
        --------------------------------
    optional inputs:
        zlib=False: -- boolean: determines if output file should be internally
                        zipped or not
        gdal_type=gdal.GDT_Float32: -- gdal data type to write
        resolution=None: -- resolution of dataset, only needed if x and y are given as upperleft coordinates
        srs=None: -- projection object (imported by osgeo.osr)
    """
    # make the geotransform
    # Give georeferences
    if hasattr(x, '__len__'):
        # x is the full axes
        xul = x[0]-(x[1]-x[0])/2
        xres = x[1]-x[0]
    else:
        # x is the top-left corner
        xul = x
        xres = resolution
    if hasattr(y, '__len__'):
        # y is the full axes
        yul = y[0]+(y[0]-y[1])/2
        yres = y[1]-y[0]
    else:
        # y is the top-left corner
        yul = y
        yres = -resolution
    geotrans = [xul, xres, 0, yul, 0, yres]

    gdal.AllRegister()
    driver1 = gdal.GetDriverByName('GTiff')
    driver2 = gdal.GetDriverByName(file_format)
    # Processing
    temp_file_name = str('{:s}.tif').format(file_name)
    logging.info(str('Writing to temporary file {:s}').format(temp_file_name))
    if zlib:
        TempDataset = driver1.Create(temp_file_name, data.shape[1],
                                     data.shape[0], 1, gdal_type,
                                     ['COMPRESS=DEFLATE'])
    else:
        TempDataset = driver1.Create(temp_file_name, data.shape[1],
                                     data.shape[0], 1, gdal_type)
    TempDataset.SetGeoTransform(geotrans)
    if srs:
        TempDataset.SetProjection(srs.ExportToWkt())
    # get rasterband entry
    TempBand = TempDataset.GetRasterBand(1)
    # fill rasterband with array
    TempBand.WriteArray(data, 0, 0)
    TempBand.FlushCache()
    TempBand.SetNoDataValue(fill_val)
    # Create data to write to correct format (supported by 'CreateCopy')
    logging.info(str('Writing to {:s}').format(file_name))
    if zlib:
        driver2.CreateCopy(file_name, TempDataset, 0, ['COMPRESS=DEFLATE'])
    else:
        driver2.CreateCopy(file_name, TempDataset, 0)
    TempDataset = None
    os.remove(temp_file_name)

def gdal_readmap(file_name, file_format, give_geotrans=False):
    """ Read geographical file into memory
    Dependencies are osgeo.gdal and numpy
    Input:
        file_name: -- string: reference path to GDAL-compatible file
        file_format: -- string: file format according to GDAL acronym
        (see http://www.gdal.org/formats_list.html)
        give_geotrans (default=False): -- return the geotrans and amount of
            cols/rows instead of x, y axis
    Output (if give_geotrans=False):
        x: -- 1D np-array: x-axis
        y: -- 1D np-array: y-axis
        data:           -- 2D np-array: raster data
        fill_val         -- float:       fill value
    Output (if give_geotrans=True):
        geotrans: -- 6-digit list with GDAL geotrans vector
        size: -- 2-digit tuple with (cols, rows)
        data:           -- 2D np-array: raster data
        fill_val         -- float:       fill value
    """
    # Open file for binary-reading
    mapFormat = gdal.GetDriverByName(file_format)
    mapFormat.Register()
    ds = gdal.Open(file_name)
    if ds is None:
        logging.warning('Could not open {:s} Shutting down').format(file_name)
        sys.exit(1)
        # Retrieve geoTransform info
    geotrans = ds.GetGeoTransform()
    originX = geotrans[0]
    originY = geotrans[3]
    resX = geotrans[1]
    resY = geotrans[5]
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x = np.linspace(originX+resX/2, originX+resX/2+resX*(cols-1), cols)
    y = np.linspace(originY+resY/2, originY+resY/2+resY*(rows-1), rows)
    # Retrieve raster
    RasterBand = ds.GetRasterBand(1)   # there's only 1 band, starting from 1
    data = RasterBand.ReadAsArray(0, 0, cols, rows)
    fill_val = RasterBand.GetNoDataValue()
    RasterBand = None
    ds = None
    if give_geotrans==True:
        return geotrans, (ds.RasterXSize, ds.RasterYSize), data, fill_val

    else:
        return x, y, data, fill_val
