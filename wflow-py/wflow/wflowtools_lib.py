# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 16:07:41 2014

@author: tollena
"""

import numpy as np
import ConfigParser
from subprocess import call
import copy
import pcraster as pcr
import sys
import logging
import logging.handlers

from osgeo import ogr
from osgeo import gdal
from osgeo.gdalconst import *
from osgeo import osr
import os
import numpy as np

Driver = ogr.GetDriverByName("ESRI Shapefile")


def setlogger(logfilename, logReference, verbose=True):
    """
    Set-up the logging system. Exit if this fails
    """
    try:
        # create logger
        logger = logging.getLogger(logReference)
        logger.setLevel(logging.DEBUG)
        ch = logging.handlers.RotatingFileHandler(
            logfilename, maxBytes=10 * 1024 * 1024, backupCount=5)
        console = logging.StreamHandler()
        console.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
        # create formatter
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        # add formatter to ch
        ch.setFormatter(formatter)
        console.setFormatter(formatter)
        # add ch to logger
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


def OpenConf(fn):
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


def configset(config, section, var, value, overwrite=False):
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
        config.set(section, var, value)
    else:
        if not config.has_option(section, var):
            config.set(section, var, value)
        else:
            if overwrite:
                config.set(section, var, value)


def get_geotransform(filename):
    ''' Return geotransform of dataset'''
    ds = gdal.Open(filename, GA_ReadOnly)
    gt = ds.GetGeoTransform()
    return gt


def get_extent(filename):
    ''' Return list of corner coordinates from a dataset'''
    ds = gdal.Open(filename, GA_ReadOnly)
    gt = ds.GetGeoTransform()
    # 'top left x', 'w-e pixel resolution', '0', 'top left y', '0', 'n-s pixel resolution (negative value)'
    nx, ny = ds.RasterXSize, ds.RasterYSize
    xmin = np.float64(gt[0])
    ymin = np.float64(gt[3]) + np.float64(ny) * np.float64(gt[5])
    xmax = np.float64(gt[0]) + np.float64(nx) * np.float64(gt[1])
    ymax = np.float64(gt[3])
    ds = None
    return xmin, ymin, xmax, ymax


def get_projection(filename):
    ds = gdal.Open(filename, GA_ReadOnly)
    WktString = ds.GetProjection()
    srs = osr.SpatialReference()
    srs.ImportFromWkt(WktString)
    ds = None
    return srs


def round_extent(extent, snap, prec):
    """Increases the extent until all sides lie on a coordinate
    divideable by 'snap'."""
    xmin, ymin, xmax, ymax = extent
    snap = float(snap)  # prevent integer division issues
    xmin = round(np.ceil(xmin / snap) * snap, prec)
    ymin = round(np.ceil(ymin / snap) * snap, prec)
    xmax = round(np.floor(xmax / snap) * snap, prec)
    ymax = round(np.floor(ymax / snap) * snap, prec)
    return xmin, ymin, xmax, ymax


def DeleteShapes(shapes):
    shapelist = list(shapes)
    for shape in shapelist:
        if os.path.exists(shape):
            Driver.DeleteDataSource(shape)
            print "shapefile deleted: " + shape


def MergeShapes(shapesin, Layer):
    for SHP in shapesin:
        if os.path.exists(SHP):
            ATT = os.path.splitext(os.path.basename(SHP))[0]
            DATA = ogr.Open(SHP)
            LYR = DATA.GetLayerByName(ATT)
            LYR.ResetReading()
            for idx, i in enumerate(range(LYR.GetFeatureCount())):
                oldfeature = LYR.GetFeature(i)
                geometry = oldfeature.geometry()
                feature = ogr.Feature(Layer.GetLayerDefn())
                feature.SetGeometry(geometry)
                feature.SetField("ID", oldfeature.GetFieldAsString(0))
                Layer.CreateFeature(feature)
            DATA.Destroy()


def readMap(fileName, fileFormat):
    """
    Read geographical file into memory
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
    resX = geotrans[1]
    resY = geotrans[5]
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x = np.linspace(originX + resX / 2, originX +
                    resX / 2 + resX * (cols - 1), cols)
    y = np.linspace(originY + resY / 2, originY +
                    resY / 2 + resY * (rows - 1), rows)
    # Retrieve raster
    RasterBand = ds.GetRasterBand(1)  # there's only 1 band, starting from 1
    data = RasterBand.ReadAsArray(0, 0, cols, rows)
    FillVal = RasterBand.GetNoDataValue()
    RasterBand = None
    ds = None
    return x, y, data, FillVal


def Reach2Nodes(SHP, EPSG, toll, storedir):
    srs = osr.SpatialReference()
    if not EPSG == None:
        srs.ImportFromEPSG(int(EPSG))
    SHP_ATT = os.path.splitext(os.path.basename(SHP))[0]
    END_SHP = storedir + SHP_ATT + "_end.shp"
    START_SHP = storedir + SHP_ATT + "_start.shp"
    CONN_SHP = storedir + SHP_ATT + "_connection.shp"

    DeleteShapes([END_SHP, START_SHP, CONN_SHP])

    fieldDef = ogr.FieldDefn("ID", ogr.OFTString)
    fieldDef.SetWidth(12)

    END_out = Driver.CreateDataSource(END_SHP)
    END_ATT = os.path.splitext(os.path.basename(END_SHP))[0]
    if not EPSG == None:
        END_LYR = END_out.CreateLayer(END_ATT, srs, geom_type=ogr.wkbPoint)
    else:
        END_LYR = END_out.CreateLayer(END_ATT, geom_type=ogr.wkbPoint)
    END_LYR.CreateField(fieldDef)

    START_out = Driver.CreateDataSource(START_SHP)
    START_ATT = os.path.splitext(os.path.basename(START_SHP))[0]
    if not EPSG == None:
        START_LYR = START_out.CreateLayer(
            START_ATT, srs, geom_type=ogr.wkbPoint)
    else:
        START_LYR = START_out.CreateLayer(START_ATT, geom_type=ogr.wkbPoint)
    START_LYR.CreateField(fieldDef)

    CONN_out = Driver.CreateDataSource(CONN_SHP)
    CONN_ATT = os.path.splitext(os.path.basename(CONN_SHP))[0]
    if not EPSG == None:
        CONN_LYR = CONN_out.CreateLayer(CONN_ATT, srs, geom_type=ogr.wkbPoint)
    else:
        CONN_LYR = CONN_out.CreateLayer(CONN_ATT, geom_type=ogr.wkbPoint)
    CONN_LYR.CreateField(fieldDef)

    StartCoord = []
    EndCoord = []

    ATT = os.path.splitext(os.path.basename(SHP))[0]
    DATA = ogr.Open(SHP)
    LYR = DATA.GetLayerByName(ATT)
    LYR.ResetReading()
    for i in range(LYR.GetFeatureCount()):
        feature = LYR.GetFeature(i)
        geometry = feature.geometry()
        StartCoord.append([geometry.GetX(0), geometry.GetY(0)])
        points = geometry.GetPointCount()
        EndCoord.append([geometry.GetX(points - 1), geometry.GetY(points - 1)])

    DATA.Destroy()

    Connections = [[np.nan, np.nan], [np.nan, np.nan]]

    for i in range(len(StartCoord)):
        if not True in np.all(np.isclose(StartCoord[i], EndCoord, rtol=0, atol=toll), axis=1):
            #point is startpoint
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(StartCoord[i][0], StartCoord[i][1])
            feature = ogr.Feature(START_LYR.GetLayerDefn())
            feature.SetGeometry(point)
            START_LYR.CreateFeature(feature)
        else:
            # point is a connection
            if not True in np.all(np.isclose(StartCoord[i], Connections, rtol=0, atol=toll), axis=1):
                Connections.append(StartCoord[i])
                point = ogr.Geometry(ogr.wkbPoint)
                point.AddPoint(StartCoord[i][0], StartCoord[i][1])
                feature = ogr.Feature(CONN_LYR.GetLayerDefn())
                feature.SetGeometry(point)
                CONN_LYR.CreateFeature(feature)
        if not True in np.all(np.isclose(EndCoord[i], StartCoord, rtol=0, atol=toll), axis=1):
            #point is end
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(EndCoord[i][0], EndCoord[i][1])
            feature = ogr.Feature(END_LYR.GetLayerDefn())
            feature.SetGeometry(point)
            END_LYR.CreateFeature(feature)
        else:
            # point is a connection
            if not True in np.all(np.isclose(EndCoord[i], Connections, rtol=0, atol=toll), axis=1):
                Connections.append(EndCoord[i])
                point = ogr.Geometry(ogr.wkbPoint)
                point.AddPoint(EndCoord[i][0], EndCoord[i][1])
                feature = ogr.Feature(CONN_LYR.GetLayerDefn())
                feature.SetGeometry(point)
                CONN_LYR.CreateFeature(feature)

    END_out.Destroy()
    START_out.Destroy()
    CONN_out.Destroy()
    return START_SHP, END_SHP, CONN_SHP


def ReachOrder(SHP, EPSG, toll, storedir):
    if not EPSG == None:
        srs.ImportFromEPSG(int(EPSG))
    SHP_ATT = os.path.splitext(os.path.basename(SHP))[0]
    ORDER_SHP = storedir + SHP_ATT + "_order1.shp"
    FileOrder = []
    FileOrder.append(ORDER_SHP)

    if os.path.exists(ORDER_SHP):
        Driver.DeleteDataSource(ORDER_SHP)

    IDField = ogr.FieldDefn("ID", ogr.OFTString)
    IDField.SetWidth(12)

    OrderField = ogr.FieldDefn("ORDER", ogr.OFTString)
    OrderField.SetWidth(12)

    StartCoord = []
    EndCoord = []

    ATT = os.path.splitext(os.path.basename(SHP))[0]
    DATA = ogr.Open(SHP)
    LYR = DATA.GetLayerByName(ATT)
    LYR.ResetReading()
    for i in range(LYR.GetFeatureCount()):
        feature = LYR.GetFeature(i)
        geometry = feature.geometry()
        StartCoord.append([geometry.GetX(0), geometry.GetY(0)])
        points = geometry.GetPointCount()
        EndCoord.append([geometry.GetX(points - 1), geometry.GetY(points - 1)])
    EndCoord_np = np.array(EndCoord)
    reaches = copy.deepcopy(i) + 1
    ReachIDs = range(reaches)
    ReachOrders = np.array([None] * len(ReachIDs))

    order = 1
    ordercoord = copy.deepcopy(EndCoord)
    tempcoord = []
    endpoints = 0
    for i in ReachIDs:
        ReachStart = StartCoord[i]
        ReachEnd = EndCoord[i]
        if not True in np.all(np.isclose(ReachStart, ordercoord, rtol=0, atol=toll), axis=1):
            ReachOrders[i] = order
            tempcoord.append(ReachEnd)
        if not True in np.all(np.isclose(ReachEnd, StartCoord, rtol=0, atol=toll), axis=1):
            endpoints += 1
    ordercoord = copy.deepcopy(tempcoord)
    order += 1
    iterations = 0

    while None in list(ReachOrders) and iterations <= 100:
        iterations += 1
        OrderMove = False
        for i in ReachIDs:
            if ReachOrders[i] == None:
                ReachStart = StartCoord[i]
                ReachEnd = EndCoord[i]
                OrderSelect = ReachOrders[np.all(np.isclose(
                    ReachStart, EndCoord_np, rtol=0, atol=toll), axis=1)]
                if not None in list(OrderSelect):
                    if all(x == list(OrderSelect)[0] for x in list(OrderSelect)) == True:
                        OrderMove = True
                        ReachOrders[i] = order
                    else:
                        ReachOrders[i] = int(np.max(OrderSelect))
        if OrderMove:
            order += 1

    if None in list(ReachOrders):
        print "Conversion of river to orders failed. Try to use a smaller tollerance"
        # sys.exit(1)

    LYR.ResetReading()

    for i in range(1, max(ReachOrders) + 1):
        order = i
        ORDER_SHP = storedir + SHP_ATT + "_order" + str(order) + ".shp"
        FileOrder.append(ORDER_SHP)
        if os.path.exists(ORDER_SHP):
            Driver.DeleteDataSource(ORDER_SHP)
        ORDER_out = Driver.CreateDataSource(ORDER_SHP)
        ORDER_ATT = os.path.splitext(os.path.basename(ORDER_SHP))[0]
        if not EPSG == None:
            ORDER_LYR = ORDER_out.CreateLayer(
                ORDER_ATT, srs, geom_type=ogr.wkbLineString)
        else:
            ORDER_LYR = ORDER_out.CreateLayer(
                ORDER_ATT, geom_type=ogr.wkbLineString)
        ORDER_LYR.CreateField(IDField)
        ORDER_LYR.CreateField(OrderField)
        for j in range(LYR.GetFeatureCount()):
            if ReachOrders[j] == order:
                orgfeature = LYR.GetFeature(j)
                geometry = orgfeature.geometry()
                feature = ogr.Feature(ORDER_LYR.GetLayerDefn())
                feature.SetGeometry(geometry)
                feature.SetField("ID", str(j))
                feature.SetField("ORDER", str(order))
                ORDER_LYR.CreateFeature(feature)
        ORDER_out.Destroy()

    DATA.Destroy()
    return FileOrder


def Burn2Tif(shapes, attribute, TIFF):
    for shape in shapes:
        shape_att = os.path.splitext(os.path.basename(shape))[0]
        os.system("gdal_rasterize -a " + str(attribute) +
                  " -l " + shape_att + " " + shape + " " + TIFF)


def ReverseMap(MAP):
    MAX = int(np.max(pcr.pcr2numpy(MAP, np.NAN)))
    REV_MAP = pcr.ordinal(pcr.ifthen(pcr.scalar(
        MAP) == pcr.scalar(-9999), pcr.scalar(0)))
    for i in range(MAX + 1):
        if i > 0:
            print i
            REV_MAP = pcr.cover(pcr.ifthen(pcr.ordinal(MAP) == pcr.ordinal(
                i), pcr.ordinal(pcr.scalar(MAX + 1) - pcr.scalar(i))), REV_MAP)
    REV_MAP = pcr.cover(REV_MAP, pcr.ordinal(MAP))
    return REV_MAP


def DeleteList(itemlist):
    for item in itemlist:
        print 'file deleted: ' + os.path.basename(item)
        os.remove(item)


def Tiff2Point(TIFF):
    DS = gdal.Open(TIFF, GA_ReadOnly)
    if DS is None:
        print 'Could not open ' + fn
        sys.exit(1)

    cols = DS.RasterXSize
    rows = DS.RasterYSize
    geotransform = DS.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]

    band = DS.GetRasterBand(1)
    NoData = band.GetNoDataValue()

    ATT = os.path.splitext(os.path.basename(TIFF))[0]
    SHP = ATT + ".shp"

    if os.path.exists(SHP):
        Driver.DeleteDataSource(SHP)

    SHP_out = Driver.CreateDataSource(SHP)
    SHP_LYR = SHP_out.CreateLayer(ATT, geom_type=ogr.wkbPoint)

    fieldDef = ogr.FieldDefn("ID", ogr.OFTInteger)
    fieldDef.SetWidth(12)
    SHP_LYR.CreateField(fieldDef)

    fieldDef = ogr.FieldDefn("X", ogr.OFTReal)
    fieldDef.SetWidth(20)
    fieldDef.SetPrecision(5)
    SHP_LYR.CreateField(fieldDef)

    fieldDef = ogr.FieldDefn("Y", ogr.OFTReal)
    fieldDef.SetWidth(20)
    fieldDef.SetPrecision(5)
    SHP_LYR.CreateField(fieldDef)

    for x in range(cols):
        for y in range(rows):
            value = band.ReadAsArray(x, y, 1, 1)
            if not value == NoData:
                xCoord = originX + (0.5 + x) * pixelWidth
                yCoord = originY + (y + 0.5) * pixelHeight
                point = ogr.Geometry(ogr.wkbPoint)
                point.AddPoint(xCoord, yCoord)
                feat_out = ogr.Feature(SHP_LYR.GetLayerDefn())
                feat_out.SetGeometry(point)
                feat_out.SetField("ID", str(int(value[0][0])))
                feat_out.SetField("X", xCoord)
                feat_out.SetField("Y", yCoord)
                SHP_LYR.CreateFeature(feat_out)

    SHP_out.Destroy()
    DS = None


def GridDef(TIFF, XML):
    DS = gdal.Open(TIFF, GA_ReadOnly)
    if DS is None:
        print 'Could not open ' + fn
        sys.exit(1)

    cols = DS.RasterXSize
    rows = DS.RasterYSize
    geotransform = DS.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    DS = None
    Grid_xml = open(XML, 'w+')
    Grid_xml.write('<regular locationId="GRID_NAME">\n')
    Grid_xml.write('\t<rows>' + str(rows) + '</rows>\n')
    Grid_xml.write('\t<columns>' + str(cols) + '</columns>\n')
    Grid_xml.write('\t<geoDatum>GEODATUM</geoDatum>\n')
    Grid_xml.write('\t<firstCellCenter>\n')
    Grid_xml.write('\t\t<x>' + str(originX + 0.5 * pixelWidth) + '</x>\n')
    Grid_xml.write('\t\t<y>' + str(originY + 0.5 * pixelHeight) + '</y>\n')
    Grid_xml.write('\t</firstCellCenter>\n')
    Grid_xml.write('\t<xCellSize>' + str(pixelWidth) + '</xCellSize>\n')
    Grid_xml.write('\t<yCellSize>' + str(pixelWidth) + '</yCellSize>\n')
    Grid_xml.write('</regular>\n')


def PCR_river2Shape(rivermap, drainmap, ordermap, lddmap, SHP_FILENAME, catchmentmap, srs=None):
    #    rivermap = riversid_map
    #    drainmap = drain_map
    #    ordermap = streamorder_map
    #    lddmap = ldd_map
    #    SHP_FILENAME = rivshp
    counter = 0.
    percentage = 0.
    file_att = os.path.splitext(os.path.basename(SHP_FILENAME))[0]
    x, y, riversid, FillVal = readMap(rivermap, 'PCRaster')
    riversid[riversid == FillVal] = -1
    x, y, strahlerorder, FillVal = readMap(ordermap, 'PCRaster')
    strahlerorder[strahlerorder == FillVal] = -1
    x, y, catchment, FillVal = readMap(catchmentmap, 'PCRaster')
    catchment[catchment == FillVal] = -1
    x, y, drain, FillVal = readMap(drainmap, 'PCRaster')
    drain[drain == FillVal] = np.nan
    x, y, ldd, FillVal = readMap(lddmap, 'PCRaster')
    xi, yi = np.meshgrid(x, y)

    # mesh of surrounding pixels
    xi_window, yi_window = np.meshgrid(range(-1, 2), range(-1, 2))
    # mesh of ldd grid values
    ldd_values = np.array([[7, 8, 9], [4, 5, 6], [1, 2, 3]])
    [iiy, iix] = np.where(riversid > 0)
    riverId = riversid[iiy, iix]
    maxRiverId = riverId.max()

    # Create new shapefile
    ogr.UseExceptions()
    ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(SHP_FILENAME)
    layer_line = ds.CreateLayer(file_att, srs, ogr.wkbLineString)

    river_ID = ogr.FieldDefn()
    river_ID.SetName('ORDER')
    river_ID.SetType(ogr.OFTInteger)
    river_ID.SetWidth(6)
    layer_line.CreateField(river_ID)

    river_ID = ogr.FieldDefn()
    river_ID.SetName('CATCHMENT')
    river_ID.SetType(ogr.OFTInteger)
    river_ID.SetWidth(6)
    layer_line.CreateField(river_ID)

    # Create a new line geometry per river segment
    for id in np.arange(1, maxRiverId + 1):
        # for id in range(25,26):
        # print 'Writing line element "' + str(id) + '"'
        y_idx, x_idx = np.where(riversid == id)
        drain_idx = drain[y_idx, x_idx]
        lat_select = yi[y_idx, x_idx]
        lon_select = xi[y_idx, x_idx]
        strahlerorder_select = strahlerorder[y_idx, x_idx]
        catchment_select = catchment[y_idx, x_idx]
        order = drain_idx.argsort()
        lat_select = lat_select[order]
        lon_select = lon_select[order]
        catchment_select = catchment_select[order]
        strahlerorder_select = strahlerorder_select[order]
        line = ogr.Geometry(type=ogr.wkbLineString)
        # add points sequentially to line segment
        for nr in range(0, len(lat_select)):
            #line_latlon.AddPoint(np.float64(lon_select[nr]), np.float64(lat_select[nr]))
            line.AddPoint(np.float64(
                lon_select[nr]), np.float64(lat_select[nr]))
        # now find the point downstream of the last pixel from the ldd, which
        # is connected with the downstream river
        try:
            xi_select = xi[y_idx[order][-1] +
                           yi_window, x_idx[order][-1] + xi_window]
            yi_select = yi[y_idx[order][-1] +
                           yi_window, x_idx[order][-1] + xi_window]
            ldd_at_pos = ldd[y_idx[order][-1], x_idx[order][-1]]
            ldd_y, ldd_x = np.where(ldd_values == ldd_at_pos)
            downstream_y = yi_select[ldd_y, ldd_x]
            downstream_x = xi_select[ldd_y, ldd_x]
            line.AddPoint(np.float64(downstream_x), np.float64(downstream_y))
        except:
            continue
            # most downstream point of segment is on the boundary of the map, so skip this step
            # print 'River segment id: %g is on boundary of the map' % id
        # Add line as a new feature to the shapefiles
        feature = ogr.Feature(feature_def=layer_line.GetLayerDefn())
        feature.SetGeometryDirectly(line)
        feature.SetField('ORDER', int(strahlerorder_select[0]))
        feature.SetField('CATCHMENT', int(catchment_select[0]))
        counter = counter + 1
        if (float(id) / float(maxRiverId)) * 100. > percentage:
            #logger.info(' ' + str(int(percentage)) + '% completed')
            percentage = percentage + 10.
        # print 'Writing polyline ' + str(id) + ' of ' + str(maxRiverId)
        layer_line.CreateFeature(feature)
        # Cleanup
        feature.Destroy()
    ds.Destroy()


def snaptomap(points, mmap):
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
    points = pcr.cover(points, 0)
    # Create unique id map of mmap cells
    unq = pcr.nominal(pcr.cover(pcr.uniqueid(
        pcr.defined(mmap)), pcr.scalar(0.0)))
    # Now fill holes in mmap map with lues indicating the closes mmap cell.
    dist_cellid = pcr.scalar(pcr.spreadzone(unq, 0, 1))
    # Get map with values at location in points with closes mmap cell
    dist_cellid = pcr.ifthenelse(points > 0, dist_cellid, 0)
    # Spread this out
    dist_fill = pcr.spreadzone(pcr.nominal(dist_cellid), 0, 1)
    # Find the new (moved) locations
    npt = pcr.uniqueid(pcr.boolean(pcr.ifthen(dist_fill == unq, unq)))
    # Now recreate the original value in the points maps
    ptcover = pcr.spreadzone(pcr.cover(points, 0), 0, 1)
    # Now get the org point value in the pt map
    nptorg = pcr.ifthen(npt > 0, ptcover)

    return nptorg


def Raster2Pol(rasterin, shapeout, srs=None, ID='ID'):
    #    rasterin = catchments_tif
    #    shapeout = catchments_shp
    #    ID = 'ID'
    sourceRaster = gdal.Open(rasterin)
    band = sourceRaster.GetRasterBand(1)
    NoData = band.GetNoDataValue()
    file_att = os.path.splitext(os.path.basename(shapeout))[0]
    driver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(shapeout):
        driver.DeleteDataSource(shapeout)
    outDatasource = driver.CreateDataSource(shapeout)
    outLayer = outDatasource.CreateLayer(file_att, srs)
    newField = ogr.FieldDefn(ID, ogr.OFTInteger)
    outLayer.CreateField(newField)
    gdal.Polygonize(band, band, outLayer, 0, [], callback=None)
    outDatasource.Destroy()
    sourceRaster = None


def windowstats(rasterin, t_rows, t_columns, t_geotransform, t_srs, resultdir, stat=np.array([50]), transform=False, logger=logging):
    """
    :param rasterin: original DEM data
    :param t_rows: number of rows in the final maps
    :param t_columns: number of columns in the final maps
    :param t_geotransform: geotransform of final maps
    :param t_srs: srs of final maps
    :param resultdir: output folder
    :param stat: percentile to compute
    :param transform:
    :param logger:
    :return:
    """

#    activate this if you want to write a window shapefile
#    if os.path.exists("windows.shp"):
#        Driver.DeleteDataSource("windows.shp")
#    window_out = Driver.CreateDataSource("windows.shp")
#    window_lyr  = window_out.CreateLayer("windows", geom_type=ogr.wkbLineString)
#    fieldDef = ogr.FieldDefn("ID", ogr.OFTString)
#    fieldDef.SetWidth(12)
#    window_lyr.CreateField(fieldDef)

    # read properties of input raster
    # print transform
    #transform = True
    if isinstance(stat, str):
        stat = np.array([stat])
    ds_in = gdal.Open(rasterin, GA_ReadOnly)
    in_trans = ds_in.GetGeoTransform()
    cellsize_in = in_trans[1]
    xorg_in = in_trans[0]
    yorg_in = in_trans[3]
    # read properties of output raster
    cellsize_out = t_geotransform[1]
    xorg = t_geotransform[0]
    yorg = t_geotransform[3]
    origin = (xorg_in, yorg_in)
    # compute statistics to new data set
    band_in = ds_in.GetRasterBand(1)
    nodata = band_in.GetNoDataValue()
    array_out = np.ones((t_rows, t_columns, len(stat)),
                        dtype=np.float) * nodata
    blocks = t_rows * t_columns
    counter = 0
    percentage = 0.
    for row in range(t_rows):
        print 'doing row ' + str(row + 1) + '/' + str(t_rows)
        for col in range(t_columns):
            counter = counter + 1
            if (float(counter) / float(blocks)) * 100. > percentage:
                logger.info(' ' + str(int(percentage)) + '% completed')
                percentage = percentage + 10.
            # determine window boundaries
            xl = xorg + (col * cellsize_out)
            xr = xorg + ((col + 1) * cellsize_out)
            yt = yorg - (row * cellsize_out)
            yb = yorg - ((row + 1) * cellsize_out)
            # print xl, xr, yt, yb
            if not transform:
                coordtl = (xl, yt)
                coordbr = (xr, yb)
            else:
                coordtl = transform.TransformPoint(xl, yt)
                coordbr = transform.TransformPoint(xr, yb)
            xOffset = int((coordtl[0] - origin[0]) / cellsize_in)
            yOffset = int((coordtl[1] - origin[1]) / -cellsize_in)
            xBlock = int((coordbr[0] - coordtl[0]) / cellsize_in + 1)
            yBlock = int((coordtl[1] - coordbr[1]) / cellsize_in + 1)
            data_block = band_in.ReadAsArray(
                xOffset, yOffset, xBlock, yBlock).astype(float)
            data_block[data_block == float(nodata)] = np.nan
            # print data_block
            #data_block = data_block.astype(int)
            # print data_block
            for idx, perc in enumerate(stat):
                if perc == 'frac':
                    array_out[row, col, idx] = np.max(
                        [(np.sum(data_block) * cellsize_in) / cellsize_out, 1])
                elif perc == 'sum':
                    array_out[row, col, idx] = np.sum(data_block)
                else:
                    array_out[row, col, idx] = np.percentile(
                        data_block, int(perc))

#    activate this if you want to write a window shapefile
#            line = ogr.Geometry(ogr.wkbLineString)
#            line.AddPoint(xl,yt)
#            line.AddPoint(xr,yt)
#            line.AddPoint(xr,yb)
#            line.AddPoint(xl,yb)
#            line.AddPoint(xl,yt)
#            feature = ogr.Feature(window_lyr.GetLayerDefn())
#            feature.SetGeometry(line)
#            feature.SetField("ID",str(counter))
#            window_lyr.CreateFeature(feature)

#    activate this if you want to write a window shapefile
#    window_out.Destroy()

    array_out[np.isnan(array_out)] = nodata
    names = []
    # write rasters
    for idx, perc in enumerate(stat):
        if perc == 100:
            name = 'max'
            # print 'computing window maximum'
            name_map = os.path.join(
                resultdir, 'wflow_dem{:s}.map'.format(name))
        elif perc == 0:
            name = 'min'
            # print 'computing window minimum'
            name_map = os.path.join(
                resultdir, 'wflow_dem{:s}.map'.format(name))
        elif perc == 'frac':
            name = 'frac'
            # print 'computing window fraction'
            name_map = os.path.join(resultdir, 'wflow_riverlength_fact.map')
        elif perc == 'sum':
            name = 'sum'
            # print 'computing window sum'
            name_map = os.path.join(resultdir, 'windowsum.map')
        else:
            logger.info('computing window {:d} percentile'.format(int(perc)))
            name_map = os.path.join(
                resultdir, 'wflow_dem{:02d}.map'.format(int(perc)))
        names.append(name)
        # name_tif = 'work\\dem_' + name + '.tif'
        ds_out = gdal.GetDriverByName('MEM').Create(
            '', t_columns, t_rows, 1, GDT_Float32)
        ds_out.SetGeoTransform(t_geotransform)
        ds_out.SetProjection(t_srs.ExportToWkt())
        band_out = ds_out.GetRasterBand(1)
        band_out.WriteArray(array_out[:, :, idx], 0, 0)
        band_out.SetNoDataValue(nodata)
        histogram = band_out.GetDefaultHistogram()
        band_out.SetDefaultHistogram(histogram[0], histogram[1], histogram[3])
        ds_in = None
        gdal.GetDriverByName('PCRaster').CreateCopy(name_map, ds_out, 0)
        ds_out = None

        # call(('gdal_translate','-of','PCRaster','-ot','Float32', name_tif,name_map))
    return names


def CreateTif(TIF, rows, columns, geotransform, srs, fill=-9999):
    ds = gdal.GetDriverByName('GTiff').Create(
        TIF, columns, rows, 1, GDT_Float32)
    ds.SetGeoTransform(geotransform)
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    array = np.ones((rows, columns), dtype=np.float) * fill
    band.WriteArray(array, 0, 0)
    band.SetNoDataValue(-9999)
    ds = None


def GetRasterTranform(rasterin, srsout):
    ds = gdal.Open(rasterin, GA_ReadOnly)
    spatialref = None
    transform = False
    EPSG = False
    spatialref = ds.GetProjection()
    if not spatialref == None:
        srsin = osr.SpatialReference()
        srsin.ImportFromWkt(spatialref)
        srsin.AutoIdentifyEPSG()
        EPSG = 'EPSG:' + srsin.GetAuthorityCode(None)
        transform = osr.CoordinateTransformation(srsin, srsout)
    return transform, EPSG
