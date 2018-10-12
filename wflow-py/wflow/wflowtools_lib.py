# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 16:07:41 2014

@author: tollena
"""

import os
import copy
import sys
from subprocess import call
import configparser
import logging
import logging.handlers

import numpy as np
import rasterio
from rasterio import warp
from osgeo import ogr
from osgeo import gdal, gdalconst
from osgeo.gdalconst import *
from osgeo import osr
import pcraster as pcr

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
            logfilename, maxBytes=10 * 1024 * 1024, backupCount=5
        )
        ch.setLevel(logging.DEBUG)
        # create formatter
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(module)s - %(levelname)s - %(message)s"
        )
        # add formatter to ch
        ch.setFormatter(formatter)
        # add ch to logger
        logger.addHandler(ch)
        logger.debug("File logging to " + logfilename)
        return logger, ch
    except IOError:
        print("ERROR: Failed to initialize logger with logfile: " + logfilename)
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
    config = configparser.ConfigParser()
    config.optionxform = str

    if os.path.exists(fn):
        config.read(fn)
    else:
        print("Cannot open config file: " + fn)
        sys.exit(1)

    return config


def configget(config, section, var, default, datatype="str"):
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
        if datatype == "int":
            ret = config.getint(section, var)
        elif datatype == "float":
            ret = config.getfloat(section, var)
        elif datatype == "boolean":
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
    """ Return geotransform of dataset"""
    ds = gdal.Open(filename, GA_ReadOnly)
    gt = ds.GetGeoTransform()
    return gt


def get_extent(filename):
    """ Return list of corner coordinates from a dataset"""
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
    xmin = round(np.floor(xmin / snap) * snap, prec)
    ymin = round(np.floor(ymin / snap) * snap, prec)
    xmax = round(np.ceil(xmax / snap) * snap, prec)
    ymax = round(np.ceil(ymax / snap) * snap, prec)
    return xmin, ymin, xmax, ymax


def DeleteShapes(shapes):
    shapelist = list(shapes)
    for shape in shapelist:
        if os.path.exists(shape):
            Driver.DeleteDataSource(shape)
            print("shapefile deleted: " + shape)


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
        print("Could not open " + fileName + ". Something went wrong!! Shutting down")
        sys.exit(1)
        # Retrieve geoTransform info
    geotrans = ds.GetGeoTransform()
    originX = geotrans[0]
    originY = geotrans[3]
    resX = geotrans[1]
    resY = geotrans[5]
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x = np.linspace(originX + resX / 2, originX + resX / 2 + resX * (cols - 1), cols)
    y = np.linspace(originY + resY / 2, originY + resY / 2 + resY * (rows - 1), rows)
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
        START_LYR = START_out.CreateLayer(START_ATT, srs, geom_type=ogr.wkbPoint)
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
        if not True in np.all(
            np.isclose(StartCoord[i], EndCoord, rtol=0, atol=toll), axis=1
        ):
            # point is startpoint
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(StartCoord[i][0], StartCoord[i][1])
            feature = ogr.Feature(START_LYR.GetLayerDefn())
            feature.SetGeometry(point)
            START_LYR.CreateFeature(feature)
        else:
            # point is a connection
            if not True in np.all(
                np.isclose(StartCoord[i], Connections, rtol=0, atol=toll), axis=1
            ):
                Connections.append(StartCoord[i])
                point = ogr.Geometry(ogr.wkbPoint)
                point.AddPoint(StartCoord[i][0], StartCoord[i][1])
                feature = ogr.Feature(CONN_LYR.GetLayerDefn())
                feature.SetGeometry(point)
                CONN_LYR.CreateFeature(feature)
        if not True in np.all(
            np.isclose(EndCoord[i], StartCoord, rtol=0, atol=toll), axis=1
        ):
            # point is end
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(EndCoord[i][0], EndCoord[i][1])
            feature = ogr.Feature(END_LYR.GetLayerDefn())
            feature.SetGeometry(point)
            END_LYR.CreateFeature(feature)
        else:
            # point is a connection
            if not True in np.all(
                np.isclose(EndCoord[i], Connections, rtol=0, atol=toll), axis=1
            ):
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
    ReachIDs = list(range(reaches))
    ReachOrders = np.array([None] * len(ReachIDs))

    order = 1
    ordercoord = copy.deepcopy(EndCoord)
    tempcoord = []
    endpoints = 0
    for i in ReachIDs:
        ReachStart = StartCoord[i]
        ReachEnd = EndCoord[i]
        if not True in np.all(
            np.isclose(ReachStart, ordercoord, rtol=0, atol=toll), axis=1
        ):
            ReachOrders[i] = order
            tempcoord.append(ReachEnd)
        if not True in np.all(
            np.isclose(ReachEnd, StartCoord, rtol=0, atol=toll), axis=1
        ):
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
                OrderSelect = ReachOrders[
                    np.all(
                        np.isclose(ReachStart, EndCoord_np, rtol=0, atol=toll), axis=1
                    )
                ]
                if not None in list(OrderSelect):
                    if (
                        all(x == list(OrderSelect)[0] for x in list(OrderSelect))
                        == True
                    ):
                        OrderMove = True
                        ReachOrders[i] = order
                    else:
                        ReachOrders[i] = int(np.max(OrderSelect))
        if OrderMove:
            order += 1

    if None in list(ReachOrders):
        print("Conversion of river to orders failed. Try to use a smaller tollerance")
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
                ORDER_ATT, srs, geom_type=ogr.wkbLineString
            )
        else:
            ORDER_LYR = ORDER_out.CreateLayer(ORDER_ATT, geom_type=ogr.wkbLineString)
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
        os.system(
            "gdal_rasterize -a "
            + str(attribute)
            + " -l "
            + shape_att
            + " "
            + shape
            + " "
            + TIFF
        )


def ReverseMap(MAP):
    MAX = int(np.max(pcr.pcr2numpy(MAP, np.NAN)))
    REV_MAP = pcr.ordinal(
        pcr.ifthen(pcr.scalar(MAP) == pcr.scalar(-9999), pcr.scalar(0))
    )
    for i in range(MAX + 1):
        if i > 0:
            print(i)
            REV_MAP = pcr.cover(
                pcr.ifthen(
                    pcr.ordinal(MAP) == pcr.ordinal(i),
                    pcr.ordinal(pcr.scalar(MAX + 1) - pcr.scalar(i)),
                ),
                REV_MAP,
            )
    REV_MAP = pcr.cover(REV_MAP, pcr.ordinal(MAP))
    return REV_MAP


def DeleteList(itemlist, logger=logging):
    for item in itemlist:
        logger.info("Deleting file: " + item)
        os.remove(item)


def Tiff2Point(TIFF):
    DS = gdal.Open(TIFF, GA_ReadOnly)
    if DS is None:
        print("Could not open " + fn)
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
        print("Could not open " + fn)
        sys.exit(1)

    cols = DS.RasterXSize
    rows = DS.RasterYSize
    geotransform = DS.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    DS = None
    Grid_xml = open(XML, "w+")
    Grid_xml.write('<regular locationId="GRID_NAME">\n')
    Grid_xml.write("\t<rows>" + str(rows) + "</rows>\n")
    Grid_xml.write("\t<columns>" + str(cols) + "</columns>\n")
    Grid_xml.write("\t<geoDatum>GEODATUM</geoDatum>\n")
    Grid_xml.write("\t<firstCellCenter>\n")
    Grid_xml.write("\t\t<x>" + str(originX + 0.5 * pixelWidth) + "</x>\n")
    Grid_xml.write("\t\t<y>" + str(originY + 0.5 * pixelHeight) + "</y>\n")
    Grid_xml.write("\t</firstCellCenter>\n")
    Grid_xml.write("\t<xCellSize>" + str(pixelWidth) + "</xCellSize>\n")
    Grid_xml.write("\t<yCellSize>" + str(pixelWidth) + "</yCellSize>\n")
    Grid_xml.write("</regular>\n")


def PCR_river2Shape(
    rivermap, drainmap, ordermap, lddmap, SHP_FILENAME, catchmentmap, srs=None
):
    #    rivermap = riversid_map
    #    drainmap = drain_map
    #    ordermap = streamorder_map
    #    lddmap = ldd_map
    #    SHP_FILENAME = rivshp
    counter = 0.0
    percentage = 0.0
    file_att = os.path.splitext(os.path.basename(SHP_FILENAME))[0]
    x, y, riversid, FillVal = readMap(rivermap, "PCRaster")
    riversid[riversid == FillVal] = -1
    x, y, strahlerorder, FillVal = readMap(ordermap, "PCRaster")
    strahlerorder[strahlerorder == FillVal] = -1
    x, y, catchment, FillVal = readMap(catchmentmap, "PCRaster")
    catchment[catchment == FillVal] = -1
    x, y, drain, FillVal = readMap(drainmap, "PCRaster")
    drain[drain == FillVal] = np.nan
    x, y, ldd, FillVal = readMap(lddmap, "PCRaster")
    xi, yi = np.meshgrid(x, y)

    # mesh of surrounding pixels
    xi_window, yi_window = np.meshgrid(list(range(-1, 2)), list(range(-1, 2)))
    # mesh of ldd grid values
    ldd_values = np.array([[7, 8, 9], [4, 5, 6], [1, 2, 3]])
    [iiy, iix] = np.where(riversid > 0)
    riverId = riversid[iiy, iix]
    maxRiverId = riverId.max()

    # Create new shapefile
    ogr.UseExceptions()
    ds = ogr.GetDriverByName("ESRI Shapefile").CreateDataSource(SHP_FILENAME)
    layer_line = ds.CreateLayer(file_att, srs, ogr.wkbLineString)

    river_ID = ogr.FieldDefn()
    river_ID.SetName("ORDER")
    river_ID.SetType(ogr.OFTInteger)
    river_ID.SetWidth(6)
    layer_line.CreateField(river_ID)

    river_ID = ogr.FieldDefn()
    river_ID.SetName("CATCHMENT")
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
            # line_latlon.AddPoint(np.float64(lon_select[nr]), np.float64(lat_select[nr]))
            line.AddPoint(np.float64(lon_select[nr]), np.float64(lat_select[nr]))
        # now find the point downstream of the last pixel from the ldd, which
        # is connected with the downstream river
        try:
            xi_select = xi[y_idx[order][-1] + yi_window, x_idx[order][-1] + xi_window]
            yi_select = yi[y_idx[order][-1] + yi_window, x_idx[order][-1] + xi_window]
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
        feature.SetField("ORDER", int(strahlerorder_select[0]))
        feature.SetField("CATCHMENT", int(catchment_select[0]))
        counter = counter + 1
        if (float(id) / float(maxRiverId)) * 100.0 > percentage:
            # logger.info(' ' + str(int(percentage)) + '% completed')
            percentage = percentage + 10.0
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
    unq = pcr.nominal(pcr.cover(pcr.uniqueid(pcr.defined(mmap)), pcr.scalar(0.0)))
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


def Raster2Pol(rasterin, shapeout, srs=None, ID="ID"):
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


def windowstats(
    rasterin,
    t_rows,
    t_columns,
    t_geotransform,
    t_srs,
    resultdir,
    stat=np.array([50]),
    transform=False,
    logger=logging,
):
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
    # transform = True
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
    array_out = np.ones((t_rows, t_columns, len(stat)), dtype=np.float) * nodata
    blocks = t_rows * t_columns
    counter = 0
    percentage = 0.0
    for row in range(t_rows):
        # print 'doing row ' + str(row + 1) + '/' + str(t_rows)
        for col in range(t_columns):
            counter = counter + 1
            if (float(counter) / float(blocks)) * 100.0 > percentage:
                logger.info(" " + str(int(percentage)) + "% completed")
                percentage = percentage + 10.0
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
            data_block = band_in.ReadAsArray(xOffset, yOffset, xBlock, yBlock)
            if data_block is None:
                continue
            data_block = data_block.astype(np.float)
            data_block[data_block == float(nodata)] = np.nan
            # print data_block
            # data_block = data_block.astype(int)
            # print data_block
            for idx, perc in enumerate(stat):
                if perc == "fact":
                    array_out[row, col, idx] = np.max(
                        [(np.sum(data_block) * cellsize_in) / cellsize_out, 1]
                    )
                elif perc == "sum":
                    array_out[row, col, idx] = np.sum(data_block)
                else:
                    array_out[row, col, idx] = np.percentile(data_block, int(perc))

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
    logger.info("writing rasters")
    for idx, perc in enumerate(stat):
        if perc == 100:
            name = "max"
            # print 'computing window maximum'
            name_map = os.path.join(resultdir, "wflow_dem{:s}.map".format(name))
            logger.info("wflow_dem{:s}.map".format(name))
        elif perc == 0:
            name = "min"
            # print 'computing window minimum'
            name_map = os.path.join(resultdir, "wflow_dem{:s}.map".format(name))
            logger.info("wflow_dem{:s}.map".format(name))
        elif perc == "fact":
            name = "fact"
            # print 'computing window fraction'
            name_map = os.path.join(resultdir, "wflow_riverlength_fact.map")
            logger.info("wflow_riverlength_fact.map")
        elif perc == "sum":
            name = "sum"
            # print 'computing window sum'
            name_map = os.path.join(resultdir, "windowsum.map")
            logger.info("wflow_dem{:s}.map".format(name))
        else:
            name_map = os.path.join(resultdir, "wflow_dem{:02d}.map".format(int(perc)))
            logger.info("wflow_dem{:02d}.map".format(int(perc)))
        names.append(name)

        # name_tif = 'work\\dem_' + name + '.tif'
        ds_out = gdal.GetDriverByName("MEM").Create(
            "", t_columns, t_rows, 1, GDT_Float32
        )
        ds_out.SetGeoTransform(t_geotransform)
        ds_out.SetProjection(t_srs.ExportToWkt())
        band_out = ds_out.GetRasterBand(1)
        band_out.WriteArray(array_out[:, :, idx], 0, 0)
        band_out.SetNoDataValue(nodata)
        histogram = band_out.GetDefaultHistogram()
        if histogram is not None:
            band_out.SetDefaultHistogram(histogram[0], histogram[1], histogram[3])
        ds_in = None
        gdal.GetDriverByName("PCRaster").CreateCopy(name_map, ds_out, 0)
        ds_out = None


def CreateTif(TIF, rows, columns, geotransform, srs, fill=-9999):
    ds = gdal.GetDriverByName("GTiff").Create(TIF, columns, rows, 1, GDT_Float32)
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
        EPSG = "EPSG:" + srsin.GetAuthorityCode(None)
        transform = osr.CoordinateTransformation(srsin, srsout)
    return transform, EPSG


def gdal_writemap(
    file_name,
    file_format,
    x,
    y,
    data,
    fill_val,
    zlib=False,
    gdal_type=gdal.GDT_Float32,
    resolution=None,
    srs=None,
    logging=logging,
    metadata=None,
):
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
        metadata=None: -- dictionary of metadata entries (key/value pairs)
    """
    # make the geotransform
    # Give georeferences
    if hasattr(x, "__len__"):
        # x is the full axes
        xul = x[0] - (x[1] - x[0]) / 2
        xres = x[1] - x[0]
    else:
        # x is the top-left corner
        xul = x
        xres = resolution
    if hasattr(y, "__len__"):
        # y is the full axes
        yul = y[0] + (y[0] - y[1]) / 2
        yres = y[1] - y[0]
    else:
        # y is the top-left corner
        yul = y
        yres = -resolution
    geotrans = [xul, xres, 0, yul, 0, yres]

    gdal.AllRegister()
    driver1 = gdal.GetDriverByName("GTiff")
    driver2 = gdal.GetDriverByName(file_format)
    # Processing
    temp_file_name = str("{:s}.tif").format(file_name)
    logging.info(str("Writing to temporary file {:s}").format(temp_file_name))
    if zlib:
        TempDataset = driver1.Create(
            temp_file_name,
            data.shape[1],
            data.shape[0],
            1,
            gdal_type,
            ["COMPRESS=DEFLATE"],
        )
    else:
        TempDataset = driver1.Create(
            temp_file_name, data.shape[1], data.shape[0], 1, gdal_type
        )
    TempDataset.SetGeoTransform(geotrans)
    if srs:
        TempDataset.SetProjection(srs.ExportToWkt())
    # get rasterband entry
    TempBand = TempDataset.GetRasterBand(1)
    # fill rasterband with array
    TempBand.WriteArray(data, 0, 0)
    TempBand.FlushCache()
    TempBand.SetNoDataValue(fill_val)
    if metadata is not None:
        TempDataset.SetMetadata(metadata)

    # Create data to write to correct format (supported by 'CreateCopy')
    logging.info(str("Writing to {:s}").format(file_name))
    if zlib:
        driver2.CreateCopy(file_name, TempDataset, 0, ["COMPRESS=DEFLATE"])
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
        logging.warning("Could not open {:s} Shutting down".format(file_name))
        sys.exit(1)
        # Retrieve geoTransform info
    geotrans = ds.GetGeoTransform()
    originX = geotrans[0]
    originY = geotrans[3]
    resX = geotrans[1]
    resY = geotrans[5]
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x = np.linspace(originX + resX / 2, originX + resX / 2 + resX * (cols - 1), cols)
    y = np.linspace(originY + resY / 2, originY + resY / 2 + resY * (rows - 1), rows)
    # Retrieve raster
    RasterBand = ds.GetRasterBand(1)  # there's only 1 band, starting from 1
    data = RasterBand.ReadAsArray(0, 0, cols, rows)
    fill_val = RasterBand.GetNoDataValue()
    RasterBand = None
    ds = None
    if give_geotrans == True:
        return geotrans, (ds.RasterXSize, ds.RasterYSize), data, fill_val

    else:
        return x, y, data, fill_val


def gdal_warp(
    src_filename,
    clone_filename,
    dst_filename,
    gdal_type=gdalconst.GDT_Float32,
    gdal_interp=gdalconst.GRA_Bilinear,
    format="GTiff",
    ds_in=None,
    override_src_proj=None,
):
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
    if not src_proj:
        # assume a WGS 1984 projection
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        src_proj = srs.ExportToWkt()
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
    dst_mem = gdal.GetDriverByName("MEM").Create("", wide, high, 1, gdal_type)
    dst_mem.SetGeoTransform(clone_geotrans)
    dst_mem.SetProjection(clone_proj)
    if not (src_nodata is None):
        dst_mem.GetRasterBand(1).SetNoDataValue(src_nodata)

    # Do the work, UUUUUUGGGGGHHHH: first make a nearest neighbour interpolation with the nodata values
    # as actual values and determine which indexes have nodata values. This is needed because there is a bug in
    # gdal.ReprojectImage, nodata values are not included and instead replaced by zeros! This is not ideal and if
    # a better solution comes up, it should be replaced.

    gdal.ReprojectImage(
        src, dst_mem, src_proj, clone_proj, gdalconst.GRA_NearestNeighbour
    )
    data = dst_mem.GetRasterBand(1).ReadAsArray(0, 0)
    idx = np.where(data == src_nodata)
    # now remove the dataset
    del data

    # now do the real transformation and replace the values that are covered by NaNs by the missing value
    if not (src_nodata is None):
        src.GetRasterBand(1).SetNoDataValue(src_nodata)

    gdal.ReprojectImage(src, dst_mem, src_proj, clone_proj, gdal_interp)
    data = dst_mem.GetRasterBand(1).ReadAsArray(0, 0)
    data[idx] = src_nodata
    dst_mem.GetRasterBand(1).WriteArray(data, 0, 0)

    if format == "MEM":
        return dst_mem
    else:
        # retrieve numpy array of interpolated values
        # write to final file in the chosen file format
        gdal.GetDriverByName(format).CreateCopy(dst_filename, dst_mem, 0)


# Check if this can fully replace the gdal_warp defined above.
# This initializes with nodata instead of 0.
def warp_like(
    input, output, like, format=None, co={}, resampling=warp.Resampling.nearest
):
    """Warp a raster to lie on top op of an existing dataset.

    This function is meant to be similar to the ``rio warp --like`` CLI,
    and uses code from the implementation. Once https://github.com/mapbox/rasterio/issues/784
    is resolved this function may no longer be required. For further information see
    https://mapbox.github.io/rasterio/cli.html#warp and
    https://mapbox.github.io/rasterio/topics/reproject.html

    Parameters
    ----------
    input : str
        Path to the input raster.
    output : str
        Path to the output raster.
    like : str or rasterio.DatasetReader
        The raster whose affine transform, size, and crs are used.
    format : str, optional
        The GDAL raster driver code of the output, see http://www.gdal.org/formats_list.html
        By default, use the same format as the input.
    co : dict, optional
        Creation options for creating the output.
    resampling : Resampling enum, optional
        Default value is ``rasterio.warp.Resampling.nearest``. For other values see
        https://mapbox.github.io/rasterio/topics/resampling.html?highlight=resampling#resampling-methods
    
    Example
    -------
    >>> warp_like('input.map', 'output.tif', 'like.tif', format='GTiff',
                  co={'COMPRESS':'DEFLATE'}, resampling=warp.Resampling.med)
    """

    dst_crs, dst_transform, dst_height, dst_width = _like(like)

    with rasterio.open(input) as src:
        out_kwargs = src.profile.copy()
        out_kwargs.update(
            {
                "crs": dst_crs,
                "transform": dst_transform,
                "width": dst_width,
                "height": dst_height,
            }
        )

        # else the format is equal to the input format
        if format is not None:
            out_kwargs["driver"] = format

        # Adjust block size if necessary.
        if "blockxsize" in out_kwargs and dst_width < out_kwargs["blockxsize"]:
            del out_kwargs["blockxsize"]
        if "blockysize" in out_kwargs and dst_height < out_kwargs["blockysize"]:
            del out_kwargs["blockysize"]

        out_kwargs.update(co)

        with rasterio.open(output, "w", **out_kwargs) as dst:
            warp.reproject(
                source=rasterio.band(src, list(range(1, src.count + 1))),
                destination=rasterio.band(dst, list(range(1, src.count + 1))),
                src_transform=src.transform,
                src_crs=src.crs,
                src_nodata=src.nodata,
                dst_transform=out_kwargs["transform"],
                dst_crs=out_kwargs["crs"],
                dst_nodata=dst.nodata,
                resampling=resampling,
            )


def _like(src):
    """Get the properties from a raster for warp_like"""
    if isinstance(src, rasterio.DatasetReader):
        return src.crs, src.transform, src.height, src.width
    else:
        with rasterio.open(src) as ds:
            return ds.crs, ds.transform, ds.height, ds.width


def ogr_burn(
    lyr,
    clone,
    burn_value,
    file_out="",
    gdal_type=gdal.GDT_Byte,
    format="MEM",
    fill_value=255,
    attribute=None,
):
    """
    ogr_burn burns polygons, points or lines from a geographical source (e.g. shapefile) onto a raster.
    Inputs:
        lyr:                Shape layer (e.g. read from ogr object) to burn
        clone:              clone file to use to define geotransform
        burn_value:         burn value
        zlib=False:         Set to True (recommended) to internally zip the
                            data
        fill_value=255      Set the fill value
        gdal_type=
        gdal.GDT_Float32:   Set the GDAL output data type.
        format='MEM':       File format (if 'MEM' is used, data is only kept in memory)
        fill_value=255:     fill value to use
        attribute=None:     alternative to burn_value, if set to attribute name, this attribute is used for burning instead of burn_value
    Output:
        The function returns a GDAL-compatible file (default = in-memory) and the numpy array raster
        TO-DO add metadata and projection information to GeoTIFF
    """
    # get geotransform
    ds_src = gdal.Open(clone, gdal.GA_ReadOnly)
    geotrans = ds_src.GetGeoTransform()
    xcount = ds_src.RasterXSize
    ycount = ds_src.RasterYSize
    # get the projection
    WktString = ds_src.GetProjection()
    srs = osr.SpatialReference()
    srs.ImportFromWkt(WktString)
    ds_src = None

    ds = gdal.GetDriverByName(format).Create(file_out, xcount, ycount, 1, gdal_type)
    ds.SetGeoTransform(geotrans)
    ds.SetProjection(srs.ExportToWkt())
    # create for target raster the same projection as for the value raster
    raster_srs = osr.SpatialReference()
    #    raster_srs.ImportFromWkt(raster.GetProjectionRef())
    #    target_ds.SetProjection(raster_srs.ExportToWkt())

    # rasterize zone polygon to raster
    if attribute is None:
        gdal.RasterizeLayer(ds, [1], lyr, burn_values=[burn_value])
    else:
        gdal.RasterizeLayer(ds, [1], lyr, options=["ATTRIBUTE={:s}".format(attribute)])
    band = ds.GetRasterBand(1)

    band.SetNoDataValue(fill_value)
    if format == "MEM":
        return ds
    else:

        band = None
        ds = None
