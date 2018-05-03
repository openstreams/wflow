# -*- coding: utf-8 -*-
"""
Created on Wed Jul 08 16:12:29 2015

@author: winsemi

$Id: create_grid.py 1126 2015-11-18 16:28:37Z winsemi $
$Date: 2015-11-18 23:28:37 +0700 (Wed, 18 Nov 2015) $
$Author: winsemi $
$Revision: 1126 $
$HeadURL: https://repos.deltares.nl/repos/Hydrology/trunk/hydro-earth/wtools/scripts/create_grid.py $
$Keywords: $

"""

# import sys packages
import sys
import os
import shutil

# import admin packages
from optparse import OptionParser

# if frozen to exe, import wflow is needed before hydrotools to run __init__
# and set the PROJ_DIR environment variable to the right path

# import general packages
import numpy as np
from osgeo import osr,ogr
from xml.etree import ElementTree
import pyproj
# import specific packages
import wflow.wflowtools_lib as wt


def parse_args():
        ### Read input arguments #####
    parser = OptionParser()
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option('-d', '--destination',
                      dest='destination', default='wflow',
                      help='Destination folder (default=./wflow)')
    parser.add_option('-l', '--logfile',
                      dest='logfilename', default='wtools_create_grid.log',
                      help='log file name')
    parser.add_option('-f', '--file', dest='inputfile',  nargs=1,
                      help='file of which extent will be read. Most logically the catchment layer\nformat: ESRI Shapefile or any gdal supported raster format (preferred GeoTiff)')
    parser.add_option('-p', '--projection',
                      dest='projection', default='EPSG:4326',
                      help='Only used if no file is provided, either of type EPSG:<####> or +proj...')
    parser.add_option('-c', '--cellsize', type='float',
                      nargs=1, dest='cellsize',
                      help='extent')
    parser.add_option('--locationid',
                      dest='locationid', default='wflow_case',
                      help='Sets the name of the locationId in the Delft-FEWS XML grid definition')
    parser.add_option('-s', '--snap',
                      dest='snap', default=False, action='store_true',
                      help='Snaps grid extents to a multiple of the resolution')
    parser.add_option('-q', '--quiet',
                      dest='verbose', default=True, action='store_false',
                      help='do not print status messages to stdout')

    (options, args) = parser.parse_args()

    print options.__dict__.items()

    if options.inputfile is None:
        parser.error('No input file (-f filename) given')
        parser.print_help()
        sys.exit(1)

    if not options.inputfile is None:
        if not os.path.exists(options.inputfile):
            parser.error('input file provided but not found, please check path')
            parser.print_help()
            sys.exit(1)

    if options.cellsize is None:
        parser.error('no cell size (-c cellsize) provided')
        parser.print_help()
        sys.exit(1)

    return options


def main(logfilename,destination,inputfile,projection,cellsize,locationid,snap=False,verbose=True):

    # open a logger, dependent on verbose print to screen or not
    logger, ch = wt.setlogger(logfilename, 'WTOOLS', verbose)

    # delete old files
    if os.path.isdir(destination):
        shutil.rmtree(destination)
    os.makedirs(destination)

    ### Get information ####
    if inputfile is not None:
        # retrieve extent from input file. Check if projection is provided
        file_ext = os.path.splitext(os.path.basename(inputfile))[1]
        if file_ext in ('.shp', '.geojson'):
            ds = ogr.Open(inputfile)
            # read the extent of the shapefile
            lyr = ds.GetLayer(0)
            extent = lyr.GetExtent()
            extent_in = [extent[0], extent[2], extent[1], extent[3]]
            # get spatial reference from shapefile
            srs = lyr.GetSpatialRef()
        else:
            # Read extent from a GDAL compatible file
            extent_in = wt.get_extent(inputfile)
            try:
                extent_in = wt.get_extent(inputfile)
            except:
                msg = 'Input file {:s} not a shape or gdal file'.format(
                    inputfile)
                wt.close_with_error(logger, ch, msg)
                sys.exit(1)

            # get spatial reference from grid file
            try:
                srs = wt.get_projection(inputfile)
            except:
                logger.warning(
                    'No projection found, assuming WGS 1984 lat long')
                srs = osr.SpatialReference()
                srs.ImportFromEPSG(4326)

        # geotransform = ds.GetGeoTransform()
        # raster_cellsize = geotransform[1]
        # ncols = ds.RasterXSize
        # nrows = ds.RasterYSize
        # extent_in = [geotransform[0],
        #              geotransform[3]-nrows*raster_cellsize,
        #              geotransform[0]+ncols*raster_cellsize,
        #              geotransform[3]]
        # # get spatial reference from grid file
        # WktString = ds.GetProjection()
        # srs = osr.SpatialReference()
        # srs.ImportFromWkt(WktString)
    else:
        lonmin, latmin, lonmax, latmax = extent
        srs_4326 = osr.SpatialReference()
        srs_4326.ImportFromEPSG(4326)
        srs = osr.SpatialReference()
        if projection is not None:
            # import projection as an srs object
            if projection.lower()[0:4] == 'epsg':
                # make a proj4 string
                srs.ImportFromEPSG(int(projection[5:]))
            elif projection.lower()[0:5] == '+proj':
                srs.ImportFromProj4(projection)
            else:
                msg = 'Projection "{:s}" is not a valid projection'.format(
                    projection)
                wt.close_with_error(logger, ch, msg)
        else:
            logger.warning('No projection found, assuming WGS 1984 lat long')
            srs.ImportFromEPSG(4326)
        xmin, ymin = pyproj.transform(pyproj.Proj(srs_4326.ExportToProj4()),
                                      pyproj.Proj(srs.ExportToProj4()), lonmin, latmin)
        xmax, ymax = pyproj.transform(pyproj.Proj(srs_4326.ExportToProj4()),
                                      pyproj.Proj(srs.ExportToProj4()), lonmax, latmax)
        # project the extent parameters to selected projection and snap to
        # selected resolution
        extent_in = [xmin, ymin, xmax, ymax]

    # srs known, extent known, prepare UTM or WGS string for grid.xml
    logger.info('Projection "{:s}" used'.format(srs.ExportToProj4()))
    if srs.IsProjected():
        utm = srs.GetUTMZone()
        if utm < 0:
            hemisphere = 'S'
        else:
            hemisphere = 'N'
        geodatum = 'UTM{:d}{:s}'.format(np.abs(utm), hemisphere)
    else:
        geodatum = 'WGS 1984'

    if snap:
        logger.info('Snapping raster')
        snap = len(str(cellsize - np.floor(cellsize))) - 2
        extent_out = wt.round_extent(extent_in, cellsize, snap)
    else:
        extent_out = extent_in
    cols = int((extent_out[2] - extent_out[0]) / cellsize)  # +2)
    rows = int((extent_out[3] - extent_out[1]) / cellsize)  # +2)
    cells = rows * cols
    xorg = extent_out[0]  # -cellsize
    yorg = extent_out[3]  # +cellsize

    # create clone raster
    print('rows: {0} cols: {1}'.format(rows, cols))

    dummy_raster = np.zeros((rows, cols)) - 9999.
    clone_file_map = os.path.abspath(
        os.path.join(destination, 'mask.map'))
    clone_file_tif = os.path.abspath(
        os.path.join(destination, 'mask.tif'))
    logger.info('Writing PCRaster clone to {:s}'.format(clone_file_map))
    wt.gdal_writemap(clone_file_map, 'PCRaster',
                      xorg, yorg, dummy_raster,
                      -9999., resolution=cellsize,
                      srs=srs)
    logger.info('Writing Geotiff clone to {:s}'.format(clone_file_tif))
    wt.gdal_writemap(clone_file_tif, 'GTiff',
                      xorg, yorg, dummy_raster,
                      -9999., resolution=cellsize,
                      zlib=True, srs=srs)

    # create grid.xml
    root = ElementTree.Element('regular', locationId=locationid)
    ElementTree.SubElement(root, 'rows').text = str(rows)
    ElementTree.SubElement(root, 'columns').text = str(cols)
    ElementTree.SubElement(root, 'geoDatum').text = geodatum
    ElementTree.SubElement(root, 'firstCellCenter')
    ElementTree.SubElement(root[3], 'x').text = str(xorg + 0.5 * cellsize)
    ElementTree.SubElement(root[3], 'y').text = str(yorg - 0.5 * cellsize)
    ElementTree.SubElement(root, 'xCellSize').text = str(cellsize)
    ElementTree.SubElement(root, 'yCellSize').text = str(cellsize)
    xml_file = os.path.abspath(os.path.join(destination, 'grid.xml'))
    logger.info('Writing Delft-FEWS grid definition to {:s}'.format(xml_file))
    gridxml = open(xml_file, 'w+')
    gridxml.write(ElementTree.tostring(root))
    gridxml.close()

    # create shape file
    Driver = ogr.GetDriverByName("ESRI Shapefile")
    shp_file = os.path.abspath(os.path.join(destination, 'mask.shp'))
    logger.info('Writing shape of clone to {:s}'.format(shp_file))
    # for encode see https://gis.stackexchange.com/a/53939
    shp_att = os.path.splitext(os.path.basename(shp_file))[0].encode('utf-8')
    shp = Driver.CreateDataSource(shp_file)
    lyr = shp.CreateLayer(shp_att, srs, geom_type=ogr.wkbPolygon)
    fieldDef = ogr.FieldDefn('ID', ogr.OFTString)
    fieldDef.SetWidth(12)
    lyr.CreateField(fieldDef)
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(xorg, yorg)
    ring.AddPoint(xorg + cols * cellsize,
                  yorg)
    ring.AddPoint(xorg + cols * cellsize,
                  yorg - rows * cellsize)
    ring.AddPoint(xorg, yorg - rows * cellsize)
    ring.AddPoint(xorg, yorg)
    ring.CloseRings
    polygon = ogr.Geometry(ogr.wkbPolygon)
    polygon.AddGeometry(ring)
    feat_out = ogr.Feature(lyr.GetLayerDefn())
    feat_out.SetGeometry(polygon)
    feat_out.SetField('ID', 'wflow_mask')
    lyr.CreateFeature(feat_out)
    shp.Destroy()
    logger.info('Model contains {:d} cells'.format(cells))
    if cells > 5000000:
        logger.warning(
            'With this amount of cells your model will run VERY slow.\nConsider a larger cell-size.\nFast models run with < 1,000,000 cells')
    elif cells > 1000000:
        logger.warning(
            'With this amount of cells your model will run slow.\nConsider a larger cell-size. Fast models run with < 1,000,000 cells')
    logger, ch = wt.closeLogger(logger, ch)
    del logger, ch


if __name__ == "__main__":
    argdict = parse_args()
    main(**vars(argdict))
