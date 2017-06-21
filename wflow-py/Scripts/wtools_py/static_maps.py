# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 12:11:43 2015

@author: winsemi

$Id: static_maps.py 1125 2015-11-18 16:28:25Z winsemi $
$Date: 2015-11-18 23:28:25 +0700 (Wed, 18 Nov 2015) $
$Author: winsemi $
$Revision: 1125 $
$HeadURL: https://repos.deltares.nl/repos/Hydrology/trunk/hydro-earth/wtools/scripts/static_maps.py $
$Keywords: $

"""

# import sys packages
from subprocess import call
import sys
import os
import shutil
import glob

# import admin packages
from optparse import OptionParser

# import general packages
import numpy as np
from osgeo import osr, gdal, ogr, gdalconst
from lxml import etree
from osgeo import gdalconst

# import specific packages
from hydrotools import gis
import pcraster as pcr
import wtools_lib

driver = ogr.GetDriverByName("ESRI Shapefile")


def clip_catchment_by_cell(cell_geom, catchment_geom):
    return cell_geom.intersection(catchment_geom)


def main():

    ### Read input arguments #####
    logfilename = 'wtools_static_maps.log'
    parser = OptionParser()
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option('-q', '--quiet',
                      dest='verbose', default=True, action='store_false',
                      help='do not print status messages to stdout')
    parser.add_option('-i', '--ini', dest='inifile', default=None,
                      help='ini file with settings for static_maps.exe')
    parser.add_option('-s', '--source',
                      dest='source', default='wflow',
                      help='Source folder containing clone (default=./wflow)')
    parser.add_option('-d', '--destination',
                      dest='destination', default='staticmaps',
                      help='Destination folder (default=./staticmaps)')
    parser.add_option('-r', '--river',
                      dest='rivshp', default=None,
                      help='river network polyline layer (ESRI Shapefile)')
    parser.add_option('-c', '--catchment',
                      dest='catchshp', default=None,
                      help='catchment polygon layer (ESRI Shapefile)')
    parser.add_option('-g', '--gauges',
                      dest='gaugeshp', default=None,
                      help='gauge point layer (ESRI Shapefile)')
    parser.add_option('-D', '--dem',
                      dest='dem_in', default=None,
                      help='digital elevation model (GeoTiff)')
    parser.add_option('-L', '--landuse',
                      dest='landuse', default=None,
                      help='land use / land cover layer (GeoTiff)')
    parser.add_option('-S', '--soiltype',
                      dest='soil', default=None,
                      help='soil type layer (GeoTiff)')
    parser.add_option('-V', '--vegetation',
                      dest='lai', default=None,
                      help='vegetation LAI layer location (containing 12 GeoTiffs <LAI00000.XXX.tif>)')
    parser.add_option('-O', '--other_maps',
                      dest='other_maps', default=None,
                      help='bracketed [] comma-separated list of paths to other maps that should be reprojected')
    parser.add_option('-C', '--clean',
                      dest='clean', default=False, action='store_true',
                      help='Clean the .xml files from static maps folder when finished')
    parser.add_option('-A', '--alltouch',
                      dest='alltouch', default=False, action='store_true',
                      help='option to burn catchments "all touching".\nUseful when catchment-size is small compared to cellsize')
    (options, args) = parser.parse_args()
    # parse other maps into an array
    options.other_maps = options.other_maps.replace(
        ' ', '').replace('[', '').replace(']', '').split(',')

    options.source = os.path.abspath(options.source)
    clone_map = os.path.join(options.source, 'mask.map')
    clone_shp = os.path.join(options.source, 'mask.shp')
    clone_prj = os.path.join(options.source, 'mask.prj')

    if None in (options.inifile,
                options.rivshp,
                options.catchshp,
                options.dem_in):
        msg = """The following files are compulsory:
        - ini file
        - DEM (raster)
        - river (shape)
        - catchment (shape)
        """
        print(msg)
        parser.print_help()
        sys.exit(1)
    if not os.path.exists(options.inifile):
        print 'path to ini file cannot be found'
        sys.exit(1)
    if not os.path.exists(options.rivshp):
        print 'path to river shape cannot be found'
        sys.exit(1)
    if not os.path.exists(options.catchshp):
        print 'path to catchment shape cannot be found'
        sys.exit(1)
    if not os.path.exists(options.dem_in):
        print 'path to DEM cannot be found'
        sys.exit(1)

    # open a logger, dependent on verbose print to screen or not
    logger, ch = wtools_lib.setlogger(logfilename, 'WTOOLS', options.verbose)

    # create directories # TODO: check if workdir is still necessary, try to
    # keep in memory as much as possible

    # delete old files (when the source and destination folder are different)
    if np.logical_and(os.path.isdir(options.destination),
                      options.destination is not options.source):
        shutil.rmtree(options.destination)
    if options.destination is not options.source:
        os.makedirs(options.destination)

    # Read mask
    if not(os.path.exists(clone_map)):
        logger.error(
            'Clone file {:s} not found. Please run create_grid first.'.format(clone_map))
        sys.exit(1)
    else:
        # set clone
        pcr.setclone(clone_map)
        # get the extent from clone.tif
        xax, yax, clone, fill_value = gis.gdal_readmap(clone_map, 'GTiff')
        trans = wtools_lib.get_geotransform(clone_map)
        extent = wtools_lib.get_extent(clone_map)
        xmin, ymin, xmax, ymax = extent
        zeros = np.zeros(clone.shape)
        ones = pcr.numpy2pcr(pcr.Scalar, np.ones(clone.shape), -9999)
        # get the projection from clone.tif
        srs = wtools_lib.get_projection(clone_map)
        unit_clone = srs.GetAttrValue('UNIT').lower()

    # READ CONFIG FILE
    # open config-file
    config = wtools_lib.OpenConf(options.inifile)

    # read settings
    snapgaugestoriver = wtools_lib.configget(config, 'settings',
                                             'snapgaugestoriver',
                                             True, datatype='boolean')
    burnalltouching = wtools_lib.configget(config, 'settings',
                                           'burncatchalltouching',
                                           True, datatype='boolean')
    burninorder = wtools_lib.configget(config, 'settings',
                                       'burncatchalltouching',
                                       False, datatype='boolean')
    verticetollerance = wtools_lib.configget(config, 'settings',
                                             'vertice_tollerance',
                                             0.0001, datatype='float')

    ''' read parameters '''
    burn_outlets = wtools_lib.configget(config, 'parameters',
                                        'burn_outlets', 10000,
                                        datatype='int')
    burn_rivers = wtools_lib.configget(config, 'parameters',
                                       'burn_rivers', 200,
                                       datatype='int')
    burn_connections = wtools_lib.configget(config, 'parameters',
                                            'burn_connections', 100,
                                            datatype='int')
    burn_gauges = wtools_lib.configget(config, 'parameters',
                                       'burn_gauges', 100,
                                       datatype='int')
    minorder = wtools_lib.configget(config, 'parameters',
                                    'riverorder_min', 3,
                                    datatype='int')
    percentiles = np.array(
        config.get('parameters', 'statisticmaps', '0, 100').replace(
            ' ', '').split(','), dtype='float')

    # read the parameters for generating a temporary very high resolution grid
    if unit_clone == 'degree':
        cellsize_hr = wtools_lib.configget(config, 'parameters',
                                           'highres_degree', 0.0005,
                                           datatype='float')
    elif (unit_clone == 'metre') or (unit_clone == 'meter'):
        cellsize_hr = wtools_lib.configget(config, 'parameters',
                                           'highres_metre', 50,
                                           datatype='float')

    cols_hr = int((float(xmax) - float(xmin)) / cellsize_hr + 2)
    rows_hr = int((float(ymax) - float(ymin)) / cellsize_hr + 2)
    hr_trans = (float(xmin), cellsize_hr, float(0),
                float(ymax), 0, -cellsize_hr)
    clone_hr = os.path.join(options.destination, 'clone_highres.tif')
    # make a highres clone as well!
    wtools_lib.CreateTif(clone_hr, rows_hr, cols_hr, hr_trans, srs, 0)

    # read staticmap locations
    catchment_map = wtools_lib.configget(config, 'staticmaps',
                                         'catchment', 'wflow_catchment.map')
    dem_map = wtools_lib.configget(config, 'staticmaps',
                                   'dem', 'wflow_dem.map')
    demmax_map = wtools_lib.configget(config, 'staticmaps',
                                      'demmax', 'wflow_demmax.map')
    demmin_map = wtools_lib.configget(config, 'staticmaps',
                                      'demmin', 'wflow_demmin.map')
    gauges_map = wtools_lib.configget(config, 'staticmaps',
                                      'gauges', 'wflow_gauges.map')
    landuse_map = wtools_lib.configget(config, 'staticmaps',
                                       'landuse', 'wflow_landuse.map')
    ldd_map = wtools_lib.configget(config, 'staticmaps',
                                   'ldd', 'wflow_ldd.map')
    river_map = wtools_lib.configget(config, 'staticmaps',
                                     'river', 'wflow_river.map')
    outlet_map = wtools_lib.configget(config, 'staticmaps',
                                      'outlet', 'wflow_outlet.map')
    riverlength_fact_map = wtools_lib.configget(config, 'staticmaps',
                                                'riverlength_fact',
                                                'wflow_riverlength_fact.map')
    soil_map = wtools_lib.configget(config, 'staticmaps',
                                    'soil', 'wflow_soil.map')
    streamorder_map = wtools_lib.configget(config, 'staticmaps',
                                           'streamorder',
                                           'wflow_streamorder.map')
    subcatch_map = wtools_lib.configget(config, 'staticmaps',
                                        'subcatch', 'wflow_subcatch.map')

    # read mask location (optional)
    masklayer = wtools_lib.configget(
        config, 'mask', 'masklayer', options.catchshp)

    # ???? empty = pcr.ifthen(ones == 0, pcr.scalar(0))

    # TODO: check if extents are correct this way
    # TODO: check what the role of missing values is in zeros and ones (l. 123
    # in old code)

    # first add a missing value to dem_in
    ds = gdal.Open(options.dem_in, gdal.GA_Update)
    RasterBand = ds.GetRasterBand(1)
    fill_val = RasterBand.GetNoDataValue()

    if fill_val is None:
        RasterBand.SetNoDataValue(-9999)
    ds = None

    # reproject to clone map: see http://stackoverflow.com/questions/10454316/how-to-project-and-resample-a-grid-to-match-another-grid-with-gdal-python
    # resample DEM
    logger.info('Resampling dem from {:s} to {:s}'.format(os.path.abspath(
        options.dem_in), os.path.join(options.destination, dem_map)))
    gis.gdal_warp(options.dem_in, clone_map, os.path.join(
        options.destination, dem_map), format='PCRaster', gdal_interp=gdalconst.GRA_Average)
    # retrieve amount of rows and columns from clone
    # TODO: make windowstats applicable to source/target with different projections. This does not work yet.
    # retrieve srs from DEM
    try:
        srs_dem = wtools_lib.get_projection(options.dem_in)
    except:
        logger.warning(
            'No projection found in DEM, assuming WGS 1984 lat long')
        srs_dem = osr.SpatialReference()
        srs_dem.ImportFromEPSG(4326)
    clone2dem_transform = osr.CoordinateTransformation(srs, srs_dem)
    # if srs.ExportToProj4() == srs_dem.ExportToProj4():
    for percentile in percentiles:
        if percentile >= 100:
            logger.info('computing window maximum')
            percentile_dem = os.path.join(
                options.destination, 'wflow_dem_max.map')
        elif percentile <= 0:
            logger.info('computing window minimum')
            percentile_dem = os.path.join(
                options.destination, 'wflow_dem_min.map')
        else:
            logger.info(
                'computing window {:d} percentile'.format(int(percentile)))
            percentile_dem = os.path.join(
                options.destination, 'wflow_dem_{:03d}.map'.format(int(percentile)))

        percentile_dem = os.path.join(
            options.destination, 'wflow_dem_{:03d}.map'.format(int(percentile)))
        stats = wtools_lib.windowstats(options.dem_in, len(yax), len(xax),
                                       trans, srs, percentile_dem, percentile, transform=clone2dem_transform, logger=logger)
#    else:
#        logger.warning('Projections of DEM and clone are different. DEM statistics for different projections is not yet implemented')

    """

    # burn in rivers
    # first convert and clip the river shapefile
    # retrieve river shape projection, if not available assume EPSG:4326
    file_att = os.path.splitext(os.path.basename(options.rivshp))[0]
    ds = ogr.Open(options.rivshp)
    lyr = ds.GetLayerByName(file_att)
    extent = lyr.GetExtent()
    extent_in = [extent[0], extent[2], extent[1], extent[3]]
    try:
        # get spatial reference from shapefile
        srs_rivshp = lyr.GetSpatialRef()
        logger.info('Projection in river shapefile is {:s}'.format(srs_rivshp.ExportToProj4()))
    except:
        logger.warning('No projection found in {:s}, assuming WGS 1984 lat-lon'.format(options.rivshp))
        srs_rivshp = osr.SpatialReference()
        srs_rivshp.ImportFromEPSG(4326)
    rivprojshp = os.path.join(options.destination, 'rivshp_proj.shp')
    logger.info('Projecting and clipping {:s} to {:s}'.format(options.rivshp, rivprojshp))
    # TODO: Line below takes a very long time to process, the bigger the shapefile, the more time. How do we deal with this?
    call(('ogr2ogr','-s_srs', srs_rivshp.ExportToProj4(),'-t_srs', srs.ExportToProj4(), '-clipsrc', '{:f}'.format(xmin), '{:f}'.format(ymin), '{:f}'.format(xmax), '{:f}'.format(ymax), rivprojshp, options.rivshp))
    """

    # TODO: BURNING!!

    # project catchment layer to projection of clone
    file_att = os.path.splitext(os.path.basename(options.catchshp))[0]
    print options.catchshp
    ds = ogr.Open(options.catchshp)
    lyr = ds.GetLayerByName(file_att)
    extent = lyr.GetExtent()
    extent_in = [extent[0], extent[2], extent[1], extent[3]]
    try:
        # get spatial reference from shapefile
        srs_catchshp = lyr.GetSpatialRef()
        logger.info('Projection in catchment shapefile is {:s}'.format(
            srs_catchshp.ExportToProj4()))
    except:
        logger.warning(
            'No projection found in {:s}, assuming WGS 1984 lat-lon'.format(options.catchshp))
        srs_catchshp = osr.SpatialReference()
        srs_catchshp.ImportFromEPSG(4326)
    catchprojshp = os.path.join(options.destination, 'catchshp_proj.shp')
    logger.info('Projecting {:s} to {:s}'.format(
        options.catchshp, catchprojshp))
    call(('ogr2ogr', '-s_srs', srs_catchshp.ExportToProj4(), '-t_srs', srs.ExportToProj4(), '-clipsrc',
          '{:f}'.format(xmin), '{:f}'.format(ymin), '{:f}'.format(xmax), '{:f}'.format(ymax), catchprojshp, options.catchshp))

    #
    logger.info('Calculating ldd')
    ldddem = pcr.readmap(os.path.join(options.destination, dem_map))
    ldd_select = pcr.lddcreate(ldddem, 1e35, 1e35, 1e35, 1e35)
    pcr.report(ldd_select, os.path.join(options.destination, 'wflow_ldd.map'))

    # compute stream order, identify river cells
    streamorder = pcr.ordinal(pcr.streamorder(ldd_select))
    river = pcr.ifthen(streamorder >= pcr.ordinal(minorder), pcr.boolean(1))
    # find the minimum value in the DEM and cover missing values with a river with this value. Effect is none!! so now left out!
    # mindem = int(np.min(pcr.pcr2numpy(pcr.ordinal(os.path.join(options.destination, dem_map)),9999999)))
    # dem_resample_map = pcr.cover(os.path.join(options.destination, dem_map), pcr.scalar(river)*0+mindem)
    # pcr.report(dem_resample_map, os.path.join(options.destination, dem_map))
    pcr.report(streamorder, os.path.join(options.destination, streamorder_map))
    pcr.report(river, os.path.join(options.destination, river_map))

    # deal with your catchments
    if options.gaugeshp == None:
        logger.info('No gauges defined, using outlets instead')
        gauges = pcr.ordinal(
            pcr.uniqueid(
                pcr.boolean(
                    pcr.ifthen(pcr.scalar(ldd_select) == 5,
                               pcr.boolean(1)
                               )
                )
            )
        )
        pcr.report(gauges, os.path.join(options.destination, gauges_map))
    # TODO: Add the gauge shape code from StaticMaps.py (line 454-489)
    # TODO: add river length map (see SticMaps.py, line 492-499)

    # report river length
    # make a high resolution empty map
    dem_hr_file = os.path.join(options.destination, 'dem_highres.tif')
    burn_hr_file = os.path.join(options.destination, 'burn_highres.tif')
    demburn_hr_file = os.path.join(options.destination, 'demburn_highres.map')
    riv_hr_file = os.path.join(options.destination, 'riv_highres.map')
    gis.gdal_warp(options.dem_in, clone_hr, dem_hr_file)
    # wtools_lib.CreateTif(riv_hr, rows_hr, cols_hr, hr_trans, srs, 0)
    file_att = os.path.splitext(os.path.basename(options.rivshp))[0]
    # open the shape layer
    ds = ogr.Open(options.rivshp)
    lyr = ds.GetLayerByName(file_att)
    gis.ogr_burn(lyr, clone_hr, -100, file_out=burn_hr_file,
                 format='GTiff', gdal_type=gdal.GDT_Float32, fill_value=0)
    # read dem and burn values and add
    xax_hr, yax_hr, burn_hr, fill = gis.gdal_readmap(burn_hr_file, 'GTiff')
    burn_hr[burn_hr == fill] = 0
    xax_hr, yax_hr, dem_hr, fill = gis.gdal_readmap(dem_hr_file, 'GTiff')
    dem_hr[dem_hr == fill] = np.nan
    demburn_hr = dem_hr + burn_hr
    demburn_hr[np.isnan(demburn_hr)] = -9999
    gis.gdal_writemap(demburn_hr_file, 'PCRaster',
                      xax_hr, yax_hr, demburn_hr, -9999.)
    pcr.setclone(demburn_hr_file)
    demburn_hr = pcr.readmap(demburn_hr_file)
    ldd_hr = pcr.lddcreate(demburn_hr, 1e35, 1e35, 1e35, 1e35)
    pcr.report(ldd_hr, os.path.join(options.destination, 'ldd_hr.map'))
    pcr.setglobaloption('unitcell')
    riv_hr = pcr.scalar(pcr.streamorder(ldd_hr) >=
                        minorder) * pcr.downstreamdist(ldd_hr)
    pcr.report(riv_hr, riv_hr_file)
    pcr.setglobaloption('unittrue')
    pcr.setclone(clone_map)
    logger.info('Computing river length')
    #riverlength = wt.windowstats(riv_hr,clone_rows,clone_columns,clone_trans,srs_clone,resultdir,'frac',clone2dem_transform)
    riverlength = wtools_lib.windowstats(riv_hr_file, len(yax), len(xax),
                                         trans, srs, os.path.join(options.destination, riverlength_fact_map), stat='fact', logger=logger)
    # TODO: nothing happends with the river lengths yet. Need to decide how to
    # use these

    # report outlet map
    pcr.report(pcr.ifthen(pcr.ordinal(ldd_select) == 5, pcr.ordinal(1)),
               os.path.join(options.destination, outlet_map))

    # report subcatchment map
    subcatchment = pcr.subcatchment(ldd_select, gauges)
    pcr.report(pcr.ordinal(subcatchment), os.path.join(
        options.destination, subcatch_map))

    # Report land use map
    if options.landuse == None:
        logger.info('No land use map used. Preparing {:s} with only ones.'.
                    format(os.path.join(options.destination, landuse_map)))
        pcr.report(pcr.nominal(ones), os.path.join(
            options.destination, landuse_map))
    else:
        logger.info('Resampling land use from {:s} to {:s}'.
                    format(os.path.abspath(options.landuse),
                           os.path.join(options.destination, os.path.abspath(landuse_map))))
        gis.gdal_warp(options.landuse,
                      clone_map,
                      os.path.join(options.destination, landuse_map),
                      format='PCRaster',
                      gdal_interp=gdalconst.GRA_Mode,
                      gdal_type=gdalconst.GDT_Int32)

    # report soil map
    if options.soil == None:
        logger.info('No soil map used. Preparing {:s} with only ones.'.
                    format(os.path.join(options.destination, soil_map)))
        pcr.report(pcr.nominal(ones), os.path.join(
            options.destination, soil_map))
    else:
        logger.info('Resampling soil from {:s} to {:s}'.
                    format(os.path.abspath(options.soil),
                           os.path.join(options.destination, os.path.abspath(soil_map))))
        gis.gdal_warp(options.soil,
                      clone_map,
                      os.path.join(options.destination, soil_map),
                      format='PCRaster',
                      gdal_interp=gdalconst.GRA_Mode,
                      gdal_type=gdalconst.GDT_Int32)

    if options.lai == None:
        logger.info('No vegetation LAI maps used. Preparing default maps {:s} with only ones.'.
                    format(os.path.join(options.destination, soil_map)))
        pcr.report(pcr.nominal(ones), os.path.join(
            options.destination, soil_map))
    else:
        dest_lai = os.path.join(options.destination, 'clim')
        os.makedirs(dest_lai)
        for month in range(12):
            lai_in = os.path.join(
                options.lai, 'LAI00000.{:03d}'.format(month + 1))
            lai_out = os.path.join(
                dest_lai, 'LAI00000.{:03d}'.format(month + 1))
            logger.info('Resampling vegetation LAI from {:s} to {:s}'.
                        format(os.path.abspath(lai_in),
                               os.path.abspath(lai_out)))
            gis.gdal_warp(lai_in,
                          clone_map,
                          lai_out,
                          format='PCRaster',
                          gdal_interp=gdalconst.GRA_Bilinear,
                          gdal_type=gdalconst.GDT_Float32)

    # report soil map
    if options.other_maps == None:
        logger.info('No other maps used. Skipping other maps.')
    else:
        logger.info('Resampling list of other maps...')
        for map_file in options.other_maps:
            map_name = os.path.split(map_file)[1]
            logger.info('Resampling a map from {:s} to {:s}'.
                        format(os.path.abspath(map_file),
                               os.path.join(options.destination, map_name)))
            gis.gdal_warp(map_file,
                          clone_map,
                          os.path.join(options.destination, map_name),
                          format='PCRaster',
                          gdal_interp=gdalconst.GRA_Mode,
                          gdal_type=gdalconst.GDT_Float32)

    if options.clean:
        wtools_lib.DeleteList(glob.glob(os.path.join(options.destination, '*.xml')),
                              logger=logger)
        wtools_lib.DeleteList(glob.glob(os.path.join(options.destination, 'clim', '*.xml')),
                              logger=logger)
        wtools_lib.DeleteList(glob.glob(os.path.join(options.destination, '*highres*')),
                              logger=logger)


if __name__ == "__main__":
    main()
