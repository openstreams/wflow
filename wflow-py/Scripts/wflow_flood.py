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

from optparse import OptionParser
import numpy as np

from hydrotools import gis
from osgeo import gdal, gdalconst
import pcraster as pcr
import netCDF4 as nc

import wflow_flood_lib as inun_lib
import pdb

def main():
    ### Read input arguments #####
    parser = OptionParser()
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option('-q', '--quiet',
                      dest='verbose', default=True, action='store_false',
                      help='do not print status messages to stdout')
    parser.add_option('-i', '--ini', dest='inifile',
                      default='hand_contour_inun.ini', nargs=1,
                      help='ini configuration file')
    parser.add_option('-f', '--flood_map',
                      nargs=1, dest='flood_map',
                      help='Flood map file (NetCDF point time series file')
    parser.add_option('-v', '--flood_variable',
                      nargs=1, dest='flood_variable',
                      default='water_level',
                      help='variable name of flood water level')
    parser.add_option('-b', '--bankfull_map',
                      dest='bankfull_map', default='',
                      help='Map containing bank full level (is subtracted from flood map, in NetCDF)')
    parser.add_option('-c', '--catchment',
                      dest='catchment_strahler', default=7, type='int',
                      help='Strahler order threshold >= are selected as catchment boundaries')
    parser.add_option('-s', '--hand_strahler',
                      dest='hand_strahler', default=7, type='int',
                      help='Strahler order threshold >= selected as riverine')
    parser.add_option('-d', '--destination',
                      dest='dest_path', default='inun',
                      help='Destination path')
    (options, args) = parser.parse_args()

    if not os.path.exists(options.inifile):
        print 'path to ini file cannot be found'
        sys.exit(1)
    options.dest_path = os.path.abspath(options.dest_path)

    # # delete old files
    # if os.path.isdir(options.dest_path):
    #     shutil.rmtree(options.dest_path)
    # os.makedirs(options.dest_path)

    # set up the logger
    flood_name = os.path.split(options.flood_map)[1].split('.')[0]
    case_name = 'inun_{:s}_hand_{:02d}_catch_{:02d}'.format(flood_name, options.hand_strahler, options.catchment_strahler)
    logfilename = os.path.join(options.dest_path, 'hand_contour_inun.log')
    logger, ch = inun_lib.setlogger(logfilename, 'HAND_INUN', options.verbose)
    logger.info('$Id: $')
    logger.info('Flood map: {:s}'.format(options.flood_map))
    logger.info('Bank full map: {:s}'.format(options.bankfull_map))
    logger.info('Destination path: {:s}'.format(options.dest_path))
    # read out ini file
    ### READ CONFIG FILE
    # open config-file
    config = inun_lib.open_conf(options.inifile)
    
    # read settings
    options.dem_file = inun_lib.configget(config, 'maps',
                                  'dem_file',
                                  True)
    options.ldd_file = inun_lib.configget(config, 'maps',
                                'ldd_file',
                                 True)
    options.stream_file = inun_lib.configget(config, 'maps',
                                'stream_file',
                                 True)
    options.riv_length_file = inun_lib.configget(config, 'maps',
                                'riv_length_file',
                                 True)
    options.riv_width_file = inun_lib.configget(config, 'maps',
                                'riv_width_file',
                                 True)
    options.file_format = inun_lib.configget(config, 'maps',
                                'file_format', 0, datatype='int')
    options.x_tile = inun_lib.configget(config, 'tiling',
                                  'x_tile', 10000, datatype='int')
    options.y_tile = inun_lib.configget(config, 'tiling',
                                  'y_tile', 10000, datatype='int')
    options.x_overlap = inun_lib.configget(config, 'tiling',
                                  'x_overlap', 1000, datatype='int')
    options.y_overlap = inun_lib.configget(config, 'tiling',
                                  'y_overlap', 1000, datatype='int')
    options.iterations = inun_lib.configget(config, 'inundation',
                                  'iterations', 20, datatype='int')
    options.initial_level = inun_lib.configget(config, 'inundation',
                                  'initial_level', 32., datatype='float')
    options.area_multiplier = inun_lib.configget(config, 'inundation',
                                  'area_multiplier', 1., datatype='float')
    logger.info('DEM file: {:s}'.format(options.dem_file))
    logger.info('LDD file: {:s}'.format(options.ldd_file))
    logger.info('Columns per tile: {:d}'.format(options.x_tile))
    logger.info('Rows per tile: {:d}'.format(options.y_tile))
    logger.info('Columns overlap: {:d}'.format(options.x_overlap))
    logger.info('Rows overlap: {:d}'.format(options.y_overlap))
    metadata_global = {}
    # add metadata from the section [metadata]
    meta_keys = config.options('metadata_global')
    for key in meta_keys:
        metadata_global[key] = config.get('metadata_global', key)
    # add a number of metadata variables that are mandatory
    metadata_global['config_file'] = os.path.abspath(options.inifile)
    metadata_var = {}
    metadata_var['units'] = 'm'
    metadata_var['standard_name'] = 'water_surface_height_above_reference_datum'
    metadata_var['long_name'] = 'Coastal flooding'
    metadata_var['comment'] = 'water_surface_reference_datum_altitude is given in file {:s}'.format(options.dem_file)
    if not os.path.exists(options.dem_file):
        logger.error('path to dem file {:s} cannot be found'.format(options.dem_file))
        sys.exit(1)
    if not os.path.exists(options.ldd_file):
        logger.error('path to ldd file {:s} cannot be found'.format(options.ldd_file))
        sys.exit(1)

    # Read extent from a GDAL compatible file
    try:
        extent = inun_lib.get_gdal_extent(options.dem_file)
    except:
        msg = 'Input file {:s} not a gdal compatible file'.format(options.dem_file)
        inun_lib.close_with_error(logger, ch, msg)
        sys.exit(1)

    try:
        x, y = inun_lib.get_gdal_axes(options.dem_file, logging=logger)
        srs = inun_lib.get_gdal_projection(options.dem_file, logging=logger)
    except:
        msg = 'Input file {:s} not a gdal compatible file'.format(options.dem_file)
        inun_lib.close_with_error(logger, ch, msg)
        sys.exit(1)

    # read history from flood file
    if options.file_format == 0:
        a = nc.Dataset(options.flood_map, 'r')
        metadata_global['history'] = 'Created by: $Id: $, boundary conditions from {:s},\nhistory: {:s}'.format(os.path.abspath(options.flood_map), a.history)
        a.close()
    else:
        metadata_global['history'] = 'Created by: $Id: $, boundary conditions from {:s},\nhistory: {:s}'.format(os.path.abspath(options.flood_map), 'PCRaster file, no history')

    # first write subcatch maps and hand maps
    ############### TODO ######
    # setup a HAND file
    dem_name = os.path.split(options.dem_file)[1].split('.')[0]
    hand_file = os.path.join(options.dest_path, '{:s}_hand_strahler_{:02d}.tif'.format(dem_name, options.hand_strahler))
    if not(os.path.isfile(hand_file)):
    # hand file does not exist yet! Generate it, otherwise skip!
        logger.info('HAND file {:s} setting up...please wait...'.format(hand_file))
        hand_file_tmp = os.path.join(options.dest_path, '{:s}_hand_strahler_{:02d}.tif.tmp'.format(dem_name, options.hand_strahler))
        ds_hand = inun_lib.prepare_gdal(hand_file_tmp, x, y, logging=logger, srs=srs)
        band_hand = ds_hand.GetRasterBand(1)

        # Open terrain data for reading
        ds_dem, rasterband_dem = inun_lib.get_gdal_rasterband(options.dem_file)
        ds_ldd, rasterband_ldd = inun_lib.get_gdal_rasterband(options.ldd_file)
        ds_stream, rasterband_stream = inun_lib.get_gdal_rasterband(options.stream_file)
        n = 0
        for x_loop in range(0, len(x), options.x_tile):
            x_start = np.maximum(x_loop, 0)
            x_end = np.minimum(x_loop + options.x_tile, len(x))
            # determine actual overlap for cutting
            for y_loop in range(0, len(y), options.y_tile):
                x_overlap_min = x_start - np.maximum(x_start - options.x_overlap, 0)
                x_overlap_max = np.minimum(x_end + options.x_overlap, len(x)) - x_end
                n += 1
                # print('tile {:001d}:'.format(n))
                y_start = np.maximum(y_loop, 0)
                y_end = np.minimum(y_loop + options.y_tile, len(y))
                y_overlap_min = y_start - np.maximum(y_start - options.y_overlap, 0)
                y_overlap_max = np.minimum(y_end + options.y_overlap, len(y)) - y_end
                # cut out DEM
                logger.debug('Computing HAND for xmin: {:d} xmax: {:d} ymin {:d} ymax {:d}'.format(x_start, x_end,y_start, y_end))
                terrain = rasterband_dem.ReadAsArray(x_start - x_overlap_min,
                                                     y_start - y_overlap_min,
                                                     (x_end + x_overlap_max) - (x_start - x_overlap_min),
                                                     (y_end + y_overlap_max) - (y_start - y_overlap_min)
                                                     )

                drainage = rasterband_ldd.ReadAsArray(x_start - x_overlap_min,
                                                     y_start - y_overlap_min,
                                                     (x_end + x_overlap_max) - (x_start - x_overlap_min),
                                                     (y_end + y_overlap_max) - (y_start - y_overlap_min)
                                                     )
                stream = rasterband_stream.ReadAsArray(x_start - x_overlap_min,
                                                       y_start - y_overlap_min,
                                                       (x_end + x_overlap_max) - (x_start - x_overlap_min),
                                                       (y_end + y_overlap_max) - (y_start - y_overlap_min)
                                                       )
                # write to temporary file
                terrain_temp_file = os.path.join(options.dest_path, 'terrain_temp.map')
                drainage_temp_file = os.path.join(options.dest_path, 'drainage_temp.map')
                stream_temp_file = os.path.join(options.dest_path, 'stream_temp.map')
                if rasterband_dem.GetNoDataValue() is not None:
                    inun_lib.gdal_writemap(terrain_temp_file, 'PCRaster',
                                      np.arange(0, terrain.shape[1]),
                                      np.arange(0, terrain.shape[0]),
                                      terrain, rasterband_dem.GetNoDataValue(),
                                      gdal_type=gdal.GDT_Float32,
                                      logging=logger)
                else:
                    # in case no nodata value is found
                    logger.warning('No nodata value found in {:s}. assuming -9999'.format(options.dem_file))
                    inun_lib.gdal_writemap(terrain_temp_file, 'PCRaster',
                                      np.arange(0, terrain.shape[1]),
                                      np.arange(0, terrain.shape[0]),
                                      terrain, -9999.,
                                      gdal_type=gdal.GDT_Float32,
                                      logging=logger)

                inun_lib.gdal_writemap(drainage_temp_file, 'PCRaster',
                                  np.arange(0, terrain.shape[1]),
                                  np.arange(0, terrain.shape[0]),
                                  drainage, rasterband_ldd.GetNoDataValue(),
                                  gdal_type=gdal.GDT_Int32,
                                  logging=logger)
                inun_lib.gdal_writemap(stream_temp_file, 'PCRaster',
                                  np.arange(0, terrain.shape[1]),
                                  np.arange(0, terrain.shape[0]),
                                  stream, rasterband_ldd.GetNoDataValue(),
                                  gdal_type=gdal.GDT_Int32,
                                  logging=logger)
                # read as pcr objects
                pcr.setclone(terrain_temp_file)
                terrain_pcr = pcr.readmap(terrain_temp_file)
                drainage_pcr = pcr.lddrepair(pcr.ldd(pcr.readmap(drainage_temp_file)))  # convert to ldd type map
                stream_pcr = pcr.scalar(pcr.readmap(stream_temp_file))  # convert to ldd type map

                # compute streams
                stream_ge, subcatch = inun_lib.subcatch_stream(drainage_pcr, stream_pcr, options.hand_strahler) # generate streams

                basin = pcr.boolean(subcatch)
                hand_pcr, dist_pcr = inun_lib.derive_HAND(terrain_pcr, drainage_pcr, 3000, rivers=pcr.boolean(stream_ge), basin=basin)
                # convert to numpy
                hand = pcr.pcr2numpy(hand_pcr, -9999.)
                # cut relevant part
                if y_overlap_max == 0:
                    y_overlap_max = -hand.shape[0]
                if x_overlap_max == 0:
                    x_overlap_max = -hand.shape[1]
                hand_cut = hand[0+y_overlap_min:-y_overlap_max, 0+x_overlap_min:-x_overlap_max]

                band_hand.WriteArray(hand_cut, x_start, y_start)
                os.unlink(terrain_temp_file)
                os.unlink(drainage_temp_file)
                band_hand.FlushCache()
        ds_dem = None
        ds_ldd = None
        ds_stream = None
        band_hand.SetNoDataValue(-9999.)
        ds_hand = None
        logger.info('Finalizing {:s}'.format(hand_file))
        # rename temporary file to final hand file
        os.rename(hand_file_tmp, hand_file)
    else:
        logger.info('HAND file {:s} already exists...skipping...'.format(hand_file))

    #####################################################################################
    #  HAND file has now been prepared, moving to flood mapping part                    #
    #####################################################################################
    # load the staticmaps needed to estimate volumes across all
    xax, yax, riv_length, fill_value = inun_lib.gdal_readmap(options.riv_length_file, 'GTiff')
    riv_length = np.ma.masked_where(riv_length==fill_value, riv_length)
    xax, yax, riv_width, fill_value = inun_lib.gdal_readmap(options.riv_width_file, 'GTiff')
    riv_width[riv_width == fill_value] = 0

    x_res = np.abs((xax[-1]-xax[0])/(len(xax)-1))
    y_res = np.abs((yax[-1]-yax[0])/(len(yax)-1))

    flood_folder = os.path.join(options.dest_path, case_name)
    flood_vol_map = os.path.join(flood_folder, '{:s}_vol.tif'.format(os.path.split(options.flood_map)[1].split('.')[0]))
    if not(os.path.isdir(flood_folder)):
        os.makedirs(flood_folder)
    inun_file_tmp = os.path.join(flood_folder, '{:s}.tif.tmp'.format(case_name))
    inun_file = os.path.join(flood_folder, '{:s}.tif'.format(case_name))
    hand_temp_file = os.path.join(flood_folder, 'hand_temp.map')
    drainage_temp_file = os.path.join(flood_folder, 'drainage_temp.map')
    stream_temp_file = os.path.join(flood_folder, 'stream_temp.map')
    flood_vol_temp_file = os.path.join(flood_folder, 'flood_warp_temp.tif')
    # load the data with river levels and compute the volumes
    if options.file_format == 0:
        # assume we need the maximum value in a NetCDF time series grid
        a = nc.Dataset(options.flood_map, 'r')
        xax = a.variables['x'][:]
        yax = a.variables['y'][:]

        flood_series = a.variables[options.flood_variable][:]
        flood = flood_series.max(axis=0)
        if yax[-1] > yax[0]:
            yax = np.flipud(yax)
            flood = np.flipud(flood)
        a.close()
    elif options.file_format == 1:
        xax, yax, flood, flood_fill_value = inun_lib.gdal_readmap(options.flood_map, 'PCRaster')
    #res_x = x[1]-x[0]
    #res_y = y[1]-y[0]

    # load the bankfull depths
    if options.bankfull_map == '':
        bankfull = np.zeros(flood.shape)
    else:
        if options.file_format == 0:
            a = nc.Dataset(options.bankfull_map, 'r')
            xax = a.variables['x'][:]
            yax = a.variables['y'][:]
            bankfull = a.variables[options.flood_variable][0, :, :]
            if yax[-1] > yax[0]:
                yax = np.flipud(yax)
                bankfull = np.flipud(bankful)
            a.close()
        elif options.file_format == 1:
            xax, yax, bankfull, bankfull_fill_value = inun_lib.gdal_readmap(options.bankfull_map, 'PCRaster')
#     flood = bankfull*2
    # res_x = 2000
    # res_y = 2000
    # subtract the bankfull water level to get flood levels (above bankfull)
    flood_vol = np.maximum(flood-bankfull, 0)
    flood_vol_m = riv_length*riv_width*flood_vol/(x_res * y_res)  # volume expressed in meters water disc (1e6 is the surface area of one wflow grid cell)
    flood_vol_m_data = flood_vol_m.data
    flood_vol_m_data[flood_vol_m.mask] = -999.
    print('Saving water layer map to {:s}'.format(flood_vol_map))
    # write to a tiff file
    inun_lib.gdal_writemap(flood_vol_map, 'GTiff', xax, yax, np.maximum(flood_vol_m_data, 0), -999.)
    ds_hand, rasterband_hand = inun_lib.get_gdal_rasterband(hand_file)
    ds_ldd, rasterband_ldd = inun_lib.get_gdal_rasterband(options.ldd_file)
    ds_stream, rasterband_stream = inun_lib.get_gdal_rasterband(options.stream_file)

    logger.info('Preparing flood map in {:s} ...please wait...'.format(inun_file))
    ds_inun = inun_lib.prepare_gdal(inun_file_tmp, x, y, logging=logger, srs=srs)
    band_inun = ds_inun.GetRasterBand(1)

    # loop over all the tiles
    n = 0
    for x_loop in range(0, len(x), options.x_tile):
        x_start = np.maximum(x_loop, 0)
        x_end = np.minimum(x_loop + options.x_tile, len(x))
        # determine actual overlap for cutting
        for y_loop in range(0, len(y), options.y_tile):
            x_overlap_min = x_start - np.maximum(x_start - options.x_overlap, 0)
            x_overlap_max = np.minimum(x_end + options.x_overlap, len(x)) - x_end
            n += 1
            # print('tile {:001d}:'.format(n))
            y_start = np.maximum(y_loop, 0)
            y_end = np.minimum(y_loop + options.y_tile, len(y))
            y_overlap_min = y_start - np.maximum(y_start - options.y_overlap, 0)
            y_overlap_max = np.minimum(y_end + options.y_overlap, len(y)) - y_end
            x_tile_ax = x[x_start - x_overlap_min:x_end + x_overlap_max]
            y_tile_ax = y[y_start - y_overlap_min:y_end + y_overlap_max]

            # cut out DEM
            logger.debug('handling xmin: {:d} xmax: {:d} ymin {:d} ymax {:d}'.format(x_start, x_end, y_start, y_end))
            hand = rasterband_hand.ReadAsArray(x_start - x_overlap_min,
                                                 y_start - y_overlap_min,
                                                 (x_end + x_overlap_max) - (x_start - x_overlap_min),
                                                 (y_end + y_overlap_max) - (y_start - y_overlap_min)
                                                 )

            drainage = rasterband_ldd.ReadAsArray(x_start - x_overlap_min,
                                                 y_start - y_overlap_min,
                                                 (x_end + x_overlap_max) - (x_start - x_overlap_min),
                                                 (y_end + y_overlap_max) - (y_start - y_overlap_min)
                                                 )
            stream = rasterband_stream.ReadAsArray(x_start - x_overlap_min,
                                                   y_start - y_overlap_min,
                                                   (x_end + x_overlap_max) - (x_start - x_overlap_min),
                                                   (y_end + y_overlap_max) - (y_start - y_overlap_min)
                                                   )
            print('len x-ax: {:d} len y-ax {:d} x-shape {:d} y-shape {:d}'.format(len(x_tile_ax), len(y_tile_ax), hand.shape[1], hand.shape[0]))
            inun_lib.gdal_writemap(hand_temp_file, 'PCRaster',
                              x_tile_ax,
                              y_tile_ax,
                              hand, rasterband_hand.GetNoDataValue(),
                              gdal_type=gdal.GDT_Float32,
                              logging=logger)
            inun_lib.gdal_writemap(drainage_temp_file, 'PCRaster',
                              x_tile_ax,
                              y_tile_ax,
                              drainage, rasterband_ldd.GetNoDataValue(),
                              gdal_type=gdal.GDT_Int32,
                              logging=logger)
            inun_lib.gdal_writemap(stream_temp_file, 'PCRaster',
                              x_tile_ax,
                              y_tile_ax,
                              stream, rasterband_stream.GetNoDataValue(),
                              gdal_type=gdal.GDT_Int32,
                              logging=logger)
            # read as pcr objects
            pcr.setclone(hand_temp_file)
            hand_pcr = pcr.readmap(hand_temp_file)
            drainage_pcr = pcr.lddrepair(pcr.ldd(pcr.readmap(drainage_temp_file)))  # convert to ldd type map
            stream_pcr = pcr.scalar(pcr.readmap(drainage_temp_file))  # convert to ldd type map
            # prepare a subcatchment map

            stream_ge, subcatch = inun_lib.subcatch_stream(drainage_pcr, stream_pcr, options.catchment_strahler) # generate subcatchments
            drainage_surf = pcr.ifthen(stream_ge > 0, pcr.accuflux(drainage_pcr, 1))  # proxy of drainage surface inaccurate at tile edges
           # compute weights for spreadzone (1/drainage_surf)
            subcatch = pcr.spreadzone(subcatch, 0, 0)

            # TODO check weighting scheme, perhaps not necessary
            # weight = 1./pcr.scalar(pcr.spreadzone(pcr.cover(pcr.ordinal(drainage_surf), 0), 0, 0))
            # subcatch_fill = pcr.scalar(pcr.spreadzone(subcatch, 0, weight))
            # # cover subcatch with subcatch_fill
            # pcr.report(weight, 'weight_{:02d}.map'.format(n))
            # pcr.report(subcatch, 'subcatch_{:02d}.map'.format(n))
            # pcr.report(pcr.nominal(subcatch_fill), 'subcatch_fill_{:02d}.map'.format(n))
            inun_lib.gdal_warp(flood_vol_map, hand_temp_file, flood_vol_temp_file, gdal_interp=gdalconst.GRA_NearestNeighbour) # ,
            x_tile_ax, y_tile_ax, flood_meter, fill_value = inun_lib.gdal_readmap(flood_vol_temp_file, 'GTiff')
            # convert meter depth to volume [m3]
            flood_vol = pcr.numpy2pcr(pcr.Scalar, flood_meter, fill_value)*((x_tile_ax[1] - x_tile_ax[0]) * (y_tile_ax[0] - y_tile_ax[1]))  # resolution of SRTM *1166400000.
            ## now we have some nice volume. Now we need to redistribute!
            inundation_pcr = inun_lib.volume_spread(drainage_pcr, hand_pcr, subcatch, flood_vol,
                                           volume_thres=0., iterations=options.iterations,
                                           area_multiplier=options.area_multiplier) # 1166400000.
            inundation = pcr.pcr2numpy(inundation_pcr, -9999.)
            # cut relevant part
            if y_overlap_max == 0:
                y_overlap_max = -inundation.shape[0]
            if x_overlap_max == 0:
                x_overlap_max = -inundation.shape[1]
            inundation_cut = inundation[0+y_overlap_min:-y_overlap_max, 0+x_overlap_min:-x_overlap_max]
            # inundation_cut
            band_inun.WriteArray(inundation_cut, x_start, y_start)
            band_inun.FlushCache()
            # clean up
            os.unlink(flood_vol_temp_file)
            os.unlink(drainage_temp_file)
            os.unlink(hand_temp_file)

            # if n == 35:
            #     band_inun.SetNoDataValue(-9999.)
            #     ds_inun = None
            #     sys.exit(0)
    os.unlink(flood_vol_map)

    logger.info('Finalizing {:s}'.format(inun_file))
    # add the metadata to the file and band
    band_inun.SetNoDataValue(-9999.)
    ds_inun.SetMetadata(metadata_global)
    band_inun.SetMetadata(metadata_var)
    ds_inun = None
    ds_hand = None
    ds_ldd = None
    # rename temporary file to final hand file
    os.rename(inun_file_tmp, inun_file)

    logger.info('Done! Thank you for using hand_contour_inun.py')
    logger, ch = inun_lib.closeLogger(logger, ch)
    del logger, ch
    sys.exit(0)


if __name__ == "__main__":
    main()

