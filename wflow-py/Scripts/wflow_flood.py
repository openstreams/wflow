#!/usr/bin/python

"""
Definition of the wflow_flood post processor.
---------------------------------------------

Performs a planar volume spreading on outputs of a wflow\_sbm|hbv|routing model run.
The module can be used to post-process model outputs into a flood map that has a 
(much) higher resolution than the model resolution. 

The routine aggregates flooded water volumes occurring across all river pixels 
across a user-defined strahler order basin scale (typically a quite small subcatchment)
and then spreads this volume over a high resolution terrain model. To ensure that the 
flood volume is not spread in a rice-field kind of way (first filling the lowest cell in
the occurring subbasin), the terrain data is first normalised to a Height-Above-Nearest-Drain
(HAND) map of the associated typically flooding rivers (to be provided by user through a catchment order)
and flooding is estimated from this HAND map.

A sequential flood mapping is performed starting from the user defined Strahler order basin scale to the highest
Strahler orders. A HAND map is derived for each river order starting from the lowest order. These maps are then
used sequentially to spread flood volumes over the high resolution terrain model starting from the lowest catchment
order provided by the user (using the corresponding HAND map) to the highest stream order.
Backwater effects from flooding of higher orders catchments to lower order catchments are taken into account by taking
the maximum flood level of both.

Preferrably a user should use from the outputs of a wflow\_routing model
because then the user can use the floodplain water level only (usually saved 
in a variable name levfp). If estimates from a HBV or SBM (or other wflow) model are used
we recommend that the user also provides a "bank-full" water level in the command line 
arguments. If not provided, wflow_flood will also spread water volumes occuring
within the banks of the river, probably leading to an overestimation of flooding.

The wflow\_sbm|hbv model must have a saved mapstacks in NetCDF containing (over-bank) water levels
The module selects the maximum occurring value in the map stacks provided (NetCDF)
and spreads this out over a high resolution terrain dataset using the high resolution
terrain, associated ldd and stream order map.

TODO:: enable selection of a time step.

Ini-file settings
-----------------

The module uses an ini file and a number of command line arguments to run. The ini-file
contains inputs that are typically the same across a number of runs with the module for
a given study area (e.g. the used DEM, LDD, and some run parameters). For instance,
a user can prepare flood maps from different flood events computed with one WFLOW model,
using the same .ini file for each event.

The .ini file sections are treated below:

::

	[HighResMaps]
	dem_file = /p/1220657-afghanistan-disaster-risk/Processed DEMs/SRTM 90m merged/BEST90m_WGS_UTM42N.tif
	ldd_file = /p/1220657-afghanistan-disaster-risk/Processed DEMs/SRTM 90m merged/LDD/ldd_SRTM0090m_WGS_UTM42N.map
	stream_file = /p/1220657-afghanistan-disaster-risk/Processed DEMs/SRTM 90m merged/stream.map
	[wflowResMaps]
	riv_length_fact_file = /p/1220657-afghanistan-disaster-risk/floodhazardsimulations/Stepf_output/river_length_fact.map
	riv_width_file = /p/1220657-afghanistan-disaster-risk/floodhazardsimulations/Stepf_output/wflow_floodplainwidth.map
	ldd_wflow = /p/1220657-afghanistan-disaster-risk/floodhazardsimulations/Stepf_output/wflow_ldd.map
	[file_settings]
	latlon = 0
	file_format = 0

The dem\_file contains a file with the high-res terrain data. It MUST be in .tif format.
This is because .tif files can contain projection information. At the moment the .tif file
must have the same projection as the WFLOW model (can be WGS84, but also any local projection 
in meters), but we intend to also facilitate projections in the future.

The ldd\_file contains the ldd, derived from the dem_file (PCRaster format)

The stream\_file contains a stream order file (made with the PCRaster stream order file) 
derived from the LDD in ldd\_file.

riv\_length\_fact\_file and riv\_width\_file contain the dimensions of the channels within the WFLOW
pixels (unit meters) and are therefore in the resolution of the WFLOW model. The riv\_length\_fact\_file is used
 to derive a riv_length by multiplying the LDD length from cell to cell within the LDD network with the
wflow\_riverlength_fact.map map, typically located in the staticmaps folder of the used WFLOW model.
The width map is also in meters, and should contain the flood plain width in case the wflow_routing 
model is used (typical name is wflow_floodplainwidth.map). If a HBV or SBM model is used, you should 
use the river width map instead (typical name wflow_riverwidth.map).

ldd\_wflow is the ldd, derived at wflow resolution (typically wflow_ldd.map)

If latlon is 0, the cell-size is given in meters (the default)

If file_format is set to 0, the flood map is expected to be given in netCDF format (in the command line after -F)
if set to 1, format is expected to be PCRaster format

::

	[metadata_global]
	source=WFLOW model XXX
	institution=Deltares
	title=fluvial flood hazard from a wflow model
	references=http://www.deltares.nl/
	Conventions=CF-1.6
	project=Afhanistan multi-peril country wide risk assessment

In the metadata\_global section the user can enter typical project details as metadata. These are
used in the outcoming .tif file. We recommend to follow Climate and Forecast conventions for typical
metadata entries (see http://cfconventions.org/). You can insert as many key/value pairs as you like.
Some examples are given above.

::

	[tiling]
	x_tile=2000
	y_tile=2000
	x_overlap=500
	y_overlap=500

When very large domains are processed, the complete rasters will not fit into memory. In this
case, the routine will break the domain into several tiles and process these separately. The
x\_tile and y\_tile parameters are used to set the tile size. If you are confident that the whole
domain will fit into memory (typically when the size is smaller than about 5,000 x 5,000 rows and
columns) then just enter a number larger than the total amount of rows and columns. The x\_overlap
and y\_overlap parameters should be large enough to prevent edge effects at the edges of each tile
where averaging subbasins are cut off from the edge. Slightly larger tiles (defined by the overlap)
are therefore processed and the edges are consequently cut off after processing one tile to get a
seamless product.

Some trial and error may be required to yield the right tile sizes and overlaps.

::

	[inundation]
	iterations=20
	initial_level=32

The inundation section contains a number of settings for the flood fill algorithm.
The number of iterations can be changed, we recommend to set it to 20 for
an accurate results. The initial\_level is the largest water level that can occur during flooding. Make 
sure it is set to a level (much) higher than anticipated to occur but not to a value close to infinity.
If you set it orders too high, the solution will not converge to a reasonable estimate.

Command line arguments
----------------------

When wflow\_flood.py is run with the -h argument, you will receive the following feedback:

::

	python wflow_flood.py -h
	Usage: wflow_flood.py [options]

	Options:
	  -h, --help            show this help message and exit
	  -q, --quiet           do not print status messages to stdout
	  -i INIFILE, --ini=INIFILE
				ini configuration file
	  -f FLOOD_MAP, --flood_map=FLOOD_MAP
				Flood map file (NetCDF point time series file
	  -v FLOOD_VARIABLE, --flood_variable=FLOOD_VARIABLE
				variable name of flood water level
	  -b BANKFULL_MAP, --bankfull_map=BANKFULL_MAP
				Map containing bank full level (is subtracted from
				flood map, in NetCDF)
	  -c CATCHMENT_STRAHLER, --catchment=CATCHMENT_STRAHLER
				Smallest Strahler order threshold over which flooding
				may occur
	  -m MAX_CATCHMENT_STRAHLER, --max_catchment=MAX_CATCHMENT_STRAHLER
				Largest Strahler order over which flooding may occur
	  -d DEST_PATH, --destination=DEST_PATH
				Destination path
	  -H HAND_FILE_PREFIX, --hand_file=HAND_FILE_PREFIX
				optional HAND file prefix of already generated HAND maps for different Strahler orders


Further explanation:

    -i = the .ini file described in the previous section

    -f = The NetCDF output time series or a GeoTIFF containing the flood event to be downscaled. In case of NetCDF, this is a typical NetCDF output file from a WFLOW model. Alternatively, you can provide a GeoTIFF file that contains the exact flood depths for a given time step user defined statistic (for example the maximum value across all time steps)

    -v = Variable within the aforementioned file that contains the depth within the flood plain (typical levfp)
    
    -b = Similar file as -f but providing the bank full water level. Can e provided in case you know that a certain water depth is blocked, or remains within banks. In cae a NetCDF is provided, the maximum values are used, alternatively, you can provide a GeoTIFF.

    -c = starting point of catchment strahler order over which flood volumes are averaged, before spreading.
         The Height-Above-Nearest-Drain maps are derived from this Strahler order on (until max Strahler order).
         NB: This is the strahler order of the high resolution stream order map
    
    -m = maximum strahler order over which flooding may occur (default value is the highest order in high res stream map)
    
    -d = path where the file is stored
    
    -H = HAND file prefix. As interim product, the module produces HAND files. This is a very time consuming process and therefore the user can also supply previously generated HAND files here (GeoTIFF format)
        The names of the HAND files should be constructed as follows: hand_prefix_{:02d}.format{hand_strahler}, so for example hand_prefix_03 for HAND map with minimum Strahler order 3. (in this case -H hand_prefix should be given)
        Maps should be made for Strahler orders from -c to -m (or maximum strahler order in the stream map)
    
Outputs
-------

wflow\_flood produces the following outputs:

Table: outputs of the wflow_flood module

+----------------------------------------------------------------------------+---------------------------------------------------------------------+
|hand\_contour\_inun.log                                                     | log file of the module, contains info and error messages            |
+----------------------------------------------------------------------------+---------------------------------------------------------------------+
|inun_<-f>\_catch\_<-c>.tif                                                  | resulting inundation map (GeoTIFF)                                  |
+----------------------------------------------------------------------------+---------------------------------------------------------------------+
|<dem_file>\_hand\_strahler_<-c>.tif                                         | HAND file based upon strahler order given with -c (only without -H  |
+----------------------------------------------------------------------------+---------------------------------------------------------------------+

Questions can be directed to hessel.winsemius@deltares.nl

$Id: wflow_flood.py $
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

from osgeo import gdal, gdalconst
import pcraster as pcr
import netCDF4 as nc

import wflow_flood_lib as inun_lib
import wflow.pcrut as pcrut
import wflow.wflow_lib as wflow_lib
import datetime as dt

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
    parser.add_option('-t', '--time',
                      dest='time', default='',
                      help='time in YYYYMMDDHHMMSS, overrides time in NetCDF input if set')
    # parser.add_option('-s', '--hand_strahler',
    #                   dest='hand_strahler', default=7, type='int',
    #                   help='Strahler order threshold >= selected as riverine')
    parser.add_option('-m', '--max_strahler',
                      dest = 'max_strahler', default=1000, type='int',
                      help='Maximum Strahler order to loop over')
    parser.add_option('-d', '--destination',
                      dest='dest_path', default='inun',
                      help='Destination path')
    parser.add_option('-H', '--hand_file_prefix',
                      dest='hand_file_prefix', default='',
                      help='optional HAND file prefix of already generated HAND files')
    (options, args) = parser.parse_args()

    if not os.path.exists(options.inifile):
        print 'path to ini file cannot be found'
        sys.exit(1)
    options.dest_path = os.path.abspath(options.dest_path)

    if not(os.path.isdir(options.dest_path)):
        os.makedirs(options.dest_path)

    # set up the logger
    flood_name = os.path.split(options.flood_map)[1].split('.')[0]
    # case_name = 'inun_{:s}_hand_{:02d}_catch_{:02d}'.format(flood_name, options.hand_strahler, options.catchment_strahler)
    case_name = 'inun_{:s}_catch_{:02d}'.format(flood_name, options.catchment_strahler)
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
    options.dem_file = inun_lib.configget(config, 'HighResMaps',
                                  'dem_file',
                                  True)
    options.ldd_file = inun_lib.configget(config, 'HighResMaps',
                                'ldd_file',
                                 True)
    options.stream_file = inun_lib.configget(config, 'HighResMaps',
                                'stream_file',
                                 True)
    options.riv_length_fact_file = inun_lib.configget(config, 'wflowResMaps',
                                'riv_length_fact_file',
                                 True)
    options.ldd_wflow = inun_lib.configget(config, 'wflowResMaps',
                                'ldd_wflow',
                                True)
    options.riv_width_file = inun_lib.configget(config, 'wflowResMaps',
                                'riv_width_file',
                                 True)
    options.file_format = inun_lib.configget(config, 'file_settings',
                                'file_format', 0, datatype='int')
    options.out_format = inun_lib.configget(config, 'file_settings',
                                'out_format', 0, datatype='int')
    options.latlon = inun_lib.configget(config, 'file_settings',
                                 'latlon', 0, datatype='int')
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
    options.flood_volume_type = inun_lib.configget(config, 'inundation',
                                  'flood_volume_type', 0, datatype='int')

    # options.area_multiplier = inun_lib.configget(config, 'inundation',
    #                               'area_multiplier', 1., datatype='float')
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
    metadata_var['long_name'] = 'flooding'
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
    # setup a HAND file for each strahler order

    max_s = inun_lib.define_max_strahler(options.stream_file, logging=logger)
    stream_max = np.minimum(max_s, options.max_strahler)

    for hand_strahler in range(options.catchment_strahler, stream_max + 1, 1):
        dem_name = os.path.split(options.dem_file)[1].split('.')[0]
        if os.path.isfile('{:s}_{:02d}.tif'.format(options.hand_file_prefix, hand_strahler)):
            hand_file = '{:s}_{:02d}.tif'.format(options.hand_file_prefix, hand_strahler)
        else:
            logger.info('No HAND files with HAND prefix were found, checking {:s}_hand_strahler_{:02d}.tif'.format(dem_name, hand_strahler))
            hand_file = os.path.join(options.dest_path, '{:s}_hand_strahler_{:02d}.tif'.format(dem_name, hand_strahler))
        if not(os.path.isfile(hand_file)):
        # hand file does not exist yet! Generate it, otherwise skip!
            logger.info('HAND file {:s} not found, start setting up...please wait...'.format(hand_file))
            hand_file_tmp = os.path.join(options.dest_path, '{:s}_hand_strahler_{:02d}.tif.tmp'.format(dem_name, hand_strahler))
            ds_hand, band_hand = inun_lib.prepare_gdal(hand_file_tmp, x, y, logging=logger, srs=srs)
            # band_hand = ds_hand.GetRasterBand(1)

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

                    #check if the highest stream order of the tile is below the hand_strahler
                    # if the highest stream order of the tile is smaller than hand_strahler, than DEM values are taken instead of HAND values.
                    max_stream_tile = inun_lib.define_max_strahler(stream_temp_file, logging=logger)
                    if max_stream_tile < hand_strahler:
                        hand_pcr = terrain_pcr
                        logger.info('For this tile, DEM values are used instead of HAND because there is no stream order larger than {:02d}'.format(hand_strahler))
                    else:
                    # compute streams
                        stream_ge, subcatch = inun_lib.subcatch_stream(drainage_pcr, hand_strahler, stream=stream_pcr) # generate streams
                        # compute basins
                        stream_ge_dummy, subcatch = inun_lib.subcatch_stream(drainage_pcr, options.catchment_strahler, stream=stream_pcr) # generate streams
                        basin = pcr.boolean(subcatch)
                        hand_pcr, dist_pcr = inun_lib.derive_HAND(terrain_pcr, drainage_pcr, 3000,
                                                                  rivers=pcr.boolean(stream_ge), basin=basin)
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
                    os.unlink(stream_temp_file)
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
    # set the clone
    pcr.setclone(options.ldd_wflow)
    # read wflow ldd as pcraster object
    ldd_pcr = pcr.readmap(options.ldd_wflow)
    xax, yax, riv_width, fill_value = inun_lib.gdal_readmap(options.riv_width_file, 'GTiff', logging=logger)

    # determine cell length in meters using ldd_pcr as clone (if latlon=True, values are converted to m2
    x_res, y_res, reallength_wflow = pcrut.detRealCellLength(pcr.scalar(ldd_pcr), not(bool(options.latlon)))
    cell_surface_wflow = pcr.pcr2numpy(x_res * y_res, 0)

    if options.flood_volume_type == 0:
        # load the staticmaps needed to estimate volumes across all
        # xax, yax, riv_length, fill_value = inun_lib.gdal_readmap(options.riv_length_file, 'GTiff', logging=logger)
        # riv_length = np.ma.masked_where(riv_length==fill_value, riv_length)
        xax, yax, riv_width, fill_value = inun_lib.gdal_readmap(options.riv_width_file, 'GTiff', logging=logger)
        riv_width[riv_width == fill_value] = 0

        # read river length factor file (multiplier)
        xax, yax, riv_length_fact, fill_value = inun_lib.gdal_readmap(options.riv_length_fact_file, 'GTiff', logging=logger)
        riv_length_fact = np.ma.masked_where(riv_length_fact==fill_value, riv_length_fact)
        drain_length = wflow_lib.detdrainlength(ldd_pcr, x_res, y_res)

        # compute river length in each cell
        riv_length = pcr.pcr2numpy(drain_length, 0) * riv_length_fact
        # riv_length_pcr = pcr.numpy2pcr(pcr.Scalar, riv_length, 0)

    flood_folder = os.path.join(options.dest_path, case_name)
    flood_vol_map = os.path.join(flood_folder, '{:s}_vol.tif'.format(os.path.split(options.flood_map)[1].split('.')[0]))
    if not(os.path.isdir(flood_folder)):
        os.makedirs(flood_folder)
    if options.out_format == 0:
        inun_file_tmp = os.path.join(flood_folder, '{:s}.tif.tmp'.format(case_name))
        inun_file = os.path.join(flood_folder, '{:s}.tif'.format(case_name))
    else:
        inun_file_tmp = os.path.join(flood_folder, '{:s}.nc.tmp'.format(case_name))
        inun_file = os.path.join(flood_folder, '{:s}.nc'.format(case_name))

    hand_temp_file = os.path.join(flood_folder, 'hand_temp.map')
    drainage_temp_file = os.path.join(flood_folder, 'drainage_temp.map')
    stream_temp_file = os.path.join(flood_folder, 'stream_temp.map')
    flood_vol_temp_file = os.path.join(flood_folder, 'flood_warp_temp.tif')
    # load the data with river levels and compute the volumes
    if options.file_format == 0:
        # assume we need the maximum value in a NetCDF time series grid
        logger.info('Reading flood from {:s} NetCDF file'.format(options.flood_map))
        a = nc.Dataset(options.flood_map, 'r')
        if options.latlon == 0:
            xax = a.variables['x'][:]
            yax = a.variables['y'][:]
        else:
            xax = a.variables['lon'][:]
            yax = a.variables['lat'][:]
        if options.time == '':
            time_list = nc.num2date(a.variables['time'][:], units = a.variables['time'].units, calendar=a.variables['time'].calendar)
            time = [time_list[len(time_list)/2]]
        else:
            time = [dt.datetime.strptime(options.time, '%Y%m%d%H%M%S')]

        flood_series = a.variables[options.flood_variable][:]
        flood_data = flood_series.max(axis=0)
        if np.ma.is_masked(flood_data):
            flood = flood_data.data
            flood[flood_data.mask] = 0
        if yax[-1] > yax[0]:
            yax = np.flipud(yax)
            flood = np.flipud(flood)
        a.close()
    elif options.file_format == 1:
        logger.info('Reading flood from {:s} PCRaster file'.format(options.flood_map))
        xax, yax, flood, flood_fill_value = inun_lib.gdal_readmap(options.flood_map, 'PCRaster', logging=logger)
        flood = np.ma.masked_equal(flood, flood_fill_value)
        if options.time == '':
            options.time = '20000101000000'
        time = [dt.datetime.strptime(options.time, '%Y%m%d%H%M%S')]

        flood[flood==flood_fill_value] = 0.
    # load the bankfull depths
    if options.bankfull_map == '':
        bankfull = np.zeros(flood.shape)
    else:
        if options.file_format == 0:
            logger.info('Reading bankfull from {:s} NetCDF file'.format(options.bankfull_map))
            a = nc.Dataset(options.bankfull_map, 'r')
            xax = a.variables['x'][:]
            yax = a.variables['y'][:]

            bankfull_series = a.variables[options.flood_variable][:]
            bankfull_data = bankfull_series.max(axis=0)
            if np.ma.is_masked(bankfull_data):
                bankfull = bankfull_data.data
                bankfull[bankfull_data.mask] = 0
            if yax[-1] > yax[0]:
                yax = np.flipud(yax)
                bankfull = np.flipud(bankfull)
            a.close()
        elif options.file_format == 1:
            logger.info('Reading bankfull from {:s} PCRaster file'.format(options.bankfull_map))
            xax, yax, bankfull, bankfull_fill_value = inun_lib.gdal_readmap(options.bankfull_map, 'PCRaster', logging=logger)
            bankfull = np.ma.masked_equal(bankfull, bankfull_fill_value)
#     flood = bankfull*2
    # res_x = 2000
    # res_y = 2000
    # subtract the bankfull water level to get flood levels (above bankfull)
    flood_vol = np.maximum(flood-bankfull, 0)
    if options.flood_volume_type == 0:
        flood_vol_m = riv_length*riv_width*flood_vol/cell_surface_wflow  # volume expressed in meters water disc
        flood_vol_m_pcr = pcr.numpy2pcr(pcr.Scalar, flood_vol_m, 0)
    else:
        flood_vol_m = flood_vol/cell_surface_wflow
    flood_vol_m_data = flood_vol_m.data
    flood_vol_m_data[flood_vol_m.mask] = -999.
    logger.info('Saving water layer map to {:s}'.format(flood_vol_map))
    # write to a tiff file
    inun_lib.gdal_writemap(flood_vol_map, 'GTiff', xax, yax, np.maximum(flood_vol_m_data, 0), -999., logging=logger)
    # this is placed later in the hand loop
    # ds_hand, rasterband_hand = inun_lib.get_gdal_rasterband(hand_file)
    ds_ldd, rasterband_ldd = inun_lib.get_gdal_rasterband(options.ldd_file)
    ds_stream, rasterband_stream = inun_lib.get_gdal_rasterband(options.stream_file)

    logger.info('Preparing flood map in {:s} ...please wait...'.format(inun_file))
    if options.out_format == 0:
        ds_inun, band_inun = inun_lib.prepare_gdal(inun_file_tmp, x, y, logging=logger, srs=srs)
        # band_inun = ds_inun.GetRasterBand(1)
    else:
        ds_inun, band_inun = inun_lib.prepare_nc(inun_file_tmp, time, x, np.flipud(y), metadata=metadata_global,
                                                 metadata_var=metadata_var, logging=logger)
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

            # stream_max = np.minimum(stream.max(), options.max_strahler)


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
            pcr.setclone(stream_temp_file)
            drainage_pcr = pcr.lddrepair(pcr.ldd(pcr.readmap(drainage_temp_file)))  # convert to ldd type map
            stream_pcr = pcr.scalar(pcr.readmap(stream_temp_file))  # convert to ldd type map

            # warp of flood volume to inundation resolution
            inun_lib.gdal_warp(flood_vol_map, stream_temp_file, flood_vol_temp_file, gdal_interp=gdalconst.GRA_NearestNeighbour) # ,
            x_tile_ax, y_tile_ax, flood_meter, fill_value = inun_lib.gdal_readmap(flood_vol_temp_file, 'GTiff', logging=logger)
            # make sure that the option unittrue is on !! (if unitcell was is used in another function)
            x_res_tile, y_res_tile, reallength = pcrut.detRealCellLength(pcr.scalar(stream_pcr), not(bool(options.latlon)))
            cell_surface_tile = pcr.pcr2numpy(x_res_tile * y_res_tile, 0)

            # convert meter depth to volume [m3]
            flood_vol = pcr.numpy2pcr(pcr.Scalar, flood_meter*cell_surface_tile, fill_value)

            # first prepare a basin map, belonging to the lowest order we are looking at
            inundation_pcr = pcr.scalar(stream_pcr) * 0
            for hand_strahler in range(options.catchment_strahler, stream_max + 1, 1):
                # hand_temp_file = os.path.join(flood_folder, 'hand_temp.map')
                if os.path.isfile(os.path.join(options.dest_path, '{:s}_hand_strahler_{:02d}.tif'.format(dem_name, hand_strahler))):
                    hand_file = os.path.join(options.dest_path, '{:s}_hand_strahler_{:02d}.tif'.format(dem_name, hand_strahler))
                else:
                    hand_file = '{:s}_{:02d}.tif'.format(options.hand_file_prefix, hand_strahler)
                ds_hand, rasterband_hand = inun_lib.get_gdal_rasterband(hand_file)
                hand = rasterband_hand.ReadAsArray(x_start - x_overlap_min,
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

                hand_pcr = pcr.readmap(hand_temp_file)

                stream_ge_hand, subcatch_hand = inun_lib.subcatch_stream(drainage_pcr, options.catchment_strahler, stream=stream_pcr)
                # stream_ge_hand, subcatch_hand = inun_lib.subcatch_stream(drainage_pcr, hand_strahler, stream=stream_pcr)
                stream_ge, subcatch = inun_lib.subcatch_stream(drainage_pcr,
                                                               options.catchment_strahler,
                                                               stream=stream_pcr,
                                                               basin=pcr.boolean(pcr.cover(subcatch_hand, 0)),
                                                               assign_existing=True,
                                                               min_strahler=hand_strahler,
                                                               max_strahler=hand_strahler) # generate subcatchments, only within basin for HAND
                flood_vol_strahler = pcr.ifthenelse(pcr.boolean(pcr.cover(subcatch, 0)), flood_vol, 0) # mask the flood volume map with the created subcatch map for strahler order = hand_strahler

                inundation_pcr_step = inun_lib.volume_spread(drainage_pcr, hand_pcr,
                                                             pcr.subcatchment(drainage_pcr, subcatch), # to make sure backwater effects can occur from higher order rivers to lower order rivers
                                                             flood_vol_strahler,
                                                             volume_thres=0.,
                                                             iterations=options.iterations,
                                                             cell_surface=pcr.numpy2pcr(pcr.Scalar, cell_surface_tile, -9999),
                                                             logging=logger,
                                                             order=hand_strahler) # 1166400000.
                # use maximum value of inundation_pcr_step and new inundation for higher strahler order
                inundation_pcr = pcr.max(inundation_pcr, inundation_pcr_step)
            inundation = pcr.pcr2numpy(inundation_pcr, -9999.)
            # cut relevant part
            if y_overlap_max == 0:
                y_overlap_max = -inundation.shape[0]
            if x_overlap_max == 0:
                x_overlap_max = -inundation.shape[1]
            inundation_cut = inundation[0+y_overlap_min:-y_overlap_max, 0+x_overlap_min:-x_overlap_max]
            # inundation_cut
            if options.out_format == 0:
                band_inun.WriteArray(inundation_cut, x_start, y_start)
                band_inun.FlushCache()
            else:
                # with netCDF, data is up-side-down.
                inun_lib.write_tile_nc(band_inun, inundation_cut, x_start, y_start)
            # clean up
            os.unlink(flood_vol_temp_file)
            os.unlink(drainage_temp_file)
            os.unlink(hand_temp_file)
            os.unlink(stream_temp_file)     #also remove temp stream file from output folder

            # if n == 35:
            #     band_inun.SetNoDataValue(-9999.)
            #     ds_inun = None
            #     sys.exit(0)
    # os.unlink(flood_vol_map)

    logger.info('Finalizing {:s}'.format(inun_file))
    # add the metadata to the file and band
    # band_inun.SetNoDataValue(-9999.)
    # ds_inun.SetMetadata(metadata_global)
    # band_inun.SetMetadata(metadata_var)
    if options.out_format == 0:
        ds_inun = None
        ds_hand = None
    else:
        ds_inun.close()

    ds_ldd = None
    # rename temporary file to final hand file
    if os.path.isfile(inun_file):
        # remove an old result if available
        os.unlink(inun_file)
    os.rename(inun_file_tmp, inun_file)

    logger.info('Done! Thank you for using hand_contour_inun.py')
    logger, ch = inun_lib.closeLogger(logger, ch)
    del logger, ch
    sys.exit(0)


if __name__ == "__main__":
    main()

