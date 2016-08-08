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
(HAND) map of the associated typically flooding rivers (to be provided by user through a strahler order)
and flooding is estimated from this HAND map. 

TODO:: perform routing from small to large scale using a sequential flood mapping from small
to large strahler orders.

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

	[maps]
	dem_file = /p/1220657-afghanistan-disaster-risk/Processed DEMs/SRTM 90m merged/BEST90m_WGS_UTM42N.tif
	ldd_file = /p/1220657-afghanistan-disaster-risk/Processed DEMs/SRTM 90m merged/LDD/ldd_SRTM0090m_WGS_UTM42N.map
	stream_file = /p/1220657-afghanistan-disaster-risk/Processed DEMs/SRTM 90m merged/stream.map
	riv_length_file = /p/1220657-afghanistan-disaster-risk/floodhazardsimulations/Stepf_output/river_length.map
	riv_width_file = /p/1220657-afghanistan-disaster-risk/floodhazardsimulations/Stepf_output/wflow_floodplainwidth.map

The dem\_file contains a file with the high-res terrain data. It MUST be in .tif format.
This is because .tif files can contain projection information. At the moment the .tif file
must have the same projection as the WFLOW model (can be WGS84, but also any local projection 
in meters), but we intend to also facilitate projections in the future.

The ldd\_file contains the ldd, derived from the dem_file (PCRaster format)

The stream\_file contains a stream order file (made with the PCRaster stream order file) 
derived from the LDD in ldd\_file.

riv\_length\_file and riv\_width\_file contain the dimensions of the channels within the WFLOW 
pixels (unit meters) and are therefore in the resolution of the WFLOW model. The user can derive
these by multiplying the LDD length from cell to cell within the LDD network with the 
wflow\_riverlength_fact.map map, typically located in the staticmaps folder of the used WFLOW model.
The width map is also in meters, and should contain the flood plain width in case the wflow_routing 
model is used (typical name is wflow_floodplainwidth.map). If a HBV or SBM model is used, you should 
use the river width map instead (typical name wflow_riverwidth.map).

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
	area_multiplier=1
	iterations=20
	initial_level=32

The inundation section contains a number of settings for the flood fill algorithm. The area_multiplier
should for the moment always be set to 1. This may change in the future of reprojection from one to 
another projection is considered. The number of iterations can be changed, we recommend to set it to 20 for
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
				Strahler order threshold >= are selected as catchment
				boundaries
	  -s HAND_STRAHLER, --hand_strahler=HAND_STRAHLER
				Strahler order threshold >= selected as riverine
	  -d DEST_PATH, --destination=DEST_PATH
				Destination path
	  -H HAND_FILE, --hand_file=HAND_FILE
				optional HAND file (already generated)


Further explanation:

    -i = the .ini file described in the previous section

    -f = The NetCDF output time series or a GeoTIFF containing the flood event to be downscaled. In case of NetCDF, this is a typical NetCDF output file from a WFLOW model. Alternatively, you can provide a GeoTIFF file that contains the exact flood depths for a given time step user defined statistic (for example the maximum value across all time steps)

    -v = Variable within the aforementioned file that contains the depth within the flood plain (typical levfp)
    
    -b = Similar file as -f but providing the bank full water level. Can e provided in case you know that a certain water depth is blocked, or remains within banks. In cae a NetCDF is provided, the maximum values are used, alternatively, you can provide a GeoTIFF.

    -c = catchment strahler order over which flood volumes are averaged, before spreading. This is the strahler order, at the resolution of the flood map (not WFLOW)
    
    -s = strahler order used to derive the Heigh-Above-Nearest-Drain, from which flood mapping originates
    
    -d = path where the file is stored
    
    -H = HAND file. As interim product, the module produces a HAND file. This is a very time consuming process and therefore the user can also supply a previously generated HAND file here (GeoTIFF format)
    
Outputs
-------

wflow\_flood produces the following outputs:

Table: outputs of the wflow_flood module

+----------------------------------------------------------------------------+---------------------------------------------------------------------+
|hand\_contour\_inun.log                                                     | log file of the module, contains info and error messages            |
+----------------------------------------------------------------------------+---------------------------------------------------------------------+
|inun_<-f>\_hand\_<-s>\_catch\_<-c>\\inun_<-f>\_hand\_<-s>\_catch\_<-c>.tif   | resulting inundation map (GeoTIFF)                                  |
+----------------------------------------------------------------------------+---------------------------------------------------------------------+
|<dem_file>\_hand\_strahler_<-s>.tif                                         | HAND file based upon strahler order given with -s (only without -H  |
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
    parser.add_option('-H', '--hand_file',
                      dest='hand_file', default='',
                      help='optional HAND file (already generated)')
    (options, args) = parser.parse_args()

    if not os.path.exists(options.inifile):
        print 'path to ini file cannot be found'
        sys.exit(1)
    options.dest_path = os.path.abspath(options.dest_path)

    if not(os.path.isdir(options.dest_path)):
        os.makedirs(options.dest_path)

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
    # setup a HAND file
    dem_name = os.path.split(options.dem_file)[1].split('.')[0]
    if os.path.isfile(options.hand_file):
        hand_file = options.hand_file
    else:
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
                pcr.report(stream_ge, 'stream_ge.map')
                pcr.report(pcr.scalar(subcatch), 'subcatch.map')
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
        logger.info('Reading flood from {:s} NetCDF file'.format(options.flood_map))
        a = nc.Dataset(options.flood_map, 'r')
        xax = a.variables['x'][:]
        yax = a.variables['y'][:]

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
        xax, yax, flood, flood_fill_value = inun_lib.gdal_readmap(options.flood_map, 'PCRaster')
        flood[flood==flood_fill_value] = 0.
    #res_x = x[1]-x[0]
    #res_y = y[1]-y[0]

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

            # TODO: loop over river order options.catchment_strahler to maximum order found in tile
            stream_ge, subcatch = inun_lib.subcatch_stream(drainage_pcr, stream_pcr, options.catchment_strahler) # generate subcatchments
            # TODO: mask subcatchments with river order lower than thres or higher than thres
            drainage_surf = pcr.ifthen(stream_ge > 0, pcr.accuflux(drainage_pcr, 1))  # proxy of drainage surface inaccurate at tile edges
           # compute weights for spreadzone (1/drainage_surf)
            subcatch = pcr.spreadzone(subcatch, 0, 0)


            # TODO: put areas outside subcatch to zero
            ## now we have some nice volume. Now we need to redistribute!
            # TODO: insert selected HAND map
            inundation_pcr = inun_lib.volume_spread(drainage_pcr, hand_pcr, subcatch, flood_vol,
                                           volume_thres=0., iterations=options.iterations,
                                           area_multiplier=options.area_multiplier) # 1166400000.
            inundation = pcr.pcr2numpy(inundation_pcr, -9999.)
            # cut relevant part
            if y_overlap_max == 0:
                y_overlap_max = -inundation.shape[0]
            if x_overlap_max == 0:
                x_overlap_max = -inundation.shape[1]
            # TODO: use maximum value of inundation_cut and new inundation for higher strahler order
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

