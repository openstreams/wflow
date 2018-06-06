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
import sys
import os
import shutil
import glob
import numpy as np

# import admin packages
import configparser
from optparse import OptionParser

# import general packages
import numpy as np
from osgeo import osr, gdal, ogr, gdalconst

# import specific packages
# if frozen to exe, import wflow is needed before hydrotools to run __init__
# and set the PROJ_DIR environment variable to the right path
import wflow
import pcraster as pcr
import wflow.wflowtools_lib as wt
import wflow.wflow_lib as tr

import fiona
from rasterio import features
import rasterio

driver = ogr.GetDriverByName("ESRI Shapefile")


def clip_catchment_by_cell(cell_geom, catchment_geom):
    return cell_geom.intersection(catchment_geom)


def parse_args():
    ### Read input arguments #####
    parser = OptionParser()
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option(
        "-q",
        "--quiet",
        dest="verbose",
        default=True,
        action="store_false",
        help="do not print status messages to stdout",
    )
    parser.add_option(
        "-i",
        "--ini",
        dest="inifile",
        default=None,
        help="ini file with settings for static_maps.exe",
    )
    parser.add_option(
        "-s",
        "--source",
        dest="source",
        default="wflow",
        help="Source folder containing clone (default=./wflow)",
    )
    parser.add_option(
        "-d",
        "--destination",
        dest="destination",
        default="staticmaps",
        help="Destination folder (default=./staticmaps)",
    )
    parser.add_option(
        "-r",
        "--river",
        dest="rivshp",
        default=None,
        help="river network polyline layer (ESRI Shapefile)",
    )
    parser.add_option(
        "-c",
        "--catchment",
        dest="catchshp",
        default=None,
        help="catchment polygon layer (ESRI Shapefile)",
    )
    parser.add_option(
        "-g",
        "--gauges",
        dest="gaugeshp",
        default=None,
        help="gauge point layer (ESRI Shapefile)",
    )
    parser.add_option(
        "-D",
        "--dem",
        dest="dem_in",
        default=None,
        help="digital elevation model (GeoTiff)",
    )
    parser.add_option(
        "-L",
        "--landuse",
        dest="landuse",
        default=None,
        help="land use / land cover layer (GeoTiff)",
    )
    parser.add_option(
        "-S", "--soiltype", dest="soil", default=None, help="soil type layer (GeoTiff)"
    )
    parser.add_option(
        "-V",
        "--vegetation",
        dest="lai",
        default=None,
        help="vegetation LAI layer location (containing 12 GeoTiffs <LAI00000.XXX.tif>)",
    )
    parser.add_option(
        "-O",
        "--other_maps",
        dest="other_maps",
        default=None,
        help="bracketed [] comma-separated list of paths to other maps that should be reprojected",
    )
    parser.add_option(
        "-C",
        "--clean",
        dest="clean",
        default=False,
        action="store_true",
        help="Clean the .xml files from static maps folder when finished",
    )
    parser.add_option(
        "-A",
        "--alltouch",
        dest="alltouch",
        default=False,
        action="store_true",
        help='option to burn catchments "all touching".\nUseful when catchment-size is small compared to cellsize',
    )
    parser.add_option(
        "-l",
        "--logfile",
        dest="logfilename",
        default="wtools_static_maps.log",
        help="log file name",
    )

    (options, args) = parser.parse_args()
    return options


def main(
    source,
    destination,
    inifile,
    dem_in,
    rivshp,
    catchshp,
    gaugeshp=None,
    landuse=None,
    soil=None,
    lai=None,
    other_maps=None,
    logfilename="wtools_static_maps.log",
    verbose=True,
    clean=True,
    alltouch=False,
    outlets=([], []),
):
    # parse other maps into an array
    if not other_maps == None:
        if type(other_maps) == str:
            print(other_maps)
            other_maps = other_maps.replace(
                ' ', '').replace('[', '').replace(']', '').split(',')



    source = os.path.abspath(source)
    clone_map = os.path.join(source, "mask.map")
    clone_shp = os.path.join(source, "mask.shp")
    clone_prj = os.path.join(source, "mask.prj")

    if None in (rivshp, catchshp, dem_in):
        msg = """The following files are compulsory:
        - DEM (raster)
        - river (shape)
        - catchment (shape)
        """
        print (msg)
        parser.print_help()
        sys.exit(1)
    if (inifile is not None) and (not os.path.exists(inifile)):
        print('path to ini file cannot be found')
        sys.exit(1)
    if not os.path.exists(rivshp):
        print('path to river shape cannot be found')
        sys.exit(1)
    if not os.path.exists(catchshp):
        print('path to catchment shape cannot be found')
        sys.exit(1)
    if not os.path.exists(dem_in):
        print('path to DEM cannot be found')
        sys.exit(1)

    # open a logger, dependent on verbose print to screen or not
    logger, ch = wt.setlogger(logfilename, "WTOOLS", verbose)

    # create directories # TODO: check if workdir is still necessary, try to
    # keep in memory as much as possible

    # delete old files (when the source and destination folder are different)
    if np.logical_and(os.path.isdir(destination), destination is not source):
        shutil.rmtree(destination)
    if destination is not source:
        os.makedirs(destination)

    # Read mask
    if not (os.path.exists(clone_map)):
        logger.error(
            "Clone file {:s} not found. Please run create_grid first.".format(clone_map)
        )
        sys.exit(1)
    else:
        # set clone
        pcr.setclone(clone_map)
        # get the extent from clone.tif
        xax, yax, clone, fill_value = wt.gdal_readmap(clone_map, "GTiff")
        trans = wt.get_geotransform(clone_map)
        extent = wt.get_extent(clone_map)
        xmin, ymin, xmax, ymax = extent
        zeros = np.zeros(clone.shape)
        ones = pcr.numpy2pcr(pcr.Scalar, np.ones(clone.shape), -9999)
        # get the projection from clone.tif
        srs = wt.get_projection(clone_map)
        unit_clone = srs.GetAttrValue("UNIT").lower()

    # READ CONFIG FILE
    # open config-file
    if inifile is None:
        config = configparser.ConfigParser()
        config.optionxform = str
    else:
        config = wt.OpenConf(inifile)

    # read settings
    snapgaugestoriver = wt.configget(
        config, "settings", "snapgaugestoriver", True, datatype="boolean"
    )
    burnalltouching = wt.configget(
        config, "settings", "burncatchalltouching", True, datatype="boolean"
    )
    burninorder = wt.configget(
        config, "settings", "burncatchalltouching", False, datatype="boolean"
    )
    verticetollerance = wt.configget(
        config, "settings", "vertice_tollerance", 0.0001, datatype="float"
    )

    """ read parameters """
    burn_outlets = wt.configget(
        config, "parameters", "burn_outlets", 10000, datatype="int"
    )
    burn_rivers = wt.configget(config, "parameters", "burn_rivers", 200, datatype="int")
    burn_connections = wt.configget(
        config, "parameters", "burn_connections", 100, datatype="int"
    )
    burn_gauges = wt.configget(config, "parameters", "burn_gauges", 100, datatype="int")
    minorder = wt.configget(config, "parameters", "riverorder_min", 3, datatype="int")
    try:
        percentiles_str = wt.configget(config, 'parameters',
            'statisticmaps', '0, 100',
            datatype='str')
        percentiles_split = percentiles_str.replace(' ', '').split(',')
        percentiles = np.array(percentiles_split, dtype='float')
    except configparser.NoOptionError:
        percentiles = [0.0, 100.0]
    # read the parameters for generating a temporary very high resolution grid
    if unit_clone == "degree":
        cellsize_hr = wt.configget(
            config, "parameters", "highres_degree", 0.0005, datatype="float"
        )
    elif (unit_clone == "metre") or (unit_clone == "meter"):
        cellsize_hr = wt.configget(
            config, "parameters", "highres_metre", 50, datatype="float"
        )

    cols_hr = int((float(xmax) - float(xmin)) / cellsize_hr + 2)
    rows_hr = int((float(ymax) - float(ymin)) / cellsize_hr + 2)
    hr_trans = (float(xmin), cellsize_hr, float(0), float(ymax), 0, -cellsize_hr)
    clone_hr = os.path.join(destination, "clone_highres.tif")
    # make a highres clone as well!
    wt.CreateTif(clone_hr, rows_hr, cols_hr, hr_trans, srs, 0)

    # read staticmap locations
    catchment_map = wt.configget(
        config, "staticmaps", "catchment", "wflow_catchment.map"
    )
    dem_map = wt.configget(config, "staticmaps", "dem", "wflow_dem.map")
    demmax_map = wt.configget(config, "staticmaps", "demmax", "wflow_demmax.map")
    demmin_map = wt.configget(config, "staticmaps", "demmin", "wflow_demmin.map")
    gauges_map = wt.configget(config, "staticmaps", "gauges", "wflow_gauges.map")
    landuse_map = wt.configget(config, "staticmaps", "landuse", "wflow_landuse.map")
    ldd_map = wt.configget(config, "staticmaps", "ldd", "wflow_ldd.map")
    river_map = wt.configget(config, "staticmaps", "river", "wflow_river.map")
    outlet_map = wt.configget(config, "staticmaps", "outlet", "wflow_outlet.map")
    riverlength_fact_map = wt.configget(
        config, "staticmaps", "riverlength_fact", "wflow_riverlength_fact.map"
    )
    soil_map = wt.configget(config, "staticmaps", "soil", "wflow_soil.map")
    streamorder_map = wt.configget(
        config, "staticmaps", "streamorder", "wflow_streamorder.map"
    )
    subcatch_map = wt.configget(config, "staticmaps", "subcatch", "wflow_subcatch.map")

    # read mask location (optional)
    masklayer = wt.configget(config, "mask", "masklayer", catchshp)

    # ???? empty = pcr.ifthen(ones == 0, pcr.scalar(0))

    # TODO: check if extents are correct this way
    # TODO: check what the role of missing values is in zeros and ones (l. 123
    # in old code)

    # first add a missing value to dem_in
    ds = gdal.Open(dem_in, gdal.GA_Update)
    RasterBand = ds.GetRasterBand(1)
    fill_val = RasterBand.GetNoDataValue()

    if fill_val is None:
        RasterBand.SetNoDataValue(-9999)
    ds = None

    # reproject to clone map: see http://stackoverflow.com/questions/10454316/how-to-project-and-resample-a-grid-to-match-another-grid-with-gdal-python
    # resample DEM
    logger.info(
        "Resampling dem from {:s} to {:s}".format(
            os.path.abspath(dem_in), os.path.join(destination, dem_map)
        )
    )
    wt.gdal_warp(
        dem_in,
        clone_map,
        os.path.join(destination, dem_map),
        format="PCRaster",
        gdal_interp=gdalconst.GRA_Average,
    )
    # retrieve amount of rows and columns from clone
    # TODO: make windowstats applicable to source/target with different projections. This does not work yet.
    # retrieve srs from DEM
    try:
        srs_dem = wt.get_projection(dem_in)
    except:
        logger.warning("No projection found in DEM, assuming WGS 1984 lat long")
        srs_dem = osr.SpatialReference()
        srs_dem.ImportFromEPSG(4326)
    clone2dem_transform = osr.CoordinateTransformation(srs, srs_dem)
    # if srs.ExportToProj4() == srs_dem.ExportToProj4():

    wt.windowstats(
        dem_in,
        len(yax),
        len(xax),
        trans,
        srs,
        destination,
        percentiles,
        transform=clone2dem_transform,
        logger=logger,
    )

    ## read catchment shape-file to create catchment map
    src = rasterio.open(clone_map)
    shapefile = fiona.open(catchshp, "r")
    catchment_shapes = [feature["geometry"] for feature in shapefile]
    image = features.rasterize(
        catchment_shapes, out_shape=src.shape, all_touched=True, transform=src.transform
    )
    catchment_domain = pcr.numpy2pcr(pcr.Ordinal, image.copy(), 0)

    ## read river shape-file and create burn layer
    shapefile = fiona.open(rivshp, "r")
    river_shapes = [feature["geometry"] for feature in shapefile]
    image = features.rasterize(
        river_shapes, out_shape=src.shape, all_touched=False, transform=src.transform
    )
    rivers = pcr.numpy2pcr(pcr.Nominal, image.copy(), 0)
    riverdem = pcr.scalar(rivers) * pcr.readmap(os.path.join(destination, dem_map))
    pcr.setglobaloption("lddin")
    riverldd = pcr.lddcreate(riverdem, 1e35, 1e35, 1e35, 1e35)

    riveroutlet = pcr.cover(pcr.ifthen(pcr.scalar(riverldd) == 5, pcr.scalar(1000)), 0)
    burn_layer = pcr.cover(
        (
            pcr.scalar(
                pcr.ifthen(pcr.streamorder(riverldd) > 1, pcr.streamorder(riverldd))
            )
            - 1
        )
        * 1000
        + riveroutlet,
        0,
    )

    outlets_x, outlets_y = outlets
    n_outlets = len(outlets_x)
    logger.info("Number of outlets: {}".format(n_outlets))
    if n_outlets >= 1:
        outlets_map_numbered = tr.points_to_map(
            pcr.scalar(0), outlets_x, outlets_y, 0.5
        )
        outlets_map = pcr.boolean(outlets_map_numbered)
        # snap outlets to closest river (max 1 cell closer to river)
        outlets_map = pcr.boolean(
            pcr.cover(tr.snaptomap(pcr.ordinal(outlets_map), rivers), 0)
        )

    ## create ldd per catchment
    logger.info("Calculating ldd")
    ldddem = pcr.scalar(clone_map)

    # per subcatchment, burn dem, then create modified dem that fits the ldd of the subcatchment
    # this ldd dem is merged over catchments, to create a global ldd that abides to the subcatchment boundaries
    for idx, shape in enumerate(catchment_shapes):
        logger.info(
            "Computing ldd for catchment "
            + str(idx + 1)
            + "/"
            + str(len(catchment_shapes))
        )
        image = features.rasterize(
            [shape], out_shape=src.shape, all_touched=True, transform=src.transform
        )
        catchment = pcr.numpy2pcr(pcr.Scalar, image.copy(), 0)
        dem_burned_catchment = (
            pcr.readmap(os.path.join(destination, dem_map))
            * pcr.scalar(catchment_domain)
            * catchment
        ) - burn_layer
        # ldddem_catchment = pcr.lddcreatedem(
        #    dem_burned_catchment, 1e35, 1e35, 1e35, 1e35)
        ldddem = pcr.cover(ldddem, dem_burned_catchment)

    pcr.report(ldddem, os.path.join(destination, "ldddem.map"))

    wflow_ldd = pcr.lddcreate(ldddem, 1e35, 1e35, 1e35, 1e35)
    if n_outlets >= 1:
        # set outlets to pit
        wflow_ldd = pcr.ifthenelse(outlets_map, pcr.ldd(5), wflow_ldd)
        wflow_ldd = pcr.lddrepair(wflow_ldd)

    pcr.report(wflow_ldd, os.path.join(destination, "wflow_ldd.map"))

    # compute stream order, identify river cells
    streamorder = pcr.ordinal(pcr.streamorder(wflow_ldd))
    river = pcr.ifthen(streamorder >= pcr.ordinal(minorder), pcr.boolean(1))
    # find the minimum value in the DEM and cover missing values with a river with this value. Effect is none!! so now left out!
    # mindem = int(np.min(pcr.pcr2numpy(pcr.ordinal(os.path.join(destination, dem_map)),9999999)))
    # dem_resample_map = pcr.cover(os.path.join(destination, dem_map), pcr.scalar(river)*0+mindem)
    # pcr.report(dem_resample_map, os.path.join(destination, dem_map))
    pcr.report(streamorder, os.path.join(destination, streamorder_map))
    pcr.report(river, os.path.join(destination, river_map))

    # deal with your catchments
    if gaugeshp == None:
        logger.info("No gauges defined, using outlets instead")
        gauges = pcr.ordinal(
            pcr.uniqueid(
                pcr.boolean(pcr.ifthen(pcr.scalar(wflow_ldd) == 5, pcr.boolean(1)))
            )
        )
        pcr.report(gauges, os.path.join(destination, gauges_map))
    # TODO: Add the gauge shape code from StaticMaps.py (line 454-489)
    # TODO: add river length map (see SticMaps.py, line 492-499)

    # since the products here (river length fraction) are not yet used
    # this is disabled for now, as it also takes a lot of computation time
    if False:
        # report river length
        # make a high resolution empty map
        dem_hr_file = os.path.join(destination, "dem_highres.tif")
        burn_hr_file = os.path.join(destination, "burn_highres.tif")
        demburn_hr_file = os.path.join(destination, "demburn_highres.map")
        riv_hr_file = os.path.join(destination, "riv_highres.map")
        wt.gdal_warp(dem_in, clone_hr, dem_hr_file)
        # wt.CreateTif(riv_hr, rows_hr, cols_hr, hr_trans, srs, 0)
        # open the shape layer
        ds = ogr.Open(rivshp)
        lyr = ds.GetLayer(0)
        wt.ogr_burn(
            lyr,
            clone_hr,
            -100,
            file_out=burn_hr_file,
            format="GTiff",
            gdal_type=gdal.GDT_Float32,
            fill_value=0,
        )
        # read dem and burn values and add
        xax_hr, yax_hr, burn_hr, fill = wt.gdal_readmap(burn_hr_file, "GTiff")
        burn_hr[burn_hr == fill] = 0
        xax_hr, yax_hr, dem_hr, fill = wt.gdal_readmap(dem_hr_file, "GTiff")
        dem_hr[dem_hr == fill] = np.nan
        demburn_hr = dem_hr + burn_hr
        demburn_hr[np.isnan(demburn_hr)] = -9999
        wt.gdal_writemap(
            demburn_hr_file, "PCRaster", xax_hr, yax_hr, demburn_hr, -9999.
        )
        pcr.setclone(demburn_hr_file)
        demburn_hr = pcr.readmap(demburn_hr_file)

        logger.info("Calculating ldd to determine river length")
        ldd_hr = pcr.lddcreate(demburn_hr, 1e35, 1e35, 1e35, 1e35)
        pcr.report(ldd_hr, os.path.join(destination, "ldd_hr.map"))
        pcr.setglobaloption("unitcell")
        riv_hr = pcr.scalar(pcr.streamorder(ldd_hr) >= minorder) * pcr.downstreamdist(
            ldd_hr
        )
        pcr.report(riv_hr, riv_hr_file)
        pcr.setglobaloption("unittrue")
        pcr.setclone(clone_map)
        logger.info("Computing river length")
        wt.windowstats(
            riv_hr_file,
            len(yax),
            len(xax),
            trans,
            srs,
            destination,
            stat="fact",
            transform=False,
            logger=logger,
        )
        # TODO: nothing happens with the river lengths yet. Need to decide how to use these

    # report outlet map
    pcr.report(
        pcr.ifthen(pcr.ordinal(wflow_ldd) == 5, pcr.ordinal(1)),
        os.path.join(destination, outlet_map),
    )

    # report subcatchment map
    subcatchment = pcr.subcatchment(wflow_ldd, gauges)
    pcr.report(pcr.ordinal(subcatchment), os.path.join(destination, subcatch_map))

    # Report land use map
    if landuse == None:
        logger.info(
            "No land use map used. Preparing {:s} with only ones.".format(
                os.path.join(destination, landuse_map)
            )
        )
        pcr.report(pcr.nominal(ones), os.path.join(destination, landuse_map))
    else:
        logger.info(
            "Resampling land use from {:s} to {:s}".format(
                os.path.abspath(landuse),
                os.path.join(destination, os.path.abspath(landuse_map)),
            )
        )
        wt.gdal_warp(
            landuse,
            clone_map,
            os.path.join(destination, landuse_map),
            format="PCRaster",
            gdal_interp=gdalconst.GRA_Mode,
            gdal_type=gdalconst.GDT_Int32,
        )

    # report soil map
    if soil == None:
        logger.info(
            "No soil map used. Preparing {:s} with only ones.".format(
                os.path.join(destination, soil_map)
            )
        )
        pcr.report(pcr.nominal(ones), os.path.join(destination, soil_map))
    else:
        logger.info(
            "Resampling soil from {:s} to {:s}".format(
                os.path.abspath(soil),
                os.path.join(destination, os.path.abspath(soil_map)),
            )
        )
        wt.gdal_warp(
            soil,
            clone_map,
            os.path.join(destination, soil_map),
            format="PCRaster",
            gdal_interp=gdalconst.GRA_Mode,
            gdal_type=gdalconst.GDT_Int32,
        )

    if lai == None:
        logger.info(
            "No vegetation LAI maps used. Preparing default maps {:s} with only ones.".format(
                os.path.join(destination, soil_map)
            )
        )
        pcr.report(pcr.nominal(ones), os.path.join(destination, soil_map))
    else:
        dest_lai = os.path.join(destination, "clim")
        os.makedirs(dest_lai)
        for month in range(12):
            lai_in = os.path.join(lai, "LAI00000.{:03d}".format(month + 1))
            lai_out = os.path.join(dest_lai, "LAI00000.{:03d}".format(month + 1))
            logger.info(
                "Resampling vegetation LAI from {:s} to {:s}".format(
                    os.path.abspath(lai_in), os.path.abspath(lai_out)
                )
            )
            wt.gdal_warp(
                lai_in,
                clone_map,
                lai_out,
                format="PCRaster",
                gdal_interp=gdalconst.GRA_Bilinear,
                gdal_type=gdalconst.GDT_Float32,
            )

    # report soil map
    if other_maps == None:
        logger.info("No other maps used. Skipping other maps.")
    else:
        logger.info("Resampling list of other maps...")
        for map_file in other_maps:
            map_name = os.path.split(map_file)[1]
            logger.info(
                "Resampling a map from {:s} to {:s}".format(
                    os.path.abspath(map_file),
                    os.path.join(
                        destination,
                        os.path.splitext(os.path.basename(map_file))[0] + ".map",
                    ),
                )
            )
            wt.gdal_warp(
                map_file,
                clone_map,
                os.path.join(
                    destination,
                    os.path.splitext(os.path.basename(map_file))[0] + ".map",
                ),
                format="PCRaster",
                gdal_interp=gdalconst.GRA_Mode,
                gdal_type=gdalconst.GDT_Float32,
            )

    if clean:
        wt.DeleteList(glob.glob(os.path.join(destination, "*.xml")), logger=logger)
        wt.DeleteList(
            glob.glob(os.path.join(destination, "clim", "*.xml")), logger=logger
        )
        wt.DeleteList(glob.glob(os.path.join(destination, "*highres*")), logger=logger)


if __name__ == "__main__":
    argdict = parse_args()
    main(**vars(argdict))
