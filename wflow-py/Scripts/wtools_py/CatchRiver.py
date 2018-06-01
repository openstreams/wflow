# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 09:22:23 2015

@author: tollena
"""

from subprocess import call
from osgeo import osr
import os
import sys
import pcraster as pcr
import getopt
import shutil
import numpy as np
from osgeo import ogr
from osgeo import gdal
from osgeo.gdalconst import *

import wflow.wflowtools_lib as wt

Driver = ogr.GetDriverByName("ESRI Shapefile")


def usage():
    print('')
    print('Usage: CatchRiver [-d dem (raster)] [-l burn_line (shape)] [-p burn_point (shape)] [-a burn_area (shape)]\n '
          '[-R riverout (shape)] [-C catchmentout (shape)]  [-O min strahler order (integer)] -B -S -K')
    print('-d   digital elevation model (GeoTiff)')
    print('-l   polylines (rivers) to be burned in the DEM (optional) (ESRI Shapefile)')
    print('-p   points (outlets) to be burned in the DEM (optional) (ESRI Shapefile)')
    print('-F   factor by which the cell-size should be scaled (default=1)')
    print('-O   minimal strahler order in river shapefile (optional, default=3) (integer)')
    print('-s   option to snap points (-p) to lines (-l) (default=no)')
    print('-R   name of output river (optional, default = river.shp) (ESRI Shapefile)')
    print('-C   name of output catchment (optional, default = catchment.shp) (ESRI Shapefile)')
    print('-B   burn value by which DEM will be lowered if burned (optional, default=1000) (integer)')
    print('-S   option to skip generation of LDD (default=no)')
    print('-K   option to keep all catchments and river networks (default=no)')
    print('-I   option to force "lddin" (default=no)')

    print('')


def removeshp(shapein, directory):
    if os.path.exists(directory + shapein):
        print(shapein + ' exists and will be deleted')
        Driver.DeleteDataSource(directory + shapein)
        shapein = directory + shapein
    else:
        shapein = directory + shapein
    if os.path.exists(directory + shapein):
        print('failed to remove ' + directory + shapein)
        counter = 1
        stopcounting = False
        shp_att = os.path.splitext(os.path.basename(shapein))[0]
        while not stopcounting:
            filename = shp_att + "(" + str(counter) + ")" + ".shp"
            if not os.path.exists(directory + filename):
                shapein = directory + filename
                stopcounting = True
            else:
                counter += 1
        print('filename used: ' + shapein)
    return shapein


# def main():
workdir = "work\\"
resultdir = "CatchRiver\\"

""" read commandline arguments """
argv = sys.argv
# argv = ['x','-d', '..\input\DEM_5M_filled.tif','-F','100','-R','river_test.shp',]


try:
    opts, args = getopt.getopt(argv[1:], "d:l:p:F:a:R:C:O:B:SsKI")
except getopt.error:
    print('error')
    usage()
    sys.exit(1)

dem_in = None
rivshp = None
catchshp = None
lineshp = None
areashp = None
pointshp = None
minorder = None
snapgaugestoriver = False
scalefactor = None
skipldd = False
EPSG = None
srs = None
keepall = False
lddin = False
burnvalue = 1000

for o, a in opts:
    if o == "-d":
        dem_in = a
    if o == "-l":
        lineshp = a
    if o == "-F":
        scalefactor = a
    if o == "-p":
        pointshp = a
    if o == "-R":
        rivshp = a
    if o == "-C":
        catchshp = a
    if o == "-O":
        minorder = int(a)
    if o == "-S":
        skipldd = True
    if o == "-s":
        snapgaugestoriver = True
    if o == "-B":
        burnvalue = float(a)
    if o == "-K":
        keepall = True
    if o == "-I":
        lddin = True

""" check if files exist """
if dem_in == None:
    if not skipldd:
        print('please provide dem')
        usage()
        sys.exit(1)
else:
    if not os.path.exists(dem_in):
        print('file ' + dem_in)
        print('Your DEM does not exist in the file-system')
        print('')
        sys.exit(1)
if not pointshp == None:
    if not os.path.exists(pointshp):
        print('file ' + pointshp)
        print('Your point-shape does not exist in the file-system')
        print('')
        sys.exit(1)
if not lineshp == None:
    if not os.path.exists(lineshp):
        print('file ' + lineshp)
        print('Your line-shape does not exist in the file-system')
        print('')
        sys.exit(1)

""" set property values """
if minorder == None:
    print('no minimum strahler order specified')
    print('default will be used: 5')
    minorder = int(5)

if burnvalue == None:
    print('no value for burning defined')
    print('default will be used: 1000 (map units)')
    print('pits will be filled till 500 (map units)')
    burnvalue = float(1000)

if rivshp == None:
    print('default name for river shape will be used: river.shp')
    rivshp = 'river.shp'

if catchshp == None:
    print('default name for river shape will be used: catchment.shp')
    catchshp = 'catchments.shp'

if not dem_in == None:
    ds = gdal.Open(dem_in)
    if ds == None:
        print('Input file specified not available or not a raster')
        sys.exit(1)
    else:
        spatialref = ds.GetProjection()
        srs = osr.SpatialReference()
        if (srs == None) or (spatialref == ''):
            print('Your DEM is not projected')
            sys.exit(1)
        print(srs)
        srs.ImportFromWkt(spatialref)
        srs.AutoIdentifyEPSG()
        EPSG = "EPSG:" + srs.GetAttrValue("AUTHORITY", 1)
        cellsize = ds.GetGeoTransform()[1]
        # transform = ds.GetGeoTransform()
    # ds = None
else:
    print('no DEM provided, no projection will be assigned to output')

""" create directories """
if not skipldd:
    if os.path.isdir(workdir):
        try:
            shutil.rmtree(workdir)
        except:
            print('cannot remove work directory')
            print('probably blocked by other process')
            sys.exit(1)
    os.makedirs(workdir)

if not os.path.isdir(resultdir):
    os.makedirs(resultdir)

rivshp = removeshp(rivshp, resultdir)
catchshp = removeshp(catchshp, resultdir)


""" convert and read DEM """
if not skipldd:
    dem_map = workdir + "dem.map"
    if not scalefactor == None:
        cellsizescaled = float(cellsize) * float(scalefactor)
        dem_scaled = workdir + "dem_scaled.tif"
        call(
            (
                "gdalwarp",
                "-overwrite",
                "-s_srs",
                EPSG,
                "-t_srs",
                EPSG,
                "-tr",
                str(cellsizescaled),
                str(-cellsizescaled),
                "-dstnodata",
                str(-9999),
                "-r",
                "cubic",
                dem_in,
                dem_scaled,
            )
        )
        dem_in = dem_scaled
    call(
        (
            "gdal_translate",
            "-of",
            "PCRaster",
            "-a_srs",
            EPSG,
            "-ot",
            "Float32",
            dem_in,
            dem_map,
        )
    )
    dem = pcr.readmap(dem_map)
    lines = dem * 0
    points = dem * 0
    # create mask (if needed)
    burndem = False
    if not (lineshp == None and areashp == None and pointshp == None):
        clone_map = workdir + "clone.map"
        clone = dem * 0
        burn = pcr.cover(dem * 0, pcr.scalar(0))
        # pcr.report(burn,'burn1.map')
        pcr.report(clone, clone_map)
        burndem = True
    # burn lines
    if not lineshp == None:
        file_att = os.path.splitext(os.path.basename(lineshp))[0]
        line_tif = workdir + "line.tif"
        line_map = workdir + "line.map"
        call(
            (
                "gdal_translate",
                "-of",
                "GTiff",
                "-a_srs",
                EPSG,
                "-ot",
                "Float32",
                clone_map,
                line_tif,
            )
        )
        call(("gdal_rasterize", "-burn", "1", "-l", file_att, lineshp, line_tif))
        call(
            (
                "gdal_translate",
                "-of",
                "PCRaster",
                "-a_srs",
                EPSG,
                "-ot",
                "Float32",
                line_tif,
                line_map,
            )
        )
        lines = pcr.scalar(pcr.readmap(line_map))
        burn = burn - (pcr.scalar(lines) * pcr.scalar(burnvalue))
        # pcr.report(burn,'burn2.map')
    # burn points
    if not pointshp == None:
        file_att = os.path.splitext(os.path.basename(pointshp))[0]
        point_tif = workdir + "point.tif"
        point_map = workdir + "point.map"
        call(
            (
                "gdal_translate",
                "-of",
                "GTiff",
                "-a_srs",
                EPSG,
                "-ot",
                "Float32",
                clone_map,
                point_tif,
            )
        )
        call(("gdal_rasterize", "-burn", "1", "-l", file_att, pointshp, point_tif))
        call(
            (
                "gdal_translate",
                "-of",
                "PCRaster",
                "-a_srs",
                EPSG,
                "-ot",
                "Float32",
                point_tif,
                point_map,
            )
        )
        points = pcr.scalar(pcr.readmap(point_map))
        if snapgaugestoriver:
            print("Snapping points to line")
            points = wt.snaptomap(pcr.ordinal(points), pcr.boolean(lines))
            points = pcr.cover(pcr.scalar(points), pcr.scalar(0))
        points = pcr.cover(points, pcr.scalar(0))
        # pcr.report(points,'points.map')
        burn = burn - (points * pcr.scalar(burnvalue) * 2)
        # pcr.report(burn,'burn3.map')

""" create ldd """
pcr.setglobaloption("lddout")
if lddin:
    pcr.setglobaloption("lddin")
ldd_map = workdir + "ldd.map"
streamorder_map = workdir + "streamorder.map"
river_map = workdir + "river.map"
catchments_map = workdir + "catchments.map"
catchments_tif = workdir + "catchments.tif"
# catchments_shp = resultdir + 'catchments.shp'

generateldd = True

if skipldd:
    print('Option -S is set')
    print('ldd will be read from ' + ldd_map)
    if os.path.exists(ldd_map):
        ldd = pcr.ldd(pcr.readmap(ldd_map))
        generateldd = False
    else:
        print('file ' + ldd_map + ' does not exist')
        print('new ldd will be generated')

if generateldd:
    print('Generating ldd...')
    if burndem:
        linescover = pcr.ifthen(lines == 1, pcr.scalar(0))
        pointscover = pcr.ifthen(pcr.scalar(points) == 1, pcr.scalar(0))
        # pcr.report(linescover,'lines.map')
        # pcr.report(pointscover,'points.map')
        dem = pcr.cover(dem, linescover, pointscover)
        # pcr.report(dem,'dem1.map')
        dem = dem + burn
        # pcr.report(dem,'dem2.map')
        ldd = pcr.lddcreate(
            dem, float("1E35"), float("1E35"), float("1E35"), float("1E35")
        )
    else:
        ldd = pcr.lddcreate(
            dem, burnvalue / 2, float("1E35"), float("1E35"), float("1E35")
        )

streamorder = pcr.ordinal(pcr.streamorder(ldd))
river = pcr.boolean(
    pcr.ifthen(
        streamorder >= int(min(np.max(pcr.pcr2numpy(streamorder, -9999)), minorder)),
        streamorder,
    )
)
outlets = pcr.ifthen(pcr.ordinal(ldd) == 5, pcr.boolean(1))
outlets = pcr.nominal(pcr.uniqueid(outlets))
catchments = pcr.nominal(pcr.catchment(ldd, outlets))

if not keepall:
    catchments = pcr.nominal(
        pcr.ifthen(
            pcr.mapmaximum(
                pcr.areatotal(pcr.scalar(catchments) * 0 + 1, pcr.nominal(catchments))
            )
            == pcr.areatotal(pcr.scalar(catchments) * 0 + 1, pcr.nominal(catchments)),
            catchments,
        )
    )

pcr.report(ldd, ldd_map)
pcr.report(streamorder, streamorder_map)
pcr.report(river, river_map)
pcr.report(catchments, catchments_map)
if not EPSG == None:
    call(
        (
            "gdal_translate",
            "-of",
            "GTiff",
            "-stats",
            "-a_srs",
            EPSG,
            "-ot",
            "Float32",
            catchments_map,
            catchments_tif,
        )
    )
else:
    call(
        (
            "gdal_translate",
            "-of",
            "GTiff",
            "-stats",
            "-ot",
            "Float32",
            catchments_map,
            catchments_tif,
        )
    )
wt.Raster2Pol(catchments_tif, catchshp, srs)

riversid_map = workdir + "riverid.map"
drain_map = workdir + "drain.map"
ldd_mask = pcr.ifthen(river, ldd)
upstream = pcr.upstream(ldd_mask, pcr.scalar(river))
downstream = pcr.downstream(ldd_mask, upstream)
# pcr.report(downstream,'downstream.map')
confluences = pcr.boolean(pcr.ifthen(downstream >= 2, pcr.boolean(1)))
# pcr.report(confluences,'confluences.map')
boundaries = pcr.boolean(pcr.ifthen(pcr.scalar(ldd_mask) == 5, pcr.boolean(1)))
catch_points = pcr.nominal(pcr.uniqueid(pcr.cover(confluences, boundaries)))
catchmentsid = pcr.nominal(pcr.subcatchment(ldd, catch_points))
drain = pcr.accuflux(ldd_mask, 1)
riversid = pcr.ifthen(river, catchmentsid)

if not keepall:
    riversid = pcr.nominal(
        pcr.ifthen(
            pcr.mapmaximum(
                pcr.areatotal(pcr.scalar(catchments) * 0 + 1, pcr.nominal(catchments))
            )
            == pcr.areatotal(pcr.scalar(catchments) * 0 + 1, pcr.nominal(catchments)),
            riversid,
        )
    )

pcr.report(riversid, riversid_map)
pcr.report(drain, drain_map)

print('converting river map-file to shape-file...')
wt.PCR_river2Shape(riversid_map, drain_map, streamorder_map,
                   ldd_map, rivshp, catchments_map, srs)
# if __name__ == "__main__":
#    main()
