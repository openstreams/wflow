# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 13:26:48 2015

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
import glob
from osgeo import ogr
from osgeo import gdal
from osgeo.gdalconst import *

import wflow.wflowtools_lib as wt


Driver = ogr.GetDriverByName("ESRI Shapefile")


def Usage():
    print('')
    print('Usage: Staticmaps [-i ini-file] [-r river] [-g gauges (shape)] [-c catchment (shape)] [-d dem (raster)] [-l landuse (raster)] [-s soiltype (raster)] -C')
    print ''
    print '-i   ini file with settings for staticmaps.exe'
    print '-r   river network polyline layer (ESRI Shapefile)'
    print '-g   gauges point layer (ESRI Shapefile)'
    print '-c   catchment polygon layer (ESRI Shapefile)'
    print '-d   digital elevation model (GeoTiff)'
    print '-l   land-use land-cover layer (GeoTiff)'
    print '-s   soil-type layer (GeoTiff)'
    print '-C   removes the .xml files from the staticmaps directory when finished'
    print '-A   option to burn catchments "all touching".'
    print 'Usuefull when catchment-size is small compared to cellsize'
    print ''
    print '!!NOTE: PREFERABLY ALL LAYERS SHOULD BE PROJECTED!!'
    print '!!IF NOT, THE PROJECTION OF MASK WILL BE USED!!'
    print ''


# def main():
clone_map = "d:\FEWS\FEWS-Kelantan\Kelantan\Modules\wflow\kelantan\mask\mask.map"
clone_shp = "d:\FEWS\FEWS-Kelantan\Kelantan\Modules\wflow\kelantan\mask\mask.shp"
clone_prj = "d:\FEWS\FEWS-Kelantan\Kelantan\Modules\wflow\kelantan\mask\mask.prj"
workdir = "work\\"
resultdir = "staticmaps\\"

''' read commandline arguments '''
argv = sys.argv
argv = ['', '-i', 'd:\FEWS\FEWS-Kelantan\Kelantan\Modules\wflow\kelantan\StaticMaps.ini', '-d', 'd:\FEWS\FEWS-Kelantan\Kelantan\Modules\data\DEM\srtm_kelantan_UTM48N_90m.tif', '-r',
        'd:\FEWS\FEWS-Kelantan\Kelantan\Modules\data\River_Network\Kelantan_Rivers_Hydrosheds_UTM48N.shp', '-c', 'd:\FEWS\FEWS-Kelantan\Kelantan\Modules\data\River_Catchment\Kelantan_Catchment_Hydrosheds_UTM48N.shp', '-C', '-A']
clone_EPSG = False


try:
    opts, args = getopt.getopt(argv[1:], 'i:g:p:r:c:d:l:s:CA')
except getopt.error:
    print 'fout'
    Usage()
    sys.exit(1)

inifile = None
rivshp = None
catchshp = None
dem_in = None
landuse = None
soiltype = None
clean = False
gaugeshp = None
alltouching = False

for o, a in opts:
    if o == '-i':
        inifile = a
    if o == '-p':
        clone_EPSG = 'EPSG:' + a
    if o == '-r':
        rivshp = a
    if o == '-c':
        catchshp = a
    if o == '-d':
        dem_in = a
    if o == '-l':
        landuse = a
    if o == '-s':
        soiltype = a
    if o == '-C':
        clean = True
    if o == '-g':
        gaugeshp = a
    if o == '-A':
        alltouching = True

if inifile == None or rivshp == None or catchshp == None or dem_in == None:
    print 'the following files are compulsory:'
    print ' - ini-file'
    print ' - DEM (raster)'
    print ' - river (shape)'
    print ' - catchment (shape)'
    Usage()
    sys.exit(1)

if landuse == None:
    print 'no raster with landuse classifications is specified. 1 class will be applied for the entire domain'

if soiltype == None:
    print 'no raster with soil classifications is specified. 1 class will be applied for the entire domain'

''' read mask '''
if not os.path.exists(clone_map):
    print 'Mask not found. Make sure the file mask\mask.map exists'
    print 'This file is usually created with the CreateGrid script'
    sys.exit(1)
else:
    pcr.setclone(clone_map)
    ds = gdal.Open(clone_map, GA_ReadOnly)
    clone_trans = ds.GetGeoTransform()
    cellsize = clone_trans[1]
    clone_rows = ds.RasterYSize
    clone_columns = ds.RasterXSize
    extent_mask = [clone_trans[0], clone_trans[3] - ds.RasterYSize *
                   cellsize, clone_trans[0] + ds.RasterXSize * cellsize, clone_trans[3]]
    xmin, ymin, xmax, ymax = map(str, extent_mask)
    ds = None
    ones = pcr.scalar(pcr.readmap(clone_map))
    zeros = ones * 0
    empty = pcr.ifthen(ones == 0, pcr.scalar(0))

''' read projection from mask.shp '''
# TODO: check how to deal with projections (add .prj to mask.shp in creategrid)
if not os.path.exists(clone_prj):
    print 'please add prj-file to mask.shp'
    sys.exit(1)
if os.path.exists(clone_shp):
    ds = ogr.Open(clone_shp)
    file_att = os.path.splitext(os.path.basename(clone_shp))[0]
    lyr = ds.GetLayerByName(file_att)
    spatialref = lyr.GetSpatialRef()
    if not spatialref == None:
        srs_clone = osr.SpatialReference()
        srs_clone.ImportFromWkt(spatialref.ExportToWkt())
        srs_clone.AutoIdentifyEPSG()
        unit_clone = False
        unit_clone = srs_clone.GetAttrValue('UNIT').lower()
        clone_EPSG = 'EPSG:' + srs_clone.GetAttrValue("AUTHORITY", 1)
        # TODO: fix hard EPSG code below
        #clone_EPSG = 'EPSG:'+ '4167'
        print 'EPSG-code is read from mask.shp: ' + clone_EPSG
        spatialref == None
if not clone_EPSG:
    print 'EPSG-code cannot be read from mask.shp'
    print 'please add prj-file to mask.shp or specify on command line'
    print 'e.g. -p EPSG:4326 (for WGS84 lat lon projection)'

ds = None
clone_EPSG_int = int(clone_EPSG[5:len(clone_EPSG)])

''' open config-file '''
config = wt.OpenConf(inifile)

''' read settings '''
snapgaugestoriver = bool(int(wt.configget(config,
                                          "settings",
                                          "snapgaugestoriver", "1")))
burnalltouching = bool(int(wt.configget(config,
                                        "settings",
                                        "burncatchalltouching",
                                        "1")))
burninorder = bool(
    int(wt.configget(config, "settings", "burncatchalltouching", "0")))
verticetollerance = float(wt.configget(
    config, "settings", "vertice_tollerance", "0.0001"))

''' read parameters '''
burn_outlets = int(wt.configget(config, "parameters", "burn_outlets", 10000))
burn_rivers = int(wt.configget(config, "parameters", "burn_rivers", 200))
burn_connections = int(wt.configget(
    config, "parameters", "burn_connections", 100))
burn_gauges = int(wt.configget(config, "parameters", "burn_gauges", 100))
minorder = int(wt.configget(config, "parameters", "riverorder_min", 3))
exec "percentile=tr.array(" + wt.configget(config, "parameters", "statisticmaps", [0, 100]) + ")"
if not unit_clone:
    print 'failed to read unit (meter or degree) from mask projection'
    unit_clone = str(wt.configget(config, "settings", "unit", 'meter'))
    print 'unit read from settings: ' + unit_clone
if unit_clone == 'degree':
    cellsize_hr = float(wt.configget(
        config, "parameters", "highres_degree", 0.0005))
elif (unit_clone == 'metre') or (unit_clone == 'meter'):
    cellsize_hr = float(wt.configget(
        config, "parameters", "highres_metre", 50))

cols_hr = int((float(xmax) - float(xmin)) / cellsize_hr + 2)
rows_hr = int((float(ymax) - float(ymin)) / cellsize_hr + 2)
hr_trans = (float(xmin), cellsize_hr, float(0), float(ymax), 0, -cellsize_hr)

''' read staticmap locations '''
catchment_map = wt.configget(
    config, "staticmaps", "catchment", "wflow_catchment.map")
dem_map = wt.configget(config, "staticmaps", "dem", "wflow_dem.map")
demmax_map = wt.configget(config, "staticmaps", "demmax", "wflow_demmax.map")
demmin_map = wt.configget(config, "staticmaps", "demmin", "wflow_demmin.map")
gauges_map = wt.configget(config, "staticmaps", "gauges", "wflow_gauges.map")
landuse_map = wt.configget(
    config, "staticmaps", "landuse", "wflow_landuse.map")
ldd_map = wt.configget(config, "staticmaps", "ldd", "wflow_ldd.map")
river_map = wt.configget(config, "staticmaps", "river", "wflow_river.map")
outlet_map = wt.configget(config, "staticmaps", "outlet", "wflow_outlet.map")
riverlength_fact_map = wt.configget(
    config, "staticmaps", "riverlength_fact", "wflow_riverlength_fact.map")
soil_map = wt.configget(config, "staticmaps", "soil", "wflow_soil.map")
streamorder_map = wt.configget(
    config, "staticmaps", "streamorder", "wflow_streamorder.map")
subcatch_map = wt.configget(
    config, "staticmaps", "subcatch", "wflow_subcatch.map")


''' read mask location (optional) '''
masklayer = wt.configget(config, "mask", "masklayer", catchshp)

''' create directories '''
if os.path.isdir(workdir):
    shutil.rmtree(workdir)
os.makedirs(workdir)

if os.path.isdir(resultdir):
    shutil.rmtree(resultdir)
os.makedirs(resultdir)

''' Preperation steps '''
zero_map = workdir + "zero.map"
zero_tif = workdir + "zero.tif"
pcr.report(zeros, zero_map)
# TODO: replace gdal_translate call
call(('gdal_translate', '-of', 'GTiff', '-a_srs',
      clone_EPSG, '-ot', 'Float32', zero_map, zero_tif))
pcr.setglobaloption("lddin")

''' resample DEM '''
dem_resample = workdir + "dem_resampled.tif"
ds = gdal.Open(dem_in, GA_ReadOnly)
band = ds.GetRasterBand(1)
nodata = band.GetNoDataValue()
proj = ds.GetGeoTransform()
cellsize_dem = proj[1]

''' read DEM projection '''
spatialref == None
spatialref = ds.GetProjection()
if not spatialref == None:
    srs = osr.SpatialReference()
    srs.ImportFromWkt(spatialref)
    srs.AutoIdentifyEPSG()
    dem_EPSG = 'EPSG:' + srs.GetAttrValue("AUTHORITY", 1)
    print 'EPSG-code is read from ' + os.path.basename(dem_in) + ': ' + dem_EPSG
    spatialref == None
    dem_EPSG_int = int(dem_EPSG[5:len(dem_EPSG)])
    srs_DEM = osr.SpatialReference()
    srs_DEM.ImportFromEPSG(dem_EPSG_int)
    clone2dem_transform = osr.CoordinateTransformation(srs_clone, srs_DEM)
else:
    dem_EPSG = clone_EPSG
    print 'No projection defined for ' + os.path.basename(dem_in)
    print 'Assumed to be the same as model projection (' + clone_EPSG + ')'

ds = None
print 'Resampling DEM...'
if nodata == None:
    call(('gdalwarp', '-overwrite', '-t_srs', clone_prj, '-te', xmin, ymin, xmax, ymax, '-tr',
          str(cellsize), str(-cellsize), '-dstnodata', str(-9999), '-r', 'cubic', dem_in, dem_resample))
else:
    call(('gdalwarp', '-overwrite', '-t_srs', clone_prj, '-te', xmin, ymin, xmax, ymax, '-tr', str(cellsize),
          str(-cellsize), '-srcnodata', str(nodata), '-dstnodata', str(nodata), '-r', 'cubic', dem_in, dem_resample))

''' create dem.map and statistic maps '''
dem_resample_map = resultdir + dem_map
call(('gdal_translate', '-of', 'PCRaster', '-a_srs', clone_EPSG,
      '-ot', 'Float32', dem_resample, dem_resample_map))
print 'Computing DEM statistics ....'
stats = wt.windowstats(dem_in, clone_rows, clone_columns,
                       clone_trans, srs_clone, resultdir, percentile)

''' burn DEM '''
ds = ogr.Open(rivshp)
file_att = os.path.splitext(os.path.basename(rivshp))[0]
lyr = ds.GetLayerByName(file_att)
spatialref = lyr.GetSpatialRef()
#    if not spatialref == None:
#        srs = osr.SpatialReference()
#        srs.ImportFromWkt(spatialref.ExportToWkt())
#        srs.AutoIdentifyEPSG()
#        rivshp_EPSG = 'EPSG:'+srs.GetAttrValue("AUTHORITY",1)
#        spatialref == None
#    else:
rivshp_EPSG = clone_EPSG
print 'No projection defined for ' + file_att + '.shp'
print 'Assumed to be the same as model projection (' + clone_EPSG + ')'

# strip rivers to nodes
xminc = str(float(xmin) + 0.5 * cellsize)
yminc = str(float(ymin) + 0.5 * cellsize)
xmaxc = str(float(xmax) - 0.5 * cellsize)
ymaxc = str(float(ymax) - 0.5 * cellsize)
if rivshp_EPSG == clone_EPSG:
    rivclipshp = workdir + 'rivshape_clip.shp'
    call(('ogr2ogr', '-s_srs', clone_EPSG, '-t_srs', clone_EPSG, '-spat', xmin,
          ymin, xmax, ymax, '-clipsrc', xminc, yminc, xmaxc, ymaxc, rivclipshp, rivshp))
else:
    rivprojshp = workdir + 'rivshape_proj.shp'
    rivclipshp = workdir + 'rivshape_clip.shp'
    call(('ogr2ogr', '-s_srs', rivshp_EPSG, '-t_srs', clone_EPSG,
          '-spat', xmin, ymin, xmax, ymax, rivprojshp, rivshp))
    call(('ogr2ogr', '-s_srs', clone_EPSG, '-t_srs', clone_EPSG, '-spat', xmin, ymin,
          xmax, ymax, '-clipsrc', xminc, yminc, xmaxc, ymaxc, rivclipshp, rivprojshp))

rivshp = rivclipshp


#### BURNING BELOW ####

# TODO: check if extraction can be done within memory and retun a burn layer
shapes = wt.Reach2Nodes(rivclipshp, clone_EPSG_int,
                        cellsize * verticetollerance, workdir)

outlets = shapes[1]
connections = shapes[2]
outlets_att = os.path.splitext(os.path.basename(outlets))[0]
connections_att = os.path.splitext(os.path.basename(connections))[0]
dem_resample_att = os.path.splitext(os.path.basename(dem_resample))[0]
connections_tif = workdir + connections_att + ".tif"
outlets_tif = workdir + outlets_att + ".tif"
# TODO: make the burning in memory
call(('gdal_translate', '-of', 'GTiff', '-a_srs', clone_EPSG,
      '-ot', 'Float32', zero_map, connections_tif))
call(('gdal_translate', '-of', 'GTiff', '-a_srs',
      clone_EPSG, '-ot', 'Float32', zero_map, outlets_tif))
call(('gdal_rasterize', '-burn', '1', '-l', outlets_att, outlets, outlets_tif))
call(('gdal_rasterize', '-burn', '1', '-l',
      connections_att, connections, connections_tif))

# convert rivers to order
rivshp_att = os.path.splitext(os.path.basename(rivshp))[0]
rivers_tif = workdir + rivshp_att + ".tif"
call(('gdal_translate', '-of', 'GTiff', '-a_srs',
      clone_EPSG, '-ot', 'Float32', zero_map, rivers_tif))
if burninorder:  # make river shape with an order attribute
    OrderSHPs = wt.ReachOrder(rivshp, clone_EPSG_int,
                              cellsize * verticetollerance, workdir)
    wt.Burn2Tif(OrderSHPs, 'order', rivers_tif)
else:
    call(('gdal_rasterize', '-burn', '1', '-l',
          rivshp_att, rivshp, rivers_tif))

# convert 2 maps
connections_map = workdir + connections_att + ".map"
rivers_map = workdir + rivshp_att + ".map"
outlets_map = workdir + outlets_att + ".map"
call(('gdal_translate', '-of', 'PCRaster', '-a_srs', clone_EPSG,
      '-ot', 'Float32', connections_tif, connections_map))
call(('gdal_translate', '-of', 'PCRaster', '-a_srs',
      clone_EPSG, '-ot', 'Float32', rivers_tif, rivers_map))
call(('gdal_translate', '-of', 'PCRaster', '-a_srs',
      clone_EPSG, '-ot', 'Float32', outlets_tif, outlets_map))

# burn the layers in DEM
outletsburn = pcr.scalar(pcr.readmap(outlets_map)) * pcr.scalar(burn_outlets)
connectionsburn = pcr.scalar(pcr.readmap(
    connections_map)) * pcr.scalar(burn_connections)
riverburn = pcr.scalar(pcr.readmap(rivers_map)) * pcr.scalar(burn_rivers)
ldddem = pcr.cover(dem_resample_map, pcr.ifthen(riverburn > 0, pcr.scalar(0)))
ldddem = ldddem - outletsburn - connectionsburn - riverburn
ldddem = pcr.cover(ldddem, pcr.scalar(0))
pcr.report(ldddem, workdir + "dem_burn.map")

''' create ldd for multi-catchments '''
ldd = pcr.ldd(empty)
# reproject catchment shape-file
ds = ogr.Open(catchshp)
file_att = os.path.splitext(os.path.basename(catchshp))[0]
lyr = ds.GetLayerByName(file_att)
spatialref = lyr.GetSpatialRef()
#    if not spatialref == None:
#        srs = osr.SpatialReference()
#        srs.ImportFromWkt(spatialref.ExportToWkt())
#        srs.AutoIdentifyEPSG()
#        catchshp_EPSG = 'EPSG:'+srs.GetAttrValue("AUTHORITY",1)
#        spatialref == None
#    else:
catchshp_EPSG = clone_EPSG
print 'No projection defined for ' + file_att + '.shp'
print 'Assumed to be the same as model projection (' + clone_EPSG + ')'

if not rivshp_EPSG == clone_EPSG:
    catchprojshp = workdir + 'catchshape_proj.shp'
    call(('ogr2ogr', '-s_srs', catchshp_EPSG,
          '-t_srs', clone_ESPG, catchprojshp, catchshp))
    catchshp = catchprojshp
ds.Destroy()

ds = ogr.Open(catchshp)
file_att = os.path.splitext(os.path.basename(catchshp))[0]
lyr = ds.GetLayerByName(file_att)

fieldDef = ogr.FieldDefn("ID", ogr.OFTString)
fieldDef.SetWidth(12)
TEMP_out = Driver.CreateDataSource(workdir + "temp.shp")
if not srs == None:
    TEMP_LYR = TEMP_out.CreateLayer("temp", srs, geom_type=ogr.wkbMultiPolygon)
else:
    TEMP_LYR = TEMP_out.CreateLayer("temp", geom_type=ogr.wkbMultiPolygon)
TEMP_LYR.CreateField(fieldDef)

for i in range(lyr.GetFeatureCount()):
    orgfeature = lyr.GetFeature(i)
    geometry = orgfeature.geometry()
    feature = ogr.Feature(TEMP_LYR.GetLayerDefn())
    feature.SetGeometry(geometry)
    feature.SetField("ID", str(i + 1))
    TEMP_LYR.CreateFeature(feature)
TEMP_out.Destroy()
ds.Destroy

# rasterize catchment map
catchments_tif = workdir + "catchments.tif"
catchments_map = workdir + "catchments.map"
call(('gdal_translate', '-of', 'GTiff', '-a_srs',
      clone_EPSG, zero_map, catchments_tif))
if alltouching:
    call(('gdal_rasterize', '-at', '-a', 'ID', '-l',
          "temp", workdir + 'temp.shp', catchments_tif))
else:
    call(('gdal_rasterize', '-a', 'ID', '-l', "temp",
          workdir + 'temp.shp', catchments_tif))
call(('gdal_translate', '-of', 'PCRaster', '-a_srs',
      clone_EPSG, catchments_tif, catchments_map))
catchments = pcr.readmap(catchments_map)
riverunique = pcr.clump(pcr.nominal(pcr.ifthen(riverburn > 0, riverburn)))
rivercatch = pcr.areamajority(pcr.ordinal(catchments), riverunique)
#catchments = pcr.cover(pcr.ordinal(rivercatch),pcr.ordinal(pcr.ifthen(catchments > 0, catchments)),pcr.ordinal(0))
catchments = pcr.cover(pcr.ifthen(catchments > 0,
                                  pcr.ordinal(catchments)),
                       pcr.ifthen(riverburn > 0,
                                  pcr.ordinal(pcr.spreadzone(pcr.nominal(catchments),
                                                             pcr.ifthen(riverburn > 0,
                                                                        pcr.scalar(1)),
                                                             1))))
rivercatch_map = workdir + "catchments_river.map"
catchclip_map = workdir + "catchments_clip.map"
pcr.report(rivercatch, rivercatch_map)
pcr.report(catchments, catchclip_map)

ds = ogr.Open(workdir + "temp.shp")
lyr = ds.GetLayerByName("temp")

print 'calculating ldd'
for i in range(lyr.GetFeatureCount()):
    feature = lyr.GetFeature(i)
    catch = int(feature.GetField("ID"))
    print "calculating ldd for catchment: " + str(i + 1) + "/" + str(lyr.GetFeatureCount()) + "...."
    ldddem_select = pcr.scalar(pcr.ifthen(
        catchments == catch, catchments)) * 0 + 1 * ldddem
    ldd_select = pcr.lddcreate(ldddem_select, float(
        "1E35"), float("1E35"), float("1E35"), float("1E35"))
    ldd = pcr.cover(ldd, ldd_select)
pcr.report(ldd, resultdir + ldd_map)
ds.Destroy()

''' report stream order, river and dem '''
streamorder = pcr.ordinal(pcr.streamorder(ldd))
river = pcr.ifthen(streamorder >= pcr.ordinal(minorder), pcr.boolean(1))
mindem = int(np.min(pcr.pcr2numpy(pcr.ordinal(dem_resample_map), 9999999)))
dem_resample_map = pcr.cover(dem_resample_map, pcr.scalar(river) * 0 + mindem)
pcr.report(dem_resample_map, resultdir + dem_map)
pcr.report(streamorder, resultdir + streamorder_map)
pcr.report(river, resultdir + river_map)

''' deal with your catchments '''
if gaugeshp == None:
    print 'No gauges defined, using outlets instead'
    gauges = pcr.ordinal(pcr.uniqueid(pcr.boolean(
        pcr.ifthen(pcr.scalar(ldd) == 5, pcr.boolean(1)))))
    pcr.report(gauges, resultdir + gauges_map)
#    ds = ogr.Open(gaugeshp)
#    file_att = os.path.splitext(os.path.basename(gaugeshp))[0]
#    lyr = ds.GetLayerByName(file_att)
#    spatialref = lyr.GetSpatialRef()
# if not spatialref == None:
##        srs = osr.SpatialReference()
# srs.ImportFromWkt(spatialref.ExportToWkt())
# srs.AutoIdentifyEPSG()
##        gaugeshp_EPSG = 'EPSG:'+srs.GetAttrValue("AUTHORITY",1)
##        spatialref == None
#    #else:
#    gaugeshp_EPSG = clone_EPSG
#    print 'No projection defined for ' + file_att + '.shp'
#    print 'Assumed to be the same as model projection (' + clone_EPSG + ')'
#
#    # reproject gauge shape if necesarry
#    if not gaugeshp_EPSG == clone_EPSG:
#        gaugeprojshp = workdir + 'gaugeshape_proj.shp'
#        call(('ogr2ogr','-s_srs',rivshp_EPSG,'-t_srs',clone_ESPG,gaugeprojshp,gaugeshp))
#        gaugeshp = gaugeprojshp
#
#    file_att = os.path.splitext(os.path.basename(gaugeshp))[0]
#    gaugestif = workdir + file_att + '.tif'
#    gaugesmap = workdir + file_att + '.map'
#    call(('gdal_translate','-of','GTiff','-a_srs',clone_EPSG,zero_map,gaugestif))
#    call(('gdal_rasterize','-burn','1','-l',file_att,gaugeshp,gaugestif))
#    call(('gdal_translate','-of','PCRaster','-a_srs',clone_EPSG,gaugestif,gaugesmap))
#    gaugelocs = pcr.readmap(gaugesmap)
#    snapgaugestoriver = True
#
#    if snapgaugestoriver:
#        print "Snapping gauges to river"
#        gauges = pcr.uniqueid(pcr.boolean(gaugelocs))
#        gauges= wt.snaptomap(pcr.ordinal(gauges),river)
#
#    gaugesmap = pcr.ifthen(gauges > 0, gauges)


''' report riverlengthfrac '''
riv_hr = workdir + 'river_highres.tif'
wt.CreateTif(riv_hr, rows_hr, cols_hr, hr_trans, srs_clone, 0)
file_att = os.path.splitext(os.path.basename(rivshp))[0]
call(('gdal_rasterize', '-burn', '1', '-l', file_att, rivshp, riv_hr))
print 'Computing river length...'
#riverlength = wt.windowstats(riv_hr,clone_rows,clone_columns,clone_trans,srs_clone,resultdir,'frac',clone2dem_transform)
riverlength = wt.windowstats(
    riv_hr, clone_rows, clone_columns, clone_trans, srs_clone, resultdir, 'frac')

''' report outlet map '''
pcr.report(pcr.ifthen(pcr.ordinal(ldd) == 5,
                      pcr.ordinal(1)), resultdir + outlet_map)

''' report  map '''
catchment = pcr.ifthen(catchments > 0, pcr.ordinal(1))
pcr.report(catchment, resultdir + catchment_map)

''' report subcatchment map '''
subcatchment = pcr.subcatchment(ldd, gauges)
pcr.report(pcr.ordinal(subcatchment), resultdir + subcatch_map)

''' report landuse map '''
if landuse == None:
    pcr.report(pcr.nominal(ones), resultdir + landuse_map)
else:
    landuse_resample = workdir + 'landuse.tif'
    landuse_map = resultdir + landuse_map
    transform = wt.GetRasterTranform(landuse, srs_clone)
    if not transform[0]:
        call(('gdalwarp', '-overwrite', '-s_srs', clone_EPSG, '-t_srs', clone_EPSG, '-te', xmin, ymin,
              xmax, ymax, '-tr', str(cellsize), str(-cellsize), '-r', 'mode', landuse, landuse_resample))
    else:
        call(('gdalwarp', '-overwrite', '-s_srs', transform[1], '-t_srs', clone_EPSG, '-te', xmin, ymin,
              xmax, ymax, '-tr', str(cellsize), str(-cellsize), '-r', 'mode', landuse, landuse_resample))
    call(('gdal_translate', '-of', 'PCRaster', '-ot',
          'Float32', landuse_resample, landuse_map))
    landuse_work = pcr.readmap(landuse_map)
    pcr.report(pcr.nominal(landuse_work), landuse_map)

''' report soil map '''
if soiltype == None:
    pcr.report(pcr.nominal(ones), resultdir + soil_map)
else:
    soiltype_resample = workdir + 'soiltype.tif'
    soil_map = resultdir + soil_map
    #transform = wt.GetRasterTranform(soiltype,srs_clone)
#        if not transform[0]:
    call(('gdalwarp', '-overwrite', '-s_srs', clone_EPSG, '-t_srs', clone_EPSG, '-te', xmin, ymin,
          xmax, ymax, '-tr', str(cellsize), str(-cellsize), '-r', 'mode', soiltype, soiltype_resample))
#        else:
#        call(('gdalwarp','-overwrite','-s_srs',transform[1],'-t_srs',clone_EPSG,'-te', xmin, ymin, xmax, ymax,'-tr',str(cellsize),str(-cellsize),'-r','mode',soiltype, soiltype_resample))
    call(('gdal_translate', '-of', 'PCRaster', '-ot',
          'Float32', soiltype_resample, soil_map))
    soiltype_work = pcr.readmap(soil_map)
    pcr.report(pcr.nominal(soiltype_work), soil_map)

if clean:
    wt.DeleteList(glob.glob(os.getcwd() + '\\' + resultdir + '/*.xml'))
# if __name__ == "__main__":
#    main()

# TODO: replace calls to absoute folders for configurable folders
# TODO: check if 'work' folder is still needed
# windowstat, replace raster write functions to gdal_writemap
