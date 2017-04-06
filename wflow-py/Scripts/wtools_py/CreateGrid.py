from subprocess import call
import os
import sys
import math
import getopt
import string
import shutil
from osgeo import osr
from osgeo import ogr
from lxml import etree
try:
    from osgeo import ogr
except ImportError:
    import ogr
try:
    from osgeo import gdal
    from osgeo.gdalconst import *
except ImportError:
    import gdal
    from gdalconst import *
import wflow.wflowtools_lib as wt

def Usage():
    print('')             
    print('Usage: CreateGrid [-f inputfile] [-c cell-size]')
    print '-f   file of which extent will be read. Most logically the catchment layer'
    print '     format: ESRI Shapefile or any gdal supported raster format (preferred GeoTiff)'
    print('')
    
def main():

    ''' parse command line arguments '''
    
    argv = sys.argv
    #argv = ['x','-f', 'Citarum.shp','-c','0.005']
       
    inputfile = None
    extent = None
    cellsize = None
    Projection = False
    try:
        opts, args = getopt.getopt(argv[1:], 'f:e::::c:')
    except getopt.error:
        print 'fout'
        Usage()
        sys.exit(1)
    
    for o, a in opts:
        if o == '-h':
         Usage()
         sys.exit()        
        elif o == '-f': inputfile = a
        elif o == '-e': extent = a
        elif o == '-c': cellsize = a
    
    ''' checks on existence if input '''
    if inputfile == None:
        if extent == None:
            print 'input file or extent need to be specified'
            Usage()
            sys.exit(1)
        else:
            extent_in = [float(x) for x in string.split(extent,',')]
    if cellsize == None:
        print 'cell-size (-c) is not specified'
        Usage()
        sys.exit(1)
    
    ''' delete existing files '''
    if os.path.isdir('mask'):
        shutil.rmtree('mask')
    os.makedirs('mask')
    
    if not inputfile == None:
        fileext = os.path.splitext(os.path.basename(inputfile))[1]
        if fileext == '.shp':
            '''ds openen als shape-file'''
            file_att = os.path.splitext(os.path.basename(inputfile))[0]
            ds = ogr.Open(inputfile)
            '''read extent'''
            lyr = ds.GetLayerByName(file_att)
            extent = lyr.GetExtent()
            extent_in = [extent[0],extent[2],extent[1],extent[3]]
            spatialref = lyr.GetSpatialRef()
            if not spatialref == None:
                srs = osr.SpatialReference()
                srs.ImportFromWkt(spatialref.ExportToWkt())
                proj = spatialref.GetUTMZone()
                Projection = True
                zone = str(abs(proj))
                if str(proj)[0] == '-':
                    hemisphere = 'S'
                else: hemisphere = 'N'
                geodatum = 'UTM' + zone + hemisphere
                if proj == 0:
                    geodatum = 'WGS 1984'
        else:
            ds = gdal.Open(inputfile)      
            if ds == None:
                print 'Input file specified not a shape-file or supported gdal raster format'
                sys.exit(1)
            '''read extent'''
            geotransform = ds.GetGeoTransform()
            raster_cellsize = geotransform[1]
            ncols = ds.RasterXSize
            nrows =ds.RasterYSize            
            extent_in = [geotransform[0],geotransform[3]-nrows*raster_cellsize,geotransform[0]+ncols*raster_cellsize,geotransform[3]]
            spatialref = ds.GetProjection()
            if not spatialref == None:
                srs = osr.SpatialReference()
                srs.ImportFromWkt(spatialref)
                if 'UTM zone' in spatialref:
                    Projection = True
                    zone = spatialref[spatialref.find('UTM zone')+9:spatialref.find('UTM zone')+11]
                    hemisphere = spatialref[spatialref.find('UTM zone')+11:spatialref.find('UTM zone')+12]
                    geodatum = 'UTM' + zone + hemisphere
                elif '"WGS 84"' in spatialref:
                    Projection = True                
                    geodatum = 'WGS 1984'
                    
    ''' warn if projection is not defined '''
    if not Projection:
        print 'No projection defined in input file:'
        print ' no projection will be assigned to mask.shp'
        print ' no projection will be assigned to geoDatum element in grid.xml'
        print ''
        geodatum = ''        
        
    ''' create mask.map '''
    snap = len(str(float(cellsize)-math.floor(float(cellsize))))-2
    extent_out = wt.round_extent(list(extent_in),float(cellsize),snap)
    columns = int((extent_out[2]-extent_out[0])/float(cellsize)+2)
    rows = int((extent_out[3]-extent_out[1])/float(cellsize)+2)
    cells = rows*columns
    xorg = extent_out[0]-float(cellsize)
    yorg = extent_out[3]+float(cellsize)
    
    ## TODO use gdal_writemap
    call(('mapattr','-s','-R',str(rows),'-C',str(columns),'-x',str(xorg),'-y',str(yorg),'-l',cellsize,'mask\mask.map'))
    print 'mask.map created'
        
    ''' create grid.xml '''
    root = etree.Element("regular", locationId='wflow_mask')
    etree.SubElement(root, 'rows').text = str(rows)
    etree.SubElement(root, 'columns').text = str(columns)
    etree.SubElement(root, 'geoDatum').text = geodatum
    etree.SubElement(root, 'firstCellCenter')
    etree.SubElement(root[3], 'x').text = str(xorg+0.5*float(cellsize))
    etree.SubElement(root[3], 'y').text = str(yorg-0.5*float(cellsize))
    etree.SubElement(root, 'xCellSize').text = str(cellsize)
    etree.SubElement(root, 'yCellSize').text = str(cellsize)
    gridxml = open('mask\grid.xml','w+')
    gridxml.write(etree.tostring(root, pretty_print=True))
    gridxml.close()
    
    print 'grid.xml created'
    
    ''' Create shape-file '''
    Driver = ogr.GetDriverByName("ESRI Shapefile")
    SHP_srs = "mask\mask.shp"
    SHP_srs_ATT = os.path.splitext(os.path.basename(SHP_srs))[0]
    SHP_srs_out = Driver.CreateDataSource(SHP_srs)
    if Projection:
        SHP_srs_LYR  = SHP_srs_out.CreateLayer(SHP_srs_ATT, srs, geom_type=ogr.wkbPolygon)
    else: SHP_srs_LYR  = SHP_srs_out.CreateLayer(SHP_srs_ATT, geom_type=ogr.wkbPolygon)
    fieldDef = ogr.FieldDefn("ID", ogr.OFTString)
    fieldDef.SetWidth(12)
    SHP_srs_LYR.CreateField(fieldDef)
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(xorg,yorg)
    ring.AddPoint(xorg + columns * float(cellsize),yorg)
    ring.AddPoint(xorg + columns * float(cellsize),yorg - rows * float(cellsize))
    ring.AddPoint(xorg,yorg - rows * float(cellsize))
    ring.AddPoint(xorg,yorg)
    ring.CloseRings
    polygon = ogr.Geometry(ogr.wkbPolygon)
    polygon.AddGeometry(ring)
    feat_out = ogr.Feature(SHP_srs_LYR.GetLayerDefn())
    feat_out.SetGeometry(polygon)
    feat_out.SetField("ID", 'wflow_mask')
    SHP_srs_LYR.CreateFeature(feat_out) 
    SHP_srs_out.Destroy()
    print 'mask.shp created'
    print ''
    print 'Number of cells: ' + str(cells)
    if cells > 1000000:
        print 'With this amount of cells your model will run VERY slow. Consider a larger cell-size'
        print 'fast models run with < 100000 cells'
    elif cells > 100000:
        print 'With this amount of cells your model will run slow. Consider a larger cell-size'
        print 'fast models run with < 100000 cells'
if __name__ == "__main__":
    main()

# TODO: remove hard folder mask from this script and staticmaps
# TODO: add extent
# TODO: add logger
# TODO: replace mapattr by gdal
# TODO: add projection code as optional argument and add projection to .shp