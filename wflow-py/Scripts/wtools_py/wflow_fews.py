
import create_grid as cg
import static_maps as sm
import glob
import os
import pcraster as pcr
import zipfile
from lxml import etree

# prameters From extent_in:

source = '.\wflow_catchment\extents\extent_in.xml'
#source = 'd:\FEWS\FEWS_Accelerator\Mekong_SA\Modules\wflow\wflow_catchment\extents\extent_in.xml'

tree = etree.parse(source)
doc = tree.getroot()
model_name = doc.attrib['locationId']
rows = float(doc[0].text)
cols = float(doc[1].text)
cellsize = float(doc[4].text)   
xmin = float(doc[3][0].text) - 0.5 * cellsize
ymin = float(doc[3][1].text) + 0.5 * cellsize
xmax = xmin + cols * cellsize
ymax = ymin + rows * cellsize 


#rename dir and change current directory
os.rename(os.path.join(os.getcwd(),'wflow_catchment'),os.path.join(os.getcwd(),model_name))
os.chdir(os.path.join(os.getcwd(),model_name))

#get base-data from GEE

'''
Martijn and Genna, do your thing!!!!

ToDo:
   - Get the .\data\catchments\catchments.shp (I think GEOJSON is also fine) from GEE
   - Get the .\data\rivers\rivers.shp (I think GEOJSON is also fine) from GEE
   - Get the .\data\dem\dem.tif from GEE
   - Get the .\data\parameters\ directory filled from GEE. *.tif is GeoTiff. LAI00000.* is PCRaster map
            .\data\parameters\FirstZoneCapacity.tif
            .\data\parameters\FirstZoneKsatVer.tif
            .\data\parameters\ FirstZoneMinCapacity.tif
            .\data\parameters\InfiltCapSoil.tif
            .\data\parameters\ M.tif
            .\data\parameters\PathFrac.tif
            .\data\parameters\thetaS.tif
            .\data\parameters\WaterFrac.tif
            .\data\parameters\wflow_landuse.tif
            .\data\parameters\wflow_soil.tif
            .\data\parameters\clim\LAI00000.001
            .\data\parameters\clim\LAI00000.002
            .\data\parameters\clim\LAI00000.003
            .\data\parameters\clim\LAI00000.004
            .\data\parameters\clim\LAI00000.005
            .\data\parameters\clim\LAI00000.006
            .\data\parameters\clim\LAI00000.007
            .\data\parameters\clim\LAI00000.008
            .\data\parameters\clim\LAI00000.009
            .\data\parameters\clim\LAI00000.010
            .\data\parameters\clim\LAI00000.011
            .\data\parameters\clim\LAI00000.012
'''

#create grid
logfilename = 'wtools_create_grid.log'
destination = '.\mask'
inputfile = '.\data\catchments\catchments.shp'
projection = 'EPSG:4326'

cg.main(logfilename,destination,inputfile,projection,cellsize,snap=True,locationid=model_name)


# create static maps
source = r'.\mask'
destination = r'.\staticmaps'
inifile = 'staticmaps.ini'
dem_in = r'.\data\dem\dem.tif'
rivshp = r'.\data\rivers\rivers.shp'
catchshp = r'.\data\catchments\catchments.shp'
lai = r'.\data\parameters\clim'
other_maps=glob.glob(os.path.join(os.getcwd(),'data\parameters','*.tif'))

sm.main(source,destination,inifile,dem_in,rivshp,catchshp,lai=lai,other_maps=other_maps)

# save default state-files in FEWS-config
state_dir = r'.\outstate'
state_files = ['CanopyStorage.map','SatWaterDepth.map','Snow.map','SnowWater.map','SurfaceRunoff.map','TSoil.map','UStoreLayerDepth_0.map','WaterLevel.map']
zip_loc = r'..\\..\\..\\Config\\ColdStateFiles\\' + model_name[:2].upper() + model_name[2:] + '_GA_Historical default.zip'

mask = pcr.readmap(os.path.join(os.getcwd(),source,'mask.map'))

zf = zipfile.ZipFile(zip_loc, mode='w')

for state_file in state_files:
    state_path = os.path.join(os.getcwd(),state_dir,state_file)
    pcr.report(pcr.cover(mask,pcr.scalar(0)),state_path)
    zf.write(state_path, state_file, compress_type=zipfile.ZIP_DEFLATED)

zf.close()   