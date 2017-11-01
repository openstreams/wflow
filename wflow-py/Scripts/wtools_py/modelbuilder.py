# currently assumes the working directory is right above the wflow cases

import os
import sys
import geojson
import requests
import wflow.create_grid as cg
import wflow.static_maps as sm
import json
import shutil
import zipfile
import tempfile
from math import sqrt
from pyproj import Geod

# get from the user:
path_settings = 'settings.json'
with open(path_settings) as f:
    d = geojson.load(f)

if not d.is_valid:
    raise AssertionError(
        '{} is not a valid GeoJSON file\n{}'.format(path_settings, d.errors()))

nfeatures = len(d['features'])

if nfeatures != 1:
    raise AssertionError(
        'Expecting 1 feature in {}, found {}'.format(path_settings, nfeatures))

region = d['features'][0]['geometry']
# assumes it is in decimal degrees, see Geod
cellsize = d['features'][0]['properties']['cellsize']
x, y = region['coordinates']
case_template = 'wflow_template'
filter_upstream_gt = 1000
crs = 'EPSG:4326'
SERVER_URL = 'http://hydro-engine.appspot.com'
case = 'wflow_demo'

g = Geod(ellps='WGS84')
# convert to meters in the center of the grid
# Earth Engine expects meters
_, _, crossdist_m = g.inv(x, y, x + cellsize, y + cellsize)
cellsize_m = sqrt(0.5 * crossdist_m ** 2)


def copycase(srccase, dstcase):
    """Set the case to a copy of the template case"""
    if os.path.isdir(dstcase):
        shutil.rmtree(dstcase)
    shutil.copytree(srccase, dstcase)


def post_data(url, data):
    """Post data using requests, throwing a clear error message if things go bad"""
    try:
        r = requests.post(url, json=data)
        r.raise_for_status()
    except requests.exceptions.HTTPError as err:
        print(err)
        print(r.text)
        sys.exit(1)
    return r


def get_data(url):
    """Get data using requests, throwing a clear error message if things go bad"""
    try:
        r = requests.get(url, stream=True)
        r.raise_for_status()
    except requests.exceptions.HTTPError as err:
        print(err)
        print(r.text)
        sys.exit(1)
    return r


def download_catchments(region, path):
    """Download a GeoJSON of the catchment upstream of `region`.
    Function copied from hydroengine, with added error reporting"""
    data = {'type': 'get_catchments', 'bounds': region, 'dissolve': True}
    r = post_data(SERVER_URL + '/get_catchments', data)
    with open(path, 'w') as f:
        f.write(r.text)


def download_rivers(region, path, filter_upstream_gt):
    """Download a GeoJSON of the rivers in the upstream catchment of `region`.
    Function copied from hydroengine, with added error reporting"""
    data = {'type': 'get_rivers', 'bounds': region}

    if filter_upstream_gt:
        data['filter_upstream_gt'] = filter_upstream_gt

    r = post_data(SERVER_URL + '/get_rivers', data)

    # download from url
    url = json.loads(r.text)['url']
    r = get_data(url)

    with open(path, 'wb') as f:
        r.raw.decode_content = True
        shutil.copyfileobj(r.raw, f)


def download_raster(region, path, variable, cell_size, crs):
    """Download a GeoTIFF raster of `variable` in the upstream catchment of `region`.
    Function copied from hydroengine, with added error reporting"""
    path_name = os.path.splitext(path)[0]

    data = {'type': 'get_raster', 'bounds': region,
            'variable': variable, 'cell_size': cell_size, 'crs': crs}

    r = post_data(SERVER_URL + '/get_raster', data)

    # download from url
    url = json.loads(r.text)['url']
    r = get_data(url)

    # download zip into a temporary file
    with tempfile.NamedTemporaryFile(delete=False) as f:
        r.raw.decode_content = True
        shutil.copyfileobj(r.raw, f)

    temp_dir = tempfile.mkdtemp()

    # unzip and rename both tfw and tif
    with zipfile.ZipFile(f.name, 'r') as zf:
        items = zf.namelist()
        zf.extractall(temp_dir)

    # move extracted files to the target path
    src_tfw = os.path.join(temp_dir, items[0])
    dst_tfw = path_name + '.tfw'
    if os.path.exists(dst_tfw):
        os.remove(dst_tfw)
    os.rename(src_tfw, dst_tfw)

    src_tif = os.path.join(temp_dir, items[1])
    if os.path.exists(path):
        os.remove(path)
    os.rename(src_tif, path)

    # clean-up
    os.rmdir(temp_dir)
    os.remove(f.name)


# start by making case an exact copy of the template
copycase(case_template, case)

# create grid
path_log = 'wtools_create_grid.log'
dir_mask = os.path.join(case, 'mask')
path_catchment = os.path.join(case, 'data/catchments/catchments.geojson')
projection = 'EPSG:4326'

download_catchments(region, path_catchment)

cg.main(path_log, dir_mask, path_catchment, projection,
        cellsize, snap=True, locationid=case)


# create static maps
dir_dest = os.path.join(case, 'staticmaps')
# use custom inifile, default high res ldd takes too long
path_inifile = os.path.join(case, 'data/staticmaps.ini')
path_dem_in = os.path.join(case, 'data/dem/dem.tif')
path_river = os.path.join(case, 'data/rivers/rivers.geojson')
path_catchment = os.path.join(case, 'data/catchments/catchments.geojson')
dir_lai = os.path.join(case, 'data/parameters/clim')

download_rivers(region, path_river, filter_upstream_gt)
download_raster(region, path_dem_in, 'dem', cellsize_m, crs)

other_maps = [
    'FirstZoneCapacity',
    'FirstZoneKsatVer',
    'FirstZoneMinCapacity',
    'InfiltCapSoil',
    'M',
    'PathFrac',
    'WaterFrac',
    'thetaS',
    'soil_type',
    'landuse'
]

# TODO rename these in hydro-engine
newnames = {
    'FirstZoneKsatVer': 'KsatVer',
    'FirstZoneMinCapacity': 'SoilMinThickness',
    'FirstZoneCapacity': 'SoilThickness',
    # 'landuse': 'wflow_landuse',  # restore after issues caused by values in this map are fixed
    'soil_type': 'wflow_soil'
}

path_other_maps = []
for param in other_maps:
    path = os.path.join(case, 'data/parameters', newnames.get(param, param) + '.tif')
    path_other_maps.append(path)


for param, path in zip(other_maps, path_other_maps):
    download_raster(region, path, param, cellsize_m, crs)

for m in range(1, 13):
    mm = str(m).zfill(2)
    path = os.path.join(dir_lai, 'LAI00000.0{}'.format(mm))
    download_raster(
        region, path, 'LAI{}'.format(mm), cellsize_m, crs)

sm.main(dir_mask, dir_dest, path_inifile, path_dem_in, path_river,
        path_catchment, lai=dir_lai, other_maps=path_other_maps)
