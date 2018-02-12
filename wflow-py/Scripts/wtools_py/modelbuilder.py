# currently assumes the working directory is right above the wflow cases

import os
import sys
import geojson
import requests
import click
import wflow.create_grid as cg
import wflow.static_maps as sm
import wflow.wflowtools_lib as wt
import json
import shutil
import zipfile
import tempfile
from math import sqrt
from pyproj import Geod
import pcraster as pcr
from osgeo import gdalconst
import rasterio
from rasterio import warp

SERVER_URL = 'http://hydro-engine.appspot.com'


@click.command()
@click.option('--geojson-path',
              help='Path to a GeoJSON file with the geometry that needs to be path of the model.')
@click.option('--cellsize',
              default=0.01, show_default=True,
              help='Desired model cell size in decimal degrees.')
@click.option('--name',
              default='wflow_case', show_default=True,
              help='Name of the WFlow case.')
@click.option('--case-template',
              default='wflow_template', show_default=True,
              help='Name of the template WFlow case.')
@click.option('--case-path',
              default='.', show_default=True,
              help='Path where both the template and created case reside.')
@click.option('--fews/--no-fews',
              default=False, show_default=True,
              help='Flag indicating whether the WFlow case is part of a Delft-FEWS setup.')
@click.option('--fews-config-path',
              default='Config', show_default=True,
              help='Path to the Delft-FEWS config directory.')
@click.option('--dem-path',
              help='Optionally provide a local or improved Digital Elevation Model (DEM) '
              'to use instead of the default global DEM.')
def build_model(geojson_path, cellsize, name, case_template, case_path, fews, fews_config_path, dem_path):
    """Prepare a simple WFlow model, anywhere, based on global datasets."""

    # force all encodings to utf-8 directly, see http://click.pocoo.org/5/python3/
    geojson_path = geojson_path.encode('utf-8')
    name= name.encode('utf-8')
    case_template = case_template.encode('utf-8')
    case_path = case_path.encode('utf-8')
    fews_config_path = fews_config_path.encode('utf-8')
    if dem_path is not None:
        dem_path = dem_path.encode('utf-8')

    # assumes it is in decimal degrees, see Geod
    region = first_geometry(geojson_path)
    x, y = region['coordinates']
    filter_upstream_gt = 1000
    crs = 'EPSG:4326'
    case = os.path.join(case_path, name)

    g = Geod(ellps='WGS84')
    # convert to meters in the center of the grid
    # Earth Engine expects meters
    _, _, crossdist_m = g.inv(x, y, x + cellsize, y + cellsize)
    cellsize_m = sqrt(0.5 * crossdist_m ** 2)

    # start by making case an exact copy of the template
    copycase(case_template, case)

    # create grid
    path_log = 'wtools_create_grid.log'
    dir_mask = os.path.join(case, 'mask')
    path_catchment = os.path.join(case, 'data/catchments/catchments.geojson')
    projection = 'EPSG:4326'

    download_catchments(region, path_catchment)

    cg.main(path_log, dir_mask, path_catchment, projection,
            cellsize, locationid=name, snap=True)

    # create static maps
    dir_dest = os.path.join(case, 'staticmaps')
    # use custom inifile, default high res ldd takes too long
    path_inifile = os.path.join(case, 'data/staticmaps.ini')
    path_dem_in = os.path.join(case, 'data/dem/dem.tif')
    path_river = os.path.join(case, 'data/rivers/rivers.geojson')
    path_catchment = os.path.join(case, 'data/catchments/catchments.geojson')
    dir_lai = os.path.join(case, 'data/parameters/clim')

    download_rivers(region, path_river, filter_upstream_gt)
    if dem_path is None:
        # download the global dem
        download_raster(region, path_dem_in, 'dem', cellsize_m, crs)
    else:
        # warp the local dem onto model grid
        mask_tif = os.path.join(dir_mask, 'mask.tif')
        wt.warp_like(dem_path, path_dem_in, mask_tif,
                     format='GTiff', co={'dtype': 'float32'}, resampling=warp.Resampling.med)

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

    if fews:
        # save default state-files in FEWS-config
        dir_state = os.path.join(case, 'outstate')
        ensure_dir_exists(dir_state)
        state_files = ['CanopyStorage.map', 'GlacierStore.map', 'ReservoirVolume.map', 'SatWaterDepth.map', 'Snow.map',
                       'SnowWater.map', 'SurfaceRunoff.map', 'SurfaceRunoffDyn.map', 'TSoil.map',
                       'UStoreLayerDepth_0.map', 'WaterLevel.map', 'WaterLevelDyn.map']
        zip_name = name + '_GA_Historical default.zip'

        zip_loc = os.path.join(fews_config_path, 'ColdStateFiles', zip_name)
        path_csf = os.path.dirname(zip_loc)
        ensure_dir_exists(path_csf)

        mask = pcr.readmap(os.path.join(dir_mask, 'mask.map'))

        with zipfile.ZipFile(zip_loc, mode='w') as zf:
            for state_file in state_files:
                state_path = os.path.join(dir_state, state_file)
                pcr.report(pcr.cover(mask, pcr.scalar(0)), state_path)
                zf.write(state_path, state_file, compress_type=zipfile.ZIP_DEFLATED)


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


def ensure_dir_exists(path_dir):
    if not os.path.exists(path_dir):
        os.makedirs(path_dir)


def first_geometry(path_geojson):
    """Provided a path to a GeoJSON file,
    check if the GeoJSON is valid and has max 1 feature,
    then return its geometry."""

    with open(path_geojson) as f:
        d = geojson.load(f)

    if not d.is_valid:
        raise AssertionError(
            '{} is not a valid GeoJSON file\n{}'.format(path_geojson, d.errors()))

    if d.type == 'FeatureCollection':
        nfeatures = len(d.features)

        if nfeatures != 1:
            raise AssertionError(
                'Expecting 1 feature in {}, found {}'.format(path_geojson, nfeatures))

        geom = d.features[0].geometry
    elif d.type == 'Feature':
        geom = d.geometry
    else:
        geom = d

    if geom.type != 'Point':
        errmsg = 'Point is the only supported geometry type, found {} in {}'
        raise AssertionError(errmsg.format(geom.type, path_geojson))
    else:
        return geom


if __name__ == '__main__':
    # use sys.argv[1:] to allow using PyCharm debugger
    # https://github.com/pallets/click/issues/536
    build_model(sys.argv[1:])
