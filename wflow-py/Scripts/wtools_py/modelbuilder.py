# currently assumes the working directory is right above the wflow cases

import json
import os
import shutil
import sys
import tempfile
import zipfile
from math import sqrt

import click
import fiona
import geojson
import numpy as np
import rasterio
import requests
from shapely.ops import unary_union
import shapely.geometry as sg
from pyproj import Geod
from rasterio import warp

import pcraster as pcr
import wflow.create_grid as cg
import wflow.static_maps as sm
import wflow.wflowtools_lib as wt
from wflow import ogr2ogr

SERVER_URL = "http://hydro-engine.appspot.com"


@click.command()
@click.option(
    "--geojson-path",
    help="Path to a GeoJSON file with the geometry that needs to be part of the model.",
)
@click.option(
    "--cellsize",
    default=0.01,
    show_default=True,
    help="Desired model cell size in decimal degrees.",
)
@click.option(
    "--model",
    type=click.Choice(["sbm", "hbv", "w3ra"]),
    default="sbm",
    show_default=True,
    help="Name of the WFlow model concept.",
)
@click.option(
    "--timestep",
    type=click.Choice(["hourly", "daily"]),
    default="daily",
    show_default=True,
    help="Model time step - daily or hourly (only for hbv).",
)
@click.option("--name", help="Name of the WFlow case.")
@click.option("--case-template", help="Name of the template WFlow case.")
@click.option(
    "--case-path",
    default=".",
    show_default=True,
    help="Path where both the template and created case reside.",
)
@click.option(
    "--fews/--no-fews",
    default=False,
    show_default=True,
    help="Flag indicating whether the WFlow case is part of a Delft-FEWS setup.",
)
@click.option(
    "--fews-config-path",
    default="Config",
    show_default=True,
    help="Path to the Delft-FEWS config directory.",
)
@click.option(
    "--dem-path",
    help="Optionally provide a local or improved Digital Elevation Model (DEM) "
    "to use instead of the default global DEM.",
)
@click.option(
    "--river-path",
    help="Optionally provide a local or improved river vector file "
    "to use instead of the default global one.",
)
@click.option(
    "--outlet-path",
    help="Optionally provide a point layer with outlets "
    "to use instead of the automatically generated one.",
)
@click.option(
    "--region-filter",
    type=click.Choice(["catchments-upstream", "catchments-intersection", "region"]),
    default="catchments-upstream",
    show_default=True,
    help="Tell hydro-engine which model area to pick, by default this "
    "is everything upstream of the provided geometry, but it is also "
    "possible to get only the current catchment (catchments-intersection), "
    "or just exactly the provided geometry (region), like your own "
    "catchment polygon.",
)
def build_model(
    geojson_path,
    cellsize,
    model,
    timestep,
    name,
    case_template,
    case_path,
    fews,
    fews_config_path,
    dem_path,
    river_path,
    outlet_path,
    region_filter,
):
    """Prepare a simple WFlow model, anywhere, based on global datasets."""

    # lists below need to stay synchronized, not sure of a better way
    [
        geojson_path,
        model,
        timestep,
        name,
        case_template,
        case_path,
        fews_config_path,
        dem_path,
        river_path,
        outlet_path,
        region_filter,
    ] = [
        encode_utf8(p)
        for p in [
            geojson_path,
            model,
            timestep,
            name,
            case_template,
            case_path,
            fews_config_path,
            dem_path,
            river_path,
            outlet_path,
            region_filter,
        ]
    ]

    # fill in the dependent defaults
    if name is None:
        name = "wflow_{}_case".format(model)
    if case_template is None:
        case_template = "wflow_{}_template".format(model)
    if model == "hbv":
        if timestep == "hourly":
            case_template = "wflow_{}_hourly_template".format(model)
        else:
            case_template = "wflow_{}_daily_template".format(model)

    # assumes it is in decimal degrees, see Geod
    case = os.path.join(case_path, name)
    path_catchment = os.path.join(case, "data/catchments/catchments.geojson")

    region = hydro_engine_geometry(geojson_path, region_filter)

    # get the centroid of the region, such that we have a point for unit conversion
    centroid = sg.shape(region).centroid
    x, y = centroid.x, centroid.y

    filter_upstream_gt = 1000
    crs = "EPSG:4326"

    g = Geod(ellps="WGS84")
    # convert to meters in the center of the grid
    # Earth Engine expects meters
    _, _, crossdist_m = g.inv(x, y, x + cellsize, y + cellsize)
    cellsize_m = sqrt(0.5 * crossdist_m ** 2)

    # start by making case an exact copy of the template
    copycase(case_template, case)

    # create folder structure for data folder
    for d in ["catchments", "dem", "rivers"]:
        dir_data = os.path.join(case, "data", d)
        ensure_dir_exists(dir_data)

    # create grid
    path_log = "wtools_create_grid.log"
    dir_mask = os.path.join(case, "mask")
    projection = "EPSG:4326"

    download_catchments(
        region, path_catchment, geojson_path, region_filter=region_filter
    )
    cg_extent = path_catchment

    cg.main(
        path_log, dir_mask, cg_extent, projection, cellsize, locationid=name, snap=True
    )
    mask_tif = os.path.join(dir_mask, "mask.tif")

    with rasterio.open(mask_tif) as ds:
        bbox = ds.bounds

    # create static maps
    dir_dest = os.path.join(case, "staticmaps")
    # use custom inifile, default high res ldd takes too long
    path_inifile = os.path.join(case, "data/staticmaps.ini")
    path_dem_in = os.path.join(case, "data/dem/dem.tif")
    dir_lai = os.path.join(case, "data/parameters/clim")

    if river_path is None:
        # download the global dataset
        river_data_path = os.path.join(case, "data/rivers/rivers.geojson")
        # raise ValueError("User must supply river_path for now, see hydro-engine#14")
        download_rivers(
            region, river_data_path, filter_upstream_gt, region_filter=region_filter
        )
    else:
        # take the local dataset, reproject and clip
        # command line equivalent of
        # ogr2ogr -t_srs EPSG:4326 -f GPKG -overwrite -clipdst xmin ymin xmax ymax rivers.gpkg rivers.shp
        river_data_path = os.path.join(case, "data/rivers/rivers.gpkg")
        ogr2ogr.main(
            [
                "",
                "-t_srs",
                "EPSG:4326",
                "-f",
                "GPKG",
                "-overwrite",
                "-clipdst",
                str(bbox.left),
                str(bbox.bottom),
                str(bbox.right),
                str(bbox.top),
                river_data_path,
                river_path,
            ]
        )

    if dem_path is None:
        # download the global dem
        download_raster(
            region, path_dem_in, "dem", cellsize_m, crs, region_filter=region_filter
        )
    else:
        # warp the local dem onto model grid
        wt.warp_like(
            dem_path,
            path_dem_in,
            mask_tif,
            format="GTiff",
            co={"dtype": "float32"},
            resampling=warp.Resampling.med,
        )

    other_maps = {
        "sbm": [
            "FirstZoneCapacity",
            "FirstZoneKsatVer",
            "FirstZoneMinCapacity",
            "InfiltCapSoil",
            "M",
            "PathFrac",
            "WaterFrac",
            "thetaS",
            "soil_type",
            "landuse",
        ],
        "hbv": [
            "BetaSeepage",
            "Cfmax",
            "CFR",
            "FC",
            "K0",
            "LP",
            "Pcorr",
            "PERC",
            "SFCF",
            "TT",
            "WHC",
        ],
    }

    # TODO rename these in hydro-engine
    newnames = {
        "FirstZoneKsatVer": "KsatVer",
        "FirstZoneMinCapacity": "SoilMinThickness",
        "FirstZoneCapacity": "SoilThickness",
        "landuse": "wflow_landuse",
        "soil_type": "wflow_soil",
    }

    # destination paths
    path_other_maps = []
    for param in other_maps[model]:
        path = os.path.join(
            case, "data/parameters", newnames.get(param, param) + ".tif"
        )
        path_other_maps.append(path)

    for param, path in zip(other_maps[model], path_other_maps):
        if model == "sbm":
            download_raster(
                region, path, param, cellsize_m, crs, region_filter=region_filter
            )
        elif model == "hbv":
            # these are not yet in the earth engine, use local paths
            if timestep == "hourly":
                path_staticmaps_global = (
                    r"p:\1209286-earth2observe\HBV-GLOBAL\staticmaps_hourly"
                )
            else:
                path_staticmaps_global = (
                    r"p:\1209286-earth2observe\HBV-GLOBAL\staticmaps"
                )
            path_in = os.path.join(path_staticmaps_global, param + ".tif")

            # warp the local staticmaps onto model grid
            wt.warp_like(
                path_in,
                path,
                mask_tif,
                format="GTiff",
                co={"dtype": "float32"},
                resampling=warp.Resampling.med,
            )

    if model == "sbm":
        ensure_dir_exists(dir_lai)
        for m in range(1, 13):
            mm = str(m).zfill(2)
            path = os.path.join(dir_lai, "LAI00000.0{}".format(mm))
            download_raster(
                region,
                path,
                "LAI{}".format(mm),
                cellsize_m,
                crs,
                region_filter=region_filter,
            )
    else:
        # TODO this creates defaults in static_maps, disable this behavior?
        # or otherwise adapt static_maps for the other models
        dir_lai = None

    # create default folder structure for running wflow
    dir_inmaps = os.path.join(case, "inmaps")
    ensure_dir_exists(dir_inmaps)
    dir_instate = os.path.join(case, "instate")
    ensure_dir_exists(dir_instate)
    for d in ["instate", "intbl", "intss", "outmaps", "outstate", "outsum", "runinfo"]:
        dir_run = os.path.join(case, "run_default", d)
        ensure_dir_exists(dir_run)

    if outlet_path is None:
        # this is for coastal catchments only, if it is not coastal and no outlets
        # are found, then it will just be the pit of the ldd
        outlets = outlets_coords(path_catchment, river_data_path)
    else:
        # take the local dataset, reproject and clip
        outlet_data_path = os.path.join(case, "data/rivers/outlets.gpkg")
        ogr2ogr.main(
            [
                "",
                "-t_srs",
                "EPSG:4326",
                "-f",
                "GPKG",
                "-overwrite",
                "-clipdst",
                str(bbox.left),
                str(bbox.bottom),
                str(bbox.right),
                str(bbox.top),
                outlet_data_path,
                outlet_path,
            ]
        )
        x = []
        y = []
        with fiona.open(outlet_data_path) as c:
            for f in c:
                coords = f["geometry"]["coordinates"]
                x.append(coords[0])
                y.append(coords[1])
            outlets_x = np.array(x)
            outlets_y = np.array(y)
        outlets = outlets_x, outlets_y

    sm.main(
        dir_mask,
        dir_dest,
        path_inifile,
        path_dem_in,
        river_data_path,
        path_catchment,
        lai=dir_lai,
        other_maps=path_other_maps,
        outlets=outlets,
    )

    if fews:
        # save default state-files in FEWS-config
        dir_state = os.path.join(case, "outstate")
        ensure_dir_exists(dir_state)
        if model == "sbm":
            state_files = [
                "CanopyStorage.map",
                "GlacierStore.map",
                "ReservoirVolume.map",
                "SatWaterDepth.map",
                "Snow.map",
                "SnowWater.map",
                "SurfaceRunoff.map",
                "SurfaceRunoffDyn.map",
                "TSoil.map",
                "UStoreLayerDepth_0.map",
                "WaterLevel.map",
                "WaterLevelDyn.map",
            ]
        elif model == "hbv":
            state_files = [
                "DrySnow.map",
                "FreeWater.map",
                "InterceptionStorage.map",
                "LowerZoneStorage.map",
                "SoilMoisture.map",
                "SurfaceRunoff.map",
                "UpperZoneStorage.map",
                "WaterLevel.map",
            ]
        zip_name = name + "_GA_Historical default.zip"

        zip_loc = os.path.join(fews_config_path, "ColdStateFiles", zip_name)
        path_csf = os.path.dirname(zip_loc)
        ensure_dir_exists(path_csf)

        mask = pcr.readmap(os.path.join(dir_mask, "mask.map"))

        with zipfile.ZipFile(zip_loc, mode="w") as zf:
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


def download_catchments(
    region, path, path_geojson, region_filter="catchments-upstream", catchment_level=6
):
    """Download a GeoJSON of the catchment upstream of `region`.
    Function copied from hydroengine, with added error reporting"""

    if region_filter == "region":
        # use the geojson as catchment
        # some checks are already done in hydro_engine_geometry
        shutil.copy2(path_geojson, path)
    else:
        data = {
            "type": "get_catchments",
            "region": region,
            "dissolve": True,
            "region_filter": region_filter,
            "catchment_level": catchment_level,
        }
        r = post_data(SERVER_URL + "/get_catchments", data)
        with open(path, "w") as f:
            f.write(r.text)


def download_rivers(
    region,
    path,
    filter_upstream_gt,
    region_filter="catchments-upstream",
    catchment_level=6,
):
    """Download a GeoJSON of the rivers in the upstream catchment of `region`.
    Function copied from hydroengine, with added error reporting"""
    data = {
        "type": "get_rivers",
        "region": region,
        "region_filter": region_filter,
        "catchment_level": catchment_level,
    }

    if filter_upstream_gt:
        data["filter_upstream_gt"] = filter_upstream_gt

    r = post_data(SERVER_URL + "/get_rivers", data)

    with open(path, "w") as f:
        f.write(r.text)


def download_raster(
    region,
    path,
    variable,
    cell_size,
    crs,
    region_filter="catchments-upstream",
    catchment_level=6,
):
    """Download a GeoTIFF raster of `variable` in the upstream catchment of `region`.
    Function copied from hydroengine, with added error reporting"""
    path_name = os.path.splitext(path)[0]
    # print region.type
    with open("db.json", "w") as f:
        f.write(geojson.dumps(region))
    # print region
    data = {
        "type": "get_raster",
        "region": region,
        "variable": variable,
        "cell_size": cell_size,
        "crs": crs,
        "region_filter": region_filter,
        "catchment_level": catchment_level,
    }

    r = post_data(SERVER_URL + "/get_raster", data)

    # download from url
    url = json.loads(r.text)["url"]
    r = get_data(url)

    # download zip into a temporary file
    with tempfile.NamedTemporaryFile(delete=False) as f:
        r.raw.decode_content = True
        shutil.copyfileobj(r.raw, f)

    temp_dir = tempfile.mkdtemp()

    # unzip and rename both tfw and tif
    with zipfile.ZipFile(f.name, "r") as zf:
        items = zf.namelist()
        zf.extractall(temp_dir)

    # move extracted files to the target path
    src_tfw = os.path.join(temp_dir, items[0])
    dst_tfw = path_name + ".tfw"
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


def hydro_engine_geometry(path_geojson, region_filter):
    """Provided a path to a GeoJSON file, check if it is valid,
    then return its geometry. Hydro-engine currently only accepts
    geometries, so we need to enforce this here."""

    with open(path_geojson) as f:
        d = geojson.load(f)

    if not d.is_valid:
        raise AssertionError(
            "{} is not a valid GeoJSON file\n{}".format(path_geojson, d.errors())
        )

    # this needs special casing since we want to be able to specify a
    # FeatureCollection of catchment polygons, though at the same time make it
    # compatible with the hydro-engine 1 geometry requirement
    polytypes = ("Polygon", "MultiPolygon")
    if region_filter == "region":
        # these polygon checks should cover all types of GeoJSON
        if d.type == "FeatureCollection":
            # this should check all features
            gtype = d.features[0].geometry.type
            if gtype not in polytypes:
                raise ValueError(
                    "Geometry type in {} is {}, needs to be a polygon".format(
                        path_geojson, gtype
                    )
                )
            # combine features into 1 geometry
            polys = []
            for fcatch in d.features:
                g = sg.shape(fcatch["geometry"])
                polys.append(g)
            # now simplify it to a rectangular polygon
            bnds = unary_union(polys).bounds
            geom = sg.mapping(sg.box(*bnds))
        elif d.type == "Feature":
            assert d.geometry.type in polytypes
            geom = d.geometry
        else:
            assert d.type in polytypes
            geom = d
        # now simplify it to a rectangular polygon
        bnds = sg.shape(geom).bounds
        geom = sg.mapping(sg.box(*bnds))
    else:
        if d.type == "FeatureCollection":
            nfeatures = len(d.features)
            if nfeatures != 1:
                raise AssertionError(
                    "Expecting 1 feature in {}, found {}".format(
                        path_geojson, nfeatures
                    )
                )
            geom = d.features[0].geometry
        elif d.type == "Feature":
            geom = d.geometry
        else:
            geom = d

    return geom


def encode_utf8(path):
    """Modify path to encode in utf-8"""
    # see http://click.pocoo.org/5/python3/
    if path is None:
        return None
    else:
        return path.encode("utf-8")


def outlets_coords(path_catchment, river_data_path):
    """Get an array of X and Y coordinates of the outlets."""

    outlets = find_outlets(path_catchment, river_data_path, max_dist=0.02)
    outlets = sg.mapping(outlets)["coordinates"]

    outlets_x = np.array([c[0] for c in outlets])
    outlets_y = np.array([c[1] for c in outlets])
    return outlets_x, outlets_y


def ne_coastal_zone(scale="medium", buffer=0.1):
    """Get a MultiPolygon of the coastal zone from Natural Earth"""

    scales = {"large": "10m", "medium": "50m", "small": "110m"}
    res = scales[scale]
    coastline_url = (
        "https://raw.githubusercontent.com/nvkelso/"
        "natural-earth-vector/master/geojson/ne_{}_coastline.geojson".format(res)
    )
    r = requests.get(coastline_url)
    r.raise_for_status()
    js = r.json()

    coasts = []

    for f in js["features"]:
        linestring = sg.shape(f["geometry"])
        polygon = linestring.buffer(buffer)
        coasts.append(polygon)

    # add an extra buffer of 0, to create a valid geometry
    return sg.MultiPolygon(coasts).buffer(0.0)


def find_outlets(catch_path, riv_path, max_dist=0.02):
    """Find outlets by seeing how close rivers end to
    coastal catchment boundaries."""
    coastal_catch = catchment_coast_boundaries(catch_path)
    rivnodes = river_ends(riv_path)

    outlets = []
    # if it is empty, all distances are 0
    if not coastal_catch.is_empty:
        for p in rivnodes:
            dist = p.distance(coastal_catch)
            if dist <= max_dist:
                outlets.append(p)
    return sg.MultiPoint(outlets)


def catchment_boundaries(catch_path):
    """Collect all catchment boundaries into a MultiLineString"""
    geoms = []
    with fiona.open(catch_path) as c:
        for f in c:
            g = sg.shape(f["geometry"]).boundary
            geoms.append(g)

    return unary_union(geoms)


def catchment_coast_boundaries(catch_path, coastal_zone=None):
    """Collect all coastal catchment boundaries into a MultiLineString"""
    catchment_boundary = catchment_boundaries(catch_path)
    if coastal_zone is None:
        coastal_zone = ne_coastal_zone()
    return catchment_boundary.intersection(coastal_zone)


def river_ends(riv_path):
    """Collect the ends of river LineStrings into a MultiPoint"""
    # the hydrosheds rivers are split into many small segments
    # so to get the true ends, we should eliminate the ends
    # that connect to other ends
    mp0 = []  # line start nodes
    mp1 = []  # line end nodes
    with fiona.open(riv_path) as c:
        for f in c:
            coords = f["geometry"]["coordinates"]
            mp0.append(sg.Point(coords[0]))
            mp1.append(sg.Point(coords[-1]))

    mp0 = sg.MultiPoint(mp0)
    mp1 = sg.MultiPoint(mp1)
    # get rid of the connecting inner nodes
    return mp0.symmetric_difference(mp1)


if __name__ == "__main__":
    # use sys.argv[1:] to allow using PyCharm debugger
    # https://github.com/pallets/click/issues/536
    build_model(sys.argv[1:])
