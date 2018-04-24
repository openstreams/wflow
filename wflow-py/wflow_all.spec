# -*- mode: python -*-

import os
import sys
import shutil
from distutils.dir_util import copy_tree, remove_tree
from pyproj import pyproj_datadir
from osgeo import gdal

# Set these for your installation

pcrasterlib = sys.argv[1]

# list identical make_wflow_exe script with --normal
# except for the wtools scripts
scriptpaths = [
    'openda_bmi/opendapy.py',
    'Scripts/bmi2runner.py',
    'Scripts/pcr2netcdf.py',
    'Scripts/wflow_prepare_step1.py',
    'Scripts/wflow_prepare_step2.py',
    'Scripts/wflow_sbm_rtc.py',
    'Scripts/wtools_py/modelbuilder.py',
    'wflow/create_grid.py',
    'wflow/static_maps.py',
    'wflow/wflow_adapt.py',
    'wflow/wflow_delwaq.py',
    'wflow/wflow_floodmap.py',
    'wflow/wflow_gr4.py',
    'wflow/wflow_hbv.py',
    'wflow/wflow_lintul.py',
    'wflow/wflow_pcrglobwb.py',
    'wflow/wflow_routing.py',
    'wflow/wflow_sbm.py',
    'wflow/wflow_sphy.py',
    'wflow/wflow_topoflex.py',
    'wflow/wflow_w3ra.py',
    'wflow/wflow_wave.py',
]


def scriptname(scriptpath):
    """Get 'wflow_hbv' from 'wflow/wflow_hbv.py'"""
    return os.path.splitext(os.path.basename(scriptpath))[0]


def do_analysis(scriptpath):
    """Run PyInstaller Analysis"""
    # note that the datas locations have to be set again in __init__.py
    # if they are to work in a bundled folder
    return Analysis([scriptpath],
                    binaries=[(pcrasterlib, '.')],
                    # TODO check if still necessary in PyInstaller 3.3 after
                    # https://github.com/pyinstaller/pyinstaller/pull/2401
                    # Though this seems more solid, submit as hook patch?
                    datas=[(gdal.GetConfigOption('GDAL_DATA'), 'gdal-data'),
                           (pyproj_datadir, 'proj-data')],
                    hiddenimports=['pywt._extensions._cwt',
                                   # in opendapy.py: importlib.import_module(sys.argv[3])
                                   # for wflow this would always be wflow.wflow_bmi
                                   'wflow.wflow_bmi',
                                   'rasterio.control',  # needed
                                   'rasterio.crs',  # needed
                                   'rasterio._shim',  # needed
                                   'rasterio.sample',  # needed
                                   'rasterio.vrt',  # needed
                                   'rasterio.coords',  # TODO test if needed
                                   'rasterio.enums',  # TODO test if needed
                                   'rasterio.env',  # TODO test if needed
                                   'rasterio.errors',  # TODO test if needed
                                   'rasterio.profiles',  # TODO test if needed
                                   'rasterio.transform',  # TODO test if needed
                                   'rasterio.vfs'  # TODO test if needed
                                   ])


def do_analysis_bare(scriptpath):
    """Run PyInstaller Analysis without extra binaries or datas"""
    # leave out the binaries and datas, only add for the first script
    # no need to copy for every script since they are later merged
    return Analysis([scriptpath])


def do_pyz(a):
    return PYZ(a.pure, a.zipped_data)


def do_exe(apyz):
    # unpack tuple created by zip
    a, pyz = apyz
    return EXE(pyz, a.scripts,
               exclude_binaries=True,
               name=scriptname(a.inputs[0]),
               upx=True,
               icon='logo.ico')


def do_collect(aexe):
    # unpack tuple created by zip
    a, exe = aexe
    return COLLECT(exe,
                   a.binaries,
                   a.zipfiles,
                   a.datas,
                   name=scriptname(a.inputs[0]),
                   upx=True)


if len(scriptpaths) == 1:
    analist = [do_analysis(scriptpaths[0])]
else:
    analist = [do_analysis(scriptpaths[0])] + \
        map(do_analysis_bare, scriptpaths[1:])

pyzlist = map(do_pyz, analist)
exelist = map(do_exe, zip(analist, pyzlist))
collist = map(do_collect, zip(analist, exelist))

# merge all individual folders in the parent 'dist' folder
# Replace with MERGE when this issue is fixed
# https://github.com/pyinstaller/pyinstaller/issues/1527
for scriptpath in scriptpaths:
    srcdir = os.path.join('dist', scriptname(scriptpath))
    # only copy if new or newer
    copy_tree(srcdir, 'dist', update=1)
    remove_tree(srcdir)
