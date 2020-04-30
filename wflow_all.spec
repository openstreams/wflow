# -*- mode: python -*-

import os
import shutil
import subprocess
import sys
from distutils.dir_util import copy_tree, remove_tree
from pathlib import Path

import pyproj
import xarray
from osgeo import gdal

# got a RecursionError: maximum recursion depth exceeded
sys.setrecursionlimit(10_000)

gdal.UseExceptions()
datas = [
    (gdal.GetConfigOption("GDAL_DATA"), "gdal-data"),
    (Path(xarray.__path__[0]) / "static", "xarray/static"),
]
pyproj_datadir = pyproj.datadir.get_data_dir()
# prevent unintentionally adding the entire workdir
if pyproj_datadir != "":
    datas.append((pyproj_datadir, "proj-data"))

print("data included in pyinstaller distribution:")
print(datas)

# list identical make_wflow_exe script with --normal
# except for the wtools scripts
scriptpaths = [
    "openda_bmi/opendapy.py",
    "Scripts/bmi2runner.py",
    "Scripts/pcr2netcdf.py",
    "Scripts/wflow_flood.py",
    "Scripts/wflow_prepare_step1.py",
    "Scripts/wflow_prepare_step2.py",
    "wflow/wflow_adapt.py",
    "wflow/wflow_delwaq.py",
    "wflow/wflow_emwaq.py",
    "wflow/wflow_sediment.py",
    "wflow/wflow_floodmap.py",
    "wflow/wflow_gr4.py",
    "wflow/wflow_hbv.py",
    "wflow/wflow_lintul.py",
    "wflow/wflow_pcrglobwb.py",
    "wflow/wflow_routing.py",
    "wflow/wflow_sbm.py",
    "wflow/wflow_sphy.py",
    "wflow/wflow_topoflex.py",
    "wflow/wflow_w3ra.py",
    "wflow/wflow_w3.py",
    "wflow/wflow_stream.py",
    "wflow/wflow_wave.py",
]


def scriptname(scriptpath):
    """Get 'wflow_hbv' from 'wflow/wflow_hbv.py'"""
    return os.path.splitext(os.path.basename(scriptpath))[0]


def do_analysis(scriptpath):
    """Run PyInstaller Analysis"""
    # note that the datas locations have to be set again in __init__.py
    # if they are to work in a bundled folder
    return Analysis(
        [scriptpath],
        # TODO check if still necessary in PyInstaller 3.3 after
        # https://github.com/pyinstaller/pyinstaller/pull/2401
        # Though this seems more solid, submit as hook patch?
        datas=datas,
        hiddenimports=[  # in opendapy.py: importlib.import_module(sys.argv[3])
            # for wflow this would always be wflow.wflow_bmi
            "wflow.wflow_bmi",
            "wflow.wflow_bmi_combined",
            "pkg_resources.py2_warn",
        ],
    )


def do_pyz(a):
    return PYZ(a.pure, a.zipped_data)


def do_exe(apyz):
    # unpack tuple created by zip
    a, pyz = apyz
    return EXE(
        pyz,
        a.scripts,
        exclude_binaries=True,
        name=scriptname(a.inputs[0]),
        upx=False,
        icon="logo.ico",
    )


def do_collect(aexe):
    # unpack tuple created by zip
    a, exe = aexe
    return COLLECT(
        exe, a.binaries, a.zipfiles, a.datas, name=scriptname(a.inputs[0]), upx=False
    )


analist = list(map(do_analysis, scriptpaths))
pyzlist = list(map(do_pyz, analist))
exelist = list(map(do_exe, zip(analist, pyzlist)))
collist = list(map(do_collect, zip(analist, exelist)))

filename = "version.txt"
versionfile = open(os.path.join("dist", filename), "w")
version = subprocess.check_output(["git", "describe", "--tags"]).strip().decode()
versionfile.write(version)
versionfile.close()


# merge all individual folders in the parent 'dist' folder
# Replace with MERGE when this issue is fixed
# https://github.com/pyinstaller/pyinstaller/issues/1527
for scriptpath in scriptpaths:
    srcdir = os.path.join("dist", scriptname(scriptpath))
    # only copy if new or newer
    copy_tree(srcdir, "dist", update=1)
    remove_tree(srcdir)
