__all__ = [
    "wflow_funcs",
    "wflow_adapt",
    "wflow_lib",
    "pcrut",
    "wf_DynamicFramework",
    "stats",
]
import os
import sys

import osgeo.gdal as gdal

if getattr(sys, "frozen", False):
    # running in a bundle
    # sys._MEIPASS is set by PyInstaller
    basedir = getattr(sys, "_MEIPASS", None)
    # support also other bundlers
    if not basedir:
        basedir = os.path.dirname(sys.executable)

    # use the included gdal-data
    gdal_data_path = os.path.join(basedir, "gdal-data")
    gdal.SetConfigOption("GDAL_DATA", gdal_data_path)

    # set environment variable instead of pyproj_datadir such
    # that child processes will inherit it
    os.environ["PROJ_DIR"] = os.path.join(basedir, "proj-data")

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions
