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
import pkg_resources

# import pyproj and gdal, to be able to set paths to shared resources
# in case we are using a frozen bundle created with pyinstaller
import pyproj
import pyproj.datadir
import osgeo.gdal as gdal

gdal.UseExceptions()

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
    os.environ["GDAL_DATA"] = gdal_data_path

    # set environment variable next to pyproj_datadir such
    # that child processes will inherit it
    proj_data_dir = os.path.join(basedir, "proj-data")
    pyproj.datadir.set_data_dir(proj_data_dir)
    os.environ["PROJ_LIB"] = proj_data_dir

try:
    from wflow.version import version

    __version__ = version
except:
    pass
