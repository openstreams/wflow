from distutils.core import setup
from bbfreeze import Freezer
from _version import *
import ctypes
import glob, os, shutil
import matplotlib


def dependencies_for_freeezing():
    import netCDF4_utils
    import requests


data_files = matplotlib.get_py2exe_datafiles()
includes = ["wflow.wflow_bmi", "wflow.wflow_w3ra", "wflow.wflow_bmi_combined"]
# include_files = glob.glob("c:\Anaconda\Lib\site-packages\zmq\backend\cython\*.pyd")

nrbits = str(ctypes.sizeof(ctypes.c_voidp) * 8)

# Extra packages needed for the bmi linkage openda:
# https://github.com/openearth/mmi-python
# https://github.com/openearth/bmi-python
# conda install cython
# conda install pyzmq
# easy_install requests
thename = "Wflow" + MVERSION + "-" + nrbits + "-wflow_kernel_openda"
f = Freezer(thename, includes=includes)
f.addScript("wflow/__init__.py")
f.addScript("openda/thrift_bmi_raster_server.py")
f.addScript("wflow/wflow_sbm.py")
f.addScript("Scripts/wflow_sbm_rtc.py")
f.addScript("wflow/wflow_hbv.py")
f.addScript("wflow/wflow_adapt.py")
f.addScript("wflow/wflow_w3ra.py")
# f.addScript("wflow/wflow_hbv_snow2.py")
f.addScript("wflow/wflow_delwaq.py")
f.addScript("wflow/wflow_wave.py")
f.addScript("wflow/wflow_gr4.py")
f.addScript("wflow/wflow_floodmap.py")
f.addScript("wflow/wflow_routing.py")
f.addScript("wflow/wflow_sphy.py", base=base)
f.addScript("Scripts/bmi2runner.py")
f.addScript("Scripts/wflow_prepare_step1.py")
# f.addScript("Scripts/area_in_out.py")
f.addScript("Scripts/wflow_prepare_step2.py")
f.addScript("Scripts/pcr2netcdf.py")
# f.addScript("wflow/wflow_fit.py") # Does not work becuse of QT
f()  # starts the freezing process


os.system("conda list" + ">" + os.path.join(thename, "packages.txt"))
# Extra data directories

ddir = "c:/pcraster4-64/lib/"
data_files.append((".", glob.glob(ddir + "/*.dll")))

gdaldata = os.getenv("GDAL_DATA")
data_files.append(("./gdal-data", glob.glob(gdaldata + "/*.*")))


print "Copying extra data files..."
for dirr in data_files:
    timake = os.path.join(thename, dirr[0])
    print timake
    if not os.path.exists(timake):
        os.makedirs(timake)
    for tocp in dirr[1]:
        shutil.copy(tocp, timake)


# test with
#    thrift_bmi_raster_server.py wflow.wflow_bmi wflowbmi_csdms 127.0.0.1 49633
