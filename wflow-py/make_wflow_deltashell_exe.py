from distutils.core import setup
from bbfreeze import Freezer
from _version import *
import ctypes
import glob,os,shutil
import matplotlib

def dependencies_for_freeezing():
    import netCDF4_utils
    import requests


data_files=matplotlib.get_py2exe_datafiles()
includes = ["zmq.backend.cython","requests","zmq.eventloop.zmqstream","pandas","matplotlib",'wflow.wflow_bmi']
#include_files = glob.glob("c:\Anaconda\Lib\site-packages\zmq\backend\cython\*.pyd")

nrbits = str(ctypes.sizeof(ctypes.c_voidp) * 8)

# Extra packages needed for the bmi linkage deltashell:
# https://github.com/openearth/mmi-python
# https://github.com/openearth/bmi-python
# conda install cython
# conda install pyzmq
# easy_install requests

f = Freezer("Wflow"+MVERSION+'-'+nrbits,includes=includes)
f.addScript("wflow/__init__.py")
f.addScript("wflow/wflow_sbm.py")
f.addScript("wflow/wflow_hbv.py")
f.addScript("wflow/wflow_adapt.py")
f.addScript("wflow/wflow_W3RA.py")
#f.addScript("wflow/wflow_hbv_snow2.py")
f.addScript("wflow/wflow_delwaq.py")
f.addScript("wflow/wflow_wave.py")
f.addScript("wflow/wflow_gr4.py")
f.addScript("wflow/wflow_floodmap.py")
f.addScript("wflow/wflow_routing.py")
#f.addScript("wflow/plottss.py")
f.addScript("Scripts/wfds_core.py")
f.addScript("Scripts/wflow_prepare_step1.py")
#f.addScript("Scripts/area_in_out.py")
f.addScript("Scripts/wflow_prepare_step2.py")
f.addScript("Scripts/pcr2netcdf.py")
#f.addScript("wflow/wflow_fit.py") # Does not work becuse of QT
f()    # starts the freezing process


ddir = "c:/pcraster4-64/lib/"
data_files.append((".", glob.glob(ddir + "/*.dll")))


shutil.copy("c:\Anaconda\Lib\site-packages\zmq\libzmq.pyd","Wflow"+MVERSION+'-'+nrbits +"/")

print "Copying extra data files..."
for dirr in data_files:
    timake = os.path.join("Wflow"+MVERSION+'-'+nrbits,dirr[0])
    print timake
    if not os.path.exists(timake):
        os.makedirs(timake)
    for tocp in dirr[1]:
        shutil.copy(tocp,timake)


