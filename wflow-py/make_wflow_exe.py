from distutils.core import setup
from bbfreeze import Freezer
from _version import *
import ctypes,glob,os,shutil
import matplotlib

def dependencies_for_freeezing():
	import netCDF4_utils 

nrbits = str(ctypes.sizeof(ctypes.c_voidp) * 8)
includes = ['wflow.wflow_bmi','wflow.wflow_w3ra','wflow.wflow_bmi_combined','bmi','bmi.wrapper']

thename = "Wflow"+MVERSION+'-'+nrbits
f = Freezer("Wflow"+MVERSION+'-'+nrbits,includes=includes)
f.addScript("wflow/__init__.py")
f.addScript("wflow/wflow_topoflex.py")
f.addScript("wflow/wflow_sbm.py")
f.addScript("Scripts/wflow_sbm_rtc.py")
f.addScript("wflow/wflow_hbv.py")
f.addScript("wflow/wflow_adapt.py")
f.addScript("wflow/wflow_w3ra.py")
#f.addScript("wflow/wflow_hbv_snow2.py")
f.addScript("wflow/wflow_delwaq.py")
f.addScript("wflow/wflow_wave.py")
f.addScript("wflow/wflow_gr4.py")
f.addScript("wflow/wflow_floodmap.py")
f.addScript("wflow/wflow_routing.py")
f.addScript("Scripts/bmi2runner.py")
f.addScript("Scripts/wflow_prepare_step1.py")
#f.addScript("Scripts/area_in_out.py")
f.addScript("Scripts/wflow_prepare_step2.py")
f.addScript("Scripts/pcr2netcdf.py")
#f.addScript("wflow/wflow_fit.py") # Does not work becuse of QT
f()    # starts the freezing process


os.system('conda list' + ">" + os.path.join(thename,'packages.txt'))
# matplolib data files
data_files=matplotlib.get_py2exe_datafiles()

# pcraster dll's
ddir = "c:/pcraster4-64/lib/"
data_files.append((".", glob.glob(ddir + "/*.dll")))

# GDAL data files
gdaldata = os.getenv("GDAL_DATA")
data_files.append(("./gdal-data", glob.glob(gdaldata + "/*.*")))

print data_files
print "Copying extra data files..."
for dirr in data_files:
    timake = os.path.join(thename ,dirr[0])
    print timake
    if not os.path.exists(timake):
        os.makedirs(timake)
    for tocp in dirr[1]:
        print tocp
        shutil.copy(tocp,timake)
