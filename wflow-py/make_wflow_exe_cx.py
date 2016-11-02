import sys
from cx_Freeze import setup, Executable, hooks
from _version import *
import ctypes,glob,os,shutil
import matplotlib
import scipy
includefiles_list=[]
scipy_path = os.path.dirname(scipy.__file__)
includefiles_list.append(scipy_path)


def load_scipy_patched(finder, module):
    """the scipy module loads items within itself in a way that causes
        problems without the entire package and a number of other subpackages
        being present."""
    finder.IncludePackage("scipy._lib")  # Changed include from scipy.lib to scipy._lib
    finder.IncludePackage("scipy.misc")

hooks.load_scipy = load_scipy_patched


nrbits = str(ctypes.sizeof(ctypes.c_voidp) * 8)
#includes = ['wflow.wflow_bmi','wflow.wflow_w3ra','wflow.wflow_bmi_combined','bmi','bmi.wrapper',"pcraster","osgeo.ogr"]

thename = "Wflow"+MVERSION+'-'+nrbits

packages = ["osgeo"]

options = {"packages": packages}
base=None

executables = [
    Executable('wflow/wflow_sbm.py', base=base),
    Executable('wflow/wflow_hbv.py', base=base)
]

setup(name='wflow',
      version='1.0',
      description='Wflow',
      options={"build_exe" : options},
      executables=executables,
      )









os.system('conda list' + ">" + os.path.join(thename,'packages.txt'))
# matplolib data files
data_files=matplotlib.get_py2exe_datafiles()

# pcraster dll's
ddir = "c:/pcraster/lib/"
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
