"""
This script makes a stand-alone 'executable' of the wflow models. It is tested using
Anaconda on windows 64 bit and ubuntu xenial 64 bit

supported targets:
--normal
--openda - includes thrift connection to openda, Make sure you have thrift installed first
--deltashell - includes bmi/mmi link to deltashell. Windows only. Make sure you have zmq, bmi and mmi
  installed. bmi and mmi can be downloaded from the openearth github repository
"""


from cx_Freeze import setup, Executable, hooks

from _version import *
import ctypes,glob,os,shutil
import matplotlib
import scipy
import sys

target = 'normal'
# Filter out wflow specific options
if "--openda" in sys.argv:
    target = 'openda'
    sys.argv.remove("--openda")
if "--normal" in sys.argv:
    target = 'normal'
    sys.argv.remove("--normal")
if "--deltashell" in sys.argv:
    target = 'deltashell'
    sys.argv.remove("--deltashell")


pdir = os.path.dirname(sys.executable) + "/"

if sys.platform == 'win32':
    # list comes from: c:\Anaconda\conda-meta\mkl-11.3.3-1.json
    MKL_files= [pdir + "Library/bin/cilkrts20.dll",
        pdir + "Library/bin/ifdlg100.dll",
        pdir + "Library/bin/libchkp.dll",
        pdir + "Library/bin/libicaf.dll",
        pdir + "Library/bin/libifcoremd.dll",
        pdir + "Library/bin/libifcoremdd.dll",
        pdir + "Library/bin/libifcorert.dll",
        pdir + "Library/bin/libifcorertd.dll",
        pdir + "Library/bin/libifportmd.dll",
        pdir + "Library/bin/libimalloc.dll",
        pdir + "Library/bin/libiomp5md.dll",
        pdir + "Library/bin/libiompstubs5md.dll",
        pdir + "Library/bin/libmmd.dll",
        pdir + "Library/bin/libmmdd.dll",
        pdir + "Library/bin/libmpx.dll",
        pdir + "Library/bin/liboffload.dll",
        pdir + "Library/bin/mkl_avx.dll",
        pdir + "Library/bin/mkl_avx2.dll",
        pdir + "Library/bin/mkl_avx512.dll",
        pdir + "Library/bin/mkl_core.dll",
        pdir + "Library/bin/mkl_def.dll",
        pdir + "Library/bin/mkl_intel_thread.dll",
        pdir + "Library/bin/mkl_mc.dll",
        pdir + "Library/bin/mkl_mc3.dll",
        pdir + "Library/bin/mkl_msg.dll",
        pdir + "Library/bin/mkl_rt.dll",
        pdir + "Library/bin/mkl_sequential.dll",
        pdir + "Library/bin/mkl_tbb_thread.dll",
        pdir + "Library/bin/mkl_vml_avx.dll",
        pdir + "Library/bin/mkl_vml_avx2.dll",
        pdir + "Library/bin/mkl_vml_avx512.dll",
        pdir + "Library/bin/mkl_vml_cmpt.dll",
        pdir + "Library/bin/mkl_vml_def.dll",
        pdir + "Library/bin/mkl_vml_mc.dll",
        pdir + "Library/bin/mkl_vml_mc2.dll",
        pdir + "Library/bin/mkl_vml_mc3.dll",
        pdir + "Library/bin/svml_dispmd.dll"]


os.system("c:\Anaconda\python mkversion_buildserver.py")



data_files=[]
scipy_path = os.path.dirname(scipy.__file__)
data_files.append(scipy_path)


def load_scipy_patched(finder, module):
    """the scipy module loads items within itself in a way that causes
        problems without the entire package and a number of other subpackages
        being present."""
    finder.IncludePackage("scipy._lib")  # Changed include from scipy.lib to scipy._lib
    finder.IncludePackage("scipy.misc")

hooks.load_scipy = load_scipy_patched

def mkdatatuples(thelist,destdir="."):
    """
    input list of input files output lis list of tuples including destination
    :param list:
    :return:
    """
    ret = []
    for item in thelist:
        destfile = os.path.join(destdir,os.path.basename(item))
        ret.append((item,destfile))
    return ret

data_files.append('packages.txt')
os.system('c:\anaconda\scripts\conda list' + ">" + os.path.join('packages.txt'))
# matplolib data files


mpl =  matplotlib.get_py2exe_datafiles()

mplfiles = []
for mpldir in mpl:
    ddir = os.path.join('mpl-data',os.path.basename(mpldir[0]))
    data_files.extend(mkdatatuples(mpldir[1],destdir=ddir))

if sys.platform == 'win32':
    # MKL files
    data_files.extend(mkdatatuples(MKL_files,destdir="."))
    # pcraster dll's
    #ddir = "c:/pcraster/lib/"
    # for teamcity
    ddir = "D:/BuildAgent/work/wflow_exe/pcraster-4.1.0_x86-64"
    data_files.extend(mkdatatuples(glob.glob(ddir + "/*.dll"),destdir='.'))

# GDAL data files
gdaldata = os.getenv("GDAL_DATA")

if gdaldata == None:
    gdaldata = "c:\Anaconda\Library\share\gdal"

data_files.extend(mkdatatuples(glob.glob(gdaldata + "/*.*"),destdir='gdal-data'))


nrbits = str(ctypes.sizeof(ctypes.c_voidp) * 8)
#includes = ['wflow.wflow_bmi','wflow.wflow_w3ra','wflow.wflow_bmi_combined','bmi','bmi.wrapper',"pcraster","osgeo.ogr"]

thename = "wflow-bin/Wflow"+MVERSION+'-'+target+'-'+sys.platform+'-'+nrbits

packages = ["osgeo"]

if target == 'openda':
    import thrift.protocol.TBinaryProtocol as TBinaryProtocol
    import thrift.transport.THttpClient as THttpClient
    import thrift.protocol.TBinaryProtocol as TBinaryProtocol
    import thrift.transport.THttpClient as THttpClient
    includes = ['wflow.wflow_bmi','wflow.wflow_w3ra','wflow.wflow_bmi_combined','lxml.etree', 'lxml._elementpath', 'gzip']
    packages.append('openda_bmi')
elif target == 'deltashell':
    import zmq.libzmq
    data_files.extend([zmq.libzmq.__file__, ])
    includes = ["zmq.backend.cython","zmq.utils.garbage","requests","zmq.eventloop.zmqstream",
                 'wflow.wflow_bmi','wflow.wflow_w3ra','wflow.wflow_bmi_combined','lxml.etree', 'lxml._elementpath', 'gzip']
    packages.append('zmq.backend.cython')
    packages.append('bmi')
    #packages.append('pkg_resources')
else:
    includes = ['wflow.wflow_bmi', 'wflow.wflow_w3ra', 'wflow.wflow_bmi_combined','lxml.etree', 'lxml._elementpath', 'gzip', 'numpy.core._methods', 'numpy.lib.format']

#  "include_msvcr": True,
options = {"includes": includes, "packages": packages,'include_files': data_files, "build_exe": thename,
            'excludes': ['collections.abc']}
base=None




if target == 'openda':
    import thrift
    executables = [
        Executable('Scripts/wtools_py/CatchRiver.py', base=base),
        Executable('wflow/create_grid.py', base=base),
        Executable('wflow/static_maps.py', base=base),
        Executable('Scripts/wtools_py/wflow_fews.py', base=base)
        Executable('Scripts/pcr2netcdf.py', base=base),
        Executable('Scripts/bmi2runner.py', base=base),
        Executable('openda_bmi/opendapy.py', base=base),
        Executable('Scripts/wflow_prepare_step2.py', base=base),
        Executable('Scripts/wflow_prepare_step1.py', base=base),
        Executable('Scripts/wflow_sbm_rtc.py', base=base),
        Executable('wflow/wflow_topoflex.py', base=base),
        Executable('wflow/wflow_sbm.py', base=base),
        Executable('wflow/wflow_adapt.py', base=base),
        Executable('wflow/wflow_w3ra.py', base=base),
        Executable('wflow/wflow_delwaq.py', base=base),
        Executable('wflow/wflow_wave.py', base=base),
        Executable('wflow/wflow_gr4.py', base=base),
        Executable('wflow/wflow_floodmap.py', base=base),
        Executable('wflow/wflow_routing.py', base=base),
        Executable('wflow/wflow_hbv.py', base=base)
    ]
elif target == 'deltashell':
    executables = [
        Executable('Scripts/wtools_py/CatchRiver.py', base=base),
        Executable('wflow/create_grid.py', base=base),
        Executable('wflow/static_maps.py', base=base),
        Executable('Scripts/wtools_py/wflow_fews.py', base=base),
        Executable('Scripts/pcr2netcdf.py', base=base),
        Executable('Scripts/bmi2runner.py', base=base),
        Executable('Scripts/wfds_core.py', base=base),
        Executable('Scripts/wflow_prepare_step2.py', base=base),
        Executable('Scripts/wflow_prepare_step1.py', base=base),
        Executable('Scripts/wflow_sbm_rtc.py', base=base),
        Executable('wflow/wflow_topoflex.py', base=base),
        Executable('wflow/wflow_routing.py', base=base),
        Executable('wflow/wflow_sbm.py', base=base),
        Executable('wflow/wflow_adapt.py', base=base),
        Executable('wflow/wflow_w3ra.py', base=base),
        Executable('wflow/wflow_delwaq.py', base=base),
        Executable('wflow/wflow_wave.py', base=base),
        Executable('wflow/wflow_gr4.py', base=base),
        Executable('wflow/wflow_floodmap.py', base=base),
        Executable('wflow/wflow_hbv.py', base=base)
    ]
else:
    executables = [
        Executable('Scripts/wtools_py/CatchRiver.py', base=base),
        Executable('wflow/create_grid.py', base=base),
        Executable('wflow/static_maps.py', base=base),
        Executable('Scripts/wtools_py/wflow_fews.py', base=base),
        Executable('Scripts/pcr2netcdf.py', base=base),
        Executable('Scripts/bmi2runner.py', base=base),
        Executable('Scripts/wflow_prepare_step2.py', base=base),
        Executable('Scripts/wflow_prepare_step1.py', base=base),
        Executable('Scripts/wflow_sbm_rtc.py', base=base),
        Executable('wflow/wflow_topoflex.py', base=base),
        Executable('wflow/wflow_sbm.py', base=base),
        Executable('wflow/wflow_routing.py', base=base),
        Executable('wflow/wflow_adapt.py', base=base),
        Executable('wflow/wflow_w3ra.py', base=base),
        Executable('wflow/wflow_delwaq.py', base=base),
        Executable('wflow/wflow_wave.py', base=base),
        Executable('wflow/wflow_gr4.py', base=base),
        Executable('wflow/wflow_floodmap.py', base=base),
        Executable('wflow/wflow_hbv.py', base=base)
    ]

setup(name='wflow',
      version=NVERSION,
      description='Wflow',
      options={"build_exe" : options},
      executables=executables,
      )
