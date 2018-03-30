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

import ctypes,glob,os,shutil
import matplotlib
import scipy
import sys
import glob
import json
import versioneer
import subprocess

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
    mkl_file = glob.glob(pdir + "conda-meta" + "/" + "mkl-[!service]*.json")[0]
    data = json.load(open(mkl_file))


MKL_files = [(pdir + s) for s in data['files']]


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
    # pcraster dll's for teamcity
    ddir = "d:/BuildAgent/work/wflow_exe/pcraster-4.1.0_x86-64/lib"
    data_files.extend(mkdatatuples(glob.glob(ddir + "/*.dll"),destdir='.'))

# GDAL data files
gdaldata = os.getenv("GDAL_DATA")

if gdaldata == None:
    gdaldata = "c:\Anaconda\Library\share\gdal"

data_files.extend(mkdatatuples(glob.glob(gdaldata + "/*.*"),destdir='gdal-data'))


nrbits = str(ctypes.sizeof(ctypes.c_voidp) * 8)
#includes = ['wflow.wflow_bmi','wflow.wflow_w3ra','wflow.wflow_bmi_combined','bmi','bmi.wrapper',"pcraster","osgeo.ogr"]

#clean wflow-bin dir
bf = os.path.join(os.getcwd(),"wflow-bin")
shutil.rmtree(bf,ignore_errors=True)

#a = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], shell=True).strip()

versions = versioneer.get_versions()
MVERSION = versions['version'].split('+')[0]

thename = "wflow-bin/Wflow-"+MVERSION+'-'+target+'-'+sys.platform+'-'+nrbits

data_files.append('\build\lib\wflow\_version.py')

packages = ["osgeo"]

if target == 'openda':
    import thrift.protocol.TBinaryProtocol as TBinaryProtocol
    import thrift.transport.THttpClient as THttpClient
    import thrift.protocol.TBinaryProtocol as TBinaryProtocol
    import thrift.transport.THttpClient as THttpClient
    includes = ['wflow.wflow_bmi','wflow.wflow_w3ra','wflow.wflow_bmi_combined','lxml.etree', 'lxml._elementpath', 'gzip', 'numpy.core._methods', 'numpy.lib.format']
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
        Executable('wflow/wflow_hbv.py', base=base),
        Executable('wflow/wflow_sphy.py', base=base),
        Executable('wflow/wflow_pcrglobwb.py', base=base)
    ]
elif target == 'deltashell':
    executables = [
        Executable('Scripts/wtools_py/CatchRiver.py', base=base),
        Executable('wflow/create_grid.py', base=base),
        Executable('wflow/static_maps.py', base=base),
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
        Executable('wflow/wflow_hbv.py', base=base),
        Executable('wflow/wflow_sphy.py', base=base),
        Executable('wflow/wflow_pcrglobwb.py', base=base)
    ]
else:
    executables = [
        Executable('Scripts/wtools_py/CatchRiver.py', base=base),
        Executable('wflow/create_grid.py', base=base),
        Executable('wflow/static_maps.py', base=base),
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
        Executable('wflow/wflow_hbv.py', base=base),
        Executable('wflow/wflow_sphy.py', base=base),
        Executable('wflow/wflow_pcrglobwb.py', base=base)
    ]


setup(name='wflow',
      version=MVERSION,
      description='Wflow',
      options={"build_exe" : options},
      executables=executables,
      )
