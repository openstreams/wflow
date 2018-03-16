import os
import datetime

import subprocess

branch = None

try:
    branch = subprocess.check_output('git rev-parse HEAD', shell=True).strip()
except:
    branch = None

if branch == None:
    branch = 'master'

vers='1'
nrversion = '1.0'

###################################
manualversion = nrversion + "." + branch + "." + vers
manualmainversion = nrversion + "." + branch
version = nrversion = '1.0' + "." + vers
###################################
a = open("_version.py","w")

build=datetime.datetime.now()

a.write("# This file is made by mkversion_buildserver.py, do not edit!")
a.write("VERSION=\"" + manualversion +  "\"\n")
a.write("MVERSION=\"" + manualmainversion +"\"\n")
a.write("NVERSION=\"" + version +"\"\n")
a.write("BUILD=\"" + str(build) +"\"\n")
a.close()

a = open("wflow/__init__.py", "w")
a.write(
    "__all__ = ['wflow_funcs','wflow_adapt','wflow_lib','pcrut','wf_DynamicFramework','stats']\n")
a.write("__version__='" + manualmainversion + "'\n")
a.write("__release__='" + manualversion + "'\n")
a.write("__versionnr__='" + version + "'\n")
a.write("__build__='" + str(build) + "'\n\n")
a.write("import os, sys\n")
a.write("import osgeo.gdal as gdal\n\n")
a.write("if getattr(sys, 'frozen', False):\n")
a.write("    # running in a bundle\n")
a.write("    # sys._MEIPASS is set by PyInstaller\n")
a.write("    basedir = getattr(sys, '_MEIPASS', None)\n")
a.write("    # support also other bundlers\n")
a.write("    if not basedir:\n")
a.write("        basedir = os.path.dirname(sys.executable)\n\n")
a.write("    # use the included gdal-data\n")
a.write("    gdal_data_path = os.path.join(basedir, 'gdal-data')\n")
a.write("    gdal.SetConfigOption('GDAL_DATA', gdal_data_path)\n\n")
a.write("    # set environment variable instead of pyproj_datadir such\n")
a.write("    # that child processes will inherit it\n")
a.write("    os.environ['PROJ_DIR'] = os.path.join(basedir, 'proj-data')\n")
a.write("    os.environ['PATH'] = basedir + ';' + os.environ['PATH']\n")
a.write("    os.environ['PYTHONPATH'] = basedir + ';' + os.environ['PYTHONPATH']\n")


print "============================================================================="
print "Now install wflow using setup.py install and regenerate the documentation...."
print "============================================================================="
