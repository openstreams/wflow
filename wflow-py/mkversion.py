import os
import datetime

import subprocess

branch = None

try:
    branch = subprocess.check_output('git rev-parse --abbrev-ref HEAD', shell=True).strip()
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

a.write("VERSION=\"" + manualversion +  "\"\n")
a.write("MVERSION=\"" + manualmainversion +"\"\n")
a.write("NVERSION=\"" + version +"\"\n")
a.write("BUILD=\"" + str(build) +"\"\n")
a.close()

a = open("wflow/__init__.py","w")
a.write("__all__ = [\"wflow_funcs\",\"wflow_adapt\",\"wflow_lib\",\"pcrut\",\"wf_DynamicFramework\",\"stats\"]\n")
a.write("__version__=\"" + manualmainversion + "\"\n")
a.write("__release__=\"" + manualversion + "\"\n")
a.write("__versionnr__=\"" + version + "\"\n")
a.write("__build__=\"" + str(build) + "\"\n")
a.write("import osgeo.gdal as gdal\n\n")
a.write("import os, sys\n")
a.write("if hasattr(sys, \"frozen\"):\n")
a.write("    _ROOT = os.path.abspath(os.path.dirname(__file__)).split(\"library.zip\")[0]\n")
a.write("    os.environ['GDAL_DATA'] = os.path.join(_ROOT,'gdal-data')\n")
a.write("else:\n")
a.write("    _ROOT = os.path.abspath(os.path.dirname(__file__))\n")

a.write("def get_data(path):\n")
a.write("    return os.path.join(_ROOT, 'data', path)\n")


print "============================================================================="
print "Now install wflow using setup.py install and regenerate the documentation...."
print "============================================================================="
