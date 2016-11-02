import os

import subprocess
branch = subprocess.check_output('git rev-parse --abbrev-ref HEAD', shell=True).strip()

vers='1'
nrversion = '1.0'

###################################
manualversion = nrversion + "." + branch + "." + vers
manualmainversion = nrversion + "." + branch
version = nrversion = '1.0' + "." + vers
###################################
a = open("_version.py","w")

a.write("VERSION=\"" + manualversion +  "\"\n")
a.write("MVERSION=\"" + manualmainversion +"\"\n")
a.write("NVERSION=\"" + version +"\"\n")

a.close()

a = open("wflow/__init__.py","w")
a.write("__all__ = [\"wflow_funcs\",\"wflow_adapt\",\"wflow_lib\",\"pcrut\",\"wf_DynamicFramework\",\"stats\"]\n")
a.write("__version__=\"" + manualmainversion + "\"\n")
a.write("__release__=\"" + manualversion + "\"\n")
a.write("__versionnr__=\"" + version + "\"\n")
a.write("import osgeo.gdal as gdal")


print "============================================================================="
print "Now install wflow using setup.py install and regenerate the documentation...."
print "============================================================================="
