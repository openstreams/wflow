__all__ = ["wflow_funcs","wflow_adapt","wflow_lib","pcrut","wf_DynamicFramework","stats"]
__version__="1.0.master"
__release__="1.0.master.1"

import sys
import os

# some check for frozen packages with bbfreeze

if hasattr(sys, "frozen"):
    _ROOT = os.path.abspath(os.path.dirname(__file__)).split("library.zip")[0]
    os.environ['GDAL_DATA'] = os.path.join(_ROOT,'gdal-data')
else:
    _ROOT = os.path.abspath(os.path.dirname(__file__))


import osgeo.gdal as gdal