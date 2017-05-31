__all__ = ["wflow_funcs","wflow_adapt","wflow_lib","pcrut","wf_DynamicFramework","stats"]
__version__="1.0.master"
__release__="1.0.master.1"
__versionnr__="1.0.1"
__build__="2017-05-23 08:34:15.293000"


import osgeo.gdal as gdal

import os, sys
if hasattr(sys, "frozen"):
    print("wflow: frozen")
    _ROOT = os.path.abspath(os.path.dirname(__file__)).split("library.zip")[0]
    os.environ['GDAL_DATA'] = os.path.join(_ROOT,'gdal-data')
    _ROOT = os.path.abspath(os.path.dirname(__file__)).split("library.zip")[0]
    os.environ['PATH'] = os.environ['PATH'] + ";" + os.path.join(_ROOT)
else:
    print("wflow: python")
    _ROOT = os.path.abspath(os.path.dirname(__file__))

def get_data(path):
    return os.path.join(_ROOT, 'data', path)
