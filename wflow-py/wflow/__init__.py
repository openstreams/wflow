__all__ = ["wflow_funcs","wflow_adapt","wflow_lib","pcrut","wf_DynamicFramework","stats"]
__version__="1.0.master"
__release__="1.0.master.1"
__versionnr__="1.0.1"
__build__="2017-06-01 15:37:12.049000"
import os, sys
if hasattr(sys, "frozen"):
    print('Frozen...')
    _ROOT = os.path.abspath(os.path.dirname(__file__)).split("library.zip")[0]
    os.environ['GDAL_DATA'] = os.path.join(_ROOT,'gdal-data')
    os.environ['PATH'] = _ROOT + ';' + os.environ['PATH']
    os.environ['PYTHONPATH'] = _ROOT + ';' + os.environ['PYTHONPATH']
    sys.path.insert(0,_ROOT)
    if _ROOT not in os.environ['PATH']:
        print('Root dir of binary not insystem path. This may cause problems...')
else:
    _ROOT = os.path.abspath(os.path.dirname(__file__))

import osgeo.gdal as gdal

def get_data(path):
    return os.path.join(_ROOT, 'data', path)
