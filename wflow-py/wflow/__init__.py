__all__ = ["wflow_funcs","wflow_adapt","wflow_lib","pcrut","wf_DynamicFramework","stats"]
__version__="1.0-fin"
__release__="1.0-RC8-fin-203-210"

import sys

if not 'sphinx' in sys.modules:
	import osgeo.gdal as gdal