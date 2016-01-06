from distutils.core import setup
from bbfreeze import Freezer
from _version import *
import ctypes

def dependencies_for_freeezing():
	import netCDF4_utils 

nrbits = str(ctypes.sizeof(ctypes.c_voidp) * 8)



f = Freezer("Wflow"+MVERSION+'-'+nrbits)
f.addScript("wflow/__init__.py")
f.addScript("wflow/wflow_sbm.py")
f.addScript("wflow/wflow_hbv.py")
f.addScript("wflow/wflow_adapt.py")
f.addScript("wflow/wflow_W3RA.py")
#f.addScript("wflow/wflow_hbv_snow2.py")
f.addScript("wflow/wflow_delwaq.py")
f.addScript("wflow/wflow_wave.py")
f.addScript("wflow/wflow_gr4.py")
f.addScript("wflow/wflow_floodmap.py")
f.addScript("wflow/wflow_routing.py")
#f.addScript("wflow/plottss.py")
f.addScript("Scripts/wfds_core.py")
f.addScript("Scripts/wflow_prepare_step1.py")
#f.addScript("Scripts/area_in_out.py")
f.addScript("Scripts/wflow_prepare_step2.py")
f.addScript("Scripts/pcr2netcdf.py")
#f.addScript("wflow/wflow_fit.py") # Does not work becuse of QT
f()    # starts the freezing process
