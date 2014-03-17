from distutils.core import setup
from bbfreeze import Freezer
import os


#vers=os.system('svnversion -n')
#print vers
#Binary dist
f = Freezer("Wflow1.0RC2")
f.addScript("wflow/__init__.py")
f.addScript("wflow/wflow_sbm.py")
f.addScript("wflow/wflow_hbv.py")
f.addScript("wflow/wflow_adapt.py")
f.addScript("wflow/wflow_delwaq.py")
f.addScript("wflow/wflow_wave.py")
f.addScript("wflow/wflow_gr4.py")
f.addScript("wflow/wflow_floodmap.py")
f.addScript("Scripts/wflow_prepare_step1.py")
f.addScript("Scripts/wflow_prepare_step2.py")
#f.addScript("wflow/wflow_fit.py") # Does not work becuse of QT
f()    # starts the freezing process
