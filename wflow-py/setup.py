import os
from _version import *

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
from setuptools import find_packages
here = os.path.abspath(os.path.dirname(__file__))

README = open(os.path.join(here, '../README.rst')).read()
TODO = open(os.path.join(here, 'TODO.txt')).read()


try:
    import pcraster
except:
    print("Could not import pcraster, make sure it is installed including the python extensions")
    print("see www.pcraster.eu")

#setup(**config)

# Source dist
setup(name='wflow',
      version= MVERSION,
      packages=['wflow'],
      package_dir={'wflow': 'wflow'},
      author='J. Schellekens',
      author_email='jaap.schellekens@deltares.nl',
      url='http://www.openstreams.nl',
      license = "GPL",
      scripts=['Scripts/pcr2netcdf.py','Scripts/tss2xml.py',
               'wflow/wflow_extract.py','wflow/wflow_sceleton.py',
               'wflow/wflow_gr4.py','wflow/plottss.py','wflow/wflow_wave.py',
               'wflow/wflow_cqf.py','wflow/wflow_floodmap.py','wflow/wflow_upscale.py',
               'wflow/wflow_fit.py','wflow/wflow_adapt.py','wflow/wflow_delwaq.py',
               'Scripts/wflow_prepare_step1.py','Scripts/wflow_prepare_step2.py',
               'wflow/wflow_sbm.py','wflow/wflow_hbv.py','wflow/wflow_W3RA.py',
               'wflow/wflow_upscale.py','wflow/wflow_routing.py'],
      description='the wflow hydrological models (part of OpenStreams)',
      )

