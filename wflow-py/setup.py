import os
from setuptools import setup
import versioneer

# since pcraster cannot be installed with pip, we check it like this
try:
    import pcraster
except:
    print(
        "Could not import pcraster, make sure it is installed including the python extensions"
    )
    print("see www.pcraster.eu")

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, "../README.rst")) as f:
    README = f.read()

setup(
    name='wflow',
    description='the wflow hydrological models (part of OpenStreams)',
    long_description=README,
    author='J. Schellekens',
    author_email='jaap.schellekens@deltares.nl',
    url='http://wflow.readthedocs.io/en/latest/',
    license="GPL",
    python_requires=">=3.6",
    install_requires=[
        'numpy',
        'gdal',
        'netCDF4',
        'pyproj',
    ],
    version=versioneer.get_version().split('+')[0],
    cmdclass=versioneer.get_cmdclass(),
    packages=['wflow', 'wflow.pcrglobwb',
              'wflow.sphy', 'wflow.wrappers.rtc'],
    package_dir={'wflow': 'wflow'},
    test_suite='UnitTests',
    zip_safe=False,
    scripts=['Scripts/pcr2netcdf.py', 'Scripts/tss2xml.py', 'Scripts/wflow_subcatch.py',
            'wflow/wflow_extract.py', 'wflow/wflow_sceleton.py',
            'wflow/wflow_gr4.py', 'wflow/plottss.py', 'wflow/wflow_wave.py', 'wflow/wflow_topoflex.py',
            'wflow/wflow_cqf.py', 'wflow/wflow_floodmap.py', 'wflow/wflow_upscale.py',
            'wflow/wflow_fit.py', 'wflow/wflow_adapt.py', 'wflow/wflow_delwaq.py',
            'Scripts/wflow_prepare_step1.py', 'Scripts/wflow_prepare_step2.py',
            'wflow/wflow_sbm.py', 'wflow/wflow_hbv.py', 'wflow/wflow_w3ra.py',
            'wflow/wflow_upscale.py', 'wflow/wflow_routing.py', 'Scripts/bmi2runner.py'],
)
