import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
from setuptools import find_packages

here = os.path.abspath(os.path.dirname(__file__))

README = open(os.path.join(here, "../README.rst")).read()
TODO = open(os.path.join(here, "TODO.txt")).read()


try:
    import osgeo
except:
    print("Could not import osgeo (gdal), make sure it is installed")

try:
    import netCDF4
except:
    print("Could not import netCDF4, make sure it is installed")

try:
    import pcraster
except:
    print(
        "Could not import pcraster, make sure it is installed including the python extensions"
    )
    print("see www.pcraster.eu")

try:
    import pyproj
except:
    print("Could not import pyproj, make sure it is installed")

try:
    import matplotlib
except:
    print("Could not import matplotlib, make sure it is installed")

try:
    import versioneer
except:
    print("Could not import versioneer, make sure it is installed")

# Source dist
setup(
    name="wflow",
    version=versioneer.get_version().split("+")[0],
    cmdclass=versioneer.get_cmdclass(),
    # version= MVERSION,
    packages=["wflow", "wflow.pcrglobwb", "wflow.sphy", "wflow.wrappers.rtc"],
    package_dir={"wflow": "wflow"},
    author="J. Schellekens",
    author_email="jaap.schellekens@deltares.nl",
    url="http://www.openstreams.nl",
    license="GPL",
    zip_safe=False,
    scripts=[
        "Scripts/pcr2netcdf.py",
        "Scripts/tss2xml.py",
        "Scripts/wflow_subcatch.py",
        "wflow/wflow_extract.py",
        "wflow/wflow_sceleton.py",
        "wflow/wflow_gr4.py",
        "wflow/plottss.py",
        "wflow/wflow_wave.py",
        "wflow/wflow_topoflex.py",
        "wflow/wflow_cqf.py",
        "wflow/wflow_floodmap.py",
        "wflow/wflow_upscale.py",
        "wflow/wflow_fit.py",
        "wflow/wflow_adapt.py",
        "wflow/wflow_delwaq.py",
        "Scripts/wflow_prepare_step1.py",
        "Scripts/wflow_prepare_step2.py",
        "wflow/wflow_sbm.py",
        "wflow/wflow_hbv.py",
        "wflow/wflow_w3ra.py",
        "wflow/wflow_upscale.py",
        "wflow/wflow_routing.py",
        "Scripts/bmi2runner.py",
    ],
    description="the wflow hydrological models (part of OpenStreams)",
)
