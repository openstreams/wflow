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
with open(os.path.join(here, "README.rst"), encoding="utf-8") as f:
    README = f.read()

setup(
    name="wflow",
    description="wflow hydrological modeling framework",
    long_description=README,
    author="J. Schellekens",
    author_email="wflow@deltares.nl",
    url="http://wflow.readthedocs.io/",
    license="GPL",
    python_requires=">=3.6",
    install_requires=[
        "numpy",
        "scipy",
        "gdal",
        "netCDF4",
        "cftime",
        "pyproj",
        "python-dateutil",
    ],
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=["wflow", "wflow.pcrglobwb", "wflow.sphy", "wflow.wrappers.rtc"],
    package_dir={"wflow": "wflow"},
    test_suite="tests",
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
    classifiers=[
        # https://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
    ],
    keywords="wflow hydrology modeling framework pcraster",
)
