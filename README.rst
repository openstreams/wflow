wflow
=====

wflow consists of a set of Python programs that can be run on the command line
and perform hydrological simulations. The models are based on the PCRaster
Python framework. In wflow this framework is extended (the wf_DynamicFramework)
so that models build using the framework can be controlled using the API.
Links to BMI, OpenMI and OpenDAP have been made.

A link to the latest version can always be found at https://github.com/openstreams/wflow

Reference documentation at:

+ https://wflow.readthedocs.io/en/latest/


Obtaining wflow
===============

Go to https://github.com/openstreams/wflow. There you can download the source or a release.
Also make sure you get the required third party models first (see below).


Which version to use
====================
The master branch can change rapidly (and break functionality without warning) so please
use one of the releases if possible. If you want to adjust things in the model(s) we
assume you should be comfortable using the master branch.


Installation
============

Install as a conda package
--------------------------

By far the easiest way to install wflow, is using the ``conda`` package manager. This
package manager comes with the Anaconda Python distribution.
``wflow`` is available in the `conda-forge <https://conda-forge.org/>`_ channel. To install
you can use the following command:

+ ``conda install -c conda-forge wflow``

If this works it will install wflow with all dependencies including Python and PCRaster,
and you skip the rest of the installation instructions.

Installing Python and PCRaster dependencies
-------------------------------------------

The main dependencies for wflow are an installation of Python 3.6+, and PCRaster 4.2+.
Only 64 bit OS/Python is supported.

*Installing Python*

For Python we recommend using the Anaconda Distribution for Python 3, which is available
for download from https://www.anaconda.com/download/. The installer gives the option to
add ``python`` to your ``PATH`` environment variable. We will assume in the instructions
below that it is available in the path, such that ``python``, ``pip``, and ``conda`` are
all available from the command line.

Note that there is no hard requirement specifically for Anaconda's Python, but often it
makes installation of required dependencies easier using the conda package manager.

*Installing pcraster*

+ If you are using conda, pcraster will be installed automatically in the section below, otherwise:
+ Download pcraster from http://pcraster.geo.uu.nl/ website (version 4.2+)
+ Follow the installation instructions at http://pcraster.geo.uu.nl/quick-start-guide/


Install as a conda environment
------------------------------

The easiest and most robust way to install wflow is by installing it in a separate
conda environment. In the root repository directory there is an ``environment.yml`` file.
This file lists all dependencies. Either use the ``environment.yml`` file from the master branch
(please note that the master branch can change rapidly and break functionality without warning),
or from one of the releases {release}.

Run this command to start installing all wflow dependencies:

+ ``conda env create -f environment.yml``

This creates a new environment with the name ``wflow``. To activate this environment in
a session, run:

+ ``activate wflow``

For the installation of wflow there are two options (from the Python Package Index (PyPI)
or from Github). To install a release of wflow from the PyPI (available from release 2018.1):

+ ``pip install wflow=={release}``

To install directly from GitHub (from the HEAD of the master branch):

+ ``pip install git+https://github.com/openstreams/wflow.git``

or from Github from a specific release:

+ ``pip install git+https://github.com/openstreams/wflow.git@{release}``

Now you should be able to start this environment's Python with ``python``, try
``import wflow`` to see if the package is installed.

More details on how to work with conda environments can be found here:
https://conda.io/docs/user-guide/tasks/manage-environments.html

If you are planning to make changes and contribute to the development of wflow, it is
best to make a git clone of the repository, and do a editable install in the location
of you clone. This will not move a copy to your Python installation directory, but
instead create a link in your Python installation pointing to the folder you installed
it from, such that any changes you make there are directly reflected in your install.

+ ``git clone https://github.com/openstreams/wflow.git``
+ ``cd wflow``
+ ``activate wflow``
+ ``pip install -e .``

Alternatively, if you want to avoid using ``git`` and simply want to test the latest
version from the ``master`` branch, you can replace the first line with downloading
a zip archive from GitHub: https://github.com/openstreams/wflow/archive/master.zip

Install using pip
-----------------

Besides the recommended conda environment setup described above, you can also install
wflow with ``pip``. For the more difficult to install Python dependencies, it is best to
use the conda package manager:

+ ``conda install numpy scipy gdal netcdf4 cftime xarray pyproj numba python-dateutil``

Then install a release {release} of wflow (available from release 2018.1) with pip:

+ ``pip install wflow=={release}``

If you want to avoid using ``conda``, an example of a PCRaster build and pip install on
Ubuntu Linux can be found in `issue #36 <https://github.com/openstreams/wflow/issues/36>`_.

Check if the installation is successful
---------------------------------------

To check it the install is successful, go to the examples directory and run the following command:

+ ``python -m wflow.wflow_sbm -C wflow_rhine_sbm -R testing``

This should run without errors.


Credits
=======

+ The stats.py script was made by Keith Cherkauer (https://engineering.purdue.edu/~cherkaue/software.htm)

+ pcraster is developed and maintained by Utrecht University (http://www.pcraster.eu)

+ netCDF4 is developed by unidata (http://unidata.github.io/netcdf4-python/)

+ GDAL is released under an X/MIT style Open Source license by the Open Source Geospatial Foundation (http://www.gdal.org).


Citation
========
See doi of the release you use. If you use a snapshot of the development (without a DOI) cite as follows:

Jaap Schellekens, Willem van Verseveld, Martijn Visser, Hessel Winsemius, Tanja Euser, Laurène Bouaziz, Christophe Thiange, Sander de Vries,
Hélène Boisgontier, Dirk Eilander, Daniel Tollenaar, Albrecht Weerts, Fedor Baart, Pieter Hazenberg, Arthur Lutz, Corine ten Velden,
Mischa Jansen, Imme Benedict, YEAR. openstreams/wflow: unstable-master. https://github.com/openstreams/wflow, obtained: DATE_OF_DOWNLOAD


Releases
--------

To check the doi of releases you use: https://doi.org/10.5281/zenodo.593510