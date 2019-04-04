wflow
=====

wflow consists of a set of Python programs that can be run on the command line
and perform hydrological simulations. The models are based on the PCRaster
Python framework. In wflow this framework is extended (the wf_DynamicFramework)
so that models build using the framework can be controlled using the API.
Links to BMI, OpenMI and OpenDAP have been made.

A link to the latest version can always be found at https://github.com/openstreams/wflow

Reference documentation at:

+ http://wflow.readthedocs.io/en/latest/


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

+ Download pcraster from http://pcraster.geo.uu.nl/ website (version 4.2+)
+ Follow the installation instructions at http://pcraster.geo.uu.nl/quick-start-guide/


Install as a conda environment
------------------------------

The easiest and most robust way to install wflow is by installing it in a separate
conda environment. In the root repository directory there is an ``environment.yml`` file.
This file lists all dependencies, except PCRaster, which must be installed manually as
described above.

Run this command to start installing wflow with all dependencies:

+ ``conda env create -f environment.yml``

This creates a new environment with the name ``wflow``. To activate this environment in
a session, run:

+ ``activate wflow``

Now you should be able to start this environment's Python with ``python``, try
``import wflow`` to see if the package is installed.

More details on how to work with conda environments can be found here:
https://conda.io/docs/user-guide/tasks/manage-environments.html


Install using pip
-----------------

Besides the recommended conda environment setup described above, you can also install
wflow with ``pip``. For the more difficult to install Python dependencies, it is best to
use the conda package manager:

+ ``conda install numpy scipy gdal netcdf4 cftime pyproj python-dateutil``

This will install the latest release of wflow:

+ ``pip install wflow``

If you are planning to make changes and contribute to the development of wflow, it is
best to make a git clone of the repository, and do a editable install in the location
of you clone. This will not move a copy to your Python installation directory, but
instead create a link in your Python installation pointing to the folder you installed
it from, such that any changes you make there are directly reflected in your install.

+ ``git clone https://github.com/openstreams/wflow.git``
+ ``cd wflow``
+ ``pip install -e .``

Alternatively, if you want to avoid using ``git`` and simply want to test the latest
version from the ``master`` branch, you can replace the first line with downloading
a zip archive from GitHub: https://github.com/openstreams/wflow/archive/master.zip


Check if the installation is successful
---------------------------------------

To check it the install is successful, go to the examples directory and run the following command:

+ ``python -m wflow.wflow_sbm -C wflow_rhine_sbm -R testing``

This should run without errors.


Linux
-----

Although you can get everything with the python packages bundled with most linux distributions
(CentOS, Ubuntu, etc) we have found that the easiest way is to install the linux version of Anaconda
and use the conda tool to install all requirements apart from pcraster which has to be installed manually.

Since version 4.2, compiled versions of PCRaster are no longer distributed, so it will
need to be built following the instructions given at http://pcraster.geo.uu.nl/getting-started/pcraster-on-linux/


OSX
---
Unfortunately there is no pcraster build for osx yet. If anybody wants to pick this up please let
the guys at pcraster.eu know!

Running wflow_sbm from docker
=============================

To run the above from the docker container, get the docker image by either installing docker and running ``docker
buid .`` in your local repository directory. Alternatively, download the image from docker hub by typing

+ docker pull ewatercycle/wflow

After obtaining the docker image, run it by typing

+ docker run -v <path-to-data>:/data ewatercycle/wflow -R testing

This will create a (root-owned) subdirectory 'testing' in your data path with the model output.

Credits
=======

+ The stats.py script was made by Keith Cherkauer (https://engineering.purdue.edu/~cherkaue/software.htm)

+ pcraster is developed and maintained by Utrecht University (http://www.pcraster.eu)

+ netCDF4 is developed by unidata (http://unidata.github.io/netcdf4-python/)

+ GDAL is released under an X/MIT style Open Source license by the Open Source Geospatial Foundation (http://www.gdal.org).


Citation
========
See doi of the release you use. If you use a snapshot of the development (without a DOI) cite as follows:

Jaap Schellekens, Willem van Verseveld, Tanja Euser, Hessel Winsemius, Christophe Thiange, Laurène Bouaziz, Daniel Tollenaar, Sander de Vries, Albrecht Weerts, YEAR. openstreams/wflow: unstable-master. https://github.com/openstreams/wflow, obtained: DATE_OF_DOWNLOAD


Recent releases
---------------

.. image:: https://zenodo.org/badge/17738134.svg
   :target: https://zenodo.org/badge/latestdoi/17738134

Jaap Schellekens, Willem van Verseveld, Tanja Euser, Hessel Winsemius, Christophe Thiange, Laurène Bouaziz, Daniel Tollenaar, Sander de Vries, 2016. openstreams/wflow: 2016.04 Test release. doi:10.5281/zenodo.167057

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.155389.svg
   :target: https://doi.org/10.5281/zenodo.155389

Jaap Schellekens, Willem van Verseveld, Tanja Euser, Hessel Winsemius, Christophe Thiange, Laurène Bouaziz, Daniel Tollenaar, Sander de Vries, 2016. openstreams/wflow: 2016.03. doi:10.5281/zenodo.155389
