============
Installation
============

The main dependencies for wflow are an installation of Python 3.6, and PCRaster 4.2+.
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
==============================

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
=================

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
=======================================

To check it the install is successful, go to the examples directory and run the following command:

+ ``python -m wflow.wflow_sbm -C wflow_rhine_sbm -R testing``

This should run without errors.

Linux
=====

Although you can get everything with the python packages bundled with most linux distributions
(CentOS, Ubuntu, etc) we have found that the easiest way is to install the linux version of Anaconda
and use the conda tool to install all requirements apart from pcraster which has to be installed manually.

Since version 4.2, compiled versions of PCRaster are no longer distributed, so it will
need to be built following the instructions given at http://pcraster.geo.uu.nl/getting-started/pcraster-on-linux/

OSX
===
Unfortunately there is no pcraster build for osx yet. If anybody wants to pick this up please let
the guys at pcraster.eu know!
