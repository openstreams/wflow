WFLOW
=====

wflow consists of a set of python programs that can be run on the command line 
and perform hydrological simulations. The models are based on the PCRaster 
python framework. In wflow this framework is extended (the wf_DynamicFramework) 
so that models build using the framework can be controlled using the API. 
Links to BMI, OpenMI and OpenDAP have been made.

A link to the latest version can always be found at http://www.openstreams.nl or
https://github.com/openstreams/wflow 

Reference documentation at:

+ http://wflow.readthedocs.org/


Obtaining wflow
===============

Goto https://github.com/openstreams/wflow. There you can download the source or a release. Also make sure
you get the required third party models first (see below). The documentation can be found at
http://wflow.readthedocs.org


Binaries for windows
====================
For windows binaries check the releases (https://github.com/openstreams/wflow/releases). These can be used 
if you do not have a python installation. However, it is recommended to install wflow as a python package (see below).

Install as a python package
===========================

Assuming you have all supporting packages install installing a new wflow
distribution entails running the setup.py script. This script follows
the general python setup.py syntax. As such running:

./setup.py install

should install the package as part of your local python installation.


in order to run wflow requires the following packages:

+ netCDF4
+ numpy
+ matplotlib
+ pcraster
+ osgeo (GDAL=1.11)
+ pyproj

The setup.py script will try to install these dependencies but it is best to make
sure you have installed and tested those before running the setup.py script.
Make sure to have 64 bit versions of all packages.

Freezer
=======
The make_wflow_exe.py script builds a binary distribution of the models.
You need the bbfreeze package installed to do this yourself.



Credits
=======

+ The stats.py script was made by Keith Cherkauer (https://engineering.purdue.edu/~cherkaue/software.htm)

+ pcraster is developed and maintained by Utrecht University (http://www.pcraster.eu)

+ netCDF4 is developed by unidata (http://unidata.github.io/netcdf4-python/)


