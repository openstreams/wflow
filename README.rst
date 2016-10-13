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
There is absolutely no guarantee the executables will work on your computer.

Install as a python package (windows quick install instructions)
================================================================

*Installing Anaconda (scientific python distribution)*

Download Anaconda for python 2.7 64 bit (Tested with anaconda2 2.5.0). From the Anaconda installer choose the following options:
+ Install to c:\Anaconda
+ Register as default python
+ Add to path

Once Anaconda is installed open a command window and install netCDF4, gdal and pyproj using the following commands:

+ Conda install netCDF4
+ Conda install gdal=1.11
+ Conda install pyproj

*Installing pcraster*

+ Download pcraster from www.pcraster.eu website (version 4.1 64 bit)
+ Extract zip to root of c: This will create c:\\pcraster-4.1.0_x86-64
+ Add c:\\pcraster-4.1.0_x86-64\\python to the PYTHONPATH environment variable

*Installing wflow itself*

Clone with git or Download the latest zip with the source code of wflow. Go to the wflow-py directory and run:

+ python setup.py install

To check it the install is successful, go to the examples directory and run the following command:

+ Python c:\\Anaconda\\scripts\\wflow_sbm.py -C wflow_rhine_sbm -T 100 -R testing

This should run without errors

Credits
=======

+ The stats.py script was made by Keith Cherkauer (https://engineering.purdue.edu/~cherkaue/software.htm)

+ pcraster is developed and maintained by Utrecht University (http://www.pcraster.eu)

+ netCDF4 is developed by unidata (http://unidata.github.io/netcdf4-python/)


