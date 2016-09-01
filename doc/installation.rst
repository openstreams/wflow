Installation
============

Introduction and targets
------------------------

The best way to install wflow is to use the setup.py install script (found in the wflow-py directory). Depending on
the way you want to use the models you can choose the install or develop targets
(e.g. python setup.py install). If you choose the develop target the installer
will make a link from the python library to the version you are working on
and ech chnage you make will end up in the library immediately. If you just want to
use the model it is best to choose the install target.


Dependencies
------------
in order to run wflow requires the following packages:

+ netCDF4
+ numpy
+ matplotlib
+ pcraster
+ osgeo
+ pyproj

The setup.py script will try to install these dependencies but it is best to make
sure you have installed and tested those before running the setup.py script.
Make sure to have 64 bit versions of all packages.

Windows
-------
We recommend installing Anaconda and install the required packages using the conda tool. We have
binary releases for windows available that allow you to run the models without installing python and all required packages:

+ https://github.com/openstreams/wflow/releases

if you want to adjust the code you should install python and all the required packages as per instructions below:

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
+ Extract zip to root of c: This will created c:\pcraster-4.1.0_x86-64
+ Add c:\pcraster-4.1.0_x86-64\python to the PYTHONPATH environment variable

*Installing wflow itself*

Clone with git or Download the latest zip with the source code of wflow. Go to the wflow-py directory and run:

+ python setup.py install

To check it the install is successfull, go to the  the the examples directory and run the following command:

+ Python c:\Anaconda\scripts\wflow_sbm.py -C wflow_rhine_hbv -T 100 -R testing

This should run without errors

Linux
-----

Although you can get everything with the python packages bundled with most linux distributions
(CentOS, Ubuntu, etc) we have found that the easiest way is to install the linux version of Anaconda
and use the conda tool to install all requirements apart from pcraster which has to be installed manually.


OSX
---
Unfortunately there is not pcraster build for osx yet. If anybody wants to pick this up please let
the guys at pcraster.eu know!
