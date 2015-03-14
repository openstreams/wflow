Installation
============

Introduction and targets
------------------------

The best way to install wflow is to use the setup.py install script. Depending on
the way you want to use the models you can choose the instal or develop targets
(e.g. python setup.py install). If you choose the develope target the installer
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

The setup.py script will try to install these dependencies but it is best to make
sure you have installed and tested those before running the setup.py script.


Windows
-------


Linux
-----

Although you can get everything with the python packages bundled with most linux distributions
(CentOS, Ubuntu, etc) we have found that the easiest wasy is to install the linux version of Anaconda
and use the conda tool to install all requirements apart from pcraster which has to ben install by hand.


OSX
---

Unfortunately there is not pcraster build for osx yet. If anybody want to pick this up please let
the guys at pcraster.eu know!