WFLOW
=====

wflow consists of a set of python programs that can be run on the command line 
and perform hydrological simulations. The models are based on the PCRaster 
python framework. In wflow this framework is extended (the wf_DynamicFramework) 
so that models build using the framework can be controlled using the API. 
Links to BMI, OpenMI and OpenDAP have been made.

A link to the latest version can always be found at http://www.openstreams.nl or
https://github.com/jaapschellekens/wflow 

Reference documentation at:

+ http://wflow.readthedocs.org/


Obtaining wflow
===============

Goto https://github.com/jaapschellekens/wflow. There you can download the source or a release. Also make sure
you get the required third party models first (see below). The documentation can be found at
http://wflow.readthedocs.org

INSTALL
=======

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
+ osgeo

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
