WFLOW
=====

wflow consists of a set of python programs that can be run on the command line 
and perform hydrological simulations. The models are based on the PCRaster 
python framework. In wflow this framework is extended (the wf_DynamicFramework) 
so that models build using the framework can be controlled using the API. 
Links to BMI, OpenMI and OpenDAP have been made.

A link to the latest version can always be found at http://www.openstreams.nl or
http://wflow.googlecode.com

INSTALL
=======

Assuming you have all supporting packages install installen a new wflow 
distribution entails running the setup.py script. This script follows
the general python setup.py syntax. As such install using

./setup.py install

should install the package as part of your local python installation.

requirements for windows:

- install python-2.7.?.msi, install in c:\python27
- install pcraster4.0 (see pcraster.eu web site)
  (we have found the Anaconda distribution to work very well in combination
   with pcraster)

Optional but highly recommended:
- install matplotlib for python 2.7
- install pyreadline
- install ipython
- install pyqt
- install spyder


Freezer
=======
The make_wflow_exe.py script builds a binary distribution of the models.
You need the bbfreeze package installed to do this yourself.



Credits
=======
the stats.py script was made by Keith Cherkauer
(https://engineering.purdue.edu/~cherkaue/software.htm)
