Install PCRaster
================


Requirements
------------
The PCRaster Python package depends on Python and Numpy:

* Python: Version 2.7 is required (http://www.python.org/download/).
* Numpy: Version 1.6 or later is fine (https://sourceforge.net/projects/numpy/files/NumPy/).


Linux
-----
Installing PCRaster on Linux involves these steps:

* Unzip the zip file containing the software
* Update the ``PATH`` and ``PYTHONPATH`` environment variables.

PCRaster can be installed anywhere you want. Typical locations are ``$HOME``, ``/opt`` and ``/usr/local``.

.. code-block:: bash

   cd  /opt
   tar zxf /tmp/pcraster-lsbcc-4_x86-64.tar.gz

In order for the PCRaster executables and the Python package to be found, the ``PATH`` and ``PYTHONPATH`` environment variables must be updated with the paths to the executables and Python package, respectively. Assuming the use of the bash shell:

.. code-block:: bash

   export PATH=/opt/pcraster-lsbcc-4_86-64/bin:$PATH
   export PYTHONPATH=/opt/pcraster-lsbcc-4_86-64/python:$PYTHONPATH

These lines can be put in $HOME/.bash_profile to have them executed each time you login.

PCRaster is now installed and ready to be used.

.. note::

   In case the software doesn't work, verify that the Linux Standard Base (LSB) 4 package is installed.

PCRaster is known to work on the folowing distributions:

* Ubuntu 12.10
* Ubuntu 12.04
* bullx Linux, which is compatible with Red Hat Enterprise Linux

PCRaster is known to work on the Dutch national supercomputer, `Cartesius`_.

.. _Cartesius: https://www.surfsara.nl/systems/cartesius

Windows
-------
Installing PCRaster on Windows involves these steps:

* Unzip the zip file containing the software
* Update the ``PATH`` and ``PYTHONPATH`` environment variables.

PCRaster can be installed anywhere you want. Typical locations are ``%PROGRAMFILES%``, ``%PROGRAMFILES(X86)%`` and ``C:\``.

In order for the PCRaster executables and the Python package to be found, the ``PATH`` and ``PYTHONPATH`` environment variables must be updated with the paths to the executables and Python package, respectively.

On Windows XP environment variables can be changed like this (http://support.microsoft.com/kb/310519):

#. Right-click My Computer, and then click Properties.
#. Click the Advanced tab.
#. Click Environment variables.
#. Click one of the following options, for either a user or a system variable:

    * Click New to add a new variable name and value.
    * Click an existing variable, and then click Edit to change its value.

On other versions of Windows a similar procedure must be folowed.

An alternative is to create a batch script that can be run before using PCRaster:

.. code-block:: bat

   rem Configure environment for use of Python.
   set python_root=C:\Python27
   set PATH=%python_root%;%PATH%
   set python_root=

   rem Configure environment for use of PCRaster.
   rem This variable is different for each version of PCRaster.
   set pcraster_version=pcraster-4.0.0-beta-20130917_x86-32_msvs-9
   set pcraster_root=%HOMEPATH%\Desktop\%pcraster_version%
   set PATH=%pcraster_root%\bin;%PATH%
   set PYTHONPATH=%pcraster_root%\python;%PYTHONPATH%
   set pcraster_root=
   set pcraster_version=

PCRaster is now installed and ready to be used.

.. note::

   In case the software doesn't work, verify that the Microsoft Visual C++ 2005 Redistributable Package is installed:

     * `Redistributable Package for PCRaster 4, 32-bit version`_
     * `Redistributable Package for PCRaster 4, 64-bit version`_

.. _Redistributable Package for PCRaster 4, 32-bit version: http://www.microsoft.com/en-us/download/details.aspx?id=3387
.. _Redistributable Package for PCRaster 4, 64-bit version: http://www.microsoft.com/en-us/download/details.aspx?id=21254


Mac OS X
--------
TODO
