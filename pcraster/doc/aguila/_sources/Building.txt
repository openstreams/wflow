.. _building:

********
Building
********

.. note::

  This page is work in progress. Lots of stuff is missing.

Prerequisites
=============
Prerequisites (http://pcraster.sourceforge.net and http://sourceforge.net/projects/pcraster):

* DevEnv-1.3
* RasterFormat-1.2
* Xsd-1.2
* Dal-1.2

Make sure these projects are configured and built.

Make sure the folowing 3rd party libraries and tools are installed (see DevEnv's `configuration/Versions.sh` for verion information):

* Boost
* Geos
* Gdal
* Qwt
* Icu
* Python
* Qt
* Xerces
* Xsd
* CMake

Checkout/configure/build::

  $ svn checkout https://pcraster.svn.sourceforge.net/svnroot/pcraster/Aguila/branches/1.1.0 Aguila-1.1.0
  $ export AGUILA=<checkout path>/Aguila-1.1.0
  $ source $AGUILA/environment/configuration/bash_profile
  $ configurecmakeproject.py release Aguila
  $ make -C $AGUILA

You may want to checkout another version of the source code.
