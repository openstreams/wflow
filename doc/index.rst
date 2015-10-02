===========================================
Welcome to wflow's |version| documentation!
===========================================

.. note::

      This documentation is for version |version| of wflow, release |release|
      This documentation was generated |today|


.. note::

    wflow is released under version 3 of the GPL

    wflow uses pcraster/python (see http://www.pcraster.eu) as it's calculation engine.

    In order to tun the python code you must install the following python packages:

    + pcraster
    + netcdf4
    + numpy
    + matplotlib
    + gdal

    Alternatively you can use one of the releases (executable) to run the models withour installing python and the
    required packages:   https://github.com/openstreams/wflow/releases


Introduction
============

This document describes the wflow distributed hydrological modelling platform.
wflow is part of the Deltares'
OpenStreams project (http://www.openstreams.nl). Wflow consists of a
set of python programs that can be run on the command line and perform
hydrological simulations. The models are based on the PCRaster python
framework. In wflow this framework is extended (the ``wf_DynamicFramework``
class) so that models build using the framework can be controlled using
the API. Links to OpenMI (www.openmi.org), BMI and OpenDA (www.openda.org)
are being tested.

The  wflow distributed hydrological model platform currently includes
the following models:

-  the wflow\_sbm  model (derived from `topog\_sbm <http://www.per.clw.csiro.au/topog/>`_ )

-  the wflow\_hbv model (a distributed version of the HBV96 model).

-  the wflow\_gr4 model (a distributed version of the gr4h/d models).

-  the wflow\_W3RA model.

-  the wflow\_routing model (a kinematic wave model that can run on the output of tone of the hydrological models).

-  the wflow\_wave model (a dynamic wave model that can run on the output of the wflow\_routing model).

-  the wflow\_floodmap model (a flood mapping model that can use the output of the wflow\_wave model or de wflow\_routing model).

The low level api and links to other frameworks allow the models to be
linked as part of larger modelling systems:


.. digraph:: Linking

    WFLOW_HBV -> WFLOWAPI;
    WFLOW_SBM -> WFLOWAPI;
    WFLOWAPI -> "PI"  [dir=both];
    "Data and Models" -> "PI";
    WFLOWAPI -> OpenMI  [dir=both];
    ModelX -> OpenMI;
    ModelY -> OpenMI;
    ModelY -> BMI;
    WFLOWAPI -> OpenDA  [dir=both];
    Calibration -> OpenDA;
    Assimilation -> OpenDA;
    WFLOWAPI [shape=square];
    OpenDA [shape=square];
    OpenMI [shape=square];
    BMI [shape=square];
    "PI" [shape=square];
    dpi=69;
    
    
.. note::

    wflow is part of the Deltares OpenStreams project
    (http://www.openstreams.nl). The OpenStreams project is a work in
    progress. Wflow functions as a toolkit for distributed hydrological
    models within OpenStreams.


The different wflow models share the same structure but are fairly
different with respect to the conceptualisation. The shared software
framework includes the basic maps (dem, landuse, soil etc) and the
hydrological routing via the kinematic wave. The Python class framework
also exposes the models as an API and is based on the PCRaster/Python
version 4.0 Beta (www.pcraster.eu).

The wflow\_sbm model maximises the use of available spatial data.
Soil depth, for example, is estimated from the DEM using a topographic
wetness index . The model is derived from the [CQFLOW]_ model that has
been applied in various countries, most notably in Central America. The
wflow\_hbv model is derived from the HBV-96 model but does not
include the routing functions, instead it uses the same kinematic wave
routine as the wflow\_sbm  model to route the water downstream.

The models are programmed python using the pcraster python extension.
As such, the structure of the model is
transparent, can be changed by other modellers easily, and the system
allows for rapid development. In order to run
the model both PCRaster 4.* and Python 2.7 are needed. At the moment
only 64 bit versions are supported.



Installation
============
.. toctree::
   :maxdepth: 2

   installation


How to use the models
=====================
.. toctree::
   :maxdepth: 2

   wflow_usage


Building a model
================
.. toctree::
   :maxdepth: 2

   wflow_building


FAQ
===
.. toctree::
   :maxdepth: 2

   faq

The wflow\_hbv model
====================
.. toctree::
   :maxdepth: 2

   wflow_hbv

The wflow\_sbm model
====================
.. toctree::
   :maxdepth: 2

   wflow_sbm

The wflow\_gr4 model
====================
.. toctree::
   :maxdepth: 2

   wflow_gr4   

The wflow\_W3RA model
=====================
.. toctree::
   :maxdepth: 2

   wflow_W3RA

The wflow\_Topoflex model
=========================
.. toctree::
   :maxdepth: 2

   wflow_topoflex

The wflow\_routing model
========================
.. toctree::
   :maxdepth: 2

   wflow_routing

The wflow\_wave model
=====================
.. toctree::
   :maxdepth: 2

   wflow_wave

The wflow\_floodmap model
=========================
.. toctree::
   :maxdepth: 2

   wflow_floodmap


The wflow Delft-FEWS adapter
============================
.. toctree::
   :maxdepth: 2

   wflow_adapt


wflow modules and libraries
===========================
.. toctree::
   :maxdepth: 2

   wflow_fit
   wflow_lib
   wflow_delwaq

   
Examples and tests
==================
.. toctree::
   :maxdepth: 2
   
   testrunner_wflowhbv
   calib_report


BMI Basic modeling interface
============================
.. toctree::
   :maxdepth: 2

   wflow_bmi

Adding a new model using the framework
======================================
.. toctree::
   :maxdepth: 2

   framework
   wf_DynamicFramework

   
Release notes
=============
.. toctree::
   :maxdepth: 2

   release_notes


OpenDA
======
.. toctree::
   :maxdepth: 2

   wflow_openda





Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


References
==========

.. [CQFLOW] Köhler, L., Mulligan, M., Schellekens, J., Schmid, S. and
    Tobón, C.: Final Technical Report DFID-FRP Project no. R7991 Hydrological
    impacts of converting tropical montane cloud forest to pasture, with
    initial reference to northern Costa Rica.,, 2006.


TODO
====

.. todolist::
