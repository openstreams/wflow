Welcome to VariableCollection's documentation!
==============================================

This document contains informations about the Python implementation of the Arrayed variables concept in oldcalc.

Note: The collection module is under development, interfaces may change!

References to oldcalc refer to the documentation of the arrayed variables concept in oldcalc, accessible online_. 

.. _online: http://pcraster.geo.uu.nl/documentation/manual_updates/ArraysInCalc.html

Contents
--------

.. toctree::
   :maxdepth: 2


   tutorial/index
   examples/index
   reference/index



Current limitations
-------------------

- parameter files: for multidimensional arrays columns holding the index type names must be separated by a single space

.. - switching variables off is not supported


Release history/Changes
-----------------------

20090901

  Changes:

  - Index and Collection constructors, object name argument removed
  - Index arguments passed to collection must be given as list
  - initialise method removed, variables initialised from parameter file now by using the ValueFromParameterTable class

  Features:

  - support for multidimensional collections added
  - support for linking variable names to external names added
  - documentation improvements


20090811

- initial version




.. Indices and tables
   ==================

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`

