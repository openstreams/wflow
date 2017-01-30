

.. index::
   single: areauniform
.. _areauniform:

***********
areauniform
***********
.. topic:: areauniform

   Value assigned to area taken from an uniform distribution

::

  Result = areauniform(areaclass)

areaclass
   spatial
   boolean,nominal,ordinal

Result
   spatial
   scalar

Operation
=========


The area to which a cell belongs is identified by areaclass: cells with corresponding values on areaclass are member of a separate area. The Result is generated with a random number generator: for each area on areaclass, a random number between 0 and 1 is taken from a uniform distribution. This value is assigned to all cells belonging to that area.  

Notes
=====


A cell with a missing value on areaclass is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Random number generators; Areas 

See Also
========
:ref:`grouppointfield`, :ref:`groupmapfield`

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Class = Class.map;
   |   initial
   |    report Result = areauniform(Class);
   |   
   | • python
   |   Class = readmap("Class.map")
   |   Result = areauniform(Class)

   ============================================== ==========================================
   Result.map                                     Class.map                                 
   .. image::  ../examples/areauniform_Result.png .. image::  ../examples/areaarea_Class.png
   ============================================== ==========================================

   | 

