

.. index::
   single: areanormal
.. _areanormal:

**********
areanormal
**********
.. topic:: areanormal

   Value assigned to an area taken from a normal distribution

::

  Result = areanormal(areaclass)

areaclass
   spatial
   boolean,nominal,ordinal

Result
   spatial
   scalar

Operation
=========


The area to which a cell belongs is identified by areaclass: cells with corresponding values on areaclass are member of a separate area. The Result is generated with a random number generator: for each area on areaclass, a random number is taken from a normal distribution with mean 0 and standard deviation 1. This value is assigned to all cells belonging to that area.  

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
   |    report Result = areanormal( Class);
   |   
   | • python
   |   Class = readmap("Class.map")
   |   Result = areanormal( Class)

   ============================================= ==========================================
   Result.map                                    Class.map                                 
   .. image::  ../examples/areanormal_Result.png .. image::  ../examples/areaarea_Class.png
   ============================================= ==========================================

   | 

