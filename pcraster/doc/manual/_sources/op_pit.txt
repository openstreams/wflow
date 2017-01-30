

.. index::
   single: pit
.. _pit:

***
pit
***
.. topic:: pit

   Unique value for each pit cell

::

  Result = pit(ldd)

ldd
   spatial
   ldd

Result
   spatial
   nominal

Operation
=========


A pit is a cell whose neighbours all have a local drain direction in
direction of the pit cell. A pit cell doesn't have a local drain direction
because all its neighbours are at a higher elevation. Additionally the
outflow cell of each catchment on the map, which is a cell at the edge of
the map is also a pit. On a local drain direction network pits have a cell
value 5. 






For every pit cell on the local drain direction network ldd an unique number starting with 1 is assigned to the corresponding cell on Result; these are cells with a value 5 on ldd. The other cells on ldd do have a local drain direction and are assigned a 0 value on Result.     

Notes
=====


A missing value on ldd is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Neighbourhood operators; local drain directions 

See Also
========
:ref:`lddmask`

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Ldd = Ldd.map;
   |   initial
   |    report Result = pit(Ldd);
   |   
   | • python
   |   Ldd = readmap("Ldd.map")
   |   Result = pit(Ldd)

   ====================================== ====================================
   Result.map                             Ldd.map                             
   .. image::  ../examples/pit_Result.png .. image::  ../examples/accu_Ldd.png
   ====================================== ====================================

   | 

