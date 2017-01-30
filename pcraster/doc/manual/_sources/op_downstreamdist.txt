

.. index::
   single: downstreamdist
.. _downstreamdist:

**************
downstreamdist
**************
.. topic:: downstreamdist

   Distance to the first cell downstream

::

  Result = downstreamdist(ldd)

ldd
   spatial
   ldd

Result
   spatial
   scalar

Options
=======
:literal:`--unittrue` or :literal:`--unitcell`

:literal:`--unittrue`
   distance is measured in true distance (default)

:literal:`--unitcell`
   distance is measured in number of cell lengths



Operation
=========


For each cell, assigns to Result the distance to the first cell downstream, where downstream cells are determined using the local drain directions on ldd. This distance is the length of one cell in case the local drain direction is to one of the right, left, top or bottom neighbouring cells or sqrt(2) multiplied by the length of one cell in case the local drain direction is to one of the 4 neighbouring cells in diagonal directions. In case a cell doesn't have a downstream cell (i.e. a pit) a zero is assigned to Result.  

Notes
=====


A cell with a missing value on ldd is assigned a missing value on Result.  

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
   |    Ldd2 = Ldd2.map;
   |   initial
   |    report Result = downstreamdist(Ldd2);
   |   
   | • python
   |   Ldd2 = readmap("Ldd2.map")
   |   Result = downstreamdist(Ldd2)

   ================================================= =====================================
   Result.map                                        Ldd2.map                             
   .. image::  ../examples/downstreamdist_Result.png .. image::  ../examples/accu_Ldd2.png
   ================================================= =====================================

   | 

