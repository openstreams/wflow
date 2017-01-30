.. index::
   single: lddrepair
.. _lddrepair:

*********
lddrepair
*********
.. topic:: lddrepair

   Reparation of unsound local drain direction map

::

  Result = lddrepair(ldd)

ldd
   spatial
   ldd

Result
   spatial
   ldd

Operation
=========


Each cell on a local drain direction map must have a pit at the end of its
downstream path. If this is not the case for one or more cells on a local
drain direction map, the map is called unsound. An unsound local drain
direction map can not be used as input expression for the operations with
local drain direction maps. 

The lddrepair operation changes the cell values on the unsound ldd in such a way that it becomes sound: all downstream paths will end in a pit cell; this adjusted ldd is saved as Result.   

The repair operation is done as follows. Two things may be the cause of
a downstream path not ending in a pit cell: a set of cells in a cycle and
cells draining to a missing value or to the outside of the map. A cycle is
a set of cells that don't drain to a pit because they drain to each other, in
a closed cycle. The smallest cycle consists of two cells with local drain
directions to each other; larger cycles may consist of several cells. First,
the cycles are removed by assigning missing values to all cells in a cycle.
Second, cells with a local drain direction to the outside of the map or to
a cell with a missing value (including cells that were in a cycle) are
assigned the ldd code of a pit cell (code: 5). Now, all downstream paths
on the local drain direction map end in a pit cell; this adjusted map is
saved as Result.

Notes
=====

A missing value on ldd is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Missing value creation 

See Also
========
:ref:`groupldd`, :ref:`formldd`

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Ldd = Ldd.map;
   |   initial
   |    report Result = lddrepair(Ldd);
   |   
   | • python
   |   Ldd = readmap("Ldd.map")
   |   Result = lddrepair(Ldd)

   ============================================ =========================================
   Result.map                                   Ldd.map                                  
   .. image::  ../examples/lddrepair_Result.png .. image::  ../examples/lddrepair_Ldd.png
   ============================================ =========================================

   | 

