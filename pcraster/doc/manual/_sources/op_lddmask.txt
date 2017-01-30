

.. index::
   single: lddmask
.. _lddmask:

*******
lddmask
*******
.. topic:: lddmask

   Local drain direction map cut into a (smaller) sound local drain direction map

::

  Result = lddmask(ldd, mask)

ldd
   spatial
   ldd

mask
   spatial
   boolean

Result
   spatial
   ldd

Operation
=========


The cell values on mask are interpreted as boolean values, where 1 is TRUE and 0 is FALSE. The part of the local drain direction map ldd which you want to cut out must totally be filled with 1 (TRUE) values on mask.   



Each cell with a mask value 0 (FALSE) is assigned a missing value on Result. Each cell with a mask value 1 (TRUE) is assigned a value which corresponds with the value on ldd, except cells with a mask value 1 that have a local drain direction on ldd towards a cell with a 0 (FALSE) on mask. These last named cells are outflow cells on the edge of the new ldd, these are assigned a cell value 5, which is the ldd code for a pit.  

Notes
=====


A cell with a missing value on mask is interpreted as a mask value 0 (FALSE) and handled in that way. In addition, a cell with a missing value on ldd is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Missing value creation 

See Also
========
:ref:`groupldd`

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Ldd = Ldd.map;
   |    Mask = Mask.map;
   |   initial
   |    report Result = lddmask(Ldd,Mask);
   |   
   | • python
   |   Ldd = readmap("Ldd.map")
   |   Mask = readmap("Mask.map")
   |   Result = lddmask(Ldd,Mask)

   ========================================== ==================================== ========================================
   Result.map                                 Ldd.map                              Mask.map                                
   .. image::  ../examples/lddmask_Result.png .. image::  ../examples/accu_Ldd.png .. image::  ../examples/lddmask_Mask.png
   ========================================== ==================================== ========================================

   | 

