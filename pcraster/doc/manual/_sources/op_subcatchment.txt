

.. index::
   single: subcatchment
.. _subcatchment:

************
subcatchment
************
.. topic:: subcatchment

   (Sub-)Catchment(s) (watershed, basin) of each one or more specified cells

::

  Result = subcatchment(ldd, points)

ldd
   spatial
   ldd

points
   spatial
   boolean, nominal, ordinal

Result
   spatial
   type of points

Operation
=========


The local drain direction for each cell is defined by ldd. For each non zero value on points its catchment is determined and all cells in its catchment are assigned this non zero value. This procedure is performed for all cells with a non zero value on points, but there is one important exception: subcatchments are also identified: if the catchment of a non zero cell on points contains another non zero cell value on points, the (smaller) catchment of the latter non zero cell is identified instead of the (larger) enclosing catchment.   



The operation is performed as follows: for each cell its downstream path
is determined which consists of the consecutively neighbouring
downstream cells on ldd. On Result each cell is assigned the non zero points cell value which is on its path and which is nearest downstream. If all cells on the downstream path of a cell have a value 0 on points a 0 is assigned to the cell on Result.  

Notes
=====


A cell with missing value on ldd is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Neighbourhood operators; local drain directions 

See Also
========
:ref:`catchment`, :ref:`lddmask`

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Ldd = Ldd.map;
   |    Points = Points.map;
   |   initial
   |    report Result = subcatchment(Ldd,Points);
   |   
   | • python
   |   Ldd = readmap("Ldd.map")
   |   Points = readmap("Points.map")
   |   Result = subcatchment(Ldd,Points)

   =============================================== ==================================== ============================================
   Result.map                                      Ldd.map                              Points.map                                  
   .. image::  ../examples/subcatchment_Result.png .. image::  ../examples/accu_Ldd.png .. image::  ../examples/catchment_Points.png
   =============================================== ==================================== ============================================

   | 

