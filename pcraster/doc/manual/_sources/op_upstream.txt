

.. index::
   single: upstream
.. _upstream:

********
upstream
********
.. topic:: upstream

   Sum of the cell values of its first upstream cell(s)

::

  Result = upstream(ldd, material)

ldd
   spatial
   ldd

material
   spatial, non spatial
   scalar

Result
   spatial
   scalar

Operation
=========


For each cell the neighbour cells that have a local drain direction on
ldd towards the cell are determined. These are cells that drain directly to the cell. On Result the cell is assigned the sum of the material values of these first upstream cells. This is done for each cell.  

Notes
=====


A cell with a missing value on ldd or material is assigned a missing value on Result. Additionally the downstream neighbour of a missing value cell on material is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Neighbourhood operators; local drain directions 

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Ldd = Ldd.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result = upstream(Ldd, Expr);
   |   
   | • python
   |   Ldd = readmap("Ldd.map")
   |   Expr = readmap("Expr.map")
   |   Result = upstream(Ldd, Expr)

   =========================================== ==================================== ===========================================
   Result.map                                  Ldd.map                              Expr.map                                   
   .. image::  ../examples/upstream_Result.png .. image::  ../examples/accu_Ldd.png .. image::  ../examples/downstream_Expr.png
   =========================================== ==================================== ===========================================

   | 

