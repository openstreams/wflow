

.. index::
   single: areaorder
.. _areaorder:

*********
areaorder
*********
.. topic:: areaorder

   Within each area ordinal numbers to cells in ascending order

::

  Result = areaorder(expression, areaclass)

expression
   spatial
   scalar, ordinal

areaclass
   spatial
   boolean, nominal, ordinal

Result
   spatial
   scalar

Operation
=========


Let :emphasis:`n` be the number of non missing value cells on expression for a specific value of areaclass.  These cell values are set in order and numbered on Result in ascending order: the cell with the smallest value on expression is assigned a 1 and the cell with the largest value is assigned a number :emphasis:`n`. Cells on expression with identical values are assigned consecutive, unique numbers; the order in which these cells are numbered is arbitrarily chosen.  

Notes
=====


A cell with missing value on expression or areaclass is not considered in the order operation; it is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Order 

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Expr = Expr.map;
   |    AreaClass = AreaClass.map;
   |   initial
   |    report Result = areaorder(
   |    Expr,
   |    AreaClass) ;
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   AreaClass = readmap("AreaClass.map")
   |   Result = areaorder(
   |    Expr,
   |    AreaClass) 

   ============================================ ========================================== ===============================================
   Result.map                                   Expr.map                                   AreaClass.map                                  
   .. image::  ../examples/areaorder_Result.png .. image::  ../examples/areaorder_Expr.png .. image::  ../examples/areaorder_AreaClass.png
   ============================================ ========================================== ===============================================

   | 

