

.. index::
   single: order
.. _order:

*****
order
*****
.. topic:: order

   Ordinal numbers to cells in ascending order

::

  Result = order(expression)

expression
   spatial
   scalar, ordinal

Result
   spatial
   scalar

Operation
=========


Let :emphasis:`n` be the number of non missing value cells on expression. These cell values are set in order and numbered on Result in ascending order: the cell with the smallest value on expression is assigned a 1 and the cell with the largest value is assigned a number :emphasis:`n`. Cells on expression with identical values are assigned consecutive, unique numbers; the order in which these cells are numbered is arbitrarily chosen.  

Notes
=====


A cell with missing value on expression is not considered in the order operation; it is assigned a missing value on Result.  

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
   |   initial
   |    report Result = order(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = order(Expr)

   ======================================== =====================================
   Result.map                               Expr.map                             
   .. image::  ../examples/order_Result.png .. image::  ../examples/succ_Expr.png
   ======================================== =====================================

   | 

