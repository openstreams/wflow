

.. index::
   single: rounddown
.. _rounddown:

*********
rounddown
*********
.. topic:: rounddown

   Rounding down of cellvalues to whole numbers

::

  Result = rounddown(expression)

expression
   spatial, non spatial
   scalar

Result
   dimension of expression
   scalar

Operation
=========


For each cell, the value on expression is rounded downwards: the next whole value that is less than or equal to the value on expression is assigned to Result.  

Notes
=====


Input values can be positive or negative. 




A cell with a missing value on expression is assigned a missing value on  Result.  

Group
=====
This operation belongs to the group of  Rounding 

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result = rounddown(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = rounddown(Expr)

   ============================================ ==========================================
   Result.map                                   Expr.map                                  
   .. image::  ../examples/rounddown_Result.png .. image::  ../examples/rounddown_Expr.png
   ============================================ ==========================================

   | 

