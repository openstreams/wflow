

.. index::
   single: roundup
.. _roundup:

*******
roundup
*******
.. topic:: roundup

   Rounding up of cellvalues to whole numbers

::

  Result = roundup(expression)

expression
   spatial, non spatial
   scalar

Result
   spatial, non spatial
   scalar

Operation
=========


For each cell, the expression value is rounded upwards: the next whole value that is greater than or equal to the value on expression is assigned to Result.  

Notes
=====


Input values can be positive or negative.
A cell with missing value on expression is assigned a missing value on Result.  

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
   |    report Result = roundup(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = roundup(Expr)

   ========================================== ==========================================
   Result.map                                 Expr.map                                  
   .. image::  ../examples/roundup_Result.png .. image::  ../examples/rounddown_Expr.png
   ========================================== ==========================================

   | 

