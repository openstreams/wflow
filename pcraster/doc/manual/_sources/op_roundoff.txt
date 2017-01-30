

.. index::
   single: roundof
.. _roundoff:

********
roundoff
********
.. topic:: roundoff

   Rounding off of cellvalues to whole numbers

::

  Result = roundoff(expression)

expression
   spatial, non spatial
   scalar

Result
   dimension of expression
   scalar

Operation
=========


For each cell, the value on expression is rounded off: the whole value whose difference with the value on expression is smallest is assigned to Result. If a value on expression is exactly halfway between two whole numbers, the number with the greatest absolute magnitude is assigned to Result.  

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
   |    report Result = roundoff(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = roundoff(Expr)

   =========================================== ==========================================
   Result.map                                  Expr.map                                  
   .. image::  ../examples/roundoff_Result.png .. image::  ../examples/rounddown_Expr.png
   =========================================== ==========================================

   | 

