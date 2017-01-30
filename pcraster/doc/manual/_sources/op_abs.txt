

.. _abs:

***
abs
***
.. index::
   single: abs
.. topic:: abs

   Absolute value

::

  Result = abs(expression)

expression
   spatial, non spatial
   scalar

Result
   dimension of expression
   scalar

Operation
=========


For each cell, calculates the absolute value of the expression cell value and assigns it to Result.  

Notes
=====


A cell with missing value on expression is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Arithmetic operators 

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result = abs(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = abs(Expr)

   ====================================== ====================================
   Result.map                             Expr.map                            
   .. image::  ../examples/abs_Result.png .. image::  ../examples/abs_Expr.png
   ====================================== ====================================

   | 

