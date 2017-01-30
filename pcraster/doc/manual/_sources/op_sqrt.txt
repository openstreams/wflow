

.. index::
   single: sqrt
.. _sqrt:

****
sqrt
****
.. topic:: sqrt

   Square root

::

  Result = sqrt(expression)

expression
   spatial, non spatial
   scalar

Result
   dimension of expression
   scalar

Operation
=========


For each cell, calculates the square root of the expression cell value and assigns it to Result.  

Notes
=====


The cell values on expression must be equal to or greater than 0. Negative values are assigned a missing value on Result.  A cell with a missing value on expression is assigned a missing value on Result.  

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
   |    report Result = sqrt(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = sqrt(Expr)

   ======================================= =====================================
   Result.map                              Expr.map                             
   .. image::  ../examples/sqrt_Result.png .. image::  ../examples/sqrt_Expr.png
   ======================================= =====================================

   | 

