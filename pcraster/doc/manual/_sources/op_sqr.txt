

.. index::
   single: sqr
.. _sqr:

***
sqr
***
.. topic:: sqr

   Square

::

  Result = sqr(expression)

expression
   spatial, non spatial
   scalar

Result
   dimension of expression
   scalar

Operation
=========


For each cell, calculates the square of the expression cell value and assigns it to Result.  

Notes
=====


A cell with a missing value on expression is assigned a missing value on Result.  

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
   |    report Result = sqr(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = sqr(Expr)

   ====================================== ====================================
   Result.map                             Expr.map                            
   .. image::  ../examples/sqr_Result.png .. image::  ../examples/sqr_Expr.png
   ====================================== ====================================

   | 

