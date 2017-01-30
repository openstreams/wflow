.. index::
   single: tan
.. _tan:

***
tan
***
.. topic:: tan

   Tangent

::

  Result = tan(expression)

expression
   spatial, non spatial
   directional, scalar

Result
   dimension of expression
   scalar

Options
=======
:literal:`--degrees` or :literal:`--radians`

:literal:`--degrees`
   if expression is a number then the unit is degrees (default)

:literal:`--radians`
   if expression is a number then the unit is radians



Operation
=========


For each cell, calculates the tangent of expression cell value and assigns it to Result.  

Notes
=====


A cell with missing value on expression is assigned a missing value on Result.  



If expression is of directional data type, a cell on expression without a direction (cell value -1) is assigned a missing value.  

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
   |    report Result = tan(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = tan(Expr)

   ====================================== ====================================
   Result.map                             Expr.map                            
   .. image::  ../examples/tan_Result.png .. image::  ../examples/tan_Expr.png
   ====================================== ====================================

   | 

