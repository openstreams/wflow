.. index::
   single: sin
.. _sin:

***
sin
***
.. topic:: sin

   Sine

::

  Result = sin(expression)

expression
   spatial, non spatial
   directional, scalar

Result
   dimension of expression
   scalar

Options
=======

if expression is a number: :literal:`--degrees` or :literal:`--radians`

:literal:`--degrees`
   direction is measured in degrees (default)

:literal:`--radians`
   direction is measured in radians



Operation
=========


For each cell, calculates the sine of the cell value on expression and assigns it to Result.  

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
   |    report Result = sin(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = sin(Expr)

   ====================================== ====================================
   Result.map                             Expr.map                            
   .. image::  ../examples/sin_Result.png .. image::  ../examples/sin_Expr.png
   ====================================== ====================================

   | 

