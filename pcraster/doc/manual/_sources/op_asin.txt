

.. index::
   single: asin
.. _asin:

****
asin
****
.. topic:: asin

   Inverse sine

::

  Result = asin(expression)

expression
   spatial, non spatial
   scalar

Result
   dimension of expression
   directional

Options
=======

if expression is a number: :literal:`--degrees` or :literal:`--radians`

:literal:`--degrees`
   direction is given in degrees (default)

:literal:`--radians`
   direction is given in radians



Operation
=========


For each cell, calculates the inverse cosine of the cell value on
expression and assigns it to Result.  

Notes
=====


The values on expression must be equal to or between -1 and 1. Cells with a value outside this range will be assigned a missing value on Result. A cell with missing value on expression is assigned a missing value on Result.  

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
   |    report Result = asin(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = asin(Expr)

   ======================================= =====================================
   Result.map                              Expr.map                             
   .. image::  ../examples/asin_Result.png .. image::  ../examples/asin_Expr.png
   ======================================= =====================================

   | 

