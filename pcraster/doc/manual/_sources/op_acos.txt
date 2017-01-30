

.. index::
   single: acos
.. _acos:

****
acos
****
.. topic:: acos

   Inverse cosine

::

  Result = acos(expression)

expression
   spatial, non spatial
   scalar

Result
   dimension of expression
   direction

Options
=======

affects Result only if expression is a number:

:literal:`--degrees`
   direction is given in degrees (default)

:literal:`--radians`
   direction is given in radians



Operation
=========


This operator calculates for each cell the inverse cosine of the cell value on
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
   |   #! --degrees
   |   binding
   |    Result = Result.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result = acos(Expr);
   |   
   | • python
   |   setglobaloption("degrees")
   |   Expr = readmap("Expr.map")
   |   Result = acos(Expr)

   ======================================= =====================================
   Result.map                              Expr.map                             
   .. image::  ../examples/acos_Result.png .. image::  ../examples/acos_Expr.png
   ======================================= =====================================

   | 

