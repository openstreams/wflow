

.. index::
   single: atan
.. _atan:

****
atan
****
.. topic:: atan

   Inverse tangent

::

  Result = atan(expression)

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


For each cell, calculates the inverse tangent of the cell value on
expression and assigns it to Result.  

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
   |    report Result = atan(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = atan(Expr)

   ======================================= =====================================
   Result.map                              Expr.map                             
   .. image::  ../examples/atan_Result.png .. image::  ../examples/atan_Expr.png
   ======================================= =====================================

   | 

