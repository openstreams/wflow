

.. index::
   single: mapminimum
.. _mapminimum:

**********
mapminimum
**********
.. topic:: mapminimum

   Minimum cell value

::

  Result = mapminimum(expression)

expression
   spatial
   ordinal, scalar

Result
   non spatial
   type of expression

Operation
=========


Determines the minimum cell value of expression cell values and assigns this value to Result (non spatial).  

Notes
=====

The value of Result is undefined if all cells of expression are missing value.  

Group
=====
This operation belongs to the group of  Map operators 

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result = mapminimum(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = mapminimum(Expr)

   ============================================= ===========================================
   Result.map                                    Expr.map                                   
   .. image::  ../examples/mapminimum_Result.png .. image::  ../examples/mapmaximum_Expr.png
   ============================================= ===========================================

   | 

