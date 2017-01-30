

.. index::
   single: mapmaximum
.. _mapmaximum:

**********
mapmaximum
**********
.. topic:: mapmaximum

   Maximum cell value

::

  Result = mapmaximum(expression)

expression
   spatial
   ordinal, scalar

Result
   non spatial
   type of expression

Operation
=========


Determines the maximum cell value of the expression cell values and assigns this value to Result (non spatial).  

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
   |    report Result = mapmaximum(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = mapmaximum(Expr)

   ============================================= ===========================================
   Result.map                                    Expr.map                                   
   .. image::  ../examples/mapmaximum_Result.png .. image::  ../examples/mapmaximum_Expr.png
   ============================================= ===========================================

   | 

