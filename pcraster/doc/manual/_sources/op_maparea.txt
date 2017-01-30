

.. index::
   single: maparea
.. _maparea:

*******
maparea
*******
.. topic:: maparea

   Total map area

::

  Result = maparea(expression)

expression
   spatial
   boolean, nominal, ordinal, scalar, directional, ldd

Result
   non spatial
   scalar

Options
=======
:literal:`--unittrue` or :literal:`--unitcell`

:literal:`--unittrue`
   area is computed in true area (default)

:literal:`--unitcell`
   area is computed in number of cells



Operation
=========

Sums the area of the non missing value cells on expression and assigns this total area to Result (non spatial). 

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
   |    report Result = maparea(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = maparea(Expr)

   ========================================== ========================================
   Result.map                                 Expr.map                                
   .. image::  ../examples/maparea_Result.png .. image::  ../examples/maparea_Expr.png
   ========================================== ========================================

   | 

