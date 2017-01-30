.. _window4total:

************
window4total
************
.. index::
   single: window4total
.. topic:: window4total

    Sum the values of the four surrounding cells.

::

   Result = window4total(expression)

expression
  scalar, spatial

Result
  scalar, spatial

Operation
=========
For each cell sums the values of the four cells which lie above, below, left and right of the cell and assigns the sum to the cell. 

Notes
=====
If there is no cell or a cell with a missing value in any of the directions it is taken as a zero.

Examples
========

#. 
   | • pcrcalc
   |   binding
   |    Result1 = Result1.map;
   |    Expr1 = Expr1.map;
   |   initial
   |    report Result1 = window4total(Expr1);
   |   
   | • python
   |   Expr1 = readmap("Expr1.map")
   |   Result1 = window4total(Expr1)

   ================================================ =====================================
   Result1.map                                      Expr1.map                            
   .. image::  ../examples/window4total_Result1.png .. image::  ../examples/max_Expr1.png
   ================================================ =====================================

   | 

