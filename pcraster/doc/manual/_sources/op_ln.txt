

.. index::
   single: ln
.. _ln:

**
ln
**
.. topic:: ln

   Natural logarithm (\ :sub:`e`)

::

  Result = ln(expression)

expression
   spatial, non spatial
   scalar

Result
   dimension of expression
   scalar

Operation
=========


For each cell, calculates the natural logarithm (\ :sub:`e`) logarithm of the cell value on expression and assigns it to the corresponding cell on Result.  

Notes
=====


The cell values on expression must be greater than 0. Any cell value outside this range is assigned a missing value on Result.  A cell with missing value on expression is assigned a missing value on Result.  

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
   |    report Result = ln(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = ln(Expr)

   ===================================== ===================================
   Result.map                            Expr.map                           
   .. image::  ../examples/ln_Result.png .. image::  ../examples/ln_Expr.png
   ===================================== ===================================

   | 

