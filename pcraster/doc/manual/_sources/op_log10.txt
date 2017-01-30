

.. index::
   single: log10
.. _log10:

*****
log10
*****
.. topic:: log10

   Log \ :sub:`10`

::

  Result = log10(expression)

expression
   spatial, non spatial
   scalar

Result
   dimension of expression
   scalar

Operation
=========


For each cell, calculates the \ :sub:`10` logarithm of the cell value on expression and assigns it to the corresponding cell on Result.  

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
   |    report Result = log10(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = log10(Expr)

   ======================================== ======================================
   Result.map                               Expr.map                              
   .. image::  ../examples/log10_Result.png .. image::  ../examples/log10_Expr.png
   ======================================== ======================================

   | 

