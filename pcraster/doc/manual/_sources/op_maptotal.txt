

.. index::
   single: maptotal
.. _maptotal:

********
maptotal
********
.. topic:: maptotal

   Sum of all cell values

::

  Result = maptotal(expression)

expression
   spatial
   scalar

Result
   non spatial
   scalar

Operation
=========

Sums all expression cell values and assigns this sum of these values  to Result (non spatial). 

Notes
=====


The value of Result is not correct if all cells of expression are missing value. In that case Result is assigned the value 0.  

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
   |    report Result = maptotal(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = maptotal(Expr)

   =========================================== ===========================================
   Result.map                                  Expr.map                                   
   .. image::  ../examples/maptotal_Result.png .. image::  ../examples/mapmaximum_Expr.png
   =========================================== ===========================================

   | 

