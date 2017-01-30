

.. index::
   single: *
.. _ster:

**
\*
**

.. topic:: \*

   Multiplication

::

  Result = expression1 * expression2

expression1
   spatial, non spatial
   scalar

expression2
   spatial, non spatial
   scalar

Result
   spatial; non spatial if expression1 and expression2 are non spatial
   scalar

Operation
=========


For each cell, the values of expression1 and expression2 are multiplied. This product is assigned to the corresponding cell on Result.  

Notes
=====

A cell with missing value on expression1 and/or expression2 is assigned a missing value on Result. 

Result \*= expression1 is an alternative notation for Result = Result \* expression1  

Group
=====
This operation belongs to the group of  Arithmetic operators 

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Expr1 = Expr1.map;
   |    Expr2 = Expr2.map;
   |   initial
   |    report Result = Expr1 * Expr2;
   |   
   | • python
   |   Expr1 = readmap("Expr1.map")
   |   Expr2 = readmap("Expr2.map")
   |   Result = Expr1 * Expr2

   ======================================= ====================================== ======================================
   Result.map                              Expr1.map                              Expr2.map                             
   .. image::  ../examples/ster_Result.png .. image::  ../examples/ster_Expr1.png .. image::  ../examples/ster_Expr2.png
   ======================================= ====================================== ======================================

   | 

