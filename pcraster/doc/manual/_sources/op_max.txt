

.. index::
   single: max
.. _max:

***
max
***
.. topic:: max

   Maximum value of multiple expressions

::

  Result = max(expression1, expression2, ..., expressionN)

expression
   spatial, non spatial
   ordinal, scalar; expression1, expression2,...expressionN must have the same data type

Result
   spatial; if all expression1, expression2,...expressionN are non spatial: non spatial
   type of expression1,expression2,...expressionN

Operation
=========


For each cell, the maximum value of expression1, expression2,...expressionN is determined and assigned to the corresponding cell on Result. As many expressions can be specified as needed.  

Notes
=====


A cell with missing value on one or more expressions results in a missing
value on the corresponding cell on Result.  

Group
=====
This operation belongs to the group of  Maximize, minimize 

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result1 = Result1.map;
   |    Expr1 = Expr1.map;
   |    Expr2 = Expr2.map;
   |   initial
   |    report Result1 = max(Expr1,Expr2);
   |   
   | • python
   |   Expr1 = readmap("Expr1.map")
   |   Expr2 = readmap("Expr2.map")
   |   Result1 = max(Expr1,Expr2)

   ======================================= ===================================== =====================================
   Result1.map                             Expr1.map                             Expr2.map                            
   .. image::  ../examples/max_Result1.png .. image::  ../examples/max_Expr1.png .. image::  ../examples/max_Expr2.png
   ======================================= ===================================== =====================================

   | 

