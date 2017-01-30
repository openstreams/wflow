

.. index::
   single: **
.. _sterster:

****
\*\*
****
.. topic:: \*\*

   Nth power of a first expression, where N is the value of a second expression

::

  Result = expression ** power

expression
   spatial, non spatial
   scalar

power
   spatial, non spatial
   scalar

Result
   spatial if expression or power is spatial, else non spatial
   scalar

Operation
=========


For each cell, raises the cell values on expression to the Nth power, where N is the cell value on power. The result of this calculation is assigned to the corresponding cell on Result.  

Notes
=====


A cell with a value 0 on expression and a value less than or equal to 0 on power is assigned a missing value on Result. Also, a cell with a value less than 0 on expression and a negative whole number on power is assigned a missing value on Result.  



A cell with a missing value on expression and/or power is assigned a missing value on Result.  

Result \*\*= power is an alternative notation for Result = Result \*\* power  

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
   |    Power = Power.map;
   |   initial
   |    report Result = Expr ** Power;
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Power = readmap("Power.map")
   |   Result = Expr ** Power

   =========================================== ========================================= ==========================================
   Result.map                                  Expr.map                                  Power.map                                 
   .. image::  ../examples/sterster_Result.png .. image::  ../examples/sterster_Expr.png .. image::  ../examples/sterster_Power.png
   =========================================== ========================================= ==========================================

   | 

