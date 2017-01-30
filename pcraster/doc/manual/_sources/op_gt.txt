

.. index::
   single: gt
.. index::
   single: >
.. _gt:

*******
gt or >
*******
.. topic:: gt or >

   Relational-greater-than operation

::

  Result = gt(expression1, expression2)

::

  Result = expression1 > expression2

expression1
   spatial, non spatial
   ordinal, scalar

expression2
   spatial, non spatial
   type of expression1

Result
   spatial; non spatial if expression1 and expression2 are non spatial
   boolean

Operation
=========


For each cell evaluates expression1 in relation to expression2. If the cell value on expression1 is greater than the value on expression2 Result has a cell value 1 (condition is TRUE) on the corresponding cell; if the cell value on expression1 equals or is less than the value on expression2 Result has a cell value 0 (condition is FALSE).   

Notes
=====


A cell with missing value on expression1 and/or expression2 results in a missing value on Result at the corresponding cell. The > sign is a alternative notation for gt.  

Group
=====
This operation belongs to the group of  Comparison operators 

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Expr1 = Expr1.map;
   |    Expr2 = Expr2.map;
   |   initial
   |    report Result = Expr1 > Expr2;
   |   
   | • python
   |   Expr1 = readmap("Expr1.map")
   |   Expr2 = readmap("Expr2.map")
   |   Result = Expr1 > Expr2

   ===================================== ==================================== ====================================
   Result.map                            Expr1.map                            Expr2.map                           
   .. image::  ../examples/gt_Result.png .. image::  ../examples/eq_Expr1.png .. image::  ../examples/eq_Expr2.png
   ===================================== ==================================== ====================================

   | 

