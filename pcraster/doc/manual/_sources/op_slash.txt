

.. index::
   single: /
.. index::
   single: div
.. _slash:

********
/ or div
********
.. topic:: / or div

   Division

::

  Result = expression1 / expression2

::

  Result = div(expression1, expression2)

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


For each cell, the value on expression1 is divided by the value on expression2. This quotient is assigned to the corresponding cell on Result.  

Notes
=====


A cell with 0 on expression2 is assigned a missing value on Result. A cell with missing value on expression1 and/or expression2 is assigned a missing value on Result.  

div is an alternative notation for /.  



Result /= expression1 is an alternative notation for Result = Result / expression1



Group
=====
This operation belongs to the group of  Arithmetic operators 

See Also
========
:ref:`idiv`, :ref:`mod`

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Expr1 = Expr1.map;
   |    Expr2 = Expr2.map;
   |   initial
   |    report Result = Expr1 / Expr2;
   |   
   | • python
   |   Expr1 = readmap("Expr1.map")
   |   Expr2 = readmap("Expr2.map")
   |   Result = Expr1 / Expr2

   ======================================== ======================================= =======================================
   Result.map                               Expr1.map                               Expr2.map                              
   .. image::  ../examples/slash_Result.png .. image::  ../examples/slash_Expr1.png .. image::  ../examples/slash_Expr2.png
   ======================================== ======================================= =======================================

   | 

