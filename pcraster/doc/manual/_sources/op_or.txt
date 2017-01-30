

.. index::
   single: or
.. index::
   single: |
.. index::
   single: pcror
.. _or:

*************
or, \|, pcror
*************
.. topic:: or, \|, pcror

   Boolean-or operation

::

  Result = expression1 or expression2 # pcrcalc
  Result = expression1 | expression2 # python
  Result = pcror(expression1,expression2) # python

expression1
   spatial, non spatial
   boolean

expression2
   spatial, non spatial
   boolean

Result
   spatial; non spatial if expression1 and expression2 are non spatial
   boolean

Operation
=========


The cell values on expression1 and expression2 are interpreted as Boolean values; where 1 is TRUE and 0 is FALSE. For each cell the Boolean or evaluation is performed: if expression1 or expression2 or both have a cell value 1 (TRUE) Result has a cell value 1 (TRUE) on the corresponding cell; if expression1 and expression2 have a cell value 0 (FALSE) Result has cell value 0 (FALSE).   

.. _tOR:

.. table:: Cross table of the or operator.

    +-----------------+-----------+------+
    |or               |expression1|Result|
    +-----------------+-----------+------+
    |expression2 True |True       |True  |
    +-----------------+-----------+------+
    |expression2 True |False      |True  |
    +-----------------+-----------+------+
    |expression2 False|True       |True  |
    +-----------------+-----------+------+
    |expression2 False|False      |False |
    +-----------------+-----------+------+

Notes
=====


A cell with missing value on expression1 and/or expression2 is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Boolean operators 

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Expr1 = Expr1.map;
   |    Expr2 = Expr2.map;
   |   initial
   |    report Result = Expr1 or Expr2;
   |   
   | • python
   |   Expr1 = readmap("Expr1.map")
   |   Expr2 = readmap("Expr2.map")
   |   Result = Expr1 | Expr2

   ===================================== ===================================== =====================================
   Result.map                            Expr1.map                             Expr2.map                            
   .. image::  ../examples/or_Result.png .. image::  ../examples/and_Expr1.png .. image::  ../examples/and_Expr2.png
   ===================================== ===================================== =====================================

   | 

