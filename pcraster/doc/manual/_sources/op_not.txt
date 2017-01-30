.. index::
   single: not
.. index::
   single: ~
.. index::
   single: pcrnot
.. _not:

**************
not, ~, pcrnot
**************
.. topic:: not, ~, pcrnot

   Boolean-not operation

::

  Result = not expression # (pcrcalc)
  Result = ~ expression # (python)
  Result = pcrnot(expression) # (python)

expression
   spatial, non spatial
   boolean

Result
   dimension of expression
   boolean

Operation
=========


The cell values on expression are interpreted as Boolean values; where 1 is TRUE and 0 is FALSE. For each cell the Boolean not evaluation is performed: if expression has a cell value 1 (TRUE) Result has a cell value 0 (FALSE) on the corresponding cell; if expression has cell value 0 (FALSE) Result has cell value 1 (TRUE).  

Notes
=====


A cell with missing value on expression is assigned a missing value on Result.   

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
   |   initial
   |    report Result = not Expr1;
   |   
   | • python
   |   Expr1 = readmap("Expr1.map")
   |   Result = ~ Expr1

   ====================================== =====================================
   Result.map                             Expr1.map                            
   .. image::  ../examples/not_Result.png .. image::  ../examples/and_Expr1.png
   ====================================== =====================================

   | 

