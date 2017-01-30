

.. index::
   single: defined
.. _defined:

*******
defined
*******
.. topic:: defined

   Boolean TRUE for non missing values and FALSE for missing values

::

  Result = defined(expression)

expression
   spatial, non spatial
   boolean, nominal, ordinal, scalar, directional, ldd

Result
   dimension of expression
   boolean

Operation
=========


For each cell on Result returns a Boolean value where 1 is TRUE and 0 is FALSE: if the cell value on expression is not a missing value a 1 (TRUE) is assigned to the corresponding cell on Result, if the cell value on expression is a missing value a 0 (FALSE) is assigned to the corresponding cell on Result.   

Notes
=====


Group
=====
This operation belongs to the group of  Missing value creation 

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result = defined(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = defined(Expr)

   ========================================== ========================================
   Result.map                                 Expr.map                                
   .. image::  ../examples/defined_Result.png .. image::  ../examples/defined_Expr.png
   ========================================== ========================================

   | 

