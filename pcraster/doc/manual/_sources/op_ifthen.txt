

.. index::
   single: if then
.. _ifthen:

*******
if then
*******
.. topic:: if then

   Return missing values if condition is not met.

::

  Result = if(condition then expression)

::

  Result = if(condition, expression)

expression
   spatial, non spatial
   boolean, nominal, ordinal, scalar, directional, ldd

condition
   spatial, non spatial
   boolean

Result
   spatial; if condition and expression are non spatial non spatial
   type of expression

Operation
=========


The cell values on condition are interpreted as Boolean values where 1 is TRUE and 0 is FALSE. For each cell, the cell value on condition determines whether the value of the corresponding cell on expression or a missing value is assigned to the corresponding cell on Result: if condition has a cell value 1 (TRUE) the value on expression is assigned to Result, if condition has a cell value 0 (FALSE) a missing value is assigned to Result.  

Notes
=====


A cell with missing value on condition and/or expression results in a missing value on Result at the corresponding cell. A comma between condition and expression in the command line is an alternative notation for then.  



If you want to cut an local drain direction map (data type ldd), use the
operator lddmask instead of if then. The operator if then allows for cutting an expression of data type ldd, but we advice to use it in very special cases only; it will result in an :ref:`unsound <SoundLDD>` ldd.  

Group
=====
This operation belongs to the group of  Conditional operators 

See Also
========
:ref:`ifthenelsepitfalls`,
:ref:`cover`, :ref:`defined`, :ref:`lddmask`

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Cond = Cond.map;
   |    Expr1 = Expr1.map;
   |   initial
   |    report Result = if(Cond then Expr1);
   |   
   | • python
   |   Cond = readmap("Cond.map")
   |   Expr1 = readmap("Expr1.map")
   |   Result = ifthen(Cond, Expr1)

   ========================================= ======================================= ========================================
   Result.map                                Cond.map                                Expr1.map                               
   .. image::  ../examples/ifthen_Result.png .. image::  ../examples/ifthen_Cond.png .. image::  ../examples/ifthen_Expr1.png
   ========================================= ======================================= ========================================

   | 

