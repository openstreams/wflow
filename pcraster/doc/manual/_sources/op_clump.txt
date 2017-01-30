

.. index::
   single: clump
.. _clump:

*****
clump
*****
.. topic:: clump

   Contiguous groups of cells with the same value ('clumps')

::

  Result = clump(expression)

expression
   spatial
   boolean, nominal, ordinal

Result
   spatial
   nominal

Options
=======
:literal:`--diagonal` or :literal:`--nondiagonal`

:literal:`--diagonal`
   cells of the same value are grouped if they are within the immediate 8-cell neighbourhoodnof each other. This includes if they are to the right or left, above or below,nor are diagonal to each other (eight nearest neighbours, default).

:literal:`--nondiagonal`
   cells of the same value are grouped only if the cells are directly to the rightnor left, or above or below each other (four nearest neighbours).



Operation
=========


Cells that have the same value on expression and are neighbours are grouped. Every group of cells satisfying these conditions is assigned a unique value on Result. Cells without neighbours with the same value on expression are also assigned a unique value on Result. The kind of connectivity needed for cells to be neighbours is specified by the option.  

Notes
=====


Cells with a missing value on expression are assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Area operators 

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result1 = Result1.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result1 = clump( Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result1 = clump( Expr)

   ========================================= ======================================
   Result1.map                               Expr.map                              
   .. image::  ../examples/clump_Result1.png .. image::  ../examples/clump_Expr.png
   ========================================= ======================================

   | 

#. 
   | • pcrcalc
   |   #! --nondiagonal
   |   binding
   |    Result2 = Result2.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result2 = clump( Expr);
   |   
   | • python
   |   setglobaloption("nondiagonal")
   |   Expr = readmap("Expr.map")
   |   Result2 = clump( Expr)

   ========================================= ======================================
   Result2.map                               Expr.map                              
   .. image::  ../examples/clump_Result2.png .. image::  ../examples/clump_Expr.png
   ========================================= ======================================

   | 

