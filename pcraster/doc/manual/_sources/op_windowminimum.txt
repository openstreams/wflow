

.. index::
   single: windowminimum
.. _windowminimum:

*************
windowminimum
*************
.. topic:: windowminimum

   Minimum value within a specified square neighbourhood

::

  Result = windowminimum(expression, windowlength)

expression
   spatial
   ordinal, scalar

windowlength
   spatial, non spatial
   scalar

Result
   spatial
   expression

Options
=======
:literal:`--unittrue` or :literal:`--unitcell`

:literal:`--unittrue`
   windowlength is measured in true length (default)

:literal:`--unitcell`
   windowlength is measured in number of cell lengths



Operation
=========


For each cell a square window with the cell in its centre is defined by
windowlength. The windowlength is the length of the window in horizontal and vertical directions. For each cell on expression, the minimum cell value within its window is determined and assigned to the corresponding cell on Result. Both cells on expression which are entirely in the window and cells which are partly in the window are considered.  

Notes
=====


The cell value on windowlength should be greater than 0, else a missing value is assigned to the corresponding cell on Result.  



A cell on windowlength with a missing value results in a missing value on Result at the corresponding cell. However, if a missing value on windowlength occurs in a cell in the window which is not the centre cell of the window the value on expression in that cell :emphasis:`is` considered for determination of the minimum cell value in the window.  

Group
=====
This operation belongs to the group of  Neigbourhood operators; window operators 

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result1 = Result1.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result1 = windowminimum( Expr, 6);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result1 = windowminimum( Expr, 6)

   ================================================= ==============================================
   Result1.map                                       Expr.map                                      
   .. image::  ../examples/windowminimum_Result1.png .. image::  ../examples/windowaverage_Expr.png
   ================================================= ==============================================

   | 

#. 
   | • pcrcalc
   |   binding
   |    Result2 = Result2.map;
   |    Expr = Expr.map;
   |    WinLen2 = WinLen2.map;
   |   initial
   |    report Result2 = windowminimum(Expr, WinLen2);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   WinLen2 = readmap("WinLen2.map")
   |   Result2 = windowminimum(Expr, WinLen2)

   ================================================= ============================================== =================================================
   Result2.map                                       Expr.map                                       WinLen2.map                                      
   .. image::  ../examples/windowminimum_Result2.png .. image::  ../examples/windowaverage_Expr.png .. image::  ../examples/windowaverage_WinLen2.png
   ================================================= ============================================== =================================================

   | 

