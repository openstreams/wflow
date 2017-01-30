

.. index::
   single: windowmajority
.. _windowmajority:

**************
windowmajority
**************
.. topic:: windowmajority

   Most occurring cell value within a specified square neighbourhood

::

  Result = windowmajority(expression, windowlength)

windowlength
   spatial, non spatial
   scalar

expression
   spatial
   boolean, nominal, ordinal

Result
   spatial
   type of expression

Options
=======
:literal:`--unittrue` or -unitcell

:literal:`--unitcell`
   windowlength is measured in true length (default)

:literal:`--unitcell`
   windowlength is measured in number of cell lengths



Operation
=========


For each cell a square window with the cell in its centre is defined by
windowlength. The windowlength is the length of the window in horizontal and vertical directions. For each cell on expression, the most often occurring cell value within its window is determined and assigned to the corresponding cell on Result. Both cells on expression which are entirely in the window and cells which are partly in the window are considered. At a cell, if two or more values occur the same (largest) number of times in its window, the windowlength at that cell is increased with two times the celllength of expression. The most often occurring cell value in this enlarged window is assigned to the cell under consideration. If this enlarged window still does not result in one most often occurring value, the window is progressively enlarged (with steps of two times the cellsize of expression) until a most often occuring cell value is found in the window.  

Notes
=====


The cell value on windowlength should be greater than 0, else a missing value is assigned to the corresponding cell on Result.   



A cell on windowlength with a missing value results in a missing value on Result at the corresponding cell. However, if a missing value on windowlength occurs which is not the centre cell of the window the value on expression in that cell :emphasis:`is` included in the computation of the average of the cell values in the window.  

Group
=====
This operation belongs to the group of  Neigbourhood operators; window operators 

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result2 = Result2.map;
   |    Expr = Expr.map;
   |    WinLen2 = WinLen2.map;
   |   initial
   |    report Result2 = windowmajority( Expr, WinLen2);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   WinLen2 = readmap("WinLen2.map")
   |   Result2 = windowmajority( Expr, WinLen2)

   ================================================== =============================================== =================================================
   Result2.map                                        Expr.map                                        WinLen2.map                                      
   .. image::  ../examples/windowmajority_Result2.png .. image::  ../examples/windowmajority_Expr.png .. image::  ../examples/windowaverage_WinLen2.png
   ================================================== =============================================== =================================================

   | 

#. 
   | • pcrcalc
   |   binding
   |    Result1 = Result1.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result1 = windowmajority( Expr, 6);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result1 = windowmajority( Expr, 6)

   ================================================== ===============================================
   Result1.map                                        Expr.map                                       
   .. image::  ../examples/windowmajority_Result1.png .. image::  ../examples/windowmajority_Expr.png
   ================================================== ===============================================

   | 

