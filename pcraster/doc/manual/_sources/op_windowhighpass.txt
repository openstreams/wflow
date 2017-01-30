

.. index::
   single: windowhighpass
.. _windowhighpass:

**************
windowhighpass
**************
.. topic:: windowhighpass

   Increases spatial frequency within a specified square neighbourhood

::

  Result = windowhighpass(expression, windowlength)

expression
   spatial
   scalar

windowlength
   spatial, non spatial
   scalar

Result
   spatial
   scalar

Options
=======
:literal:`--unittrue` or :literal:`--unitcell`

:literal:`--unittrue`
   windowlength is measured in true length (default)

:literal:`--unitcell`
   windowlength is measured in number of cell lengths



Operation
=========


For each cell :emphasis:`c` its windowhighpass value is computed as follows. A square window with the cell :emphasis:`c` in its centre is defined by windowlength. The windowlength is the length of the window in horizontal and vertical directions. For each cell :emphasis:`i` which is entirely or partly in the window and which is not the centre cell :emphasis:`c` the fraction of the cell in the window is determined. This is the area of the part of the cell in the window divided by the total area of a cell. Let :emphasis:`fraction(i)` be this fraction; let expression(i) be the expression value of a surrounding cell :emphasis:`i` and expression(c) the expression value of the centre cell :emphasis:`c`. The windowhighpass filter value on the centre cell :emphasis:`c` is calculated according to:  

.. math::

  windowhighpass(c) = 2 * expression(c) * { \sum^n_{i=1} fraction(i) } -
                          { \sum^n_{i=1} ( fraction(i) * expression(i) ) }

where n is the number of cells i surrounding the centre cell C in the
window.
This computation is performed for all cells: for each cell the
highpass(C) value is assigned to the corresponding cell on Result.  

Notes
=====


The cell value on windowlength should be greater than 0, else a missing value is assigned to the corresponding cell on Result.  


A cell on windowlength with a missing value results in a missing value on Result at the corresponding cell. However, if a missing value on windowlength occurs in a cell which is not the centre cell of the window the expression value in that cell :emphasis:`is` included in the computation of the highpass value in the window. 

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
   |    report Result2 = windowhighpass( Expr, WinLen2);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   WinLen2 = readmap("WinLen2.map")
   |   Result2 = windowhighpass( Expr, WinLen2)

   ================================================== ============================================== =================================================
   Result2.map                                        Expr.map                                       WinLen2.map                                      
   .. image::  ../examples/windowhighpass_Result2.png .. image::  ../examples/windowaverage_Expr.png .. image::  ../examples/windowaverage_WinLen2.png
   ================================================== ============================================== =================================================

   | 

#. 
   | • pcrcalc
   |   binding
   |    Result1 = Result1.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result1 = windowhighpass( Expr, 6);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result1 = windowhighpass( Expr, 6)

   ================================================== ==============================================
   Result1.map                                        Expr.map                                      
   .. image::  ../examples/windowhighpass_Result1.png .. image::  ../examples/windowaverage_Expr.png
   ================================================== ==============================================

   | 

