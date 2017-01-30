

.. index::
   single: ycoordinate
.. _ycoordinate:

***********
ycoordinate
***********
.. topic:: ycoordinate

   Y-coordinate of each Boolean TRUE cell

::

  Result = ycoordinate(expression)

expression
   spatial, non spatial
   boolean

Result
   spatial
   scalar

Options
=======
:literal:`--unittrue` or :literal:`--unitcell`

:literal:`--unittrue`
   coordinates are expressed in true distance coordinates (default)

:literal:`--unitcell`
   coordinates are expressed in unit cell lengths, where thenminimum y coordinate of a cell centre is 0.5 and the maximum y coordinate of a cell centre isny - 0.5, where y is the number of rows of cells.




assignment of coordinates

:literal:`--coorcentre`
   for each cell, the coordinate of the cell centre is assigned (default)

:literal:`--coorul`
   for each cell, the coordinate of the upper left corner is assigned

:literal:`--coorlr`
   for each cell, the coordinate of the lower right corner is assigned



Operation
=========


For each cell that has a value 1 (TRUE) on expression assigns the y coordinate of the cell to the cell on Result. Cells with a value 0 (FALSE) on expression are assigned a missing value.   

Notes
=====


A cell with a missing value on expression is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Coordinates, unique ID's 

Examples
========
#. 
   | • pcrcalc
   |   #! --coorcentre
   |   binding
   |    Result = Result.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result = ycoordinate(Expr);
   |   
   | • python
   |   setglobaloption("coorcentre")
   |   Expr = readmap("Expr.map")
   |   Result = ycoordinate(Expr)

   ============================================== ============================================
   Result.map                                     Expr.map                                    
   .. image::  ../examples/ycoordinate_Result.png .. image::  ../examples/xcoordinate_Expr.png
   ============================================== ============================================

   | 

