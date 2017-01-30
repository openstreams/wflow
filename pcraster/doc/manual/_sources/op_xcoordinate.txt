

.. index::
   single: xcoordinate
.. _xcoordinate:

***********
xcoordinate
***********
.. topic:: xcoordinate

   X-coordinate of each Boolean TRUE cell

::

  Result = xcoordinate(expression)

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
   coordinates are expressed in unit cell lengths, where the minimumnx coordinate of a cell centre is 0.5 and the maximum x coordinate of a cell centre is xn- 0.5, where x is the number of columns of cells.




assignment of coordinates

:literal:`--coorcentre`
   for each cell, the coordinate of the cell centre is assigned (default)

:literal:`--coorul`
   for each cell, the coordinate of the upper left corner is assigned

:literal:`--coorlr`
   for each cell, the coordinate of the lower right corner is assigned



Operation
=========


For each cell that has a value 1 (TRUE) on expression assigns the x coordinate of the cell to the cell on Result. Cells with a value 0 (FALSE) on expression are assigned a missing value.   

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
   |   #! --coorlr
   |   binding
   |    Result = Result.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result = xcoordinate(Expr);
   |   
   | • python
   |   setglobaloption("coorlr")
   |   Expr = readmap("Expr.map")
   |   Result = xcoordinate(Expr)

   ============================================== ============================================
   Result.map                                     Expr.map                                    
   .. image::  ../examples/xcoordinate_Result.png .. image::  ../examples/xcoordinate_Expr.png
   ============================================== ============================================

   | 

