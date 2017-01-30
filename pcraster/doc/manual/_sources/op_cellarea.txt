

.. index::
   single: cellarea
.. _cellarea:

********
cellarea
********
.. topic:: cellarea

   Area of one cell

::

  Result = cellarea()

Result
   non spatial
   scalar

Options
=======
:literal:`--unittrue` or :literal:`--unitcell`:

:literal:`--unittrue`
   area is computed in true area; Result is true area of a cell (default)

:literal:`--unitcell`
   area is computed in number of cells; Result is 1



Operation
=========


Calculates the area represented by one cell and assigns this cell area
to Result (non spatial).  

Notes
=====


Group
=====
This operation belongs to the group of  Map operators 

Examples
========
#. 
   | • pcrcalc
   |   #! --unitcell
   |   binding
   |    Result2 = Result2.map;
   |    Expr = Expr.map;
   |   areamap
   |    Expr.map;
   |   initial
   |    report Result2 = cellarea();
   |   
   | • python
   |   setglobaloption("unitcell")
   |   Expr = readmap("Expr.map")
   |   setclone("Expr.map");
   |   Result2 = cellarea()

   ============================================ =========================================
   Result2.map                                  Expr.map                                 
   .. image::  ../examples/cellarea_Result2.png .. image::  ../examples/cellarea_Expr.png
   ============================================ =========================================

   | 

#. 
   | • pcrcalc
   |   #! --unittrue
   |   binding
   |    Result1 = Result1.map;
   |    Expr = Expr.map;
   |   areamap
   |    Expr.map;
   |   initial
   |    report Result1 = cellarea();
   |   
   | • python
   |   setglobaloption("unittrue")
   |   Expr = readmap("Expr.map")
   |   setclone("Expr.map");
   |   Result1 = cellarea()

   ============================================ =========================================
   Result1.map                                  Expr.map                                 
   .. image::  ../examples/cellarea_Result1.png .. image::  ../examples/cellarea_Expr.png
   ============================================ =========================================

   | 

