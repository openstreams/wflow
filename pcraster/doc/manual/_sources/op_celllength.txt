

.. index::
   single: celllength
.. _celllength:

**********
celllength
**********
.. topic:: celllength

   Horizontal and vertical length of a cell

::

  Result = celllength()

Result
   non spatial
   scalar

Options
=======
:literal:`--unittrue` or :literal:`--unitcell`

:literal:`--unittrue`
   length of cells is computed in true length of cells; Result is true length of a cell (default)

:literal:`--unitcell`
   length of cells is computed in unit cell length; Result is 1



Operation
=========


Calculates the length (horizontal or vertical) of one cell and assigns this
cell length to Result (non spatial).  

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
   |    report Result2 = celllength();
   |   
   | • python
   |   setglobaloption("unitcell")
   |   Expr = readmap("Expr.map")
   |   setclone("Expr.map");
   |   Result2 = celllength()

   ============================================== =========================================
   Result2.map                                    Expr.map                                 
   .. image::  ../examples/celllength_Result2.png .. image::  ../examples/cellarea_Expr.png
   ============================================== =========================================

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
   |    report Result1 = celllength();
   |   
   | • python
   |   setglobaloption("unittrue")
   |   Expr = readmap("Expr.map")
   |   setclone("Expr.map");
   |   Result1 = celllength()

   ============================================== =========================================
   Result1.map                                    Expr.map                                 
   .. image::  ../examples/celllength_Result1.png .. image::  ../examples/cellarea_Expr.png
   ============================================== =========================================

   | 

