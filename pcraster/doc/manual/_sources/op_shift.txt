.. _shift:

************
shift,shift0
************
.. index::
   single: shift
.. index::
   single: shift0
.. topic:: shift, shift0

   Shifts the value of each cell a specified number of cells in the assigned direction.

::

   Result = shift(expression, northing, westing)

::


  Result = shift0(expression, northing, westing)

expression
  boolean, nominal, ordinal, scalar, directional, ldd, spatial

northing
  scalar, non-spatial

westing
  scalar, non-spatial

Result
  datatype of expression, spatial

Operation
=========

Each cell of expression is shifted a number of cells equal to the value of northing to the north (up) if northing is positive, or to the south (down) if northing is negative and each cell of expression is shifted a number of cells equal to the value of westing to the west (left) if westing is positive, or to the east (right) if westing is negative. The values of northing and southing are rounded down to the nearest integer. The shift and shift0 operation will result in a number equal to the value of northing of rows being lost on the north side of the map and created on the south side of the map. Similarly a number of columns equal to the value of westing will be lost on the west side of the map and created on the east side of the map. The newly created rows or columns will have missing values when using shift and will be assigned 0 when using shift0.

Notes
=====

When using shift0 missing values that already exist in expression will also be assigned 0.

Examples
========

#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result = shift(Expr, 1,1);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = shift(Expr, 1,1)

   ======================================== ==============================================
   Result.map                               Expr.map                                      
   .. image::  ../examples/shift_Result.png .. image::  ../examples/windowaverage_Expr.png
   ======================================== ==============================================

   | 

#. 
   | • pcrcalc
   |   binding
   |    Result2 = Result2.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result2 = shift(Expr, -1,-1);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result2 = shift(Expr, -1,-1)

   ========================================= ==============================================
   Result2.map                               Expr.map                                      
   .. image::  ../examples/shift_Result2.png .. image::  ../examples/windowaverage_Expr.png
   ========================================= ==============================================

   | 

#. 
   | • pcrcalc
   |   binding
   |    Result1 = Result1.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result1 = shift0(Expr, -1,-1);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result1 = shift0(Expr, -1,-1)

   ========================================== ==============================================
   Result1.map                                Expr.map                                      
   .. image::  ../examples/shift0_Result1.png .. image::  ../examples/windowaverage_Expr.png
   ========================================== ==============================================

   | 

