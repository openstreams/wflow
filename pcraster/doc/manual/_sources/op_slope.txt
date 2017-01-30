

.. index::
   single: slope
.. _slope:

*****
slope
*****
.. topic:: slope

   Slope of cells using a digital elevation model

::

  Result = slope(dem)

dem
   spatial
   scalar

Result
   spatial
   scalar

Options
=======
:literal:`--unittrue` or :literal:`--unitcell` (see also notes)

:literal:`--unittrue`
   horizontal and vertical scale is measured in true distance; values on dem are interpreted as real heights (default).

:literal:`--unitcell`
   horizontal and vertical scale is measured in number of cell lengths; values on dem are interpreted as number of cell lengths.



Operation
=========


For each cell, calculates the slope on basis of the elevation dem of its eight nearest neighbours in a 3 x 3 cell window. The third-order finite difference method is used, proposed by :ref:`Horn <horn81>`, also used by :ref:`Skidmore <skidmore89>`. The slope on Result is given in dZ/dX, which is the increase in height (vertical direction dZ) per distance in horizontal direction (dX). This result value is often reffered to as a percentage. Thus if slope returns a value of 0.12, one says a slope value of 12 %.

Notes
=====


Always set the option :literal:`--unittrue`; the option :literal:`--unitcell` is only used in very very special cases. In addition, note that for a correct calculation of the slope the scales for the horizontal distance on your map and the vertical distance (height) on dem must be the same.  For instance, choose for both distances metres.    



If a cell has a missing value on dem, a missing value is assigned to Result, in any case.   



For each cell, the slope is calculated using its 8 neighbours in a 3 x 3 cells
window. Elevation in all these cells must be known, else the finite
difference method can not be performed. It may occur that one of these
values is unknown: this is the case if a surrounding cell is a missing value
or if the centre cell is at the edge of the map resulting in the absence of
some surrounding cells. If this occurs the slope operator uses a built in neighbourhood interpolator to fill these missing values or absent cells in; this will make calculation of the slope for the centre cell still possible. For each missing value cell or absent cell, the elevation is determined by taking the average elevation of non missing value cells in a 3 x 3 cell window, with the missing value cell or absent cell in the centre of the window.   

Group
=====
This operation belongs to the group of  Derivatives of elevation maps 

See Also
========
:ref:`nodirection`, :ref:`aspect`

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Dem = Dem.map;
   |   initial
   |    report Result = slope(Dem);
   |   
   | • python
   |   Dem = readmap("Dem.map")
   |   Result = slope(Dem)

   ======================================== =====================================
   Result.map                               Dem.map                              
   .. image::  ../examples/slope_Result.png .. image::  ../examples/slope_Dem.png
   ======================================== =====================================

   | 

