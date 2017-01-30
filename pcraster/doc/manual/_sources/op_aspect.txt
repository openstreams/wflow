

.. index::
   single: aspect
.. _aspect:

******
aspect
******
.. topic:: aspect

   Aspects of a map using a digital elevation model

::

  Result = aspect(dem)

dem
   spatial
   scalar

Result
   spatial
   directional

Operation
=========


For each cell, calculates the direction of maximum rate of change in
elevation (aspect) on basis of the elevation of its 8 neighbours in a 3
x 3 cells window. The third-order finite difference method is used,
proposed by :ref:`Horn <horn81>`, also used by :ref:`Skidmore <skidmore89>`. Result is of a directional data type, with aspect values assigned clockwise; aspect to the top of the map is taken as zero aspect. When using aguila the aspect will be expressed in degrees or radians (depending on the global setting :literal:`--degrees` (default) or :literal:`--radians`).  

Notes
=====


If a cell has a missing value on dem, a missing value is assigned to Result, in any case.   




For each cell, the aspect is calculated using its 8 neighbours in a 3 x 3 cells
window. Elevation in all these cells must be known, else the finite
difference method can not be performed. It may occur that one of these
values is unknown: this is the case if a surrounding cell is a missing value
or if the centre cell is at the edge of the map resulting in the absence of
some surrounding cells. If this occurs the aspect operator uses a built in neighbourhood interpolator to fill these missing values or absent cells in; this will make calculation of the slope for the centre cell still possible. For each missing value cell or absent cell, the elevation is determined by taking the average elevation of non missing value cells in a 3 x 3 cell window, with the missing value cell or absent cell in the centre of the window.    

Group
=====
This operation belongs to the group of  Derivatives of elevation maps 

See Also
========
:ref:`nodirection`

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Dem = Dem.map;
   |   initial
   |    report Result = aspect(Dem);
   |   
   | • python
   |   Dem = readmap("Dem.map")
   |   Result = aspect(Dem)

   ========================================= =====================================
   Result.map                                Dem.map                              
   .. image::  ../examples/aspect_Result.png .. image::  ../examples/slope_Dem.png
   ========================================= =====================================

   | 

