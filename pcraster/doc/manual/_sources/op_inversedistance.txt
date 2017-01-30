.. _inversedistance:

***************
inversedistance
***************
.. index::
   single: inversedistance
.. index::
   single: interpolation
.. topic:: inversedistance

    Interpolate values

::

  Result = inversedistance( mask, points, idp,radius,maxNr)

Result
 scalar, spatial

mask
 boolean, spatial or non-spatial

points
  scalar, spatial

idp
 scalar, spatial or non-spatial

radius
  scalar, spatial or non-spatial
maxNr
  scalar, spatial or non-spatial

Options
=======

:emphasis:`--unittrue`,:emphasis:`--unitcell`

*--unittrue*
  radius is measured in true distance (default).
*--unitcell* 
  radius is measured in number of cell lengths. 

Operation
=========

For each cell, having a value 1 (TRUE) in mask, a value is interpolated using weighted average of all non missing value points on points. The weighted average method is an inverse distance scheme with a weight of d-idp.

radius allows to select only the points at a distance less or equal to the cell a value is computed. If radius is 0 or less then there is no restriction on the distance towards the cell; all points are candidate to contribute to the computation.

maxNr allows to select only the maxNr closest points as the maximum number of points used in the computation. Less points are used if if no maxNr points can be found with in the radius. If maxNr is 0 or less then all points are candidate to contribute to the computation.

A missing value is generated if the the combination of radius and maxNr setting yield no points at all.

Inversedistance is typically used to interpolate sparse point samples to a continuous surface of values. An idp value of 2 is commonly used to achieve a smooth surface. A Result with idp value of 2 is also much faster to compute then another idp value.

A common use of inversedistance is in combination with timeinput, to create a value surface each timestep. This setup will also handle the case where timeseries have missing values (1e31): in that case a point is not defined in that timestep and will not take part in the computation of the weighted average.

Notes
=====
The execution time of this function increases linear with the number of non-missing values on points.
Using all points in the interpolation by setting both radius and maxNr to 0 will gives the fasted execution time.

Examples
========

#. ::

    binding
     # a map with columns id's MV's elsewhere
     pointId = pointid.map;
     # a map with 1 (TRUE) for area of interest
     mask    = mask.map;
     # a timeseries possible with missing data (1e31) 
     inputTss= input.tss;
     ....
    dynamic
     points = timeinputscalar(inputTss, pointId);
     # use all points if radius is 0 and maxNr is 0
     result = inversedistance(mask, points, 2,0,0);



#. Since radius can be spatial we may specify a different radius for each cell. In combination with spread() we can enlarge the radius in such a way that at least 1 point is always selected in areas where no sampling points are near. Note that we add one celllength to adjust for cell roundings.

 ::

  binding
   # a map with columns id's MV's elsewhere
   pointId = pointid.map;
   # a map with 1 (TRUE) for area of interest
   mask    = mask.map;
   # a timeseries possible with missing data (1e31) 
   inputTss= input.tss;
   ....
  dynamic
   points = timeinputscalar(inputTss, pointId);
   spreadSurface = cover(points,0);
   distanceToNearestPoint = spread(spreadSurface,0,1)+celllength();
   # 5000 is the default radius
   radius = max(spreadSurface,5000);
   # use at max 5 points within radius
   result = inversedistance(mask, points, 2,radius,5);


