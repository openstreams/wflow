.. _horizontan:

**********
horizontan
**********
.. index::
   single: horizontan
.. topic:: horizontan

   Calculates the maximum tangent of the angles of neighbouring cells in the direction of the sun.


::

   result = horizontan(dem, viewAngle)

dem
  scalar; spatial

viewAngle
  directional; spatial, non-spatial

result
  scalar; spatial


Operation
=========
To determine whether a cell receives direct radiation, a critical angle for each cell is required.
When the solar angle is larger than the critical angle, the cell receives direct solar radiation.
This function calculates the maximum tangent of the angles of neighbouring cells in the direction of the sun (``viewAngle``, the solar azimuth).


Notes
=====
Since the function can only be executed if there are neighbouring cells, no value can be given to cells at the edge of the ``dem`` map.
They become missing values.
