.. _spatial:

*******
spatial
*******
.. index::
   single: spatial
.. topic:: spatial

   Conversion of a non-spatial value to a spatial data type.


::

   result = spatial(expression)

expression
  spatial, non-spatial; boolean, nominal, ordinal, scalar, directional, ldd

result
  spatial


Operation
=========

If ``expression`` is a non-spatial value, ``spatial`` converts the value to a spatial PCRaster map of the given data type.


Group
=====
This operation belongs to the group of Conversion and assignment

See also
========
:ref:`secdatbasemaptype`