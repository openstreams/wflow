.. _influencesimplegauss:

********************
influencesimplegauss
********************
.. index::
   single: influencesimplegauss
.. topic:: influencesimplegauss

   Simple unweighted gaussian decrease of influence from sources.


::

   Result = influencesimplegauss(sources, range, threshold)

sources
  scalar, spatial
range
  scalar, spatial, non-spatial

threshold
 scalar, spatial, non-spatial


Operation
=========

Calculates the unweighted gaussian decrease of the influence of sources where range is the range at which the influence is 35%. It is calculated according to the formula:

Result = sources * exp(-distance/range)

In the formula distance is the distance from source cell to the cell which is being calculated.
The value threshold is the amount above which the influence must lie before it is calculated. If it is below threshold a value of 0 will be returned.

Notes
=====
The influence on cells with missing values in sources is not calculated.
