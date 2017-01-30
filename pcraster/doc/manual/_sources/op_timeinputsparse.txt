.. _timeinputsparse:

***************
timeinputsparse
***************
.. index::
   single: timeinputsparse
.. topic:: timeinputsparse

   Returns a map for each timestep where map-timesteps combinations can be missing.

::

  Result = timeinputsparse(SuffixOfMap)


SuffixOfMap
 boolean, nominal, ordinal, scalar, directional, ldd, spatial

Result
 type of SuffixMap, spatial

Operation
=========

timeinputsparse performs the same as :ref:`timeinput` with the difference that timeinputsparse
use the following algoritm if an input map is missing:

1. If the first input map is available at timestep A, then use this first map for the timestep 1 to A.
2. If an input map is missing for timestep B, but a map is available for a prior timestep then use the prior map.

Example
=======

A dynamic modelling script with 7 timesteps contains the function:

 |  rainFall = timeinputsparse(rain);

But the filesystem only contains:

 | rain0000.003
 | rain0000.006

These 2 maps are assigned to rainFaill as follows:

  =========== ============
  in timestep map
  =========== ============
  1           rain0000.003
  2           rain0000.003
  3           rain0000.003
  4           rain0000.003
  5           rain0000.003
  6           rain0000.006
  7           rain0000.006
  =========== ============
