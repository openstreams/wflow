.. _timeinputmodulo:

***************
timeinputmodulo
***************
.. index::
   single: timeinputmodulo
.. topic:: timeinputmodulo

   Returns a map for each timestep using a modulo index into a shorter map stack.

::

   Result = timeinputmodulo(SuffixMap, highestTimestepAvail)

SuffixMap
  boolean, nominal, ordinal, scalar, directional, ldd, spatial, non spatial
highestTimestepAvail
  ordinal, non spatial
Result
  type of SuffixMap, spatial

Operation
=========
Works like the :ref:`timeinput` except that in timesteps larger than highestTimestepAvail a stack item I is read
according to the formula:

   m = time() mod highestTimestepAvail

   I = if (m eq 0, then highestTimestepAvail else m)

Notes
=====
highestTimestepAvail must be a single integer and cannot be a computed variable.

Example
=======

A dynamic modelling script with 7 timesteps contains the function:

 |  rainFall = timeinputmodulo(rain,3);

The maps assigned to rainFaill are as follows:

  =========== ============
  in timestep map
  =========== ============
  1           rain0000.001
  2           rain0000.002
  3           rain0000.003
  4           rain0000.001
  5           rain0000.002
  6           rain0000.003
  7           rain0000.001
  =========== ============
