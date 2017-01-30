.. _lookupmapstack:

**************
lookupmapstack
**************
.. index::
   single: lookupmapstack
.. topic:: lookupmapstack

   Reads a variable assigned map number from a map stack.

::

   pcrcalc Result = lookupmapstack(SuffixMap, index)

SuffixMap
  boolean, nominal, ordinal, scalar, directional, ldd, spatial

index
  ordinal, non spatial

Result type of SuffixMap
  spatial


Operation
=========
lookupmapstack read a stack item of map stack. The index argument is
computed run time, and can therefor be the result of complicated
calculations such as a lookupordinal operation. The function does
require that stack item with index 1 exists, even if it is not used at
runtime. Stack item 1 is used to check the datatype of the resulting
maps. The function may cause runtime errors if then stack item for the
computed index are not available.



Example
=======

 | lookupmapstack(mapStack,3)

will read mapStack.003
