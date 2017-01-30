

.. index::
   single: extentofview
.. _extentofview:

************
extentofview
************
.. topic:: extentofview

   Total length of the lines in a number of directions from the cell under consideration to the first cell with a different value.

::

  Result = extentofview(classes, nrdirections)

nrdirections
   non spatial
   scalar

classes
   spatial
   boolean, nominal, ordinal

Result
   spatial
   scalar

Options
=======
:literal:`--unittrue` or :literal:`--unitcell`

:literal:`--unittrue`
   distance is measured in true distance (default)

:literal:`--unitcell`
   distance is measured in number of cell lengths



Operation
=========
For each cell and for each direction this function determines the distance untill a cell with a different value than the current cell value is encountered. For each cell these distances are summed. To calculate the average distance, divide the result by the number of distances. This average extent of view can be used as an indicator for the openness of the landscape, for example.

Group
=====
This operation belongs to the group of  Neigbourhood operator; operators for visibility analysis 

Examples
========
