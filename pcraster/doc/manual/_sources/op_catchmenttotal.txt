

.. index::
   single: catchmenttotal
.. _catchmenttotal:

**************
catchmenttotal
**************
.. topic:: catchmenttotal

   Total catchment for the entire upstream area

::

  Result = catchmenttotal(amount, ldd)

ldd
   spatial
   ldd

amount
   spatial, non spatial
   scalar

Result
   spatial
   scalar

Operation
=========

This operation is identical to accuflux, except that accuflux does not accept negative values. Catchmenttotal calculates for each cell the accumulated amount of material that flows out
of the cell into its neighbouring downstream cell. This accumulated amount
is the amount of material in the cell itself plus the amount of material in
upstream cells of the cell. For each cell, the following procedure is performed: using the local
drain direction network on ldd, the catchment of a cell its outflow is determined which is made up the cell itself and all cells that drain to the cell (i.e. which are in upstream direction of the cell). The material values of all cells in the catchment are summed and send to the cell on Resultflux. This value is the amount of material which accumulates during transport in downstream direction to the outflow of the cell. 

Examples
========
