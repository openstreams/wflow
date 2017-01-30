

.. index::
   single: accuflux

.. _accuflux:

********
accuflux
********
.. topic:: accuflux

   Accumulated material flowing into downstream cell

::

  Resultflux = accuflux(ldd, material)

ldd
   spatial
   ldd

material
   spatial, non spatial
   scalar

Resultflux
   spatial
   scalar

Operation
=========


This operation calculates for each cell the accumulated amount of material that flows out
of the cell into its neighbouring downstream cell. This accumulated amount
is the amount of material in the cell itself plus the amount of material in
upstream cells of the cell. For each cell, the following procedure is performed: using the local
drain direction network on ldd, the catchment of a cell its outflow is determined which is made up the cell itself and all cells that drain to the cell (i.e. which are in upstream direction of the cell). The material values of all cells in the catchment are summed and send to the cell on Resultflux. This value is the amount of material which accumulates during transport in downstream direction to the outflow of the cell.  

Notes
=====


The values on material must be equal to or larger than zero. A cell with missing value on material is assigned a missing value on Resultflux. Additionally, all its downstream cells are assigned a missing value.  

Group
=====
This operation belongs to the group of  Neighbourhood operators; local drain directions 

See Also
========
:ref:`secstatneightr`, :ref:`lddmask`

Examples
========
