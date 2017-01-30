DRN package
^^^^^^^^^^^
The DRN package specifies the use of drains in a simulation. To enable the drain package use

.. code-block:: c

   res = mf::setDrain(elevation, conductance, LAYER);

where

elevation [:math:`L`]
   is the name of a spatial, scalar PCRaster map containing drain elevation values;

conductance [:math:`L^2T^-1`]
   is the name of a spatial, scalar PCRaster map containing drain conductance values and

LAYER [:math:`-`]
   is the layer number the map values will be assigned to.

Water seepage into a drain is calculated if the conductance value in a cell is larger than zero.
