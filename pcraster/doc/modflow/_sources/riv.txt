RIV package
^^^^^^^^^^^
The RIV package specifies the contribution or drainage of the rivers to/from the aquifer. To enable the river package use

.. code-block:: c

   res = mf::setRiver(head, bottom, conductance, LAYER);

where

head [:math:`L`]
   is the name of a spatial, scalar PCRaster map containing the head values;

bottom [:math:`L`]
   is the name of a spatial, scalar PCRaster map containing bottom values;

conductance [:math:`L^2T^-1`]
   is the name of a spatial, scalar PCRaster map containing the conductance values. The cell values must contain either 0 (no river) or the conductance value of the corresponding river cell;

LAYER [:math:`-`]
   is the layer number the map values will be assigned to.
