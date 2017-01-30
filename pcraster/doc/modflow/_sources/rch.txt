RCH package
^^^^^^^^^^^
The RCH package specifies the areally distributed recharge. To enable the recharge package either use

.. code-block:: c

   res = mf::setRecharge(recharge, NRCHOP);

where

recharge [:math:`LT^-1`]
   is the name of a spatial, scalar PCRaster map containing the recharge flux and

NRCHOP [:math:`-`]
   is the recharge option code (1 - recharge is only to the top grid layer, 3 - recharge is applied to the highest active cell in each vertical column),

or to indicate the layer number where recharge is applied to use

.. code-block:: c

   res = mf::setIndicatedRecharge(recharge, layer);

where

recharge [:math:`LT^-1`]
   is the name of a spatial, scalar PCRaster map containing the recharge flux and

layer [:math:`-`]
   is the name of a spatial, nominal PCRaster map containing layer number variables defining the layer in each vertical column where recharge is applied.
