WEL package
^^^^^^^^^^^
The WEL package specifies the the use of wells in a simulation.  To enable the well package use

.. code-block:: c

   res = mf::setWell(rates, LAYER);

where

rates
   is the name of a spatial, scalar PCRaster map containing the well values. The cell values must contain either 0 (no well), positive (injection) or negative (pumping) volumetric rates and

LAYER
   is the layer number the map values will be assigned to.

Only layer containing wells need to be assigned, the remaining layer will be set to 0.
