BAS package
^^^^^^^^^^^
The BAS package specifies the input for the basic package. It is necessary both to set the boundary and the starting head values for each layer. If a command is not given for a layer the cells are set to the default values of 0 (inactive cell) and a initial head value of 0.

The operation for the boundary specification is

.. code-block:: c

   res = mf::setBoundary(boundaryValues, LAYER);

where

boundaryValues
   is the name of a spatial, nominal PCRaster map (-1: cell has constant head, 0: cell is inactive, 1: cell is active) and

LAYER
   is the Modflow layer number the map values will be assigned to.

The operation for the specification of the starting or initial head values

.. code-block:: c

   res = mf::setInitialHead(initialValues, LAYER);

where

initialValues
   is the name of a spatial, scalar PCRaster map containing the head values and

LAYER
   is the layer number the map values will be assigned to.

Optional operation
~~~~~~~~~~~~~~~~~~
The options for the BAS package can be specified with

.. code-block:: c

   res = mf::setNoFlowHead(VALUE);

where

VALUE
   is the scalar, non-spatial value of the head to be assigned to all no flow cells (HNOFLO).

If this operation is not used the value will be set to a default value of -999.99.
