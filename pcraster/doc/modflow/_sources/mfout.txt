Retrieving Modflow's output
---------------------------
Afterwards the specification of the grid Modflow can be called with

.. code-block:: c

   res = mf::run();

After a successful Modflow run the results of a stress period can be reported with the commands described in this section. If the run fails the end of the Modflow list file is displayed, the execution of the script ends.

The head and boundary values are retrieved automatically and must not be set again for the next stress period. Applying the operations to layer specified as quasi-3D confining beds will result in an error.

.. code-block:: c

   head.map = mf::getHeads(LAYER);

``head.map`` is a scalar PCRaster map containing the resulting head values of the layer ``layer``. Cells converted to dry obtain the data type missing value.

.. code-block:: c

   river.map = mf::getRiverLeakage(LAYER);

``river.map`` is a scalar PCRaster map containing the resulting river cell-by-cell flow values (in [:math:`L^3T^-1`]) of the layer ``layer``.

.. code-block:: c

   rch.map = mf::getRecharge(LAYER);

``rch.map`` is a scalar PCRaster map containing the resulting recharge cell-by-cell flow values (in [:math:`L^3T^-1]`) of the layer ``layer``.

.. code-block:: c

   drn.map = mf::getDrain(LAYER);

``drn.map`` is a scalar PCRaster map containing the resulting drain cell-by-cell flow values (in [:math:`L^3T^-1]`) of the layer ``layer``.




.. code-block:: c

   st.map = mf::getStorage(LAYER);

``st.map`` is a scalar PCRaster map containing the resulting cell-by-cell storage terms of the layer ``layer``.


.. code-block:: c

   ch.map = mf::getConstantHead(LAYER);

``ch.map`` is a scalar PCRaster map containing the resulting cell-by-cell constant head flow terms of the layer ``layer``.


.. code-block:: c

   rf.map = mf::getRightFace(LAYER);

``rf.map`` is a scalar PCRaster map containing the resulting internal cell-by-cell flows (right) of the layer ``layer``.


.. code-block:: c

   ff.map = mf::getFrontFace(LAYER);

``ff.map`` is a scalar PCRaster map containing the resulting internal cell-by-cell flows (front) of the layer ``layer``.


.. code-block:: c

   lf.map = mf::getLowerFace(LAYER);

``lf.map`` is a scalar PCRaster map containing the resulting internal cell-by-cell flows (lower) of the layer ``layer``.




