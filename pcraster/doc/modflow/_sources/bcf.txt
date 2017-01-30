BCF package
^^^^^^^^^^^
The BCF package specifies the input for the block-centered flow package.

.. note::

   The PCRasterModflow extension calculates --- based on the user input --- values for transmissivities and the vertical leakance and passes these to Modflow. See the descriptions below for more detail.

.. tip::

   Note that Modflow uses square cells to calculate volumes. In case the PCRaster model uses projected cells the user should modify the inputs first before these are passed to Modflow; this is done by dividing the input volume (calculated with projected cell) by the Modflow cell area (square). The user should apply this to VCONT, Storage, and recharge.

The operation for the specification of the conductivity values is

.. code-block:: c

   res = mf::setConductivity(LAYTYPE, hConductivity, vConductivity, LAYER);

where

LAYTYPE
   contains the combined code for the method of computing interblock conductance (left digit; 0 or blank: harmonic mean (Modflow-88), 1: arithmetic mean, 2: logarithmic mean, 3: arithmetic mean of saturated thickness and logarithmic mean hydraulic conductivity) and the layer type (LAYCON, right digit; 0: confined, 1: unconfined, 2: confined/unconfined (transmissivity constant), 3: confined/unconfined (transmissivity varies));

hConductivity
   is the name of a spatial, scalar PCRaster map containing the horizontal conductivity values. hConductivity will be used for
   layer type 0 or 2  to calculate the transmissivity along rows as :math:`TRAN = hConductivity * thickness`, for layer type 1 or 3 as :math:`HY = hConductivity`

vConductivity
   is the name of a spatial, scalar PCRaster map containing the vertical conductivity values. VCOND for a layer l will be calculated for each cell either as :math:`VCOND=\frac{1}{\frac{0.5 thick(l)}{vConductivity(l)} + \frac{0.5thick(l+1)}{vConductivity(l+1)}}` or as :math:`VCOND=\frac{1}{\frac{0.5thick(l)}{vConductivity(l)} + \frac{thick(cb)}{vConductivity(cb)} + \frac{0.5thick(l+1)}{vConductivity(l+1)} + }` when there is a Quasi-3D confining bed cb.

LAYER
   is the layer number the map values will be assigned to.

Annotation: For calculations in Modflow only vertical hydraulic conductivity values for quasi-3D confining beds are used. Horizontal conductivity values must also be specified for those layer due to technical reasons.

Transient simulations
~~~~~~~~~~~~~~~~~~~~~
If the simulation state is set to transient the specification of the storage coefficients is required. The operation is

.. code-block:: c

   res = mf::setStorage(primary, secondary, LAYER);

where

primary
   is the name of a spatial, scalar PCRaster map containing the Sf1 values.

secondary
   is the name of a spatial, scalar PCRaster map containing the Sf2 values and

LAYER
   is the layer number the map values will be assigned to.

Wetting capability
~~~~~~~~~~~~~~~~~~
To activate the wetting capability the following operations must be used. The non-spatial parameters are set with

.. code-block:: c

   res = mf::setWettingParameter(WETFCT, IWETIT, IHDWET);

where

WETFCT
   is the factor that is included in the calculation of the head that is initially established at a cell when that cell is converted from dry to wet;

IWETIT
   is the iteration interval for attempting to wet cells. Wetting is attempted every IWETIT iteration. This applies to outer iterations and not inner iterations. If IWETIT is 0, the value is changed to 1;

IHDWET
   flag that determines which equation is used to define the initial head at cells that become wet (0: :math:`h = BOT + WETFCT(h_n - BOT)`, 1: :math:`h = BOT + WETFCT(THRESH)`).

The wetting threshold and flag values are specified with

.. code-block:: c

   res = mf::setWetting(map, LAYER);

where

map
   is the name of a spatial, scalar PCRaster map holding WETDRY (WETDRY < 0, only the cell below a dry cell can cause the cell to become wet; WETDRY > 0, the cell below a dry cell and the four horizontally adjacent cells can cause a cell to become wet; WETDRY is 0, the cell cannot be wetted; absolute value of WETDRY is the wetting threshold) and

LAYER
   is the layer number the map values will be assigned to.

Optional operations
~~~~~~~~~~~~~~~~~~~
The head value that is assigned to cells that are converted to dry during a simulation (HDRY) can be specified with

.. code-block:: c

   res = mf::setDryHead(VALUE);

where

VALUE
   is the scalar, non-spatial head value.

If this operation is not used the value will be set to a default value of -999.9.

The variable containing the horizontal anisotropy factor (TPRY) can be specified with

.. code-block:: c

   res = mf::setHorizontalAnisotropy(VALUE, LAYER);

where

VALUE
   is the scalar, non-spatial horizontal anisotropy value.

LAYER
   is the layer number the value will be assigned to.

If this operation is not used for a specific layer the value will be set to a default value of 1.0 (isotropic conditions).
