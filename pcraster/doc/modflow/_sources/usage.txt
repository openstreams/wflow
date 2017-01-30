Usage in PCRCalc
================
This chapter describes the usage of the PCRasterModflow extension. The first step is to create one single extension object in the script's inital section:

.. code-block:: c

   initial
   ...
   # Construct Modflow extension object
   object mf = PCRasterModflow::initialise();


This will initialise the data structures used in the extension. All operations described in the package sections will refer to the object ``mf``.

.. warning::

   Operations specifying Modflow input data currently return a boolean map containing false values, this behaviour may change in the future. It is strongly recommended to keep the coding style as shown below (i.e. ``res = mf::operation();``) and advised not to use the result in oncomming operations.

The next step is the grid specification using the operations described in the DIS package. This must be done before any other package is defined. Afterwards packages can be defined in arbitrary order. For a Modflow simulation at least the DIS, BAS and BCF package must be specified.

The DIS, BAS, BCF and a solver package must be set in the initial section of a script. Stress packages (RIV, DRN, RCH and WEL) can be activated and modified in the dynamic section.

.. note::

   One timestep in PCRaster represents one stress period in Modflow.

   In the following non-spatial arguments are written capitalised.

.. toctree::
   :maxdepth: 1

   mfin
   mfout
