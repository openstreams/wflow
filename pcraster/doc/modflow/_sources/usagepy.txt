Usage in Python
===============
This chapter describes the usage of the PCRasterModflow extension within Python scrips. Using the extension requires the specification of the geographic extents and the creation of an extension object:

.. code-block:: c

   from pcraster import *

   setclone("clone.map")
   mf = initialise(clone())

This will initialise the data structures used in the extension. All operations will operate on the object ``mf``.

The differences of Python operations to the PCRcalc syntax is minimal and therefor only brief explained: The Python extension operations have the same name and take the same arguments as the PCRcalc operations. Setting for example boundary values in PCRCalc is done by

.. code-block:: c

   res = mf::setBoundary(boundaryValues, LAYER);

and in Python by

.. code-block:: c

   mf.setBoundary(boundaryValues, LAYER)

where ``boundaryValues`` is either a string holding a filename to a spatial, scalar PCRaster map or a variable holding an object returned by the readmap operation.

See the Python example script for further operations.
