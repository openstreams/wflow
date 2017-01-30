Differences between PCRcalc and PCRaster Python
-----------------------------------------------

..
   CW on first read I thought where are the examples?
   Also most of what is in Syntax difference could be included
   when discussing the Additional operations.

The majority of PCRcalc operations can be easily adjusted to PCRaster Python operations. Nevertheless, both modelling languages differ in a few areas. The following sections explain these differences and give sample applications. The differences in the operators (e.g +,-,...) are already listed in :ref:`operators`.

Setting geographical attributes of a model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To use the PCRaster operations you need to define the model area and resolution. While this is done in PCRcalc with the ``areamap`` section, PCRasterPython does not make use of sections. Apply the ``setclone`` operation instead:

.. code-block:: python

   setclone(map)

Sets the clone.

map
   Filename of clone map.

Thus, the PCRcalc areamap section

.. code-block:: c

   areamap clone.map;


would become in a PCRaster Python script

.. code-block:: python

   setclone("clone.map")

.. note::

   The first operation in a script must be the ``setclone`` operation.

   Like in PCRcalc scripts it is not possible to mix maps with different extents, i.e. varying cellsize, number of rows or columns.

Reading input maps
^^^^^^^^^^^^^^^^^^
An input map is identified by a Python string containing the pathname to the disk dataset. The function ``readmap(pathname)`` reads a map.

* PCRaster operations with a function syntax can accept the pathname.

  .. code-block:: python

     Dem = readmap("dem.map")
     Slope1 = slope(Dem)
     # or direct by pathname
     Slope2 = slope("dem.map")

* PCRaster operations with an operator syntax do not accept a pathname.

  .. code-block:: python

     # will yield "dem.mapdem.mapdem.map"
     stringTimes3 = "dem.map" * 3
     # will yield a Spatial object
     DemTime3 = readmap("dem.map") * 3

* When using the PCRaster Python Framework one usually uses the class method ``readmap(self, filename)`` instead of the function ``readmap(pathname)``.

  .. code-block:: python

     # ... denotes removed fragments
     ...
     class RunoffModel(object):
     ...
       def initial(self):
         # read static data
         self.dem = self.readmap("dem.map")
       def dynamic(self):
         # reads rain0000.001, rain0000.002, etc.
         self.rainFall = self.readmap("rain")

     dynModel = DynamicFramework(RunoffModel, 50)
     dynModel.run()

Writing output maps
^^^^^^^^^^^^^^^^^^^
To store maps use:

.. code-block:: python

   report(map, filename)

Writes a map to a file.

map
   Map you want to write.

filename
   Filename to use.

Note the differences between PCRcalc and Python.

.. code-block:: c

   # PCRcalc script code
   binding
       gradient = output.map;
       dem = dem.map;
   initial
       report gradient = slope(dem);

.. code-block:: python

   # Python code
   gradient = slope("dem.map")
   report(gradient, "gradient.map")

The ``report`` operation only writes spatial data to disk. For writing non spatial data as timeseries in the DynamicFramework use the TimeoutTimeseries construct.

Accessing cell values of a map
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
With PCRaster Python you can query cell values from a map:

.. code-block:: python

   cellvalue(map, row, col)

Returns a cell value from a map.

map
   Map you want to query.

row
   Row index of a cell in the map. Indices range from [1, number-of-rows].

col
   Col index of a cell in the map. Indices range from [1, number-of-cols].

Returns a tuple with two elements: the first is the cell value, the second is a boolean value which shows whether the first element, is valid or not. If the second element is False, the cell contains a missing value.

.. code-block:: python

   cellvalue(map, index)

Returns a cell value from a map.

map
   Map you want to query.

index
   Linear index of a cell in the map. Indices range from [1, number-of-cells].

Returns a tuple with two elements: the first is the cell value, the second is a boolean value which shows whether the first element, is valid or not. If the second element is False, the cell contains a missing value.

Global options
^^^^^^^^^^^^^^
To set the global options of a script use

.. code-block:: python

   setglobaloption(name)

Sets the global option name. The option name must not contain the leading dashes as used on the command line of pcrcalc.

.. code-block:: python

   # Python example:
   setglobaloption("unitcell")

.. code-block:: c

   # The pcrcalc equivalent:
   pcrcalc --unitcell -f model.mod

Writing time series
^^^^^^^^^^^^^^^^^^^
PCRaster Python does not provide a ``timeoutput`` operation like PCRcalc does. Instead, a separate class in the modelling framework is handling the PCRcalc timeoutput style timeseries.

Therefore, the PCRcalc

.. code-block:: c

   binding
       outlet=samples.map;

   dynamic
       ...
       Q = ...
       report runoff.tss = timeoutput(outlet , Q);

will become in your Python model script

.. code-block:: python

   from pcraster.framework import *

   class RunoffModel(object):

       def initial(self):
           ...
           self.runoffTss=TimeoutputTimeseries("runoff", self, "samples.map",
               noHeader=False)

       def dynamic(self):
           ...
           runoff = ...
           self.runoffTss.sample(runoff)

   myModel = RunoffModel("clone.map")
   dynModelFw = DynamicFramework(myModel, endTimeStep=28, firstTimestep=1)
   dynModelFw.run()

In the ``initial`` section of the model class you create a member variable ``self.runoffTss`` holding the ``TimeoutputTimeseries`` object. The output is written to the file ``runoff.tss`` (for the DynamicFramework in the current working directory, for the MonteCarloFramework it will store the file into the corresponding sample subdirectories). ``samples.map`` is either a boolean, nominal or ordinal map holding the sample locations. By default the header section is written to the timeseries file.

In the ``dynamic`` section ``self.runoffTss.sample(runOff)`` samples the values of the expression (here ``runoff``) at the given locations for the current timestep. Note that for sequenced calls of ``sample()`` the values of the last call are sampled.

See also the example script ``runoff.py`` in the deterministic subdirectory of the workspace/framework directory.

Converting to and from NumPy arrays
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. note::

   The conversion functions are no longer provided by default by importing the PCRaster module, instead use: ``from pcraster.numpy import *``.

You can convert PCRaster maps to NumPy arrays and vice versa with the following conversion functions:

.. code-block:: python

   numpy2pcr(dataType, array, mv)

Converts entities from NumPy to PCRaster.

dataType
   Either Boolean, Nominal, Ordinal, Scalar, Directional or Ldd.
array
   Array you want to convert.
mv
   Value that identifies a missing value in array.

Returns a map.

.. code-block:: python

   pcr2numpy(map, mv)

Converts entities from PCRaster to NumPy.

map
   Map you want to convert.

mv
   Value to use in result array cells as a missing value.

Returns an array.

The following example script demonstrates the usage of the conversion functions:

.. code-block:: python

   import numpy
   from pcraster import *
   from pcraster.numpy import *

   setclone("clone.map")

   a = numpy.array([[12,5,21],[9,7,3],[20,8,2],[5,6,-3]])

   # conver a to a PCRaster Python map
   # with the value 20 indicating a 'missing value' cell
   n2p = numpy2pcr(Nominal, a, 20)
   print "cellvalue:", cellvalue(n2p, 2, 3)[0]

   # write it to a PCRaster map
   report(n2p, "n2p.map")
   # read the PCRaster map back
   p2n = readmap("n2p.map")
   # print it as a numpy array
   # missing value is replaced by 9999
   print pcr2numpy(p2n, 9999)

Execution of the script will result in the following output:

.. code-block:: python

   cellvalue: 3
   [[  12    5   21]
    [   9    7    3]
    [9999    8    2]
    [   5    6   -3]]

PCRcalc operations returning multiple results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In PCRcalc some operations can be combined and return two result values like:

.. code-block:: c

   # PCRcalc
   state,flux = accufractionstate,accufractionflux(ldd.map,material.map,0.5);

This syntax is not supported in PCRaster Python scripts. Use two separate operations instead, although this will double the execution time:

.. code-block:: python

   # Python
   state = accufractionstate("ldd.map", "material.map", 0.5)
   flux = accufractionflux("ldd.map", "material.map", 0.5)

Boolean operations
^^^^^^^^^^^^^^^^^^
Do not use PCRaster objects in context of Python boolean operations.

   "In the context of Boolean operations, and also when expressions are used by control flow statements, the following values are interpreted as false: None, numeric zero of all types, empty sequences (strings, tuples and lists), and empty mappings (dictionaries). All other values are interpreted as true."

   -- Python Reference Manual [http://docs.python.org/ref].

This means that PCRaster objects will always be interpreted as true when used in the above mentioned cases.

When the PCRaster (field) objects are used in a Python boolean context you will recieve the error::

   The truth value for PCRaster spatial data types is ambiguous.
   See the section Boolean operations in the PCRaster Python manual.

Most likely you used one of the following Python constructs:

.. code-block:: python

   booleanMap = readmap("dump.map")
   booleanMap2 = readmap("household.map")

   if booleanMap:
     selection = map1
   # instead use
   selection = ifthen(booleanMap, map1)

   if booleanMap:
     something = map1
   else:
     something = map2
   # instead use
   something = ifthenelse(booleanMap, map1, map2)

   dumpAndHouseHold = booleanMap and booleanMap2
   # instead use
   dumpAndHouseHold = booleanMap & booleanMap2

   dumpOrHouseHold = booleanMap or booleanMap2
   # instead use
   dumpOrHouseHold = booleanMap | booleanMap2

   dumpXorHouseHold = booleanMap xor booleanMap2
   # instead use
   dumpXorHouseHold = booleanMap ^ booleanMap2

   notDump = not booleanMap
   # instead use
   notDump = ~ booleanMap
