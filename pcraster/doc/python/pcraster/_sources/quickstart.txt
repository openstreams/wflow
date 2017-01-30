PCRaster Python
---------------
Quickstart
^^^^^^^^^^
The Python example below is included in the workspace/Demo directory of the PCRaster distribution and shows the general usage of PCRaster functions within a Python script. To run the script change to the demo directory and execute ``python example.py``. The syntax of the PCRaster Python functions equates the syntax of the PCRaster functions. The results of the functions are written to disk with the ``report`` function.

The following sections give a detailed introduction about how to build environmental models in the Python scripting language.

*Python script example.py demonstrating the use of PCRaster functions.*

.. literalinclude:: ../../../data/demo/example.py
   :language: python


Introduction
^^^^^^^^^^^^
You can start using PCRaster Python by importing the main module in your Python script:

.. code-block:: python

   # Python
   import pcraster

This will put all PCRaster functions into the pcraster namespace. Here is an example of how you can use the ``slope`` function to calculate the slope of a digital elevation model:

.. code-block:: python

   # Python
   import pcraster

   gradient = pcraster.slope("dem.map")
   pcraster.report(gradient, "gradient.map")

This is equivalent to the PCRcalc script

.. code-block:: c

   # PCRcalc
   report gradient.map = slope(dem.map);

After importing the ``pcraster`` module the ``slope`` function is called with a filename as its argument. The function will open the raster file, read its values and calculate the slope. The resulting raster is returned and assigned to the variable ``gradient``.

Both, PCRcalc and PCRaster Python operations use exactly the same algorithm. If you compare the Python and the PCRcalc code you will see a minimal difference. For information about ``slope`` function or any of the other functions of PCRaster you can look them up in the PCRaster reference manual.

..
   TODO Not entirely true: The difference with PCRaster is that you can use filenames as arguments to functions whereas in PCRaster you should bind the filenames to variable names in the binding section of your PCRaster script.

The resulting gradient calculated above can be used as input to another function like this:

.. code-block:: python

   # Python
   import pcraster

   gradient = pcraster.slope("dem.map")
   smoothGradient = pcraster.windowaverage(gradient, 3)
   pcraster.report(smoothGradient, "smoothGradient.map")

By combining functions like this environmental models can be created.

In the examples given so far we had to explicitly state the module which the PCRaster functions are part of (``pcraster``). We can do better by importing all symbols from the ``pcraster`` module into the current scope. This way our script will become shorter and easier to read:

.. code-block:: python

   # Python
   from pcraster import *

   gradient = slope("dem.map")
   smoothGradient = windowaverage(gradient, 3)
   report(smoothGradient, "smoothGradient.map")

.. note::

   From now on the ``import`` statement will be discarded from the examples.

.. warning::

   By importing all symbols from a module there is an increasing chance that symbols with the same name end op in the same namespace. Keep this in mind when, for example, the abs function seems to work but gives an unexpected result: which function got called, math.abs() or pcraster.abs()?

.. _operators:

Operators
^^^^^^^^^
Besides functions operators can be used in expressions. Most boolean, comparison and arithmetic operators of PCRcalc have corresponding operators defined in PCRaster Python. Here's an example of the use of an operator:

.. code-block:: python

   # Python
   gradient = slope("dem.map")
   steepSlopes = gradient > 20

Below is a table of all operators defined in PCRaster Python and their their equivalents in PCRcalc. For more information about operators see the PCRaster manual pages and the Python operator module.


*List of PCRaster Python Operators*

========================= ========================= ===================
PCRaster Python operators PCRaster Python functions PCRcalc equivalents
========================= ========================= ===================
a & b                                               a and b
a | b                                               a or b
a ^ b                                               a xor b
~a                                                  not a
a != b                    ne(a, b)                  a ne b, a != b
a == b                    eq(a, b)                  a eq b, a == b
a > b                     gt(a, b)                  a gt b, a > b
a >= b                    ge(a, b)                  a ge b, a >= b
a < b                     lt(a, b)                  a lt b, a < b
a <= b                    le(a, b)                  a le b, a <= b
a * b                     mul(a, b)                 a * b
a / b                     div(a, b)                 a / b, a div b
a // b                    floordiv(a, b)            a idiv b
a ** b                    pow(a, b)                 a ** b
a % b                     mod(a, b)                 a mod b
a + b                     add(a, b)                 a + b
a - b                     sub(a, b)                 a - b
-a                                                  -a
+a                                                  +a
========================= ========================= ===================
