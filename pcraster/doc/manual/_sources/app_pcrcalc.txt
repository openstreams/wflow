
.. index:: pcrcalc; usage

.. _pcrcalc:

*******
pcrcalc
*******
.. topic:: pcrcalc

   Execute a map algebra statement or script file

::

  pcrcalc [options] MapAlgebraStatement

::

  pcrcalc [options] -f scriptFile


+-------+-------------------------------------------------------------------------------------+
|Options|Description                                                                          |
+=======+=====================================================================================+
|-s #   |set seed (integer > 0) for random generator                                          |
|       |Default is based on current time.                                                    |
+-------+-------------------------------------------------------------------------------------+
|-b f   |:ref:`Overrule script bindings <externalBindingFile>`                                |
+-------+-------------------------------------------------------------------------------------+
|-1     |:ref:`Update timeseries files at end of each timestep <pcrcalcoption1>`.             |
+-------+-------------------------------------------------------------------------------------+
|-r f   |Set run directory                                                                    |
+-------+-------------------------------------------------------------------------------------+
|-d f   |:ref:`Debug mode <debugmode>`, check MV creation on assignment                       |
|       |comparing against clone or areamap boolean                                           |
|       |mask.                                                                                |
+-------+-------------------------------------------------------------------------------------+
|-c     |Strict Case significant filename check                                               |
|       |(Unix portability)                                                                   |
+-------+-------------------------------------------------------------------------------------+
|-p     |Print :ref:`profile information <pcrcalcprofile>`                                    |
+-------+-------------------------------------------------------------------------------------+
|-m     |:ref:`Optimize with areamap MV compression <mvcompression>`                          |
+-------+-------------------------------------------------------------------------------------+
|-l     |:ref:`Use less memory <uselessmemory>` but more temporary disk                       |
|       |storage                                                                              |
+-------+-------------------------------------------------------------------------------------+
| -t    | test :ref:`argument substitution <secsubstitution>`.                                |
+-------+-------------------------------------------------------------------------------------+


.. index:: binding; external binding file

.. index:: external binding file

.. _externalBindingFile:

Overrule script bindings
========================

Usage:

 | pcrcalc -f script.mod -b binding.txt

All bindings specified in an external binding file have
priority over binding section found in the script itself. As
such one can overwrite binding from the script.

The syntax is identical to the binding section of a model
script including the ability to add comments, but the exernal
binding file has also the the following features: 

- The semi-colon (;) is optional 
- multiple definitions are allowed, the last definition is used.

Example:

::

   # timestep (seconds)
   DtSec=86400;

   # number of 1st day in model run
   StartTime=182;

.. _pcrcalcoption1:

Saving timeseries at each timestep
==================================

To optimize file access times pcrcalc only writes data to timeseries every 128 timesteps and at the end of executing the model.  For inspection of timeseries output during a model run it might be easier to update the timeseries at each time step. This will incur a runtime penalty.

To update the timeseries files at end of each timestep use the option --1.

.. _debugmode:

Debug mode
==========
A frequent problem when writing scripts is that sometimes the results have
unwanted missing values due to invalid ranges on some functions, such as
division by zero or taking the square root of a negative value.
The --d option of pcrcalc can help in finding the places where this
occurs. But only if you follow the *convention* that at each creation point of a
map, the map must have non missing values at each location where the boolean
areamap (or clone map) has a true (1) value. For example, the script test.mod:

::

  areamap mask.map;
  timer 1 100 1;
  dynamic
   Result.map = (VarA.map*VarB.map)/SomeZero.map;

If the script is called like pcrcalc -d debug.map -f test.mod then at each
timestep when Result.map receives a new map computed, the new map is checked if
it has any missing values where the areamap (mask.map) is true (1). If so,
pcrcalc terminates immediately with an error message, and writes the map
debug.map. The :ref:`error message <pcrcalcerrormessage>` gives the exact script location
of the statement that fails the areamap check. The map debug.map may contains the following values:

MV
  where a missing value is in the areamap

0
  where a 0 value is in the areamap
1 
  where a 1 value is in the areamap and the checked computation does not have a missing value
2 
  where a 1 value is in the areamap and the checked computation does have a missing value

With this technique you find may find the points where additional work is
needed. In the example above, we know that the formula only applies to cases
where SomeZero.map is larger than 0 and other cases should have the value 0:

::

  Result.map = if (SomeZero.map gt 0 then (VarA.map*VarB.map)/SomeZero.map else 0);


Another solution is to cover the generated missing values:

::

   Result.map = cover((VarA.map*VarB.map)/SomeZero.map, 0);

This kind of debugging is of limited use if one applies operations that
generate missing values on purpose (such as if then without the else clause) in
this way:

::

    Temp.map = if(VarA.map eq VarB.map then VarC.map);
    Result.map = cover(Temp.map, 0);

If the -d option is used then the creation of Temp.map will create an error
message. Rewriting the operation fixes this:

::

     Temp.map = if(VarA.map eq VarB.map then VarC.map else 0);

.. index:: missing value compression
.. index:: optimization; missing value compression

.. _mvcompression:

Missing value compression
=========================

This optimisation decreases the computation time and RAM memory
by a fraction almost equal to the number of cells in the non
Missing Value area divided by the number cells in the
rectangular grid. This fraction is called the mvFraction. To
compute this fraction of a particular dataset, divide the
number of MV cells by the total number of cells of the areamap
of the dataset. For example:

::

   pcrcalc define.map=defined(area.map)
   table --unitcell define.map define.tbl
   type/cat defined.tbl
    0 1860
    1 3140

For this particular dataset the mvFraction is
1860/(1860+3140)=0.37. This means the model can run about 37%
faster and using 37% less memory with most type of models.

The areamap has the value true(1) on each cell where the model
is active. No computation will be done on other cells. This
will reduce the execution time and memory requirements roughly
by the percentage of non true(1) cells on the areamap for most
models. A slight overhead is added by the number of global
operating functions, e.g. spread, accuflux and reading and
writing maps. The number of global functions used in a model is
in most models small compared to the number of point
operations, e.g. \*,/,+,-. When in doubt one should measure the
execution time with and without the --m option.

To enable this feature, start pcrcalc with the --m option:

::

   pcrcalc --m -f model.mod

When the --m option is active the areamap model section of the
model (or the --clone setting) must define a map with a cell
value equal to 0 or MV for the cells that should should be
excluded.

Important differences
---------------------

The --m option will mask out indifferent of the type of
operation. The following command will yield different maps.
coverWithM.map will still have MV's outside the true-defined
area.map area, coverNoM.map will not.

::

   pcrcalc coverNoM.map = cover(area.map,0);
   pcrcalc --m --clone area.map coverWithM.map = cover(area.map,0);

Be careful with ldd maps. Masking out parts of an ldd map may
result in an unsound ldd, with flow paths broken at unexpected
places. If possible use the ldd map of the model as areamap. Or
repair the ldd explicitly using :ref:`lddmask`.

::

  binding
   InputLdd = ldd.map;
   Area     = area.map;
  areamap Area;
  initial
   Ldd = lddmask(InputLdd,Area);

.. index:: missing value compression
.. index:: optimization; use less memory

.. _uselessmemory:

Use less memory
===============

The --l option of pcrcalc will turn pcrcalc into a disk based
computing system. Normally pcrcalc assumes it can keep all data
used more than once in memory. This is not the case with
extremely complicated models and/or datasets with large maps
due to large memory demands. If the model runs out of memory,
try running with the --l option. This will increase the
computing time significant.

The --l option will create a temporary directory where it reads
and writes temporary data. The name of the directory starts
with the pcrcalcSwap (e.g. pcrcalcSwap1, pcrcalcSwap2). Under
normal circumstances this directory is removed after the model
run. The location of the directory is different under windows
and linux:

 * windows: the TMP or TEMP environment variable (see SDK GetTempPath() for details)
 * linux: the TMP environment variable or the current directory if not set

The --l option will generates error messages prior to model
execution if different model parameters are bound to the same
external symbol. For example, below A1 and A2 are bound to the
same map: a.map, D and d.map are also bound to the same map:
d.map.

::

  binding
   D  = d.map;
   A1 = a.map;
   A2 = a.map;
  initial
   B  = A1*A2*D*d.map;

.. index:: memory usage; profile
.. _pcrcalcprofile:

profiling and tuning
====================

Memory demands of a model can now be measured by supplying the
--p option to pcrcalc. This will print the maximum number of
bytes per cell (bpc) needed for the model. Minimizing the
number of spatial parameters, pre computed static parameters,
feedback parameters will decrease this number.

Calculations for the total memory demand of a model can be done
by multiplying the maximum bpc by the number of cells used in
the raster. Without the --m option the number of cells in the
raster is equal to the number of cells in the raster, with the
--m option it is equal the number of non Missing Value cells.
For the --m option 8 bytes should be added to the maximum bpc as
printed by the --p option. Some sample calculations:

 | Number of total   cells = 500000 (500 rows by 1000 columns)
 | Number of defined cells = 314000
 | maximum bpc (--p)        = 533

Memory demand

 | no --m: 500000 * 533     = 266500000 bytes =  254 Mb
 | with --m: 314000 * (533+8) = 169874000 bytes =  162 Mb

The memory demand computed is only the memory needed for the
model maps. Memory needed for the program executable,
timeseries and tables are not included. In most cases these
additional resources can be neglected or estimated as:

 * pcrcalc program: about 10 Mb.
 * timeseries and tables: The total file size.

..
   The total amount of memory pcrcalc can allocate for a model run
   depends on the Operating System:
     * linux and windows (32 bit): 2 Gb minus the resources
       already allocated by the Operating System itself and other
       processes.
     * windows (64 bit): 4 Gb.
     * linux (64 bit): ? Gb.
   If pcrcalc is compiled as a native 64 bit application on
   windows and/or linux the total amount of memory can reach about
   8200 Gb.
