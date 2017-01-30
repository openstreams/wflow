.. _quickStart:

***********
Quick start
***********

Displaying attributes
=====================
An attribute that is stored in a single file can be displayed by simply specifying the file name on the command line (see section :ref:`programOptions`). Basic commands for maps and time series are:

::

  # Show a map.
  aguila dem.map

.. image:: Images/dem_Map.png

::

  # Show time series.
  aguila rain.tss

.. image:: Images/rain_TimeGraph.png

::

  # Show some maps and a time series.
  aguila dem.map soil.map rain.tss

  # Show a map in 2.5D.
  aguila --drapeView dem.map

Each attribute is displayed in a separate window, here called a view (see section :ref:`views`). To display multiple attributes in a single view insert a ``+``-symbol between the names, with spaces surrounding the ``+``-symbol.

..
  KDJ: This actually does not work!
  # Multiple data items as drape and buildg.map only in the value matrix
  aguila --drapeView dem.map + soil.map --valueOnly buildg.map

::

  # Multiple data items as drape.
  aguila --drapeView dem.map + soil.map

.. image:: Images/dem+soil_Drape.png

Spatio-temporal attributes stored in the PCRaster raster format, is not stored in a single file but as a sequence of numbered rasters. This is true for most raster formats. To show such a spatio-temporal data set, two notations are allowed:

::

  # Old style: show discharge and rainfall map stack.
  aguila dis00000.001+36 rainfall.001+36

  # New style: show discharge and rainfall map stack.
  aguila --timesteps [1,36] dis rainfall

That is all information you need to display spatio-temporal raster attributes and timeseries. This manual will further instruct you how to create more sophisticated and combined views of data.

Required information
====================
There are three types of information Aguila has to know about before it can visualise your data: 1) the name(s) of the attributes you want to visualise, 2) the dimensions of the attributes in which it is defined (see section :ref:`dimensions`), and 3) the views you want to use.

Of course, the name of the data is required. Telling Aguila about the dimensions of the data may be required or not, depending on the data format (see :ref:`dataFormats`). For each type of data there is a default view which may be fine or has to be chosen otherwise by configuring specific start up options. See :ref:`programOptions` for information about using program options to pass information to Aguila.

Once Aguila has started visualising data it maintains links between data properties and the views, and between the views. For example, if two spatial data sets are visualised with 2D map views and in one of them the map shown is zoomed into, the other view will also zoom into the same area.

.. _programOptions:

Program options
===============
Aguila can be configured by passing it program options on the command line. For example, the command

::

  aguila dem.map

starts Aguila with the name of a dataset to visualise. If Aguila can read the data it will show it in the default view. Type

::

  aguila --help

for a list of all command line options supported.

Since typing long commands with a lot of options can be tedious, some program options can also be read from a configuration file. For example, the command

::

  aguila --config settings.cfg

tells Aguila to look for program options in the file settings.cfg.

Aguila can also be started from an XML specification file, which allows for more fine grained control over the Aguila configuration::

  aguila --xml settings.xml

See section :ref:`xmlStartupConfiguration` for details.

Command line options and configuration file options can be used at the same time, so

::

  aguila -f settings.cfg dem.map

starts Aguila with options on the command line and in a configuration file. Some options will be combined when given more than once and others will be overridden (this is documented in the tables below or obvious). When an option is overridden because it is given both on the command line and in a configuration file, the one on the command line is given precedence. One way to use configuration files is to put common options in them which can be used together with the command line (see also the examples_ below).

Some options take a range or a set of values as an argument. The syntax for a *range* of values is ``[first, last, step]`` where ``first`` is the first value, ``last`` is the last value and ``step`` is the interval between individual values between first and last. For example, the range ``[1,6,2]`` consists of the values ``1``, ``3``, ``5`` (note that ``6`` is not used in this case because ``5 + step`` equals ``7``, which is considered outside the range). The step is optional and the default value is dependent on the kind of range values. The range of timesteps ``[1,6]`` consists of the values ``1``, ``2``, ``3``, ``4``, ``5``, ``6``. The syntax for a *set* of values is ``{value1, value2, ..., valuen}``.

.. warning::

  The range and set notations can lead to suprising effects when used in a Unix (or Cygwin) shell (eg: bash). If you use such a shell, you must escape the ``{``, ``}``, ``[`` and ``]`` characters, because these have special meaning to the shell interpreter. Or you may quote the whole argument::

    aguila --scenarios \{a,b,c\} concentration
    aguila --scenarios "{a, b, c}" concentration

  When quoting the argument, you can use spaces between the values, otherwise you cannot.

Some options take the name of a data set as an argument. Aguila supports certain naming conventions which depend on the format the data is stored in. For example, a table stored in an ASCII column file is named by its filename, but the same table stored by a database management system might be named as ``myname(mypasswd)@mydbmsserver:mydatabase/mytable`` or just ``mydatabase/mytable``. For more information about these naming conventions see section :ref:`dataFormats`.

Command line options which can also occur in a configuration file must be named by the long option name, without the leading double dash. So, while ``-n`` and ``--scenarios`` are equivalent when used on the command line, only ``scenarios`` can be used in a configuration file. Furthermore, the value of an option given on the command line is put immediately after the option name or character, optionally separated by white space. In a configuration file the option name and value are separated by an equals sign, optionally surrounded by white space. See also section `Examples`_.

.. table:: Command line options

  ================================= ============================================
  Option                            Description
  ================================= ============================================
  -f [ --config ] arg               Read options from the configuration file
                                    named ``arg``.
  -x [ --xml ] arg                  Read options from the xml file named
                                    ``arg``, see section
                                    :ref:`xmlStartupConfiguration`.
  -h [ --help ]                     Show the command synopsis and exit.
  -l [ --lock ] arg                 Create a lock file named ``arg``. If the
                                    file does not already exists it
                                    is created. When Aguila exits the
                                    file is deleted again. This can be
                                    useful if Aguila is started by another
                                    application which wants to be able to
                                    check whether Aguila is still running.
  --license                         Show the software license and exit.
  -v [ --version ]                  Show the software version and exit.
  -2 [ --mapView ] arg              Show attribute named ``arg`` in a 2D map
                                    view.  ``arg`` can contain the names
                                    of more than one attribute. When the
                                    names of two attributes are separated
                                    by ``whitespace + whitespace``,
                                    they are stacked on top of each
                                    other in the same view. Otherwise
                                    each attribute is visualised in its
                                    own view.
  -3 [ --drapeView ] arg            Show attribute named ``arg`` in a 3D map
                                    view. ``arg`` can contain the names
                                    of more than one attribute. When the
                                    names of two attributes are separated
                                    by ``whitespace + whitespace``,
                                    they are stacked on top of each
                                    other in the same view. Otherwise
                                    each attribute is visualised in its
                                    own view. The first attribute is
                                    used for the height values. This
                                    attribute must contain scalar values.
  -t [ --timeGraphView ] arg        Show attribute named ``arg`` in a time graph
                                    view. ``arg`` can contain the names of
                                    more than one attribute. The optional
                                    selection specification in the attribute
                                    name must contain 2 column numbers of
                                    which the first one is regarded as the
                                    time step column and the second as the
                                    attribute column.
  -p [ --probabilityGraphView ] arg Show attribute named ``arg`` in a
                                    probability graph view. ``arg``
                                    can contain the names of more than
                                    one attribute.
  --valueOnly arg                   Show attribute named ``arg`` only
                                    in the value cursor matrix. ``arg``
                                    can contain the names of more than
                                    one data set.
  ================================= ============================================

.. table:: Command line and configuration file options

  ============================== ===============================================
  Option                         Description
  ============================== ===============================================
  -n [ --scenarios ] arg         Configures the scenario dimension using the set
                                 of scenarios in ``arg``. Multiple scenarios
                                 options are merged into one scenario
                                 dimension.
  -s [ --timesteps ] arg         Configures the time dimension using the range
                                 or set of time steps in ``arg``. Time steps
                                 must be larger than ``0``. Multiple time steps
                                 options are merged into one time dimension.
  -q [ --quantiles] arg          Configures the cumulative probability dimension
                                 using the range or set of quantiles in
                                 ``arg``. Quantiles must be larger than ``0``
                                 and smaller than ``1``. Multiple quantiles
                                 options are merged into one cumulative
                                 probability dimension.
  --cursorValueMonitorFile arg   Tells Aguila to append an ``aguilaCursorValue``
                                 element to the value monitor
                                 file named ``arg`` each time ``Save`` is
                                 pressed in the Cursor Value Window. On
                                 start up, ``arg`` is created with 0
                                 ``aguilaCursorValue`` sub-elements. The file
                                 is written in XML conforming to the Aguila
                                 XML Schema.
  -m [ --multi ] arg             When visualising scenarios of a spatial
                                 attribute, this option can be used to
                                 put all scenario's side by side in one 2D
                                 map view. ``arg`` should be formatted as
                                 ``<number or rows>x<number of columns>``, eg:
                                 ``2x3``.
  ============================== ===============================================

Examples
========
These examples assume ``dem``, ``ldd`` and ``erosion`` are valid names of raster attributes and ``discharge`` is a valid name of a time series attribute. These attributes are presented in a format supported by Aguila. Note that attributes in different formats can be combined.

2D raster
---------
Visualise a raster in 2D. Default view for rasters is 2D map.

::

  aguila dem

See also section :ref:`mapView`.

2D rasters on top of each other
-------------------------------
Stack rasters on each other. Spaces around the ``+``-sign.

::

  aguila dem + ldd

See also section :ref:`mapView`.

2.5D raster draped
------------------
2.5D is also possible.

::

  aguila -3 dem + ldd

See also section :ref:`drapeView`.

2D raster stack
---------------
Raster attribute might be temporal. ``dem00000.001+100`` is deprecated. Separate the name of the dataset from the dimension information.

::

  aguila --timesteps [1,100] dem

See also sections :ref:`mapView`, :ref:`timeGraphView`.
 
Time series
-----------
Visualise all time series in one time series plot.

::

  aguila discharge

See also section :ref:`timeGraphView`.

Time series selection
---------------------
Select some time series (the fifth and seventh columns) from discharge.

::

  aguila discharge{1,5} + discharge{1,7}

See also section :ref:`timeGraphView`.

2D raster draped and time series
--------------------------------
Combine rasters and time series data.

::

  aguila dem + ldd discharge

See also sections :ref:`mapView`,  :ref:`timeGraphView`.

Scenarios of temporal quantiles
-------------------------------
Select some more dimensions. For each scenario one 2D map view of the median value of ``erosion``.

::

  aguila --scenarios {simple,complex} --quantiles [0.01,0.99] --timesteps [1,100] erosion

This example assumes that two erosion models where created, a simple one and a more complex one. Each of these models was used in a Monte Carlo analysis resulting in a distribution of erosion outcomes. For each timestep the distribution in erosion outcomes was summarised by a range of percentiles, for example the ``0.01``, ``0.05``, ``0.1``, ``0.25``, ``0.5``, ``0.75``, ``0.9``, ``0.95``, ``0.99`` percentiles. Aguila starts by presenting the median value of ``erosion`` for each scenario as a 2D map view. Given this attribute it is possible to get a time series plot and a cumulative distribution plot for each location. All views can be animated. Aguila interpolates the data for percentiles that are not provided.

See also sections :ref:`mapView`,  :ref:`timeGraphView`, :ref:`cumulativeProbabilitiesView`.

.. _environmentVariables:

Environment variables
=====================
Aguila has support for many data set formats. While searching for data, Aguila tries each format driver in turn to see whether it is able to read the data. If you only use a fixed set of formats for your data sets, you can decrease startup time by limiting the number of potential data set formats Aguila considers. For this, define the environment variable called ``PCRASTER_DAL_FORMATS``. Its contents should be a comma separated list of format names. Most of them correspond to the names used in the GDAL and OGR data I/O libraries used by Aguila (:ref:`dataSetTypes`).

Example::

  $ export PCRASTER_DAL_FORMATS="CSF, HDF4, WCS, ESRI Shapefile"

.. note::

  Although GDAL comes with a PCRaster driver built in, Aguila makes use of its own PCRaster driver and skips GDAL's one. The name of the PCRaster raster format driver is ``CSF``.

Unknown format names are skipped. If none of the format names are known, Aguila will use all formats it knows of.

.. warning::

  Whenever Aguila is unable to read your attributes, make sure the setting of PCRASTER_DAL_FORMATS is not excluding a format driver required for reading the attribute.

