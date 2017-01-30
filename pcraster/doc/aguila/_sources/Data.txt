.. _data:

****
Data
****

.. _dataSets:

Data sets
=========
Attributes are stored in data sets. Depending on the data set formats, data sets are stored in a single file in the local or in a remote filesystem, or in multiple files, or in tables in a database, or are part of a single file containing multiple data sets, for example. Whatever is needed to store the values of a *single* attribute comprises the data set.

Attribute values in a data set describe the variation of a single attribute (e.g.: height, concentration, temperature). Dimensions are used to structure data set values (see section :ref:`dimensions`). Unique combinations of such dimensions give rise to different data set types (see section :ref:`dataSetTypes`). And various data formats are used to store data set values (see section :ref:`dataFormats`).

.. _dimensions:

Dimensions
==========
Dimensions structure data set values. Attribute values (potentially) vary along dimensions. For example, a raster contains spatially varying attribute values. Similar, a time series contains temporally varying attribute values.

A data set lacking attribute values for a certain dimension is taken to be constant along that dimension. For example, a single raster contains values that vary in space, but not in time. Its values are constant in time because the data set contains no information about the temporal variation.

The next sections describe the dimension that are supported by Aguila.

.. _space:

Space
-----
Aguila supports visualisation of both raster and feature data. Spatial raster data supported by Aguila is formatted as a raster with cells of equal width and height.

Aguila offers 2D and 2.5D views for rasters and 2D views for features.

.. _time:

Time
----
Temporal data supported is data which is defined at discrete integer time steps. Sophisticated schemes to link these discrete time steps to real time are only possible by using the XML interface described in section :ref:`xmlStartupConfiguration`.

Aguila offers time series views for temporal data. Additionally, it adds support for iteration through time (animation) to every view which is able to visualise temporal data, for example the 2D map view.

.. _uncertainty:

Uncertainty
-----------
If, instead of fixed variables, the available data are random variables (e.g. resulting from a Monte Carlo analysis or statistical modelling), Aguila can view them when they are encoded as (cumulative) probability distribution functions. The encoding is done by supplying maps that have for each location (and optionally for each time step) the attribute value corresponding to a fixed cumulative probability level, a value in ``[0,1]``. The simplest example would be if all variables followed a uniform distribution, in which case we only need to specify the map with the minimum (distribution value ``0``) and maximum (distribution value ``1``); Aguila uses linear interpolation for values in between. More complex distributions functions can be specified with maps for a set of distribution values. E.g. for a normal distribution one could specify the maps corresponding to probability distribution values ``0.01``, ``0.05``, ``0.1``, ``0.25``, ``0.5``, ``0.75``, ``0.9``, ``0.95`` and ``0.99``. Aguila assumes then straight lines in between these values.

Distribution values do not need to be symmetric; if much emphasis is put on the upper tail, the lower tail may be discretized by very few values (like, say, ``0.01``, ``0.5``, ``0.9``, ``0.91``, ``0.92``, ``0.93``, ...).

The interface of Aguila assumes that all distribution values provided are part of a regular sequence. Both examples above are subsets of the sequence ``0.01``, ``0.02``, ``0.03``, ..., ``0.99``. Minimum, maximum, and step size need to be specified for a sequence. Aguila will read those values present, and will linearly interpolate inbetween them.

.. _scenarios:

Scenarios
---------
Different scenarios contain data values for the same attributes which should be visualised at the same time, and with identical topological constraints. When data is loaded from different scenarios, Aguila will create views for each scenario. Scenarios are explicitly supported by the :ref:`2DMultimapView`. This view shows spatial attributes from different scenarios in a single window.

.. note::

  When scenarios are used, we still speak of a single attribute.

.. _dataSetTypes:

Data set types
==============
As described above, Aguila supports data sets with various kinds of dimensions. Every unique combination of dimensions gives rise to a data set type.

For some types of data sets, there are convenient data formats available, for example raster formats. For other data set types, we had to come up with some conventions to be able to store the values. There is, for example, no convenient data format availabe for storing uncertain spatio-temporal values.

This section describes the various types of data sets, and the next section describes the formats for storing the attribute values and some naming conventions used. The grey nodes in the figure below show data set types are not implemented. They are drawn and mentioned in the text for completeness.

..
  TODO Add a label and refer to it from the text.
  TODO Add links from the nodes to the sections.
  TODO Color the non-implemented nodes differently and remove them from the text.

.. graphviz::

  digraph DataSetTypes {
    constantValues[
           label="Constant values",
           color="grey55",
           fontcolor="grey55"]
    uncertainValues[
           label="Uncertain values",
           color="grey55",
           fontcolor="grey55"]
    uncertainTemporalValues[
           label="Uncertain temporal values",
           color="grey55",
           fontcolor="grey55"]
    spatialValues[
           label="Spatial values"]
    temporalValues[
           label="Temporal values"]
    spatioTemporalValues[
           label="Spatio-temporal values"]
    uncertainSpatialValues[
           label="Uncertain spatial values"]
    uncertainSpatioTemporalValues[
           label="Uncertain spatio-temporal values"]

    constantValues -> spatialValues[
           label="+space"];
    constantValues -> temporalValues[
           label="+time"];
    constantValues -> uncertainValues[
           label="+uncertainty",
           color="grey55",
           fontcolor="grey55"];
    spatialValues -> spatioTemporalValues[
           label="+time"];
    spatialValues -> uncertainSpatialValues[
           label="+uncertainty"];
    spatioTemporalValues -> uncertainSpatioTemporalValues[
           label="+uncertainty"];
    temporalValues -> uncertainTemporalValues[
           label="+uncertainty",
           color="grey55",
           fontcolor="grey55"];
    temporalValues -> spatioTemporalValues[
           label="+space"];
  }

* `Constant values`_
* `Uncertain values`_
* `Temporal values`_
* `Uncertain temporal values`_
* `Spatial values`_
* `Uncertain spatial values`_
* `Spatio-temporal values`_
* `Uncertain spatio-temporal values`_
* `Scenario values`_

.. _constantValues:

Constant values
---------------
The simplest data set type is the constant. A constant is just a single value, its value does not vary along one of the data dimensions.

Examples:

* An empirical constant.
* A boundary value, not to be exceeded.

.. Note::

  Not implemented.

.. _uncertainValues:

Uncertain values
----------------
Uncertainty adds information about the probability distribution of values that may represent the attribute. Note that there is no spatial or temporal variation in this case.

.. Note::

  Not implemented.

.. _temporalValues:

Temporal values
---------------
Temporal values vary in time, but they have no spatial variation or even spatial reference.

Examples:

* Rainfall for a relatively small area and assumed te be constant for the whole area

Aguila has support for reading temporal values from many `table data formats`_.

..
  For temporal values the folowing fields must be present in the table:
  * ``date``
  * <attribute name>

.. _uncertainTemporalValues:

Uncertain temporal values
-------------------------
TODO

.. Note::

  Not implemented.

.. _spatialValues:

Spatial values
--------------
Spatial values vary in space, but not in time.

Examples:

* Raster with elevation values.
* Feature layer with population density per administrative area.

Raster, feature and vector data sets contain spatial values. An exception to this rule is that Aguila can also visualise only the geometry of a feature layer.

Vector attributes are attributes with a direction and a magnitude. They are currently raster based, meaning that two rasters with magnitudes in x- and y-directins are used to store the values.

Aguila has support for reading spatial data from many `raster data formats`_, `feature layer formats`_ and `vector data formats`.

.. _uncertainSpatialValues:

Uncertain spatial values
------------------------
Uncertain spatial values are `spatial values`_ for which the probability distribution is known.

Currently, Aguila supports reading quantiles of spatial values, for example a range of the ``0.01``, ``0.05``, ``0.10``, ``0.25``, ``0.50``, ``0.75``, ``0.90``, ``0.95``, ``0.99`` quantile levels.

.. _spatioTemporalValues:

Spatio-temporal values
----------------------
Spatio-temporal values are `spatial values`_ which also vary in time.

Currently, Aguila supports reading time steps of spatial values, for example a range of the steps 1 through 1000.

.. _uncertainSpatioTemporalValues:

Uncertain spatio-temporal values
--------------------------------
Uncertain spatio-temporal values are `spatio-temporal values`_ for which the probability distribution is known.

Currently, Aguila supports reading quantiles of spatial values, for example a range of the ``0.01``, ``0.05``, ``0.10``, ``0.25``, ``0.50``, ``0.75``, ``0.90``, ``0.95``, ``0.99`` quantile levels.

.. _scenarioValues:

Scenario values
---------------
TODO

.. _dataFormats:

Data formats
============
Exactly which data formats Aguila supports depends on how Aguila is compiled. Potentially, Aguila can support many formats. The distributed Aguila supports about 50 raster formats (eg: PCRaster raster format, Hdf4, GeoTIFF, GML, ESRI's binary grid), a dozen feature formats (eg: ESRI's Shapefile) and a few table formats (eg: text, ODBC, Sqlite). Exactly which formats are supported can be deduced by looking at the script__ used to build the support libraries [#]_.

__ http://pcraster.svn.sourceforge.net/viewvc/pcraster/devenv/trunk/scripts/buildPcrTeamExtern.sh?view=markup

The data format used to store the data set determines the way the database must be configured. A data format might support data with some dimensions, like raster data with two spatial dimensions, but lack support for other dimensions, like the time dimension. Support for missing dimensions can be added to formats using naming conventions for file names, or by storing this information in an attribute table, for example. Most raster formats do not support more dimensions than the two space dimensions, but using naming conventions, information about the time dimension can be added easily. This results, for example, in a name like dem_1.map for the first time step. If a raster format supports the time dimension, such a naming convention is not needed and the application can read the information about this dimension from the data source.

.. [#] Eventually Aguila will be able to list the formats built in.

.. _rasterDataFormats:

Raster data formats
-------------------
The raster formats supported are PCRasters CSF 2.0 raster file format and the formats supported by the `Geospatial Data Abstraction Library`__ (GDAL). Information about these formats, including free (conversion) software and manuals can be found at the PCRaster website and the GDAL website.

__ http://www.gdal.org

For formats with default file name extensions there is no need to supply an extension.

Aguila uses the value scale property to decide how to visualise an attribute. In the CSF raster file format this value scale is stored in the data set. In data sets read using GDAL, this information may or may not be available. If a meta data item named ``PCRASTER_VALUESCALE`` is available, it is used to determine the value scale. If it is not available, the `color interpretation`__ of the data set is used to determine the value scale. Lastly, the type of the attribute values is used to determine the value scale. Integral types result in a nominal value scale and floating point types result in a scalar value scale to be used.

__ http://www.gdal.org/gdal_datamodel.html

The folowing values for the ``PCRASTER_VALUESCALE`` meta data item are recognized:

- ``VS_BOOLEAN``: boolean attribute
- ``VS_NOMINAL``: nominal attribute
- ``VS_ORDINAL``: ordinal attribute
- ``VS_SCALAR``: scalar attribute
- ``VS_DIRECTION``: directional attribute
- ``VS_LDD``: ldd attribute

The folowing color interpretation values are recognized:

- ``GCI_GrayIndex``: scalar attribute

.. _rasterDataFormatsNamingConventions:

Naming conventions
^^^^^^^^^^^^^^^^^^
Raster formats store spatial data in seperate files. Each file contains spatial attribute values for one map. To add information about other dimensions Aguila supports the folowing naming conventions:

.. table:: Raster file format conventions

  ========================= ====================================================
  Dimension                 Convention
  ========================= ====================================================
  Scenarios                 Data for different scenarios must be stored in sub
                            directories. When a scenarios dimension
                            is configured Aguila searches its data in
                            directories named after the scenarios.
  Cumulative probabilities  The different quantile levels are reflected in the
                            filename. The filenaming rule is
                            ``name_<level>{.extension}``. The quantile level
                            must be a floating point value between ``0`` and
                            ``1.0``. The file name extension is optional.
  Time                      Time steps are reflected in the filename. The
                            filenaming rule is ``name_<step>{.extension}``. The
                            time step must be a positive integral value
                            greater than ``0``. The file name extension is
                            optional.
  ========================= ====================================================

This means that, for example, a raster data source with scenarios ``simple`` and ``complex``, and temporal quantile levels has rasters named ``simple/co2_1_0.001.map`` and ``complex/co2_100_0.5.map``.

The PCRaster convention for naming spatio-temporal raster data is also supported. PCRaster uses the 8.3 DOS convention where each member of a stack is named by its name, possibly 0's and a time step number. So, a PCRaster model might output ``dem00000.001``, ``dem00000.002``, etc. When the time dimension is set up right, the name (``dem`` in this case) can be used to name such a stack. Also supported is the (depricated) convention of naming a stack by the first member, folowed by a ``+``-sign and the last time step of the stack, ``dem00000.001+1000`` for example.

.. _featureLayerFormats:

Feature layer formats
---------------------
The feature layer formats supported are the formats supported by the `OGR Simple Feature Library`__ (OGR). Information about these formats, including free (conversion) software and manuals can be found at the OGR website.

__ http://www.gdal.org/ogr

Naming conventions
^^^^^^^^^^^^^^^^^^
Aguila visualises attributes, not whole data sources. Feature data sources contain one or more feature layers and each layer contains zero or more attributes. The folowing naming rule must be used to tell Aguila which attribute to visualize::

  datasource/layer{/attribute}

When the attribute is not provided, Aguila will show the geometry of the layer. For example, given a ESRI Shapefile with information about countries, the folowing command will draw the country borders::

  aguila countries.shp/countries

And the folowing command will draw each country's population::

  aguila countries.shp/countries/population

Known file name extensions are optional and can be left out.

.. note::

  In the case of the Shapefile feature format, directories can be used as data sources too. In that case the names of the individual shape files (without extension) must be used as layer names.

Each feature layer attribute contains information about the spatial variation of the attribute. To add information about other dimensions, Aguila supports the use of an external attribute table in a database. The attribute table must be named after the feature layer and contain a field called ``fid`` (which is short for feature id). Aguila will join the internal attribute table of the feature layer data source to the external table using the feature id fields present in both tables.

See :ref:`tableDataFormats` for more information about setting up the external attribute table.

.. note::

   Currently, only SQLite databases have been used for storing external attribute tables for feature layers. In principle, all other table formats should also be supported. Also, since the name of the database and table are the same as the datasource and feature layer, it is not possible to pass server name and/or user credentials.

.. _vectorDataFormats:

Vector data formats
-------------------
Two raster datasets are used to hold the magnitude values for the x- and y-direction of a vector attribute. Because of this, all formats described in section :ref:`rasterDataFormats` can be used to store vector magnitude values.

Naming conventions
^^^^^^^^^^^^^^^^^^
There is a simple naming scheme for vector datasets::

  name_(x|y){.extension}

So, the values for a vector attribute called ``flow`` are stored in two raster datasets called ``flow_x`` and ``flow_y``. Default extensions are optional.

Naming conventions for vector datasets with more dimensions are the same as those described in section :ref:`rasterDataFormatsNamingConventions`.

.. _tableDataFormats:

Table data formats
------------------
Table file formats supported are the text and geoEAS formatted files and formats supported by the `QtSql module`__.

__ http://doc.trolltech.com/qtsql.html

Naming conventions
^^^^^^^^^^^^^^^^^^^
Naming a table data set depends on the format used to store the table in. For file formats the name of the data set is the name of the file it is stored in. For tables served by database management systems different conventions apply: ``myname(mypasswd)@myserver:mydatabase/mytable``. Some of the elements of this naming conventions might be optional, depending on the configuration of your database server. The server might, for example, grant read access to everybody, which means the account information in the name is redundant: ``myserver:mydatabase/mytable``. The database server might be your production machine, which means the server name is redundant: ``mydatabase/mytable``.

.. Depending on the application, extra information about the columns to select from the table can be given in the name of the dataset. To select the first and fifth column from a dataset you can use: ``mytable{1,5}``. The selection specification consists of a list of column numbers, seperated by a comma and optional whitespace, and surrounded by curly braces.

The table must contain a field named after the attribute. The values of this field will be visualized.

.. TODO How about table formats without a header. Default take the first two columns? The default should make PCRaster's timeseries files useable. Get rid of the selection stuff in the above paragraph(?).

.. TODO pick info from mail. Here, only info for joining the attribute. Info about the table layout can be described in the next section.

Attribute variation for one or more dimensions must be stored using a primary key consisting of one or more of the folowing field names: ``scenario``, ``date``, ``quantile``. For each unique combination of these, an attribute value can be stored in the table.

.. note::

  In case the table is used as an external attribute table for a feature layer, the primary key must also include the feature id field (``fid``).

