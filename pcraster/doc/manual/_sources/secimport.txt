

.. _secimport:

**********************************************************
How to Import or Export Data, Display Maps, Global Options
**********************************************************


.. _secimportintro:

Introduction
============

This Chapter contains information about how to in- and export data from PCRaster.
It discusses how to convert data from other
GIS to one of the kinds of data used in PCRaster (:ref:`secimportmapcrea`) and how to export data from PCRaster (:ref:`secimportmapexpo`).  Other database management operations are also described, such as how to cut or join together maps and how to change location attributes or the data type of maps (:ref:`secimportmaptype`). Furthermore it will be explained how to attach a legend to a map (:ref:`secimportmapleg`) and how to display PCRaster maps on the computer screen and print a map (:ref:`secimportmapdis`). 


The other sorts of data used in PCRaster (tables, time series
and point data column files) are ascii formatted. Unlike PCRaster
maps, these can easily be created, analyzed and edited with other
software packages.  :ref:`secimportlook`  and :ref:`secimporttime`  explain the creation of tables and time series, respectively. Creation and conversion of point data column files is described in :ref:`secimportpoint`.  The above mentioned sections do not describe the formats of these sort of data. The format of PCRaster maps (:ref:`secdatbasemap`),  tables (:ref:`secdatbaselookform`), time series (:ref:`secdatbasetimeform`) and point data column files (:ref:`secdatbasepointform`) are described succesively. 

:ref:`secimportopt` gives a list of :emphasis:`global
options` used in PCRaster. Global options are set only once (for instance at the beginning of a project) and will affect all operations to which they are relevant. :emphasis:`Local options` which are used each time a operation is done, must always be redefined. 



.. _secimportmap:

PCRaster maps: database management
==================================


.. _secimportmapcrea:

Creation of a PCRaster map; data import
---------------------------------------

.. index::
   single: clone map

If you start a project and want to use PCRaster maps for analyzing your data, first make an empty clone map to define the geographical and cartographical location attributes and set this map as global clone map with the global option :ref:`--clone <GOClone>`.  To make the clone map, use the operator :ref:`mapattr`; the location attributes can be entered using the menu given by mapattr.  This map can be used as a clone (mask map) for all data sets you want to import. A data type and its cell representation can also be attached to the clone map using mapattr. This data type and cell representation will by default be assigned to the data that are imported, but the type of the imported data can also be specified during the import. 

.. _LocAtChoice:

.. index::
   single: location attributes; choice of

The choice of geographical location attributes must be based upon the spatial characteristics of the data set you want to import to the PCRaster map format and upon whether x,y coordinates are attached to the data in the data set or not. There are two possibilities for data import:  

.. _ColomnInput:

.. index::
   single: point data column files; conversion to a map

1) import of point data with x,y coordinates using either a :emphasis:`column file in simplified Geo-EAS format` or :emphasis:`plain column file format` 

.. _AngleAdv:

.. index::
   single: angle

The data file you want to import contains x and y coordinates with data values. In this case the location attributes of the PCRaster map must be chosen in correspondence with the spatial distribution of the data given by the x and y coordinates. If the data are regularly spaced on a rectangular grid, you probably want location attributes that match the set of x and y coordinates of the data set. If the data are irregularly spaced we advise to choose a map size of the smallest bounding rectangle (or somewhat larger) that comprises the study area, as shown in :ref:`fig3.5`. In this Figure a rectangle has been chosen that is rotated with respect to the real coordinate system (positive angle). Quite often, it is better not to rotate the map, and to choose an unrotated (smallest) rectangle.  Rotation may result in a map that better fits the shape of the study area. But remember that rotation has an important effect which may be very clumsy: rotation will not only result in a rotated map but of course also in a grid of cells that is rotated with respect to the real coordinate system. As a result, cells that are in one row on the PCRaster map will not have the same y coordinates; the same holds for cells in one column and their x coordinates.   


Two sorts of column files with x and y coordinates may be imported:
a :emphasis:`column file in simplified Geo-EAS format` or a :emphasis:`plain column file format`. These formats are used in PCRaster for representation of point data, especially in the gstat module for geostatistical analysis. :ref:`secdatbasepointform` describes these formats. Point data are imported to PCRaster map format with the PCRaster operator :ref:`col2map <col2map>`. 


2) import of data without x and y coordinates (ascii formatted) 



The data do not contain x and y coordinates: the asciifile with your
data contains a sequence of cell values, without coordinates. In this
case, the number of rows and columns of the PCRaster map must correspond
exactly with the number of rows and columns of the file you want to import
or else the import of data will result in nonsense.  Among other data,
maps from the ARC/INFO or Genamap GIS packages are imported in this
way. It is done with the PCRaster operator :ref:`asc2map <asc2map>`. 



.. _secimportmapexpo:

Data export from a PCRaster map
-------------------------------

Data can be exported from a PCRaster map using one of the operators
:ref:`map2col <map2col>` and :ref:`map2asc <map2asc>`. 

.. _DataExpFrMap:

.. index::
   single: conversion; from map to point data column file

The operator map2col exports data to an ascii column file in simplified Geo-EAS format or a plain column file format. Both contain x,y coordinates and data values. These kind of data are also used in the PCRaster package for representation of point data, see for formats :ref:`secdatbasepointform`. The plain column file format can easily be imported in spreadsheet, database management or word processing programs. It is also used in the gstat module of PCRaster. 

.. _DataExpFrMAsc:

.. index::
   single: conversion; from map to ascii file without x,y

The operator map2asc exports to an ascii file which will contain only data values, without x and y coordinates.   This operator is used if you want to export data to the ARC/INFO package. 



.. _secimportmaptype:

Cutting or joining maps; changing geographical location attributes or data types
--------------------------------------------------------------------------------

The spatial characteristics of the spatial data in PCRaster map format can
be changed with both the operator :ref:`resample <resample>` and :ref:`mapattr <mapattr>`.  Note that these two operations have a totally different result:  

.. _ResampleAMap:

.. index::
   single: cutting a map

If you want to resample your data in a PCRaster map to a new map with different geographical location attributes use the operator resample. First create a new map with the location attributes you wish (this is done with mapattr). For instance it may have geographical location attributes that define an area that only partly covers the old map, with a somewhat smaller or larger cell size. Now the resample operator can be used to resample your old data to the new map: for each cell of the new map a new cell value is computed on basis of the cell values on the old map that come into the cell on the new map. 


The resample operator can also be used to join two maps together. The maps that are joined together will be resampled and may have different location attributes: for instance they may overlap or may not overlap or may have different cell sizes.   


The geographical location attributes of a map can be changed using
the operator mapattr. Using this operator will not result in resampling of the data: each cell of the new map will contain a value that corresponds with the value on that cell of the old map. For instance halving the cell width of a map that consists of 50 x 50 cells of width 10 m. results in a (smaller) map of 50 x 50 cells of width 5 m., containing values that are taken directly from the old map. So changing the location attributes with mapattr will result in a new location of your data with respect to the real world coordinate system. Maybe this sounds silly but you may want to change the geographical location attributes after a map has been made, especially if you made an error in the specification of the location attributes with mapattr. 

.. _ConDataType:

.. index::
   single: conversion; between data types

Conversion between data types is done using one of the conversion of data type operators (:ref:`boolean <boolean>`,  :ref:`nominal <nominal>`,  :ref:`ordinal <ordinal>`,  :ref:`scalar <scalar>`, :ref:`directional <directional>` and  :ref:`ldd <ldd>`. These operators change the data type of the input map to the data type that corresponds with the name of the operator. Conversion is only possible if it results in a map that has some meaning with the new data type attached to it. 



.. _secimportmapleg:

Attaching a legend to a PCRaster map
------------------------------------

The operator :ref:`legend <legend>` attaches a legend to a PCRaster map. 



.. _secimportmapdis:

Screen display, hard copy output of PCRaster maps
-------------------------------------------------
:emphasis:`hard copy output: not yet included in software`


Visualisation of PCRaster maps and time series is done with aguila. 



.. _secimportlook:

Tables: database management
===========================


.. _secimportlookintro:

Introduction
------------
:ref:`secimportlookcrea` describes how to create or edit a table. 



.. _secimportlookcrea:

Creating and editing tables
---------------------------

By default, the PCRaster package uses column tables. For defining
relations between two maps of boolean, nominal, ordinal or ldd
data type it is sometimes better to use matrices instead of column
tables. If you want to use the matrix tables, set the global option
:literal:`--matrixtable`.  This option can be set for one separate operation or as general global option. How this is done, will be described later on in this chapter (:ref:`secimportopt`). If a relationship is specified between more than two maps in the matrix, the relation can not be described by a matrix table: PCRaster will automatically use a column table. 


A table for the PCRaster operator 
:ref:`lookup <lookup>` for creating a new map on basis of a table) or an input table for the operator :ref:`table <table>` for counting the number of cells that match the key) can be made with a text editor or alternatively with spreadsheet or word processing programs by exporting your table as a file in ascii text format. Additionally, the table operator itself generates a correctly formatted table, which can be used (possibly with a few edits) as input table for the lookup operator. 



.. _secimporttime:

Time series: database management
================================


.. _secimporttimeintro:

Introduction
------------

The format of time series has been described in an earlier chapter (:ref:`secdatbasetimeform`). Underneath (:ref:`secimporttimecrea`) database management with time series is described. 



.. _secimporttimecrea:

Creating and editing time series
--------------------------------

If you have time series data in a spread sheet program, database
management program or a package for statistics you can create a
PCRaster time series file by exporting your data as ascii formatted
text. A PCRaster time series file must have the lay-out as described in
:ref:`secdatbasetimeform`.  We advice to export the time series data in such a way that the resulting ascii formatted file will have a lay-out that looks like the lay-out used in PCRaster. Minor changes (for instance adding a header) can be made with a text editor. 


Additionally, a timeseries file can be created by reporting a time series
file to the database during a run of a dynamic model (see 
:ref:`secseqinreport`). 


You can import a PCRaster time series file to an other software package by
importing the time series as ascii text.




.. _secimportpoint:

Point data column files: database management
============================================


.. _secimportpointintro:

Introduction
------------

The format of point data column files has been explained in the previous chapter (:ref:`secdatbasepointform`). Underneath (:ref:`secimportpointcrea`) the database management with point data column filess described. 



.. _secimportpointcrea:

Creating point data column files, conversion to/from PCRaster maps
------------------------------------------------------------------

If you have spatial data in a spread sheet program, database management
program or a package for statistics you can create a point data column
file by exporting your data as ascii formatted text. Before exporting,
we advise you to put the x and y coordinates in the first and second
column respectively. The third and following columns may contain
the data. After exporting as ascii text, you can check and edit the
ascii column file with a text editor. If you want to convert it to a
PCRaster map it must have the simplified Geo-EAS format or the plain
column file format, which is described in the previous chapter (:ref:`secdatbasepointform`). Conversion is done with the col2map operator (:ref:`secimportmapcrea`) about importing to a PCRaster map or see the operator :ref:`col2map <col2map>`. Interpolation of regular or irregular spaced point data to a PCRaster map can be done with the gstat module. 

.. _PoiCreFrMap:

.. index::
   single: point data column files; conversion from a map

Additionally, you can create a point data column file by exporting data from a PCRaster map, with the map2col operator. The point data column file will be in plain column file format or in simplified Geo-EAS format. See the :ref:`map2col operator <map2col>`. 


.. index::
   single: options; global versus local
.. index::
   single: global options

.. _secimportoptintro:
.. _secimportopt:

Global options
==============

All :ref:`applications <applications>` have options to define the exact, detailed functionality of the operation.
Two types of options are used in PCRaster. :emphasis:`Local options` begin with a single \- (dash) are 1 character long and do mean
some different depending on the application. Local options are documented in the application reference documentation.

`Global options` begin with a double \- (dash) are 1 word long do have the same meaning across all applications.

If you want to set a different global option for only one operation,you can specify a global option in the command line. This is done
in the same way as it is done for a local option, by typing the
:literal:`--globaloption` after pcrcalc. If a non-pcrcalc operator is applied (for instance table) the global option is typed after the operator. The general global option setting described above is overruled only for the execution of the operation that is given in the command line. 

In addition there are two other methods to set global options, who are explained next:

- using environment variable PCROPTIONS
- in the first comment line of a pcrcalc script

.. index::
   single: PCROPTIONS; environment variable

using environment variable PCROPTIONS
-------------------------------------

If you use the MS-Windows version of PCRaster, you can set general global
options by typing after the command prompt:

 | set PCROPTIONS = globaloption1 globaloption2 ...globaloptionN

where globaloption1 is one of the global options, which, unlike a local option, is preceded by :emphasis:`two` - characters. If you use the UNIX version of PCRaster, global options are set by typing after the UNIX prompt:

 | PCROPTIONS=globaloption1 globaloption2 ...globaloptionN
 | export PCROPTIONS

Note that in UNIX no spaces are typed on either side of the = sign. For instance:

.. parsed-literal::

 PCROPTIONS='--clone CloneStudyArea.map --lddin'; export PCROPTIONS 
 
After applying PCROPTIONS the global options have the setting as specified or, if they are not specified, the default values. This set of options is used until a new set of options is specified with PCROPTIONS. If you set the options again with PCROPTIONS, options which are not specified that time are always set to default, no old settings are taken. 

.. index::
   single: global options; within script

within script
-------------

Global options that apply to a script can be stored in the script itself instead of typing it on the command line when the script is called. The first two characters of the first line of the script should read #!. After #! the global options can be specified. For example

.. parsed-literal::

  #! --lddfill --radians
  ..etc.

.. note::

  The #! should exactly be the first positions on the first line. Leading (empty) lines or spaces will change the line to a normal comment, and no global options from that line are set. No error messages are printed in such a case.

The remainder of this discussion on global options in the script applies only to those who use pcrcalc from a UNIX-like shell. The syntax of this first line is modelled after the alternative interpreter syntax of UNIX shells. If the first 2 characters of this first line are #! then all words on that line that can be recognized as global options, are parsed as global options, overriding global options at the command line. For example, example.mod:

.. parsed-literal::

  #! /usr/local/bin/pcrcalc -F --lddfill
  ..etc.

To be fully compatible with the UNIX alternative interpreter notation, one needs to specify the -F flag on the first line if and only if one wants to put global options on the first line. Note that pcrcalc will always scan the first line of a script, even if the script is not executed by itself, in other words:

 | pcrcalc -f example.mod
 | is equal to
 | chmod +x example.mod
 | example.mod

In both cases the --lddfill option will be activated. Since pcrcalc scans the first line itself, these feature is also enabled on NON-UNIX OS's (MSDOS).

-F is only neccessary and can only be used if one wants to specify global options. If the goal is to only describe pcrcalc as the interpreter of the script then use -f. This is due to a limitation in the syntax of an alternative interpreter, only one option can be specified. For example, example2.mod:

.. parsed-literal::

  #! /usr/local/bin/pcrcalc -f
  ..etc.




.. _secimportoptover:

Overview of global options
--------------------------
.. _LocAtrOp:

.. index::
   single: clone map

.. _GOClone:

:emphasis:`global options related to location attributes:`

:literal:`--clone` :emphasis:`CloneMap`
   The CloneMap is a PCRaster map that must have the location attributes of the maps you want to use during a project. It may be an empty map made at the start of a project using the operator :ref:`mapattr <mapattr>`, see also :ref:`secimportmapcrea`. Alternatively you may specify as CloneMap an existing PCRraster map containing data.

:literal:`--unittrue` (default) or :literal:`--unitcell` 
   This option specifies the units used for the coordinates and sizes ofnthe cells. It is of importance to the operations that make calculationsnwith distances or areas in the map or to operations that import ornexport coordinates of cells.  nnDefault, with the option :literal:`--unittrue`, PCRaster uses true distances and the true coordinate system of the map. These are given during creation of your map or the clone map of your map with :ref:`mapattr <mapattr>` and :ref:`secdatbasemaphead` about the location attributes of a map. The cell length is defined by the real length of a cell. The x coordinates are real distance coordinates and increase from left to right, starting with the x coordinate at the left edge of your map; the y coordinates increase from top to bottom, starting with the y coordinate at the top edge of your map, or from bottom to top, starting with the y coordinate at the bottom edge of the map (depends on the projection you have chosen).  PCRaster uses a sort of matrix coordinates if the option :literal:`--unitcell` is set.  This option is seldom used. Both the real coordinate system and the real cell length of the maps (given with the mapattr operator) are totally ignored. In all operations a cell length of 1 is used. If coordinates are exported or imported, the top left corner of the map is assigned the x,y co-ordinates (0,0); x increases from left to right and y increases always from top to bottom. As a result the centre of the top left cell of a PCRaster map has always :literal:`--unitcell` x,y coordinates (0.5,0.5).

:literal:`--coorcentre` (default) or :literal:`--coorul` or :literal:`--coorlr` 
   This option gives the coordinate position that is linked to a cell. If coordinates of cells are exported, you can choose to export for each cell the x,y coordinates of the centre of the cell (:literal:`--coorcentre`, default), the upper left corner of the cell (:literal:`--coorul`) or the lower right corner of the cell (:literal:`--coorlr`).  If point data with x,y coordinates are imported to a PCRaster map the option determines the assignment of point values that have x,y coordinates exactly at the edges of cells. This is only relevant to the operator :ref:`col2map <col2map>`.

:emphasis:`global options for defining the sort of directional data type:`

:literal:`--degrees` (default) or :literal:`--radians`   
   Default, the program interprets directional data as degrees (domain [0,360>). If the option :literal:`--radians` is set, directional data are interpreted and displayed as radians (domain [0,2pi>).

:emphasis:`global options for defining cell representations:`

.. _GOrepres:

.. index::
   single: cell; representation, global options

These options define the representation of cell values used for storage and processing of data. The settings are only applied to the maps that are created: if you change the cell representation settings during a project, only the new maps you generate are assigned the new cell representation. The old maps created before you changed the settings will keep the 'old' cell representation.   

scalar and directional data type:  
:literal:`--single` (default) or :literal:`--double`

Default, the cell representation is :emphasis:`single real` (:literal:`--single`). If you set :literal:`--double` it will result in :emphasis:`double real` representation of cell values. See also scalar data type, :ref:`formscalar` and directional data type, :ref:`formdirectional`.

:ref:`nominal <formnominal>` and ordinal data type: :literal:`--small` (default) / :literal:`--large`

Default, nominal and ordinal data are represented by :emphasis:`smallninteger` cell representation. If you set :literal:`--large`, :emphasis:`large integer` is used. See also and ordinal data type, :ref:`formordinal`.

:emphasis:`global option for specifying the format of tables`  

:literal:`--columntable` (default) or :literal:`--matrixtable`

Default, the operators :ref:`lookup <lookup>`  and :ref:`table <table>`  use column tables (setting :literal:`--columntable`).  If the option :literal:`--matrixtable` is set, a matrix table is used. For a description of these formats, see  :ref:`secdatbaselookform`.

:emphasis:`global options related to generation of a local drain direction map`  

:literal:`--lddout` (default) or :literal:`--lddin`

This option determines whether small catchments at the edge of the map are or are not considered as pits.    See the operators :ref:`lddcreate <lddcreate>`  and :ref:`lddcreatedem <lddcreate>`.

:literal:`--lddfill` (default) or :literal:`--lddcut`

This option determines the way the elevation model is modified in the corenof pits. See the operator :ref:`lddcreatedem.  <lddcreatedem>`

:emphasis:`global option for specifying the neighbourhood of cells`  

:literal:`--diagonal` (default) or :literal:`--nondiagonal`
This options specifies the neighbourhood of cells in the clump operation, see the :ref:`clump operator <clump>`.

global option for defining the sort of message printed during execution of an operation  

:literal:`--noprogress` (default), :literal:`--progress` or :literal:`--nothing`

This option affects the message printed on the screen during execution of an operation.  Default, the name, copyright and version of the software is printed  and error messages if errors occur (:literal:`--noprogress`). The option :literal:`--progress` will result in additional information printed during execution, such as 'busy with row x', etc. If you set the option :literal:`--nothing`, nothing is printed on the screen, except error messages.

.. index:: noheader
.. index:: time series; noheader
 
.. _noheader:

global option for defining the format of timeseries output.
:literal:`--noheader`

Only use this option if the output format is
not required to be compatible with PCRaster software
components because the created timeseries can not be read back as 
a timeinput series nor can be vizualized by Aguila. 

The option will affect timeseries as such that:

  #. No header is written.

  #. :ref:`Selective reports <selectivereports>` for timeseries are enabled.

  #. no MV data is prepended when start time is larger than 1.

  #. Timeseries are only created at the end of the model run and not partial written when running the model. This means the
     -1 option is discarded.
