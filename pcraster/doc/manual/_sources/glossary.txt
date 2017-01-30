########
Glossary
########
.. glossary::
   :sorted:


angle
   Location attribute: the angle is the angle between the horizontal direction on the PCRaster map and the x axis of the real world coordinate system.

area
   A fundamental spatial unit of a PCRaster map consisting of one cell or a set of cells.

areamap section
   The section in a dynamic modelling script that defines the location attributes of the maps used in the model.

attribute
   Property of a geographic object or location.

boolean data type
   data type for attributes that only may have a value TRUE (cellvalue 1) or FALSE (cell value 0).

Cartographic Model
   A model that computes new attribute values from those attribute values already present. It represents one static change in the property of cells.

cartographic modelling script
   A script of a Cartographic Model.

cell
   The basic spatial element of a PCRaster map.

cell length
   Location attribute. The length of cells on a PCRaster map which is the same in horizontal and vertical direction. Measured in the distance unit of the real coordinate system.

cell property
   The combined information at one cell location stored in one PCRaster map or several PCRaster maps.

cell representation
   The method of storage used in the computer for cell values of PCRaster maps.

column
   A column of cells in vertical direction on a PCRaster map.

database
   A collection of interrelated information, usually stored on a harddisk. The PCRaster database includes data stored as binary PCRaster maps, and ascii formatted tables, time series and point data column files.

data type
   Data description attached to a PCRaster map. It defines the scale and the domain of the data stored in the map and as result the behaviour with respect to operations performed upon the map.

Dynamic Model
   A model that computes new attribute values as a function of attribute changes over time. It represents a change in the property of cells over time.

directional data type
   Data type for continuous data with a direction.

dynamic modelling script
   Script of a Dynamic Model.

dynamic section
   The iterative section in a dynamic modelling script which contains calc operations that are consecutively performed at each timestep.

expression
   A PCRaster map or an operation resulting in a PCRaster map.

initial section
   Cartographic Modelling script; the section containing the consecutively performed calc operations which describe the Cartographic Model. Dynamic Modelling script; the section which defines the map values of maps used in the dynamic section at the start of a model run.

input
   The set of one or more PCRaster maps, tables, time series or point data column files which are used in an operation to generate the result of the operation.

iterative section
   iterative section: a section in a dynamic modelling script that describes the temporal change in map values.

keyword
   A primitive word used in a script written in the PCRaster modelling language. It is a word which has a special meaning in the PCRaster modelling language. It is typed in lower case.

larger integer
   large integer: Optional cell representation for nominal and ordinal data type; cell values are stored in computer as INT4.

ldd data type
   Data type for maps that represent a local drain direction network.

local drain direction network
   PCRaster map representation of flow paths in a landscape. Each cell contains a pointer towards its downstream neighbour or no pointer in case it is a pit.

location attributes
   The information entities attached to a PCRaster map which give together all spatial properties of the map.

map
   The collection of digital information about a part of the earth's surface. The kind of maps used in PCRaster is the PCRaster map.

module
   A piece of the PCRaster package with a distinct functionality.

nominal data type
   Data type for classified data without order.

number of columns
   Location attribute. The number of columns in a PCRaster map.

number of rows
   Location attribute. The number of rows in a PCRaster map.

operation
   One static manipulation of one or more database components with one operator or several operators nested in one operation resulting in one or more database components (PCRaster map, table, time series, point data column file).

operator
   A primitive PCRaster 'function' of Map Algebra, Cartographic Modelling, Dynamic Modelling and GIS. It calculates a result on basis of one or more inputs. Both the result and the inputs may be a PCRaster map, table, time series or point data column file.

ordinal data type
   Data type for classified data with order.

outlet point
   The pit cell at the end of a downstream path from a cell in a local drain direction network.

PCRaster database
   The database of PCRaster; see database.

PCRaster map
   One of the kind of data in the PCRaster database. Contains spatial data of one attribute encoded in the form of a regular grid of cells covering an area. Binary format.

PCRaster modelling language
   The computer language provided by PCRaster for building Cartographic or Dynamic Models using a script.

pit
   A cell in a local drain direction network that only has neighbours pointing towards it and no neighbours at lower or equal elevation that it can point to.

point data column file
   One of the kind of data in the PCRaster database. Contains ascii formatted point data (x,y coordinates with attribute value(s)).

projection
   Location attribute. The projection of the real world co- ordinate system assigned to the PCRaster map. It is an x,y field (also used in basic mathematics). The x coordinates increase from left to right. The y coordinates increase from top to bottom or from bottom to top.

result
   The set of one or more PCRaster maps, tables, time series or point data column files that are generated by an operation.

row
   PCRaster map; a series of cells in horizontal direction.

scalar data type
   data type for continuous data that do not represent a direction.

script
   An ascii formatted computer programme of a Cartographic or Dynamic Model written in the PCRaster modelling language.

section
   A separate part of a script identified by a section keyword.

single real
   Default cell representation for scalar and directional data type; cell values are stored as REAL4 in computer.

section keyword
   A keyword that identifies the start of a section.

small integer
   Cell representation for boolean and ldd data type, default for nominal, ordinal and ldd data type; cell values are stored as UINT1 in the computer.

statement
   One line in a section of a script terminated with a semi colon (;).

table
   One of the kind of data in the PCRaster database. Contains relations between PCRaster maps. Ascii formatted.

timeseries
   One of the kind of data used in the PCRaster database. Contains a time series of (aggregated) cell values. Ascii formatted.

variable
   A PCRaster map, table, time series or point data column file in a Cartographic Model or Dynamic Model. Unlike a keyword its name is chosen by the model builder and starts with an upper case.


