

.. _map2col:

*******
map2col
*******
.. index::
   single: map2col
.. topic:: map2col

   Converts from PCRaster map format to column file format

::

  map2col [options] PCRmap1 PCRmap2....PCRmapN columnfile

PCRmap
 boolean, nominal, ordinal, scalar, directional, ldd
 spatial
 The maps must have the same projection, the other location attributes and the data types may be different between the maps.

columnfile
  asciifile

Options
=======
:literal:`--unittrue` or :literal:`--unitcell`

:literal:`--unittrue`
   coordinates in columnfile are interpreted as real distance (default)

:literal:`--unitcell`
   coordinates in columnfile are interpreted as distance in number of cell lengths

 


coordinate positions:

:literal:`--coorcentre`
   cell values on columnfile are assigned the coordinates of the centres of the cells on PCRmap1 (default).

:literal:`--coorul`
   cell values on columnfile are assigned the coordinates of the upper left corner of the cells on PCRmap1.

:literal:`--coorlr`
   cell values on columnfile are assigned the coordinates of the lower right corner of the cells on PCRmap1.

 


columnnumbers

-x :emphasis:`columnnumberx`
   :emphasis:`columnnumberx`  is the column number of the x coordinate in columnfile; in append mode: the colum number of the x coordinate in the :emphasis:`inputcolumnfile` (default 1).

-y :emphasis:`columnnumbery`
   :emphasis:`columnnumbery` is the column number of the y coordinate in columnfile; in append mode: the colum number of the y coordinate in the :emphasis:`inputcolumnfile` (default 2).

 


-m :emphasis:`nodatavalue` 


   :emphasis:`nodatavalue` is the value in :emphasis:`columnfile` which is converted to a missing value on PCRresult. It can be one ascii character (letters, figures, symbols) or a string of ascii charaters. For instance: -m -99.98 or -m ! or -m j5w. Default, if this option is not set, 1e31 is recognized as a missing value.

 


-M



   Default, if PCRmap1 has a missing value in a cell, the cell is not saved in columnfile. If the option -M is set these cells :emphasis:`are` saved in columnfile. They are assigned the :emphasis:`nodatavalue`.

 


-s :emphasis:`seperator` 


   By default, whitespace (one or more tabs, spaces) is recognized as separator between the values of a row in the columnfile. If the values are separated by a different separator, you can specify it with the option. The :emphasis:`separator` can be one of the ascii characters (always one). In that case, col2map recognizes the specified separator with or without whitespace as separator. For instance, if the values in columnfile are separated by a ; character followed by 5 spaces, specify -s ; in the command line (you do not need to specify the whitespace characters).

 


specifying format of columnfile 

-p
   columnfile is in plain format without header (default)

-g
   columnfile is in simplified :ref:`Geo-EAS format <secdatbasepointform>`.

 


row wise or column wise output 

-r
   row wise output (default). The cell values in the first row onnPCRmap1, PCRmap2,...PCRmapn will be at the top of the columnfile, underneath the cells in the second row, etc.

-c
   column wise output. The cell values in the first column onnPCRmap1, PCRmap2,...PCRmapn will be at the top of the columnfile, underneath the cells in the second column, etc.

 


-a :emphasis:`inputcolumnfile` 


   append mode. The inputcolumnfile is the name of the column file to which the data are appended. The file with the appended data is saved as columnfile. See operation in append mode below.

 


-f :emphasis:`C-typeformat` 
 This option is used to specify the sort of format which is assigned to the cell values in columnfile. The format determines the maximal and minimal cell value which can be converted, the precision the data can be saved in columnfile and the number of positions used for each cell value. The default format that is used depends on the data type of PCRmap. For boolean, nominal and ordinal maps, containing only whole values, the smallest possible number of positions is used for each cell value field in columnfile, taking into account all cell values and the number of positions needed for the missing value. For instance, a nominal map with nominal values between 12 and 19 and a no data value -999 (given by the option -m :emphasis:`nodatavalue`) is written to columnfile using four positions for each cell value (resulting from the number of positions needed for the :emphasis:`nodatavalue`). For instance, the value 12 is printed as a cell value field made up of 2 spaces followed by 12.

Maps of scalar and directional data type are always printed in the :emphasis:`C-type format` 11.6g, also used in the C programming language: the precision the data are saved on columnfile corresponds with the precision they are available on PCRmap, with the restriction that the maximal number of significant figures which can be saved on columnfile is six per cell value. The maximum and minimum value which can be saved is 10\ :sup:`99` and -10\ :sup:`99` respectively; a notation with base\ :sub:`10` exponents is used if the value is larger than 10\ :sup:`6` or smaller than -10\ :sup:`6`.  

Examples::

   PCRmap         columnfile
   -894.41000     -894.41
   -5674935         -5.67494e+06
   453628190.6       4.53628e+08
   0.000000000031    3.1e-11
   -0.02000012      -0.0200001
   -1.0200001       -1.02

If you want to prevent the usage of base 10 exponents for scalar or directional data use the C-type format f and specify -f :emphasis:`a`\ .\ :emphasis:`d`\ f, where :emphasis:`a` and :emphasis:`d`  must be whole numbers equal to or larger than 0. The value :emphasis:`d` is the number of decimal figures which will be used for each cell value, the value :emphasis:`a` is the minimal total number of positions used for each value; if more positions are needed (large values), more positions are used.

Examples::

    C-typeformat    PCRmap        columnfile
    5.6f            1234.1981     1234.198100
    5.3f            1234.1981     1234.198
    5.0f            1234.1981     1234
    3.1f            1234.1981     1234.2
    15.6f           1234.1981     1234.1981
    2.1f            1289128932.75 1289128932.7

You can also specify other C-type formats, see for description of these formats a C programming language standard work.

 

Operation
=========

Default operation (no append mode)
----------------------------------

Operation if PCRmap1, PCRmap2,...PCRmapn have corresponding location attributes:

In most cases only one PCRmap1 is given in the command line or several maps PCRmap1, PCRmap2,...PCRmapn are given which have the same location attributes. In these cases, the operation is performed  as follows. The PCRaster expression(s) PCRmap1, PCRmap2,...PCRmapn are converted to columnfile, which is an ascii file in column format. Each line in this columnfile represents one cell on PCRmap1, PCRmanp2,...PCRmapn. The x and y coordinates of the cells will be in the column numbers specified by the options -x :emphasis:`columnumberx` and -y :emphasis:`columnnumbery`. The cell values of PCRmap1 will be in the first 'empty' column, the values of PCRmap2 in the next column etc. For instance if you set -x 2 -y 3, values of PCRmap1 are written in the first column, values of PCRmap2 in the fourth, values of PCRmap3 in the fifth etc. 


Operation if PCRmap1, PCRmap2,...PCRmapn have different location attributes:

If more than one PCRmap is given in the command line and the given maps have different location attributes (with the exception of the projection which must be the same) the operation is performed in a somewhat different way. Only the x, y coordinates of the first map PCRmap1 are printed in columnfile. Real world coordinates or cell coordinates are printed, as specified by the option :literal:`--unittrue`, :literal:`--unitcell`; the coordinate position that is printed is specified by the option :literal:`--coorul`, :literal:`--coorlr`, :literal:`--coorcentre`. The cell values are printed as follows: first, each x, y coordinate pair is supplemented with its cell value of PCRmap1. Than, each line in columnfile is supplemented with the cell values of the remaining maps PCRmap2,...PCRmapn. No new lines are appended for these maps. For each of these maps and each line, the cell value is printed of the cell which has a real world location that corresponds with the real world location of the PCRmap1 cell that is already printed on that line. The real world location corresponds if the real world x,y coordinate of the PCRmap1 cell comes into the cell of the PCRmap2,...PCRmapn under consideration (the x,y coordinate of upper left corner, lower right corner or centre of each PCRmap1 cell is used, as specified by the :literal:`--coorcentre`, :literal:`--coorlr`, :literal:`--coorcentre` option).  A line on the columnfile that represents a PCRmap1 cell with a real world x,y coordinate that does not come into a cell on the PCRmap2,...PCRmapn under consideration is assigned a missing value in the appended field. 

operation in append mode
------------------------

Data can also be appended to an existing column file :emphasis:`inputcolumnfile`. This columnfile may be a plain column file without a header or a column file in simplified Geo-EAS format. These formats don't need to be specified, the map2col operator will detect the format of inputcolumnfile itself. For the append mode, the option -a :emphasis:`inputcolumnfile` is used. The :emphasis:`inputcolumnfile` is the name of the column file to which the data are appended. The file with the appended data is saved as columnfile. The data are appended as follows: each line (record) of the :emphasis:`inputcolumnfile` is supplemented with the PCRmap1, PCRmap2,...PCRmapn values of the cell in which the x,y coordinates of the line (record) are. On each line, the values of PCRmap1, PCRmap2,...PCRmapn will be typed in the order they are specified in the command line (i.e. PCRmap1 values are printed in the first column after the columns in the :emphasis:`inputfile`, PCRmap2 in the second column, etc.). 


The append mode results in the appending of columns only, no lines will
be appended: a cell on PCRmap1, PCRmap2,...PCRmapn without a x,y coordinate in the :emphasis:`inputcolumnfile` that comes into the cell will :emphasis:`not` be saved in a new line (record). A line (record) on :emphasis:`inputcolumnfile` with a x,y coordinate that does not come into a cell on PCRmap1, PCRmap2,...PCRmapn is assigned a missing value in the appended column(s). 


In append mode the options -M, -r, -c, -p and -g must not be used. The
other options can be used as normal, but the default values will be
appropriate in almost any case; the flags -m, -f and -s will only affect the
columns which are appended. The options :literal:`--coorcentre`, :literal:`--coorul` and :literal:`--coorlr` have the following meaning when used in append mode:   

:literal:`--coorcentre` (default) and :literal:`--coorul`
  lines (records) with coordinates in :emphasis:`inputcolumnfile` that are exactly at the upper or left edge of a cell are supplemented with the cell value of that cell, records with coordinates at the lower or right edge are supplemented with the value of a neighbouring cell.

   

:literal:`--coorlr`
  lines (records) with coordinates in :emphasis:`inputcolumnfile` that are exactly at the lower or right edge of a cell are supplemented with the cell value of that cell, records with coordinates at the upper or left edge are supplemented with the value of a neighbouring cell.

 

Notes
=====


Group
=====
This operation belongs to the group of  Creation of PCRaster maps

See Also
========
:ref:`map2asc`

