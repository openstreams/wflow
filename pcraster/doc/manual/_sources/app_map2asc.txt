

.. _map2asc:

*******
map2asc
*******
.. index::
   single: map2asc
.. topic:: map2asc

   Converts from PCRaster maps format to ascii file format

::

  map2asc [options] PCRmap asciifile

PCRmap
   spatial
   boolean, nominal, ordinal, scalar, directional, ldd

asciifile
   formatted asciifile
   None

Options
=======

Options can be given for assigning a special layout to asciifile. These options are described in the operation section. Other options are: 


-m :emphasis:`nodatavalue` 


   :emphasis:`nodatavalue` may be one ascii character (letters, figures, symbols) or a string of ascii charaters. For instance: -m -999, -m ? or -m 8kf.f. The missing values on the map will be assigned the :emphasis:`nondatavalue` on the asciifile

 


-s :emphasis:`seperator` 


   :emphasis:`separator` may be one ascii character (letters, figures, symbols, space, tab) or a string of ascii characters. A space or tab is specified using "\" \"", for instance: -s "\" \"" is a space as separator. The cell values which are on one row on PCRmap will be separated by :emphasis:`separator` on asciifile. Default, if -s is not set, map2asc prints one space :emphasis:`before` each cell value field. Note that if you want to convert the asciifile back to PCRaster map format with the asc2map operator, it should contain whitespace characters only or whitespace characters with only :emphasis:`one` non whitespace character as separator.

 


-f C :emphasis:`-typeformat` 


   This option is used to specify the sort of format which is assigned to thencell values in asciifile. See the discussion of -f in map2col.

 

Operation
=========

simple conversion
-----------------

The data in PCRaster map format on PCRmap are converted to ascii format and saved as asciifile. Default a simple conversion is performed. This is a rowwise conversion to an asciifile without header. The rows on PCRmap are scanned from left to right, starting with the top row and ending with the bottom row. Each time a value is scanned it is added to asciifile, starting with a new line on asciifile if a new row on PCRmap is scanned. The values in one row on PCRmap will be on one line on asciifile, with a separator defined by the option - s (default one or more spaces, see above). The lines on asciifile always end in a cell value.  

conversion to ARC/INFO input format
-----------------------------------

In ARC/INFO ascii data with a special lay-out containing a header with location attributes and the :emphasis:`nodatavalue` can be imported using the ARC/INFO command asciigrid. An asciifile with the lay-out needed for this command can be created with map2asc specifying the option -a. Additionally only the option -m can be specified. 

conversion to other ascii file lay-outs
---------------------------------------

Two options can be set to specify asciifiles with other lay-outs: 

-c
   Default a rowwise output is performed: lines on asciifile correspond with rows on PCRmap. If this option is set the output will be columnwise: the first column on PCRmap (from top to bottom) is printed as the first line (from left to right) on asciifile, the second column as the second line etc. The number of columns on PCRmap will correspond with the number of rows on asciifile.

  

-n
   :emphasis:`numberofcellsonline` must be a whole value larger or equal to 1. Default, each row of cell values on PCRmap is saved as one line in asciifile. This option can be used to print a different number of cell values (the value of :emphasis:`numberofcellsonline`) on one line in asciifile.

  

Notes
=====


Group
=====
This operation belongs to the group of  Data export from PCRaster map

See Also
========
:ref:`map2col`

