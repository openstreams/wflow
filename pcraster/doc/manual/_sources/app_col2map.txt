
.. index::
   single: col2map

.. _col2map:

*******
col2map
*******
.. topic:: col2map

   Converts from column file format to PCRaster map format

::

  col2map [options] columnfile PCRresult

columnfile
   asciifile

PCRresult
   spatial
   specified by data type option; if data type option is not set: data type of :emphasis:`PCRclone`

Options
=======


:literal:`--clone`:emphasis:`PCRclone`
   :emphasis:`PCRclone` is taken as clonemap. If you have set a global clonemap as :ref:`global option <GOClone>`, you don't need to set clone in the command line: the clonemap you have set as global option is taken as clonemap. If you have not set a global clonemap or if you want to use a different clonemap than the global clonemap, you must specify the clonemap in the command line with the clone option.



:literal:`--unittrue` or :literal:`--unitcell`

:literal:`--unittrue`
   coordinates in columnfile are interpreted as real distance (default)

:literal:`--unitcell`
   coordinates in columnfile are interpreted as distance in number of cell lengths

 


-B, -N, -O, -S, -D and -L 



   This data type option specifies the type options which is assigned to PCRresult (respectively boolean, nominal, ordinal, scalar, directional, ldd). If the option is not set, PCRresult is assigned the data type of :emphasis:`PCRclone` or the global clone. The data in columnfile must be in the domain of the data type which is assigned to PCRresult. For description of these domains see the description of the different :ref:`data types <secdatbasemaptype>`.

 

if option -D is set; --degrees of --radians

:literal:`--degrees`
   values on columnfile are interpreted as degrees (default)

:literal:`--radians`
   values on columnfile are interpreted as radians

 


-m :emphasis:`nodatavalue` 


   :emphasis:`nodatavalue` is the value in :emphasis:`columnfile` which is converted to a missing value on PCRresult. It can be one ascii character (letters, figures, symbols) or a string of ascii characters. For instance: -m -99.89 or -m j5w. Default, if this option is not set, 1e31 is recognized as a missing value.

 


-s :emphasis:`separator`


   By default, whitespace (one or more tabs, spaces) is recognized as separator between the values of a row in the columnfile. If the values are separated by a different separator, you can specify it with the option. The :emphasis:`separator` can be one of the ascii characters (always one). In that case, col2map recognizes the specified separator with or without whitespace as separator. For instance, if the values in columnfile are separated by a ; character followed by 5 spaces, specify -s ; in the command line (you do not need to specify the whitespace characters).

 


columnnumbers


-x :emphasis:`columnnumberx` 
   is the column number of the x coordinate in columnfile (default 1)


-y :emphasis:`columnnumbery` 
   :emphasis:`columnnumbery` is the column number of the y coordinate in columnfile (default 2)


-v :emphasis:`columnnumberv` 
   :emphasis:`columnnumberv` is the column number of the cell values in columnfile (default 3)

 


Each cell on PCRresult is assigned the cell value on columnfile which has x,y coordinates that define a point in that cell";" for assignment of values in columnfile which have x,y coordinates at the edges of cells on PCRresult, the following options are used: 

:literal:`--coorcentre`, :literal:`--coorul` or :literal:`--coorlr`

:literal:`--coorcentre` (default) or :literal:`--coorul` 
   values in columnfile that have x,y coordinates at the upper and left margins of a cell come into that cell, values at the bottom and right margins come into neighbouring cells. So, cell values with x, y coordinates at vertexes of cells come into the cell at the lower right side of the vertex.

:literal:`--coorlr`
   values in columnfile that have x, y coordinates at the bottom and right margins of a cell come into that cell, values at the upper and left margins come into neighbouring cells. So, cell values with x, y coordinates at vertexes of cells come into the cell at the upper left side of the vertex.

 


Options to specify which value is assigned if two or more values in
columnfile are found which all come into the same cell on PCRresult: 


-a, -h, -l, -H, -M, -t

-a
   average value of the values found within the cell is assigned (default for scalar and directional data; for directional data and assignment of records without a direction, see notes)

-h
   highest score: most occuring value found for the cell is assigned; if two values are found the same (largest) number of times, the highest value of these values is assigned, this is called a majority conflict (default for boolean, nominal, ordinal and ldd data)

-l
   lowest score: least occurring value found for the cell is assigned (option for nominal, ordinal, boolean, ldd data); if two values are found the same (smallest) number of times, the smallest value of these values is assigned, this is called a minority conflict.

-H
   highest value found for the cell is assigned (option for scalar or ordinal data)

-M
   lowest value found for the cell is assigned (option for scalar or ordinal data)

-t
   total (sum) of the columnfile values is assigned (option for scalar data)

 

Operation
=========

The columnfile is converted to PCRresult, which is an expression in PCRaster map format with the location attributes of :emphasis:`PCRclone`. The columnfile must be in the format described in :ref:`secdatbasepointform`.

For each cell on PCRresult the operator searches in columnfile for records that have x,y co-ordinates that come into that cell on PCRresult. If one single record is found, the value of this record is assigned to the cell, if several records are found, the value which is assigned is specified by the option (-a, -h, -l, -H or -M). A cell on PCRresult without a value on columnfile that falls into the cell is assigned a missing value on PCRresult. 

Notes
=====


Directional data: If the option -a (average, default) is set, and both records
without a direction (value -1) and records with a direction come into a cell
(a so called direction conflict), the records without a direction are
discarded and the cell value is computed from the records containing a
direction only. Thus a cell is assigned a no direction value (value -1) only
if all records for that cell don't have a direction.
Using col2map for generating a PCRresult of data type ldd is quite risky: probably it will result in a ldd which is unsound. If you do want to create a PCRresult of data type ldd use the operator lddrepair afterwards. This operator will modify the ldd in such a way that it will be sound, see the operator lddrepair. 

Group
=====
This operation belongs to the group of  Creation of PCRaster maps

See Also
========
:ref:`asc2map`

Examples
========
#. ::

      col2map --clone mapclone.map -S -m mv -v 4 ColFile1.txt Result1.map

   ==================================================== =========================================== ============================================
   `ColFile1.txt`                                       `Result1.map`                               `mapclone.map`                      
   .. literalinclude:: ../examples/col2map_ColFile1.txt .. image::  ../examples/col2map_Result1.png .. image::  ../examples/mapattr_mapclone.png
   ==================================================== =========================================== ============================================

#. ::

      col2map --clone mapclone.map -O -m mv -x 2 -y 3 -v 6 --coorlr -H ColFile2.txt Result2.map

   ==================================================== =========================================== ============================================
   `ColFile2.txt`                                       `Result2.map`                               `mapclone.map`                      
   .. literalinclude:: ../examples/col2map_ColFile2.txt .. image::  ../examples/col2map_Result2.png .. image::  ../examples/mapattr_mapclone.png
   ==================================================== =========================================== ============================================

