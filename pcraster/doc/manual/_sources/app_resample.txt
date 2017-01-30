

.. _resample:

********
resample
********
.. index::
   single: resample
.. topic:: resample

   Cuts one map or joins together several maps by resampling to the cells of the result map.

::

   resample [options] Map1 Map2....PCRmapN Result

Map1-N
  boolean, nominal, ordinal, scalar, directional (must have the same data type)
  spatial

Result
  data type of Map1, Map2,...MapN
  spatial


Options
=======

general options (additional options for specifying the location attributes of
Result are given in the description section): 


cell value assignment


-m
   for each cell on Result, the maximum value of the Map1, Map2,...Mapn cells that cover the Result cell is taken (see also description).

 


subpixellength


-e
   :emphasis:`subpixellength`: :emphasis:`subpixellength` must be a number equal to or between 0 and 2.5 (default 2.5). This options gives the subpixel length as percentage of the cell length on Result. Subpixels are used for calculating the percentage of a Result cell that is covered by a Map1, Map2,...Mapn cell (see also description). If -e 0 is set, the highest accuracy possible (smallest subpixel length) is taken.

 


missing value assignment


-p
   :emphasis:`percentmv`: :emphasis:`percentmv` must be a number equal to or between 0 and 100, default is 0. This option specifies the assignment of missing values. For each cell on Result, if the percentage of the cell that is covered by non missing value cells on Map1, Map2,...Mapn is less than :emphasis:`percentmv`, a missing value is assigned to the corresponding cell on Result.

-k 
   set the minimum and maximum of the result map to the minimum and maximum of the input maps.

Operation
=========

With the resample operator one or more input maps Map1, Map2,...Mapn are pasted in a new map Result. The location attributes will be changed to the location attributes of Result. The way this is done depends on some options: you can specify the location attributes of Result with a clone map and paste in that map, or determine the location attributes of Result on basis of the input expressions (cookie cutter or cell size modifier); these options will be described below. 


If several maps Map1, Map2,...Mapn are specified they must have the same data type and projection. The location attributes x\ :sub:`UL`,y\ :sub:`UL` coordinate, number of rows and columns and cell length may be different.  The angles of Map1, Map2,...Mapn may be different only if you use resample for pasting in a clone map, else these must be the same. If you specify more than one Map1, Map2,...Mapn, these maps may have any spatial location with respect to each other: they may overlap, may be adjoining or they may be separated in space. 


Almost in any case, the separate cells on Map1, Map2,...Mapn will not exactly overlap the separate cells on Result. So the raster data on Map1, Map2,...Mapn must be resampled to the raster of Result. For each cell on Result this is done as follows: for each cell on Map1, Map2,...Mapn which is partly or entirely in the cell on Result the area of the Result cell covered by that cell is calculated. This is done by laying down a fine raster of subpixels over the Result cell (default 40 x 40 subpixels per Result cell), and counting the number of subpixels covered by each Map1, Map2,...Mapn cell. These areas are used for assignment of the Result cell value: if the data type of Map1, Map2,...Mapn is scalar or directional an area weighted average of the Map1, Map2,...Mapn cell values is taken, where the weights are the numbers of subpixels covered by the cells. If the data type is boolean, nominal or ordinal the value is taken of the Map1, Map2,...Mapn cell which covers the largest area in the Result cell. If two or more cells both cover the same largest area, the maximum value of these cells is assigned. The maximum value is chosen in any case, if the option -m is set: the areas covered will be totally ignored. 


The subpixel length is specified by the option -e :emphasis:`subpixellength`, with a default length of 2.5 % of the Result cell length, which results in 40 x 40 subpixels per Result cell: smaller subpixels will reduce the error made in the computation of the areas, but the time needed to perform the operation will increase. 


As above said, the location attributes of the Result map can be specified in three ways: 


with a clone map



specifying clonemap

   :literal:`--clone`:emphasis:`Clonemap`: :emphasis:`Clonemap` is taken as clone. If a global clonemap is set as a :ref:`global option  <GOClone>`, the option can be omitted, the global clone map is taken as clone map. If the clone map is not set as a global option or if you want to use a different clone map than the global clone, you must specify the clone map in the command line with the option.

 


This functionality of resample is performed if no other options are used than the general options described at the start of the resample text. The clone map must be given as described above. Map1, Map2,...Mapn will be pasted in Result which has the location attributes of the clone map. The clone map and each map Map1, Map2,...Mapn must have the same projection. The other location attributes may be different. 


on basis of Map1, Map2,...Mapn (cookie cutter) 


specifying border around map(s)

   -b :emphasis:`borderwidth`: The smallest rectangle around cells (including missing value cells) is determined. Result will cover an area of this size plus borders or minus borders around this rectangle, where :emphasis:`borderwidth` is the width of the border. A positive :emphasis:`borderwidth` results in a larger map than the rectangle, a negative value in a smaller map. If -b 0 is specified Result will have (approximately) the size of the rectangle.

 




   -c :emphasis:`borderwidth`: idem, the smallest rectangle around non missing value cells is determined.






map expansion or contraction

   -x: if the area covered as defined by -b :emphasis:`borderwidth` or -c :emphasis:`borderwidth` contains a fractional number of rows and columns of cells on Result the number of rows and columns is rounded off upwards: the map is expanded (default).






   -a: if the area covered as defined by -b :emphasis:`borderwidth` or -c :emphasis:`borderwidth` contains a fractional number of rows and columns of cells on Result the number of rows and columns is rounded off downwards: the map is contracted.




This functionality of resample (cookie cutter) generates a Result with location attributes determined on basis of Map1, Map2,...Mapn. One of the options -b :emphasis:`borderwidth` or -c :emphasis:`borderwidth` must be specified and additionally -x or -a and the general options (described at the start of the resample text) may be given (optional). 


If more than one input map Map1, Map2,...Mapn is given, these must have the same projection and angle; the remaining location attributes may be different. Result will have the same projection and angle as the input maps; the cell size is taken from the first input map (Map1). The x\ :sub:`UL`,y\ :sub:`UL` coordinates and the number of rows and columns are calculated as follows: first the operations related to the options -b :emphasis:`bordersize` or -c :emphasis:`bordersize` are performed: the smallest rectangle around the edges of the input maps is determined, including or excluding missing values. The rectangle is enlarged or reduced by adding or removing a border at all sides of the map. This new rectangle is the approximate size of the Result, its top left vertex is the x\ :sub:`UL`,y\ :sub:`UL` coordinate of Result. Rows and columns of cells are laid down in the rectangle, starting at x\ :sub:`UL`, y\ :sub:`UL`. If the number of columns or rows needed to fill up the rectangle is a fractional number the rectangle is somewhat (always less than one cel length) expanded or contracted at the right and bottom sides until a whole number of rows and columns of cells fits into the rectangle. This number of rows and columns is assigned to Result.  Expansion or contraction is specified with -x (default) or -a, respectively. 


to modify cell length



celllength

   -r :emphasis:`celllength`: :emphasis:`celllength` is the cell length which is assigned to Result

 



:literal:`--unittrue` or :literal:`--unitcell` 
   :literal:`--unittrue`: :emphasis:`cellength` in the option -r is real distance (default)


   :literal:`--unitcell`: :emphasis:`cellength` in the option -r is distance in unit cell lengths




map expansion or contraction



-x: if the area covered by the smallest rectangle around the input maps
contains a fractional number of rows and columns of Result cells the number of rows and columns is rounded off upwards: the map is expanded (default). 
   


   -a: if the area covered by the smallest rectangle around the input mapsncontains a fractional number of rows and columns of Result cells the number of rows and columns is rounded off upwards: the map is contracted.

 


This functionality of resample is meant for changing the cell size of the first input map. No clone map must be given. The option -r :emphasis:`cellength` must be set, additionally you can specify :literal:`--unittrue` or  :literal:`--unitcell`, -x or -a or the general options described at the top of the resample text. 


It is quite unlikely that you want to specify more than one map, so
first the operation with one map is explained. Result will have the projection, angle, x\ :sub:`UL`, y\ :sub:`UL` coordinate of the input map Map1. The cell length of the input map is changed according to the option -r :emphasis:`cellength` and this length is assigned to Result. The area covered by the input map is filled up with cells of the new cell size, starting at x\ :sub:`UL`,y\ :sub:`UL`. If this results in a fractional number of rows and columns the map is somewhat (less than one new cell length) expanded (default) or contracted until a whole number of columns and rows is reached. This number of rows and collumns is assigned to Result. 


If more than one input map is given the operation performed corresponds
with the operation as a cookie cutter (described above), but you can
:emphasis:`not` use the options -b and -c: no borders can be specified. Result will approximately have the size of the smallest rectangle around cells (including missing value cells) on the input maps, x\ :sub:`UL`,y\ :sub:`UL` will be the top left vertex of the rectangle. 

See Also
========
:ref:`Import map types <secimportmaptype>`

