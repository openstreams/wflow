

.. _mapattr:

*******
mapattr
*******
.. index::
   single: mapattr
.. topic:: mapattr

   Create a new PCRaster map, change or display location attributes of existing PCRaster map

::

  mapattr [options] PCRmap1 PCRmap2....PCRmapN

PCRmap
  boolean, nominal, ordinal, scalar, directional, ldd
  spatial


Options
=======

Options that apply to all modes of operations
---------------------------------------------


Set default attribute values with existing map


--clone
   :literal:`--clone PCRclone`. :emphasis:`PCRclone` is taken as clone map. This option can be set in the command line or as a :ref:`global option <GOClone>`. If a map attribute is not specified by the user (in the menu or with options) the attribute value of :emphasis:`PCRclone` is assigned to the map that is generated. If  :emphasis:`clone` is not set the following defaults are used: data type: boolean; cell representation: small integer; projection: y increases from bottom to top, upper left corner (:emphasis:`xUL`, :emphasis:`yUL`): (0,0), cell length 1 and angle 0.

 


Print information for all maps



-p
   gives information for all maps PCRmap1, PCRmap2,...PCRmapn (from top to bottom): n number of rows, number of columns, cell length, data type, cell representation, projection, angle in degrees, upper left corner (:emphasis:`xUL`, :emphasis:`yUL`), minimum value, maximum value, PCRaster version number, file_id (i.e. the number of the file, used internally by PCRaster), native/not native (specifies whether the byte order of the file is native ('y' is printed) or not native ('n' is printed): big endian or little endian), attribute table (specifies whether the file contains additional attributes ('y' is printed) like for example a legend or a colour palette or not ('n').

 

Options that apply to operation with the menu
---------------------------------------------


Change location attributes of an existing map



-e

   starts the menu allowing to change the location attributes projection, upper left corner (:emphasis:`xUL`, :emphasis:`yUL`), cell length and angle of the existing map PCRmap1. The other map attributes are fixed. Only one map (PCRmap1) can be specified. Do not use -e in combination with other options.

 


Copy location attributes of a map to other maps



-c

   when using this option two or more PCRmap1, PCRmap2,...PCRmapn must be specified. These maps must have corresponding number of rows and corresponding number of columns. The location attributes cell length, projection, angle, upper left corner (:emphasis:`xUL`, :emphasis:`yUL`), of the first map PCRmap1 are copied to the other map(s) PCRmap2,...PCRmapn. The other map attributes are not changed. Do not use -c in combination with other options.

 

Options that result in operation without the menu
-------------------------------------------------

It is possible to use mapattr without invoking the menu by setting the `-s` option.
If you do not invoke the menu you must define all attributes with the following options else the default values (or the values of :emphasis:`PCRclone`, see the  :ref:`global clone option <GOClone>`\ ) are assigned.


:literal:`-s`
   set attributes or create map without menu


 


-R :emphasis:`NumberOfRows` 


   Specifies the number of rows in the map that is generated. :emphasis:`NumberOfRows` must be a whole positive value.

 


-C :emphasis:`NumberOfColumns` 


   Specifies the number of columns in the map that is generated. :emphasis:`NumberOfColumns` must be a whole positive value.

 


-B,-N,-O,-S,-D or -L



   Specifies the data type of the map that is generated (respectively boolean, nominal, ordinal, scalar, directional, ldd). Default: boolean,nor if :literal:`--clone` is set, the data type of the clone map PCRclone.

 

 --small or --large


   In most case, the default cell representation will be sufficient. If you want, you can specify the cell representations:


   Nominal and ordinal data types

:literal:`--small`
   cell values are represented by small integer cell representation (default)

:literal:`--large`
   cell values are represented by large integer cell representation

 


-P yb2t (or -P yt2b)



   Specifies the projection: n-P yb2t, y increases from bottom to top (or -P yt2b, y increases from top to bottom). Default: yb2t or if :literal:`--clone` is set the projection of :emphasis:`PCRclone`. This option is an historical error, it should always be the default.

 


-x :emphasis:`XCorULC` 


   Specifies the x-coordinate of the upper left corner of the map that is generated. :emphasis:`XCorULC` is a real value. Default: 0.0 or if :literal:`--clone` is set, the value of the clone map :emphasis:`PCRclone`.

 


-y :emphasis:`YCorULC` 


   Specifies the y-coordinate of the upper left corner of the map that is generated. :emphasis:`YCorULC` is a real value. Default: 0.0 or if :literal:`--clone` is set, the value of the clone map :emphasis:`PCRclone`.

 


-l :emphasis:`CellLength` 


   Specifies the celllength of the map that is generated. :emphasis:`CellLength` is a real value. Default: 1 or if :literal:`--clone` is set, the value of the clone map :emphasis:`PCRclone`.

 


-a :emphasis:`Angle` 


   Specifies the angle of the map that is generated. :emphasis:`Angle` is the angle in degrees between -90 and 90 degrees. Default: 0.0 or if :literal:`--clone` is set, the value of the clone map :emphasis:`PCRclone`.

 

Operation
=========

The mapattr application generates a new PCRaster map with map attributes specified by the user, changes the :ref:`location <secdatbasemaphead>` attributes of an existing map or prints map attribute information.

Operation with the menu
-----------------------
A new map is generated by specifying one input file PCRmap1 and not setting options. The application invokes a menu. In the menu the location attributes can be entered of the new map PCRmap1 that is created. You can scroll through the menu with the arrow up or arrow down keys or the keys listed at the bottom of the menu. An menu item is entered by pressing <Enter> and typing the value followed by <Enter>. The menu items 'data type' and 'projection'  are filled in by selecting one of the options with the arrow left or arrow right keys (or one of the keys listed at the bottom of the menu) instead of typing the entry. Quit the menu by pressing 'q'. The program asks you whether the map must be created. You can answer by pressing 'Y' (Yes, create the map), 'N' (No, do not create the map and leave the menu) or by pressing <Escape which effects that you can resume editing.


Location attributes of one map can be changed by specifying one existing
map PCRmap1 and the option -e. Note that this is not meant for cutting or resampling the map which is done with the application resample. See also the section on :ref:`the import map type <secimportmaptype>`.


Location attributes of the first, existing, PCRmap1 are copied to the existing PCRmap2,...PCRmapn with the option -c. 

 
Map attribute information of multiple PCRmap1-n is printed with the option -p.  

Operation without the menu
--------------------------

The map attributes can also be entered without the menu, by setting the -s options.
All modes of operation possible with the menu can also be executed in this command line/no menu mode. So you can also use the options -e or -c for changing map attributes of an existing map and copying map attributes from one map to other maps respectively.


If one of the map attributes is not specified with an option in the
command line, mapattr assigns the default value or, if a clone map is specified with the option :literal:`--clone`, the value of the :emphasis:`PCRclone`. See also the :literal:`--clone` option in the option list at the top of the mapattr description. 

See Also
========
:ref:`resample`

Examples
========

Generation of a new map that is stored under the filename
'mask.map'; invokes the menu:

.. parsed-literal::

    mapattr mask.map 

Changing the location attributes of an existing map :emphasis:`clone.map`; invokes the menu:  

.. parsed-literal::

  mapattr -e clone.map

Copying the location attributes of clone1.map to dem.map and ldd.map; does not invoke the menu: 

.. parsed-literal::

  mapattr -c clone1.map dem.map ldd.map

Printing (on the screen) the map info for clone1.map, dem.map and ldd.map: 

.. parsed-literal::

  mapattr -p clone1.map dem.map ldd.map


Generation of a new map that is stored under the filename mask2.map; does not invoke the menu (here, the option -a is not set so the angle of the map is assigned the default value 0): 

.. parsed-literal::

  mapattr -s -R 19 -C 68 -B -P yb2t -x 12 -y -14.02 -l 0.8 mask2.map
