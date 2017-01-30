

.. _legend:

******
legend
******
.. index::
   single: legend
.. topic:: legend

   Attaches a legend to or changes the legend of one or more maps.

::

  legend [options] PCRmap1 PCRmap2.....PCRmapN


PCRmap1-N
  boolean, nominal, ordinal; they must have the same data type

Options
=======

By default, legend starts a menu. Options that apply to this default menu mode are -l and -h. 


enter also labels for cell values with a smaller value than which occurs on PCRmap11....PCRmapN

-l :emphasis:`minvalues`
   the menu will allow you to fill in labels for wholencell values equal to or between the (whole) value :emphasis:`minvalue` and the maximum value on the PCRmap1,PCRmap2,...PCRmapN.

 


enter also labels for cell values with a higher value than which occurs on PCRmap1....PCRmapN

-h :emphasis:`maxvalues`
   the menu will allow you to fill in labels for wholencell values equal to or between the minimum value on thenPCRmap1,PCRmap2,...PCRmapn and the (whole) number :emphasis:`maxvalue`.

 


The operator can also be used in a mode without the menu. The options
that invoke this mode are described in the second part of the operation
section. 


Operation
=========

The one or more maps PCRmap1,PCRmap2,...PCRmapn must have the same data type (boolean,nominal or ordinal), the location attributes do not need to correspond. The PCRmap1,PCRmap2,...PCRmapn may have (corresponding or different) legends attached or may not have legends attached. The legend operator will overwrite or assign in both cases the same legend to all PCRmap1,PCRmap2,...PCRmapn.  

Operation with the menu
-----------------------

This is the default action when no options or the options -l or -h are set. The operator invokes a menu which is used to change or enter the legend labels for the cell values that are found on one or more PCRmap1,PCRmap2,...PCRmapn. In addition the title of the legend can be assigned. The keys that are supported by the menu are given at the bottom of the menu.   

Operation without the menu
--------------------------
The operator can also be used without the menu. This is done by setting one of the following options. 


copy the legend of the first map to the other maps
   -c: More than one PCRmap1,PCRmap2,...PCRmapN are specified. The legend of the first map PCRmap1 is copied to the other maps PCRmap2,...PCRmapN.

 


store the legend labels in an ascii formatted legend file
   -w :emphasis:`outputlegendfile`: One or more PCRmap1,PCRmap2,...PCRmapN may be specifed. The cell values of these maps with the labels are stored in the ascii formatted :emphasis:`outputlegendfile`. For each cell value, the label is stored of the PCRmapI (I is 1...N) with the lowest I that contains the cell value under consideration. The layout of the legend file :emphasis:`outputlegendfile` is given below.

read the legend labels for the maps from an ascii formatted legend file
   -f :emphasis:`inputlegendfile`: One or more PCRmap1,PCRmap2,...PCRmapN may be specified. The labels given in the :emphasis:`inputlegendfile` are assigned to the legend of these maps. A cell value on an input map that does not occur in the :emphasis:`inputlegendfile` is assigned a label. The layout of the legend file :emphasis:`inputlegendfile` is given below.


The general layout of a legend file is as follows. The ascii formatted
file constists of two columns, separated by one or more whitespace
characters (space(s), tab(s)). The cell values are in the first column,
the labels are in the second column. The first row contains in the first
column the field -0 and in the second column the title of the legend.
The following rows contain in the first column the cell value and in the
second column the label for that cell value. 

When the legend file is the output from the legend operator, the title row (if a title occurs in one or more of the input maps) is the first row. The cell values with the labels are written in rising order from top to bottom, starting with row two.  

When the legend file is used as input file, the row with the title may
be omitted. If it is given it :emphasis:`must` be the first row of the file. The following rows contain the cell values with the labels. The order of these rows does not matter.  


An example of an input legend file is: 


::

   -0  landuse
    4  arable land
    1  woodland
    3  buildings
    2  lake
    2  river


