

.. index::
   single: areaarea
.. _areaarea:

********
areaarea
********
.. topic:: areaarea

   The area of the area to which a cell belongs

::

  Result = areaarea(areaclass)

areaclass
   spatial
   boolean, nominal, ordinal

Result
   spatial
   scalar

Options
=======
:literal:`--unittrue` or :literal:`--unitcell`

:literal:`--unittrue`
   area is computed in true area represented by cells (default)

:literal:`--unitcell`
   area is computed in number of cells



Operation
=========


The class to which a cell belongs is identified by areaclass: all cells with corresponding values on areaclass are grouped and together they form one class. For each separate class the total area that is represented by the cells belonging to that class is calculated (cell value times the total number of classes).  This value is assigned to all cells belonging to that class. This is done for all classes and saved as Result.  

Notes
=====


A cell with missing value on areaclass is assigned a missing value on Result at the corresponding cell.  

Group
=====
This operation belongs to the group of  Area operators 

See Also
========
:ref:`secstatar`

Examples
========
#. 
   | • pcrcalc
   |   #! --unittrue
   |   binding
   |    Result = Result.map;
   |    Class = Class.map;
   |   initial
   |    report Result = areaarea( Class);
   |   
   | • python
   |   setglobaloption("unittrue")
   |   Class = readmap("Class.map")
   |   Result = areaarea( Class)

   =========================================== ==========================================
   Result.map                                  Class.map                                 
   .. image::  ../examples/areaarea_Result.png .. image::  ../examples/areaarea_Class.png
   =========================================== ==========================================

   | 

