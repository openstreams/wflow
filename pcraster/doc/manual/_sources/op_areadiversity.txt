

.. index::
   single: areadiversity
.. _areadiversity:

*************
areadiversity
*************
.. topic:: areadiversity

   Number of unique cell values within an area

::

  Result = areadiversity(expression, areaclass)

expression
   spatial
   boolean, nominal, ordinal

areaclass
   spatial
   boolean, nominal, ordinal

Result
   spatial
   scalar

Operation
=========
areaclass Identifies the class to which a cell belongs: cells with corresponding values on areaclass are member of a separate class. For each separate class the number of unique cell values on expression is counted. This number is assigned to all cells belonging to that class. This is done for all classes and saved as Result.  

Notes
=====


A cell with missing value on areaclass will result in a missing value on Result at the corresponding cell.  

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
   |   binding
   |    Result = Result.map;
   |    Class = Class.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result = areadiversity( Expr, Class);
   |   
   | • python
   |   Class = readmap("Class.map")
   |   Expr = readmap("Expr.map")
   |   Result = areadiversity( Expr, Class)

   ================================================ ========================================== ==============================================
   Result.map                                       Class.map                                  Expr.map                                      
   .. image::  ../examples/areadiversity_Result.png .. image::  ../examples/areaarea_Class.png .. image::  ../examples/areadiversity_Expr.png
   ================================================ ========================================== ==============================================

   | 

