

.. index::
   single: areaaverage
.. _areaaverage:

***********
areaaverage
***********
.. topic:: areaaverage

   Average cell value of within an area

::

  Result = areaaverage(expression, areaclass)

expression
   spatial
   scalar

areaclass
   spatial
   boolean, nominal, ordinal

Result
   spatial
   scalar

Operation
=========
areaclass Identifies the class to which a cell belongs: cells with corresponding values on areaclass together form a separate class. For each separate class the expression values of the cells belonging to that class are averaged. This average value is assigned to all cells belonging to that class. This is done for all classes and saved as Result.  

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
   |    report Result = areaaverage( Expr, Class);
   |   
   | • python
   |   Class = readmap("Class.map")
   |   Expr = readmap("Expr.map")
   |   Result = areaaverage( Expr, Class)

   ============================================== ========================================== ============================================
   Result.map                                     Class.map                                  Expr.map                                    
   .. image::  ../examples/areaaverage_Result.png .. image::  ../examples/areaarea_Class.png .. image::  ../examples/areamaximum_Expr.png
   ============================================== ========================================== ============================================

   | 

