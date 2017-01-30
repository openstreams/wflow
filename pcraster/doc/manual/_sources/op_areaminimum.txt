

.. index::
   single: areaminimum
.. _areaminimum:

***********
areaminimum
***********
.. topic:: areaminimum

   Minimum cell value within an area

::

  Result = areaminimum(expression, areaclass)

expression
   spatial
   ordinal, scalar

areaclass
   spatial
   boolean, nominal, ordinal

Result
   spatial
   type of expression

Operation
=========
areaclass Identifies the class to which a cell belongs: cells with corresponding values on areaclass are member of a separate class. For each separate class the minimum expression value of the cells belonging to that class is determined. This value is assigned to all cells belonging to that class. This is done for all classes and saved as Result.  

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
   |   binding
   |    Result = Result.map;
   |    Class = Class.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result = areaminimum( Expr, Class);
   |   
   | • python
   |   Class = readmap("Class.map")
   |   Expr = readmap("Expr.map")
   |   Result = areaminimum( Expr, Class)

   ============================================== ========================================== ============================================
   Result.map                                     Class.map                                  Expr.map                                    
   .. image::  ../examples/areaminimum_Result.png .. image::  ../examples/areaarea_Class.png .. image::  ../examples/areamaximum_Expr.png
   ============================================== ========================================== ============================================

   | 

