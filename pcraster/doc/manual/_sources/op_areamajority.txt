

.. index::
   single: areamajority
.. _areamajority:

************
areamajority
************
.. topic:: areamajority

   Most often occurring cell value within an area

::

  Result = areamajority(expression, areaclass)

areaclass
   spatial
   boolean, nominal, ordinal

expression
   spatial
   boolean, nominal, ordinal

Result
   spatial
   type of expression

Operation
=========
areaclass Identifies the class to which a cell belongs: cells with corresponding values on areaclass are member of a separate class. For each separate class the most often occurring cell value on expression is determined. This value is assigned to all cells belonging to that class. This is done for all classes and saved as Result. If two values both occur the same number of times, the largest value of these values is assigned.  

Notes
=====


A cell on areaclass with missing value will result in a missing value on Result at the corresponding cell.  

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
   |    report Result = areamajority( Expr, Class);
   |   
   | • python
   |   Class = readmap("Class.map")
   |   Expr = readmap("Expr.map")
   |   Result = areamajority( Expr, Class)

   =============================================== ========================================== =============================================
   Result.map                                      Class.map                                  Expr.map                                     
   .. image::  ../examples/areamajority_Result.png .. image::  ../examples/areaarea_Class.png .. image::  ../examples/areamajority_Expr.png
   =============================================== ========================================== =============================================

   | 

