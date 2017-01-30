

.. index::
   single: areatotal
.. _areatotal:

*********
areatotal
*********
.. topic:: areatotal

   Sum of cell values within an area

::

  Result = areatotal(expression, areaclass)

areaclass
   spatial
   boolean, nominal, ordinal

expression
   spatial
   scalar

Result
   spatial
   scalar

Operation
=========
areaclass Identifies the class to which a cell belongs: cells with corresponding values on areaclass are member of a separate class. For each separate class the expression values of the cells belonging to that class are summed. This sum is assigned to all cells belonging to that class. This is done for all classes and saved as Result.  

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
   |    report Result = areatotal( Expr, Class);
   |   
   | • python
   |   Class = readmap("Class.map")
   |   Expr = readmap("Expr.map")
   |   Result = areatotal( Expr, Class)

   ============================================ ========================================== ============================================
   Result.map                                   Class.map                                  Expr.map                                    
   .. image::  ../examples/areatotal_Result.png .. image::  ../examples/areaarea_Class.png .. image::  ../examples/areamaximum_Expr.png
   ============================================ ========================================== ============================================

   | 

