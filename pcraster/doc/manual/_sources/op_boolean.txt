

.. index::
   single: boolean
.. _boolean:

*******
boolean
*******
.. topic:: boolean

   Conversion data type to boolean data type

::

  Result = boolean(expression)

expression
   spatial, non spatial
   nominal, ordinal, scalar, directional, ldd

Result
   dimension of expression
   boolean

Operation
=========


If expression is a PCRaster map or a calculation resulting in a PCRaster map of one  of the data types nominal, ordinal, scalar, directional or ldd, converts to a Boolean data type with Boolean cell values, where 1 is TRUE and 0 is FALSE. Each expression cell value not equal to 0 is assigned a 1 (TRUE) on Result; each expression cell value equal to 0 is assigned a 0 (FALSE) on Result. Or it generates a map of boolean data type with one constant value. 



If expression has no PCRaster data type, boolean generates a boolean Result. This is the case if expression is a number. This number must be in the domain of the boolean map type, i.e. 0 or 1. Result will be a map with the same location attributes as the :ref:`global clone map <GOClone>` ; all cells will have the value of expression.  

Notes
=====


A cell with missing value on expression is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Conversion and assignment 

See Also
========
:ref:`secdatbasemaptype`, :ref:`formboolean`:ref:`DataTyToNumb`

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result = boolean(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = boolean(Expr)

   ========================================== ========================================
   Result.map                                 Expr.map                                
   .. image::  ../examples/boolean_Result.png .. image::  ../examples/boolean_Expr.png
   ========================================== ========================================

   | 

