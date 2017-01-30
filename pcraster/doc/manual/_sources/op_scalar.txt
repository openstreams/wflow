

.. index::
   single: scalar
.. _scalar:

******
scalar
******
.. topic:: scalar

   Data conversion to the scalar data type

::

  Result = scalar(expression)

expression
   spatial, non spatial
   boolean, nominal, ordinal, directional, ldd

Result
   dimension of expression
   scalar

Operation
=========



If expression is a PCRaster map or a calculation resulting in a PCRaster map, it is converted: the cell values of expression are assigned without change to the corresponding cells on Result.  Or it generates a map of scalar data type with one constant value.  



If expression has no PCRaster data type, a Result with data type scalar is generated. This is the case if expression is a number or a calculation with numbers. Result will be a map with the same location attributes as the :ref:`global clone map <GOClone>`; all cells will have the value of expression.   

Notes
=====


A cell with missing value on expression is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Conversion and assignment 

See Also
========
:ref:`secdatbasemaptype`, :ref:`formscalar`, :ref:`DataTyToNumb`

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Expr1 = Expr1.map;
   |   areamap
   |    Expr1.map;
   |   initial
   |    report Result = scalar(1);
   |   
   | • python
   |   Expr1 = readmap("Expr1.map")
   |   setclone("Expr1.map");
   |   Result = scalar(1)

   ========================================= =====================================
   Result.map                                Expr1.map                            
   .. image::  ../examples/scalar_Result.png .. image::  ../examples/and_Expr1.png
   ========================================= =====================================

   | 

