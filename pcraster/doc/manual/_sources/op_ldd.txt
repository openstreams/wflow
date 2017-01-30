

.. index::
   single: ldd
.. _ldd:

***
ldd
***
.. topic:: ldd

   Data conversion from specific data types to local drain direction data type

::

  Result = ldd(expression)

expression
   spatial, non spatial
   nominal, ordinal, directional

Result
   spatial, non spatial
   ldd

Options
=======

conversion from directional data type: :literal:`--degrees` or :literal:`--radians`

:literal:`--degrees`
   values on expression are interpreted as degrees (default)

:literal:`--radians`
   values on expression are interpreted as radians



Operation
=========


If the expression is a PCRaster map or a calculation resulting in a PCRaster map,  it is converted. If the  expression has a data type nominal or ordinal, only the data type is changed; the cell values on Result correspond with the values on expression. For this conversion it is required that the cell values (or directions) on expression are in the domain of the ldd data type, i.e. a whole number from 1 up to and including 9. The values resemble the layout of the numeric key pad of your computer.  If expression has a directional data type, the circular directional scale of expression is converted to the local drain direction codes of the ldd data type as follows: the local drain direction codes are interpreted as real clockwise directions where a local drain direction to the top of the map (ldd code 8) is 0 degrees. Each directional cell value on expression is assigned the ldd code of the local drain direction which is closest to the direction given by expression. For instance, assuming the option :literal:`--degrees` is set, all expression values equal to 22.5 or between 22.5 and 67.5 (i.e values in the range [22.5,67.5>) are assigned a ldd code 9 on Result.
  Or it generates a map of local drain direction data type with one constant value.  



If expression has no PCRaster data type, a Result with data type ldd is generated. This is the case if expression is a number. The value of expression must be in the domain of the ldd data type, i.e. a whole number from 1 up to and including 9. Result will be a map with the same location attributes as the :ref:`global clone map <GOClone>`; all cells will have the value of expression.  



Notes
=====


A missing value on expression is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Conversion and assignment 

See Also
========
:ref:`secdatbasemaptype`, :ref:`formldd`, :ref:`DataTyToNumb`

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result = ldd(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = ldd(Expr)

   ====================================== ====================================
   Result.map                             Expr.map                            
   .. image::  ../examples/ldd_Result.png .. image::  ../examples/ldd_Expr.png
   ====================================== ====================================

   | 

