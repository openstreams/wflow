

.. index::
   single: nominal
.. _nominal:

*******
nominal
*******
.. topic:: nominal

   Data conversion data type nominal data type

::

  Result = nominal(expression)

expression
   spatial, non spatial
   boolean, ordinal, scalar, directional, ldd

Result
   dimension of expression
   nominal

Operation
=========


If expression is a PCRaster map or a calculation resulting in a PCRaster map, it is converted: if expression is of one of the data types boolean, ordinal or ldd, the cell values on expression are assigned without change to the corresponding cells on Result; if expression is of data type scalar or direction, the values on expression are truncated.
 Or the operation generates a map of nominal data type with one constant value.  



If expression is not a PCRaster map, a nominal Result is generated. This is the case if expression is a number. This number must be in the domain of the nominal map type, i.e. a whole value. Result will be a map with the same location attributes as the  :ref:`global clone map <GOClone>`; all cells will have the value of expression.   

Notes
=====


A cell with missing value on expression is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Conversion and assignment 

See Also
========
:ref:`secdatbasemaptype`, :ref:`formnominal`, :ref:`DataTyToNumb`

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result = nominal(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = nominal(Expr)

   ========================================== ==========================================
   Result.map                                 Expr.map                                  
   .. image::  ../examples/nominal_Result.png .. image::  ../examples/rounddown_Expr.png
   ========================================== ==========================================

   | 

