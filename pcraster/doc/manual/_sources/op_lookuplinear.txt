.. _lookuplinear:

************
lookuplinear
************
.. index::
   single: lookuplinear
.. index::
   single: interpolation
.. topic:: lookuplinear

   Assigns table key values with possible interpolation between key values.

::

   Result = lookuplinear(table,expression)

table
  ascii text table
expression
  scalar, spatial, non-spatial
Result
  scalar, spatial


Operation
=========

The value of expression is compared to the values of the search key (first column) of the table and is assigned the value of the corresponding record (second column) on the same row if the expression matches the value of the search key. If the value of expression lies between two values of the search key of the table a linear interpolation is executed between the two corresponding records. The interpolation is executed according to the formula:

Result = ((expression – lowerKey) / (upperKey – lowerKey)) * (lowerRecord – upperRecord) + lowerRecord

Notes
=====

If the cell value of the expression lies outside the range of the search key a missing value is assigned.
If a cell of the expression has a missing value a missing value is assigned to the result.
The search key values of the table should be in ascending order or missing values are assigned to all cells of Result.

Examples
========

#. 
   | • pcrcalc
   |   binding
   |    Result1 = Result1.map;
   |    Table = Table.txt;
   |    Expr = Expr.map;
   |   initial
   |    report Result1 = lookuplinear(Table,Expr);
   |   
   | • python
   |   Table = "Table.txt"
   |   Expr = readmap("Expr.map")
   |   Result1 = lookuplinear(Table,Expr)

   ================================================ ====================================================== =============================================
   Result1.map                                      Table.txt                                              Expr.map                                     
   .. image::  ../examples/lookuplinear_Result1.png .. literalinclude:: ../examples/lookuplinear_Table.txt .. image::  ../examples/lookuplinear_Expr.png
   ================================================ ====================================================== =============================================

   | 

#. 
   | • pcrcalc
   |   binding
   |    Result2 = Result2.map;
   |    Table2 = Table2.txt;
   |    Expr2 = Expr2.map;
   |   initial
   |    report Result2 = lookuplinear(Table2,Expr2);
   |   
   | • python
   |   Table2 = "Table2.txt"
   |   Expr2 = readmap("Expr2.map")
   |   Result2 = lookuplinear(Table2,Expr2)

   ================================================ ======================================================= ==============================================
   Result2.map                                      Table2.txt                                              Expr2.map                                     
   .. image::  ../examples/lookuplinear_Result2.png .. literalinclude:: ../examples/lookuplinear_Table2.txt .. image::  ../examples/lookuplinear_Expr2.png
   ================================================ ======================================================= ==============================================

   | 

