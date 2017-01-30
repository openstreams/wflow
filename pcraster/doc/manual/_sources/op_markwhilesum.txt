.. _markwhilesum:

*****************************
markwhilesumle,markwhilesumge
*****************************
.. index::
   single: markwhilesumle
.. index::
   single: markwhilesumge
.. topic:: markwhilesumle,markwhilesumge

    Marks each cell in specified order until the sum reaches a specified limit.

::

   Result = markwhilesumle(expression1, expression2, limit)

::

   Result = markwhilesumge(expression1, expression2, limit)

expression1
  ordinal, spatial

expression2
  scalar, spatial, non spatial

limit
  scalar, non spatial

Result
  boolean, spatial


Operation
=========
Each cell of expression2 is summed in the order specified by the ascending ordinal values in expression1 while the sum of the cells in expression2 is less or equal than limit for markwhilesumle, or greater or equal than limit for markwhilesumge. Result is a boolean map where every cell that matches the condition is marked 1 (True) or 0 (False) otherwise.

Notes
=====
Missing values in expression1 or expression2 will return missing values in Result.
expression1 should have an ascending order with no multiple instances of the same number.

Examples
========

#. 
   | • pcrcalc
   |   binding
   |    Result2 = Result2.map;
   |    Expr1 = Expr1.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result2 = markwhilesumge(Expr1, scalar(Expr), 40);
   |   
   | • python
   |   Expr1 = readmap("Expr1.map")
   |   Expr = readmap("Expr.map")
   |   Result2 = markwhilesumge(Expr1, scalar(Expr), 40)

   ================================================ ============================================== ====================================
   Result2.map                                      Expr1.map                                      Expr.map                            
   .. image::  ../examples/markwhilesum_Result2.png .. image::  ../examples/markwhilesum_Expr1.png .. image::  ../examples/cos_Expr.png
   ================================================ ============================================== ====================================

   | 

#. 
   | • pcrcalc
   |   binding
   |    Result1 = Result1.map;
   |    Expr1 = Expr1.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result1 = markwhilesumle(Expr1, scalar(Expr), 40);
   |   
   | • python
   |   Expr1 = readmap("Expr1.map")
   |   Expr = readmap("Expr.map")
   |   Result1 = markwhilesumle(Expr1, scalar(Expr), 40)

   ================================================ ============================================== ====================================
   Result1.map                                      Expr1.map                                      Expr.map                            
   .. image::  ../examples/markwhilesum_Result1.png .. image::  ../examples/markwhilesum_Expr1.png .. image::  ../examples/cos_Expr.png
   ================================================ ============================================== ====================================

   | 

