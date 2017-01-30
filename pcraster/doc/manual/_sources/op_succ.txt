

.. index::
   single: succ
.. _succ:

****
succ
****
.. topic:: succ

   Ordinal number of the next higher ordinal class

::

  Result = succ(expression)

expression
   spatial
   ordinal

Result
   spatial
   ordinal

Operation
=========


The result of the operation depends on wheter expression has a legend or not. If expression has a legend, the legend determines the domain of expression: the domain consists of the ordinal numbers linked to the classes in the legend. This domain with these ordinal classes are also assigned to Result. Cells on Result may have values in this domain. For each expression cell value the first higher ordinal number which is in the domain is determined. This is assigned to the corresponding cell on Result.   



If expression does not have a legend an ordinal number is assigned to Result which is the ordinal number on expression plus 1, on a cell- by-cell basis.  

Notes
=====


A cell on expression with missing value is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Order 

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result1 = Result1.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result1 = succ(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result1 = succ(Expr)

   ======================================== =====================================
   Result1.map                              Expr.map                             
   .. image::  ../examples/succ_Result1.png .. image::  ../examples/succ_Expr.png
   ======================================== =====================================

   | 

