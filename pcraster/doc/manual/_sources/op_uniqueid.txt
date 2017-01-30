

.. index::
   single: uniqueid
.. _uniqueid:

********
uniqueid
********
.. topic:: uniqueid

   Unique whole value for each Boolean TRUE cell

::

  Result = uniqueid(expression)

expression
   spatial, non spatial
   boolean

Result
   spatial
   scalar

Operation
=========


For each cell that has a value 1 (TRUE) on expression assigns a unique whole positive value to Result, starting with 1. Cells that have a value 0 (FALSE) on expression are assigned a value 0.  

Notes
=====


A cell with a missing value on expression is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Coordinates, unique ID's 

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result = uniqueid(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = uniqueid(Expr)

   =========================================== =========================================
   Result.map                                  Expr.map                                 
   .. image::  ../examples/uniqueid_Result.png .. image::  ../examples/uniqueid_Expr.png
   =========================================== =========================================

   | 

