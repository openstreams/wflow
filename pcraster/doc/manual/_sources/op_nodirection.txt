

.. index::
   single: nodirection
.. _nodirection:

***********
nodirection
***********
.. topic:: nodirection

   Expression of directional data type

::

  Result = nodirection(expression)

expression
   spatial, non spatial
   directional

Result
   dimension of expression
   boolean

Operation
=========


The expression cell values represent directions, with values equal to 0 or between 0 and 360 (or 2 :emphasis:`pi` if the global option :literal:`--radians` is set). A cell with no direction (a flat area) has a value -1. Result is a Boolean expression where 1 is TRUE for cells with no direction and 0 is FALSE for cells with a direction.   



Each cell without a direction on expression is assigned a 1 (Boolean TRUE) on Result. Cells that have a direction on expression are assigned a 0 (Boolean FALSE) on Result.  

Notes
=====


Group
=====
This operation belongs to the group of  Missing value creation 

See Also
========
:ref:`aspect`, :ref:`slope`

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result = nodirection(Expr);
   |   
   | • python
   |   Expr = readmap("Expr.map")
   |   Result = nodirection(Expr)

   ============================================== ============================================
   Result.map                                     Expr.map                                    
   .. image::  ../examples/nodirection_Result.png .. image::  ../examples/nodirection_Expr.png
   ============================================== ============================================

   | 

