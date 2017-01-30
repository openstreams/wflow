

.. index::
   single: downstream
.. _downstream:

**********
downstream
**********
.. topic:: downstream

   Cell gets value of the neighbouring downstream cell

::

  Result = downstream(ldd, expression)

expression
   spatial
   boolean, nominal, ordinal, scalar, directional, ldd

ldd
   spatial
   ldd

Result
   spatial
   type of expression

Operation
=========


For each cell, assigns to Result the expression value of the neighbouring downstream cell, where downstream cells are determined using the local drain directions on ldd. In case a cell doesn't have a downstream cell (i.e. a pit) its own expression value is assigned to Result.  

Notes
=====


A cell with a missing value on ldd and/or expression is assigned a missing value on Result. Its upstream neighbours are assigned a missing value also.  

Group
=====
This operation belongs to the group of  Neighbourhood operators; local drain directions 

See Also
========
:ref:`lddmask`

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Ldd = Ldd.map;
   |    Expr = Expr.map;
   |   initial
   |    report Result = downstream( Ldd, Expr);
   |   
   | • python
   |   Ldd = readmap("Ldd.map")
   |   Expr = readmap("Expr.map")
   |   Result = downstream( Ldd, Expr)

   ============================================= ==================================== ===========================================
   Result.map                                    Ldd.map                              Expr.map                                   
   .. image::  ../examples/downstream_Result.png .. image::  ../examples/accu_Ldd.png .. image::  ../examples/downstream_Expr.png
   ============================================= ==================================== ===========================================

   | 

