

.. index::
   single: streamorder
.. _streamorder:

***********
streamorder
***********
.. topic:: streamorder

   Stream order index of all cells on a local drain direction network

::

  Result = streamorder(ldd)

ldd
   spatial
   ldd

Result
   spatial
   ordinal

Operation
=========


The classification of stream networks was originally developed by Horton,
and modified by :ref:`Strahler <strahler64>`. Following the scheme of Strahler, streamorder designates an order 1 to the smallest channels, which are the cells with no upstream cells connected to that cell. Where two channels of order 1 join, a channel of order 2 results downstream. In general, where two channels of order i join, a channel of order i+1 results.   

Group
=====
This operation belongs to the group of  Neighbourhood operators; local drain directions 

See Also
========
:ref:`accuflux`

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Ldd = Ldd.map;
   |   initial
   |    report Result = streamorder(Ldd);
   |   
   | • python
   |   Ldd = readmap("Ldd.map")
   |   Result = streamorder(Ldd)

   ============================================== ====================================
   Result.map                                     Ldd.map                             
   .. image::  ../examples/streamorder_Result.png .. image::  ../examples/accu_Ldd.png
   ============================================== ====================================

   | 

