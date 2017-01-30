

.. index::
   single: argorder
.. index::
   single: argorderwithid
.. _argorder:

***********************
argorder,argorderwithid
***********************
.. topic:: argorder

   identify highest value by argument order

::

  Result = argorder(chances, chances, chances)

::

  Result = argorderwithid(chances, id, chances, id, chances, id)

id
   nonspatial
   ordinal

chances
   spatial
   scalar

Result
   spatial
   ordinal

Operation
=========


Assign to Result the argument number of  highest chances value. For the argorder function this argument number is a value between 1 and n. The argorderwithid function will assign the value of id_i when chances_i is the highest value. 

Notes
=====


A cell with missing value on chances  is not considered in the operation; it is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Order 

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Chances11 = Chances11.map;
   |    Chances12 = Chances12.map;
   |   initial
   |    report Result = argorder(
   |    Chances11,
   |    Chances12);
   |   
   | • python
   |   Chances11 = readmap("Chances11.map")
   |   Chances12 = readmap("Chances12.map")
   |   Result = argorder(
   |    Chances11,
   |    Chances12)

   =========================================== =============================================================== ===============================================================
   Result.map                                  Chances11.map                                                   Chances12.map                                                  
   .. image::  ../examples/argorder_Result.png .. image::  ../examples/argorderwithidarealimited_Chances11.png .. image::  ../examples/argorderwithidarealimited_Chances12.png
   =========================================== =============================================================== ===============================================================

   | 

