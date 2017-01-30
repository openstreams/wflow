

.. index::
   single: argorderarealimited
.. index::
   single: argorderwithidarealimited
.. _argorderarealimited:

*********************************************
argorderarealimited,argorderwithidarealimited
*********************************************
.. topic:: argorderarealimited

   identify highest value by argument order with a limit per argument

::

  Result = argorderarealimited(chances, areaLimit, chances, areaLimit, chances, areaLimit)

::

  Result = argorderwithidarealimited(chances, id, areaLimit, chances, id, areaLimit, chances, id, areaLimit)

areaLimit
   nonspatial
   scalar

chances
   spatial
   scalar

id
   nonspatial
   ordinal

Result
   spatial
   ordinal

Operation
=========


Assign to Result the argument number of  highest chances value. If the assignments for  argument i equals areaLimit_1, argument i can no longer be assigned and the next highest is assigned. If all assignments are exhausted a 0 value is assigned. For the argorderarealimited function this argument number is a value between 1 and n. The argorderwithidarealimited function will assign the value of id_i when i is assigned. 

Notes
=====


A cell with missing value on chances is not considered in the operation; it is assigned a missing value on Result. An id_i value of 0 will lead to unclear results.  

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
   |    report Result = argorderarealimited(
   |    Chances11,
   |     2, # limit in pixels
   |    Chances12,
   |     2  # limit in pixels
   |     )
   |     ;
   |   
   | • python
   |   Chances11 = readmap("Chances11.map")
   |   Chances12 = readmap("Chances12.map")
   |   Result = argorderarealimited(
   |    Chances11,
   |     2, # limit in pixels
   |    Chances12,
   |     2  # limit in pixels
   |     )
   |     

   ====================================================== =============================================================== ===============================================================
   Result.map                                             Chances11.map                                                   Chances12.map                                                  
   .. image::  ../examples/argorderarealimited_Result.png .. image::  ../examples/argorderwithidarealimited_Chances11.png .. image::  ../examples/argorderwithidarealimited_Chances12.png
   ====================================================== =============================================================== ===============================================================

   | 

