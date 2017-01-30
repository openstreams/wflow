

.. index::
   single: argorderaddarealimited
.. index::
   single: argorderwithidarealimited

.. _argorderaddarealimited:

***************************************************
argorderaddarealimited,argorderwithidaddarealimited
***************************************************
.. topic:: argorderaddarealimited

   variation on argorder

::

  Result = argorderaddarealimited(currentId, chances, areaLimit, chances, areaLimit, chances, areaLimit)

::

  Result = argorderwithidaddarealimited(currentId, chances, id, areaLimit, chances, id, areaLimit, chances, id, areaLimit)

chances
   spatial
   scalar

currentId
   ordinal
   spatial

areaLimit
   nonspatial
   scalar

id
   nonspatial
   ordinal

Result
   spatial
   ordinal

Operation
=========

Operation is like the :ref:`argorderarealimited` variant except that the areaLimit argument has been replaced by the areaAdded argument.
The operation will extend the area with the size of areadAdded on base of "maximum likelyhood".


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
   |    CurrentId = CurrentId.map;
   |    Chances1 = Chances1.map;
   |    Chances2 = Chances2.map;
   |   initial
   |    report Result = argorderaddarealimited(
   |    CurrentId,
   |    Chances1,
   |     1,# add area
   |    Chances2,
   |     4  # limit in pixels
   |     )
   |     ;
   |   
   | • python
   |   CurrentId = readmap("CurrentId.map")
   |   Chances1 = readmap("Chances1.map")
   |   Chances2 = readmap("Chances2.map")
   |   Result = argorderaddarealimited(
   |    CurrentId,
   |    Chances1,
   |     1,# add area
   |    Chances2,
   |     4  # limit in pixels
   |     )
   |     

   ========================================================= ============================================================ =========================================================== ===========================================================
   Result.map                                                CurrentId.map                                                Chances1.map                                                Chances2.map                                               
   .. image::  ../examples/argorderaddarealimited_Result.png .. image::  ../examples/argorderaddarealimited_CurrentId.png .. image::  ../examples/argorderaddarealimited_Chances1.png .. image::  ../examples/argorderaddarealimited_Chances2.png
   ========================================================= ============================================================ =========================================================== ===========================================================

   | 

