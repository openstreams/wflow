

.. index::
   single: lddcreatedem
.. _lddcreatedem:

************
lddcreatedem
************
.. topic:: lddcreatedem

   Modified digital elevation model

::

  Result = lddcreatedem(elevation, outflowdepth, corevolume, corearea, catchmentprecipitation)

catchmentprecipitation
   spatial, non spatial
   scalar

corearea
   spatial, non spatial
   scalar

elevation
   spatial
   scalar

outflowdepth
   spatial, non spatial
   scalar

corevolume
   spatial, non spatial
   scalar

Result
   spatial
   scalar

Options
=======

pit removing at edges of the map

:literal:`--lddout`
   small catchments at the edge of the map are not considered as potentially being affected by the pit removing process: pits which are at the edge of the map are not removed. These pits remain in the map as outflow points of these small catchmentsn(default).

:literal:`--lddin`
   pits at the edge of the map (outflow points of a catchment) are removed like the other pits (if they have core dimensions smaller than the pit dimension thresholds).nOn the result their original catchment cells (including the pit cell) will drainnin another catchment.



:literal:`--unittrue` (default) or :literal:`--unitcell`

:literal:`--unittrue`
   elevation, outflowdepth and catchmentprecipitation is measured in true length, corearea in true area and corevolume in true volume. Units used for elevation and horizontal distance in x and y direction must be the same (default).

:literal:`--unitcell`
   elevation, outflowdepth and catchmentprecipitation is measured in number of cell lengths, corearea in number of cells and corevolume in number of 3D blocks with edges of one cell length.




assignment of elevation in pits

:literal:`--lddfill`
   for each pit which is removed, cells in the area which was formerly the corenof the pit are assigned an elevation equal to the overflow level of the pit coren(default).

:literal:`--lddcut`
   for each pit which is removed, cells on the path between the pit and the outflowncell are assigned the elevation of the pit cell. The elevation of the other cellsnin the core of the pit is not changed.



Operation
=========


This operation corresponds with the local drain direction maker
lddcreate with the difference that a modified digital elevation model is created instead of a local drain direction map. The modified digital elevation model fits  the local drain direction map generated on the basis of the original digital elevation model. 'Not real' cores are removed from the local drain direction map. Additionally an extra option needed for creation of the modified digital elevation model can be specified.   



The expressions used for the pit removing process outflowdepth, corevolume, corearea, catchmentprecipitation and the options - -unitcell/unittrue and :literal:`--lddout`/:literal:`--lddin` have exactly the same meaning and are used in the same way as with the lddcreate operation. So, before you start making modified dem's using lddcreatedem we advise you to read and study the description of the lddcreate command first.   



First, a local drain direction map is generated internally, using
elevation. This is done after the manner of the lddcreate operator, but it is not saved as an expression. Second, the original digital elevation model elevation is modified in such a way that it fits this local drain direction map. This modified digital elevation model is saved as Result. The cell values on Result correspond with the values on elevation, with the exception that the elevation of cells in cores of pit cells is changed. This is done for cores of pit cells which are removed only; the elevation in cores of pits which are not removed remains unaffected. The way elevation values in cores are changed is specified with the :literal:`--lddfill` and :literal:`--lddcut` options. Setting the option :literal:`--lddfill` the elevation of cells in a core is increased until the overflow level is reached. This can be compared with fluviatile or lacustrine sedimentation in the core depression until a maximum sedimentation level is reached: the level of the core pass which is at the lowest elevation. The option :literal:`--lddcut` does not fill the core but reduces the elevation of the cells on the path between the pit cell and the overflow cell. This can be compared with digging a canal in the core between the pit and the pass with the lowest elevation; the canal bottom is at the elevation of the pit.  

Notes
=====


A cell with missing value on one or more of the input expressions is
totally ignored during operation of lddcreatedem; it is assigned a missing value on Result.  



Here, a somewhat generalized description of pit removing and reversal of local
drain directions is given. For a detailed description see Van Deursen, 1995.



Group
=====
This operation belongs to the group of  Derivatives of elevation maps 

Examples
========
#. 
   | • pcrcalc
   |   #! --lddcut
   |   binding
   |    Result2 = Result2.map;
   |    Dem = Dem.map;
   |   initial
   |    report Result2 = lddcreatedem(Dem,999999,9999999,9999999,9999999);
   |   
   | • python
   |   setglobaloption("lddcut")
   |   Dem = readmap("Dem.map")
   |   Result2 = lddcreatedem(Dem,999999,9999999,9999999,9999999)

   ================================================ =========================================
   Result2.map                                      Dem.map                                  
   .. image::  ../examples/lddcreatedem_Result2.png .. image::  ../examples/lddcreate_Dem.png
   ================================================ =========================================

   | 

#. 
   | • pcrcalc
   |   binding
   |    Result1 = Result1.map;
   |    Dem = Dem.map;
   |   initial
   |    report Result1 = lddcreatedem(Dem,9999999,9999999,9999999,9999999);
   |   
   | • python
   |   Dem = readmap("Dem.map")
   |   Result1 = lddcreatedem(Dem,9999999,9999999,9999999,9999999)

   ================================================ =========================================
   Result1.map                                      Dem.map                                  
   .. image::  ../examples/lddcreatedem_Result1.png .. image::  ../examples/lddcreate_Dem.png
   ================================================ =========================================

   | 

