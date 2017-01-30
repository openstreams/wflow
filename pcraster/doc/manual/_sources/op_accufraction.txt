

.. index::
   single: accufractionstate
.. index::
   single: accufractionflux
.. _accufraction:

***********************************
accufractionflux, accufractionstate
***********************************
.. topic:: accufractionflux, accufractionstate

   Fractional material transport downstream over local drain direction network

::

  Resultflux = accufractionflux(ldd, material, transportfraction)

::

  Resultstate = accufractionstate(ldd, material, transportfraction)

::

  Resultflux Resultstate = accufractionflux, accufractionstate(ldd, material, transportfraction)

transportfraction
   spatial, non spatial
   scalar

ldd
   spatial
   ldd

material
   spatial, non spatial
   scalar

Resultflux
   spatial
   scalar

Resultstate
   spatial
   scalar

Operation
=========


These operations describe the accumulation of material in a drainage network
with transport of a certain fraction. The remaining material is withdrawn
from the stream. The operators enable the description of phenomena such
as loss of a certain percentage of organic matter over a river stretch.  



For each cell, accufractionflux assigns the amount of material which is transported out of the cell, accufractionstate assigns the amount which is stored in the cell.  Both perform the same function of accumulation of material with a transport fraction, the only difference between the operators is the sort of result that is saved: accufractionstate yields storages of material in cells, accufractionflux yields fluxes of material out of cells.   



For each cell, the amount of material input, for instance the amount of
rain, is given by material. This is transported in downstream direction through the consecutively neighbouring downstream cells, following the local drain directions on ldd. Each time material moves through a cell a certain amount is stored in the cell. These storages are saved as Resultstate, if the accufractionstate operator is used. The remaining material is transported out of the cell, these amounts of outflow from each cell into its neighbouring downstream cell are the result of the accufractionflux operator, they are saved as Resultflux.   



The function can be described by flow of material through a set of linked
systems, where a cell represents a system. The flow starts at the
cells/systems at the watershed boundaries (defined by ldd) and ends at a pit cell. The systems are linked by the local drain directions on ldd; these define the path of flow through the set of cells/systems. Each time a system is passed, the amount of flow changes.   



For a cell/system somewhere on the map, the flow of material is described
by a system. The inflow of the cell is the
sum of the outflow amounts of its upstream neighbours. This inflow
amount is added to the material value in the cell itself. This amount of material is potentially available for transport out of the cell. The amount actually transported is this amount multiplied by the transportfraction value of the cell. The remaining material is stored in the cell. Since transportfraction is a fraction it must contain values equal to or between 0 and 1 ([0,1]). (If transportfraction is 0 nothing will be transported, if it is 1 all material will be transported).   



For each cell, the amount of material which is transported to its
downstream neighbour (or out of the map if the cell is a pit cell) is saved
as Resultflux (use the operator accufractionflux); the amount of material which is stored to the cell is saved as Resultstate (use accufractionstate)  

Notes
=====


The values on material must be equal to or larger than zero. The values on transportfraction must be equal to or between 0 and 1.   



A cell with missing value on material and/or transportfraction is assigned a missing value on Resultflux or Resultstate. Additionally, all its downstream cells are assigned a missing value. The local drain direction network on ldd must be :ref:`sound <SoundLDD>`.  

Group
=====
This operation belongs to the group of  Neighbourhood operators; local drain directions 

See Also
========
:ref:`secstatneightr`, :ref:`lddmask`

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    State2 = State2.map;
   |    Flux2 = Flux2.map;
   |    Ldd = Ldd.map;
   |    Material = Material.map;
   |    TransFra = TransFra.map;
   |   initial
   |    report State2, Flux2 = accufractionstate,accufractionflux(Ldd,Material,TransFra);
   |   
   | • python
   |   Ldd = readmap("Ldd.map")
   |   Material = readmap("Material.map")
   |   TransFra = readmap("TransFra.map")
   |   
   |   State2 = accufractionstate(Ldd,Material,TransFra)
   |   Flux2 = accufractionflux(Ldd,Material,TransFra)

   =============================================== ============================================== ==================================== ================================================= =================================================
   State2.map                                      Flux2.map                                      Ldd.map                              Material.map                                      TransFra.map                                     
   .. image::  ../examples/accufraction_State2.png .. image::  ../examples/accufraction_Flux2.png .. image::  ../examples/accu_Ldd.png .. image::  ../examples/accufraction_Material.png .. image::  ../examples/accufraction_TransFra.png
   =============================================== ============================================== ==================================== ================================================= =================================================

   | 

#. 
   | • pcrcalc
   |   binding
   |    State1 = State1.map;
   |    Flux1 = Flux1.map;
   |    Ldd = Ldd.map;
   |    Material = Material.map;
   |   initial
   |    report State1, Flux1 = accufractionstate,accufractionflux(Ldd,Material,0.5);
   |   
   | • python
   |   Ldd = readmap("Ldd.map")
   |   Material = readmap("Material.map")
   |   
   |   State1 = accufractionstate(Ldd,Material,0.5)
   |   Flux1 = accufractionflux(Ldd,Material,0.5)
   |   

   =============================================== ============================================== ==================================== =========================================
   State1.map                                      Flux1.map                                      Ldd.map                              Material.map                             
   .. image::  ../examples/accufraction_State1.png .. image::  ../examples/accufraction_Flux1.png .. image::  ../examples/accu_Ldd.png .. image::  ../examples/accu_Material.png
   =============================================== ============================================== ==================================== =========================================

   | 

