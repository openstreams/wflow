wflow_delwaq Module
===================


The wflow\_delwaq module provides a set of functions to construct a 
delwaq pointer file from a PCRaster local drainage network.
A command-line interface is provide that allows you to 
create a delwaq model that can be linked to a wflow model.

.. figure:: _images/ldd-rhine.png
	:width: 640px


The script sets-up a one-layer model (representing the kinematic wave
reservoir). Water is labeled according to the area and flux where it
enters the kinematic wave reservoir.

For the script to work a run of the wflow model must be available and a
template directory in which the delwaq model is created should also be
available. These are indicated by the -C -R and -D command line
options. The -R and -C options indicate the wflow case and wflow run
directories while the -D option indicates the delwaq template
directory.

The template used is shown below:

::

	debug/
	fixed/
	fixed/B2_numsettings.inc
	fixed/B4_dispersion.inc
	fixed/B4_dispx.inc
	fixed/B9_Hisvar.inc
	fixed/B9_Mapvar.inc
	includes_deltashell/
	includes_flow/
	run.bat
	dlwqlib.dll
	libcta.dll
	libiomp5md.dll
	waq_plugin_wasteload.dll
	delwaq1.exe
	delwaq2.exe
	deltashell.inp

The debug, includes_flow, and includes_deltashell directories are
filled by the script. After that delwaq1.exe and delwaq2.exe programs
may be run (the run.bat file shows how this is done)::

	delwaq1.exe deltashell.inp -np
	delwaq2.exe deltashell.inp
	
the script sets up delwaq such that the result for the wflow gauges locations
are stored in the deltashell.his file.



How the script works
--------------------

The pointer file for delwaq is made using the following information:

#. The wflow\_ldd.map files is used to create the internal flow network, 
   it defines the segments and how water flows between the segments
#. The number of inflows into each segment determines is taken 
   from the sources mapstacks (-S option). Together these sources 
   should include all the water that enters the kinematic wave 
   reservoir. These are the red and green arrows in the figure below
#. The delwaq network is generated for the area define in the
   wflow_catchment map. The included area is define by all cells were the
   catchment id in a cel is larger than 1.


.. figure:: _images/delwaq-exchanges.png
	:width: 640px

	Figure: How exchanges and inflows are connected


Within the includes_flow directory the following files are made:

- volume.dat - volumes (N+1) noseg
- flow.dat - flows (N). Contents is noq
- area.dat - N timesteps. Content is noq
- surface.dat - surface area of the water per segment (N+1), noseq
- length.dat - One timestep only (constant). Content is two times noq

Here nrseg is the number of segments (taken from the non-missing grid
cell in the wflow_ldd.map file) and noq is the number of exchanges which
is calculated as the number of segments plus number the of inflows
(in each segment) times the number of segments

Delwaq expects volumes to be instantanious values at the start of a
timestes while flows are integrated between tow timesteps. For volumes
N+1 timesteps are needed, for flows N timesteps.  The figure below
demonstrates this principle for N=4.

.. digraph:: Flows
	
	node[shape=record,width=.1,height=.1];
	node0 [label="{Time|Volume|Flow integrated}"];
	node1 [label="{T=0|Volume=0|Flow=0 to 1}"];
	node2 [label="{T=1|Volume=1|Flow=1 to 2}"];
	node3 [label="{T=2|Volume=2|Flow=2 to 3}"];
	node4 [label="{T=3|Volume=3|Flow=3 to 4}"];
	node5 [label="{T=4|Volume=4|May be zero}"];
	node1 -> node2
	node2 -> node3
	node3 -> node4
	node4 -> node5

The volume.dat file is filled with N+1 steps of volumes of the wflow
kinematic wave reservoir. To obtain the needed lag between the flows and
the volumes the volumes are taken from the kinematic wave reservoir one
timestep back (OldKinWaveVolume). 

The flow.dat files is filled as follows. For each timestep internal
flows (within the kinematic wave reservoir, i.e. flows from segment to
segment) are written first (blue in the layout above). Next the flow
into each segment are written. Depending on how many inflow types are
given to the script (sources). For one type, one set of flows is
written, if there are two types two sets etc (green and red in the
layout above).


Very simple example:
--------------------

The following very simple example demonstrated how the pointer file is
created.  First the pcraster ldd:

.. figure:: _images/ldd-mini.png

The resulting network consist of 10 points:

.. figure:: _images/ldd-res.png


As can be seen both 9 and 10 are bottom points. The
generated pointer file is shown below::

        ;Written by dw_WritePointer
	;nr of pointers is:  20
         1          3          0          0
         2          4          0          0
         3          5          0          0
         4          6          0          0
         5          7          0          0
         6          8          0          0
         7          9          0          0
         8         10          0          0
         9         -1          0          0
        10         -2          0          0
        -3          1          0          0
        -4          2          0          0
        -5          3          0          0
        -6          4          0          0
        -7          5          0          0
        -8          6          0          0
        -9          7          0          0
       -10          8          0          0
       -11          9          0          0
       -12         10          0          0



Case study for Malaysia and Singapore
-------------------------------------

To estimate load of different nutrients to Johor strait two wflow\_sbm models
have been setup. Next these models where linked to delwaq as follows:

#. A delwaq segment network similar to the wflow D8 ldd was made
#. The volumes in the delwaq segment are taken from the wflow\_sbm 
   kinematic wave volumes
#. For each segment two sources (inflows) are constructed, fast and slow 
   each representing different runoff compartments from the wflow model. 
   Fast represents SOF [#]_, HOF [#]_ and SSSF [#]_ while Slow represent 
   groundwater flow.
#. Next the flow types are combined with the available land-use classes. 
   As such a Luclass times flowtypes matrix of constituents is made. 
   Each constituent (e.g. Slow flow of LU class 1) is traced throughout 
   the system. All constituents are conservative and have a 
   concentration of 1 as they flow in each segement.
#. To check for consistency an Initial water type and a Check water 
   type are introduced. The Initial water will leave the system gradually 
   after a cold start, the Check water type is added to each flow 
   component and should be 1 at each location in the system (Mass Balance 
   Check).
  
.. [#] SOF: Saturation Overland Flow
.. [#] HOF: Hortonian Overland Flow (or infiltration excess Overland Flow)
.. [#] SSSF: SubSurface Storm Flow. Rapid lateral flow through the 
             top part of the soil profile.

The above results in a system in which the different flow types (including
the LU type where they have been generated) can be traced throughout the
system. Each each gauge location the discharge and the flow components
that make up the discharge are reported.

.. figure:: _images/sing_24.png
	:width: 800px
	
	Figure: Discharge and flow types for a small Singapore catchment. The
	Singapore catchment are dominated by fast flow types but during the end
	of the dry periods the slow flow types start to rise in importance.

.. figure:: _images/malay_7.png
	:width: 800px

	Figure: Discharge, flow types and resulting total P for a catchment in 
	Malaysia.

By assuming each flow type is an end-member in a mixing model we can add
fixed concentration of real parameters to the flow fractions and multiply
those with the concentrations of the end-membesrt modelled concentration
at the gauge locations can be obtained for each timestep.

.. figure:: _images/sbm_flowtypes.png

    Figure: Flow types in the topog\_sbm models used in the
    Singapore/Malaysia case. HOF = Hortonian or Infiltration excess
    overland flow. SOF = Saturation overland flow, GW = exfiltrating
    groundwater. Note that the subcell representation of saturated
    areas means that both SOF and GW can occur before a cell is
    completely saturated.

The figure above shows the flow types in the models used in Singapore and 
Malaysia. Groundwater flow (both from completely saturated cell and 
subcell groundwater flow) makes up the Slow flow that is fed into the delwaq 
model while SOF and HOF make up the Fast flow to the delwaq model. In addition
the water is also labelled according to the landuse type of the cell that 
it flows out off.


The whole procedure was setup in a Delft-FEWS configuration that can run 
the following steps operationally:

.. digraph:: Flows
	
    node[shape=record,width=.1,height=.1];
    "Pre-Process P, T and PET data to match the model grid" -> "Run the hydrological Model"
    "Run the hydrological Model" -> " Save all flows per cell"
    " Save all flows per cell" -> "Feed flows per LU type and flow type to delwaq"
    "Feed flows per LU type and flow type to delwaq" -> "Obtain flow fraction per LU and flow type at gauge locations"
    "Obtain flow fraction per LU and flow type at gauge locations" -> "Multiply constituent concentration per LU and flow type with fraction"
    "Multiply constituent concentration per LU and flow type with fraction" -> "Sum all fraction concentrations to obtain total concentration at gauge locations"



wflow\_delwaq module documentation
----------------------------------

.. automodule:: wflow_delwaq
    :members:
    :undoc-members:
