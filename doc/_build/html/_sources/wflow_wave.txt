The wflow_wave Model
====================


.. warning::

    The documentation of this model is incomplete

Introduction
------------
An *experimental* implementation of the full dynamic wave equations has been 
implemented. The current implementation is fairly unstable and *very* slow.

However, in flat of tidal areas and areas that flood the dynamic wave can provide
much better results. The plot below is from the Rio Mamore in Bolivia in the
lower partes of the river with extensive wetlands that flood nearly each year.


.. plot:: plots/kin-dyn.py


Dependencies
------------
This module is setup to be run in *an existing case and runid* of a wflow\_sbm
or wflow\_hbv model. In order for wflow\_wave to run they must have saved
discharge and waterlevel for each timesteps. This output will be used as a forcing
for the dynamic wave module. The wflow\_wave module will also use the existing ldd
and DEM


Configuration
-------------


It needs anumber of settings in the ini file. The default name for the file
is wflow\_wave.ini. it is also possible to insert this section in the
wflow\_sbm or wflow\_hbv ini file and point to that file.

See below for an example: 

::

    [inputmapstacks]
    # Name of the mapstack with discharge (output from the hydrological model)
    Q = run
    # Name of the mapstack with waterlevel (output from the hydrological model)
    H = lev

    [dynamicwave]
    
    # Number of timeslices per dynamic wave substep
    TsliceDyn=100
    
    # number of substeps for the dynamic wave with respect to the model timesteps
    dynsubsteps=24
    
    # map with level boundary points
    wflow_hboun = staticmaps/wflow_outlet.map
    
    # Optional river map for the dynamic wave that must be the same size or smaller as that of the
    # kinematic wave
    wflow_dynriver = staticmaps/wflow_dynriver.map
    
    # a fixed water level for each non-zero point in the wflow_hboun map 
    # level > 0.0 use that level
    # level == 0.0 use supplied timeseries (see levelTss)
    # level < 0.0 use upstream water level
    fixedLevel = 3.0

    # if this is set the program will try to keep the volume at the pits at
    # a constant value
    lowerflowbound = 1
    
    # instead of a fixed level a tss file with levels for each timesteps and each 
    # non-zero value in the wflow_hboun map
    #levelTss=intss/Hboun.tss
    
    # If set to 1 the program will try to optimise the timestep
    # Experimental, mintimestep is the smallest to be used
    AdaptiveTimeStepping = 1
    mintimestep =0.1


A description of the implementation of the dynamicwave is given on the 
`pcraster website <http://pcraster.geo.uu.nl/documentation/PCRaster/html/op_dynamicwave.html>`_.
 
In addition to the settings in the ini file you need to give the model additional maps
or lookuptables in the staticmaps or intbl directories:


- ChannelDepth.[map|tbl] - depth of the main channel in metres
- ChannelRoughness.[map|tbl] - 1/n n = manning roughness coefficient (default = 1/0.03)
- ChannelForm.[map|tbl] - form of the channel (default = 1.0)
- FloodplainWidth.[map|tbl] - width of the floodplain in metres (default = 0)


The following variables for the dynamicwave function are set as follows and are taken from the
hydrological model run by default:

- ChannelBottomLevel - Taken from the dem (wflow_dem.map in the staticmaps dir)
- ChannelLength - Taken from the length in the kinematic wave (DCL.map from the outsum dir)
- ChannelBottomWidth - taken from the wflow_riverwidth map from the outsum dir





wflow_wave module documentation
-------------------------------

.. automodule:: wflow_wave
    :members:
    :undoc-members:
    :show-inheritance:
