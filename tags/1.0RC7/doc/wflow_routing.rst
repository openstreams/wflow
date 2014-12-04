THe wflow_routing Model
=======================


.. warning::

    The documentation of this model is incomplete

Introduction
------------
A


Dependencies
------------
T


Configuration
-------------


It needs a number of settings in the ini file. The default name for the file
is wflow\_routing.ini. it is also possible to insert this section in the
wflow\_sbm or wflow\_hbv ini file and point to that file.

See below for an example: 

::

    [inputmapstacks]
    # Name of the mapstack with specific discharge (mm/timestep output from the hydrological model)
    IW= inmaps/IW





A description of the implementation of the kinematic wave is given on the pcraster website
 
In addition to the settings in the ini file you need to give the model additional maps
or lookuptables in the staticmaps or intbl directories:





wflow_routing module documentation
----------------------------------

.. automodule:: wflow_routing
    :members:
    :undoc-members:
    :show-inheritance:
