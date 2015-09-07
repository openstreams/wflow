The wflow_topoflex Model
========================


Introduction
------------

Topoflex applies the concept of using different model configurations for different hydrological response units (HRUs)
. These HRUs can for a large part be derived from topography ([savenije]); however, additional data such as land use
and geology can further optimize the selection of hydrological response units. In contrast to other models, topoflex
generally uses a minimum amount of HRUs, which are defined manually, i.e. 2-5 depending on the size of and the
variablity within the  catchment. Nevertheless, the model code is written such that it can handle an infinte number
of classes. The individual HRUs are treated as parallel model structures and only connected via the groundwater
reservoir and the stream network.

The model code is written in a modular setup: different conceptualisations can be selected for each reservoir for
each HRU. The code currently contains some reservoir conceptualisations (see below), but new conceptualisations can
easily be added.

Examples of the application of topoflex can be found in [gharari], [gao] and [euser].

Figure 1 shows a possible model conceptualistion: one that used two HRUs (wetland (W) and hillslope (H), adapted from [euser])

.. figure:: _images/topoflex.jpg
   :width: 600px

   Schematic view of the relevant components of the topoflex model

Limitations
~~~~~~~~~~~

- Using a set of HRUs introduces additional complexity (structural and computational) in the model. Therefore,
calibration of the model should be carried out carefully (see below for some tips and tricks that might help).

- The selection and deliniation of the HRUs is a relatively subjective exercise: different data sources and
preferably some expert knowledge migth help to construct meaningful HRUs.


Spatial resolution
------------------

The model uses a grid based on the available input or required output data. The combination with the HRUs has to be
done in the preparation step: for each cell the percentage of each class needs to be determined and stored in a
staticmap.

The cell sizes do not have to have the same size: the size of individual cells can be stored in a staticmap, which is
 used to calculate the contribution of each cell to the generated discharge.


Input data
----------

The required input data consists of timeseries of precipitation, temperature potential evaporation. In case the
Jarvis equations are used to determine transpiration more input variable are required, such as humidity, radiation,
wind speed and LAI.


Different HRUs
--------------

Currently there are conceptualisations for three main HRUs: wetland, hillslope and (agricultural) plateau. These
conceptualisations are simply a set of reservoirs, which match the perception of different landscape elements (in
western Europe). Dependent on the area or interest the HRUs can be adapted, was well as the set of reservoirs per HRU.

wetland
~~~~~~~

Wetlands are areas close to the stream, where the groundwater levels are assumed to be shallow and assumed to rise
quickly during an event. The dominant runoff process in wetlands is assumed to be saturation overland flow, leading
to quick and sharp responses to rainfall.

hillslope
~~~~~~~~~

Hillslopes are steep areas in an catchment, generally covered with forest. The dominant runoff process is assumed to
be quick subsurface flow: the hillslopes mainly contribute to the discharge during the winter period.

(agricultural) plateau
~~~~~~~~~~~~~~~~~~~~~~

Plateus are flat areas high above the stream, thus with deep ground water levels. Depending on the specific
conditions, the dominant runoff processes are ground water recharge, quick subsurface flow and hortonian overland
flow. The latter is especially important in agricultural areas.

Other modelling options
-----------------------
.. Jarvis equations
.. snow and frozen ground
.. deep infiltration losses

Routing of generated discharge
------------------------------

The routing of generated discharge is based on the average velocity of water through the river, which is currently
set to 1 m/s. For each cell the average distance to the outlet is calculated and multiplied with the selected flow
velocity to determine the delay at the outlet.

There are currently two options to apply this routing:

1. only calculating the delay relevant for the discharge at the outlet
2. calculating the delay (and thus discharge) over the stream network. This option is mainly relevant for calculations with a finer grid

 
Calibrating the  wflow_topoflex model
-------------------------------------

Including more HRUs in the model leads to an increase in parameters. To make it still feasible to calibrate the
model, a set of constraints is introduced: parameter and process constraints. These constraints are assumed relations
 between parameters and fluxes of different HRUs and prevent the selection of unreaslistic parameters. The
 constraints are an important part of the perceptual model, but are not (yet) included in de wflow code. Below some
 examples of constraints are given, more examples of constraints can be found in [gharari], [gao] and [euser].

Parameter constraints
~~~~~~~~~~~~~~~~~~~~~

Parameter constraints are relations between parameters of different HRUs, for example the root zone storage capacity ( S_{u,max}), which is assumed to be larger on hillslopes than in wetlands. As in the latter groundwater levels quicky rise during a rain event, reducing the root zone storage capacity. Parameter constraints are calculated before the model runs.


Process constraints
~~~~~~~~~~~~~~~~~~~

Process constraints are comparable with parameter constraints, but focus on fluxes from different HRUs, for example
the fast response from the wetlands is assumed to be larger than the fast response of the hillslopes in the summer
period. As on the hillslopes generally more storage is required before a runoff generation threshold is exceeded.
Process constraints are calculated after the model runs.


References
----------
.. [euser]     Euser, T., Hrachowitz, M., Winsemius, H. C. and Savenije, H. H. G.: The effect of forcing and landscape distribution on performance and consistency of model structures. Hydrol. Process., doi: 10.1002/hyp.10445, 2015.	
.. [gao]       Gao, H., Hrachowitz, M., Fenicia, F., Gharari, S., and Savenije, H. H. G.: Testing the realism of a topography-driven model (FLEX-Topo) in the nested catchments of the Upper Heihe, China, Hydrol. Earth Syst. Sci., 18, 1895-1915, doi:10.5194/hess-18-1895-2014, 2014.
.. [gharari]   Gharari, S., Hrachowitz, M., Fenicia, F., Gao, H., and Savenije, H. H. G.: Using expert knowledge to increase realism in environmental system models can dramatically reduce the need for calibration, Hydrol. Earth Syst. Sci., 18, 4839-4859, doi:10.5194/hess-18-4839-2014, 2014.
.. [savenije]  Savenije, H. H. G.: HESS Opinions "Topography driven conceptual modelling (FLEX-Topo)", Hydrol. Earth Syst. Sci., 14, 2681-2692, doi:10.5194/hess-14-2681-2010, 2010.

example ini file
----------------
The .ini file below shows the available options

.. literalinclude:: _download/wflow_topoflex_example.ini

An example ini file be found :download:`here. <_download/wflow_topoflex_example.ini>`


wflow_topoflex module documentation
------------------------------

.. automodule:: wflow_topoflex
    :members:
    :undoc-members:
    :show-inheritance:

    .. autoattribute:: usage


