===========================================
Welcome to wflow's documentation!
===========================================

.. note::

      This documentation was generated |today|

      Latest version documentation (development):

      http://wflow.readthedocs.org/en/latest/

      Latest release (stable) version documentation

      http://wflow.readthedocs.org/en/stable/


.. note::

    wflow is released under version 3 of the GPL

    wflow uses PCRaster/Python (see http://www.pcraster.eu) as it's calculation engine.


Introduction
============

This document describes the wflow distributed hydrological modelling platform.
wflow is part of the Deltares'
OpenStreams project (http://www.openstreams.nl). Wflow consists of a
set of python programs that can be run on the command line and perform
hydrological simulations. The models are based on the PCRaster python
framework (www.pcraster.eu). In wflow this framework is extended (the ``wf_DynamicFramework``
class) so that models build using the framework can be controlled using
the API. Links to BMI and OpenDA (www.openda.org)
have been established. All code is available at github (https://github.com/openstreams/wflow/)
and distributed under the GPL version 3.0.

The  wflow distributed hydrological model platform currently includes
the following models:

-  the wflow\_sbm  model (derived from `topog\_sbm <https://csdms.colorado.edu/wiki/Model:TOPOG>`_ )

-  the wflow\_hbv model (a distributed version of the HBV96 model).

-  the wflow\_gr4 model (a distributed version of the gr4h/d models).

-  the wflow\_W3RA and wflow\_w3 models (implementations and adaptations of the Australian Water Resources Assessment Landscape model (AWRA-L))

-  the wflow\_topoflex model (a distributed version of the FLEX-Topo model)

-  the wflow\_pcrglobwb model (`PCR-GLOBWB <http://globalhydrology.nl/models/pcr-globwb-2-0/>`_ (PCRaster Global Water Balance, v2.1.0_beta_1))

-  the wflow\_sphy model (`SPHY <https://www.futurewater.eu/methods/sphy/>`_ (Spatial Processes in HYdrology, version 2.1))

-  the wflow\_stream model (`STREAM <https://www.ivm.vu.nl/en/projects/Projects/spatial-analysis/stream/index.aspx>`_ 
   (Spatial Tools for River Basins and Environment and Analysis of Management Options))

-  the wflow\_routing model (a kinematic wave model that can run on the output of one of the hydrological models
   optionally including a floodplain for more realistic simulations in areas that flood).

-  the wflow\_wave model (a dynamic wave model that can run on the output of the wflow\_routing model).

-  the wflow\_floodmap model (a flood mapping model that can use the output of the wflow\_wave model or de wflow\_routing model).

-  the wflow\_sediment model (an experimental erosion and sediment dynamics model that uses the output of the wflow\_sbm model).

-  the wflow\_lintul model (rice crop growth model `LINTUL <https://edepot.wur.nl/461276>`_ (Light Interception and Utilization))

The low level api and links to other frameworks allow the models to be
linked as part of larger modelling systems:


.. digraph:: Linking

    WFLOW_HBV -> WFLOWAPI;
    WFLOW_SBM -> WFLOWAPI;
    WFLOWAPI -> "PI"  [dir=both];
    "Data and Models" -> "PI";
    WFLOWAPI -> OpenMI  [dir=both];
    ModelX -> OpenMI;
    ModelY -> OpenMI;
    ModelY -> BMI;
    BMI -> OpenDA  [dir=both];
    WFLOWAPI -> BMI  [dir=both];
    Calibration -> OpenDA;
    Assimilation -> OpenDA;
    WFLOWAPI [shape=square];
    OpenDA [shape=square];
    OpenMI [shape=square];
    BMI [shape=square];
    "PI" [shape=square];
    dpi=69;


.. note::

    wflow is part of the Deltares OpenStreams project
    (http://www.openstreams.nl). The OpenStreams project is a work in
    progress. Wflow functions as a toolkit for distributed hydrological
    models within OpenStreams.


.. note::

   As part of the eartH2Observe project global dataset of forcing data has been compiled that can also be used with the
   wflow models. A set of tools is available that can work with wflow (the wflow_dem.map file) to extract data from the server and downscale
   these for your wflow model. Check https://github.com/earth2observe/downscaling-tools for the tools. A description
   of the project can be found at http://www.earth2observe.eu and the data server can be access via http://wci.earth2observe.eu


The different wflow models share the same structure but are fairly
different with respect to the conceptualisation. The shared software
framework includes the basic maps (dem, landuse, soil etc) and the
hydrological routing via the kinematic wave. The Python class framework
also exposes the models as an API and is based on PCRaster/Python.

The wflow\_sbm model maximises the use of available spatial data.
Soil depth, for example, is estimated from the DEM using a topographic
wetness index . The model is derived from the CQflow model (Köhler et al., 2006) that has
been applied in various countries, most notably in Central America. The
wflow\_hbv model is derived from the HBV-96 model but does not
include the routing functions, instead it uses the same kinematic wave
routine as the wflow\_sbm  model to route the water downstream.

The models are programmed in Python using the PCRaster Python extension.
As such, the structure of the model is
transparent, can be changed by other modellers easily, and the system
allows for rapid development.

.. toctree::
   :hidden:
   :maxdepth: 2

   installation.rst
   wflow_usage.rst
   wflow_building.rst
   faq.rst
   index_models.rst
   wf_DynamicFramework.rst
   wflow_adapt.rst
   index_modules_libs.rst
   index_bmi.rst
   wflow_modelbuilder.rst
   Release notes <release_notes.rst>
   OpenDA <wflow_openda.rst>



References
==========

- Köhler, L., Mulligan, M., Schellekens, J., Schmid, S. and Tobón, C.: Final Technical Report DFID-FRP Project no. R7991 Hydrological
  impacts of converting tropical montane cloud forest to pasture, withinitial reference to northern Costa Rica. 2006.


Papers/reports using wflow
==========================

Arnal, L., 2014. An intercomparison of flood forecasting models for the Meuse River basin (MSc Thesis). Vrije Universiteit,
Amsterdam.

Azadeh Karami Fard, 2015. Modeling runoff of an Ethiopian catchment with WFLOW (MSc thesis). Vrije Universiteit,
Amsterdam.

de Boer-Euser, T., Bouaziz, L., De Niel, J., Brauer, C., Dewals, B., Drogue, G., Fenicia, F., Grelier, B., Nossent, J.,
Pereira, F., Savenije, H., Thirel, G., Willems, P., 2017. Looking beyond general metrics for model comparison – lessons
from an international model intercomparison study. Hydrol. Earth Syst. Sci. 21, 423–440. doi:10.5194/hess-21-423-2017

Emerton, R.E., Stephens, E.M., Pappenberger, F., Pagano, T.C., Weerts, A.H., Wood, A.W., Salamon, P., Brown, J.D.,
Hjerdt, N., Donnelly, C., Baugh, C.A., Cloke, H.L., 2016. Continental and global scale flood forecasting systems. WIREs
Water 3, 391–418. doi:10.1002/wat2.1137

Hally, A., Caumont, O., Garrote, L., Richard, E., Weerts, A., Delogu, F., Fiori, E., Rebora, N., Parodi, A., Mihalović,
A., Ivković, M., Dekić, L., van Verseveld, W., Nuissier, O., Ducrocq, V., D’Agostino, D., Galizia, A., Danovaro, E.,
Clematis, A., 2015. Hydrometeorological multi-model ensemble simulations of the 4 November 2011 flash flood event in
Genoa, Italy, in the framework of the DRIHM project. Nat. Hazards Earth Syst. Sci. 15, 537–555.
doi:10.5194/nhess-15-537-2015

Hassaballah, K., Mohamed, Y., Uhlenbrook, S., Biro, K., 2017. Analysis of streamflow response to land use land cover
changes using satellite data and hydrological modelling: case study of Dinder and Rahad tributaries of the Blue Nile.
Hydrol. Earth Syst. Sci. Discuss. 2017, 1–22. doi:10.5194/hess-2017-128

Jeuken, A., Bouaziz, L., Corzo, G., Alfonso, L., 2016. Analyzing Needs for Climate Change Adaptation in the Magdalena
River Basin in Colombia, in: Filho, W.L., Musa, H., Cavan, G., O’Hare, P., Seixas, J. (Eds.), Climate Change Adaptation,
Resilience and Hazards, Climate Change Management. Springer International Publishing, pp. 329–344.

López López, P., Wanders, N., Schellekens, J., Renzullo, L.J., Sutanudjaja, E.H., Bierkens, M.F.P., 2016. Improved
large-scale hydrological modelling through the assimilation of streamflow and downscaled satellite soil moisture
observations. Hydrol. Earth Syst. Sci. 20, 3059–3076. doi:10.5194/hess-20-3059-2016

Maat, W.H., 2015. Simulating discharges and forecasting floods using a conceptual rainfall-runoff model for the Bolivian
Mamoré basin (MSc thesis). University of Twente, Enschede.

Research paper: HYDROLOGIC MODELING OF PRINCIPAL SUB-BASINS OF THE MAGDALENA-CAUCA LARGE BASIN USING WFLOW MODEL [WWW
Document], n.d. . ResearchGate. URL
https://www.researchgate.net/publication/280293861_HYDROLOGIC_MODELING_OF_PRINCIPAL_SUB-BASINS_OF_THE_MAGDALENA-CAUCA_LARGE_BASIN_USING_WFLOW_MODEL
(accessed 4.4.17).

Tangdamrongsub, N., Steele-Dunne, S.C., Gunter, B.C., Ditmar, P.G., Weerts, A.H., 2015. Data assimilation of GRACE
terrestrial water storage estimates into a regional hydrological model of the Rhine River basin. Hydrol. Earth Syst.
Sci. 19, 2079–2100. doi:10.5194/hess-19-2079-2015

Tretjakova, D., 2015. Investigating the effect of using fully-distributed model and data assimilation on the performance
of hydrological forecasting in the Karasu catchment, Turkey (MSc thesis). Wageningen University.

Wang, X., Zhang, J., Babovic, V., 2016. Improving real-time forecasting of water quality indicators with combination of
process-based models and data assimilation technique. Ecological Indicators 66, 428–439.
doi:10.1016/j.ecolind.2016.02.016



TODO
====

.. todolist::