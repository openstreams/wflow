The wflow_floodmap model
========================


Introduction
------------
The wflow\_floodmap module can generate flood maps from output
of a wflow\_sbm|hbv| or wflow\_routing model.

At the moment there are two approaches for flood mapping
1. wflow_floodmap.py -- this is a regular wflow-type model, running at the same
resolution as the wflow model used to establish flood maps. The benefit is that
it produces dynamic flood maps. The down side is that the results are in the
same resolution as the original model. The method is also not volume conservative
as it only does a planar inundation and bound to lead to large overestimations
in very flat areas.

2. wflow_flood.py (see Scripts folder). This is a postprocessor to results of a wflow
model and transforms low-resolution model results into a high-resolution flood map
using a (possibly much) higher resolution terrain dataset as input. We recommend
to retrieve high resolution terrain from the Bare-Earth SRTM dataset by Bristol University.
See https://data.bris.ac.uk/data/dataset/10tv0p32gizt01nh9edcjzd6wa

Method
------
PM

Configuration
-------------
PM


	


Description of the wflow_floodmap model
---------------------------------------

.. automodule:: wflow_floodmap
	:members:
    
Description of the wflow_flood post processor
---------------------------------------------

.. automodule:: wflow_flood
	:members:










    


