============================
Using the wflow modelbuilder
============================

wflow modelbuilder tutorial
===========================

.. note::

	**UNDER CONSTRUCTION**

The wflow modelbuilder is a new tool with which you can set up a wflow
model in a few simple steps. The default setup of the wflow modelbuilder
is fully based on global data sets. The modelbuilder uses a set of tools
called hydro-engine (https://github.com/openearth/hydro-engine), which
are built on top of Google Earth Engine (https://earthengine.google.com/).

Installation
------------

The wflow modelbuilder is part of the standard wflow distribution.

You can download the wflow source code from
https://github.com/openstreams/wflow (follow the wflow installation
instructions). In due time the wflow modelbuilder will also be available
as an executable and you will only have to download the wflow executables.

Set up a wflow model
--------------------

Once you have downloaded and installed wflow, the modelbuilder is
available in this location:

::

	<your_wflow_folder>/Scripts/wtools_py/modelbuilder.py

You can run the modelbuilder script from the command line or in a batch
file. The script uses the settings.json file for the location of the
model. To run the modelbuilder script with python, use the following
command:

::

	python modelbuilder.py --geojson-path settings.json

The settings.json file must be present. The modelbuilder comes with a
default settings.json file. Its contents look like this:

::

	{
	  "type": "FeatureCollection",
	  "features": [
		{
		  "type": "Feature",
		  "properties": {},
		  "geometry": {
			"type": "Point",
			"coordinates": [
			  7.466239929199218,
			  50.31565429419649
			]
		  }
		}
	  ]
	}

In this default settings.json a set of point coordinates is specified;
the coordinates in the example are for the Moselle catchment in Germany,
a tributary of the Rhine river. To make your own settings.json file, you
can use the website `geojson.io <http://geojson.io>`__. This is a
website that can be used to easily select a location and create your own
settings.json file. Alternatively, you can change the coordinates in the default settings.json file. This settings.json must be
present for the modelbuilder to know the location of the model.

Furthermore you can specify the following options:

::

	--geojson-path		Path to a GeoJSON file with the geometry that needs to be path of the model

	--cellsize		Desired model cell size in decimal degrees (default=0.01)

	--name			Name of the wflow case (default=wflow_<modelconcept>_case)

	--model			Name of the wflow model concept (options: sbm, hbv, w3ra)(default=sbm)

	--timestep		Model time step for hbv (options: hourly, daily) (default=daily)

	--case-template		Name of the template wflow case (default=wflow_<modelconcept>_template)

	--case-path		Path where both the template and created case reside(default is the current directory)

	--fews/--no-fews	Flag indicating whether the wflow case is part of a Delft-FEWS setup (default=no-fews)

	--fews-config-path	Path to the Delft-FEWS config directory (to save default FEWS states) (default=Config)

	--dem-path		Path to a local DEM if available

	--river-path		Path to a local river shapefile if available
	
	--region-filter		Tell hydro-engine which model area to pick, by default this is everything upstream of the provided geometry, but it is also possible to get only the current catchment (catchments-intersection), or just exactly the provided geometry (region), like your own catchment polygon (options: catchments-upstream, catchments-intersection, region)(default=catchments-upstream)

Example:

::

	python modelbuilder.py --geojson-path settings.json --name wflow_moselle --cellsize 0.01

Run this command from the command line or in a batch file, and you will
have your model.

The generated model structure looks like this:

::

    data\
    inmaps\
    instate\
    intbl\
    mask\
    run_default\
    staticmaps\
    wflow_sbm.ini

To run the wlfow model, you need the staticmaps and intbl directories
and the wflow_sbm.ini file. Also the inmaps and the instate directories
are needed to run the model, but these are not filled yet. By default,
results of your model run are stored in the run_default directory, and
this directory including all its subfolders is required if you run the
model within FEWS.

In the mask folder you will find the mask that is used to clip the
model, and the grid definition in FEWS format (in grid.xml), which you
can copy-paste into the Grids.xml file in your FEWS configuration. In
the data folder you will find the data that was used to generate the
model, after clipping it from the global data: geojson files for the
catchments and rivers, and raster files for the DEM and the parameter
maps.

The wflow_sbm.ini file is the file with configuration settings that is
needed to run the wflow-sbm model. This is an example file – please
change the settings in the ini file according to your specific model
setup (see :ref:`ini-file`).

Model data
----------

Where does the data come from? This default setup of the wflow
modelbuilder is fully based on global data sets. Below you find the
specifications of the global data sets used.

Catchment delineation
~~~~~~~~~~~~~~~~~~~~~

The clipping of the global maps is done based on the model area. The
model area is based on the HydroBASINS subcatchments, level 9
(http://hydrosheds.org/page/hydrobasins). The modelbuilder determines
within which HydroBASINS subcatchment the coordinates are located that
you specified in the settings.json file, and queries all upstream
catchments as a single or multiple polygons. These subcatchments
together define the area of your model. The data sets described below
are all clipped based on this area.

Rivers
~~~~~~

For the river network, the HydroSHEDS drainage network is queried as
polylines (http://hydrosheds.org/).

Optionally, a local or improved river vector file (shapefile, geojson,
etc.) can be provided to the modelbuilder with the option ``--river-path``.
If a local river vector file is specified, this will be used instead of
the default global river file.

DEM
~~~

For the elevation data the digital elevation model (DEM) used is SRTM
v4, 30m (https://www2.jpl.nasa.gov/srtm/)

Optionally, a local or improved Digital Elevation Model (DEM) can be
provided to the modelbuilder with the option ``--dem-path``. If a local DEM
is specified, this will be used instead of the default global DEM.

Land use
~~~~~~~~

For land use the 0.5 km MODIS-based Global Land Cover Climatology map by
the USGS Land Cover Institute (LCI) is used
(https://landcover.usgs.gov/global_climatology.php). This land cover
dataset consists of 17 different classes for land cover types. The
legend for this land cover map is also provided in the template case
(and copied to your wflow model) in data/parameters/lulegend.txt

LAI
~~~

LAI (Leaf Area Index) maps for the wflow-sbm model are stored in the
staticmaps/clim directory. These are twelve maps with monthly average
LAI, based on combined AVHRR and MODIS data, derived from Liu et al. 2012 [Liu2012]_, calculated as averages over 1981-2011.

Soil type
~~~~~~~~~

A soil map indicating major soil texture types is also downloaded with
the modelbuilder (wflow_soil.map), which is derived from the Harmonized
World Soil Database (HWSD) (FAO et al. 2009 [FAO2009]_). The legend for
this soil dataset is also provided in the template case in
data/parameters/wflow_soil.csv. In the current setup with global data,
this soil map is not used, since all soil-based parameters are specified
as rasters. It can however be useful if you want to differentiate
parameters in the intbl directory based on soil type, or if you want add
more parameters as .tbl files.

Model parameters
~~~~~~~~~~~~~~~~

Parameters linked to LAI:

-  Specific leaf storage: determined from Liu 1998 [Liu1998]_
-  Storage on the woody part of the vegetation (branch and trunk
   storage): determined from Liu 1998 [Liu1998]_
-  Extinction coefficient: Van Dijk & Bruijnzeel 2001 [VanDijk2001]_

Parameters linked to soil and land use:

-  Parameters provided as maps in the staticmaps directory: based on Dai et al. 2013 [Dai2013]_ and Shangguan et al. 2014 [Shangguan2014]_ 
-  Other parameters provided as intbl files: the parameters that are not
   specified as rasters, are given in the intbl directory as .tbl files,
   which can be linked to either land use, soil type or subcatchment
   (see :ref:`Input-Parameters`). For these parameters
   a default value or values have been established.

It is important to note that with the modelbuilder setup you can easily
generate a functioning model, including the model structure and all the
rasters and other files you need, resampled to your model resolution.
However, this results by no means in a calibrated model. The parameter
maps and tables are a best first estimate based on global datasets, but
most likely need tweaking for application in a regional- or local-scale
model.

Current limitations
-------------------

At the moment it is only possible to set up a model with the
modelbuilder in the WGS84 coordinate system (EPSG:4326).

References
----------

.. [Dai2013] Dai, Y., W. Shangguan, Q. Duan, B. Liu, S. Fu, G. Niu, 2013. Development of a China Dataset of Soil Hydraulic Parameters Using Pedotransfer Functions for Land Surface Modeling. Journal of Hydrometeorology, 14:869-887.

.. [VanDijk2001] Dijk, A.I.J.M. van and L.A. Bruijnzeel (2001), Modelling rainfall interception by vegetation of variable density using an adapted analytical model. Part 1. Model description. Journal of Hydrology 247, 230-238.

.. [FAO2009] FAO/IIASA/ISRIC/ISS-CAS/JRC, 2009. Harmonized World Soil Database (version 1.1). FAO, Rome, Italy and IIASA, Laxenburg, Austria.

.. [Liu1998] Liu, S. (1998), Estimation of rainfall storage capacity in the canopies of cypress wetlands and slash pine uplands in North-Central Florida. Journal of Hydrology 207, 32-41.

.. [Liu2012] Liu, Y., R. Liu, and J. M. Chen (2012), Retrospective retrieval of long-term consistent global leaf area index (1981–2011) from combined AVHRR and MODIS data. J. Geophys. Res., 117, G04003, doi:10.1029/2012JG002084.

.. [Shangguan2014] Shangguan, W., Dai, Y., Duan, Q., Liu, B. and Yuan, H., 2014. A Global Soil Data Set for Earth System Modeling. Journal of Advances in Modeling Earth Systems, 6: 249-263.
