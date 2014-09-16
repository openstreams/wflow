Calibration of the wflow\_sbm model for the Rhine catchment using EOBS data
===========================================================================

DEM and landuse data
--------------------

The digital elevation model used was the SRTM 90x90 m resolution dataset. The 
DEM was used to determine the river network and the altitude of each cell using the
following steps:

# The 90x90 m DEM was used to determine the local drainage network (D8) and pits were removed from the network
# The network determined above was used to fix the river network when resampling the DEM to a 900x900 m resolution DEM.
# [pm] resampling procedure

Corina land-use map and reclassification [PM]




Forcing data
------------

The wflow\sbm model was calibrated for the Rhine basin using daily temperature
and precipitation data from the EOBS dataset [Haylock]_ (version 5.0). 
Discharge was taken from the CHR daily dataset ([gorgen]_). Potential
evapotranspiration was determined using and adjusted version Hargraeves' method
([weiland]_):


.. math::

	E_{pot} = 0.0031  (T + 17.78)  abs(T_{max} - T_{min})^{0.5}  R
	
In which 0.0031 is a calibration factor and :math:`R` is the clear sky solar
radiation expressed in mm/day.
	
Within the wflow\_sbm model the temperature is corrected for altitude
differences within one EOBS cell. For each wflow\_sbm 1x1km cell the temperature
is corrected by using the differences between the altitude in that cell and the
average altitude within the matching EOBS cell and applying a lapse rate.

Gauge data from CHR

Calibration procedure
---------------------



Calibration was performed in two steps. First an initial gues of all parameters
was made based on land-use type. Secondly  a full matrix search of a number of parameters 
was performed and best perfroming sets have been chosen per subcatchment (based on NS and bias).

The plots below show the calibration results for the period 1985 -- 1995 for the
following stations::

	Rhein-Basel, Rheinhalle
	Kalkhoven
	Rockenau
	Kaub
	Koeln
	Lobith
	Raunheim
	Cochem
	Andernach
	Maxau
	Schermbeck
	Menden
	Hattingen

The title of the plots also shows the model efficiency (Nash
and Sutcliffe). 

Original landuse
----------------

.. plot:: plots/calibplot.py

Modified landuse
----------------

.. plot:: plots/calibplot_lunew.py



.. [Haylock] Haylock, M.R., N. Hofstra, A.M.G. Klein Tank, E.J. Klok, P.D.
Jones, M. New. 2008: A European daily high-resolution gridded dataset of surface
temperature and precipitation. J. Geophys. Res (Atmospheres), 113, D20119,
doi:10.1029/2008JD10201"


.. [gorgen] Görgen, K., Beersma, J., Brahmer, G., Buiteveld, H., Carambia, M.,
de Keizer, O., Krahe, P., Nilson, E., Lammersen, R., Perrin, C. and Volken, D.
(2010) Assessment of Climate Change Impacts on Discharge in the Rhine River
Basin: Results of the RheinBlick2050 Project, CHR report, I-23, 229 pp.,
Lelystad, ISBN 978-90-70980-35-1



.. [weiland] Sperna Weiland, F. C.,  C. Tisseuil, H. H. Dürr, M. Vrac, and L. P.
H. Van Beek, “Selecting the optimal method to calculate daily global reference
potential evaporation from CFSR reanalysis data for application in a
hydrological model study,” Hydrology and Earth System Sciences, vol. 16, pp.
983–1000, 2012.
