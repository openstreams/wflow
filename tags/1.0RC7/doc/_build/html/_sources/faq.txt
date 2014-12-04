Questions and answers
=====================

Questions
---------

[1]_ The discharge in the timeseries output gives weird numbers (1E31) what is going wrong?

[2]_ How do a setup a wflow model?

[3]_ Why do I have missing values in my model output?

[4]_ wflow stops and complains about types not matching

[5]_ wflow complains about missing initial state maps

[6]_ in some areas the mass balance error seems large

Answers
-------

.. [1] *The discharge in the timeseries output gives weird numbers (1E31) what is going wrong?*
    The 1E31 values indicate missing values. This probably means that at least one
    of the cells in the part upstreasm of the discharge points has a missing value. 
    Missing values are routed downstream so any missing values upstreams of a discharge
    will cause the discharge to eventually become a missing value. To resolve this check the following:

    - Check if the .tbl files are correct (do they cover all values in the landuse soil and subcatchment maps)
    - check for missing values in the input maps
    - check of model parameters are within the working range: e.g. you have set a parameter (e.g. the canopy gap fraction in the interception model > 1) to an unrealistic value
    - check all maps in the runId/outsum directory so see at which stage the missing values starts
    - the soil/landuse/catchment maps does not cover the whole domain

    .. note::
		note that missing values in upstreams cells are routed down and will eventually make
		all downstreams values missing. Check the maps in the runid/outsum directory to see if the tbl files are correct


.. [2] *How do a setup a wflow model?*
    First read the section on :ref:`Setting-up a-new-model`. Next check one of the supplied example models

.. [3] *Why do I have missing values in my model output?*
    See question [1]_


.. [4] *wflow stops and complains about types not matching*
	The underlying pcraster framework is very picky about data types. As such the maps must all be of the
	expected type. e.g. your landuse map MUST be nominal. See the pcraster documentation at pcraster.eu
	for more information

     .. note::
          If you create maps with qgis (or gdal) specify the right output type (e.g. Float32 for scalar maps)

.. [5] *wflow complains about missing initial state maps*
    run the model with the -I option first and copy the resulting files in runid/outstate back to the instate directory

.. [6] *in some areas the mass balance error seems large*
   The simple explicit solution of most models can cuase this, especially when parameter values are outside
   the nomally used range and with large timsteps. For example, setting the soil depth to zero will usually cuase large errors. The solution is usually to check the parameters throughout the model.