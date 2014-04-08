Questions and answers
=====================

Questions
---------

[1]_ The discharge in the timeseries output gives weird numbers (1E31) what is going wrong?

[2]_ How do a setup a wflow model?

[3]_ Why do I have missing values in my model output?

Answers
-------

.. [1] *The discharge in the timeseries output gives weird numbers (1E31) what is going wrong?*
    The 1E31 values indicate missing values. This probably means that at least one
    of the cells in the part upstreasm of the discharge points has a missing value. 
    Missing values are routed downstream so any missing values upstreams of a discharge
    will cause the discharge to eventually become a missing value. To resolve this check the following:

    - Check if the .tbl files are correct ( do that cover all values in the landuse soil and subcatchment maps)
    - check for missing values in the input maps
    - check of model parameters are within the working range
    - check all maps in the runId/outsum directory so see at which stage the missing values starts

.. [2] *How do a setup a wflow model?*
    First read the section on :ref:`Setting-up a-new-model`. Next check one of the supplied example models

.. [3] *Why do I have missing values in my model output?*
    See question [1]_
