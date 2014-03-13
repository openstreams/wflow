Questions and answers
=====================

Questions
---------

[1]_ The discharge in the timeseries output gives weird numbers (!E31) what is going wrong?

[2]_ How do a setup a wflow model?



Answers
-------

.. [1] The 1E31 values indicate missing values. This probably means that at least one
    of the cells in the part upstreasm of the discharge points has a missing value. 
    Missing values are routed downstreams so any missing values upstreams of a discharge
    will cause the discharge to evetually become a missing value. To resolve this check the following:

    - Check if the .tbl files are correct ( do that cover all va;ues in the landuse soil and subcatchment maps)
    - check for missing values in the input maps
    - check of model parameters are within the working range
    - check all maps in the runId/outsum directory so see at which stage the missing values starts

.. [2] First read the section on :ref:`Setting-up a-new-model`. Next check one of the supplied example models