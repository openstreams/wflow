run wflow\_hbv via API 
======================

Settings in the ini file
------------------------

In the ini file example below several variables are configured to be available via the
API.  For most settings this only defines
what the API will expose to the outside world. However, if you specify 0 (input) 
as a role for one of the forcing variables the ``wf_readmap`` function will no longer
read maps from disk for that variable but will return the contents of that 
variable via the API.


Settings in the API section
---------------------------

The API section specifies variables that are exposed via the ere. Use the following 
convention:

::

    variable_name_in_model=variable_role,variable_unit


::


    role: 0 = input (to the model)
          1 = is output (from the model)
          2 = input/output (state information)
          3 = model parameter
    unit: 0 = mm/timestep
          1 = m^3/sec
          2 = m
          3 = degree Celcius
          4 = mm
          5 = -


Use role 0 for input maps to the model (those that are normally read from disk), role 1
for outputs, role 2 for state variables and role3 for model parameters.


.. include:: ../examples/wflow_rhine_hbv/wflow_hbv_mem.ini
	:literal:

.. automodule:: testrunner_wflowhbv

.. literalinclude:: ../wflow-py/wflow/testrunner_wflowhbv.py


