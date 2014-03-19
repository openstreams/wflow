The wflow_gr4 model
===================


.. warning::

    The documentation is incomplete

Introduction
------------
An *experimental* implementation of the gr4 model. It is based on the hourly (gr4h) version

Dependencies
------------
[PM]

Configuration
-------------


The model  needs a number of settings in the ini file. The default name for the ini file
is wflow\_gr4.ini. 

See below for an example: 

::

	[model]

	Tslice=1
	# Maximum upstream distance to update the flow in metres


	[gr4]
	dt = 1
	B = 0.9
	D = 1.25
	X4 = 32.83
	# X1,X2 and X3 are given as .tbl files or maps

	[layout]
	# if set to zero the cell-size is given in lat/long (the default)
	sizeinmetres = 1

	[outputmaps]
	# Add maps here

	# List all timeseries in tss format to save in this section. Timeseries are
	# produced as averages per subcatchment. The discharge (run) timeseries 
	# is always saved (as samples at the gauge location)s.
	[outputtss]
	self.S_X1=S_X1
	self.R_X3=R_X3
	self.Pr=Pr
	self.Q=Q 





wflow_gr4 module documentation
------------------------------

.. automodule:: wflow_gr4
    :members:
    :undoc-members:
    :show-inheritance:
