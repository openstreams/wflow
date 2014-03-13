The wflow_fit module
====================


Introduction
------------
The wflow\_fit module provides simple automated least square fitting 
for the wflow models. It uses the scipy.optimize function to perform
the fitting.

The program works mu multipling the fit parameter with a factor and optimise this
factor. To get the new optimised parameters for your model you have to 
multiply your original parameters with the optimised factor. You
can specify measured and simulated Q pairs to use and which area of the model
you wan to adjust for each Simulated/Measured pair

In order to use the fit module you must have a:

- A working wflow model
- a tss file with measured discharge
- an [fit] section in the ini file




The ini file
------------
To be able to use the fit module you must add a [fit] section to the
.ini file of the wflow model you want to fit.


::

	[fit]
          # The parameters are name parameter_0 to parameter_n
	parameter_0 = M
	parameter_1 = RootingDepth
          # Q specifies the tss file with measure discharge data
          # the path is relative to the case directory
	Q = testing.tss
          # The columns in the measured Q you want to fit to
	ColMeas = [1,5]
          # The columns in the measured Q you want to fit
	ColSim = [1,5]
          # Number of warup timesteps. This are not used in fitting
	WarmUpSteps = 1
          # The map defining the areas you want to adjust
	areamap=staticmaps/wflow_catchment.map
          # The areas you want to adjust for each Qmeas/Qsim combination
	areacode=[1,5]
	

	
	
Fitting results
---------------

Results are saved in the wflow_fit.res file in the case/runid directory. In addition,
the program saves a graph of modelled and observed data in the file fit.png and maps
of the original and fitted parameters are also saved.

If you specify the -U option the resulting maps are saved in the staticmaps directory
after each steps. As such, next steps (if you calibrate multiple subcatchments/areas)
also include the results of the previous steps. Note that this will overwrite your 
maps if you already have those!



How to fit
----------

Although wflow\_sbm has a fairly large number of parameters most should not be 
fitted automatically. The parameters that are most suited for fitting are:

- M
- FirstZoneKsatVer
- RunoffGeneratingGWPerc (if this is switched on. It is usually best to first 
  setup the model without this parameter!)
- RootingDepth

It is recommended to only fit one or two parameters at one time.

The wflow\_rhine\_sbm example can be used to test the fitting procedure. 

::

    wflow_fit.py -M wflow\_sbm -T 300 -C wflow\rhine\_sbm

Description of the python module
================================

.. automodule:: wflow_fit
	:members:
    



    
        









    


