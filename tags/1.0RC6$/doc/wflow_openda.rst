Linking wflow to OpenDA
=======================

.. warning:: 

	This is experimental and incomplete



There are lot of folders that need to be set in the PATH and PYTHONPATH environment variables for this to work. Also you need the latest wflow python scripts, since Jaap fixed a bug a few days ago. These are also in Albrecht's postbox::

\\hafilerk\homes\Deltabox\Postbox\Weerts, Albrecht\vanArno\testWflow\wflow_bin\


You need to set the following folders in the following environment variables:

PATH:
::

 1. folder with jep.dll file (e.g. f:\testWflow\openda_bin\win32_ifort\)
 2. folder with PCRaster dll files (e.g. C:\PCRaster\apps).
 3. folder with CPython installation (e.g. C:\Python25)

PYTHONPATH:
::

 1. folder with WFLOW python scripts (e.g. f:\testWflow\wflow_bin\)
 2. folder with PCRaster python scripts (e.g. C:\PCRaster\Python).


 Furthermore here is the description of the arguments that are configured in the bbStochModelConfig.xml file:
 6 arguments expected:

::

 1: piRunFile: pi run info file path (relative to working dir)."
 2: name of the python module that contains the model to use (e.g. wflow_hbv)."
 3: caseDirectory: case directory path (relative to working dir)."
 4: templateRunId: name of the runId within the current case (relative to the caseDirectory)."
 5: configFileName: name of the model configuration file (relative to the caseDirectory)."
 6: cloneMapFileName: name of the map file that describes the catchment (relative to the staticmaps folder in the caseDirectory)."


This should make it run without any errors, however Albrecht and I are still working to get the OpenDA code working correctly, so we might run into some more problems to solve.

.. note::
  See also the chapter on configuring the models. There are a number of settings in the ini file (the API section) that are 
  important. See [run wflow\_hbv via API]_