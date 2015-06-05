Using the .ini file
===================


A number of settings of the framework can be set in the ini file for each model.
The settings are explained din the section below.


Settings in the run section
===========================

Information for the current run can be given in the run section. Here the start
and end-time of the run as well as the timestep can be given. Alternatively a link
to a Delft-FEWS runinfo.xml file can be given. An example is shown below.

::

    [run]
    # either a runinfo file or a start and end-time are required
    #runinfofile=runinfo.xml
    starttime= 1995-01-31 00:00:00
    endtime= 1995-02-28 00:00:00
    # required, base timestep of the model
    timestepsecs = 86400



Settings in the framework section
=================================

The in/output file formats can be specified in the framework section. At present only pcraster mapstacks and netcdf are
available fot input. See the supplied pcr2netcdf.py script for information on the layout of the netcdf files.
If netcdf files are used the name of the mapstack is used as the standardname in the netcdf file.

::

    [framework]
    # outputformat for the *dynamic* mapstacks (not the states and summary maps)
    # 1: pcraster
    # 2: numpy
    # 3: matlab
    outputformat=1

    # netcdfoutput requires also outputformat = 1 (default) and additionally the name of the file
    netcdfoutput = outmaps.nc
    netcdfwritebuffer=100
    netcdfinput= inmaps.nc

    # Provide a lot of debug info
    # debug=1




Settings in the API section
---------------------------

In the ini file example below several variables are configured to be available via the
API.  For most settings this only defines
what the API will expose to the outside world. However, if you specify 0 (input)
as a role for one of the forcing variables the ``wf_readmap`` function will no longer
read maps from disk for that variable but will return the contents of that
variable via the API.

The API section specifies variables that are exposed via the api. Use the following
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

example:
::

    [API]
    FreeWater=2,4
    SoilMoisture=2,4
    UpperZoneStorage=2,4
    LowerZoneStorage=2,4
    InterceptionStorage=2,4
    SurfaceRunoff=2,1
    WaterLevel=2,2
    DrySnow=2,4
    Percolation=1,0
    ForecQ_qmec=0,1
    PERC=3,5
    FC=3,4
    # Below are the forcing variables. By putting these here you MUST
    # supply them via the API, if not these will default to 0.0
    #P=0,0
    PET=0,0
    TEMP=0,3



Settings in the modelparameters section
---------------------------------------

Most of the time this section is not needed as this will mostly be configured
in the python code by the model developer. However, in some case this section can be used
alter the model for example force the model to read RootingDepth from an external data source.

The format of entries in this section is as follows::

    name=stack,type,default,verbose,[lookupmap_1],[lookupmap_2],lookupmap_n]


+ name - Name of the parameter (internal variable, without the self.)
+ stack - Name of the mapstack (representation on disk or in mem) relative to case
+ type - Type of parameter (default = static)
+ default - Default value if map/tbl is not present
+ verbose - If set to 1 (True) the maps a log message will be generated is a default values is used
             instead of a map
+ lookupmap - map(s) [0 to n] to be used in the lookuptable in the case the type is tbl


Possible parameter types (the second option)are:

+ staticmap: Read at startup from map
+ statictbl: [deprecated] Read at startup from tbl, fallback to map (need Landuse, Soil and TopoId (subcatch) maps)!
+ tbl: Read at startup from tbl and ar runtime fallback to map, lookuptable maps define here.
+ timeseries: read map for each timestep
+ monthlyclim: read a map corresponding to the current month (12 maps in total)
+ dailyclim: read a map corresponding to the current day of the year (366 maps in total)
+ hourlyclim: [not implemented yet] read a map corresponding to the current hour of the day (24 maps in total)


Example::

    [modelparameters]
    RootingDepth=monthlyclim/ROOTS,monthlyclim,75,0
    # Force the model to read monthly climatology of P
    Precipitation=inmaps/P,monthlyclim,0.0,1


Example::

    [modelparameters]
    RootingDepth=monthlyclim/ROOT,monthyclim,100,1
    Sl=inmaps/clim/LCtoSpecificLeafStorage.tbl,tbl,0.5,1,inmaps/clim/LC.map
    Kext=inmaps/clim/LCtoSpecificLeafStorage.tbl,tbl,0.5,1,inmaps/clim/LC.map
    Swood=inmaps/clim/LCtoBranchTrunkStorage.tbl,tbl,0.5,1,inmaps/clim/LC.map
    LAI=inmaps/clim/LAI,monthlyclim,1.0,1


Settings in the summary_* sections
----------------------------------

By adding variable in one or several of these sectiosn the framework will save
these variables to disk (using the value at the end, sum, min, max or avg) at the end of a run.

the available sections are:

+ summary - Saves the actual value of the variable
+ summary_avg - Saves the average value over all timesteps of the variable
+ summary_sum - Saves the sum over all timesteps of the variable
+ summary_min - Saves the minimum value over all timesteps of the variable
+ summary_max - Saves the maximum value over all timesteps of the variable

All maps are saved in the outsum directory of the current runid.

Example::

    [summary]
    self.MaxLeakage=MaxLeakage.map
    # Save and average these per LU type

    [summary_sum]
    self.Precipitation=Sumprecip.map

    [summary_max]
    self.Precipitation=maxprecip.map

    [summary_min]
    self.Temperature=mintemp.map

    [summary_avg]
    self.Precipitation=avgprecip.map


wf_DynamicFramework Module
==========================

.. automodule:: wf_DynamicFramework
    :members: wf_DynamicFramework
    :undoc-members:
