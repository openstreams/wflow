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







wf_DynamicFramework Module
==========================

.. automodule:: wf_DynamicFramework
    :members: wf_DynamicFramework
    :undoc-members:
