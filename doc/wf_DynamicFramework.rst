Using the .ini file
===================


A number of settings of the framework can be set in the ini file for each model.
The settings are explained in the section below.


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


If this section is not present and a runinfo.xml is also not used you will need
to specify the number of timesteps using the -T option on the command line (for most models).


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


Use role 0 for input maps to the model (forcing data that are normally read from disk) only, role 1
for outputs, role 2 for state variables and role 3 for model parameters. The units may be choose freely and
be strings also.

example:
::

    [API]
    FreeWater=2,4
    SoilMoisture=2,4
    UpperZoneStorage=2,4
    LowerZoneStorage=2,4
    InterceptionStorage=2,4
    SurfaceRunoff=2,m^3/sec
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
Not all models support this. You can check if the model you uses support this by looking for the
wf_updateparameters() function in de model code.

The format of entries in this section is as follows:

::

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
+ tbl: Read at startup from tbl lookuptable maps define here.
+ tblts: Lookup tables for each timestep, including initial section
+ tblsparse: Lookup tables for each timestep, including initial section. Fills in missing
  timestep using previous timestep
+ timeseries: read map for each timestep
+ monthlyclim: read a map corresponding to the current month (12 maps in total)
+ dailyclim: read a map corresponding to the current day of the year (366 maps in total)
+ hourlyclim: [not implemented yet] read a map corresponding to the current hour of the day (24 maps in total)
+ tss: read a tss file and link to lookupmap (only one allowed) a map using timeinputscalar

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



Settings in the variable_change_timestep/once section
-----------------------------------------------------
In the two sections "variable_change_timestep" and "variable_change_once" you can set
operations on parameters and variable that are executed at the start of each timestep or once in the initialisation
of the model respectively. What you specify here should be valid python code and include variable that exists
in the model you are using. This only works if the actual model you are using includes the wf_multparameters() function.
At the moment wflow\_hbv, wflow\_sbm, wflow\_w3ra and wflow\_routing include this.
See below for a configuration example. Some models may also support this via the -P and -p command-line options.

::

    [variable_change_timestep]
    self.Precipitation = self.Precipitation * 1.2
    # Mutiplies the precipitation input times 1.2 every timestep

    [variable_change_once]
    self.PathFrac = self.PathFrac * 1.1
    # Increases the paved area fraction of the model by 10 %



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



Settings in the outputtss/outputcsv sections
--------------------------------------------
[outputcsv_0-n]
[outputtss_0-n]

Number of sections to define output timeseries in csv format. Each section
should at lears contain one samplemap item and one or more variables to save.
The samplemap is the map that determines how the timeseries are averaged/sampled. The function
key specifies how the data is sample: average(default), minimum, maximum, total, majority.

All other items are variable=filename pairs. The filename is given relative
to the case directory.

Example:

::

    [outputcsv_0]
    samplemap=staticmaps/wflow_subcatch.map
    self.SurfaceRunoffMM=Qsubcatch_avg.csv
    function=average
    # average is the default

    [outputcsv_1]
    samplemap=staticmaps/wflow_gauges.map
    self.SurfaceRunoffMM=Qgauge.csv
    self.WaterLevel=Hgauge.csv

    [outputtss_0]
    samplemap=staticmaps/wflow_landuse.map
    self.SurfaceRunoffMM=Qlu.tss
    function=total



In the above example the discharge of this model (self.SurfaceRunoffMM) is
saved as an average per subcatchment, a sample at the gauge locations and as
an average per landuse.



[run] section: The use of date and time
---------------------------------------

Available options in the [run] section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The run section can contain information about the model timesteps, the date/time range,
how to initialize the model and how interpret the forcing data.

::

    [run]
    # either a runinfo file or a start and end-time are required
    #runinfo=runinfo.xml
    starttime=2016-01-01 00:00:00
    endtime=2016-01-04 00:00:00
    # Base timestep of the model in seconds, default = 86400
    timestepsecs = 86400
    #start model with cold state
    reinit=1
    # Default behaviour: steps
    runlengthdetermination=steps

Data and time and timesteps
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The original pcraster framework has no notion of date and time, only timesteps that are used to propagate a
model forward. However, to be able to support the BMI and netcdf files dat and time functionality
has been inserted into the framework.

As most of the forcing in hydrological modelling is an accumulation the dat/time of an input time series is
assumed to represent the total (or average) of that variable one timestep length back in time from the timestamp.

For example, if the forcing data has the following four timestamps the model will run four timesteps. The first timesteps
will be to propagate the state from T0 to T1 (in that case the date/time of the state files is assumed to be T0).
As such the state going into the model should be valid for the T0, see graph below:


.. graphviz::

    digraph G {
        node [shape=box]

        subgraph cluster_1 {
	    	node [style=filled];
		    " "  -> T1 -> T2 -> T3 -> T4;
		    label = "Model steps";
		    color=blue
	    }

        subgraph cluster_2 {
	    	node [style=filled];
		    "  " ->"1 Jan 2016" -> "2 Jan 2016" -> "3 Jan 2016" -> "4 jan 2016";
		    label = "Forcing data";
		    color=blue
	    }


        "Model state: T0 31 Dec 2015";


    }

Here the first column shows the model steps and the second column the timestamp of the input/output data for that timestep
(a - means there is nothing for that step). The first empty-row is regarded as the timestamp of the initial conditions. The
first forcing data is used to propagate the model from T0 (the state, 31 Dec 2015) to T1 (1 Jan 2016) whcih will also
be the first output the model writes..

The above corresponds to the following date/time in rhe [run] section:

::

    [run]
    # either a runinfo file or a start and end-time are required
    #runinfo=runinfo.xml
    starttime=2016-01-01 00:00:00
    endtime=2016-01-04 00:00:00
    # required, base timestep of the model
    timestepsecs = 86400
    #start model with cold state
    reinit=1
    # Default behaviour: steps
    runlengthdetermination=steps



The above shows the default behaviour of the framework. For each data point in the input forcing a model steps is performed.
The 'runlengthdetermination' variable in the run section can also be set to 'intervals'. In that case the the number
 of steps is determined from the number of intervals in the forcing data. Hence, the following run will be performed:


.. graphviz::

    digraph G {
        node [shape=box]

        subgraph cluster_1 {
	    	node [style=filled];
		    " "  -> T1 -> T2 -> T3;
		    label = "Model steps";
		    color=blue
	    }

        subgraph cluster_2 {
	    	node [style=filled];
		    "1 Jan 2016" -> "2 Jan 2016" -> "3 Jan 2016" -> "4 jan 2016";
		    label = "Forcing data";
		    color=blue
	    }


        "Model state: T0 01 Jan 2016";


    }


::

    [run]
    # either a runinfo file or a start and end-time are required
    #runinfo=runinfo.xml
    starttime=2016-01-01 00:00:00
    endtime=2016-01-04 00:00:00
    # required, base timestep of the model
    timestepsecs = 86400
    #start model with cold state
    reinit=1
    # Default behaviour: steps
    runlengthdetermination=intervals



In this case the forcing data has the same input as in the previous case (4 timestamps) but in this case the first
timestamp is regarded as the time of the initial conditions. As such, one timestep less (now three in total) is performed and the
forcing data from 01 Jan 2016 is NOT used as it is regarded as the initial state. The first output point
will be  02 Jan 2016 generated from the forcing marked as 2 Jan 2016.


The same applies when the start and end time of the model run are supplied via the bmi interface.




wf_DynamicFramework Module
==========================

.. automodule:: wf_DynamicFramework
    :members: wf_DynamicFramework
    :undoc-members:
