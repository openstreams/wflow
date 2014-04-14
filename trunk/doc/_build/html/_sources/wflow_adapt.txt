.. _wflow_adapt:



wflow_adapt Module
==================


Introduction
------------

wflow_adapt is an adapter that links wflow to Delft-FEWS 
(http://publicwiki.deltares.nl/display/FEWSDOC/Home). it is typically run from the
Delft-FEWS general adapter.  


Linking wflow models to Delft-FEWS
----------------------------------


To run the model from Delft-FEWS the following actions need to be
performed:

-  wflow_[sbm|hbv].py needs to be run with the -F option where the argument refers
   to a Delft-FEWS runinfo.xml file

-  you need to specify fewsrun=1 in the model section of the .ini file

-  The postadapter (wflow\_adapt.py) needs to be run after the wflow run

Because DELFT-FEWS exports the mapstacks beginning at 0 and pcraster
expects them to start at 1 you will need to add a delay of one
timesstep to mapstack timeseries exported to wflow. This will mean
the first timestep (.000) is empty but that one will be ignored by
wflow anyway.

The adapter also tries to convert the log messages of the model into Delft-FEWS diagnostics
XML format.

- Casename\runid\wflow.log is converted to wflow_diag.xml

- Also the adapter log files is converted to wflow_adapt_diag.xml


Command line arguments:



An example of executing wflow from the Delft-FEWS general adapter
is shown below:

   
::

    <executeActivities>
    <executeActivity> 
    <description>Run wflow</description> 
    <command><executable>bin-wflow\wflow_sbm.exe</executable></command>
        <arguments>
        <argument>-C</argument>
        <argument>rhine</argument>
        <argument>-F</argument>
        <argument>rhine/inmaps/runinfo.xml</argument>
        <argument>-f</argument>
    </arguments> 
    <timeOut>7200000</timeOut> 
    </executeActivity> 
    <executeActivity> 
    <description>Run wflow post</description> 
    <command> <executable>bin-wflow\wflow_adapt.exe</executable> </command> <arguments> 
        <argument>-M</argument> 
        <argument>Post</argument> 
        <argument>-s</argument> 
        <argument>rhine/instate/state.xml</argument> 
        <argument>-o</argument> 
        <argument>rhine/instate/outstate.xml</argument> 
        <argument>-r</argument> 
        <argument>rhine/inmaps/runinfo.xml</argument> 
        <argument>-w</argument> 
        <argument>./</argument> 
        <argument>-C</argument> 
        <argument>rhine</argument> 
        <argument>-I</argument>
        <argument>wflow_sbm.ini</argument> 
        <argument>-T</argument>
        <argument>86400</argument>
    </arguments> 
    <timeOut>1200000</timeOut>
    <overrulingDiagnosticFile>wflow_diag.xml</overrulingDiagnosticFile>
    </executeActivity> 
    </executeActivities> 



The wflow_adapt module can also be used by other programs to convert .tss files to
pi-xml vv. Below the API documentation of the module is given.




Module function documentation
-----------------------------

.. automodule:: wflow_adapt
    :members:
