============================
The wflow Delft-FEWS adapter
============================
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

-  The runinfo.xml file should be specified in the [run]section of the ini file
-  The use of netcdf input and output should be switched on
-  The postadapter (wflow\_adapt.py) needs to be run after the wflow run


The postadapter also converts the log messages of the model into Delft-FEWS diagnostics
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
        <argument>-w</argument> 
        <argument>./</argument> 
        <argument>-C</argument> 
        <argument>rhine</argument> 
        <argument>-I</argument>
        <argument>wflow_sbm.ini</argument>
    </arguments> 
    <timeOut>1200000</timeOut>
    <overrulingDiagnosticFile>wflow_diag.xml</overrulingDiagnosticFile>
    </executeActivity> 
    </executeActivities> 



The wflow_adapt module can also be used by other programs to convert .tss files to
pi-xml vv. Below the API documentation of the module is given.


In the above example the state files belonging to the model should be configed as per  below in the General Adapter XML. In
Fews the read and write locations are as viewed from the model's point of view:

::

        <stateLocation>
            <readLocation>WaterLevel.map</readLocation>
            <writeLocation>../run_default/outstate/WaterLevel.map</writeLocation>
        </stateLocation>
        # Repeat for all state variables


Module function documentation
-----------------------------

.. automodule:: wflow_adapt
    :members:
