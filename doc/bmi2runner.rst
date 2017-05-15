bmi2runner
==========

Introduction
------------

bmi2runner.py is a simple script that runs two or more wflow modules connection via the BMI
interface (the combined version). A configfile
is used to control which models to start as well as the exchange of data between the models.

The config file contains a list of models configured with the name of the wflow modules to run
and the ini file that is used by the model. Furthermore, in the exchanges section
the data flows from model to model are configured.
::

    [models]
    # module name = path to config of module relative to the dir of this ini file
    wflow_sbm=wflow_sbm/wflow_sbm_comb.ini
    wflow_routing=wflow_routing/wflow_routing_comb.ini

    [exchanges]
    # From_model/var -> To_model/var
    wflow_sbm@InwaterMM=wflow_routing@IW

To setup a combined model  you should first configure and setup the individual models. They can
be setup in separate case directories or they can be merged in one case directory.
Each model should have it's own config/ini file. The following principles apply
when using the bmi2runner script:

+ the models are executed in the order they are listed in the models section
+ the variables are get/set in the order they appear in the exchanges section
+ the script runs explicitly, no iteration is performed

Example
-------

In the examples directory the file bmirunner.ini is present.
You can use this to run a combined wflow_sbm/wflow_routing
model. Start this up using the following command (in the examples dir):

::

    bmi2runner.py -c bmirunner.ini


A second example runs wflow_sbm next wflow_routing followed by
the wflow_floodmap module:

::

     bmi2runner.py -c bmirunner-floodmap.ini

The contents of the ini file (bmirunner-floodmap.ini) is given below:

::

    [models]
    wflow_sbm=wflow_rhine_sbm/wflow_sbm.ini
    wflow_routing=wflow_routing/wflow_routing_BMI.ini
    wflow_floodmap=wflow_routing/wflow_floodmap_BMI.ini

    [exchanges]
    # From_model.var -> To_model.var
    wflow_sbm@InwaterMM=wflow_routing@IW
    wflow_routing@WaterLevel=wflow_floodmap@H

In this case the floodmap module uses the same directory as the routing module (but a different config file).



bmi2runner script
-----------------

.. automodule:: bmi2runner
    :members: bmi2runner
    :undoc-members:

