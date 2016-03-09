bmi2runner
==========

Introduction
------------

bmi2runner.py is a simple script that runs two or more wflow modules connection via the BMI. A configfile
is used to control the exchange of data between the models.

The config file contains a list of models configured as the name of the wflow module and the ini file that is used
by the model. Furthermore, in the exchanges section the data flows from model to model are configures.
::

    [models]
    wflow_sbm=wflow_sbm@wflow_sbm_comb.ini
    wflow_routing=wflow_routing@wflow_routing_comb.ini

    [exchanges]
    # From_model/var -> To_model/var
    wflow_sbm@InwaterMM=wflow_routing@IW

To setup a combined model the you shoudl first configure and setup the individual models. They can be setup in separate
case directories or they can be setup in one case directory. Each model should have it's onw config/ini files.

+ the models are executed in the order they are listed in the models section
+ the variables are get/set in the order they appear in the exchanges section
+ the script runs explicitly, no iteration is performed

bmi2runner script
=================

.. automodule:: bmi2runner
    :members: bmi2runner
    :undoc-members:

