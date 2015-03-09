The wflow_bmi interface
=======================


Introduction
------------
In order to simplify conversion of an existing model to a reusable,
plug-and-play model component, CSDMS has developed a simple interface
called the Basic Model Interface or BMI that model developers are asked to implement.
Recall that in this context an interface is a named set of functions with prescribed
function names, argument types and return types. The BMI functions make the model
self-describing and fully controllable by a modeling framework.

See also: http://csdms.colorado.edu/wiki/BMI_Description

This is the first implementation of the BMI for the wflow pcraster/python models

Configuration
-------------

Mapping of long_var_name to model variables not yet implemented. The long_var_name
should be model for names for now


wflow_bmi module documentation
------------------------------

.. automodule:: wflow_bmi
    :members:
    :undoc-members:
    :show-inheritance:
