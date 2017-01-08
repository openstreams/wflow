The wflow_bmi_combined interface
================================


Introduction
------------
In order to simplify conversion of an existing model to a reusable,
plug-and-play model component, CSDMS has developed a simple interface
called the Basic Model Interface or BMI that model developers are asked to implement.
Recall that in this context an interface is a named set of functions with prescribed
function names, argument types and return types. The BMI functions make the model
self-describing and fully controllable by a modeling framework.

See also: http://csdms.colorado.edu/wiki/BMI_Description

The wflow\_bmi\_combined module implements a class that connects 2 or more python bmi modules
and exports those to the outside as a single bmi model.

+ A @ character is used to separate the module from the variable. For example, the variable
  Flow in module wflow_hbv becomes wflow_hbv@Flow in the combined model (it was Flow in the single interface)
+ The individual models can run in a separate case dir or can be combined into one directory (and share maps
  that are identical)


The bmi2runner.py script can be used to run a set of combined models, it is documented seperately.

wflow_bmi_combined module documentation
---------------------------------------

.. automodule:: wflow_bmi_combined
    :members:
    :undoc-members:
    :show-inheritance:
