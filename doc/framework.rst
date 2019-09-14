===============================================
Adding a new model yourself using the framework
===============================================

Using the framework
===================

This section only gives a brief description of the framework focusing on the extensions
made for OpenStreams. A full description of the current version of the framework can be found
at http://www.pcraster.eu.

In order to build a dynamic model you will needs to define a model class and add several methods
to the class to describe the model behaviour. The easiest way to get started is to copy
and modify the ``wflow_sceleton.py`` example model. You can also use the other models
for inspiration.

In order to facilitate reusing data between models the data is stored in the following
directory tree:

.. digraph:: file_system

   //rankdir=LR;
   size="8,11";
   "Case" -> "inmaps";
   "Case" -> "instate";
   "Case" -> "intbl";
   "Case" -> "intss";
   "Case" -> "outstate";
   "Case" -> "Run";
   "Case" -> "staticmaps";
      "Run" -> " intbl";
      "Run" -> "outmaps";
      "Run" -> " outstate";
      "Run" -> "outsum";
      "Run" -> "runinfo";

Although it is possible to deviate from this layout it is highly recommended to 
adhere to this if you build your own model. Also make sure you use an ini file to 
specify model settings instead of putting those in the python code.

A basic sceleton of a model is given below:

.. automodule:: wflow_sceleton
    :members: WflowModel, main
    :private-members:
    :special-members:
    :member-order: bysource



Annotated source code for the above
-----------------------------------

.. include:: ../wflow/wflow_sceleton.py
	:literal:





