

.. index::
   single: exp
.. _exp:

***
exp
***
.. topic:: exp

   Base\ :sub:`e` exponential

::

  Result = exp(power)

power
   spatial, non spatial
   scalar

Result
   dimension of power
   scalar

Operation
=========


For each cell, raises e to the Nth power, where N is the cell value
on power. The result of this calculation is assigned to the corresponding cell on Result.  

Notes
=====


A cell a with missing value on power is assigned a missing value on Result.  

Group
=====
This operation belongs to the group of  Arithmetic operators 

Examples
========
#. 
   | • pcrcalc
   |   binding
   |    Result = Result.map;
   |    Power = Power.map;
   |   initial
   |    report Result = exp(Power);
   |   
   | • python
   |   Power = readmap("Power.map")
   |   Result = exp(Power)

   ====================================== =====================================
   Result.map                             Power.map                            
   .. image::  ../examples/exp_Result.png .. image::  ../examples/exp_Power.png
   ====================================== =====================================

   | 

