

.. index::
   single: time
.. _time:

****
time
****
.. topic:: time

   Timestep

::

  Result = time()

Result
   non spatial
   scalar

Operation
=========


This operation is used in the iterative sections (dynamic, storage and
transport sections) of a dynamic model script only. For each timestep in a
model run, the operator assigns to Result the time :emphasis:`t` at the timestep :emphasis:`i` under consideration. The time at the first timestep (i = 1) is t(start), the time at the second timestep (i = 2) is t(start) + timeslice, at the third time step (i = 3) is t(start) + 2 x timeslice, etc. The time dimension in a model (i.a. t(start), timeslice) is defined in the timer section of a :ref:`Dynamic Modelling Script <secseqscrtime>`.   

Group
=====
This operation belongs to the group of  Time operators 

Examples
========
