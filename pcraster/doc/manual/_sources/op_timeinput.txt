

.. index::
   single: timeinput
.. _timeinput:

*********
timeinput
*********
.. topic:: timeinput

   Set of output maps per timestep with an extension that refers to the time at the timestep

::

  Result = timeinput(prefixOfMap)

prefixOfMap
   spatial, non spatial
   boolean, nominal, ordinal, scalar, directional, ldd

Result
   spatial
   type of prefixOfMap

Operation
=========


This operation is used in the iterative sections (dynamic, storage and
transport sections) of a dynamic model script only. For each timestep
timeinput assigns to Result one map of a set of sequential maps given by prefixOfMap.  

prefixOfMap refers to a set of maps that all have a filename starting with the prefix prefixOfMap. Additionally these maps contain a unique time extension Ext in their file name, directly after prefixOfMap. For each map, this time extension refers to the time at a timestep in a model run: each timestep, Result is assigned the prefixOfMapExt with the time extension Ext that corresponds with the time :emphasis:`t(i)` at the timestep :emphasis:`i` under consideration.  



For each timestep :emphasis:`i` a prefixOfMapExt must be available in the PCRaster database with a time extension Ext corresponding with the time :emphasis:`t(i)` at the timestep under consideration. These maps referred to as prefixOfMap must have the following filenames. The filenames consist of 8 characters, a dot and 3 characters. This is in accordance with the ordinary rules for filenames in DOS. Each filename starts with the name of the prefix prefixOfMap. This prefix name may be maximal 8 characters long. All the remaining character positions must be used for the time extension, which is unique and different for each map. The time extension must be a whole value; the map will be assigned to Result at timestep :emphasis:`i` with time :emphasis:`t(i)` = Ext.  



Two examples are given to illustrate the use of the operator. Imagine a
model with startime = 4, endtime 10 and a timeslice of 2. As a result, this
model consists of 4 timesteps at time 4, 6, 8, 10. During a model run, the
operation timeinput(Rain) queries for maps in the PCRaster database with filenames: Rain0000.004, Rain0000.006, Rain0000.008 and Rain0000.010. These maps are assigned to Result at the sequential timesteps at time 4, 6, 8, 10.  



In a model with starttime 990, endtime 1010 and a timeslice of 10, the
operation timeinput(Water) will query for Water000.990, Water001.000 and Water001.010.  

Notes
=====


Maps that are reported in a model with the :ref:`report <secseqinreport>` keyword are stored in the database with a filename format that corresponds with the format needed for the timeinput operation. So, one model may be used to generate a set of prefix maps; these maps can be used in another model as input maps to a timeinput operation.  



A stack of maps generated and stored in the database during a model run
by the report keyword cannot be imported during the same model run with
the timeinput operator.   

Group
=====
This operation belongs to the group of  Time operators 

Examples
========
