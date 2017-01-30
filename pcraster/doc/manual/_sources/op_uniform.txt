

.. index::
   single: uniform
.. _uniform:

*******
uniform
*******
.. topic:: uniform

   Boolean TRUE cell gets value from an uniform distribution

::

  Result = uniform(expression)

expression
   spatial, non spatial
   boolean

Result
   spatial
   scalar

Operation
=========


A random generator is used to generate the Result: for each cell that has a value 1 (TRUE) on expression, a value is taken from a uniform distribution between 0 and 1 and assigned to the cell on Result. Cells that have a value 0 (FALSE) on expression area assigned a missing value.   

Notes
=====


Group
=====
This operation belongs to the group of  Random number generators; Cells 

See Also
========
:ref:`groupareafield`, :ref:`groupmapfield`

Examples
========
