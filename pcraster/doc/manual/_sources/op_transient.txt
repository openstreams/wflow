.. _transient:

*********
transient
*********
.. index::
   single: transient
.. topic:: transient

   Simulates transient groundwater flow according to the implicit finite difference method.


::

   result = transient(elevation, recharge, transmissivity, flowCondition, storageCoefficient, timestep, tolerance)

elevation
  scalar; spatial ([L])
recharge
  scalar; spatial, non-spatial ([L  T-1])
transmissivity
  scalar; spatial, non-spatial ([L2  T-1])
flowCondition
  nominal; spatial, non-spatial
storageCoefficient
  scalar; spatial, non-spatial ([L3  L-3])
timestep
  scalar; non-spatial ([T])
tolerance
  scalar; non-spatial

result
 scalar; spatial


Operation
=========
This function calculates the groundwater head of an aquifer system according to the Gauss-Seidel iteration and Crank-Nicolson method (:ref:`Wang & Anderson <WangAnderson95>`).
The ``flowCondition`` indicates whether cells are either inactive (0), active (1) or have a constant head (2).
The ``tolerance`` specifies the maximum difference between the current elevation and the new elevation.
