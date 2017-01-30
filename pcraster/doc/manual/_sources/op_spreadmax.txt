

.. index::
   single: spreadmax
.. _spreadmax:

*********
spreadmax
*********
.. topic:: spreadmax

   Total friction of the shortest accumulated friction path over a map with friction values from a source cell to cell under consideration

::

  Result = spreadmax(expression)

expression
   spatial, non spatial
   nominal, ordinal, scalar, directional, ldd

Result
   dimension of expression
   boolean

Operation
=========


Identical to spread and spreadzone but with a fourth parameter, a maximum spread distance. Areas that are not 
reached are given the value for MV for spreadmax and 0 for spreadmaxzone. 



Examples
========
