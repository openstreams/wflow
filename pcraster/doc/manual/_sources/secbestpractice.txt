.. _secbestpractice:

*******************
best practice guide
*******************

.. index:: if then (else); pittfalls

.. _ifthenelsepitfalls:

if then (else) pitfalls
=======================


The use of boolean true/false data in models is a difficult
subject for many users. We encounter a lot of very complicated expressions in user
models that are unnecessary and error prune. 

Very often fragments are used like:

::

   nearGroundWater = if ( GroundWaterLevel < 20 then 1 else 0);

or

::

   nearGroundWater = boolean(if ( GroundWaterLevel < 20 then 1 else 0));

This is identical to the much shorter and faster:

::

    nearGroundWater = GroundWaterLevel < 20;

The first expression of the :ref:`if-then-else <ifthenelse>` is already a boolean expression
consisting of 0 (false) and 1 (true) and possible missing values, so there is
no need to write the values of 1 and 0 explicit. This happens often when using
any of the :ref:`groupbool` or :ref:`groupcomp`.

More dangerous is the common believe that a statement like:

::

  nearGroundWater = if ( GroundWaterLevel < 20 then 1); # WRONG!

will result in a true/false map, as in the statements given above. This is
wrong: :ref:`if without the else <ifthen>` will create missing values not 0 (false).

People also tend to think in terms of "if this then that happens, if not then
if this the case then etc. etc.". This results in statements like:

::

  A      = aDepth < 20;
  B      = bDepth < 20;
  result = boolean(if( A eq 1, 1, if(B eq 1, 1, 0)));

The same is better expressed as:

::

  result = aDepth < 20 or bDepth < 20;

To summarize in a few guidelines:

- if-then and if-then-else are two totally different operations:

   - if-then is only to create missing values in maps.
   - if-then-else is to make a selection between 2 values.

- in pcrcalc: Use the syntax if( then else ) and if( then ) instead of if( , , ) and if( , ) to see better which type of operation you use.

- Avoid patterns like if( expression then 0 else 1) or if( expression then
  1 else 0) in your model. In most cases these patterns can be rewritten.
  If these patterns are nested like in if (.. (if (... try rephrasing in
  terms of :ref:`and`, :ref:`or`, :ref:`xor` and ref:`not`.
