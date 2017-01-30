.. _tutorial:

Tutorial
********

The tutorial gives a short introduction into the usage of the arrayed variables in Python.


Importing the module
--------------------

To use the arrayed variables in your Python script use::

  from PCRaster.Collection import *


Creating collection index types
-------------------------------

In oldcalc, index types are created as follows::

  # oldcalc
  PlantSpecies = [
        Species1,
        Species2,
        Species3];

In Python, an index type is created as::

  # Python
  PlantSpecies = Index(["Species1", "Species2"])

.. As a convention, the name argument of the index must equal the variable name.

The indices are given as a list of names.

Creating variable collections
-----------------------------

In oldcalc, arrayed variables are created as follows::

  # oldcalc
  Variable1[Species1] = ...;
  Variable1[Species2] = ...;

In Python, you create a variable collection object by using the index types::

  # Python
  PlantSpecies = Index(["Species1", "Species2"])
  Variable1 = VariableCollection([PlantSpecies], value=None)

This will initialise the indices with the default value ``None``. Afterwards, you can assign values as::

  Variable1[PlantSpecies.Species1] = ...
  Variable1[PlantSpecies.Species2] = ...

As an alternative, you can use an external parameter file to initialise the values (collection values not specified in the external file will be set to the default value ``None``) like::

  # Python
  Variable1 = VariableCollection([PlantSpecies], value=ValueFromParameterTable("Variable1", "values.tbl", dataType=Scalar))

where ``Variable1`` is the name of the collection variable referred to in the external file, ``values.tbl`` the filename of the external file and ``Scalar`` the PCRaster data type of the collection. In case of initialising collection indices with PCRaster field objects, i.e. reading maps from disk, set the ``dataType`` argument to ``None``.

Creating multidimensional collections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Multidimensional collections are created in a similar way as one-dimensional arrays::

  # Python
  PlantSpecies = Index(["Plant1", "Plant2"])
  Herbivores = Index(["Species1", "Species2", "Species3"])

  Interaction = VariableCollection([Herbivores, PlantSpecies], value=None)

This will create a variable ``Interaction`` holding a 2-dimensional collection. You can set values like::

  Interaction[Herbivores.Species1, PlantSpecies.Plant1] = scalar(2.5)
  # repeat initialisation for the remaining Herbivore-PlantSpecies combinations

Again, the variable can be initialised from a parameter file insetead::

  Interaction = VariableCollection([Herbivores, PlantSpecies], value=ValueFromParameterTable("Interaction", "interaction.tbl", dataType=Scalar))

The parameter file looks similar to::

  #interaction.tbl
  Interaction  Species1 Plant1 0.2
  Interaction  Species1 Plant2 0.35
  Interaction  Species2 Plant1 0.33
  Interaction  Species2 Plant2 0.35
  Interaction  Species3 Plant1 0.4
  Interaction  Species3 Plant2 0.5



Iterating over the collection
-----------------------------

While in oldcalc a ``foreach`` construct iterates over the array::

  # oldcalc
  foreach species in PlantSpecies {
    Variable1[species] = ...;
    Variable2[species] = ...;
  }

In Python there are two options to iterate over the collection. First, you can use the ``for in`` construct and iterate over the index types::

  for species in PlantSpecies:
    Variable1[species] = ...
    Variable2[species] = ...


Iterating over multidimensional collections can be done with nested ``for in`` statements::

  for herb in Herbivores:
    for plant in Plants:
      Interaction[herb, plant] = ...

As a second option you can iterate over the collection::

  for herb, plant in Interaction:
    Interaction[herb, plant] = ...

Using ``None`` values in a collection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not necessarily all values in a collection need to be specified. E.g. in our Herbivore-PlantSpecies interaction example, some interactions between species do not appear. If so, you can assign ``None`` to that variable::

  Interaction[Herbivores.Species3, PlantSpecies.Plant2] = None

When iterating over index types, it is now required to test if values are set::

  for herb in Herbivores:
    for plant in Plants:
      if not Interaction[herb, plant] is None:
        Interaction[herb, plant] = ...

The collection iterator only operates on values being not ``None``, therefore it is not necessary to change the code::

  for herb, plant in Interaction:
    Interaction[herb, plant] = ...

Reporting timeseries
--------------------

To report timeseries for each collection index the variable collection class needs to be initialised with the ``ValueTimeoutputTimeseries`` class::

  InteractionTss = VariableCollection([Herbivores, PlantSpecies], value=ValueTimeoutputTimeseries("Interaction", self, idMap="locations.map", noHeader=False)))

Afterwards you can sample values for each collection index::

  for herb, plant in Interaction:
    InteractionTss[herb, plant].sample(Interaction[herb, plant])

The resulting output files will be named like ``Interaction-Species3-Plant2.tss``.

Linking variable names to external names
----------------------------------------

In the case you want to use variable names in your model script different from the ones in the parameter files you can link those external names to variable names while you initialise the ``Index`` class::

  Herbs = Index(["ES=EvergreenShrub", "SPG=ShortPerennialGrass", "Cal=Calluna"])
  GrowthFactor = VariableCollection([Herbs], value=ValueFromParameterTable("GrowthFactor", "parameterFile.tbl", Scalar))

To assign external names use the format ``Modelname = Externalname``. Here, you create in your model script a collection variable ``GrowthFactor`` with the three indices ``Herbs.ES``, ``Herbs.SPG`` and ``Herbs.Cal``. In the parameter file you can refer to the external names ``EvergreenShrub``, ``ShortPerennialGrass`` and ``Calluna``::

  GrowthFactor  EvergreenShrub       0.3
  GrowthFactor  ShortPerennialGrass  0.37
  GrowthFactor  Calluna              0.29
