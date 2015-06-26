#!/usr/bin/env python
# -*- coding: utf-8 -*-
import inspect
import re

class Index(object):
  """
  This class introduces the index type

  Parameters:

  values
    list of array-index names. The names must be CamelCased

  Example:

  PlantSpecies = Index([ "Species1", "Species2", "Species3" ])

  This introduces the index type 'PlantSpecies'. PlantSpecies (the type)
  has three names of type array-index: Species1, Species2 and Species3.
  """
  def __init__(self, values):
    if len(values) == 0:
      raise AttributeError("Error in initialisation of class Index: no array indices provided")

    externalNames = {}
    counter = 0
    for val in values:
      if val.count("=") == 1:
        modelName, sep, externalName = val.partition("=")
        values[counter] = modelName.strip()
        externalNames[modelName.strip()] = externalName.strip()
      elif val.count("=") > 1:
        msg = "Error in initialisation of class Index: format of %s does not match Modelname = Externalname" % (val)
        raise AttributeError(msg)

      counter += 1

    self.__dict__["_values"] = values
    self.__dict__["_externalNames"] = externalNames


    for v in values:
      if not self.__dict__.get(v) == None:
        raise AttributeError("Error in initialisation of class Index: array indices must be unique, %s already used" % (v))
      self.__dict__[v] = v


  def __setattr__(self, name, value):
    raise ValueError("Modification of an Index attribute not permitted")


  def __delattr__(self, value):
    raise ValueError("Removal of an Index attribute not permitted")


  def __getitem__(self, index):
    return self._values[index]


  def _getExternalName(self, variable):
    """
    returns name if var has an external name; None if not
    """
    return self.__dict__["_externalNames"].get(variable)


  def __repr__(self):
    msg = ""
    for v in self.__dict__["_values"]:
      msg += v + " "
    return msg
