# -*- coding: utf-8 -*-
import staticBase



class StaticModel(staticBase.StaticBase):

  def __init__(self):
    """
    Base class for classes implementing the
    :ref:`Static Model Concept <staticModelConcept>`.
    """
    staticBase.StaticBase.__init__(self)

  def initial(self):
    """
    Reimplemented from StaticBase.
    """
    print("Implement 'initial' method")

  def report(self,
    variable,
    name,
    style=1):
    """
    Store map data to disk

    `variable`
      object containing the PCRaster map data

    `name`
      name used as filename. Use filename without extension. File extension
      for deterministic static models is ".map" and will be appended
      automatically.

    .. todo::

      `style` argument is not used.
    """
    self._reportNew(variable, name)

  def readmap(self,
    name,
    style=1):
    """
    Read map data from disk.

    `name`
      name used as filename. Use filename without extension. File extension
      for deterministic static models is ".map" and will be appended
      automatically.

    .. todo::

      `style` argument is not used.
    """
    return self._readmapNew(name)

