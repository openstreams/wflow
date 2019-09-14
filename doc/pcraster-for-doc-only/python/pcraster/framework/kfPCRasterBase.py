# -*- coding: utf-8 -*-
from . import kfBase



class EnKfModel(kfBase.EnKfBase):
  def __init__(self):
    kfBase.EnKfBase.__init__(self)

  ## \brief Storing map data to disk.
  #
  # \param variable object containing the PCRaster map data
  # \param name name used as filename. Use filename without extension. File extension for deterministic
  # static models is ".map" and will be appended automatically.
  def report(self, variable, name, style=1):
    self._reportNew(variable, name)

  ## \brief Reads map data from disk.
  #
  # \param name name used as filename. Use filename without extension. File extension for deterministic
  # static models is ".map" and will be appended automatically.
  def readmap(self, name, style=1):
    return self._readmapNew(name)

