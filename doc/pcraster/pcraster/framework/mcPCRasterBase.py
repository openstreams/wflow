# -*- coding: utf-8 -*-
import pcraster
import frameworkBase
import mcBase



class MonteCarloModel(mcBase.MonteCarloBase):

  def __init__(self):
    mcBase.MonteCarloBase.__init__(self)

  def premcloop(self):
    print "premcloop not implemented"

  def postmcloop(self):
    print "postmcloop not implemented"

  def report(self,
    variable,
    name):
    """
    Report map data to disk.

    Standard extension is "map" in initial and timestep in the dynamic section.
    Output directory will be sample directory.
    """
    self._reportNew(variable, name)

  def readmap(self,
    name):
    """
    Read sample data from disk.

    Returns the map of the current time step from the current sample directory.
    """
    return self._readmapNew(name)

  def readDeterministic(self,
    name):
    """
    Read deterministic data from disk.

    Returns the map of the file with current time step, from the current
    working directory.
    """
    if self._inPremc() or self._inPostmc() or self._inInitial():
      newName = name + ".map"
    else:
      newName = frameworkBase.generateNameT(name, self.currentTimeStep())

    return pcraster.readmap(newName)

