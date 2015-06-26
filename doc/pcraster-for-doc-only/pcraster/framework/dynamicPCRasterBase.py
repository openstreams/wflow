# -*- coding: utf-8 -*-
import _pcraster as pcr
import dynamicBase



class DynamicModel(dynamicBase.DynamicBase):

  def __init__(self):
    dynamicBase.DynamicBase.__init__(self)
    self.silentModelOutput = False

  def initial(self):
    print "Implement 'initial' method"

  def dynamic(self):
    print  "Implement 'dynamic' method"

  def _silentModelOutput(self):
    return self.silentModelOutput

  def timeSteps(self):
    """
    Return a list of time steps configured
    """
    return range(self.firstTimeStep(), self.nrTimeSteps() + 1)

  def nrTimeSteps(self):
    """
    Return the number of time steps
    """
    assert self._d_nrTimeSteps
    return self._d_nrTimeSteps

  def currentTimeStep(self):
    """
    Return the current time step in the range from firstTimeStep to nrTimeSteps.
    """
    assert self.currentStep >= 0
    return self.currentStep

  def firstTimeStep(self):
    """
    Return first timestep of a model.
    """
    assert self._d_firstTimeStep
    return self._d_firstTimeStep

  def report(self,
    variable,
    name):
    """
    Storing map data to disk

    `variable`
      Variable containing the PCRaster map data

    `name`
      Name used as filename. Use a filename with less than eight
      characters and without extension. File extension for dynamic models
      is ".map" in the initial section and the 8.3 style format name in
      the dynamic section. File extensions will be appended automatically.
    """
    self._reportNew(variable, name)

  def readmap(self, name, style=1):
    """
    Read map data from disk.

    `name`
      Name used as filename. Use filename with less than eight characters
      and without extension. File extension for dynamic models is ".map"
      in initial section and the 8.3 style format name in the dynamic
      section. File extensions will be appended automatically.

    .. todo::

      `style` argument is not used.
    """
    return self._readmapNew(name)

  def _setNrTimeSteps(self,
    timeSteps):
    """
    Configure the number of time steps.

    In addition to the setting the number of timesteps we need to pass
    the value to the PCRaster runtime engine.
    """
    dynamicBase.DynamicBase._setNrTimeSteps(self, timeSteps)

    pcr._rte().setNrTimeSteps(timeSteps)

  def _setCurrentTimeStep(self,
    step):
    """
    Set the current time step.

    In addition to the setting the current timestep within the framework,
    we need to pass the value to the PCRaster runtime engine.
    """
    dynamicBase.DynamicBase._setCurrentTimeStep(self, step)

    pcr._rte().setCurrentTimeStep(step)

