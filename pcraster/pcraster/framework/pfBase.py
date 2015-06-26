# -*- coding: utf-8 -*-



class ParticleFilterBase(object):
  def __init__(self):
    if self.__class__ is ParticleFilterBase:
      raise NotImplementedError

    self._d_filterPeriod = 0
    self._d_inFilterPeriod = False
    self._d_filterTimesteps = []
    self._d_inUpdateWeight = False

  def setDebug(self):
    msg = "Class needs to implement 'setDebug' method"
    raise NotImplementedError(msg)

  def filterPeriod(self):
    """
    Return the current filter period.
    """
    #assert hasattr(self._userModel(), "_d_filterPeriod")
    return self._d_filterPeriod + 1
  #def filterPeriod(self):
    #msg = "Class needs to implement 'filterPeriod' method"
    #raise NotImplementedError(msg)

  def filterTimesteps(self):
    msg = "Class needs to implement 'filterTimesteps' method"
    raise NotImplementedError(msg)

  def particleWeight(self):
    msg = "Class needs to implement 'particleWeight' method"
    raise NotImplementedError(msg)

  def reportState(self):
    msg = "Class needs to implement 'reportState' method"
    raise NotImplementedError(msg)

  def readState(self):
    msg = "Class needs to implement 'readState' method"
    raise NotImplementedError(msg)

