# -*- coding: utf-8 -*-



class DynamicBase(object):
  """
  Base class for DynamicModel framework class.

  .. todo::

    KDJ: What is the reason this class exists? Can it be merged with
    DynamicModel?
  """

  def __init__(self):
    if self.__class__ is DynamicBase:
      raise NotImplementedError

    self._d_nrTimeSteps = 0
    self.currentStep = 0
    self._d_firstTimeStep = 1
    self.inTimeStep = False
    self.inInitial = False
    self.inDynamic = False

  def setDebug(self):
    msg = "Class needs to implement 'setDebug' method"
    raise NotImplementedError(msg)

  def initial(self):
    """  """
    msg = "Class needs to implement 'initial' method"
    raise NotImplementedError(msg)

  def dynamic(self):
    """  """
    msg = "Class needs to implement 'dynamic' method"
    raise NotImplementedError(msg)

  def timeSteps(self):
    """  """
    msg = "Class needs to implement 'timeSteps' method"
    raise NotImplementedError(msg)

  def nrTimeSteps(self):
    """  """
    msg = "Class needs to implement 'nrTimeSteps' method"
    raise NotImplementedError(msg)

  def firstTimeStep(self):
    """
    Return the first timestep that is executed.
    """
    msg = "Class needs to implement 'firstTimeStep' method"
    raise NotImplementedError(msg)

  def setQuiet(self,
    quiet=True):
    """
    Disables the progress display of timesteps.
    """
    self.silentModelOutput = quiet

  def currentTimeStep(self):
    """  """
    msg = "Class needs to implement 'currentTimeStep' method"
    raise NotImplementedError(msg)

  def _inDynamic(self):
    return self.inDynamic

  def _inInitial(self):
    return self.inInitial


  def _setInInitial(self, value):
    assert isinstance(value, bool)
    self.inInitial = value

  def _setInDynamic(self, value):
    assert isinstance(value, bool)
    self.inDynamic = value





  def _inTimeStep(self):
    """
    Returns whether a time step is currently executing.
    """
    #if hasattr(self._userModel(), "_d_inTimeStep"):
    #  return self._userModel()._d_inTimeStep
    #else:
    #  return False
    return self.inTimeStep

  def _setInTimeStep(self, value):
    assert isinstance(value, bool)
    self.inTimeStep = value


  def _setFirstTimeStep(self, firstTimeStep):

    if not isinstance(firstTimeStep, int):
      msg = "first timestep argument (%s) of DynamicFramework must be of type int" % (type(firstTimeStep))
      raise AttributeError(msg)

    if firstTimeStep <= 0:
      msg = "first timestep argument (%s) of DynamicFramework must be > 0" % (firstTimeStep)
      raise AttributeError(msg)

    if firstTimeStep > self.nrTimeSteps():
      msg = "first timestep argument (%s) of DynamicFramework must be smaller than given last timestep (%s)" % (firstTimeStep, self.nrTimeSteps())
      raise AttributeError(msg)

    self._d_firstTimeStep = firstTimeStep

  def _setCurrentTimeStep(self, step):

    if step <= 0:
      msg = "Current timestep must be > 0"
      raise AttributeError(msg)

    if step > self.nrTimeSteps():
      msg = "Current timestep must be <= %d (nrTimeSteps)"
      raise AttributeError(msg)

    self.currentStep = step

  def _setNrTimeSteps(self, lastTimeStep):
    """
    Set the number of time steps to run.
    """
    if not isinstance(lastTimeStep, int):
      msg = "last timestep argument (%s) of DynamicFramework must be of type int" % (type(lastTimeStep))
      raise AttributeError(msg)

    if lastTimeStep <= 0:
      msg = "last timestep argument (%s) of DynamicFramework must be > 0" % (lastTimeStep)
      raise AttributeError(msg)

    self._d_nrTimeSteps = lastTimeStep

