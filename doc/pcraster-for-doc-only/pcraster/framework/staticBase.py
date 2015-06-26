# -*- coding: utf-8 -*-



class StaticBase(object):
  """
  Base class for StaticModel framework class.

  .. todo::

    KDJ: What is the reason this class exists? Can it be merged with
    StaticModel?
  """

  def __init__(self):
    if self.__class__ is StaticBase:
      raise NotImplementedError

    self.inInitial = False

  def initial(self):
    """  """
    msg = "Class needs to implement 'initial' method"
    raise NotImplementedError(msg)

  def setDebug(self):
    """  """
    msg = "Class needs to implement 'setDebug' method"
    raise NotImplementedError(msg)

  def _inInitial(self):
    return self.inInitial

  def _setInInitial(self,
    value):
    assert isinstance(value, bool)
    self.inInitial = value

