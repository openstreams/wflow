# -*- coding: utf-8 -*-

from . import frameworkBase


class EnKfBase(object):
  def __init__(self):
    if self.__class__ is EnKfBase:
      raise NotImplementedError

  def initial(self):
    msg = "Class needs to implement 'initial' method"
    raise NotImplementedError(msg)

  def setDebug(self):
    msg = "Class needs to implement 'setDebug' method"
    raise NotImplementedError(msg)

  def setState(self):
    msg = "Class needs to implement 'setState' method"
    raise NotImplementedError(msg)

  def setObservations(self):
    msg = "Class needs to implement 'setObservations' method"
    raise NotImplementedError(msg)

