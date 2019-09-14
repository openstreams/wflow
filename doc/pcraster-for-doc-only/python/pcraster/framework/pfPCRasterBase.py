# -*- coding: utf-8 -*-
import os
import re
import pcraster
from . import frameworkBase
from . import pfBase



class ParticleFilterModel(pfBase.ParticleFilterBase):

  def __init__(self):
    pfBase.ParticleFilterBase.__init__(self)

  def reportState(self, variable, variableName):
    """
    Report a map into the state variable directory.
    """
    sample = str(self.currentSampleNumber())
    if re.search(".map", variableName):
      filename = variableName
    else:
      filename = frameworkBase.generateNameT(variableName,
        self.currentTimeStep())
    name = os.path.join(sample, "stateVar", filename)
    pcraster.report(variable, name)

  def readState(self, variableName):
    """
    Read a state variable map.
    """
    sample = str(self.currentSampleNumber())
    if re.search(".map", variableName):
      filename = variableName
    else:
      timestep = self.firstTimeStep() - 1
      filename = frameworkBase.generateNameT(variableName, timestep)
    name = os.path.join(sample, "stateVar", filename)
    return pcraster.readmap(name)

