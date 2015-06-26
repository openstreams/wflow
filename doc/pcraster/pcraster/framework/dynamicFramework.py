# -*- coding: utf-8 -*-
import sys
import frameworkBase


class DynamicFramework(frameworkBase.FrameworkBase):
  """
  Framework class for dynamic models.

  `userModel`
    Instance that models the :ref:`Dynamic Model Concept <dynamicModelConcept>`.

  `lastTimeStep`
    Last timestep to run.

  `firstTimestep`
    Sets the starting timestep of the model (optional, default is 1).
  """

  def __init__(self,
    userModel,
    lastTimeStep=0,
    firstTimestep=1):
    frameworkBase.FrameworkBase.__init__(self)

    self._d_model = userModel
    self._testRequirements()

    # fttb
    self._addMethodToClass(self._readmapNew)
    self._addMethodToClass(self._reportNew)

    try:
      self._userModel()._setNrTimeSteps(lastTimeStep)
      self._d_firstTimestep = firstTimestep
      self._userModel()._setFirstTimeStep(self._d_firstTimestep)
    except Exception, msg:
      sys.stderr.write('Error: %s\n' % str(msg))
      sys.exit(1)

  def _userModel(self):
    """
    Return the model instance provided by the user.
    """
    return self._d_model

  def run(self):
    """
    Run the dynamic user model.

    .. todo::

      This method depends on the filter frameworks concept. Shouldn't its run
      method call _runSuspend()?
    """
    self._atStartOfScript()
    if(hasattr(self._userModel(), "resume")):
      if self._userModel().firstTimeStep() == 1:
        self._runInitial()
      else:
        self._runResume()
    else:
      self._runInitial()

    self._runDynamic()

    # Only execute this section while running filter frameworks.
    if hasattr(self._userModel(), "suspend") and \
      hasattr(self._userModel(), "filterPeriod"):
      self._runSuspend()

    return 0

  def _testRequirements(self):
    """
    Test whether the user model models the
    :ref:`Dynamic Model Concept <dynamicModelConcept>`.
    """
    if hasattr(self._userModel(), "_userModel"):
      msg = "The _userModel method is deprecated and obsolete"
      self.showWarning(msg)

    if( not hasattr(self._userModel(), "dynamic") and not hasattr(self._userModel(), "run")):
      msg = "Cannot run dynamic framework: Implement dynamic method"
      raise frameworkBase.FrameworkError(msg)


    if not hasattr(self._userModel(), "initial"):
      if self._debug():
        self.showWarning("No initial section defined.")

  def setQuiet(self,
    quiet=True):
    """
    Disables the progress display of timesteps.
    """
    self._d_quietProgressDots = quiet

