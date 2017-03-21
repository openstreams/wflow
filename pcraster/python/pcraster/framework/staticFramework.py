# -*- coding: utf-8 -*-
import frameworkBase



class StaticFramework(frameworkBase.FrameworkBase):
  """
  Framework class for static models.

  `userModel`
    Instance that models the :ref:`Static Model Concept <staticModelConcept>`.
  """

  def __init__(self,
    userModel):
    frameworkBase.FrameworkBase.__init__(self)
    self._d_model = userModel
    self._testRequirements()

    self._addMethodToClass(self._readmapNew)
    self._addMethodToClass(self._reportNew)

  def run(self):
    """
    Re-implemented from ShellScript.

    Run the execute or initial section of the user model.
    """
    self._atStartOfScript()
    # Execute initial section.
    self._runInitial()

    return 0

  def _userModel(self):
    """
    Return the model instance provided by the user.
    """
    return self._d_model

  def _testRequirements(self):
    """
    Test whether the user model models the
    :ref:`Static Model Concept <staticModelConcept>`.
    """
    if hasattr(self._userModel(), "_userModel"):
      msg = "The _userModel method is deprecated and obsolete"
      self.showWarning(msg)

    if not hasattr(self._userModel(), "initial") \
      and not hasattr(self._userModel(), "run"):
      msg = "Cannot run static framework: Implement either an initial or a run method in the user class"
      raise frameworkBase.FrameworkError(msg)

