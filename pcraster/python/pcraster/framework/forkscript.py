"""Functionality for applications with sub processes."""

import copy, os, signal, time, sys



class ChildProcess:
  """Class to represent child processes."""

  def __init__(self, name, pid):
    """Constructor.

    name -- Name of the process.

    pid -- Process id of the process.

    Upon creation the start time of the process will be set."""
    self.d_name = name
    self.d_pid = pid
    self.d_startTime = time.time()
    self.d_endTime = None

  def stopped(self, result):
    """Sets the end time and the result of the process.

    result -- Result.

    Should be called after the process stopped executing."""
    self.d_endTime = time.time()
    self.d_result = result

  def pid(self):
    """Returns the process id."""
    return self.d_pid

  def name(self):
    """Returns the name."""
    return self.d_name

  def exitedWithZeroStatus(self):
    """Returns whether or not the process exited with status 0."""
    return os.WIFEXITED(self.d_result) and os.WEXITSTATUS(self.d_result) == 0

  def duration(self):
    """Returns the duration of the process as an amount of seconds."""
    return self.d_endTime - self.d_startTime

  def durationAsString(self):
    """Returns the duration of the process as a string."""
    minsec = divmod((self.duration()), 60)
    return "%d min %d sec" % (int(minsec[0]), int(minsec[1]))

  def message(self):
    """Returns a message with result information of the process.

    Should only be called after the process stopped executing."""
    message = ""
    if os.WIFSTOPPED(self.d_result):
      message = "Child %s with process id %d stopped with signal %d" \
            % (self.d_name, self.d_pid, os.WSTOPSIG(self.d_result))
    elif os.WIFSIGNALED(self.d_result):
      message = "Child %s with process id %d exited due to signal %d" % \
            (self.d_name, self.d_pid, os.WTERMSIG(self.d_result))
    elif os.WIFEXITED(self.d_result):
      message = "Child %s with process id %d exited with parameter %d" % \
            (self.d_name, self.d_pid, os.WEXITSTATUS(self.d_result))
    assert message
    message += "\nDuration (indication): %s" % (self.durationAsString())
    return message



class ForkScript:

  def __init__(self):
    self.d_parentProcess = True
    self.d_childProcesses = []
    self.d_nrForkedChildProcesses = 0
    self.d_nrFinishedChildProcesses = 0
    self.d_nrFailedChildProcesses = 0
    self.d_childProcessDurations = []
    self.d_decayFactor = 0.5

  def __del__(self):
    if self.isParentProcess():
      for process in self.childProcesses():
        os.kill(process.pid(), signal.SIGKILL)

  def isParentProcess(self):
    return self.d_parentProcess

  def isChildProcess(self):
    return not self.isParentProcess()

  def childProcesses(self):
    return self.d_childProcesses

  def nrChildProcesses(self):
    return len(self.d_childProcesses)

  def nrForkedChildProcesses(self):
    return self.d_nrForkedChildProcesses

  def nrFinishedChildProcesses(self):
    return self.d_nrFinishedChildProcesses

  def nrFailedChildProcesses(self):
    return self.d_nrFailedChildProcesses

  def setSleepTimeDecayFactor(self,
         decayFactor):
    self.d_decayFactor = decayFactor

  def fork(self, name):
    pid = os.fork()
    assert pid != -1
    if pid == 0:
      self.d_parentProcess = False
      # set seeds different for each child
      if "random" in sys.modules:
        import random
        random.seed()
      if "numpy.random" in sys.modules:
        import numpy.random
        numpy.random.seed()
      if "pcraster" in sys.modules:
        import pcraster
        pcraster.setrandomseed(0)
      # # for backwards compatibility
      # if "pcraster" in sys.modules:
      #   import PCRasterPython
      #   PCRasterPython.setrandomseed(0)
    else:
      self.d_childProcesses.append(ChildProcess(name, pid))
      self.d_nrForkedChildProcesses += 1
    return pid

  def handleFinishedChildProcess(self, process, result):
    self.d_childProcesses.remove(process)
    process.stopped(result[1])
    if process.exitedWithZeroStatus():
      self.d_childProcessDurations.append(process.duration())
    else:
      self.d_nrFailedChildProcesses += 1
    self.d_nrFinishedChildProcesses += 1

  def handleFinishedChildProcesses(self):
    """Performs some administratives in case child processes have finished."""
    assert self.isParentProcess()
    finishedChildProcesses = []
    i = 0
    while i < len(self.childProcesses()):
      process = self.childProcesses()[i]
      result = os.waitpid(process.pid(), os.WNOHANG)
      if result[0] > 0:
        self.handleFinishedChildProcess(process, result)
        finishedChildProcesses.append(process)
        i -= 1
      i += 1
    return finishedChildProcesses

  def waitForChildProcessesToFinish(self):
    """Waits for the current child processes to finish.

    Also performs some administratives.
    """
    assert self.isParentProcess()
    finishedChildProcesses = []
    while self.childProcesses():
      process = self.childProcesses()[0]
      result = os.waitpid(process.pid(), 0)
      assert result[0] > 0
      self.handleFinishedChildProcess(process, result)
      finishedChildProcesses.append(process)
    return finishedChildProcesses

  def averageChildProcessDuration(self, lastNrSamplesToUse):
    assert self.isParentProcess()
    assert len(self.d_childProcessDurations)
    nrDurations = len(self.d_childProcessDurations)
    nrDurationsToUse = min(nrDurations, lastNrSamplesToUse)
    duration = 0.0
    for i in range(0, nrDurationsToUse):
      duration += self.d_childProcessDurations[nrDurations - 1 - i]
    return duration / nrDurationsToUse

  def sleepTime(self):
    assert self.isParentProcess()
    defaultDuration = 60
    totalDuration = defaultDuration
    if self.d_childProcessDurations:
      totalDuration = self.d_decayFactor * self.averageChildProcessDuration(3)
    totalDuration = min(totalDuration, defaultDuration)
    return totalDuration

