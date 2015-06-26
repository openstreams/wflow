"""
Contains general utilities that don't belong in one of the *utils.py modules.

@todo: This module should be empty. In general utilities must be classified and put in a thematic *utils.py module. I (CW) vote for systemutils.py; they are all wrappers around os,process and sys.
"""
import os
import subprocess
import sys


def callWithRaise(command):
  """
     Executes the command like call() does but then raise OSError it it fails.
  """
  status = call(command)
  if status:
     raise OSError("Failure on cmd: "+command)
  return status

def call(command):
  """
  Executes the command and returns the result. Output of the command will be printed to stdout and/or stderr as it becomes available. If you need to process these outputs use the execute function.

  @param command: Command to execute.
  @return: status
  """
  # if sys.platform != "win32":
  #   command = command.split()

  return subprocess.call(command, shell=True)

def execute(
         command,
         input=None):
  """
  Executes the command and returns the results.

  Calling this function can result in three states (assuming the command passed in adheres to the conventions):
    1. Command exists, usage is good: status will be 0, outLines will contain the output and errLines will be empty.
    2. Command exists, wrong usage: status will be 2, outLines will be empty and errLines will contain the error message.
    3. Command does not exist: status will be 1, outLines and errLines will be empty and an OSError exception will be thrown.

  @param command: Command to execute.
  @return: status, outLines, errLines.
  @rtype: Tuple with the three elements.
  @raise OSError: In case the command doesn't exist for example.
  """
  # if sys.platform != "win32":
  #   command = command.split()

  noneOrPipe = None

  if input:
    noneOrPipe = subprocess.PIPE

  child = subprocess.Popen(command,
         shell=True,
         stdin=noneOrPipe,
         stdout=subprocess.PIPE,
         stderr=subprocess.PIPE)
  out, err = child.communicate(input)
  status = child.returncode

  return status, out, err

def nativePath(
         name):
  """
  Converts the path name passed in to a native representation.

  @param name: Path name.
  @return: Converted path name.

  The path is normalized.

  TODO base this on the characteristics of the python interperter, not on
       platform. If cygwin's python is used, the path should be different
       than when window's python is used. Does platform take that into account?
  """
  result = name

  if sys.platform == "win32":
    # Also converts cygwin links to windows paths, so python for windows
    # can find the directory.
    result = execute("cygpath -m \"%s\"" % (name))[1].rstrip()

  return os.path.normpath(result)

def environmentVariableAsNativePath(environmentVariable):
  """
  Assumes environmentVariable to contain a path and
  returns the value as nativePath

  On Windows the value is converted according to the Windows pathnaming conventions.
  """
  assert(os.environ.has_key(environmentVariable)), environmentVariable
  return nativePath(os.environ[environmentVariable])

def duration(
         startTime,
         endTime):
  """
  Returns a formatted string with the duration between endTime and startTime.

  Use time.time() to create startTime and endTime.

  @param startTime: A date before endTime.
  @param endTime: A date after startTime.
  """
  minsec = divmod((endTime - startTime), 60)
  return "%s min %s sec" % (str(int((minsec[0]))), str(int((minsec[1]))))



  # def systemIsOccupied(self):
  #   print "deprecated, update your script"
  #   return system.loadPerCPU() > 150.0

