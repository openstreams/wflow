## \file
# Base class for shell scripts.



import os.path
import sys
import time
import traceback



## Base class for shell scripts.
#
# Once created the run function can be called to start the script. This
# function calls the _parseOptions and _run functions which should be
# reimplemented in the specializations (the default _run does nothing but
# returning 0).
#
# \todo      Improve message printing, using different streams
class ShellScript(object):

  ## Constructor.
  #
  # Stores the argument vector. The first element is used as the name of the
  # shell script. The remaining arguments are stored in a list which can be
  # accessed using the arguments() property.
  #
  # \param     argv List of commandline arguments.
  # \param     optionParser Parser for the command line options.
  #
  # If \a optionParser is given, the results of parsing the command line
  # arguments can be obtained by calling the options() and arguments()
  # properties.
  #
  # The \a optionParser argument should be a optparse.OptionParser object.
  # This object can be configured in a subclass before this constructor is
  # called. This way subclasses can steer the parsing of command line
  # arguments. Example:
  #
  # \code
  # import optparse
  # ...
  # def __init__(self,
  #      argv,
  #      options=[]):
  #   # Create and configure an option parser.
  #   parser = optparse.OptionParser("usage: %prog [options] arguments",
  #      option_list=options)
  #   parser.add_option("-f", "--source",
  #      dest="source", help="read from SOURCE")
  #   parser.add_option("-v", "--verbose",
  #      action="store_true", dest="verbose", default=False,
  #      help="print more status messages")
  #
  #   # Call the constructor of the base class with our parser.
  #   shellscript.ShellScript.__init__(self, argv, parser)
  #
  # def _run(self):
  #   # Access parsed options.
  #   if self.options.source:
  #     # ...
  # \endcode
  #
  # See $DEVENV/Templates/sources/script.py for an empty script that you can
  # use to kick-start the implementation of a new shell script.
  def __init__(self,
         argv,
         optionParser=None):
    self._name = os.path.basename(argv[0])
    self._options = None
    self._arguments = argv[1:]
    if optionParser:
      (self._options, self._arguments) = optionParser.parse_args(self.arguments)
    self._startTime = time.time()

  def _getName(self):
    return self._name

  def _getArguments(self):
    return self._arguments

  def _getOptions(self):
    return self._options

  ## Name property, read-only.
  name = property(_getName)

  ## Arguments property, read-only.
  #  Note: usualy equals sys.argv[1:], so 1 less than sys.argv
  arguments = property(_getArguments)

  ## Options property, read-only.
  options = property(_getOptions)

  def _formatMessage(self,
         message,
         prefix=""):
    messages = message.split("\n")
    for i in range(0, len(messages)):
      if prefix:
        messages[i] = "%s: %s" % (prefix, messages[i])
    return "%s\n" % ("\n".join(messages)).encode("utf-8")

  ## Prints a message on stdout.
  #
  # \param     message Message to print.
  #
  # - Message is stripped.
  # - A newline is appended.
  def showMessage(self, message):
    sys.stdout.write(self._formatMessage(message))

  ## Prints a warning message on stdout.
  #
  # \param     message Message to print.
  #
  # - Message is prepended by "Warning: ".
  # - Message is stripped.
  # - A newline is appended.
  def showWarning(self, message):
    sys.stdout.write(self._formatMessage(message, "Warning"))

  ## Prints an error message on stderr.
  #
  # \param     message Message to print.
  #
  # - Message is prepended by "Error: ".
  # - Message is stripped.
  # - A newline is appended.
  def showError(self, message):
    sys.stderr.write(self._formatMessage(message, "Error"))

  ## Prints a message on stdout.
  #
  # \param     message Message to print.
  #
  # The stream is flushed, so the output is imediately visible.
  def write(self, message):
    sys.stdout.write(message.strip())
    sys.stdout.flush()

  ## Wraps the call to _parseOptions() and _run() by exception handling code.
  #
  # Returns the result of _run() if it is not None, otherwise 0 is returned.
  def run(self):
    result = 1

    try:
      self._parseOptions()
      result = self._run()
    except KeyboardInterrupt:
      self.showMessage("Interrupted")
    except AssertionError:
      # TODO provide stream to write on.
      traceback.print_exc()
    except Exception:
      # self.showError(str(exception))
      traceback.print_exc()
    except:
      self.showError("Unexpected error: %s" % (sys.exc_info()[0]))

    # Functions return None if nothing is explicitly returned.
    if result is None:
      result = 0

    return result

  def _parseOptions(self):
    pass

  def _run(self):
    return 0




  ##
  # \deprecated Use arguments().
  def argv(self):
    print("deprecated, update your script")
    return [self.name] + self.arguments

  ##
  # \deprecated Use len(arguments()).
  def nrArguments(self):
    print("deprecated, update your script")
    return len(self.argv())

  ##
  # \deprecated Use name().
  def appName(self):
    print("deprecated, update your script")
    return self.name

  ##
  # \deprecated Not crucial for a basic shellscript class.
  def duration(self):
    print("deprecated, update your script")
    self._endTime = time.time()
    return "%s" % (utils.duration(self._startTime, self._endTime))

  ##
  # \deprecated Not crucial for a basic shellscript class.
  def printDuration(self):
    print("deprecated, update your script")
    self._endTime = time.time()
    self.showMessage("duration %s: %s" % (self.name,
         utils.duration(self._startTime, self._endTime)))

  ##
  # \deprecated Not crucial for a basic shellscript class.
  def preRun(self):
    print("deprecated, update your script")

  ##
  # \deprecated Not crucial for a basic shellscript class.
  def postRun(self):
    print("deprecated, update your script")

  ##
  # \deprecated Use showMessage().
  def message(self, msg):
    print("deprecated, update your script")
    self.showMessage(msg)

  ##
  # \deprecated Use showWarning().
  def warning(self, msg):
    print("deprecated, update your script")
    self.showWarning(msg)

  ##
  # \deprecated Use showError().
  def error(self, msg):
    print("deprecated, update your script")
    self.showError(msg)

