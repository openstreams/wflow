import math, os, os.path, re
import stat
import string
import sys
import time, types, warnings

warnings.filterwarnings(u"ignore",
         u"tempnam is a potential security risk to your program",
         RuntimeWarning)



Error = 'Exception raised in pcr.py library'

class Exception(Exception):
  def __init__(self, message):
    Exception.__init__(self, message)

# copy all output printed also to file
# like unix tee(1)
# also calls flush at each print
class StdoutTee:
  def __init__(self, fileName,openMode):
      self.d_orgStdout = sys.stdout
      self.d_toFile    = open(fileName,openMode)
      sys.stdout = self
  def write(self, string):
      self.d_toFile.write(string)
      self.d_orgStdout.write(string)
      self.flush()
  def flush(self):
      self.d_toFile.flush()
      self.d_orgStdout.flush()


# Some constants.
EXEMODE  = 0o755
DIRMODE  = 0o755
FILEMODE = 0o444



def message(msg):
  for m in string.split(msg, '\n'):
    if len(string.strip(m)) > 0:
      sys.stdout.write(m + '\n')
      sys.stdout.flush()

def warning(msg):
  for m in string.split(msg, '\n'):
    if len(string.strip(m)) > 0:
      sys.stdout.write('warning: ' + m + '\n')
      sys.stdout.flush()

def error(msg):
  for m in string.split(msg, '\n'):
    if len(string.strip(m)) > 0:
      sys.stderr.write('error: ' + m + '\n')
      sys.stderr.flush()

def write(msg):
   sys.stdout.write(msg)
   sys.stdout.flush()

def testNonEmptyFile(p):
    if os.path.getsize(p) == 0:
      msg = "%s: is 0 bytes" % (p)
      raise Exception(msg)

# CW BTW, there is os.path.walk(..) with
#  similar functionality (see end of file)
def processFiles(dir, func):
  for f in os.listdir(dir):
    p = os.path.join(dir,f)
    # somehow some links do not exists
    # try without this exist claus .. in your
    # home directory
    if not (os.path.exists(p)):
      continue
    mode = os.stat(p)[stat.ST_MODE]
    if stat.S_ISDIR(mode):
      processFiles(p, func)
    elif stat.S_ISREG(mode):
      func(p)

def processDirs(dir, func):
  for f in os.listdir(dir):
    mode = os.stat(os.path.join(dir, f))[stat.ST_MODE]
    if stat.S_ISDIR(mode):
      processFiles(os.path.join(dir, f), func)
    func(os.path.join(dir, f))

def dateString():
  return time.strftime('%y%m%d', time.localtime(time.time()))

def timeString():
  return time.strftime('%y%b%d_%H%M%S', time.localtime(time.time()))

def createArchive(dir):
  currentPathName = os.getcwd()
  newPathName = os.path.split(dir)[0]

  os.chdir(newPathName)
  if regex.match('.*x.*', sys.platform) > 0:
    archiveName = os.path.join(currentPathName, os.path.basename(dir) + \
                   '.tar.gz')
    command = 'tar zcf ' + archiveName + ' ' + os.path.basename(dir) + \
                   ' --preserve'
  else:
    archiveName = os.path.join(currentPathName, os.path.basename(dir) + \
                   '.zip')
    command = 'pkzip -add -rec ' + archiveName + ' ' + os.path.basename(dir)

  os.system(command)
  os.chdir(currentPathName)

# substitute a word with different casings
# caseType can be a list of chars: c(orrect case) u(pper case) l(ower case)
def replaceCase(caseTypes, str, wordToReplace, replaceWithWord):
  newStr = str;
  for i in caseTypes.lower():
   if i == 'c':
     newStr = newStr.replace(wordToReplace,replaceWithWord)
   elif i == 'u':
     newStr = newStr.replace(wordToReplace.upper(),replaceWithWord.upper())
   elif i == 'l':
     newStr = newStr.replace(wordToReplace.lower(),replaceWithWord.lower())
   else:
     assert i in "ulc"

  return newStr

# remove all trailing and leading space and
# change all in between space to a single space
# change all occurences of 1 or more newlines
# to nl string, default is a newline
def normalizeText(str,nl="\n"):
       lines = string.split(str,"\n")
       l   = []
       for i in lines:
         str = " ".join(string.split(i))
         if (str != ""):
           l.append(str)
       return nl.join(l)

# Returns 1 if str contains whitespace characters.
def containsWhiteSpace(str):
  # split default splits at any whitespace character.
  return len(str) > 0 and \
    (str[0].isspace() or str[len(str) - 1].isspace() or len(str.split()) > 1)

# no args return value of PCRTREE environment variable
# with arg: return os.path.join(PCRTREE,arg0,arg1,...)
def pcrtree( *path ):
  def getEnv():
    assert("PCRTREE2" in os.environ)
    return os.environ['PCRTREE2']
    # if os.name == 'posix':
    #   assert(os.environ.has_key('PCRTREE'))
    #   return os.environ['PCRTREE']
    # elif os.name == 'nt':
    #   assert(os.environ.has_key('PCRTREE_NTPATH'))
    #   return os.environ['PCRTREE_NTPATH']
    # else:
    #   assert false, 'os not supported yet'
  p=getEnv()
  for d in path:
    p = os.path.join(p,d)
  return p

# returns value of DEVENV environment variable
def devenv( *path ):
  def getEnv():
    assert("DEVENV" in os.environ)
    return os.environ['DEVENV']
  p=getEnv()
  return p

def pcrtreeRelativePath(absPath):
 """ path of absPath relative to PCRTREE env. variable:\n
     PCRTREE=/home/cees/pcrtree
     absPath=/home/cees/pcrtree/libs/pcrme/calc_test.cc
     will return libs/pcrme/calc_test.cc

    *Under Windows the return path is also ALWAYS with "/"!
    *Under Cygwin all output options are accepted for absPath
     --windows c:\d\pcrtree
     --mixed   c:/d/pcrtree
     --unix    /cygdrive/c/d/pcrtree
 """
 assert os.path.isabs(absPath)
 if absPath.find("/cygdrive")==0:
   absPath =executeOneLine("cygpath --mixed \""+absPath+"\"")
 # +1 for intervening slash
 return absPath[len(pcrtree())+1:].replace("\\","/")

def python2():
  return sys.version_info[0] == 2

def python24():
  return python2() and sys.version_info[1] == 5

if not python24():

  def execute(cmd, bufferSize = -1):
    pipe = os.popen("%s 2>&1" % (cmd), "r", bufferSize)
    output = pipe.readlines()

    if pipe.close():
      if not output:
        output = "Command '%s': error while executing" % (cmd)

      raise Exception("".join(output))

    return output

else:

  import subprocess

  # execute with no data needed on stdin
  # TODO: obsoleted, use subprocess module directly (2.4 and up)
  class ExecNoStdin:
    def __init__(self, cmd):
      p = subprocess.Popen(string.split(cmd),
                shell=True, bufsize=-1,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
      self.d_stdoutLines = p.stdout.readlines()
      self.d_stderrLines = p.stderr.readlines()
      self.d_exitStatus  = p.wait()

    def stdoutEmpty(self):
     return len(self.d_stdoutLines)==0
    def stderrEmpty(self):
     return len(self.d_stderrLines)==0
    def stdoutLines(self):
     return self.d_stdoutLines
    def stdoutText(self):
     return " ".join(self.d_stdoutLines)
    def stderrLines(self):
     return self.d_stderrLines
    def stderrText(self):
     return " ".join(self.d_stderrLines)
    def exitStatus(self):
     return self.d_exitStatus

  # wrapper around subprocess.Popen that
  #  - use temp files for stdout and stderr, due to
  #    PIPE problems on windows
  #  - executes cmd and waits for termination
  #  - throws pcr.Exception in case of OSError
  #  - sets shell to True (CW dunno if that's handy)
  #  - returns a class with the following attributes:
  #     - stdout
  #     - stderr   both a list with lines (\\n terminated) empty list if no
  #                output on that stream.
  #     - exitcode exit code of wait()
  def popenWrapper(cmd, bufferSize=0):

   class PopenWrapperReturn:
    def __init__(self,cmd,bufferSize):
     class Stream:
       def __init__(self):
        self.fileName = os.tempnam()
       def createFile(self):
        return file(self.fileName, "wb")
       def close(self):
        if os.path.exists(self.fileName):
         os.remove(self.fileName)
       def lines(self):
        if os.path.exists(self.fileName) and os.path.getsize(self.fileName):
          return file(self.fileName, "r").readlines()
        return []

     sOut = Stream()
     sErr = Stream()

     self.exitcode = -999
     self.stderr = self.stdout = []

     try:
       try:
         # linux requires string as 1st arg!
         if type(cmd) == list:
          cmdStr=" ".join(cmd)
         else:
          cmdStr=cmd
         process = subprocess.Popen(cmdStr,
              shell=True, bufsize=bufferSize,
              stdout=sOut.createFile(), stderr=sErr.createFile() )
         self.exitcode = process.wait()
       except OSError as exception:
         sOut.close()
         sErr.close()
         raise Exception(str(exception))
     finally:
       self.stdout=sOut.lines()
       self.stderr=sErr.lines()

     sOut.close()
     sErr.close()
   return PopenWrapperReturn(cmd,bufferSize)

  # expecte one line of stdout and strip it
  def executeOneLine(cmd):
    return string.strip(execute(cmd)[0])

  def execute(cmd, bufferSize=0):
    """ execute command and
         - return stdout as list of strings
         - raise  pcr.Exception on non-0 exit code or any other error condition
         - cache stderr, if pcr.Exception is raised then str() on this
           exception yield the stderr output
    """

    r = popenWrapper(cmd,bufferSize)
    if r.exitcode != 0 and r.exitcode != 1100:
      assert len(r.stderr)
      raise Exception("".join(r.stderr))
    return r.stdout

# concatenate lists
def listConc(src):
 dest = []
 for i in src:
  for j in i:
   dest.append(j)
 return dest

# HANDIG-> if isinstance(item, ListType):

def listAppend(dest,*listsToAppend):
  for l in listsToAppend:
    for j in l:
     dest.append(j)
  return dest

# Joins all pathnames in the list pathNames and returns the result.
# The result is normalised.
def pathJoin(pathNames):
  pathName = ''
  for n in pathNames:
    pathName = os.path.join(pathName, n)
  return os.path.normpath(pathName)

# Converts pathnames to urls.
def pathNameToUrl(pathName):
  url = pathName

  if os.sep != '/':
    url = url.replace(os.sep, '/')

  return url

# Tests if a file can be opened for reading.
def testOpenForReading(pathName):
  if not os.path.exists(pathName):
    raise Exception('File \'%s\': Does not exist' % (pathName))
  mode = os.stat(pathName)[stat.ST_MODE]
  # if not stat.S_ISREG(mode):
  #   raise Exception('File \'%s\': Not a regular file' % (pathName))
  if not (mode & stat.S_IRUSR or mode & stat.S_IRGRP or mode & stat.S_IROTH):
    raise Exception('File \'%s\': No permission to read' % (pathName))

# Tests if a file can be opened for writing.
def testOpenForWriting(pathName):
  if os.path.exists(pathName):
    mode = os.stat(pathName)[stat.ST_MODE]
    # if not stat.S_ISREG(mode):
    #   raise Exception('File \'%s\': Not a regular file' % (pathName))
    if not (mode & stat.S_IWUSR or mode & stat.S_IWGRP or mode & stat.S_IWOTH):
      raise Exception('File \'%s\': No permission to write' % (pathName))

# Tests if a file doesn't already exist.
def testFileDoesNotExist(pathName):
  if os.path.exists(pathName):
    raise Exception('File \'%s\': Already exists' % (pathName))

# Search a file in a file tree. Start at the bottom. If found, this function
# returns the file name, else an empty string.
def searchFileUpTree(directory, fileToFind):
  fileNotFound = 1
  fn = os.path.join(directory, fileToFind)
  while fileNotFound and not len(fn) == 0:
    if os.path.exists(fn):
      fileNotFound = 0
      break;
    fn = os.path.split(fn)[0]
    # Test if we can go up the tree.
    if len(fn) > 0 and fn != os.sep:
      fn = os.path.join(os.path.split(fn)[0], fileToFind)

  return fn

# Returns true (not 0) if DOM node node has one child with name name.
def hasOnlyChild(node, name):
  elements = node.getElementsByTagName(name)
  return len(elements) == 1

# Returns the only child of DOM node node with name name.
def getOnlyChild(node, name):
  elements = node.getElementsByTagName(name)
  assert len(elements) == 1
  child = elements[0]
  assert child.nodeName == name
  return child

# Returns path with all environment variables expanded.
def expandEnvVars(path):
  # Search and replace variables.
  pattern = '\$[A-Z]+'
  expression = re.compile(pattern)

  start = 0;
  end = len(path)

  match = re.search(expression, path)
  while match != None:
    variable = match.group()
    if not variable[1:] in os.environ:
      start = start + len(variable)
    else:
      value = os.environ[variable[1:]]
      path = path.replace(variable, value)
      start = start + len(value)
    match = re.search(expression, path[start:end])

  return path

# Searches for executable fileName in all standard locations.
def findExecutable(fileName):
  tryName = os.path.join('/usr/local/bin', fileName)
  if os.path.exists(tryName):
    return tryName

  tryName = os.path.join(expandEnvVars('$HOME'), fileName)
  if os.path.exists(tryName):
    return tryName

  return ''

# Replaces all occurences of certain strings in file with name fileName.
# Returns new string.
def replaceFromFile(fileName, searchStrings, replaceStrings):
  testOpenForReading(fileName)
  result = ''
  for line in open(fileName, "r").readlines():
    result += replaceFromString(line, searchStrings, replaceStrings)
  return result

# Replaces all occurences of certain strings in sourceString. Returns new
# string.
def replaceFromString(sourceString, searchStrings, replaceStrings):
  assert len(searchStrings) == len(replaceStrings)
  result = sourceString
  for i in range(len(searchStrings)):
    result = result.replace(searchStrings[i], replaceStrings[i])
  return result

# Removes duplicate occurences in list. Returns new list.
def removeDuplicates(oldList):
  newList = []
  for item in oldList:
    if not item in newList:
      newList.append(item)
  return newList

# Checks if address contains a valid email address.
def isValidEmailAddress(address):
  if address.find("@") == -1:
    return 0
  return 1


# Encoding is in utf-8.
def createMailMessage(fromAddress, toAddress, subject, text):

  from email.MIMEText import MIMEText
  message = MIMEText(text.encode("utf-8"), "plain", "utf-8")
  message["Subject"] = subject
  message["From"] = fromAddress
  message["To"] = toAddress
  return message.as_string()


def log(fileName, msg):

  testOpenForWriting(fileName)
  open(fileName, 'a').write(msg)



# quite expensive way of what findDirs also does, maybe substite
def discover(dirsWithMakefile, path, names):
  # [:] iter over copy, so we changes name in loop
  for n in names[:] :
    # remove all, except those added again
    names.remove(n)
    p = os.path.join(path,n)
    # seems to happen?
    # permissions?
    if not (os.path.exists(p)):
      continue
    mode = os.stat(p)[stat.ST_MODE]
    # skip links
    if stat.S_ISLNK(mode):
      continue
    #  if stat.S_ISDIR(mode):
    # else :
# os.path.walk(os.environ['PCRTREE']+"/apps", discover, dirsWithMakefiles)

# Returns the relative path from the directory to the top.
def relPathUp(directory):
  directory = os.path.normpath(string.strip(directory))
  if not directory or directory == "." or directory.find("..") != -1:
    return directory

  path = "."
  dirNames = string.split(directory, os.sep)
  for name in dirNames:
    if name and name != "..":
      path = os.path.join(path, "..")
  return path

# Creates a directory with a specific name in a temp location and returns the
# name of the created directory.
def createTempDirectory(name):
  directoryName = os.path.join(os.tmpnam(), name)
  os.makedirs(directoryName)
  return directoryName

# filter for
#   Microsoft cl compiler to use with vim
#   pcrtree unittest in "cygdrive"-fmt
# input list of lines (ending newline optional)
# output list of reformat lines with newline char at end
def clVimFilter(msgs):
 output = []
 for l in msgs:
   l = string.strip(l)
   # empty string
   if not len(l):
     continue
   # scons start
   if not l.find("scons: Entering directory"):
     continue
   # compilation unit printed (e.g. main.c on single line)
   if len(l.split())==1:
     singleWordSplit = l.split(".")
     if len(singleWordSplit) >= 2 and singleWordSplit[1] in [ 'c','cc','cpp','res']:
      continue

   def reformat(l):
     def absCygFile(l,afterFile):
         fileName    = l[0:afterFile]
         fileName    = l[0:afterFile]
         if not os.path.isabs(fileName):
           fileName = pcrtree(fileName)
         return executeOneLine("cygpath --unix \""+fileName+"\"")

     afterFile = l.find("(")
     e = l.find(") :")
     if afterFile != -1 and e != -1:
       lineNoStr = l[afterFile+1:e]
       if lineNoStr.isdigit():
          # like col2map.c(85) : warning C4244: etc.
          #  sep counter:    123
         fileName=absCygFile(l,afterFile)
         msg         = l[e+4:]
         return "%s:%s::%s" % (fileName,lineNoStr,msg)
     afterFile = l.find(":")
     if afterFile != -1:
        afterLineNr= l.find(":",afterFile+1)
        if afterLineNr != -1:
          lineNoStr = l[afterFile+1:afterLineNr]
          if lineNoStr.isdigit():
            # like col2map_col2maptest.cc:101: TODO: etc.
            fileName    = absCygFile(l,afterFile)
            msg         = l[afterLineNr+1:]
            return "%s:%s::%s" % (fileName,lineNoStr,msg)
     # as is
     return l

   # OK :cl :cn works
   # :cw gives me all messages in window
   output.append(reformat(l))

 return "\n".join(output)
