import os, string, re
import pcr

def nrCPUs():
  nrCPUs = 0
  pattern = "^processor[ \t]*:[ \t][0-9]+$"
  matchObject = re.compile(pattern)
  for line in file("/proc/cpuinfo", "r").readlines():
    if matchObject.match(line):
      nrCPUs += 1
  assert nrCPUs > 0
  return nrCPUs

# Returns average system load for last minute as a float. A load of 100 means
# one processor is fully occupied.
def load():
  return 100.0 * float(os.getloadavg()[0])

def loadPerCPU():
  return load() / nrCPUs()
