"""
This is the PCRaster multicore package.


PCRaster is a collection of tools and software libraries tailored to the
construction of spatio-temporal models. Information about the development
of PCRaster and its application in environmental modelling can be found
at http://www.pcraster.eu

The multicore module provides alternative, multi-threaded implementations
of several PCRaster operations. The following operators and operations
are able to exploit multiple processors: ::

  + - * /
  < <= > >= == !=
  **
  & | ! ^
  defined, cover
  ifthen, ifthenelse
  max, min
  boolean, nominal, ordinal, scalar
  sqr, sqrt, abs, fac, ln, log10
  rounddown, roundup, roundoff
  cos, sin, tan, acos, asin, atan
  slope, window4total, windowtotal, windowaverage
  mapmaximum, mapminimum

To set or query the number of worker threads use: ::

  set_nr_worker_threads, nr_worker_threads
"""
import os
import sys

# On Windows prepend the path to our dlls to the PATH environment variable.
# Otherwise our dlls won't be found when our Python extensions are loaded
# by Python.
if sys.platform == "win32":
    path_ = os.environ["PATH"]
    pcraster_installation_root = os.path.abspath(os.path.join(
        os.path.dirname(__file__), "..", "..", ".."))
    pcraster_dll_pathname = os.path.join(pcraster_installation_root, "lib")
    if os.path.exists(pcraster_dll_pathname):
        os.environ["PATH"] = pcraster_dll_pathname + os.pathsep + path_

try:
    from ._operations import *
    from . import _operators

    from .__about__ import (
        __version__, __author__, __uri__, __license__, __copyright__
    )
finally:
    if sys.platform == "win32":
        os.environ["PATH"] = path_
