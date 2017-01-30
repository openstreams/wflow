import os
import sys


# On Windows prepend the path to our dlls to the PATH environment variable.
# Otherwise our dlls won't be found when our Python extensions are loaded
# by Python.
if sys.platform == "win32":
    path_ = os.environ["PATH"]
    pcraster_installation_root = os.path.abspath(os.path.join(
        os.path.dirname(__file__), "..", ".."))
    pcraster_dll_pathname = os.path.join(pcraster_installation_root, "lib")
    if os.path.exists(pcraster_dll_pathname):
        os.environ["PATH"] = pcraster_dll_pathname + os.pathsep + path_

try:
    from operations import *
    import operators
    from _pcraster import *
    from _pcraster_modflow import *
    from aguila import *
    from numpy_operations import *
finally:
    if sys.platform == "win32":
        os.environ["PATH"] = path_
