import os
import sys
import platform
from distutils.version import StrictVersion
import warnings


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
    try:
        from ._pcraster import *
        from ._pcraster_modflow import *
    except ImportError as error:
        # Test matching Python versions
        major, minor, patchlevel = platform.python_version_tuple()
        major = int(major)

        if major != 3:
            msg = "The 'pcraster' module was built for Python 3.6, the version used is {}.{}".format(
                major, minor)
            raise ImportError(msg)

        # Test matching bitness
        nr_bits = 64 if sys.maxsize > 2**32 else 32

        if nr_bits != 64:
            msg = "The 'pcraster' module was built for 64-bit, the Python version used is {}-bit".format(
                nr_bits)
            raise ImportError(msg)

        # Too old Python 2 versions might not work
        if major == 2:
            if StrictVersion("{}.{}".format(major, minor)) < StrictVersion("2.7"):
                msg = "The 'pcraster' module was built for Python 3.6, the version used is {}.{}".format(
                    major, minor)
                warnings.warn(msg, RuntimeWarning)

        msg = ""

        # VS2015 and runtime components related issues
        if sys.platform == "win32":
            if StrictVersion("{}.{}".format(major, minor)) != StrictVersion("3.6"):
                msg += "The 'pcraster' module was built for Python 3.6, the version used is {}.{}".format(
                    major, minor)
                raise ImportError(msg)

            msg += "The 'Microsoft Visual C++ 2015 Redistributable Update 3' is required to run PCRaster, available at:\nhttps://www.microsoft.com/en-us/download/details.aspx?id=53840\n"

        # Something else went wrong...
        msg += "{}\n".format(error)
        msg += "In case you cannot solve the issue please consult the PCRaster mailing list at:\nhttps://lists.geo.uu.nl/mailman/listinfo/pcraster-info\n"
        raise ImportError(msg)

    from .operations import *
    from . import operators
    from .aguila import *
    from .numpy_operations import *

    try:
        _var = "PCRASTER_NR_WORKER_THREADS"
        _workers = int(os.environ[_var])
        from .multicore import *
        set_nr_worker_threads(_workers)
    except ValueError as err:
        msg = "Could not obtain the number of worker threads. The environment variable {} is set to '{}'; use an integer value instead.".format(
            _var, os.environ[_var])
        raise RuntimeError(msg)
    except KeyError:
        pass  # No multicore module use intended
    except ImportError:
        pass  # Multicore module was not installed

finally:
    if sys.platform == "win32":
        os.environ["PATH"] = path_
