# try:
#     import numpy
# except:
#     raise ImportError, "Import of module 'numpy' failed."


# _minimal_numpy_version = "1.2.1"
# if numpy.__version__[:5] < _minimal_numpy_version:
#     raise ImportError, "NumPy needs to be of version %s or higher, %s is installed" % (_minimal_numpy_version, numpy.__version__[:5])


import numpy
import _pcraster


def pcr2numpy(
        map,
        mv):
    """
    Convert entities from PCRaster to NumPy.

    map -- Map you want to convert.
    mv -- Value to use in the result array cells as a missing value.

    Returns an array.
    """
    return _pcraster.pcr2numpy(_pcraster.clone(), map, mv)


def numpy2pcr(
        dataType,
        array,
        mv):
    """
    Convert entities from NumPy to PCRaster.

    dataType -- Either Boolean, Nominal, Ordinal, Scalar, Directional or
                Ldd.
    array -- Array you want to convert.
    mv -- Value that identifies a missing value in the array.

    Returns a map.
    """
    return _pcraster.numpy2pcr(_pcraster.clone(), dataType, array, mv)


def pcr_as_numpy(
        map):
    """
    Reference PCRaster maps from NumPy arrays.

    map -- Map to reference.

    Returns an array.
    """
    return _pcraster.pcr_as_numpy(_pcraster.clone(), map)
