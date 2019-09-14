# This file was generated during the PCRaster build.
# Do not edit!
import sys
import math
import pcraster.operations
import pcraster as _pcraster
from . import _pcraster_multicore as _pmc

if sys.version_info.major == 2:
    import __builtin__ as builtins
else:
    import builtins


def set_nr_worker_threads(arg1):
    _pmc.set_nr_worker_threads(arg1)


def nr_worker_threads():
    return _pmc.nr_worker_threads()


# just get a spatial or nonspatial
# assign a valuescale if necessary
# correct SV of an operation should be checked in the C++ code
def _get_field(arg):
    if isinstance(arg, _pcraster.Field):
        return(arg)
    elif isinstance(arg, str):
        return(_pcraster.readmap(arg))
    elif isinstance(arg, bool):
        return(_pmc._newNonSpatialBoolean(arg))
    elif isinstance(arg, int):  # Not in python-3.5 or isinstance(arg, long):
        return(_pmc._newNonSpatialNominal(arg))
    elif isinstance(arg, float):
        return(_pmc._newNonSpatialScalar(arg))
    else:
        msg = "conversion of argument with type '{}' to PCRaster not possible".format(type(arg).__name__)
        raise RuntimeError(msg)


def _all_python_pod(*arguments):
    # Helper function used to determine wheter we can call
    # builtin functions or the overloaded PCRaster ones
    for argument in arguments:
        if isinstance(argument, _pcraster.Field):
            return False
    return True



def _sort(*args):
    # Places PCRaster spatial fields at the start of the list
    # min/max need spatials as first argument
    values = []
    for item in args:
        if isinstance(item, _pcraster.Field):
            if item.isSpatial():
                values.insert(0, item)
            else:
                values.append(item)
        else:
            values.append(item)
    return values

def _oversized_kernel(radius):
    """ Check whether the desired kernel exceeds the raster extent """
    # ToDo: this error handling should be handled by Fern...
    # determine exceedance based on array/cell units
    maximum_window = min(_pcraster.clone().nrRows(), _pcraster.clone().nrCols())
    maximum_radius = int(math.floor((maximum_window - 1) / 2.0))

    if radius > maximum_radius:
        return True
    else:
        return False

def _radius_in_cells(window_length):
    """ returns kernel radius suitable for Fern, 0 otherwise """
    kernel_radius = 0
    kernel_length_cells = 0

    # Determine kernel width in number of cells
    if _pmc._unittrue():
        cellSize = _pcraster.cellvalue(_pcraster.celllength(), 1)[0]

        if float(window_length) % float(cellSize) == 0:
            # Window length 'matches' a number of cells
            kernel_length_cells = int(window_length / cellSize)
    else:
        kernel_length_cells = int(window_length)

    # Radius must be an odd number in Fern...
    if (kernel_length_cells % 2) != 0 and kernel_length_cells > 2:
        kernel_radius = int(math.floor(kernel_length_cells / 2.0))

    # Kernel must not exceed raster extent
    if _oversized_kernel(kernel_radius):
      return 0
    else:
      return kernel_radius

def windowtotal(arg1, arg2):
    # Field as window length argument: use traditional PCRaster implementation
    if isinstance(arg2, _pcraster.Field):
      return pcraster.operations.windowtotal(arg1, arg2)

    radius = _radius_in_cells(arg2)

    if radius > 0:
      try:
        arg1 = _get_field(arg1)
        return _pmc.windowtotal(arg1, int(radius))
      except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore windowtotal: {}\n".format(str(exception)))

    # Otherwise,  fall back to the traditional PCRaster implementation
    return pcraster.operations.windowtotal(arg1, arg2)

def windowaverage(arg1, arg2):
    # Field as window length argument: use traditional PCRaster implementation
    if isinstance(arg2, _pcraster.Field):
      return pcraster.operations.windowaverage(arg1, arg2)

    radius = _radius_in_cells(arg2)

    if radius > 0:
      # Multicore window must be of size odd number of cells
      try:
        arg1 = _get_field(arg1)
        return _pmc.windowaverage(arg1, int(radius))
      except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore windowaverage: {}".format(str(exception)))

    # Otherwise,  fall back to the traditional PCRaster implementation
    return pcraster.operations.windowaverage(arg1, arg2)

def cover(*args):
    try:
        fields = list(args)
        for i, field in enumerate(fields):
            fields[i] = _get_field(field)
        if fields[0].dataType() == _pcraster.VALUESCALE.Ldd:
            return _pcraster.operations.cover(*args)
        if fields[0].dataType() == _pcraster.VALUESCALE.Directional:
            return _pcraster.operations.cover(*args)
        return _pmc.cover(fields)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore cover: {}".format(str(exception)))

def max(*args):
    # Use the builtin Python max in case of non-PCRaster data types
    if _all_python_pod(*args):
        return builtins.max(args)
    try:
        fields = _sort(*args)
        all_nonspatial = True
        for i, field in enumerate(fields):
            fields[i] = _get_field(field)
            if fields[i].isSpatial():
                all_nonspatial = False
        if all_nonspatial == False:
            return _pmc.maximum(fields)
        else:
            return _pcraster.operations.max(*args)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore max: {}".format(str(exception)))


def min(*args):
    # Use the builtin Python max in case of non-PCRaster data types
    if _all_python_pod(*args):
        return builtins.min(args)
    try:
        fields = _sort(*args)
        all_nonspatial = True
        for i, field in enumerate(fields):
            fields[i] = _get_field(field)
            if fields[i].isSpatial():
                all_nonspatial = False
        if all_nonspatial == False:
            return _pmc.minimum(fields)
        else:
            return _pcraster.operations.min(*args)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore min: {}".format(str(exception)))


def ifthen(arg1, arg2):

    # to mimic the PCRaster model engine behaviour
    if not isinstance(arg2, _pcraster.Field):
        raise RuntimeError("Use a conversion function to pick a data type")

    try:
        arg1 = _get_field(arg1)
        arg2 = _get_field(arg2)
        return _pmc.ifthen(arg1, arg2)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore ifthen: {}".format(str(exception)))


def ifthenelse(arg1, arg2, arg3):

    # to mimic the PCRaster model engine behaviour
    if not ( isinstance(arg2, _pcraster.Field) or isinstance(arg3, _pcraster.Field) ):
        raise RuntimeError("Use a conversion function to pick a data type")
    try:
        arg1 = _get_field(arg1)
        arg2 = _get_field(arg2)
        arg3 = _get_field(arg3)
        return _pmc.ifthenelse(arg1, arg2, arg3)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore ifthenelse: {}".format(str(exception)))




def boolean(arg1):

    # Fallback to model_engine
    if isinstance(arg1, _pcraster.Field):
        if arg1.dataType() ==  _pcraster.VALUESCALE.Directional:
            return pcraster.operations.boolean(arg1)

    try:
        arg1 = _get_field(arg1)
        return _pmc.boolean(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore boolean: {}".format(str(exception)))


def nominal(arg1):

    # Fallback to model_engine
    if isinstance(arg1, _pcraster.Field):
        if arg1.dataType() ==  _pcraster.VALUESCALE.Directional:
            return pcraster.operations.nominal(arg1)

    try:
        arg1 = _get_field(arg1)
        return _pmc.nominal(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore nominal: {}".format(str(exception)))


def ordinal(arg1):

    # Fallback to model_engine
    if isinstance(arg1, _pcraster.Field):
        if arg1.dataType() ==  _pcraster.VALUESCALE.Directional:
            return pcraster.operations.ordinal(arg1)

    try:
        arg1 = _get_field(arg1)
        return _pmc.ordinal(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore ordinal: {}".format(str(exception)))


def scalar(arg1):

    # Fallback to model_engine
    if isinstance(arg1, _pcraster.Field):
        if arg1.dataType() ==  _pcraster.VALUESCALE.Directional:
            return pcraster.operations.scalar(arg1)
    try:
        arg1 = _get_field(arg1)
        return _pmc.scalar(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore scalar: {}".format(str(exception)))




# Overload operator function names from pcraster module

def pcrand(arg1, arg2):
    return _and(arg1, arg2)

def pcror(arg1, arg2):
    return _or(arg1, arg2)

def pcrxor(arg1, arg2):
    return _xor(arg1, arg2)

def pcrnot(arg1):
    return _not(arg1)

def pcrne(arg1, arg2):
    return unequal(arg1, arg2)

def pcreq(arg1, arg2):
    return equal(arg1, arg2)

def pcrgt(arg1, arg2):
    return greater(arg1, arg2)

def pcrge(arg1, arg2):
    return greater_equal(arg1, arg2)

def pcrlt(arg1, arg2):
    return less(arg1, arg2)

def pcrle(arg1, arg2):
    return less_equal(arg1, arg2)


# these ones below here can be generated?!

def slope(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.slope(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore slope: {}".format(str(exception)))


def equal(arg1, arg2):
    try:
        arg1 = _get_field(arg1)
        arg2 = _get_field(arg2)
        return _pmc.equal(arg1, arg2)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore equal: {}".format(str(exception)))


def unequal(arg1, arg2):
    try:
        arg1 = _get_field(arg1)
        arg2 = _get_field(arg2)
        return _pmc.unequal(arg1, arg2)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore unequal: {}".format(str(exception)))


def greater(arg1, arg2):
    try:
        arg1 = _get_field(arg1)
        arg2 = _get_field(arg2)
        return _pmc.greater(arg1, arg2)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore greater: {}".format(str(exception)))


def greater_equal(arg1, arg2):
    try:
        arg1 = _get_field(arg1)
        arg2 = _get_field(arg2)
        return _pmc.greater_equal(arg1, arg2)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore greater_equal: {}".format(str(exception)))


def less(arg1, arg2):
    try:
        arg1 = _get_field(arg1)
        arg2 = _get_field(arg2)
        return _pmc.less(arg1, arg2)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore less: {}".format(str(exception)))


def less_equal(arg1, arg2):
    try:
        arg1 = _get_field(arg1)
        arg2 = _get_field(arg2)
        return _pmc.less_equal(arg1, arg2)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore less_equal: {}".format(str(exception)))


def _and(arg1, arg2):
    try:
        arg1 = _get_field(arg1)
        arg2 = _get_field(arg2)
        return _pmc._and(arg1, arg2)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore and: {}".format(str(exception)))


def _or(arg1, arg2):
    try:
        arg1 = _get_field(arg1)
        arg2 = _get_field(arg2)
        return _pmc._or(arg1, arg2)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore or: {}".format(str(exception)))


def _xor(arg1, arg2):
    try:
        arg1 = _get_field(arg1)
        arg2 = _get_field(arg2)
        return _pmc._xor(arg1, arg2)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore xor: {}".format(str(exception)))


def _not(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc._not(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore not: {}".format(str(exception)))


def defined(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.defined(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore defined: {}".format(str(exception)))


def add(arg1, arg2):
    try:
        arg1 = _get_field(arg1)
        arg2 = _get_field(arg2)
        return _pmc.add(arg1, arg2)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore add: {}".format(str(exception)))


def sub(arg1, arg2):
    try:
        arg1 = _get_field(arg1)
        arg2 = _get_field(arg2)
        return _pmc.sub(arg1, arg2)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore sub: {}".format(str(exception)))


def mul(arg1, arg2):
    try:
        arg1 = _get_field(arg1)
        arg2 = _get_field(arg2)
        return _pmc.mul(arg1, arg2)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore mul: {}".format(str(exception)))


def div(arg1, arg2):
    try:
        arg1 = _get_field(arg1)
        arg2 = _get_field(arg2)
        return _pmc.div(arg1, arg2)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore div: {}".format(str(exception)))


def sqrt(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.sqrt(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore sqrt: {}".format(str(exception)))


def power(arg1, arg2):
    try:
        arg1 = _get_field(arg1)
        arg2 = _get_field(arg2)
        return _pmc.power(arg1, arg2)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore power: {}".format(str(exception)))


def abs(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.abs(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore abs: {}".format(str(exception)))


def cos(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.cos(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore cos: {}".format(str(exception)))


def sin(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.sin(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore sin: {}".format(str(exception)))


def tan(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.tan(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore tan: {}".format(str(exception)))


def sqr(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.sqr(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore sqr: {}".format(str(exception)))



def mapmaximum(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.mapmaximum(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore mapmaximum: {}".format(str(exception)))


def mapminimum(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.mapminimum(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore mapminimum: {}".format(str(exception)))




def acos(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.acos(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore acos: {}".format(str(exception)))


def asin(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.asin(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore asin: {}".format(str(exception)))


def atan(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.atan(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore atan: {}".format(str(exception)))


def fac(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.fac(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore fac: {}".format(str(exception)))


def window4total(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.window4total(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore window4total: {}".format(str(exception)))


def ln(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.ln(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore ln: {}".format(str(exception)))


def log10(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.log10(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore log10: {}".format(str(exception)))


def rounddown(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.rounddown(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore rounddown: {}".format(str(exception)))


def roundup(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.roundup(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore roundup: {}".format(str(exception)))


def roundoff(arg1):
    try:
        arg1 = _get_field(arg1)
        return _pmc.roundoff(arg1)
    except RuntimeError as exception:
        raise RuntimeError("pcraster.multicore roundoff: {}".format(str(exception)))



# to be generated...

abs.__doc__ = _pmc.abs.__doc__
acos.__doc__ = _pmc.acos.__doc__
add.__doc__ = _pmc.add.__doc__
asin.__doc__ = _pmc.asin.__doc__
atan.__doc__ = _pmc.atan.__doc__
boolean.__doc__ = _pmc.boolean.__doc__
cos.__doc__ = _pmc.cos.__doc__
cover.__doc__ = _pmc.cover.__doc__
defined.__doc__ = _pmc.defined.__doc__
div.__doc__ = _pmc.div.__doc__
equal.__doc__ = _pmc.equal.__doc__
fac.__doc__ = _pmc.fac.__doc__
greater.__doc__ = _pmc.greater.__doc__
greater_equal.__doc__ = _pmc.greater_equal.__doc__
ifthen.__doc__ = _pmc.ifthen.__doc__
ifthenelse.__doc__ = _pmc.ifthenelse.__doc__
less.__doc__ = _pmc.less.__doc__
less_equal.__doc__ = _pmc.less_equal.__doc__
ln.__doc__ = _pmc.ln.__doc__
log10.__doc__ = _pmc.log10.__doc__
max.__doc__ = _pmc.maximum.__doc__
min.__doc__ = _pmc.minimum.__doc__
mul.__doc__ = _pmc.mul.__doc__
nominal.__doc__ = _pmc.nominal.__doc__
nr_worker_threads.__doc__ = _pmc.nr_worker_threads.__doc__
ordinal.__doc__ = _pmc.ordinal.__doc__
power.__doc__ = _pmc.power.__doc__
rounddown.__doc__ = _pmc.rounddown.__doc__
roundoff.__doc__ = _pmc.roundoff.__doc__
roundup.__doc__ = _pmc.roundup.__doc__
scalar.__doc__ = _pmc.scalar.__doc__
set_nr_worker_threads.__doc__ = _pmc.set_nr_worker_threads.__doc__
sin.__doc__ = _pmc.sin.__doc__
sqr.__doc__ = _pmc.sqr.__doc__
sqrt.__doc__ = _pmc.sqrt.__doc__
sub.__doc__ = _pmc.sub.__doc__
tan.__doc__ = _pmc.tan.__doc__
unequal.__doc__ = _pmc.unequal.__doc__
# Focal operations
slope.__doc__ = _pmc.slope.__doc__
window4total.__doc__ = _pmc.window4total.__doc__
windowtotal.__doc__ = _pmc.windowtotal.__doc__
windowaverage.__doc__ = _pmc.windowaverage.__doc__
#  Operations on full map extent
mapmaximum.__doc__ = _pmc.mapmaximum.__doc__
mapminimum.__doc__ = _pmc.mapminimum.__doc__
