# This file was generated during the PCRaster build.
# Do not edit!
try:
    import pcraster._pcraster as _pcraster
except:
    # we rely on a correct PCRaster setup
    raise ImportError("Failed to import the module '_pcraster'")

from . import _operations as pmcop


def pcrmcAnd(self, field):
    return pmcop._and(self, field)


def pcrmcRAnd(self, number):
    return pmcop._and(number, self)


def pcrmcOr(self, field):
    return pmcop._or(self, field)


def pcrmcROr(self, number):
    return pmcop._or(number, self)


def pcrmcXOr(self, field):
    return pmcop._xor(self, field)


def pcrmcNot(self):
    return pmcop._not(self)


def pcrmcEQ(self, field):
    return pmcop.equal(self, field)


def pcrmcNE(self, field):
    return pmcop.unequal(self, field)


def pcrmcGT(self, field):
    return pmcop.greater(self, field)


def pcrmcGE(self, field):
    return pmcop.greater_equal(self, field)


def pcrmcLT(self, field):
    return pmcop.less(self, field)


def pcrmcLE(self, field):
    return pmcop.less_equal(self, field)


def pcrmcAdd(self, field):
    return pmcop.add(self, field)


def pcrmcRAdd(self, number):
    return pmcop.add(number, self)


def pcrmcSub(self, field):
    return pmcop.sub(self, field)


def pcrmcRSub(self, number):
    return pmcop.sub(number, self)


def pcrmcMul(self, field):
    return pmcop.mul(self, field)


def pcrmcRMul(self, number):
    return pmcop.mul(number, self)


def pcrmcDiv(self, field):
    return pmcop.div(self, field)


def pcrmcRDiv(self, number):
    return pmcop.div(number, self)


def pcrmcPow(self, variable):
    return pmcop.power(self, variable)


def pcrmcRPow(self, number):
    return pmcop.power(number, self)


# Some syntactic sugar to still allow the use of operators on Field objects.
# These operators should overwrite the ones from the PCRaster module.
_pcraster.Field.__and__ = pcrmcAnd      # &
_pcraster.Field.__rand__ = pcrmcRAnd

_pcraster.Field.__or__ = pcrmcOr        # |
_pcraster.Field.__ror__ = pcrmcROr

_pcraster.Field.__invert__ = pcrmcNot   # ~
_pcraster.Field.__xor__ = pcrmcXOr      # ^

_pcraster.Field.__ne__ = pcrmcNE
_pcraster.Field.__eq__ = pcrmcEQ

_pcraster.Field.__gt__ = pcrmcGT
_pcraster.Field.__ge__ = pcrmcGE

_pcraster.Field.__lt__ = pcrmcLT
_pcraster.Field.__le__ = pcrmcLE

_pcraster.Field.__add__ = pcrmcAdd
_pcraster.Field.__radd__ = pcrmcRAdd

_pcraster.Field.__sub__ = pcrmcSub
_pcraster.Field.__rsub__ = pcrmcRSub

_pcraster.Field.__mul__ = pcrmcMul
_pcraster.Field.__rmul__ = pcrmcRMul

_pcraster.Field.__div__ = pcrmcDiv
_pcraster.Field.__rdiv__ = pcrmcRDiv

_pcraster.Field.__pow__ = pcrmcPow
_pcraster.Field.__rpow__ = pcrmcRPow
