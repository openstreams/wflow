import sys
import pcraster
from . import operations, _pcraster

def pcrAnd(self, field):
  return operations.pcrand(self, field)

def pcrRAnd(self, number):
  return operations.pcrand(number, self)

def pcrOr(self, field):
  return operations.pcror(self, field)

def pcrROr(self, number):
  return operations.pcror(number, self)

def pcrXOr(self, field):
  return operations.pcrxor(self, field)

def pcrNot(self):
  return operations.pcrnot(self)

def pcrNE(self, field):
  return operations.pcrne(self, field)

def pcrEQ(self, field):
  return operations.pcreq(self, field)

def pcrGT(self, field):
  return operations.pcrgt(self, field)

def pcrGE(self, field):
  return operations.pcrge(self, field)

def pcrLT(self, field):
  return operations.pcrlt(self, field)

def pcrLE(self, field):
  return operations.pcrle(self, field)

def pcrMul(self, field):
  return operations.pcrmul(self, field)

def pcrRMul(self, number):
  return operations.pcrmul(number, self)

def pcrDiv(self, field):
  return operations.pcrfdiv(self, field)

def pcrRDiv(self, number):
  return operations.pcrfdiv(number, self)

def pcrFloorDiv(self, field):
  return operations.pcridiv(self, field)

def pcrRFloorDiv(self, number):
  return operations.pcridiv(number, self)

def pcrPow(self, field):
  return operations.pcrpow(self, field)

def pcrRPow(self, number):
  return operations.pcrpow(number, self)

def pcrMod(self, field):
  return operations.pcrmod(self, field)

def pcrRMod(self, number):
  return operations.pcrmod(number, self)

def pcrAdd(self, field):
  return operations.pcrbadd(self, field)

def pcrRAdd(self, number):
  return operations.pcrbadd(number, self)

def pcrSub(self, field):
  return operations.pcrbmin(self, field)

def pcrRSub(self, number):
  return operations.pcrbmin(number, self)

def pcrNeg(self):
  return operations.pcrumin(self)

def pcrPos(self):
  return operations.pcruadd(self)

def _bool(self):
  if self.isSpatial():
    raise Exception(
         "The truth value for PCRaster spatial data types is ambiguous. "
         "See the section Boolean operations in the PCRaster Python manual.")

  result = None
  value, isValid = pcraster.cellvalue(self, 0)

  if isValid:
    result = bool(value)

  return result

def _int(self):
  if self.isSpatial():
    raise Exception(
         "The integer value for PCRaster spatial data types is ambiguous.")

  result = None
  value, isValid = pcraster.cellvalue(self, 0)

  if isValid:
    result = int(value)

  return result

def _float(self):
  if self.isSpatial():
    raise Exception(
         "The floating point value for PCRaster spatial data types is ambiguous.")

  result = None
  value, isValid = pcraster.cellvalue(self, 0)

  if isValid:
    result = float(value)

  return result


# Some syntactic sugar to allow the use of operators on Field objects.
_pcraster.Field.__and__      = pcrAnd      # &
_pcraster.Field.__rand__     = pcrRAnd
_pcraster.Field.__or__       = pcrOr       # |
_pcraster.Field.__ror__      = pcrROr
_pcraster.Field.__xor__      = pcrXOr      # ^
_pcraster.Field.__invert__   = pcrNot      # ~

_pcraster.Field.__ne__       = pcrNE
_pcraster.Field.__eq__       = pcrEQ
_pcraster.Field.__gt__       = pcrGT
_pcraster.Field.__ge__       = pcrGE
_pcraster.Field.__lt__       = pcrLT
_pcraster.Field.__le__       = pcrLE

_pcraster.Field.__mul__      = pcrMul
_pcraster.Field.__rmul__     = pcrRMul
_pcraster.Field.__floordiv__ = pcrFloorDiv
_pcraster.Field.__rfloordiv__ = pcrRFloorDiv
_pcraster.Field.__pow__      = pcrPow
_pcraster.Field.__rpow__     = pcrRPow
_pcraster.Field.__mod__      = pcrMod
_pcraster.Field.__rmod__     = pcrRMod
_pcraster.Field.__add__      = pcrAdd
_pcraster.Field.__radd__     = pcrRAdd
_pcraster.Field.__sub__      = pcrSub
_pcraster.Field.__rsub__     = pcrRSub

_pcraster.Field.__neg__      = pcrNeg
_pcraster.Field.__pos__      = pcrPos

_pcraster.Field.__bool__ = _bool
_pcraster.Field.__int__ = _int
_pcraster.Field.__float__ = _float


if sys.version_info[0] < 3:
  _pcraster.Field.__div__  = pcrDiv
  _pcraster.Field.__rdiv__     = pcrRDiv
  _pcraster.Field.__nonzero__  = _bool # pcrNonzero
else:
  _pcraster.Field.__truediv__  = pcrDiv
  _pcraster.Field.__rtruediv__ = pcrRDiv
  # Instead of __nonzero__, Python 3 calls __bool__


