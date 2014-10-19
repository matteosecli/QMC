# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.8
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.



from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_lib', [dirname(__file__)])
        except ImportError:
            import _lib
            return _lib
        if fp is not None:
            try:
                _mod = imp.load_module('_lib', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _lib = swig_import_helper()
    del swig_import_helper
else:
    import _lib
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


ZERO = _lib.ZERO
class TID(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, TID, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, TID, name)
    __repr__ = _swig_repr
    __swig_setmethods__["tick"] = _lib.TID_tick_set
    __swig_getmethods__["tick"] = _lib.TID_tick_get
    if _newclass:tick = _swig_property(_lib.TID_tick_get, _lib.TID_tick_set)
    __swig_setmethods__["sec"] = _lib.TID_sec_set
    __swig_getmethods__["sec"] = _lib.TID_sec_get
    if _newclass:sec = _swig_property(_lib.TID_sec_get, _lib.TID_sec_set)
    __swig_setmethods__["min"] = _lib.TID_min_set
    __swig_getmethods__["min"] = _lib.TID_min_get
    if _newclass:min = _swig_property(_lib.TID_min_get, _lib.TID_min_set)
    __swig_setmethods__["hour"] = _lib.TID_hour_set
    __swig_getmethods__["hour"] = _lib.TID_hour_get
    if _newclass:hour = _swig_property(_lib.TID_hour_get, _lib.TID_hour_set)
    def __init__(self): 
        this = _lib.new_TID()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _lib.delete_TID
    __del__ = lambda self : None;
TID_swigregister = _lib.TID_swigregister
TID_swigregister(TID)
cvar = _lib.cvar


def time_step(*args):
  return _lib.time_step(*args)
time_step = _lib.time_step

def matrix(*args):
  return _lib.matrix(*args)
matrix = _lib.matrix

def free_matrix(*args):
  return _lib.free_matrix(*args)
free_matrix = _lib.free_matrix

def rk4(*args):
  return _lib.rk4(*args)
rk4 = _lib.rk4

def ludcmp(*args):
  return _lib.ludcmp(*args)
ludcmp = _lib.ludcmp

def lubksb(*args):
  return _lib.lubksb(*args)
lubksb = _lib.lubksb

def tqli(*args):
  return _lib.tqli(*args)
tqli = _lib.tqli

def tred2(*args):
  return _lib.tred2(*args)
tred2 = _lib.tred2

def pythag(*args):
  return _lib.pythag(*args)
pythag = _lib.pythag

def gauleg(*args):
  return _lib.gauleg(*args)
gauleg = _lib.gauleg

def jacobi(*args):
  return _lib.jacobi(*args)
jacobi = _lib.jacobi

def rectangle_rule(*args):
  return _lib.rectangle_rule(*args)
rectangle_rule = _lib.rectangle_rule

def trapezoidal_rule(*args):
  return _lib.trapezoidal_rule(*args)
trapezoidal_rule = _lib.trapezoidal_rule

def spline(*args):
  return _lib.spline(*args)
spline = _lib.spline

def splint(*args):
  return _lib.splint(*args)
splint = _lib.splint

def polint(*args):
  return _lib.polint(*args)
polint = _lib.polint

def rtbis(*args):
  return _lib.rtbis(*args)
rtbis = _lib.rtbis

def rtsec(*args):
  return _lib.rtsec(*args)
rtsec = _lib.rtsec

def rtnewt(*args):
  return _lib.rtnewt(*args)
rtnewt = _lib.rtnewt

def zbrent(*args):
  return _lib.zbrent(*args)
zbrent = _lib.zbrent

def ran0(*args):
  return _lib.ran0(*args)
ran0 = _lib.ran0

def ran1(*args):
  return _lib.ran1(*args)
ran1 = _lib.ran1

def ran2(*args):
  return _lib.ran2(*args)
ran2 = _lib.ran2

def ran3(*args):
  return _lib.ran3(*args)
ran3 = _lib.ran3
# This file is compatible with both classic and new-style classes.


