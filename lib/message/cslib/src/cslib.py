# ------------------------------------------------------------------------
#   CSlib - Client/server library for code coupling
#   http://cslib.sandia.gov, Sandia National Laboratories
#   Steve Plimpton, sjplimp@sandia.gov
#
#   Copyright 2018 National Technology & Engineering Solutions of
#   Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
#   NTESS, the U.S. Government retains certain rights in this software.
#   This software is distributed under the modified Berkeley Software
#   Distribution (BSD) License.
#
#   See the README file in the top-level CSlib directory.
# -------------------------------------------------------------------------

# Python wrapper on CSlib library via ctypes

# ctypes and Numpy data types:
# 32-bit int = c_int = np.intc = np.int32
# 64-bit int = c_longlong = np.int64
# 32-bit floating point = c_float = np.float32
# 64-bit floating point = c_double = np.float = np.float64

import sys,traceback
from ctypes import *

# Numpy and mpi4py packages may not exist

try:
  import numpy as np
  numpyflag = 1
except:
  numpyflag = 0

try:
  from mpi4py import MPI
  mpi4pyflag = 1
except:
  mpi4pyflag = 0

# wrapper class

class CSlib:

  # instantiate CSlib thru its C-interface
  
  def __init__(self,csflag,mode,ptr,comm):

    # load libcslib.so
    
    try:
      if comm: self.lib = CDLL("libcsmpi.so",RTLD_GLOBAL)
      else: self.lib = CDLL("libcsnompi.so",RTLD_GLOBAL)
    except:
      etype,value,tb = sys.exc_info()
      traceback.print_exception(etype,value,tb)
      raise OSError,"Could not load CSlib dynamic library"

    # define ctypes API for each library method

    self.lib.cslib_open.argtypes = [c_int,c_char_p,c_void_p,c_void_p,
                                    POINTER(c_void_p)]
    self.lib.cslib_open.restype = None

    self.lib.cslib_close.argtypes = [c_void_p]
    self.lib.cslib_close.restype = None

    self.lib.cslib_send.argtypes = [c_void_p,c_int,c_int]
    self.lib.cslib_send.restype = None
    
    self.lib.cslib_pack_int.argtypes = [c_void_p,c_int,c_int]
    self.lib.cslib_pack_int.restype = None
    
    self.lib.cslib_pack_int64.argtypes = [c_void_p,c_int,c_longlong]
    self.lib.cslib_pack_int64.restype = None
    
    self.lib.cslib_pack_float.argtypes = [c_void_p,c_int,c_float]
    self.lib.cslib_pack_float.restype = None
    
    self.lib.cslib_pack_double.argtypes = [c_void_p,c_int,c_double]
    self.lib.cslib_pack_double.restype = None
    
    self.lib.cslib_pack_string.argtypes = [c_void_p,c_int,c_char_p]
    self.lib.cslib_pack_string.restype = None
    
    self.lib.cslib_pack.argtypes = [c_void_p,c_int,c_int,c_int,c_void_p]
    self.lib.cslib_pack.restype = None
    
    self.lib.cslib_pack_parallel.argtypes = [c_void_p,c_int,c_int,c_int,
                                             POINTER(c_int),c_int,c_void_p]
    self.lib.cslib_pack_parallel.restype = None
    
    self.lib.cslib_recv.argtypes = [c_void_p,POINTER(c_int),
                                    POINTER(POINTER(c_int)),
                                    POINTER(POINTER(c_int)),
                                    POINTER(POINTER(c_int))]
    self.lib.cslib_recv.restype = c_int
    
    self.lib.cslib_unpack_int.argtypes = [c_void_p,c_int]
    self.lib.cslib_unpack_int.restype = c_int

    self.lib.cslib_unpack_int64.argtypes = [c_void_p,c_int]
    self.lib.cslib_unpack_int64.restype = c_longlong

    self.lib.cslib_unpack_float.argtypes = [c_void_p,c_int]
    self.lib.cslib_unpack_float.restype = c_float

    self.lib.cslib_unpack_double.argtypes = [c_void_p,c_int]
    self.lib.cslib_unpack_double.restype = c_double

    self.lib.cslib_unpack_string.argtypes = [c_void_p,c_int]
    self.lib.cslib_unpack_string.restype = c_char_p

    # override return in unpack()
    self.lib.cslib_unpack.argtypes = [c_void_p,c_int]
    self.lib.cslib_unpack.restype = c_void_p

    self.lib.cslib_unpack_data.argtypes = [c_void_p,c_int,c_void_p]
    self.lib.cslib_unpack_data.restype = None

    # override last arg in unpack_parallel()
    self.lib.cslib_unpack_parallel.argtypes = [c_void_p,c_int,c_int,
                                               POINTER(c_int),c_int,c_void_p]
    self.lib.cslib_unpack_parallel.restype = None

    self.lib.cslib_extract.argtypes = [c_void_p,c_int]
    self.lib.cslib_extract.restype = c_int

    # create an instance of CSlib with or w/out MPI communicator

    self.cs = c_void_p()
    
    if not comm:
      self.lib.cslib_open(csflag,mode,ptr,None,byref(self.cs))
    elif not mpi4pyflag:
      print "Cannot pass MPI communicator to CSlib w/out mpi4py package"
      sys.exit()
    else:
      address = MPI._addressof(comm)
      comm_ptr = c_void_p(address)
      if mode == "mpi/one":
        address = MPI._addressof(ptr)
        ptrcopy = c_void_p(address)
      else: ptrcopy = ptr
      self.lib.cslib_open(csflag,mode,ptrcopy,comm_ptr,byref(self.cs))

  # destroy instance of CSlib
  
  def __del__(self):
    if self.cs: self.lib.cslib_close(self.cs)

  def close(self):
    self.lib.cslib_close(self.cs)
    self.lib = None

  # send a message
  
  def send(self,msgID,nfield):
    self.nfield = nfield
    self.lib.cslib_send(self.cs,msgID,nfield)

  # pack one field of message
  
  def pack_int(self,id,value):
    self.lib.cslib_pack_int(self.cs,id,value)

  def pack_int64(self,id,value):
    self.lib.cslib_pack_int64(self.cs,id,value)

  def pack_float(self,id,value):
    self.lib.cslib_pack_float(self.cs,id,value)

  def pack_double(self,id,value):
    self.lib.cslib_pack_double(self.cs,id,value)

  def pack_string(self,id,value):
    self.lib.cslib_pack_string(self.cs,id,value)

  def pack(self,id,ftype,flen,data):
    cdata = self.data_convert(ftype,flen,data)
    self.lib.cslib_pack(self.cs,id,ftype,flen,cdata)

  def pack_parallel(self,id,ftype,nlocal,ids,nper,data):
    cids = self.data_convert(1,nlocal,ids)
    cdata = self.data_convert(ftype,nper*nlocal,data)
    self.lib.cslib_pack_parallel(self.cs,id,ftype,nlocal,cids,nper,cdata)

  # convert input data to a ctypes vector to pass to CSlib
  
  def data_convert(self,ftype,flen,data):
       
    # tflag = type of data
    # tflag = 1 if data is list or tuple
    # tflag = 2 if data is Numpy array
    # tflag = 3 if data is ctypes vector
    # same usage of tflag as in unpack function
    
    txttype = str(type(data))
    if "numpy" in txttype: tflag = 2
    elif "c_" in txttype: tflag = 3
    else: tflag = 1
    
    # create ctypes vector out of data to pass to lib
    # cdata = ctypes vector to return
    # NOTE: error check on ftype and tflag everywhere, also flen
    
    if ftype == 1:
      if tflag == 1: cdata = (flen * c_int)(*data)
      elif tflag == 2: cdata = data.ctypes.data_as(POINTER(c_int))
      elif tflag == 3: cdata = data
    elif ftype == 2:
      if tflag == 1: cdata = (flen * c_longlong)(*data)
      elif tflag == 2: cdata = data.ctypes.data_as(POINTER(c_longlong))
      elif tflag == 3: cdata = data
    elif ftype == 3:
      if tflag == 1: cdata = (flen * c_float)(*data)
      elif tflag == 2: cdata = data.ctypes.data_as(POINTER(c_float))
      elif tflag == 3: cdata = data
    elif ftype == 4:
      if tflag == 1: cdata = (flen * c_double)(*data)
      elif tflag == 2: cdata = data.ctypes.data_as(POINTER(c_double))
      elif tflag == 3: cdata = data

    return cdata

  # receive a message
  
  def recv(self):
    self.lib.cslib_recv.restype = c_int
    nfield = c_int()
    fieldID = POINTER(c_int)()
    fieldtype = POINTER(c_int)()
    fieldlen = POINTER(c_int)()
    msgID = self.lib.cslib_recv(self.cs,byref(nfield),
                                byref(fieldID),byref(fieldtype),byref(fieldlen))

    # copy returned C args to native Python int and lists
    # store them in class so unpack() methods can access the info
    
    self.nfield = nfield = nfield.value
    self.fieldID = fieldID[:nfield]
    self.fieldtype = fieldtype[:nfield]
    self.fieldlen = fieldlen[:nfield]
    
    return msgID,self.nfield,self.fieldID,self.fieldtype,self.fieldlen

  # unpack one field of message
  # tflag = type of data to return
  # 3 = ctypes vector is default, since no conversion required
  
  def unpack_int(self,id):
    return self.lib.cslib_unpack_int(self.cs,id)

  def unpack_int64(self,id):
    return self.lib.cslib_unpack_int64(self.cs,id)

  def unpack_float(self,id):
    return self.lib.cslib_unpack_float(self.cs,id)

  def unpack_double(self,id):
    return self.lib.cslib_unpack_double(self.cs,id)

  def unpack_string(self,id):
    return self.lib.cslib_unpack_string(self.cs,id)

  def unpack(self,id,tflag=3):
    index = self.fieldID.index(id)

    # reset data type of return so can morph by tflag
    # cannot do this for the generic c_void_p returned by CSlib
    
    if self.fieldtype[index] == 1:
      self.lib.cslib_unpack.restype = POINTER(c_int)
    elif self.fieldtype[index] == 2:
      self.lib.cslib_unpack.restype = POINTER(c_longlong)
    elif self.fieldtype[index] == 3:
      self.lib.cslib_unpack.restype = POINTER(c_float)
    elif self.fieldtype[index] == 4:
      self.lib.cslib_unpack.restype = POINTER(c_double)
    #elif self.fieldtype[index] == 5:
    #  self.lib.cslib_unpack.restype = POINTER(c_char)

    cdata = self.lib.cslib_unpack(self.cs,id)

    # tflag = user-requested type of data to return
    # tflag = 1 to return data as list
    # tflag = 2 to return data as Numpy array
    # tflag = 3 to return data as ctypes vector
    # same usage of tflag as in pack functions
    # tflag = 2,3 should NOT perform a data copy
    
    if tflag == 1:
      data = cdata[:self.fieldlen[index]]
    elif tflag == 2:
      if numpyflag == 0:
        print "Cannot return Numpy array w/out numpy package"
        sys.exit()
      data = np.ctypeslib.as_array(cdata,shape=(self.fieldlen[index],))
    elif tflag == 3:
      data = cdata
      
    return data

  # handle data array like pack() or unpack_parallel() ??
  
  def unpack_data(self,id,tflag=3):
    index = self.fieldID.index(id)

  # unpack one field of message in parallel
  # tflag = type of data to return
  # 3 = ctypes vector is default, since no conversion required
  # NOTE: allow direct use of user array (e.g. Numpy), if user provides data arg?
  #       as opposed to creating this cdata
  #       does that make any performance difference ?
  #       e.g. should we allow CSlib to populate an existing Numpy array's memory
  
  def unpack_parallel(self,id,nlocal,ids,nper,tflag=3):
    cids = self.data_convert(1,nlocal,ids)

    # allocate memory for the returned data
    # pass cdata ptr to the memory to CSlib unpack_parallel()
    # this resets data type of last unpack_parallel() arg
    
    index = self.fieldID.index(id)
    if self.fieldtype[index] == 1: cdata = (nper*nlocal * c_int)()
    elif self.fieldtype[index] == 2: cdata = (nlocal*nper * c_longlong)()
    elif self.fieldtype[index] == 3: cdata = (nlocal*nper * c_float)()
    elif self.fieldtype[index] == 4: cdata = (nlocal*nper * c_double)()
    #elif self.fieldtype[index] == 5: cdata = (nlocal*nper * c_char)()

    self.lib.cslib_unpack_parallel(self.cs,id,nlocal,cids,nper,cdata)

    # tflag = user-requested type of data to return
    # tflag = 1 to return data as list
    # tflag = 2 to return data as Numpy array
    # tflag = 3 to return data as ctypes vector
    # same usage of tflag as in pack functions
    
    if tflag == 1:
      data = cdata[:nper*nlocal]
    elif tflag == 2:
      if numpyflag == 0:
        print "Cannot return Numpy array w/out numpy package"
        sys.exit()
      # NOTE: next line gives ctypes warning for fieldtype = 2 = 64-bit int
      # not sure why, reported as bug between ctypes and Numpy here:
      # https://stackoverflow.com/questions/4964101/pep-3118-
      #         warning-when-using-ctypes-array-as-numpy-array
      # but why not same warning when just using unpack() ??
      # in Python these lines give same warning:
      # >>> import ctypes,numpy
      # >>> a = (10 * ctypes.c_longlong)()
      # >>> b = numpy.ctypeslib.as_array(a)
      data = np.ctypeslib.as_array(cdata,shape=(nlocal*nper,))
    elif tflag == 3:
      data = cdata
      
    return data

  # extract a library value
  
  def extract(self,flag):
   return self.lib.cslib_extract(self.cs,flag)
