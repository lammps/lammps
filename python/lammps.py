# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   http://lammps.sandia.gov, Sandia National Laboratories
#   Steve Plimpton, sjplimp@sandia.gov
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under 
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
# -------------------------------------------------------------------------

# Python wrapper on LAMMPS library via ctypes

import types
from ctypes import *

LMPINT = 0
LMPDOUBLE = 1
LMPIPTR = 2
LMPDPTR = 3
LMPDPTRPTR = 4

class lammps:
  def __init__(self,name="",cmdlineargs=None):

    # load liblmp.so by default
    # if name = "g++", load liblmp_g++.so
    
    try:
      if not name: self.lib = CDLL("liblmp.so")
      else: self.lib = CDLL("liblmp_%s.so" % name)
    except:
      raise OSError,"Could not load LAMMPS dynamic library"

    # create an instance of LAMMPS
    # don't know how to pass an MPI communicator from PyPar
    # no_mpi call lets LAMMPS use MPI_COMM_WORLD
    # cargs = array of C strings from args
    
    if cmdlineargs:
      cmdlineargs.insert(0,"lammps.py")
      narg = len(cmdlineargs)
      cargs = (c_char_p*narg)(*cmdlineargs)
      self.lmp = c_void_p()
      self.lib.lammps_open_no_mpi(narg,cargs,byref(self.lmp))
    else:
      self.lmp = c_void_p()
      self.lib.lammps_open_no_mpi(0,None,byref(self.lmp))
      # could use just this if LAMMPS lib interface supported it
      # self.lmp = self.lib.lammps_open_no_mpi(0,None)

  def __del__(self):
    if self.lmp: self.lib.lammps_close(self.lmp)

  def close(self):
    self.lib.lammps_close(self.lmp)
    self.lmp = None

  def file(self,file):
    self.lib.lammps_file(self.lmp,file)

  def command(self,cmd):
    self.lib.lammps_command(self.lmp,cmd)

  def extract_global(self,name,type):
    if type == LMPDOUBLE:
      self.lib.lammps_extract_global.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_global(self.lmp,name)
      return ptr[0]
    if type == LMPINT:
      self.lib.lammps_extract_global.restype = POINTER(c_int)
      ptr = self.lib.lammps_extract_global(self.lmp,name)
      return ptr[0]
    return None

  def extract_atom(self,name,type):
    if type == LMPDPTRPTR:
      self.lib.lammps_extract_atom.restype = POINTER(POINTER(c_double))
      ptr = self.lib.lammps_extract_atom(self.lmp,name)
      return ptr
    if type == LMPDPTR:
      self.lib.lammps_extract_atom.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_atom(self.lmp,name)
      return ptr
    if type == LMPIPTR:
      self.lib.lammps_extract_atom.restype = POINTER(c_int)
      ptr = self.lib.lammps_extract_atom(self.lmp,name)
      return ptr
    return None

  def extract_compute(self,id,style,type):
    if type == 0:
      if style > 0: return None
      self.lib.lammps_extract_compute.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
      return ptr[0]
    elif type == 1:
      self.lib.lammps_extract_compute.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
      return ptr
    elif type == 2:
      self.lib.lammps_extract_compute.restype = POINTER(POINTER(c_double))
      ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
      return ptr
    return None

  # in case of global datum, free memory for 1 double via lammps_free()
  # double was allocated by library interface function
  
  def extract_fix(self,id,style,type,i=0,j=0):
    if type == 0:
      if style > 0: return None
      self.lib.lammps_extract_fix.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_bix(self.lmp,id,style,type,i,j)
      result = ptr[0]
      self.lib.lammps_free(ptr)
      return result
    elif type == 1:
      self.lib.lammps_extract_fix.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,i,j)
      return ptr
    elif type == 2:
      self.lib.lammps_extract_fix.restype = POINTER(POINTER(c_double))
      ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,i,j)
      return ptr
    return None

  # free memory for 1 double or 1 vector of doubles via lammps_free()
  # for vector, must copy nlocal returned values to local c_double vector
  # memory was allocated by library interface function
  
  def extract_variable(self,name,group,type):
    if type == 0:
      self.lib.lammps_extract_variable.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_variable(self.lmp,name,group)
      result = ptr[0]
      self.lib.lammps_free(ptr)
      return result
    if type == 1:
      self.lib.lammps_extract_global.restype = POINTER(c_int)
      nlocalptr = self.lib.lammps_extract_global(self.lmp,"nlocal")
      nlocal = nlocalptr[0]
      result = (c_double*nlocal)()
      self.lib.lammps_extract_variable.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_variable(self.lmp,name,group)
      for i in xrange(nlocal): result[i] = ptr[i]
      self.lib.lammps_free(ptr)
      return result
    return None

  def get_natoms(self):
    return self.lib.lammps_get_natoms(self.lmp)

  def get_coords(self):
    nlen = 3 * self.lib.lammps_get_natoms(self.lmp)
    coords = (c_double*nlen)()
    self.lib.lammps_get_coords(self.lmp,coords)
    return coords

  # assume coords is an array of c_double, as created by get_coords()
  # could check if it is some other Python object and create c_double array?
  # constructor for c_double array can take an arg to use to fill it?
  
  def put_coords(self,coords):
    self.lib.lammps_put_coords(self.lmp,coords)
