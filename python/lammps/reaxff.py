# -------------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   https://www.lammps.org/ Sandia National Laboratories
#   LAMMPS Development team: developers@lammps.org
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
# -------------------------------------------------------------------------

###########################################################################
#   Written by Mitch Murphy alphataubio at gmail
#   MINIMAL changes needed for FitSNAP-ReaxFF
###########################################################################

from lammps.core import lammps
from ctypes import c_void_p, c_int, c_double

# -------------------------------------------------------------------------

class lammps_reaxff(lammps):

  def __init__(self,name='',cmdargs=None,ptr=None,comm=None):

    super().__init__(name=name,cmdargs=cmdargs,ptr=ptr,comm=comm)

    #if has_package('REAXFF')

    self.lib.lammps_set_reaxff_atm_parameter.argtypes = [c_void_p,c_int,c_int,c_double]
    self.lib.lammps_set_reaxff_bnd_parameter.argtypes = [c_void_p,c_int,c_int,c_int,c_double]
    self.lib.lammps_set_reaxff_ofd_parameter.argtypes = [c_void_p,c_int,c_int,c_int,c_double]
    self.lib.lammps_set_reaxff_ang_parameter.argtypes = [c_void_p,c_int,c_int,c_int,c_int,c_double]
    self.lib.lammps_set_reaxff_tor_parameter.argtypes = [c_void_p,c_int,c_int,c_int,c_int,c_int,c_double]
    self.lib.lammps_set_reaxff_hbd_parameter.argtypes = [c_void_p,c_int,c_int,c_int,c_double]
    self.lib.lammps_set_reaxff_atm_parameter.restype = None
    self.lib.lammps_set_reaxff_bnd_parameter.restype = None
    self.lib.lammps_set_reaxff_ofd_parameter.restype = None
    self.lib.lammps_set_reaxff_ang_parameter.restype = None
    self.lib.lammps_set_reaxff_tor_parameter.restype = None
    self.lib.lammps_set_reaxff_hbd_parameter.restype = None


  def set_reaxff_parameters(self, parameters):

    #with ExceptionCheck(self):

      # ATM
      for p in parameters[0]:
        self.lib.lammps_set_reaxff_atm_parameter(self.lmp,p[0],p[1],p[2])

      # BND
      for p in parameters[1]:
        self.lib.lammps_set_reaxff_bnd_parameter(self.lmp,p[0],p[1],p[2],p[3])

      # OFD
      for p in parameters[2]:
        self.lib.lammps_set_reaxff_ofd_parameter(self.lmp,p[0],p[1],p[2],p[3])

      # ANG
      for p in parameters[3]:
        self.lib.lammps_set_reaxff_ang_parameter(self.lmp,p[0],p[1],p[2],p[3],p[4])

      # TOR
      for p in parameters[4]:
        self.lib.lammps_set_reaxff_tor_parameter(self.lmp,p[0],p[1],p[2],p[3],p[4],p[5])

      # HBD
      for p in parameters[5]:
        self.lib.lammps_set_reaxff_hbd_parameter(self.lmp,p[0],p[1],p[2],p[3],p[4])

