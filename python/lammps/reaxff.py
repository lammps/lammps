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

from lammps.core import lammps, ExceptionCheck
from ctypes import c_void_p, c_int, c_double, POINTER

# -------------------------------------------------------------------------

class lammps_reaxff(lammps):

  # -------------------------------------------------------------------------

  def __init__(self,name='',cmdargs=None,ptr=None,comm=None):

    super().__init__(name=name,cmdargs=cmdargs,ptr=ptr,comm=comm)

    #if has_package('REAXFF')

    self.lib.lammps_set_reaxff_gen_parameter.argtypes = [c_void_p,c_int,c_double]
    self.lib.lammps_set_reaxff_atm_parameter.argtypes = [c_void_p,c_int,c_int,c_double]
    self.lib.lammps_set_reaxff_bnd_parameter.argtypes = [c_void_p,c_int,c_int,c_int,c_double]
    self.lib.lammps_set_reaxff_ofd_parameter.argtypes = [c_void_p,c_int,c_int,c_int,c_double]
    self.lib.lammps_set_reaxff_ang_parameter.argtypes = [c_void_p,c_int,c_int,c_int,c_int,c_double]
    self.lib.lammps_set_reaxff_tor_parameter.argtypes = [c_void_p,c_int,c_int,c_int,c_int,c_int,c_double]
    self.lib.lammps_set_reaxff_hbd_parameter.argtypes = [c_void_p,c_int,c_int,c_int,c_int,c_double]
    self.lib.lammps_set_reaxff_gen_parameter.restype = None
    self.lib.lammps_set_reaxff_atm_parameter.restype = None
    self.lib.lammps_set_reaxff_bnd_parameter.restype = None
    self.lib.lammps_set_reaxff_ofd_parameter.restype = None
    self.lib.lammps_set_reaxff_ang_parameter.restype = None
    self.lib.lammps_set_reaxff_tor_parameter.restype = None
    self.lib.lammps_set_reaxff_hbd_parameter.restype = None

    self.lib.lammps_reset_box_single.argtypes = \
      [c_void_p,POINTER(c_double),POINTER(c_double),c_double,c_double,c_double]
    self.lib.lammps_reset_box_single.restype = None

  # -------------------------------------------------------------------------

  def set_reaxff_parameters(self, parameters, values):

    with ExceptionCheck(self):

      for p, v in zip(parameters, values):

        p_block_index = p[0]

        # GEN
        if p_block_index == 0:
          self.lib.lammps_set_reaxff_gen_parameter(self.lmp,p[1],v)

        # ATM
        if p_block_index == 1:
          self.lib.lammps_set_reaxff_atm_parameter(self.lmp,p[1],p[2],v)

        # BND
        if p_block_index == 2:
          self.lib.lammps_set_reaxff_bnd_parameter(self.lmp,p[1],p[2],p[3],v)

        # OFD
        if p_block_index == 3:
          self.lib.lammps_set_reaxff_ofd_parameter(self.lmp,p[1],p[2],p[3],v)

        # ANG
        if p_block_index == 4:
          self.lib.lammps_set_reaxff_ang_parameter(self.lmp,p[1],p[2],p[3],p[4],v)

        # TOR
        if p_block_index == 5:
          self.lib.lammps_set_reaxff_tor_parameter(self.lmp,p[1],p[2],p[3],p[4],p[5],v)

        # HBD
        if p_block_index == 6:
          self.lib.lammps_set_reaxff_hbd_parameter(self.lmp,p[1],p[2],p[3],p[4],v)

  # -------------------------------------------------------------------------

  def reset_box_single(self,boxlo,boxhi,xy,yz,xz):
    """Reset simulation box parameters

    This is a wrapper around the :cpp:func:`lammps_reset_box_single` function
    of the C-library reaxff interface.

    :param boxlo: new lower box boundaries
    :type boxlo: list of 3 floating point numbers
    :param boxhi: new upper box boundaries
    :type boxhi: list of 3 floating point numbers
    :param xy: xy tilt factor
    :type xy: float
    :param yz: yz tilt factor
    :type yz: float
    :param xz: xz tilt factor
    :type xz: float
    """
    cboxlo = (3*c_double)(*boxlo)
    cboxhi = (3*c_double)(*boxhi)
    with ExceptionCheck(self):
      self.lib.lammps_reset_box_single(self.lmp,cboxlo,cboxhi,xy,yz,xz)
