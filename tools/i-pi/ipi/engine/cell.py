"""Contains the classes which deal with the system box.

Copyright (C) 2013, Joshua More and Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.


Used for implementing the minimum image convention.

Classes:
   Cell: Base cell class with the generic methods and attributes.
"""

__all__ = ['Cell']

import numpy as np
from ipi.utils.depend import *
from ipi.utils.mathtools import *
from ipi.utils import units


class Cell(dobject):
   """Base class to represent the simulation cell in a periodic system.

   This class has the base attributes required for either flexible or
   isotropic cell dynamics. Uses an upper triangular lattice vector matrix to
   represent the cell.

   Depend objects:
      h: An array giving the lattice vector matrix.
      ih: An array giving the inverse of the lattice vector matrix.
      V: The volume of the cell.
   """

   def __init__(self, h=None):
      """Initializes base cell class.

      Args:
         h: Optional array giving the initial lattice vector matrix. The
            reference cell matrix is set equal to this. Must be an upper
            triangular 3*3 matrix. Defaults to a 3*3 zeroes matrix.
      """

      if h is None:
         #h = np.identity(3,float)
         h = np.zeros((3,3), float)
      dset(self,"h",depend_array(name = 'h', value = h) )

      dset(self,"ih",
         depend_array(name = "ih", value = np.zeros((3,3),float),
            func=self.get_ih, dependencies=[dget(self,"h")]) )
      dset(self,"V",
         depend_value(name = 'V', func=self.get_volume,
            dependencies=[dget(self,"h")]) )

   def get_ih(self):
      """Inverts the lattice vector matrix."""

      return invert_ut3x3(self.h)

   def get_volume(self):
      """Calculates the volume of the system box."""

      return det_ut3x3(self.h)

   def apply_pbc(self, atom):
      """Uses the minimum image convention to return a particle to the
         unit cell.

      Args:
         atom: An Atom object.

      Returns:
         An array giving the position of the image that is inside the
         system box.
      """

      s = np.dot(self.ih,atom.q)


      for i in range(3):
         s[i] = s[i] - round(s[i])

      return np.dot(self.h,s)

   def array_pbc(self, pos):
      """Uses the minimum image convention to return a list of particles to the
         unit cell.

      Args:
         atom: An Atom object.

      Returns:
         An array giving the position of the image that is inside the
         system box.
      """

      s = depstrip(pos).copy()
      s.shape = (len(pos)/3,3)

      s = np.dot(depstrip(self.ih),s.T)
      s = s - np.round(s)

      s = np.dot(depstrip(self.h),s).T

      pos[:] = s.reshape((len(s)*3))

   def minimum_distance(self, atom1, atom2):
      """Takes two atoms and tries to find the smallest vector between two
      images.

      This is only rigorously accurate in the case of a cubic cell,
      but gives the correct results as long as the cut-off radius is defined
      as smaller than the smallest width between parallel faces even for
      triclinic cells.

      Args:
         atom1: An Atom object.
         atom2: An Atom object.

      Returns:
         An array giving the minimum distance between the positions of atoms
         atom1 and atom2 in the minimum image convention.
      """

      s = np.dot(self.ih,atom1.q-atom2.q)
      for i in range(3):
         s[i] -= round(s[i])
      return np.dot(self.h, s)
