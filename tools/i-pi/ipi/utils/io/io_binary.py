"""Contains the functions used to print the trajectories and read input
configurations (or even full status dump) as unformatted binary.

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


Functions:
   print_bin: Prints an atomic configuration.
"""

__all__ = ['print_bin']

import os
import numpy as np
import math, sys
from ipi.utils.depend import depstrip

def print_bin(atoms, cell, filedesc = sys.stdout, title=""):
   """Prints the centroid configurations, into a binary file.

   Args:
      beads: An atoms object giving the centroid positions.
      cell: A cell object giving the system box.
      filedesc: An open writable file object. Defaults to standard output.
      title: This gives a string to be appended to the comment line.
   """

   buff = filedesc # .buffer
   cell.h.tofile(buff)
   nat = np.asarray([atoms.natoms])
   nat.tofile(buff)
   atoms.names.tofile(buff)
   atoms.q.tofile(buff)

