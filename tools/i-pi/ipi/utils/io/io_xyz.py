"""Contains the functions used to print the trajectories and read input
configurations with xyz formatting.

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
   print_xyz_path: Prints all the bead configurations.
   print_xyz: Prints the centroid configurations.
   read_xyz: Reads the cell parameters and atom configurations from a xyz file.
"""

__all__ = ['print_xyz_path', 'print_xyz', 'read_xyz', 'iter_xyz']

import numpy as np
import math, sys
import ipi.utils.mathtools as mt
from ipi.utils.depend import depstrip
from ipi.engine.atoms import Atoms
from ipi.utils.units import *

def print_xyz_path(beads, cell, filedesc = sys.stdout):
   """Prints all the bead configurations, into a xyz formatted file.

   Prints all the replicas for each time step separately, rather than all at
   once.

   Args:
      beads: A beads object giving the bead positions.
      cell: A cell object giving the system box.
      filedesc: An open writable file object. Defaults to standard output.
   """

   a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h)

   natoms = beads.natoms
   nbeads = beads.nbeads
   for j in range(nbeads):
      filedesc.write("%d\n# bead: %d CELL(abcABC): %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f \n" % (natoms, j, a, b, c, alpha, beta, gamma))
      for i in range(natoms):
         qs = depstrip(beads.q)
         lab = depstrip(beads.names)
         filedesc.write("%8s %12.5e %12.5e %12.5e\n" % (lab[i], qs[j][3*i], qs[j][3*i+1], qs[j][3*i+2]))

def print_xyz(atoms, cell, filedesc = sys.stdout, title=""):
   """Prints the centroid configurations, into a xyz formatted file.

   Args:
      atoms: An atoms object giving the centroid positions.
      cell: A cell object giving the system box.
      filedesc: An open writable file object. Defaults to standard output.
      title: This gives a string to be appended to the comment line.
   """

   a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h)

   natoms = atoms.natoms
   filedesc.write("%d\n# CELL(abcABC): %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %s\n" % ( natoms, a, b, c, alpha, beta, gamma, title))
   # direct access to avoid unnecessary slow-down
   qs = depstrip(atoms.q)
   lab = depstrip(atoms.names)
   for i in range(natoms):
      filedesc.write("%8s %12.5e %12.5e %12.5e\n" % (lab[i], qs[3*i], qs[3*i+1], qs[3*i+2]))

def read_xyz(filedesc):
   """Takes a xyz-style file and creates an Atoms object.

   Args:
      filedesc: An open readable file object from a xyz formatted file.

   Returns:
      An Atoms object with the appropriate atom labels, masses and positions.
   """

   natoms = filedesc.readline()
   if natoms == "":
      raise EOFError("The file descriptor hit EOF.")
   natoms = int(natoms)
   comment = filedesc.readline()

   qatoms = []
   names = []
   masses = []
   iat = 0
   while (iat < natoms):
      body = filedesc.readline()
      if body.strip() == "":
         break
      body = body.split()
      name = body[0]
      names.append(name)
      masses.append(Elements.mass(name))
      x = float(body[1])
      y = float(body[2])
      z = float(body[3])
      qatoms.append(x)
      qatoms.append(y)
      qatoms.append(z)
      iat += 1

   if natoms != len(names):
      raise ValueError("The number of atom records does not match the header of the xyz file.")

   atoms = Atoms(natoms)
#   for i in range(natoms):
#      nat = atoms[i]
#      nat.q = qatoms[i]
#      nat.name = names[i]
#      nat.m = Elements.mass(names[i])
   atoms.q = np.asarray(qatoms)
   atoms.names = np.asarray(names, dtype='|S4')
   atoms.m = np.asarray(masses)

   return atoms

def iter_xyz(filedesc):
   """Takes a xyz-style file and yields one Atoms object after another.

   Args:
      filedesc: An open readable file object from a xyz formatted file.

   Returns:
      Generator over the xyz trajectory, that yields
      Atoms objects with the appropriate atom labels, masses and positions.
   """

   try:
      while 1:
         atoms = read_xyz(filedesc)
         yield atoms
   except EOFError:
      pass
