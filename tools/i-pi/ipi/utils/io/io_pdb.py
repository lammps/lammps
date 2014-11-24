"""Contains the functions used to print the trajectories and read input
configurations with pdb formatting.

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
   print_pdb_path: Prints all the bead configurations, and shows the ring
      polymer connectivity.
   print_pdb: Prints the centroid configurations.
   read_pdb: Reads the cell parameters and atom configurations from a pdb file.
"""

__all__ = ['print_pdb_path', 'print_pdb', 'read_pdb']

import numpy as np
import sys
import ipi.utils.mathtools as mt
from ipi.utils.depend  import depstrip
from ipi.engine.cell   import Cell
from ipi.engine.atoms  import Atoms
from ipi.utils.units   import *

def print_pdb_path(beads, cell, filedesc = sys.stdout):
   """Prints all the bead configurations, into a pdb formatted file.

   Prints the ring polymer springs as well as the bead positions using the
   CONECT command. Also prints the cell parameters in standard pdb form. Note
   that the angles are in degrees.

   Args:
      beads: A beads object giving the bead positions.
      cell: A cell object giving the system box.
      filedesc: An open writable file object. Defaults to standard output.
   """

   a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h)

   z = 1 #What even is this parameter?
   filedesc.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s%4i\n" % (a, b, c, alpha, beta, gamma, " P 1        ", z))

   natoms = beads.natoms
   nbeads = beads.nbeads
   for j in range(nbeads):
      for i in range(natoms):
         qs = depstrip(beads.q)
         lab = depstrip(beads.names)
         filedesc.write("ATOM  %5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2i\n" % (j*natoms+i+1, lab[i],' ','  1',' ',1,' ', qs[j][3*i], qs[j][3*i+1], qs[j][3*i+2],0.0,0.0,'  ',0))

   if nbeads > 1:
      for i in range(natoms):
         filedesc.write("CONECT%5i%5i\n" % (i+1, (nbeads-1)*natoms+i+1))
      for j in range(nbeads-1):
         for i in range(natoms):
            filedesc.write("CONECT%5i%5i\n" % (j*natoms+i+1, (j+1)*natoms+i+1))

   filedesc.write("END\n")

def print_pdb(atoms, cell, filedesc = sys.stdout, title=""):
   """Prints the atom configurations, into a pdb formatted file.

   Also prints the cell parameters in standard pdb form. Note
   that the angles are in degrees.

   Args:
      atoms: An atoms object giving the atom positions.
      cell: A cell object giving the system box.
      filedesc: An open writable file object. Defaults to standard output.
      title: An optional string of max. 70 characters.
   """


   if title != "" :
      filedesc.write("TITLE   %70s\n" % (title))

   a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h)

   z = 1
   filedesc.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s%4i\n" % (a, b, c, alpha, beta, gamma, " P 1        ", z))

   natoms = atoms.natoms
   qs = depstrip(atoms.q)
   lab = depstrip(atoms.names)
   for i in range(natoms):
      filedesc.write("ATOM  %5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2i\n" % (i+1, lab[i], ' ', '  1', ' ', 1, ' ', qs[3*i], qs[3*i+1], qs[3*i+2], 0.0, 0.0, '  ', 0))

   filedesc.write("END\n")

def read_pdb(filedesc):
   """Takes a pdb-style file and creates an Atoms and Cell object.

   Args:
      filedesc: An open readable file object from a pdb formatted file.

   Returns:
      An Atoms object with the appropriate atom labels, masses and positions,
      and a Cell object with the appropriate cell dimensions and an estimate
      of a reasonable cell mass.
   """

   header = filedesc.readline()
   if "TITLE" in header: header = filedesc.readline()   # skip the comment field
   if header == "":
      raise EOFError("End of file or empty header in PDB file")

   a = float(header[6:15])
   b = float(header[15:24])
   c = float(header[24:33])
   alpha = float(header[33:40])
   beta = float(header[40:47])
   gamma = float(header[47:54])
   alpha *= np.pi/180.0
   beta *= np.pi/180.0
   gamma *= np.pi/180.0
   h = mt.abc2h(a, b, c, alpha, beta, gamma)
   cell = Cell(h)

   natoms = 0
   body = filedesc.readline()
   qatoms = []
   names = []
   masses = []
   while (body.strip() != "" and body.strip() != "END"):
      natoms += 1
      name = body[12:16].strip()
      names.append(name)
      masses.append(Elements.mass(name))
      x = float(body[31:39])
      y = float(body[39:47])
      z = float(body[47:55])
      qatoms.append(x)
      qatoms.append(y)
      qatoms.append(z)

      body = filedesc.readline()

   atoms = Atoms(natoms)
   atoms.q = np.asarray(qatoms)
   atoms.names = np.asarray(names,dtype='|S4')
   atoms.m = np.asarray(masses)

   return atoms, cell

def iter_pdb(filedesc):
   """Takes a pdb-style file and yields one Atoms, Cell tuple after another.

   Args:
      filedesc: An open readable file object from a pdb formatted file.

   Returns:
      Generator over the pdb trajectory, that yields
      (Atoms, Cell) tuple with the appropriate atom labels, masses and positions.
   """

   try:
      while 1:
         atoms, cell = read_pdb(filedesc)
         yield atoms, cell
   except EOFError:
      pass
