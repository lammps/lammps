"""Contains the classes which deal with all the beads.

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


Used for holding information about the beads, including their positions, masses
momenta and kinetic energy. Has different objects for the position and normal
mode representations, and has a special centroid atoms object for when the
centroid coordinate is required.

Classes:
   Beads: Class with methods dealing with all the beads.
"""

__all__ = ['Beads']

import numpy as np
from ipi.utils.depend import *
from ipi.engine.atoms import Atoms
from ipi.utils import units

class Beads(dobject):
   """Storage for the beads positions and velocities.

   Everything is stored as (nbeads,3*natoms) sized contiguous arrays,
   and a convenience-access to each replica of the system is provided through a
   list of Atoms objects. Contains arrays of both the normal mode representation
   and the position representation, and various sized arrays for the atom
   labels and masses. Also contains the potential and force between
   neighbouring replicas.

   Attributes:
      natoms: The number of atoms.
      nbeads: The number of beads.
      _blist: A list of Atoms objects for each replica of the system. Each
         replica is assumed to have the same mass and atom label.
      centroid: An atoms object giving the centroid coordinate of the beads.

   Depend objects:
      names: An array giving the atom names.
      m: An array giving the atom masses.
      m3: An array giving the mass associated with each degree of freedom.
      sm3: An array giving the square root of m3.
      q: An array giving all the bead positions.
      p: An array giving all the bead momenta.
      qc: An array giving the centroid positions. Depends on qnm.
      pc: An array giving the centroid momenta. Depends on pnm.
      vpath: The spring potential between the beads, divided by omegan**2.
         Depends on q.
      fpath: The spring force between the beads, divided by omegan**2.
         Depends on q.
      kins: A list of the kinetic energy of each replica.
      kin: The total kinetic energy of the system. Note that this is not the
         same as the estimate of the kinetic energy of the system, which is
         contained in the properties module.
      kstress: The total kinetic stress tensor for the system.
      rg: An array giving the radius of gyration of each atom.
   """

   def __init__(self, natoms, nbeads):
      """Initializes Beads.

      Args:
         natoms: Number of atoms.
         nbeads: Number of beads.
      """

      self.resize(natoms, nbeads)

   def resize(self, natoms, nbeads):
      """Creates all the data arrays needed in the simulation.

      Effectively initializes the whole Beads object, according to the
      specified number of atoms and beads. Is also used, as the name suggests,
      to resize the data to a new number of beads when this is necessary, for
      example in initialization from a simulation with a different number of
      beads.

      Also creates, or recreates, the dependency network, as this requires
      the data arrays to be created for it to work.

      Args:
         natoms: The number of atoms.
         nbeads: The number of beads.
      """

      self.natoms = natoms
      self.nbeads = nbeads

      dset(self,"names",
         depend_array(name="names",value=np.zeros(natoms, np.dtype('|S6'))) )

      # atom masses, and mass-related arrays
      dset(self,"m",depend_array(name="m",value=np.zeros(natoms, float)) )   # this is the prototype mass array (just one independent of bead n)
      dset(self,"m3",
         depend_array(name="m3",value=np.zeros((nbeads,3*natoms), float),    # this is m conveniently replicated to be (nb,3*nat)
            func=self.mtom3, dependencies=[dget(self,"m")]))
      dset(self,"sm3",
         depend_array(name="sm3",value=np.zeros((nbeads,3*natoms), float),   # this is just the square root of m3
            func=self.m3tosm3, dependencies=[dget(self,"m3")]))

      # positions and momenta. bead representation, base storage used everywhere
      dset(self,"q",
         depend_array(name="q",value=np.zeros((nbeads,3*natoms), float)) )
      dset(self,"p",
         depend_array(name="p",value=np.zeros((nbeads,3*natoms), float)) )

      # position and momentum of the centroid
      dset(self,"qc",
         depend_array(name="qc",value=np.zeros(3*natoms, float),
            func=self.get_qc, dependencies=[dget(self,"q")] ) )
      dset(self,"pc",
         depend_array(name="pc",value=np.zeros(3*natoms, float),
            func=self.get_pc, dependencies=[dget(self,"p")] ) )

      # create proxies to access the centroid and the individual beads as Atoms objects
      self.centroid = Atoms(natoms, _prebind=(self.qc, self.pc, self.m, self.names))
      self._blist = [Atoms(natoms, _prebind=( self.q[i,:], self.p[i,:], self.m,  self.names )) for i in range(nbeads) ]

      # path springs potential and force
      dset(self,"vpath",
         depend_value(name="vpath", func=self.get_vpath,
            dependencies=[dget(self,"q")]))
      dset(self,"fpath",
         depend_array(name="fpath", value=np.zeros((nbeads,3*natoms), float),
            func=self.get_fpath, dependencies=[dget(self,"q")]))

      # kinetic energies of thhe beads, and total (classical) kinetic stress tensor
      dset(self,"kins",
         depend_array(name="kins",value=np.zeros(nbeads, float),
            func=self.kin_gather,
               dependencies=[dget(b,"kin") for b in self._blist]))
      dset(self,"kin",
         depend_value(name="kin", func=self.get_kin,
            dependencies=[dget(self,"kins")]))
      dset(self,"kstress",
         depend_array(name="kstress",value=np.zeros((3,3), float),
            func=self.get_kstress,
               dependencies=[dget(b,"kstress") for b in self._blist]))

   def copy(self):
      """Creates a new beads object from the original.

      Returns:
         A Beads object with the same q, p, m and names arrays as the original.
      """

      newbd = Beads(self.natoms, self.nbeads)
      newbd.q[:] = self.q
      newbd.p[:] = self.p
      newbd.m[:] = self.m
      newbd.names[:] = self.names
      return newbd


   def m3tosm3(self):
      """Takes the mass array and returns the square rooted mass array."""

      return np.sqrt(depstrip(self.m3))

   def mtom3(self):
      """Takes the mass array for each bead and returns one with an element
      for each degree of freedom.

      Returns:
         An array of size (nbeads,3*natoms), with each element corresponding
         to the mass associated with the appropriate degree of freedom in q.
      """

      m3 = np.zeros((self.nbeads,3*self.natoms),float)
      m3[:,0:3*self.natoms:3] = self.m
      m3[:,1:3*self.natoms:3] = m3[:,0:3*self.natoms:3]
      m3[:,2:3*self.natoms:3] = m3[:,0:3*self.natoms:3]
      return m3

   def get_qc(self):
      """Gets the centroid coordinates."""

      return np.dot(np.ones(self.nbeads,float),depstrip(self.q))/float(self.nbeads)

   def get_pc(self):
      """Gets the centroid momenta."""

      return np.dot(np.ones(self.nbeads,float),depstrip(self.p))/float(self.nbeads)

   def kin_gather(self):
      """Gets the kinetic energy for all the replicas.

      Returns:
         A list of the kinetic energy for each system.
      """

      return np.array([b.kin for b in self._blist])

   def get_kin(self):
      """Gets the total kinetic energy of all the replicas.

      Note that this does not correspond to the total kinetic energy estimate
      for the system.

      Returns:
         The sum of the kinetic energy of each replica.
      """

      return self.kins.sum()

   def get_kstress(self):
      """Calculates the total kinetic stress tensor of all the replicas.

      Note that this does not correspond to the quantum kinetic stress tensor
      estimate for the system.

      Returns:
         The sum of the kinetic stress tensor of each replica.
      """

      ks = np.zeros((3,3),float)
      for b in range(self.nbeads):
         ks += self[b].kstress
      return ks

   def get_vpath(self):
      """Calculates the spring potential between the replicas.

      Note that this is actually the harmonic potential without being
      multiplied by the factor omegan**2, which is only available in the
      ensemble as the temperature is required to calculate it.
      """

      epath = 0.0
      q = depstrip(self.q)
      m = depstrip(self.m3)[0]
      for b in range(self.nbeads):
         if b > 0:
            dq = q[b,:] - q[b-1,:]
         else:
            dq = q[b,:] - q[self.nbeads-1,:]
         epath += np.dot(dq, m*dq)
      return epath*0.5

   def get_fpath(self):
      """Calculates the spring force between the replicas.

      Note that this is actually the harmonic force without being
      multiplied by the factor omegan**2, which is only available in the
      ensemble as the temperature is required to calculate it.
      """

      nbeads = self.nbeads
      natoms = self.natoms
      f = np.zeros((nbeads,3*natoms),float)

      q = depstrip(self.q)
      m = depstrip(self.m3)[0]
      for b in range(nbeads):
         if b > 0:
            dq = q[b,:] - q[b-1,:]
         else:
            dq = q[b,:] - q[self.nbeads-1,:]
         dq *= m
         f[b] -= dq
         if b > 0:
            f[b-1] += dq
         else:
            f[nbeads-1] += dq
      return f

   # A set of functions to access individual beads as Atoms objects
   def __len__(self):
      """Length function.

      This is called whenever the standard function len(beads) is used.

      Returns:
         The number of beads.
      """

      return self.nbeads

   def __getitem__(self,index):
      """Overwrites standard getting function.

      This is called whenever the standard function beads[index] is used.
      Returns an Atoms object with the appropriate position and momenta arrays.

      Args:
         index: The index of the replica of the system to be accessed.

      Returns:
         The replica of the system given by the index.
      """

      return self._blist[index]

   def __setitem__(self,index,value):
      """Overwrites standard setting function.

      This is called whenever the standard function beads[index]=value is used.
      Changes the position and momenta of the appropriate slice of the global
      position and momentum arrays to those given by value.

      Args:
         index: The replica of the system to be changed.
         value: The Atoms object that holds the new values.
      """

      self._blist[index].p[:] = value.p
      self._blist[index].q[:] = value.q
      self._blist[index].m[:] = value.m
      self._blist[index].names[:] = value.names
