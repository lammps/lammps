"""Contains the classes that connect the driver to the python code.

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


Communicates with the driver code, obtaining the force, virial and potential.
Deals with creating the jobs that will be sent to the driver, and
returning the results to the python code.

Classes:
   ForceField: Base forcefield class with the generic methods and attributes.
   FFSocket: Deals with a single replica of the system
   ForceBeads: Deals with the parallelization of the force calculation over
      different beads.
   Forces: Deals with the parallelizatoin of the force calculation over
      different forcefields.
"""

__all__ = ['ForceField', 'ForceBeads', 'Forces', 'FFSocket']

import numpy as np
import time
from ipi.utils.softexit import softexit
from ipi.utils.messages import verbosity, warning
from ipi.utils.depend import *
from ipi.utils.nmtransform import nm_rescale
from ipi.interfaces.sockets import InterfaceSocket
from ipi.engine.beads import Beads

class ForceField(dobject):
   """Base forcefield class.

   Gives the standard methods and quantities needed in all the forcefield
   classes.

   Attributes:
      atoms: An Atoms object containing all the atom positions.
      cell: A Cell object containing the system box.

   Depend objects:
      ufvx: A list of the form [pot, f, vir]. These quantities are calculated
         all at one time by the driver, so are collected together. Each separate
         object is then taken from the list. Depends on the atom positions and
         the system box.
      extra: A string containing some formatted output returned by the client. Depends on ufvx.
      pot: A float giving the potential energy of the system. Depends on ufvx.
      f: An array containing all the components of the force. Depends on ufvx.
      fx: A slice of f containing only the x components of the forces.
      fy: A slice of f containing only the y components of the forces.
      fz: A slice of f containing only the z components of the forces.
      vir: An array containing the components of the virial tensor in upper
         triangular form, not divided by the volume. Depends on ufvx.
   """

   def __init__(self):
      """Initializes ForceField."""

      # ufvx is a list [ u, f, vir, extra ]  which stores the results of the force
      #calculation
      dset(self,"ufvx", depend_value(name="ufvx", func=self.get_all))

   def copy(self):
      """Creates a deep copy without the bound objects.

      Used in ForceBeads to create a ForceField for each replica of the system.

      Returns:
         A ForceField object without atoms or cell attributes.
      """

      return type(self)(self.nbeads, self.weight)

   def bind(self, atoms, cell):
      """Binds atoms and cell to the forcefield.

      This takes an atoms object and a cell object and makes them members of
      the forcefield. It also then creates the objects that will hold the data
      that the driver returns and the dependency network.

      Args:
         atoms: The Atoms object from which the atom positions are taken.
         cell: The Cell object from which the system box is taken.
      """

      # stores a reference to the atoms and cell we are computing forces for
      self.atoms = atoms
      self.cell = cell

      # ufv depends on the atomic positions and on the cell
      dget(self,"ufvx").add_dependency(dget(self.atoms,"q"))
      dget(self,"ufvx").add_dependency(dget(self.cell,"h"))

      # potential and virial are to be extracted very simply from ufv
      dset(self,"pot",
         depend_value(name="pot", func=self.get_pot,
            dependencies=[dget(self,"ufvx")]))

      dset(self,"vir",
         depend_array(name="vir", value=np.zeros((3,3),float),func=self.get_vir,
            dependencies=[dget(self,"ufvx")]))

      # NB: the force requires a bit more work, to define shortcuts to xyz
      # slices without calculating the force at this point.
      fbase = np.zeros(atoms.natoms*3, float)
      dset(self,"f",
         depend_array(name="f", value=fbase, func=self.get_f,
             dependencies=[dget(self,"ufvx")]))

      dset(self,"extra",
         depend_value(name="extra", func=self.get_extra,
            dependencies=[dget(self,"ufvx")]))

      dset(self,"fx", depend_array(name="fx", value=fbase[0:3*atoms.natoms:3]))
      dset(self,"fy", depend_array(name="fy", value=fbase[1:3*atoms.natoms:3]))
      dset(self,"fz", depend_array(name="fz", value=fbase[2:3*atoms.natoms:3]))
      depcopy(self,"f", self,"fx")
      depcopy(self,"f", self,"fy")
      depcopy(self,"f", self,"fz")

   def queue(self):
      """Dummy queueing method."""

      pass

   def stop(self):
      """Dummy queueing method."""

      pass

   def run(self):
      """Dummy queueing method."""

      pass

   def get_all(self):
      """Dummy driver routine.

      Returns:
         A list of the form [potential, force, virial] where the potential
         and all components of the force and virial have been set to zero.
      """

      return [0.0, np.zeros(3*self.atoms.natoms), np.zeros((3,3),float), ""]

   def get_pot(self):
      """Calls get_all routine of forcefield to update potential.

      Returns:
         Potential energy.
      """

      return self.ufvx[0]

   def get_f(self):
      """Calls get_all routine of forcefield to update force.

      Returns:
         An array containing all the components of the force.
      """

      return depstrip(self.ufvx[1])

   def get_vir(self):
      """Calls get_all routine of forcefield to update virial.

      Returns:
         An array containing the virial in upper triangular form, not divided
         by the volume.
      """

      vir = depstrip(self.ufvx[2])
      vir[1,0] = 0.0
      vir[2,0:2] = 0.0
      return vir

   def get_extra(self):
      """Calls get_all routine of forcefield to update potential.

      Returns:
         A string containing all formatted additional output that the
         client might have produced.
      """

      return self.ufvx[3]


class FFSocket(ForceField):
   """Interface between the PIMD code and the socket for a single replica.

   Deals with an individual replica of the system, obtaining the potential
   force and virial appropriate to this system. Deals with the distribution of
   jobs to the interface.

   Attributes:
      parameters: A dictionary of the parameters used by the driver. Of the
         form {'name': value}.
      socket: The interface object which contains the socket through which
         communication between the forcefield and the driver is done.
      request: During the force calculation step this holds a dictionary
         containing the relevant data for determining the progress of the step.
         Of the form {'atoms': atoms, 'cell': cell, 'pars': parameters,
                      'status': status, 'result': result, 'id': bead id,
                      'start': starting time}.
   """

   def __init__(self, pars=None, interface=None):
      """Initializes FFSocket.

      Args:
         pars: Optional dictionary, giving the parameters needed by the driver.
         interface: Optional Interface object, which contains the socket.
      """

      # a socket to the communication library is created or linked
      super(FFSocket,self).__init__()
      if interface is None:
         self.socket = InterfaceSocket()
      else:
         self.socket = interface

      if pars is None:
         self.pars = {}
      else:
         self.pars = pars
      self.request = None

   def bind(self, atoms, cell):
      """Pass on the binding request from ForceBeads.

      Also makes sure to set the socket's softexit.

      Args:
         atoms: Atoms object from which the bead positions are taken.
         cell: Cell object from which the system box is taken.
      """

      super(FFSocket,self).bind(atoms, cell)

   def copy(self):
      """Creates a deep copy without the bound objects.

      Used in ForceBeads to create a FFSocket for each replica of the system.

      Returns:
         A FFSocket object without atoms or cell attributes.
      """

      # does not copy the bound objects
      # (i.e., the returned forcefield must be bound before use)
      return type(self)(self.pars, self.socket)

   def get_all(self):
      """Driver routine.

      When one of the force, potential or virial are called, this sends the
      atoms and cell to the driver through the interface, requesting that the
      driver does the calculation. This then waits until the driver is finished,
      and then returns the ufvx list.

      Returns:
         A list of the form [potential, force, virial, extra].
      """

      # this is converting the distribution library requests into [ u, f, v ]  lists
      if self.request is None:
         self.request = self.socket.queue(self.atoms, self.cell, pars=self.pars, reqid=-1)
      while self.request["status"] != "Done":
         if self.request["status"] == "Exit":
            break
         time.sleep(self.socket.latency)
      if self.request["status"] == "Exit":
         softexit.trigger(" @Force: Requested returned a Exit status")

      # data has been collected, so the request can be released and a slot
      #freed up for new calculations
      self.socket.release(self.request)
      result = self.request["result"]
      self.request = None

      return result

   def queue(self, reqid=-1):
      """Sends the job to the interface queue directly.

      Allows the ForceBeads object to ask for the ufvx list of each replica
      directly without going through the get_all function. This allows
      all the jobs to be sent at once, allowing them to be parallelized.

      Args:
         reqid: An optional integer that indentifies requests of the same type,
            e.g. the bead index.
      """

      if self.request is None and dget(self,"ufvx").tainted():
         self.request = self.socket.queue(self.atoms, self.cell, pars=self.pars, reqid=reqid)

   def run(self):
      """Makes the socket start looking for driver codes.

      Tells the interface code to start the thread that looks for
      connection from the driver codes in a loop. Until this point no
      jobs can be queued.
      """

      if not self.socket.started():
         self.socket.start_thread()

   def stop(self):
      """Makes the socket stop looking for driver codes.

      Tells the interface code to stop the thread that looks for
      connection from the driver codes in a loop. After this point no
      jobs can be queued.
      """

      if self.socket.started():
         self.socket.end_thread()


class ForceBeads(dobject):
   """Class that gathers the forces for each replica together.

   Deals with splitting the bead representation into
   separate replicas, and collecting the data from each replica.

   Attributes:
      natoms: An integer giving the number of atoms.
      nbeads: An integer giving the number of beads.
      f_model: A model used to create the forcefield objects for each replica
         of the system.
      _forces: A list of the forcefield objects for all the replicas.
      weight: A float that will be used to weight the contribution of this
         forcefield to the total force.

   Depend objects:
      f: An array containing the components of the force. Depends on each
         replica's ufvx list.
      pots: A list containing the potential energy for each system replica.
         Depends on each replica's ufvx list.
      virs: A list containing the virial tensor for each system replica.
         Depends on each replica's ufvx list.
      pot: The sum of the potential energy of the replicas.
      vir: The sum of the virial tensor of the replicas.
      extras: Strings containing some formatted output returned by the client.
         Depends on each replica's ufvx list.
   """

   def __init__(self, model, nbeads=0, weight=1.0):
      """Initializes ForceBeads

      Args:
         model: A model to be used to create the forcefield objects for all
            the replicas of the system.
         nbeads: The number of replicas.
         weight: A relative weight to be given to the values obtained with this
            forcefield. When the contribution of all the forcefields is
            combined to give a total force, the contribution of this forcefield
            will be weighted by this factor.
      """

      self.f_model = model
      self.nbeads = nbeads
      self.weight = weight

   def copy(self):
      """Creates a deep copy without the bound objects.

      Used so that we can create multiple Forces objects from the same
      Forcebeads model, without binding a particular ForceBeads object twice.

      Returns:
         A ForceBeads object without beads or cell attributes.
      """

      # does not copy the bound objects (i.e., the returned forcefield must be bound before use)
      return type(self)(self.f_model, self.nbeads, self.weight)


   def bind(self, beads, cell):
      """Binds beads, cell and force to the forcefield.

      Takes the beads, cell objects and makes them members of the forcefield.
      Also takes the force object and copies it once for each replica of the
      system, then binds each replica to one of the copies so that the force
      calculation can be parallelized. Creates the objects that will
      hold the data that the driver returns and the dependency network.

      Args:
         beads: Beads object from which the bead positions are taken.
         cell: Cell object from which the system box is taken.
      """

      # stores a copy of the number of atoms and of beads
      #!TODO! make them read-only properties
      self.natoms = beads.natoms
      if (self.nbeads != beads.nbeads):
         raise ValueError("Binding together a Beads and a ForceBeads objects with different numbers of beads")

      # creates an array of force objects, which are bound to the beads
      #and the cell
      self._forces = [];
      for b in range(self.nbeads):
         new_force = self.f_model.copy()
         new_force.bind(beads[b], cell)
         self._forces.append(new_force)

      # f is a big array which assembles the forces on individual beads
      dset(self,"f",
         depend_array(name="f",value=np.zeros((self.nbeads,3*self.natoms)),
            func=self.f_gather,
               dependencies=[dget(self._forces[b],"f") for b in range(self.nbeads)]))

      # collection of pots and virs from individual beads
      dset(self,"pots",
         depend_array(name="pots", value=np.zeros(self.nbeads,float),
            func=self.pot_gather,
               dependencies=[dget(self._forces[b],"pot") for b in range(self.nbeads)]))
      dset(self,"virs",
         depend_array(name="virs", value=np.zeros((self.nbeads,3,3),float),
            func=self.vir_gather,
               dependencies=[dget(self._forces[b],"vir") for b in range(self.nbeads)]))
      dset(self,"extras",
         depend_value(name="extras", value=np.zeros(self.nbeads,float),
            func=self.extra_gather,
               dependencies=[dget(self._forces[b],"extra") for b in range(self.nbeads)]))

      # total potential and total virial
      dset(self,"pot",
         depend_value(name="pot", func=(lambda: self.pots.sum()),
            dependencies=[dget(self,"pots")]))
      dset(self,"vir",
         depend_array(name="vir", func=self.get_vir, value=np.zeros((3,3)),
            dependencies=[dget(self,"virs")]))

   def run(self):
      """Makes the socket start looking for driver codes.

      Tells the interface code to start the thread that looks for
      connection from the driver codes in a loop. Until this point no
      jobs can be queued.
      """

      for b in range(self.nbeads):
         self._forces[b].run()

   def stop(self):
      """Makes the socket stop looking for driver codes.

      Tells the interface code to stop the thread that looks for
      connection from the driver codes in a loop. After this point no
      jobs can be queued.
      """

      for b in range(self.nbeads):
         self._forces[b].stop()

   def queue(self):
      """Submits all the required force calculations to the interface."""

      # this should be called in functions which access u,v,f for ALL the beads,
      # before accessing them. it is basically pre-queueing so that the
      # distributed-computing magic can work
      for b in range(self.nbeads):
         self._forces[b].queue(reqid=b)

   def pot_gather(self):
      """Obtains the potential energy for each replica.

      Returns:
         A list of the potential energy of each replica of the system.
      """

      self.queue()
      return np.array([b.pot for b in self._forces], float)

   def extra_gather(self):
      """Obtains the potential energy for each replica.

      Returns:
         A list of the potential energy of each replica of the system.
      """

      self.queue()
      return [b.extra for b in self._forces]

   def vir_gather(self):
      """Obtains the virial for each replica.

      Returns:
         A list of the virial of each replica of the system.
      """

      self.queue()
      return np.array([b.vir for b in self._forces], float)

   def f_gather(self):
      """Obtains the force vector for each replica.

      Returns:
         An array with all the components of the force. Row i gives the force
         array for replica i of the system.
      """

      newf = np.zeros((self.nbeads,3*self.natoms),float)

      self.queue()
      for b in range(self.nbeads):
         newf[b] = depstrip(self._forces[b].f)

      return newf

      #serial
#      for b in range(self.nbeads): newf[b]=self._forces[b].f
      # threaded
#      bthreads=[]
#      print "starting threads"
#      for b in range(self.nbeads):
#         thread=threading.Thread(target=self._getbead, args=(b,newf,))
#         thread.start()
#         bthreads.append(thread)

#      print "waiting threads"
#      for b in range(self.nbeads): bthreads[b].join()
#      print "threads joined in"

   def get_vir(self):
      """Sums the virial of each replica.

      Not the actual system virial, as it has not been divided by either the
      number of beads or the cell volume.

      Returns:
          Virial sum.
      """

      vir = np.zeros((3,3))
      for v in depstrip(self.virs):
         vir += v
      return vir

   def __len__(self):
      """Length function.

      This is called whenever the standard function len(forcebeads) is used.

      Returns:
         The number of beads.
      """

      return self.nbeads

   def __getitem__(self,index):
      """Overwrites standard getting function.

      This is called whenever the standard function forcebeads[index] is used.
      Returns the force on bead index.

      Args:
         index: The index of the replica of the system to be accessed.

      Returns:
         The forces acting on the replica of the system given by the index.
      """

      return self._forces[index]


class Forces(dobject):
   """Class that gathers all the forces together.

   Collects many forcefield instances and parallelizes getting the forces
   in a PIMD environment.

   Attributes:
      natoms: An integer giving the number of atoms.
      nbeads: An integer giving the number of beads.
      nforces: An integer giving the number of ForceBeads objects.
      mforces: A list of all the forcefield objects.
      mbeads: A list of all the beads objects. Some of these may be contracted
         ring polymers, with a smaller number of beads than of the simulation.
      mweights: A list of the weights of all the forcefields.
      mrpc: A list of the objects containing the functions required to
         contract the ring polymers of the different forcefields.

   Depend objects:
      f: An array containing the components of the force. Depends on each
         replica's ufvx list.
      pots: A list containing the potential energy for each system replica.
         Depends on each replica's ufvx list.
      virs: A list containing the virial tensor for each system replica.
         Depends on each replica's ufvx list.
      extras: A list containing the "extra" strings for each replica.
      pot: The sum of the potential energy of the replicas.
      vir: The sum of the virial tensor of the replicas.
   """

   def bind(self, beads, cell, flist):

      self.natoms = beads.natoms
      self.nbeads = beads.nbeads
      self.nforces = len(flist)

      # flist should be a list of tuples containing ( "name", forcebeads)
      self.mforces = []
      self.mbeads = []
      self.mweights = []
      self.mrpc = []

      # a "function factory" to generate functions to automatically update
      #contracted paths
      def make_rpc(rpc, beads):
         return lambda: rpc.b1tob2(depstrip(beads.q))

      # creates new force objects, possibly acting on contracted path
      #representations
      for (ftype, fbeads) in flist:

         # creates an automatically-updated contracted beads object
         newb = fbeads.nbeads
         newforce = fbeads.copy()
         newweight = fbeads.weight

         # if the number of beads for this force component is unspecified,
         #assume full force evaluation
         if newb == 0:
            newb = beads.nbeads
            newforce.nbeads = newb

         newbeads = Beads(beads.natoms, newb)
         newrpc = nm_rescale(beads.nbeads, newb)

         dget(newbeads,"q")._func = make_rpc(newrpc, beads)
         for b in newbeads:
            # must update also indirect access to the beads coordinates
            dget(b,"q")._func = dget(newbeads,"q")._func

         # makes newbeads.q depend from beads.q
         dget(beads,"q").add_dependant(dget(newbeads,"q"))

         #now we create a new forcebeads which is bound to newbeads!
         newforce.bind(newbeads, cell)

         #adds information we will later need to the appropriate lists.
         self.mweights.append(newweight)
         self.mbeads.append(newbeads)
         self.mforces.append(newforce)
         self.mrpc.append(newrpc)

      #now must expose an interface that gives overall forces
      dset(self,"f",
         depend_array(name="f",value=np.zeros((self.nbeads,3*self.natoms)),
            func=self.f_combine,
               dependencies=[dget(ff, "f") for ff in self.mforces] ) )

      # collection of pots and virs from individual ff objects
      dset(self,"pots",
         depend_array(name="pots", value=np.zeros(self.nbeads,float),
            func=self.pot_combine,
               dependencies=[dget(ff, "pots") for ff in self.mforces]) )

      # must take care of the virials!
      dset(self,"virs",
         depend_array(name="virs", value=np.zeros((self.nbeads,3,3),float),
            func=self.vir_combine,
               dependencies=[dget(ff, "virs") for ff in self.mforces]) )

      dset(self,"extras",
         depend_value(name="extras", value=np.zeros(self.nbeads,float),
            func=self.extra_combine,
               dependencies=[dget(ff, "extras") for ff in self.mforces]))

      # total potential and total virial
      dset(self,"pot",
         depend_value(name="pot", func=(lambda: self.pots.sum()),
            dependencies=[dget(self,"pots")]))
      dset(self,"vir",
         depend_array(name="vir", func=self.get_vir, value=np.zeros((3,3)),
            dependencies=[dget(self,"virs")]))

   def run(self):
      """Makes the socket start looking for driver codes.

      Tells the interface code to start the thread that looks for
      connection from the driver codes in a loop. Until this point no
      jobs can be queued.
      """

      for ff in self.mforces:
         ff.run()

   def stop(self):
      """Makes the socket stop looking for driver codes.

      Tells the interface code to stop the thread that looks for
      connection from the driver codes in a loop. After this point no
      jobs can be queued.
      """

      for ff in self.mforces:
         ff.stop()

   def queue(self):
      """Submits all the required force calculations to the forcefields."""

      for ff in self.mforces:
         ff.queue()

   def get_vir(self):
      """Sums the virial of each forcefield.

      Not the actual system virial.

      Returns:
          Virial sum.
      """

      vir = np.zeros((3,3))
      for v in depstrip(self.virs):
         vir += v
      return vir

   def f_combine(self):
      """Obtains the total force vector."""

      self.queue()
      rf = np.zeros((self.nbeads,3*self.natoms),float)
      for k in range(self.nforces):
         # "expand" to the total number of beads the forces from the
         #contracted one
         rf += self.mweights[k]*self.mrpc[k].b2tob1(depstrip(self.mforces[k].f))
      return rf

   def pot_combine(self):
      """Obtains the potential energy for each forcefield."""

      self.queue()
      rp = np.zeros(self.nbeads,float)
      for k in range(self.nforces):
         # "expand" to the total number of beads the potentials from the
         #contracted one
         rp += self.mweights[k]*self.mrpc[k].b2tob1(self.mforces[k].pots)
      return rp

   def extra_combine(self):
      """Obtains the potential energy for each forcefield."""

      self.queue()
      rp = [ "" for b in range(self.nbeads) ]
      for k in range(self.nforces):
         # "expand" to the total number of beads the potentials from the
         #contracted one
         for b in range(self.nbeads):
            rp[b] += self.mforces[k].extras[b]
      return rp

   def vir_combine(self):
      """Obtains the virial tensor for each forcefield."""

      self.queue()
      rp = np.zeros((self.nbeads,3,3),float)
      for k in range(self.nforces):
         virs = depstrip(self.mforces[k].virs)
         # "expand" to the total number of beads the virials from the
         #contracted one, element by element
         for i in range(3):
            for j in range(3):
               rp[:,i,j] += self.mweights[k]*self.mrpc[k].b2tob1(virs[:,i,j])
      return rp
