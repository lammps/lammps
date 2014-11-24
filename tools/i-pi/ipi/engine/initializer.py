"""Contains the classes that are used to initialize data in the simulation.

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


These classes can either be used to restart a simulation with some different
data or used to start a calculation. Any data given in these classes will
overwrite data given elsewhere.

Classes:
   Initializer: Holds the functions that are required to initialize objects in
      the code. Data can be initialized from a file, or according to a
      particular parameter. An example of the former would be initializing
      the configurations from a xyz file, an example of the latter would be
      initializing the velocities according to the physical temperature.
   InitBase: Simple class that reads data from a string or file.
   InitIndexed: The same as init base, but can also optionally hold
      information about which atom or bead to initialize from.

Functions:
   init_xyz: Reads beads data from a xyz file.
   init_pdb: Reads beads and cell data from a pdb file.
   init_chk: Reads beads, cell and thermostat data from a checkpoint file.
   init_beads: Initializes a beads object from an Initializer object.
   init_vector: Initializes a vector from an Initializer object.
   set_vector: Initializes a vector from another vector.
"""

import numpy as np

from ipi.engine.beads import Beads
from ipi.engine.cell import Cell
from ipi.engine.normalmodes import NormalModes
from ipi.engine.ensembles import Ensemble
from ipi.utils.io.io_xyz import read_xyz
from ipi.utils.io.io_pdb import read_pdb
from ipi.utils.io.io_xml import xml_parse_file
from ipi.utils.depend import dobject
from ipi.utils.units import Constants, unit_to_internal
from ipi.utils.nmtransform import nm_rescale
from ipi.utils.messages import verbosity, warning, info

__all__ = ['Initializer', 'InitBase', 'InitIndexed']

class InitBase(dobject):
   """Base class for initializer objects.

   Attributes:
      value: A duck-typed stored value.
      mode: A string that determines how the value is to be interpreted.
      units: A string giving which unit the value is in.
   """

   def __init__(self, value="", mode="", units="", **others):
      """Initializes InitFile.

      Args:
         value: A string which specifies what value to initialize the
            simulation property to.
         mode: A string specifiying what style of initialization should be
            used to read the data.
         units: A string giving which unit the value is in.
      """

      self.value = value
      self.mode = mode
      self.units = units

      for (o, v) in others.items():
         self.__dict__[o] = v


class InitIndexed(InitBase):
   """Class to initialize objects which can be set for a particular bead.

   Attributes:
      index: Which atom to initialize the value of.
      bead: Which bead to initialize the value of.
   """

   def __init__(self, value="", mode="", units="", index=-1, bead=-1):
      """Initializes InitFile.

      Args:
         value: A string which specifies what value to initialize the
            simulation property to.
         mode: A string specifiying what style of initialization should be
            used to read the data.
         units: A string giving which unit the value is in.
         index: Which atom to initialize the value of.
         bead: Which bead to initialize the value of.
      """

      super(InitIndexed,self).__init__(value=value, mode=mode, units=units, index=index, bead=bead)


def init_xyz(filename):
   """Reads an xyz file and returns the data contained in it.

   Args:
      filename: A string giving the name of the xyz file to be read from.

   Returns:
      A list of Atoms objects as read from each frame of the xyz file.
   """

   rfile = open(filename,"r")
   ratoms = []
   while True:
   #while loop, so that more than one configuration can be given
   #so multiple beads can be initialized at once.
      try:
         myatoms = read_xyz(rfile)
      except EOFError:
         break
      ratoms.append(myatoms)
   return ratoms

def init_pdb(filename):
   """Reads an pdb file and returns the data contained in it.

   Args:
      filename: A string giving the name of the pdb file to be read from.

   Returns:
      A list of Atoms objects as read from each frame of the pdb file, and
      a Cell object as read from the final pdb frame.
   """

   rfile = open(filename,"r")
   ratoms = []
   while True:
   #while loop, so that more than one configuration can be given
   #so multiple beads can be initialized at once.
      try:
         myatoms, rcell  = read_pdb(rfile)
      except EOFError:
         break
      ratoms.append(myatoms)
   return ( ratoms, rcell ) # if multiple frames, the last cell is returned

def init_chk(filename):
   """Reads an checkpoint file and returns the data contained in it.

   Args:
      filename: A string giving the name of the checkpoint file to be read from.

   Returns:
      A Beads object, Cell object and Thermostat object as read from the
      checkpoint file.
   """

   # reads configuration from a checkpoint file
   rfile = open(filename,"r")
   xmlchk = xml_parse_file(rfile) # Parses the file.

   from ipi.inputs.simulation import InputSimulation
   simchk = InputSimulation()
   simchk.parse(xmlchk.fields[0][1])
   rcell = simchk.cell.fetch()
   rbeads = simchk.beads.fetch()
   rthermo = simchk.ensemble.thermostat.fetch()

   return (rbeads, rcell, rthermo)

def init_beads(iif, nbeads):
   """A file to initialize a beads object from an appropriate initializer
   object.

   Args:
      iif: An Initializer object which has information on the bead positions.
      nbeads: The number of beads.

   Raises:
      ValueError: If called using an Initializer object with a 'manual' mode.
   """

   mode = iif.mode; value = iif.value
   if mode == "xyz" or mode == "pdb":
      if mode == "xyz": ratoms = init_xyz(value)
      if mode == "pdb": ratoms = init_pdb(value)[0]
      rbeads = Beads(ratoms[0].natoms,len(ratoms))
      for i in range(len(ratoms)): rbeads[i] = ratoms[i]
   elif mode == "chk":
      rbeads = init_chk(value)[0]
   elif mode == "manual":
      raise ValueError("Cannot initialize manually a whole beads object.")

   return rbeads

def init_vector(iif, nbeads, momenta=False):
   """A file to initialize a vector from an appropriate initializer
   object.

   Args:
      iif: An Initializer object specifying the value of a vector.
      nbeads: The number of beads.
      momenta: If bead momenta rather than positions are being initialized
         from a checkpoint file, this is set to True.
   """

   mode = iif.mode; value = iif.value
   if mode == "xyz" or mode == "pdb":
      rq = init_beads(iif, nbeads).q
   elif mode == "chk":
      if momenta: rq = init_beads(iif, nbeads).p
      else:       rq = init_beads(iif, nbeads).q
   elif mode == "manual":
      rq = value

   # determines the size of the input data
   if mode == "manual":
      if iif.bead >= 0: # if there is a bead specifier then we return a single bead slice
         nbeads = 1
      natoms = len(rq)/nbeads/3
      rq.shape = (nbeads,3*natoms)

   return rq

def set_vector(iif, dq, rq):
   """A file to initialize a vector from an another vector.

   If the first dimension is different, i.e. the two vectors correspond
   to a different number of beads, then the ring polymer contraction/expansion
   is used to rescale the original vector to the one used in the simulation,
   as described in the paper T. E. Markland and D. E. Manolopoulos, J. Chem.
   Phys. 129, 024105, (2008).

   Args:
      iif: An Initializer object specifying the value of a vector.
      dq: The vector to be initialized.
      rq: The vector to initialize from.
   """

   (nbeads, natoms) = rq.shape; natoms /= 3
   (dbeads, datoms) = dq.shape; datoms /= 3

   # Check that indices make sense
   if iif.index < 0 and natoms != datoms:
      raise ValueError("Initialization tries to mix up structures with different atom numbers.")
   if iif.index >= datoms:
      raise ValueError("Cannot initialize single atom as atom index %d is larger than the number of atoms" % iif.index)
   if iif.bead >= dbeads:
      raise ValueError("Cannot initialize single bead as bead index %d is larger than the number of beads" % iif.bead)

   if iif.bead < 0:   # we are initializing the path
      res = nm_rescale(nbeads,dbeads)  # path rescaler
      if nbeads != dbeads:
         info(" # Initialize is rescaling from %5d beads to %5d beads" % (nbeads, dbeads), verbosity.low)
      if iif.index < 0:
         dq[:] = res.b1tob2(rq)
      else: # we are initializing a specific atom
         dq[:,3*iif.index:3*(iif.index+1)] = res.b1tob2(rq)
   else:  # we are initializing a specific bead
      if iif.index < 0:
         dq[iif.bead] = rq
      else:
         dq[iif.bead,3*iif.index:3*(iif.index+1)] = rq

class Initializer(dobject):
   """Class that deals with the initialization of data.

   This can either be used to initialize the atom positions and the cell data
   from a file, or to initialize them from a beads, atoms or cell object.

   Currently, we use a ring polymer contraction scheme to create a new beads
   object from one given in initialize if they have different numbers of beads,
   as described in the paper T. E. Markland and D. E. Manolopoulos, J. Chem.
   Phys. 129, 024105, (2008). If the new beads object has more beads than
   the beads object it was initialized from, we set the higher ring polymer
   normal modes to zero.

   Attributes:
      queue: A list of things to initialize. Each member of the list is a tuple
         of the form ('type', 'object'), where 'type' specifies what kind of
         initialization is being done, and 'object' gives the data to
         initialize it from.
   """

   def __init__(self, nbeads=0, queue=None):
      """Initializes Initializer.

      Arguments:
         nbeads: The number of beads that we need in the simulation. Not
            necessarily the same as the number of beads of the objects we are
            initializing the data from.
         queue: A list of things to initialize. Each member of the list is a
            tuple of the form ('type', 'object'), where 'type' specifies what
            kind of initialization is being done, and 'object' gives the data to
            initialize it from.
      """

      self.nbeads = nbeads

      if queue is None:
         self.queue = []
      else:
         self.queue = queue

   def init_stage1(self, simul):
      """Initializes the simulation -- first stage.

      Takes a simulation object, and uses all the data in the initialization
      queue to fill up the beads and cell data needed to run the simulation.

      Args:
         simul: A simulation object to be initialized.

      Raises:
         ValueError: Raised if there is a problem with the initialization,
            if something that should have been has not been, or if the objects
            that have been specified are not compatible with each other.
      """

      if simul.beads.nbeads == 0:
         fpos = fmom = fmass = flab = fcell = False   # we don't have an explicitly defined beads object yet
      else:
         fpos = fmom = fmass = flab = fcell = True

      for (k,v) in self.queue:
         info(" # Initializer (stage 1) parsing " + str(k) + " object.", verbosity.high)

         if k == "cell":
            if fcell :
               warning("Overwriting previous cell parameters", verbosity.medium)
            if v.mode == "pdb":
               rh = init_pdb(v.value)[1].h
            elif v.mode == "chk":
               rh = init_chk(v.value)[1].h
            else:
               rh = v.value.reshape((3,3))
            rh *= unit_to_internal("length",v.units,1.0)

            simul.cell.h = rh
            if simul.cell.V == 0.0:
               ValueError("Cell provided has zero volume")

            fcell = True
         elif k == "masses":
            if simul.beads.nbeads == 0:
               raise ValueError("Cannot initialize the masses before the size of the system is known")
            if fmass:
               warning("Overwriting previous atomic masses", verbosity.medium)
            if v.mode == "manual":
               rm = v.value
            else:
               rm = init_beads(v, self.nbeads).m
            rm *= unit_to_internal("mass",v.units,1.0)

            if v.bead < 0:   # we are initializing the path
               if (fmom and fmass):
                  warning("Rescaling momenta to make up for changed mass", verbosity.medium)
                  simul.beads.p /= simul.beads.sm3   # go to mass-scaled momenta, that are mass-invariant
               if v.index < 0:
                  simul.beads.m = rm
               else: # we are initializing a specific atom
                  simul.beads.m[v.index:v.index+1] = rm
               if (fmom and fmass):  # finishes correcting the momenta
                  simul.beads.p *= simul.beads.sm3  # back to normal momenta
            else:
               raise ValueError("Cannot change the mass of a single bead")
            fmass = True

         elif k == "labels":
            if simul.beads.nbeads == 0:
               raise ValueError("Cannot initialize the labels before the size of the system is known")
            if flab:
               warning("Overwriting previous atomic labels", verbosity.medium)
            if v.mode == "manual":
               rn = v.value
            else:
               rn = init_beads(v, self.nbeads).names

            if v.bead < 0:   # we are initializing the path
               if v.index < 0:
                  simul.beads.names = rn
               else: # we are initializing a specific atom
                  simul.beads.names[v.index:v.index+1] = rn
            else:
               raise ValueError("Cannot change the label of a single bead")
            flab = True

         elif k == "positions":
            if fpos:
               warning("Overwriting previous atomic positions", verbosity.medium)
            # read the atomic positions as a vector
            rq = init_vector(v, self.nbeads)
            rq *= unit_to_internal("length",v.units,1.0)
            (nbeads, natoms) = rq.shape;   natoms /= 3

            # check if we must initialize the simulation beads
            if simul.beads.nbeads == 0:
               if v.index >= 0:
                  raise ValueError("Cannot initialize single atoms before the size of the system is known")
               simul.beads.resize(natoms,self.nbeads)

            set_vector(v, simul.beads.q, rq)
            fpos = True

         elif (k == "velocities" or k == "momenta") and v.mode == "thermal" :   # intercept here thermal initialization, so we don't need to check further down
            if fmom:
               warning("Overwriting previous atomic momenta", verbosity.medium)
            if simul.beads.natoms == 0:
               raise ValueError("Cannot initialize momenta before the size of the system is known.")
            if not fmass:
               raise ValueError("Trying to resample velocities before having masses.")

            rtemp = v.value * unit_to_internal("temperature",v.units,1.0)
            if rtemp <= 0:
               warning("Using the simulation temperature to resample velocities", verbosity.low)
               rtemp = simul.ensemble.temp
            else:
               info(" # Resampling velocities at temperature %s %s" % (v.value, v.units), verbosity.low)

            # pull together a mock initialization to get NM masses right
            #without too much code duplication
            if v.bead >= 0:
               raise ValueError("Cannot thermalize a single bead")
            if v.index >= 0:
               rnatoms = 1
            else:
               rnatoms = simul.beads.natoms
            rbeads = Beads(rnatoms, simul.beads.nbeads)
            if v.index < 0:
               rbeads.m[:] = simul.beads.m
            else:
               rbeads.m[:] = simul.beads.m[v.index]
            rnm = NormalModes(mode=simul.nm.mode, transform_method=simul.nm.transform_method, freqs=simul.nm.nm_freqs)
            rens = Ensemble(dt=simul.ensemble.dt, temp=simul.ensemble.temp)
            rnm.bind(rbeads,rens)
            # then we exploit the sync magic to do a complicated initialization
            # in the NM representation
            # with (possibly) shifted-frequencies NM
            rnm.pnm = simul.prng.gvec((rbeads.nbeads,3*rbeads.natoms))*np.sqrt(rnm.dynm3)*np.sqrt(rbeads.nbeads*rtemp*Constants.kb)

            if v.index < 0:
               simul.beads.p = rbeads.p
            else:
               simul.beads.p[:,3*v.index:3*(v.index+1)] = rbeads.p
            fmom = True

         elif k == "momenta":
            if fmom:
               warning("Overwriting previous atomic momenta", verbosity.medium)
            # read the atomic momenta as a vector
            rp = init_vector(v, self.nbeads, momenta = True)
            rp *= unit_to_internal("momentum",v.units,1.0)
            (nbeads, natoms) = rp.shape;   natoms /= 3

            # checks if we must initialize the simulation beads
            if simul.beads.nbeads == 0:
               if v.index >= 0 :
                  raise ValueError("Cannot initialize single atoms before the size of the system is known")
               simul.beads.resize(natoms,self.nbeads)

            rp *= np.sqrt(self.nbeads/nbeads)
            set_vector(v, simul.beads.p, rp)
            fmom = True

         elif k == "velocities":
            if fmom:
               warning("Overwriting previous atomic momenta", verbosity.medium)
            # read the atomic velocities as a vector
            rv = init_vector(v, self.nbeads)
            rv *= unit_to_internal("velocity",v.units,1.0)
            (nbeads, natoms) = rv.shape;   natoms /= 3

            # checks if we must initialize the simulation beads
            if simul.beads.nbeads == 0 or not fmass:
               ValueError("Cannot initialize velocities before the masses of the atoms are known")
               simul.beads.resize(natoms,self.nbeads)

            warning("Initializing from velocities uses the previously defined masses -- not the masses inferred from the file -- to build momenta", verbosity.low)
            if v.index >= 0:
               rv *= simul.beads.m[v.index]
            elif v.bead >= 0:
               rv *= simul.beads.m3[0]
            else:
               rv *= simul.beads.m3
            rv *= np.sqrt(self.nbeads/nbeads)
            set_vector(v, simul.beads.p, rv)
            fmom = True
         elif k == "thermostat": pass   # thermostats must be initialised in a second stage

      if simul.beads.natoms == 0:
         raise ValueError("Initializer could not initialize the atomic positions")
      if simul.cell.V == 0:
         raise ValueError("Initializer could not initialize the cell")
      for i in range(simul.beads.natoms):
         if simul.beads.m[i] <= 0:
            raise ValueError("Initializer could not initialize the masses")
         if simul.beads.names[i] == "":
            raise ValueError("Initializer could not initialize the atom labels")
      if not fmom:
         warning("Momenta not specified in initialize. Will start with zero velocity if they are not specified in beads.", verbosity.low)

   def init_stage2(self, simul):
      """Initializes the simulation -- second stage.

      Takes a simulation object which has been fully generated,
      and restarts additional information such as the thermostat internal state.

      Args:
         simul: A simulation object to be initialized.

      Raises:
         ValueError: Raised if there is a problem with the initialization,
            if something that should have been has not been, or if the objects
            that have been specified are not compatible with each other.
      """

      for (k,v) in self.queue:
         info(" # Initializer (stage 2) parsing " + str(k) + " object.", verbosity.high)

         if k == "gle":
            # read thermostat parameters from file
            if not ( hasattr(simul.ensemble, "thermostat") ):
               raise ValueError("Ensemble does not have a thermostat to initialize")
            if not ( hasattr(simul.ensemble.thermostat, "s") ):
               raise ValueError("There is nothing to initialize in non-GLE thermostats")
            ssimul = simul.ensemble.thermostat.s
            if v.mode == "manual":
               sinput = v.value.copy()
               if (sinput.size() != ssimul.size() ):
                  raise ValueError("Size mismatch in thermostat initialization data")
               sinput.shape = ssimul.shape
            elif v.mode == "chk":
               rthermo = init_chk(v.value)[2]
               if not hasattr(rthermo,"s"):
                  raise ValueError("Checkpoint file does not contain usable thermostat data")
               sinput = rthermo.s.copy()
               if sinput.shape != ssimul.shape :
                  raise ValueError("Shape mismatch in thermostat initialization data")

            # if all the preliminary checks are good, we can initialize the s's
            ssimul[:] = sinput
