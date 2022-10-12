"""Contains the classes that deal with the different dynamics required in
different types of ensembles.

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


Holds the algorithms required for normal mode propagators, and the objects to
do the constant temperature and pressure algorithms. Also calculates the
appropriate conserved energy quantity for the ensemble of choice.

Classes:
   Ensemble: Base ensemble class with generic methods and attributes.
   NVEEnsemble: Deals with constant energy dynamics.
   NVTEnsemble: Deals with constant temperature dynamics.
   NPTEnsemble: Deals with constant pressure dynamics.
   ReplayEnsemble: Takes a trajectory, and simply sets the atom positions to
      match it, rather than doing dynamics. In this way new properties can
      be calculated on an old simulation, without having to rerun it from
      scratch.
"""

__all__ = ['Ensemble', 'NVEEnsemble', 'NVTEnsemble', 'NPTEnsemble', 'ReplayEnsemble']

import numpy as np
import time

from ipi.utils.depend import *
from ipi.utils import units
from ipi.utils.softexit import softexit
from ipi.utils.io.io_xyz import read_xyz
from ipi.utils.io.io_pdb import read_pdb
from ipi.utils.io.io_xml import xml_parse_file
from ipi.utils.units import Constants, unit_to_internal
from ipi.inputs.thermostats import InputThermo
from ipi.inputs.barostats import InputBaro
from ipi.engine.thermostats import *
from ipi.engine.barostats import *


class Ensemble(dobject):
   """Base ensemble class.

   Gives the standard methods and attributes needed in all the
   ensemble classes.

   Attributes:
      beads: A beads object giving the atoms positions.
      cell: A cell object giving the system box.
      forces: A forces object giving the virial and the forces acting on
         each bead.
      prng: A random number generator object.
      nm: An object which does the normal modes transformation.
      fixcom: A boolean which decides whether the centre of mass
         motion will be constrained or not.

   Depend objects:
      econs: The conserved energy quantity appropriate to the given
         ensemble. Depends on the various energy terms which make it up,
         which are different depending on the ensemble.
      temp: The system temperature.
      dt: The timestep for the algorithms.
      ntemp: The simulation temperature. Will be nbeads times higher than
         the system temperature as PIMD calculations are done at this
         effective classical temperature.
   """

   def __init__(self, dt, temp, fixcom=False):
      """Initializes Ensemble.

      Args:
         dt: The timestep of the simulation algorithms.
         temp: The temperature.
         fixcom: An optional boolean which decides whether the centre of mass
            motion will be constrained or not. Defaults to False.
      """

      dset(self, "econs", depend_value(name='econs', func=self.get_econs))
      dset(self, "temp",  depend_value(name='temp',  value=temp))
      dset(self, "dt",    depend_value(name='dt',    value=dt))
      self.fixcom = fixcom


   def bind(self, beads, nm, cell, bforce, prng):
      """Binds beads, cell, bforce and prng to the ensemble.

      This takes a beads object, a cell object, a forcefield object and a
      random number generator object and makes them members of the ensemble.
      It also then creates the objects that will hold the data needed in the
      ensemble algorithms and the dependency network. Note that the conserved
      quantity is defined in the init, but as each ensemble has a different
      conserved quantity the dependencies are defined in bind.

      Args:
         beads: The beads object from which the bead positions are taken.
         nm: A normal modes object used to do the normal modes transformation.
         cell: The cell object from which the system box is taken.
         bforce: The forcefield object from which the force and virial are
            taken.
         prng: The random number generator object which controls random number
            generation.
      """

      # store local references to the different bits of the simulation
      self.beads = beads
      self.cell = cell
      self.forces = bforce
      self.prng = prng
      self.nm = nm

      # n times the temperature
      dset(self,"ntemp", depend_value(name='ntemp',func=self.get_ntemp,
         dependencies=[dget(self,"temp")]))

      # dependencies of the conserved quantity
      dget(self,"econs").add_dependency(dget(self.beads, "kin"))
      dget(self,"econs").add_dependency(dget(self.forces, "pot"))
      dget(self,"econs").add_dependency(dget(self.beads, "vpath"))


   def get_ntemp(self):
      """Returns the PI simulation temperature (P times the physical T)."""

      return self.temp*self.beads.nbeads


   def pstep(self):
      """Dummy momenta propagator which does nothing."""

      pass

   def qcstep(self):
      """Dummy centroid position propagator which does nothing."""

      pass

   def step(self):
      """Dummy simulation time step which does nothing."""

      pass

   def get_econs(self):
      """Calculates the conserved energy quantity for constant energy
      ensembles.
      """

      return self.beads.vpath*self.nm.omegan2 + self.nm.kin + self.forces.pot


class NVEEnsemble(Ensemble):
   """Ensemble object for constant energy simulations.

   Has the relevant conserved quantity and normal mode propagator for the
   constant energy ensemble. Note that a temperature of some kind must be
   defined so that the spring potential can be calculated.

   Attributes:
      ptime: The time taken in updating the velocities.
      qtime: The time taken in updating the positions.
      ttime: The time taken in applying the thermostat steps.

   Depend objects:
      econs: Conserved energy quantity. Depends on the bead kinetic and
         potential energy, and the spring potential energy.
   """

   def __init__(self, dt, temp, fixcom=False):
      """Initializes NVEEnsemble.

      Args:
         dt: The simulation timestep.
         temp: The system temperature.
         fixcom: An optional boolean which decides whether the centre of mass
            motion will be constrained or not. Defaults to False.
      """

      super(NVEEnsemble,self).__init__(dt=dt,temp=temp, fixcom=fixcom)

   def rmcom(self):
      """This removes the centre of mass contribution to the kinetic energy.

      Calculates the centre of mass momenta, then removes the mass weighted
      contribution from each atom. If the ensemble defines a thermostat, then
      the contribution to the conserved quantity due to this subtraction is
      added to the thermostat heat energy, as it is assumed that the centre of
      mass motion is due to the thermostat.

      If there is a choice of thermostats, the thermostat
      connected to the centroid is chosen.
      """

      if (self.fixcom):
         pcom = np.zeros(3,float);

         na3 = self.beads.natoms*3
         nb = self.beads.nbeads
         p = depstrip(self.beads.p)
         m = depstrip(self.beads.m3)[:,0:na3:3]
         M = self.beads[0].M

         for i in range(3):
            pcom[i] = p[:,i:na3:3].sum()

         if hasattr(self,"thermostat"):
            if hasattr(self.thermostat, "_thermos"):
               self.thermostat._thermos[0].ethermo += np.dot(pcom,pcom)/(2.0*M*nb)
            else:
               self.thermostat.ethermo += np.dot(pcom,pcom)/(2.0*M*nb)

         # subtracts COM _velocity_
         pcom *= 1.0/(nb*M)
         for i in range(3):
            self.beads.p[:,i:na3:3] -= m*pcom[i]

   def pstep(self):
      """Velocity Verlet momenta propagator."""

      self.beads.p += depstrip(self.forces.f)*(self.dt*0.5)

   def qcstep(self):
      """Velocity Verlet centroid position propagator."""

      self.nm.qnm[0,:] += depstrip(self.nm.pnm)[0,:]/depstrip(self.beads.m3)[0]*self.dt

   def step(self):
      """Does one simulation time step."""

      self.ptime = -time.time()
      self.pstep()
      self.ptime += time.time()

      self.qtime = -time.time()
      self.qcstep()

      self.nm.free_qstep()
      self.qtime += time.time()

      self.ptime -= time.time()
      self.pstep()
      self.ptime += time.time()

      self.ttime = -time.time()
      self.rmcom()
      self.ttime += time.time()


class NVTEnsemble(NVEEnsemble):
   """Ensemble object for constant temperature simulations.

   Has the relevant conserved quantity and normal mode propagator for the
   constant temperature ensemble. Contains a thermostat object containing the
   algorithms to keep the temperature constant.

   Attributes:
      thermostat: A thermostat object to keep the temperature constant.

   Depend objects:
      econs: Conserved energy quantity. Depends on the bead kinetic and
         potential energy, the spring potential energy and the heat
         transferred to the thermostat.
   """

   def __init__(self, dt, temp, thermostat=None, fixcom=False):
      """Initializes NVTEnsemble.

      Args:
         dt: The simulation timestep.
         temp: The system temperature.
         thermostat: A thermostat object to keep the temperature constant.
            Defaults to Thermostat()
         fixcom: An optional boolean which decides whether the centre of mass
            motion will be constrained or not. Defaults to False.
      """

      super(NVTEnsemble,self).__init__(dt=dt,temp=temp, fixcom=fixcom)

      if thermostat is None:
         self.thermostat = Thermostat()
      else:
         self.thermostat = thermostat

   def bind(self, beads, nm, cell, bforce, prng):
      """Binds beads, cell, bforce and prng to the ensemble.

      This takes a beads object, a cell object, a forcefield object and a
      random number generator object and makes them members of the ensemble.
      It also then creates the objects that will hold the data needed in the
      ensemble algorithms and the dependency network. Also note that the
      thermostat timestep and temperature are defined relative to the system
      temperature, and the the thermostat temperature is held at the
      higher simulation temperature, as is appropriate.

      Args:
         beads: The beads object from which the bead positions are taken.
         nm: A normal modes object used to do the normal modes transformation.
         cell: The cell object from which the system box is taken.
         bforce: The forcefield object from which the force and virial are
            taken.
         prng: The random number generator object which controls random number
            generation.
      """

      super(NVTEnsemble,self).bind(beads, nm, cell, bforce, prng)
      fixdof = None
      if self.fixcom:
         fixdof = 3

      # first makes sure that the thermostat has the correct temperature, then proceed with binding it.
      deppipe(self,"ntemp", self.thermostat,"temp")
      deppipe(self,"dt", self.thermostat, "dt")

      #decides whether the thermostat will work in the normal mode or
      #the bead representation.
      if isinstance(self.thermostat,ThermoNMGLE) or isinstance(self.thermostat,ThermoNMGLEG) or isinstance(self.thermostat,ThermoPILE_L) or isinstance(self.thermostat,ThermoPILE_G):
         self.thermostat.bind(nm=self.nm,prng=prng,fixdof=fixdof )
      else:
         self.thermostat.bind(beads=self.beads,prng=prng, fixdof=fixdof)

      dget(self,"econs").add_dependency(dget(self.thermostat, "ethermo"))

   def step(self):
      """Does one simulation time step."""

      self.ttime = -time.time()
      self.thermostat.step()
      self.rmcom()
      self.ttime += time.time()

      self.ptime = -time.time()
      self.pstep()
      self.ptime += time.time()

      self.qtime = -time.time()
      self.qcstep()
      self.nm.free_qstep()
      self.qtime += time.time()

      self.ptime -= time.time()
      self.pstep()
      self.ptime += time.time()

      self.ttime -= time.time()
      self.thermostat.step()
      self.rmcom()
      self.ttime += time.time()

   def get_econs(self):
      """Calculates the conserved energy quantity for constant temperature
      ensemble.
      """

      return NVEEnsemble.get_econs(self) + self.thermostat.ethermo


class NPTEnsemble(NVTEnsemble):
   """Ensemble object for constant pressure simulations.

   Has the relevant conserved quantity and normal mode propagator for the
   constant pressure ensemble. Contains a thermostat object containing the
   algorithms to keep the temperature constant, and a barostat to keep the
   pressure constant.

   Attributes:
      barostat: A barostat object to keep the pressure constant.

   Depend objects:
      econs: Conserved energy quantity. Depends on the bead and cell kinetic
         and potential energy, the spring potential energy, the heat
         transferred to the beads and cell thermostat, the temperature and
         the cell volume.
      pext: External pressure.
   """

   def __init__(self, dt, temp, pext, thermostat=None, barostat=None, fixcom=False):
      """Initializes NPTEnsemble.

      Args:
         dt: The simulation timestep.
         temp: The system temperature.
         pext: The external pressure.
         thermostat: A thermostat object to keep the temperature constant.
            Defaults to Thermostat().
         barostat: A barostat object to keep the pressure constant.
            Defaults to Barostat().
         fixcom: An optional boolean which decides whether the centre of mass
            motion will be constrained or not. Defaults to False.
      """

      super(NPTEnsemble,self).__init__(dt, temp, thermostat, fixcom=fixcom)
      if barostat == None:
         self.barostat = Barostat()
      else:
         self.barostat = barostat

      dset(self,"pext",depend_value(name='pext'))
      if not pext is None:
         self.pext = pext
      else: self.pext = 0.0


   def bind(self, beads, nm, cell, bforce, prng):
      """Binds beads, cell, bforce and prng to the ensemble.

      This takes a beads object, a cell object, a forcefield object and a
      random number generator object and makes them members of the ensemble.
      It also then creates the objects that will hold the data needed in the
      ensemble algorithms and the dependency network. Also note that the cell
      thermostat timesteps and temperatures are defined relative to the system
      temperature, and the the thermostat temperatures are held at the
      higher simulation temperature, as is appropriate.

      Args:
         beads: The beads object from which the bead positions are taken.
         nm: A normal modes object used to do the normal modes transformation.
         cell: The cell object from which the system box is taken.
         bforce: The forcefield object from which the force and virial are
            taken.
         prng: The random number generator object which controls random number
            generation.
      """


      fixdof = None
      if self.fixcom:
         fixdof = 3

      super(NPTEnsemble,self).bind(beads, nm, cell, bforce, prng)
      self.barostat.bind(beads, nm, cell, bforce, prng=prng, fixdof=fixdof)


      deppipe(self,"ntemp", self.barostat, "temp")
      deppipe(self,"dt", self.barostat, "dt")
      deppipe(self,"pext", self.barostat, "pext")
      dget(self,"econs").add_dependency(dget(self.barostat, "ebaro"))

   def get_econs(self):
      """Calculates the conserved energy quantity for the constant pressure
      ensemble.
      """

      return NVTEnsemble.get_econs(self) + self.barostat.ebaro

   def step(self):
      """NPT time step.

      Note that the barostat only propagates the centroid coordinates. If this
      approximation is made a centroid virial pressure and stress estimator can
      be defined, so this gives the best statistical convergence. This is
      allowed as the normal mode propagation is approximately unaffected
      by volume fluctuations as long as the system box is much larger than
      the radius of gyration of the ring polymers.
      """

      self.ttime = -time.time()
      self.thermostat.step()
      self.barostat.thermostat.step()
      self.rmcom()
      self.ttime += time.time()

      self.ptime = -time.time()
      self.barostat.pstep()
      self.ptime += time.time()

      self.qtime = -time.time()
      self.barostat.qcstep()
      self.nm.free_qstep()
      self.qtime += time.time()

      self.ptime -= time.time()
      self.barostat.pstep()
      self.ptime += time.time()

      self.ttime -= time.time()
      self.barostat.thermostat.step()
      self.thermostat.step()
      self.rmcom()
      self.ttime += time.time()


class ReplayEnsemble(Ensemble):
   """Ensemble object that just loads snapshots from an external file in sequence.

   Has the relevant conserved quantity and normal mode propagator for the
   constant energy ensemble. Note that a temperature of some kind must be
   defined so that the spring potential can be calculated.

   Attributes:
      intraj: The input trajectory file.
      ptime: The time taken in updating the velocities.
      qtime: The time taken in updating the positions.
      ttime: The time taken in applying the thermostat steps.

   Depend objects:
      econs: Conserved energy quantity. Depends on the bead kinetic and
         potential energy, and the spring potential energy.
   """

   def __init__(self, dt, temp, fixcom=False, intraj=None):
      """Initializes ReplayEnsemble.

      Args:
         dt: The simulation timestep.
         temp: The system temperature.
         fixcom: An optional boolean which decides whether the centre of mass
            motion will be constrained or not. Defaults to False.
         intraj: The input trajectory file.
      """

      super(ReplayEnsemble,self).__init__(dt=dt,temp=temp,fixcom=fixcom)
      if intraj == None:
         raise ValueError("Must provide an initialized InitFile object to read trajectory from")
      self.intraj = intraj
      if intraj.mode == "manual":
         raise ValueError("Replay can only read from PDB or XYZ files -- or a single frame from a CHK file")
      self.rfile = open(self.intraj.value,"r")

   def step(self):
      """Does one simulation time step."""

      self.ptime = self.ttime = 0
      self.qtime = -time.time()

      try:
         if (self.intraj.mode == "xyz"):
            for b in self.beads:
               myatoms = read_xyz(self.rfile)
               myatoms.q *= unit_to_internal("length",self.intraj.units,1.0)
               b.q[:] = myatoms.q
         elif (self.intraj.mode == "pdb"):
            for b in self.beads:
               myatoms, mycell = read_pdb(self.rfile)
               myatoms.q *= unit_to_internal("length",self.intraj.units,1.0)
               mycell.h  *= unit_to_internal("length",self.intraj.units,1.0)
               b.q[:] = myatoms.q
            self.cell.h[:] = mycell.h
         elif (self.intraj.mode == "chk" or self.intraj.mode == "checkpoint"):
            # reads configuration from a checkpoint file
            xmlchk = xml_parse_file(self.rfile) # Parses the file.

            from ipi.inputs.simulation import InputSimulation
            simchk = InputSimulation()
            simchk.parse(xmlchk.fields[0][1])
            mycell = simchk.cell.fetch()
            mybeads = simchk.beads.fetch()
            self.cell.h[:] = mycell.h
            self.beads.q[:] = mybeads.q
            softexit.trigger(" # Read single checkpoint")
      except EOFError:
         softexit.trigger(" # Finished reading re-run trajectory")

      self.qtime += time.time()


