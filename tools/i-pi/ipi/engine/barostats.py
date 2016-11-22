"""Contains the classes that deal with constant pressure dynamics.

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


Contains the algorithms which propagate the position and momenta steps in the
constant pressure ensemble. Holds the properties directly related to
these ensembles, such as the internal and external pressure and stress.

Classes:
   Barostat: Base barostat class with the generic methods and attributes.
   BaroBZP: Generates dynamics with a stochastic barostat -- see
            Ceriotti, More, Manolopoulos, Comp. Phys. Comm. 2013 for
            implementation details.
"""

# NB: this file also contains a 'BaroMHT' class, that follows more closely the
# Martyna, Hughes, Tuckerman implementation of a PIMD barostat. However it is so
# close to the BZP implementation that we disabled it for the sake of simplicity
# BaroMHT: Generates dynamics according to the method of G. Martyna, A.
# Hughes and M. Tuckerman, J. Chem. Phys., 110, 3275.

__all__ = ['Barostat', 'BaroBZP']

import numpy as np
from ipi.utils.depend import *
from ipi.utils.units import *
from ipi.utils.mathtools import eigensystem_ut3x3, invert_ut3x3, exp_ut3x3, det_ut3x3
from ipi.inputs.thermostats import InputThermo
from ipi.engine.thermostats import Thermostat

class Barostat(dobject):
   """Base barostat class.

   Gives the standard methods and attributes needed in all the barostat classes.

   Attributes:
      beads: A beads object giving the atoms positions
      cell: A cell object giving the system box.
      forces: A forces object giving the virial and the forces acting on
         each bead.
      nm: An object to do the normal mode transformation.
      thermostat: A thermostat coupled to the barostat degrees of freedom.
      mdof: The number of atomic degrees of freedom

   Depend objects:
      dt: The time step used in the algorithms. Depends on the simulation dt.
      temp: The (classical) simulation temperature. Higher than the physical
         temperature by a factor of the number of beads.
      tau: The timescale associated with the piston
      pext: The external pressure
      ebaro: The conserved quantity associated with the barostat.
      pot: The potential energy associated with the barostat.
      kstress: The system kinetic stress tensor.
      stress: The system stress tensor.
      press: The system pressure.
   """

   def __init__(self, dt=None, temp=None, pext=None, tau=None, ebaro=None, thermostat=None):
      """Initialises base barostat class.

      Note that the external stress and the external pressure are synchronized.
      This makes most sense going from the stress to the pressure, but if you
      must go in the other direction the stress is assumed to be isotropic.

      Args:
         dt: Optional float giving the time step for the algorithms. Defaults
            to the simulation dt.
         temp: Optional float giving the temperature for the thermostat.
            Defaults to the simulation temp.
         pext: Optional float giving the external pressure.
         tau: Optional float giving the time scale associated with the barostat.
         ebaro: Optional float giving the conserved quantity already stored
            in the barostat initially. Used on restart.
         thermostat: The thermostat connected to the barostat degree of freedom.
      """

      dset(self,"dt",depend_value(name='dt'))
      if not dt is None:
         self.dt = dt
      else: self.dt = 1.0

      dset(self, "temp", depend_value(name="temp"))
      if not temp is None:
         self.temp = temp
      else: self.temp = 1.0

      dset(self,"tau",depend_value(name='tau'))
      if not tau is None:
         self.tau = tau
      else: self.tau = 1.0

      dset(self,"pext",depend_value(name='pext'))
      if not pext is None:
         self.pext = pext
      else: self.pext = 0.0

      dset(self,"ebaro",depend_value(name='ebaro'))
      if not ebaro is None:
         self.ebaro = ebaro
      else: self.ebaro = 0.0

      if thermostat is None:
         thermostat = Thermostat()
      self.thermostat = thermostat

      # pipes timestep and temperature to the thermostat
      deppipe(self,"dt", self.thermostat, "dt")
      deppipe(self, "temp", self.thermostat,"temp")


   def bind(self, beads, nm, cell, forces, prng=None, fixdof=None):
      """Binds beads, cell and forces to the barostat.

      This takes a beads object, a cell object and a forcefield object and
      makes them members of the barostat. It also then creates the objects that
      will hold the data needed in the barostat algorithms and the dependency
      network.

      Args:
         beads: The beads object from which the bead positions are taken.
         nm: The normal modes propagator object
         cell: The cell object from which the system box is taken.
         forces: The forcefield object from which the force and virial are
            taken.
         prng: The parent PRNG to bind the thermostat to
         fixdof: The number of blocked degrees of freedom.
      """

      self.beads = beads
      self.cell = cell
      self.forces = forces
      self.nm = nm

      dset(self,"pot",
         depend_value(name='pot', func=self.get_pot,
            dependencies=[ dget(cell,"V"), dget(self,"pext") ]))
      dset(self,"kstress",
         depend_value(name='kstress', func=self.get_kstress,
            dependencies=[ dget(beads,"q"), dget(beads,"qc"), dget(beads,"pc"), dget(forces,"f") ]))
      dset(self,"stress",
         depend_value(name='stress', func=self.get_stress,
            dependencies=[ dget(self,"kstress"), dget(cell,"V"), dget(forces,"vir") ]))
      dset(self,"press",
         depend_value(name='press', func=self.get_press,
            dependencies=[ dget(self,"stress") ]))

      if fixdof is None:
         self.mdof = float(self.beads.natoms)*3.0
      else:
         self.mdof = float(self.beads.natoms)*3.0 - float(fixdof)

   def get_pot(self):
      """Calculates the elastic strain energy of the cell."""

      # NOTE: since there are nbeads replicas of the unit cell, the enthalpy contains a nbeads factor
      return self.cell.V*self.pext*self.beads.nbeads

   def get_kstress(self):
      """Calculates the quantum centroid virial kinetic stress tensor
      estimator.
      """

      kst = np.zeros((3,3),float)
      q = depstrip(self.beads.q)
      qc = depstrip(self.beads.qc)
      pc = depstrip(self.beads.pc)
      m = depstrip(self.beads.m)
      na3 = 3*self.beads.natoms
      fall = depstrip(self.forces.f)

      for b in range(self.beads.nbeads):
         for i in range(3):
            for j in range(i,3):
               kst[i,j] -= np.dot(q[b,i:na3:3] - qc[i:na3:3],
                  fall[b,j:na3:3])

      # NOTE: In order to have a well-defined conserved quantity, the Nf kT term in the
      # diagonal stress estimator must be taken from the centroid kinetic energy.
      for i in range(3):
         kst[i,i] += np.dot(pc[i:na3:3],pc[i:na3:3]/m) *self.beads.nbeads

      return kst

   def get_stress(self):
      """Calculates the internal stress tensor."""

      return (self.kstress + self.forces.vir)/self.cell.V

   def get_press(self):
      """Calculates the internal pressure."""

      return np.trace(self.stress)/3.0

   def pstep(self):
      """Dummy momenta propagator step."""

      pass

   def qcstep(self):
      """Dummy centroid position propagator step."""

      pass


class BaroBZP(Barostat):
   """Bussi-Zykova-Parrinello barostat class.

   Just extends the standard class adding finite-dt propagators for the barostat
   velocities, positions, piston.

   Depend objects:
      p: The momentum associated with the volume degree of freedom.
      m: The mass associated with the volume degree of freedom.
   """

   def __init__(self, dt=None, temp=None, pext=None, tau=None, ebaro=None, thermostat=None, p=None):
      """Initializes BZP barostat.

      Args:
         dt: Optional float giving the time step for the algorithms. Defaults
            to the simulation dt.
         temp: Optional float giving the temperature for the thermostat.
            Defaults to the simulation temp.
         pext: Optional float giving the external pressure.
         tau: Optional float giving the time scale associated with the barostat.
         ebaro: Optional float giving the conserved quantity already stored
            in the barostat initially. Used on restart.
         thermostat: The thermostat connected to the barostat degree of freedom.
         p: Optional initial volume conjugate momentum. Defaults to 0.
      """


      super(BaroBZP, self).__init__(dt, temp, pext, tau, ebaro, thermostat)

      dset(self,"p", depend_array(name='p', value=np.atleast_1d(0.0)))

      if not p is None:
         self.p = np.asarray([p])
      else:
         self.p = 0.0

   def bind(self, beads, nm, cell, forces, prng=None, fixdof=None):
      """Binds beads, cell and forces to the barostat.

      This takes a beads object, a cell object and a forcefield object and
      makes them members of the barostat. It also then creates the objects that
      will hold the data needed in the barostat algorithms and the dependency
      network.

      Args:
         beads: The beads object from which the bead positions are taken.
         nm: The normal modes propagator object
         cell: The cell object from which the system box is taken.
         forces: The forcefield object from which the force and virial are
            taken.
         prng: The parent PRNG to bind the thermostat to
         fixdof: The number of blocked degrees of freedom.
      """

      super(BaroBZP, self).bind(beads, nm, cell, forces, prng, fixdof)

      # obtain the thermostat mass from the given time constant
      # note that the barostat temperature is nbeads times the physical T
      dset(self,"m", depend_array(name='m', value=np.atleast_1d(0.0),
         func=(lambda:np.asarray([self.tau**2*3*self.beads.natoms*Constants.kb*self.temp])),
            dependencies=[ dget(self,"tau"), dget(self,"temp") ] ))

      # binds the thermostat to the piston degrees of freedom
      self.thermostat.bind(pm=[ self.p, self.m ], prng=prng)

      dset(self,"kin",depend_value(name='kin',
         func=(lambda:0.5*self.p[0]**2/self.m[0]),
            dependencies= [dget(self,"p"), dget(self,"m")] ) )

      # the barostat energy must be computed from bits & pieces (overwrite the default)
      dset(self, "ebaro", depend_value(name='ebaro', func=self.get_ebaro,
         dependencies=[ dget(self, "kin"), dget(self, "pot"),
            dget(self.cell, "V"), dget(self, "temp"),
               dget(self.thermostat,"ethermo")] ))

   def get_ebaro(self):
      """Calculates the barostat conserved quantity."""

      return self.thermostat.ethermo + self.kin + self.pot - np.log(self.cell.V)*Constants.kb*self.temp


   def pstep(self):
      """Propagates the momenta for half a time step."""

      dthalf = self.dt*0.5
      dthalf2 = dthalf**2
      dthalf3 = dthalf**3/3.0

      # This differs from the BZP thermostat in that it uses just one kT in the propagator.
      # This leads to an ensemble equaivalent to Martyna-Hughes-Tuckermann for both fixed and moving COM
      # Anyway, it is a small correction so whatever.
      self.p += dthalf*3.0*( self.cell.V* ( self.press - self.beads.nbeads*self.pext ) +
                Constants.kb*self.temp )

      fc = np.sum(depstrip(self.forces.f),0)/self.beads.nbeads
      m = depstrip(self.beads.m3)[0]
      pc = depstrip(self.beads.pc)

      # I am not 100% sure, but these higher-order terms come from integrating the pressure virial term,
      # so they should need to be multiplied by nbeads to be consistent with the equations of motion in the PI context
      # again, these are tiny tiny terms so whatever.
      self.p += (dthalf2*np.dot(pc,fc/m) + dthalf3*np.dot(fc,fc/m)) * self.beads.nbeads

      self.beads.p += depstrip(self.forces.f)*dthalf

   def qcstep(self):
      """Propagates the centroid position and momentum and the volume."""

      v = self.p[0]/self.m[0]
      expq, expp = (np.exp(v*self.dt), np.exp(-v*self.dt))

      m = depstrip(self.beads.m3)[0]

      self.nm.qnm[0,:] *= expq
      self.nm.qnm[0,:] += ((expq-expp)/(2.0*v))* (depstrip(self.nm.pnm)[0,:]/m)
      self.nm.pnm[0,:] *= expp

      self.cell.h *= expq


class BaroMHT(Barostat):
   """Martyna-Hughes-Tuckerman barostat class.

   Just extends the standard class adding finite-dt propagators for the barostat
   velocities, positions, piston.

   Depend objects:
      p: The momentum associated with the volume degree of freedom.
      m: The mass associated with the volume degree of freedom.
   """

   def __init__(self, dt=None, temp=None, pext=None, tau=None, ebaro=None, thermostat=None, p=None):
      """Initializes MHT barostat.

      Args:
         dt: Optional float giving the time step for the algorithms. Defaults
            to the simulation dt.
         temp: Optional float giving the temperature for the thermostat.
            Defaults to the simulation temp.
         pext: Optional float giving the external pressure.
         tau: Optional float giving the time scale associated with the barostat.
         ebaro: Optional float giving the conserved quantity already stored
            in the barostat initially. Used on restart.
         thermostat: The thermostat connected to the barostat degree of freedom.
         p: Optional initial volume conjugate momentum. Defaults to 0.
      """

      super(BaroMHT, self).__init__(dt, temp, pext, tau, ebaro, thermostat)

      dset(self,"p", depend_array(name='p', value=np.atleast_1d(0.0)))

      if not p is None:
         self.p = np.asarray([p])
      else:
         self.p = 0.0

   def bind(self, beads, nm, cell, forces, prng=None, fixdof=None):
      """Binds beads, cell and forces to the barostat.

      This takes a beads object, a cell object and a forcefield object and
      makes them members of the barostat. It also then creates the objects that
      will hold the data needed in the barostat algorithms and the dependency
      network.

      Args:
         beads: The beads object from which the bead positions are taken.
         nm: The normal modes propagator object
         cell: The cell object from which the system box is taken.
         forces: The forcefield object from which the force and virial are
            taken.
         prng: The parent PRNG to bind the thermostat to
         fixdof: The number of blocked degrees of freedom.
      """

      super(BaroMHT, self).bind(beads, nm, cell, forces, prng, fixdof)

      # obtain the thermostat mass from the given time constant
      # note that the barostat temperature is nbeads times the physical T
      dset(self,"m", depend_array(name='m', value=np.atleast_1d(0.0),
         func=(lambda:np.asarray([self.tau**2*3*self.beads.natoms*Constants.kb*self.temp])),
            dependencies=[ dget(self,"tau"), dget(self,"temp") ] ))

      # binds the thermostat to the piston degrees of freedom
      self.thermostat.bind(pm=[ self.p, self.m ], prng=prng)

      dset(self,"kin",depend_value(name='kin',
         func=(lambda:0.5*self.p[0]**2/self.m[0]),
            dependencies=[dget(self,"p"), dget(self,"m")] ) )

      # the barostat energy must be computed from bits & pieces (overwrite the default)
      dset(self, "ebaro", depend_value(name='ebaro', func=self.get_ebaro,
         dependencies=[ dget(self, "kin"), dget(self, "pot"),
            dget(self.cell, "V"), dget(self, "temp"),
               dget(self.thermostat,"ethermo")]))

   def get_ebaro(self):
      """Calculates the barostat conserved quantity."""

      return self.thermostat.ethermo + self.kin + self.pot

   def pstep(self):
      """Propagates the momenta for half a time step."""

      dthalf = self.dt*0.5
      dthalf2 = dthalf**2
      dthalf3 = dthalf**3/3.0

      fc = np.sum(depstrip(self.forces.f),0)/float(self.beads.nbeads)
      m = depstrip(self.beads.m3)[0]
      pc = depstrip(self.beads.pc)

      self.p += dthalf*3.0*( self.cell.V* ( self.press - self.beads.nbeads*self.pext ) +
                float(self.beads.nbeads)/self.mdof*np.dot(pc,pc/m) )

      self.beads.p += depstrip(self.forces.f)*dthalf

   def qcstep(self):
      """Propagates the centroid position and momentum and the volume."""

      v = self.p[0]/self.m[0]
      adof = (1 + 3.0/self.mdof)
      expq, expp = (np.exp(v*self.dt), np.exp( -v*self.dt * adof  ) )

      m = depstrip(self.beads.m3)[0]

      self.nm.qnm[0,:] *= expq
      self.nm.qnm[0,:] += ((expq-expp)/(v*(1+adof)) *
                    (depstrip(self.nm.pnm)[0,:])/m)
      self.nm.pnm[0,:] *= expp

      self.cell.h *= expq
