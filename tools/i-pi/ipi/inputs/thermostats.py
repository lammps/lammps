"""Deals with creating the thermostats class.

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


Chooses between the different possible thermostat options and creates the
appropriate thermostat object, with suitable parameters.

Classes:
   InputThermo: Deals with creating the thermostat object from a file, and
      writing the checkpoints.
"""

__all__ = ['InputThermo']

import numpy as np
from ipi.utils.depend   import *
from ipi.utils.inputvalue  import *
from ipi.engine.thermostats import *

class InputThermo(Input):
   """Thermostat input class.

   Handles generating the appropriate thermostat class from the xml input file,
   and generating the xml checkpoiunt tags and data from an instance of the
   object.

   Attributes:
      mode: An optional string giving the type of the thermostat used. Defaults
         to 'langevin'.

   Fields:
      ethermo: An optional float giving the amount of heat energy transferred
         to the bath. Defaults to 0.0.
      tau: An optional float giving the damping time scale. Defaults to 1.0.
      pile_scale: Scaling for the PILE damping relative to the critical damping.
      A: An optional array of floats giving the drift matrix. Defaults to 0.0.
      C: An optional array of floats giving the static covariance matrix.
         Defaults to 0.0.
      s: An optional array of floats giving the additional momentum-scaled
         momenta in GLE. Defaults to 0.0.
   """

   attribs = { "mode": (InputAttribute, { "dtype"   : str,
                                      "options" : [ "", "langevin", "svr", "pile_l", "pile_g", "gle", "nm_gle", "nm_gle_g" ],
                                      "help"    : "The style of thermostatting. 'langevin' specifies a white noise langevin equation to be attached to the cartesian representation of the momenta. 'svr' attaches a velocity rescaling thermostat to the cartesian representation of the momenta. Both 'pile_l' and 'pile_g' attaches a white noise langevin thermostat to the normal mode representation, with 'pile_l' attaching a local langevin thermostat to the centroid mode and 'pile_g' instead attaching a global velocity rescaling thermostat. 'gle' attaches a colored noise langevin thermostat to the cartesian representation of the momenta, 'nm_gle' attaches a colored noise langevin thermostat to the normal mode representation of the momenta and a langevin thermostat to the centroid and 'nm_gle_g' attaches a gle thermostat to the normal modes and a svr thermostat to the centroid."
                                         }) }
   fields = { "ethermo" : (InputValue, {  "dtype"     : float,
                                          "default"   : 0.0,
                                          "help"      : "The initial value of the thermostat energy. Used when the simulation is restarted to guarantee continuity of the conserved quantity.",
                                          "dimension" : "energy" }),
            "tau" : (InputValue, {  "dtype"     : float,
                                    "default"   : 0.0,
                                    "help"      : "The friction coefficient for white noise thermostats.",
                                    "dimension" : "time" }),
            "pile_scale" : (InputValue, { "dtype" : float,
                                    "default"   : 1.0,
                                    "help"      : "Scaling for the PILE damping relative to the critical damping."} ),
            "A" : (InputArray, {    "dtype"     : float,
                                    "default"   : input_default(factory=np.zeros, args = (0,)),
                                    "help"      : "The friction matrix for GLE thermostats.",
                                    "dimension" : "frequency" }),
            "C" : (InputArray, {    "dtype"     : float,
                                    "default"   : input_default(factory=np.zeros, args = (0,)),
                                    "help"      : "The covariance matrix for GLE thermostats.",
                                    "dimension" : "temperature" }),
            "s" : (InputArray, {    "dtype"     : float,
                                    "default"   : input_default(factory=np.zeros, args = (0,)),
                                    "help"      : "Input values for the additional momenta in GLE.",
                                    "dimension" : "ms-momentum" })
             }

   default_help = "Simulates an external heat bath to keep the velocity distribution at the correct temperature."
   default_label = "THERMOSTATS"

   def store(self, thermo):
      """Takes a thermostat instance and stores a minimal representation of it.

      Args:
         thermo: A thermostat object.

      Raises:
         TypeError: Raised if the thermostat is not a recognized type.
      """

      super(InputThermo,self).store(thermo)
      if type(thermo) is ThermoLangevin:
         self.mode.store("langevin")
         self.tau.store(thermo.tau)
      elif type(thermo) is ThermoSVR:
         self.mode.store("svr")
         self.tau.store(thermo.tau)
      elif type(thermo) is ThermoPILE_L:
         self.mode.store("pile_l")
         self.tau.store(thermo.tau)
         self.pile_scale.store(thermo.pilescale)
      elif type(thermo) is ThermoPILE_G:
         self.mode.store("pile_g")
         self.tau.store(thermo.tau)
         self.pile_scale.store(thermo.pilescale)
      elif type(thermo) is ThermoGLE:
         self.mode.store("gle")
         self.A.store(thermo.A)
         if dget(thermo,"C")._func is None:
            self.C.store(thermo.C)
         self.s.store(thermo.s)
      elif type(thermo) is ThermoNMGLE:
         self.mode.store("nm_gle")
         self.A.store(thermo.A)
         if dget(thermo,"C")._func is None:
            self.C.store(thermo.C)
         self.s.store(thermo.s)
      elif type(thermo) is ThermoNMGLEG:
         self.mode.store("nm_gle_g")
         self.A.store(thermo.A)
         self.tau.store(thermo.tau)
         if dget(thermo,"C")._func is None:
            self.C.store(thermo.C)
         self.s.store(thermo.s)
      elif type(thermo) is Thermostat:
         self.mode.store("")
      else:
         raise TypeError("Unknown thermostat mode " + type(thermo).__name__)
      self.ethermo.store(thermo.ethermo)

   def fetch(self):
      """Creates a thermostat object.

      Returns:
         A thermostat object of the appropriate type and with the appropriate
         parameters given the attributes of the InputThermo object.

      Raises:
         TypeError: Raised if the thermostat type is not a recognized option.
      """

      super(InputThermo,self).fetch()
      if self.mode.fetch() == "langevin":
         thermo = ThermoLangevin(tau=self.tau.fetch())
      elif self.mode.fetch() == "svr":
         thermo = ThermoSVR(tau=self.tau.fetch())
      elif self.mode.fetch() == "pile_l":
         thermo = ThermoPILE_L(tau=self.tau.fetch(), scale=self.pile_scale.fetch())
      elif self.mode.fetch() == "pile_g":
         thermo = ThermoPILE_G(tau=self.tau.fetch(), scale=self.pile_scale.fetch())
      elif self.mode.fetch() == "gle":
         rC = self.C.fetch()
         if len(rC) == 0:
            rC = None
         thermo = ThermoGLE(A=self.A.fetch(),C=rC)
         thermo.s = self.s.fetch()
      elif self.mode.fetch() == "nm_gle":
         rC = self.C.fetch()
         if len(rC) == 0:
            rC = None
         thermo = ThermoNMGLE(A=self.A.fetch(),C=rC)
         thermo.s = self.s.fetch()
      elif self.mode.fetch() == "nm_gle_g":
         rC = self.C.fetch()
         if len(rC) == 0:
            rC = None
         thermo = ThermoNMGLEG(A=self.A.fetch(),C=rC, tau=self.tau.fetch())
         thermo.s = self.s.fetch()
      elif self.mode.fetch() == "" :
         thermo=Thermostat()
      else:
         raise TypeError("Invalid thermostat mode " + self.mode.fetch())

      thermo.ethermo = self.ethermo.fetch()

      return thermo

   def check(self):
      """Checks that the parameter arrays represents a valid thermostat."""

      super(InputThermo,self).check()

      if self.mode.fetch() in ["langevin", "svr", "pile_l", "pile_g", "nm_gle_g"]:
         if self.tau.fetch() <= 0:
            raise ValueError("The thermostat friction coefficient must be set to a positive value")
      if self.mode.fetch() in ["gle", "nm_gle", "nm_gle_g"]:
         pass  # PERHAPS DO CHECKS THAT MATRICES SATISFY REASONABLE CONDITIONS (POSITIVE-DEFINITENESS, ETC)
