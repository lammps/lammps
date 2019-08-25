"""Deals with creating the ensembles class.

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


Classes:
   InputEnsemble: Deals with creating the Ensemble object from a file, and
      writing the checkpoints.
"""

import numpy as np
import ipi.engine.thermostats
import ipi.engine.initializer
import ipi.engine.barostats
from ipi.engine.ensembles import *
from ipi.utils.inputvalue import *
from ipi.inputs.barostats import *
from ipi.inputs.thermostats import *
from ipi.inputs.initializer import *
from ipi.utils.units import *

__all__ = ['InputEnsemble']

class InputEnsemble(Input):
   """Ensemble input class.

   Handles generating the appropriate ensemble class from the xml input file,
   and generating the xml checkpoint tags and data from an instance of the
   object.

   Attributes:
      mode: An optional string giving the mode of ensemble to be simulated.
         Defaults to 'unknown'.

   Fields:
      thermostat: The thermostat to be used for constant temperature dynamics.
      barostat: The barostat to be used for constant pressure or stress
         dynamics.
      timestep: An optional float giving the size of the timestep in atomic
         units. Defaults to 1.0.
      temperature: An optional float giving the temperature in Kelvin. Defaults
         to 1.0.
      pressure: An optional float giving the external pressure in atomic units.
         Defaults to 1.0.
      fixcom: An optional boolean which decides whether the centre of mass
         motion will be constrained or not. Defaults to False.
      replay_file: An optional string that gives an input file name to get
         a trajectory to be re-run.
   """

   attribs={"mode"  : (InputAttribute, {"dtype"   : str,
                                    "help"    : "The ensemble that will be sampled during the simulation. 'replay' means that a simulation is restarted from a previous simulation.",
                                    "options" : ['nve', 'nvt', 'npt', 'replay']}) }
   fields={"thermostat" : (InputThermo, {"default"   : input_default(factory=ipi.engine.thermostats.Thermostat),
                                         "help"      : "The thermostat for the atoms, keeps the atom velocity distribution at the correct temperature."} ),
           "barostat" : (InputBaro, {"default"       : input_default(factory=ipi.engine.barostats.Barostat),
                                     "help"          : InputBaro.default_help}),
           "timestep": (InputValue, {"dtype"         : float,
                                     "default"       : 1.0,
                                     "help"          : "The time step.",
                                     "dimension"     : "time"}),
           "temperature" : (InputValue, {"dtype"     : float,
                                         "default"   : 1.0,
                                         "help"      : "The temperature of the system.",
                                         "dimension" : "temperature"}),
           "pressure" : (InputValue, {"dtype"        : float,
                                      "default"      : 1.0,
                                      "help"         : "The external pressure.",
                                      "dimension"    : "pressure"}),
           "fixcom": (InputValue, {"dtype"           : bool,
                                   "default"         : True,
                                   "help"            : "This describes whether the centre of mass of the particles is fixed."}),
           "replay_file": (InputInitFile, {"default" : input_default(factory=ipi.engine.initializer.InitBase),
                           "help"            : "This describes the location to read a trajectory file from."})
         }

   default_help = "Holds all the information that is ensemble specific, such as the temperature and the external pressure, and the thermostats and barostats that control it."
   default_label = "ENSEMBLE"

   def store(self, ens):
      """Takes an ensemble instance and stores a minimal representation of it.

      Args:
         ens: An ensemble object.
      """

      super(InputEnsemble,self).store(ens)
      if type(ens) is ReplayEnsemble:
         self.mode.store("rerun")
         tens = 0
      elif type(ens) is NVEEnsemble:
         self.mode.store("nve")
         tens = 1
      elif type(ens) is NVTEnsemble:
         self.mode.store("nvt")
         tens = 2
      elif type(ens) is NPTEnsemble:
         self.mode.store("npt")
         tens = 3

      self.timestep.store(ens.dt)
      self.temperature.store(ens.temp)

      if tens == 0:
         self.replay_file.store(ens.intraj)
      if tens > 1:
         self.thermostat.store(ens.thermostat)
         self.fixcom.store(ens.fixcom)
      if tens > 2:
         self.barostat.store(ens.barostat)
      if tens == 3:
         self.pressure.store(ens.pext)


   def fetch(self):
      """Creates an ensemble object.

      Returns:
         An ensemble object of the appropriate mode and with the appropriate
         objects given the attributes of the InputEnsemble object.
      """

      super(InputEnsemble,self).fetch()

      if self.mode.fetch() == "nve" :
         ens = NVEEnsemble(dt=self.timestep.fetch(),
            temp=self.temperature.fetch(), fixcom=self.fixcom.fetch())
      elif self.mode.fetch() == "nvt" :
         ens = NVTEnsemble(dt=self.timestep.fetch(),
            temp=self.temperature.fetch(), thermostat=self.thermostat.fetch(), fixcom=self.fixcom.fetch())
      elif self.mode.fetch() == "npt" :
         ens = NPTEnsemble(dt=self.timestep.fetch(),
            temp=self.temperature.fetch(), thermostat=self.thermostat.fetch(), fixcom=self.fixcom.fetch(),
                  pext=self.pressure.fetch(), barostat=self.barostat.fetch() )
      elif self.mode.fetch() == "replay":
         ens = ReplayEnsemble(dt=self.timestep.fetch(),
            temp=self.temperature.fetch(),fixcom=False,intraj=self.replay_file.fetch() )
      else:
         raise ValueError("'" + self.mode.fetch() + "' is not a supported ensemble mode.")

      return ens

   def check(self):
      """Function that deals with optional arguments.

      Makes sure that if the ensemble requires a thermostat or barostat that
      they have been defined by the user and not given the default values.
      """

      super(InputEnsemble,self).check()
      if self.mode.fetch() == "nvt":
         if self.thermostat._explicit == False:
            raise ValueError("No thermostat tag supplied for NVT simulation")
      if self.mode.fetch() == "npt":
         if self.thermostat._explicit == False:
            raise ValueError("No thermostat tag supplied for NPT simulation")
         if self.barostat._explicit == False:
            raise ValueError("No barostat tag supplied for NPT simulation")
         if self.barostat.thermostat._explicit == False:
            raise ValueError("No thermostat tag supplied in barostat for NPT simulation")

      if self.timestep.fetch() <= 0:
         raise ValueError("Non-positive timestep specified.")
      if self.temperature.fetch() <= 0:
            raise ValueError("Non-positive temperature specified.")

      if self.mode.fetch() == "npt":
         if not self.pressure._explicit:
            raise ValueError("Pressure should be supplied for constant pressure simulation")
      if self.mode.fetch() == "npt" or self.mode.fetch() == "nvt":
         if not self.temperature._explicit:
            raise ValueError("Temperature should be supplied for constant temperature simulation")
