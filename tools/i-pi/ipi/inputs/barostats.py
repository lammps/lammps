"""Deals with creating the barostat class.

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
   InputBaro: Deals with creating the Barostat object from a file, and
      writing the checkpoints.
"""

import numpy as np
import ipi.engine.thermostats
from ipi.engine.barostats import *
from ipi.utils.inputvalue import *
from ipi.inputs.thermostats import *

__all__ = ['InputBaro']

class InputBaro(Input):
   """Barostat input class.

   Handles generating the appropriate barostat class from the xml input file,
   and generating the xml checkpoint tags and data from an
   instance of the object.

   Attributes:
      mode: An optional string giving the type of barostat used. Defaults to
         'rigid'.

   Fields:
      thermostat: A thermostat object giving the cell thermostat.
      tau: The time constant associated with the dynamics of the piston.
      p: The conjugate momentum to the volume degree of freedom.
   """

   attribs={ "mode": (InputAttribute, {"dtype"    : str,
                                   "default" : "dummy",
                                   "help"     : """The type of barostat.  Currently, only a 'isotropic' barostat is implemented, that combines
                                   ideas from the Bussi-Zykova-Parrinello barostat for classical MD with ideas from the
                                   Martyna-Hughes-Tuckerman centroid barostat for PIMD; see Ceriotti, More, Manolopoulos, Comp. Phys. Comm. 2013 for
                                   implementation details.""",
                                   "options"  : ["dummy", "isotropic"]}) }
   fields={ "thermostat": (InputThermo, {"default" : input_default(factory=ipi.engine.thermostats.Thermostat),
                                         "help"    : "The thermostat for the cell. Keeps the cell velocity distribution at the correct temperature. Note that the 'pile_l', 'pile_g', 'nm_gle' and 'nm_gle_g' options will not work for this thermostat."}),
            "tau": (InputValue, {"default" : 1.0,
                                  "dtype" : float,
                                  "dimension" : "time",
                                  "help"    : "The time constant associated with the dynamics of the piston."}),
            "p": (InputArray, {  "dtype"     : float,
                                 "default"   : input_default(factory=np.zeros, args = (0,)),
                                 "help"      : "Momentum (or momenta) of the piston.",
                                 "dimension" : "momentum" })
           }

   default_help = "Simulates an external pressure bath."
   default_label = "BAROSTAT"

   def store(self, baro):
      """Takes a barostat instance and stores a minimal representation of it.

      Args:
         baro: A barostat object.
      """

      super(InputBaro,self).store(baro)
      self.thermostat.store(baro.thermostat)
      self.tau.store(baro.tau)
      if type(baro) is BaroBZP:
         self.mode.store("isotropic")
         self.p.store(baro.p)
      elif type(baro) is Barostat:
         self.mode.store("dummy")
      else:
         raise TypeError("The type " + type(baro).__name__ + " is not a valid barostat type")


   def fetch(self):
      """Creates a barostat object.

      Returns:
         A barostat object of the appropriate type and with the appropriate
         thermostat given the attributes of the InputBaro object.
      """

      super(InputBaro,self).fetch()
      if self.mode.fetch() == "isotropic":
         baro = BaroBZP(thermostat=self.thermostat.fetch(), tau=self.tau.fetch())
         if self.p._explicit: baro.p = self.p.fetch()
      elif self.mode.fetch() == "dummy":
         baro = Barostat(thermostat=self.thermostat.fetch(), tau=self.tau.fetch())
      else:
         raise ValueError(self.mode.fetch() + " is not a valid mode of barostat")

      return baro
