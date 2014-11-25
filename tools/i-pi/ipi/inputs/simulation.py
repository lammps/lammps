"""Deals with creating the simulation class.

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
   InputSimulation: Deals with creating the Simulation object from a file, and
      writing the checkpoints.
"""

__all__ = ['InputSimulation']

import numpy as np
import os.path, sys
import ipi.engine.simulation
from ipi.utils.depend import *
from ipi.utils.inputvalue import *
from ipi.utils.units  import *
from ipi.utils.prng   import *
from ipi.utils.io     import *
from ipi.utils.io.io_xml import *
from ipi.utils.messages import verbosity
from ipi.inputs.forces import InputForces
from ipi.inputs.prng import InputRandom
from ipi.inputs.initializer import InputInitializer
from ipi.inputs.beads import InputBeads
from ipi.inputs.cell import InputCell
from ipi.inputs.ensembles import InputEnsemble
from ipi.inputs.outputs import InputOutputs
from ipi.inputs.normalmodes import InputNormalModes
from ipi.engine.normalmodes import NormalModes
from ipi.engine.atoms import Atoms
from ipi.engine.beads import Beads
from ipi.engine.cell import Cell
from ipi.engine.initializer import Initializer

class InputSimulation(Input):
   """Simulation input class.

   Handles generating the appropriate forcefield class from the xml input file,
   and generating the xml checkpoint tags and data from an instance of the
   object.

   Attributes:
      verbosity: A string saying how much should be output to standard output.

   Fields:
      force: A restart force instance. Used as a model for all the replicas.
      ensemble: A restart ensemble instance.
      beads: A restart beads instance.
      normal_modes: Setup of normal mode integrator.
      cell: A restart cell instance.
      output: A list of the required outputs.
      prng: A random number generator object.
      step: An integer giving the current simulation step. Defaults to 0.
      total_steps: The total number of steps. Defaults to 1000
      total_time:  The wall clock time limit. Defaults to 0 (no limit).
      initialize: An array of strings giving all the quantities that should
         be output.
   """

   fields = { "forces" :   (InputForces,    { "help"  : InputForces.default_help }),
             "ensemble": (InputEnsemble, { "help"  : InputEnsemble.default_help } ),
             "prng" :    (InputRandom,   { "help"  : InputRandom.default_help,
                                         "default" : input_default(factory=Random)} ),
             "initialize" : (InputInitializer, { "help" : InputInitializer.default_help,
                                                "default" : input_default(factory=Initializer) } ),
             "beads" :   (InputBeads, { "help"     : InputBeads.default_help,
                                        "default"  : input_default(factory=Beads, kwargs={'natoms': 0, 'nbeads': 0}) } ),
             "normal_modes" :   (InputNormalModes, { "help"     : InputNormalModes.default_help,
                                        "default"  : input_default(factory=NormalModes, kwargs={'mode': "rpmd"}) } ),
             "cell" :    (InputCell,   { "help"    : InputCell.default_help,
                                        "default"  : input_default(factory=Cell) }),
             "output" :  (InputOutputs, { "help"   : InputOutputs.default_help,
                                          "default": input_default(factory=InputOutputs.make_default)  }),
             "step" :       ( InputValue, { "dtype"    : int,
                                            "default"  : 0,
                                            "help"     : "The current simulation time step." }),
             "total_steps": ( InputValue, { "dtype"    : int,
                                            "default"  : 1000,
                                            "help"     : "The total number of steps that will be done. If 'step' is equal to or greater than 'total_steps', then the simulation will finish." }),
             "total_time" :       ( InputValue, { "dtype"    : float,
                                            "default"  : 0,
                                            "help"     : "The maximum wall clock time (in seconds)." }),
                                             }

   attribs = { "verbosity" : (InputAttribute, { "dtype"   : str,
                                      "default" : "low",
                                      "options" : [ "quiet", "low", "medium", "high", "debug" ],
                                      "help"    : "The level of output on stdout."
                                         })
             }

   default_help = "This is the top level class that deals with the running of the simulation, including holding the simulation specific properties such as the time step and outputting the data."
   default_label = "SIMULATION"

   def store(self, simul):
      """Takes a simulation instance and stores a minimal representation of it.

      Args:
         simul: A simulation object.
      """

      super(InputSimulation,self).store()
      self.forces.store(simul.flist)
      self.ensemble.store(simul.ensemble)

      self.beads.store(simul.beads)

      self.normal_modes.store(simul.nm)
      self.cell.store(simul.cell)
      self.prng.store(simul.prng)
      self.step.store(simul.step)
      self.total_steps.store(simul.tsteps)
      self.total_time.store(simul.ttime)
      self.output.store(simul.outputs)

      # this we pick from the messages class. kind of a "global" but it seems to
      # be the best way to pass around the (global) information on the level of output.
      if verbosity.debug:
         self.verbosity.store("debug")
      elif verbosity.high:
         self.verbosity.store("high")
      elif verbosity.medium:
         self.verbosity.store("medium")
      elif verbosity.low:
         self.verbosity.store("low")
      elif verbosity.quiet:
         self.verbosity.store("quiet")
      else:
         raise ValueError("Invalid verbosity level")

   def fetch(self):
      """Creates a simulation object.

      Returns:
         A simulation object of the appropriate type and with the appropriate
         properties and other objects given the attributes of the
         InputSimulation object.

      Raises:
         TypeError: Raised if one of the file types in the stride keyword
            is incorrect.
      """

      super(InputSimulation,self).fetch()

      # small hack: initialize here the verbosity level -- we really assume to have
      # just one simulation object
      verbosity.level=self.verbosity.fetch()

      # this creates a simulation object which gathers all the little bits
      #TODO use named arguments since this list is a bit too long...
      rsim = ipi.engine.simulation.Simulation(self.beads.fetch(), self.cell.fetch(),
               self.forces.fetch(), self.ensemble.fetch(), self.prng.fetch(),
                  self.output.fetch(), self.normal_modes.fetch(),
                     self.initialize.fetch(), self.step.fetch(),
                        tsteps=self.total_steps.fetch(),
                           ttime=self.total_time.fetch())

      # this does all of the piping between the components of the simulation
      rsim.bind()

      return rsim

   def check(self):
      """Function that deals with optional arguments.

      Deals with the difference between classical and PI dynamics. If there is
      no beads argument, the bead positions are generated from the atoms, with
      the necklace being fixed at the atom position. Similarly, if no nbeads
      argument is specified a classical simulation is done.

      Raises:
         TypeError: Raised if no beads or atoms attribute is defined.
      """

      super(InputSimulation,self).check()
      if self.total_steps.fetch() <= self.step.fetch():
         raise ValueError("Current step greater than total steps, no dynamics will be done.")
