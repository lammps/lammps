"""Classes to deal with output of simulation data.

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


Holds classes to deal with the output of different properties, trajectories
and the restart files.

Classes:
   PropertyOutput: Deals with outputting properties.
   TrajectoryOutput: Deals with outputting trajectories.
   CheckpointOutput: Deals with outputting restart files.
"""

import os
import numpy as np
import ipi.inputs.simulation
from ipi.utils.depend import *
from ipi.utils.io.io_xml import *
from ipi.engine.properties import getkey

__all__ = [ 'PropertyOutput', 'TrajectoryOutput', 'CheckpointOutput' ]

class PropertyOutput(dobject):
   """Class dealing with outputting a set of properties to file.

   Does not do any calculation, just manages opening a file, getting data
   from a Properties object and outputting with the desired stride.

   Attributes:
      filename: The name of the file to output to.
      outlist: A list of the properties to be output.
      stride: The number of steps that should be taken between outputting the
         data to file.
      flush: How often we should flush to disk.
      nout: Number of steps since data was last flushed.
      out: The output stream on which to output the properties.
      simul: The simulation object to get the data to be output from.
   """


   def __init__(self, filename="out", stride=1, flush=1, outlist=None):
      """Initializes a property output stream opening the corresponding
      file name.

      Also writes out headers.

      Args:
         filename: A string giving the name of the file to be output to.
         stride: An integer giving how many steps should be taken between
            outputting the data to file.
         flush: Number of writes to file between flushing data.
         outlist: A list of all the properties that should be output.
      """

      if outlist is None:
         outlist = np.zeros(0,np.dtype('|S1024'))
      self.filename = filename
      self.outlist = np.asarray(outlist,np.dtype('|S1024'))
      self.stride = stride
      self.flush = flush
      self.nout = 0
      self.out = None

   def bind(self, simul):
      """Binds output proxy to simulation object.

      Args:
         simul: A simulation object to be bound.
      """

      self.simul = simul

      # Checks as soon as possible if some asked-for properties are
      # missing or misspelled
      for what in self.outlist:
         key = getkey(what)
         if not key in self.simul.properties.property_dict.keys():
            print "Computable properties list: ", self.simul.properties.property_dict.keys()
            raise KeyError(key + " is not a recognized property")

      self.open_stream()

   def open_stream(self):
      """Opens the output stream."""

      try:
         self.out = open(self.filename, "a")
      except:
         raise ValueError("Could not open file " + self.filename + " for output")

      # print nice header if information is available on the properties
      if (self.simul.step == 0) :
         icol = 1
         for what in self.outlist:
            ohead = "# "
            key = getkey(what)
            prop = self.simul.properties.property_dict[key]

            if "size" in prop and prop["size"] > 1:
               ohead += "cols.  %3d-%-3d" % ( icol, icol+prop["size"] - 1 )
               icol += prop["size"]
            else:
               ohead += "column %3d    " % ( icol )
               icol += 1
            ohead += " --> %s " % (what)
            if "help" in prop:
               ohead += ": " + prop["help"]
            self.out.write(ohead + "\n")

   def close_stream():
      """Closes the output stream."""

      self.out.close()

   def write(self):
      """Outputs the required properties of the system.

      Note that properties are outputted using the same format as for the
      output to the xml checkpoint files, as specified in io_xml.

      Raises:
         KeyError: Raised if one of the properties specified in the output list
            are not contained in the property_dict member of properties.
      """

      if not (self.simul.step + 1) % self.stride == 0:
         return
      self.out.write("  ")
      for what in self.outlist:
         try:
            quantity = self.simul.properties[what]
         except KeyError:
            raise KeyError(what + " is not a recognized property")
         if not hasattr(quantity,"__len__") :
            self.out.write(write_type(float, quantity) + "   ")
         else:
            for el in quantity:
               self.out.write(write_type(float, el) + " ")

      self.out.write("\n")

      self.nout += 1
      if self.flush > 0 and self.nout >= self.flush :
         self.out.flush()
         os.fsync(self.out)  # we REALLY want to print out! pretty please OS let us do it.
         self.nout = 0


class TrajectoryOutput(dobject):
   """Class dealing with outputting atom-based properties as a
   trajectory file.

   Does not do any calculation, just manages opening a file, getting data
   from a Trajectories object and outputting with the desired stride.

   Attributes:
      filename: The (base) name of the file to output to.
      format: The format of the trajectory file to be created.
      what: The trajectory that needs to be output.
      stride: The number of steps that should be taken between outputting the
         data to file.
      out: The output stream on which to output the trajectories.
      flush: How often we should flush to disk.
      nout: Number of steps since data was last flushed.
      ibead: Index of the replica to print the trajectory of.
      cell_units: The units that the cell parameters are given in.
      simul: The simulation object to get the data to be output from.
   """

   def __init__(self, filename="out", stride=1, flush=1, what="", format="xyz", cell_units="atomic_unit", ibead=-1):
      """ Initializes a property output stream opening the corresponding
      file name.

      Also writes out headers.

      Args:
         filename: A string giving the name of the file to be output to.
         stride: An integer giving how many steps should be taken between
            outputting the data to file.
         flush: How often we should flush to disk
         what: A string specifying what trajectory should be output.
         format: A string specifying the type of trajectory file to be created.
         cell_units: A string specifying the units that the cell parameters are
            given in.
         ibead: If positive, prints out only the selected bead. If negative, prints out one file per bead.
      """

      self.filename = filename
      self.what = what
      self.stride = stride
      self.flush = flush
      self.ibead = ibead
      self.format = format
      self.cell_units = cell_units
      self.out = None
      self.nout = 0

   def bind(self, simul):
      """Binds output proxy to simulation object.

      Args:
         simul: A simulation object to be bound.
      """

      self.simul = simul

      # Checks as soon as possible if some asked-for trajs are missing or misspelled
      key = getkey(self.what)
      if not key in self.simul.trajs.traj_dict.keys():
         print "Computable trajectories list: ", self.simul.trajs.traj_dict.keys()
         raise KeyError(key + " is not a recognized output trajectory")

      self.open_stream()

   def open_stream(self):
      """Opens the output stream(s)."""

      if getkey(self.what) in [ "positions", "velocities", "forces", "extras" ]:
         # must write out trajectories for each bead, so must create b streams
         self.out = []
         for b in range(self.simul.beads.nbeads):
            # zero-padded bead number
            padb = ( ("%0" + str(int(1 + np.floor(np.log(self.simul.beads.nbeads)/np.log(10)))) + "d") % (b) )
            try:
               if (self.ibead < 0 or self.ibead == b):
                  if getkey(self.what) == "extras":
                     self.out.append(open(self.filename + "_" + padb, "a"))
                  else:
                     self.out.append(open(self.filename + "_" + padb + "." + self.format, "a"))
               else:
                  self.out.append(None) # creates null outputs if a
                                        # single bead output is chosen
            except:
               raise ValueError("Could not open file " + self.filename + "_" + padb + "." + self.format + " for output")
      else:
         try:
            self.out = ( open(self.filename + "." + self.format, "a") )
         except:
            raise ValueError("Could not open file " + self.filename + "." + self.format + " for output")

   def close_stream():
      """Closes the output stream."""

      if hasattr(self.out, "__getitem__"):
         for o in self.out:
            o.close()
      else:
         self.out.close()

   def write(self):
      """Writes out the required trajectories."""

      if not (self.simul.step + 1) % self.stride == 0:
         return

      doflush = False
      self.nout += 1
      if self.flush > 0 and self.nout >= self.flush :
         doflush = True
         self.nout = 0

      # quick-and-dirty way to check if a trajectory is "global" or per-bead
      # Checks to see if there is a list of files or just a single file.
      if hasattr(self.out, "__getitem__"):
         if self.ibead < 0:
            for b in range(len(self.out)):
               self.simul.trajs.print_traj(self.what, self.out[b], b, format=self.format, cell_units=self.cell_units, flush=doflush)
         elif self.ibead < len(self.out):
            self.simul.trajs.print_traj(self.what, self.out[self.ibead], self.ibead, format=self.format, cell_units=self.cell_units, flush=doflush)
         else:
            raise ValueError("Selected bead index " + str(self.ibead) + " does not exist for trajectory " + self.what)
      else:
         self.simul.trajs.print_traj(self.what, self.out, b=0, format=self.format, cell_units=self.cell_units, flush=doflush)


class CheckpointOutput(dobject):
   """Class dealing with outputting checkpoints.

   Saves the complete status of the simulation at regular intervals.

   Attributes:
      filename: The (base) name of the file to output to.
      step: the number of times a checkpoint has been written out.
      stride: The number of steps that should be taken between outputting the
         data to file.
      overwrite: If True, the checkpoint file is overwritten at each output.
         If False, will output to 'filename_step'. Note that no check is done
         on whether 'filename_step' exists already.
      simul: The simulation object to get the data to be output from.
      status: An input simulation object used to write out the checkpoint file.
   """


   def __init__(self, filename="restart", stride=1000, overwrite=True, step=0):
      """Initializes a checkpoint output proxy.

      Args:
         filename: A string giving the name of the file to be output to.
         stride: An integer giving how many steps should be taken between
            outputting the data to file.
         overwrite: If True, the checkpoint file is overwritten at each output.
            If False, will output to 'filename_step'. Note that no check is done
            on whether 'filename_step' exists already.
         step: The number of checkpoint files that have been created so far.
      """

      self.filename = filename
      self.step = step
      self.stride = stride
      self.overwrite = overwrite

   def bind(self, simul):
      """Binds output proxy to simulation object.

      Args:
         simul: A simulation object to be bound.
      """

      self.simul = simul
      self.status = ipi.inputs.simulation.InputSimulation()
      self.status.store(simul)

   def store(self):
      """Stores the current simulation status.

      Used so that, if halfway through a step a kill signal is received,
      we can output a checkpoint file corresponding to the beginning of the
      current step, which is the last time that both the velocities and
      positions would have been consistent.
      """

      self.status.store(self.simul)

   def write(self, store=True):
      """Writes out the required trajectories.

      Used for both the checkpoint files and the soft-exit restart file.
      We have slightly different behaviour for these two different types of
      checkpoint file, as the soft-exit files have their store() function
      called automatically, and we do not want this to be updated as the
      status of the simulation after a soft-exit call is unlikely to be in
      a consistent state. On the other hand, the standard checkpoint files
      are not automatically updated in this way, and we must manually store the
      current state of the system before writing them.

      Args:
         store: A boolean saying whether the state of the system should be
            stored before writing the checkpoint file.
      """

      if not (self.simul.step + 1) % self.stride == 0:
         return

      if self.overwrite:
         filename = self.filename
      else:
         filename = self.filename + "_" + str(self.step)

      if store:
         self.step += 1    # advances the step counter before saving, so next time the correct index will be loaded.
         self.store()
      check_file = open(filename, "w")
      check_file.write(self.status.write(name="simulation"))
      check_file.close()
