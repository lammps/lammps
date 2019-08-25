"""Deals with creating the output objects.

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
   InputOutputs: Creates a list of all the output objects.
   InputProperties: Deals with property output.
   InputTrajectory: Deals with trajectory output.
   InputCheckpoint: Deals with restart file output.
"""
import numpy as np
from copy import copy
import ipi.engine.outputs
from ipi.utils.depend import *
from ipi.utils.inputvalue import *
from ipi.engine.properties import getkey

__all__=['InputOutputs', 'InputProperties', 'InputTrajectory',
         'InputCheckpoint']

class InputProperties(InputArray):
   """Simple input class to describe output for properties.

   Storage class for PropertyOutput.

   Attributes:
      filename: The name of the file to output to.
      stride: The number of steps that should be taken between outputting the
         data to file.
      flush: An integer describing how often the output streams are flushed,
         so that it doesn't wait for the buffer to fill before outputting to
         file.
   """

   default_help = """This class deals with the output of properties to one file. Between each property tag there should be an array of strings, each of which specifies one property to be output."""
   default_label = "PROPERTIES"

   attribs = copy(InputArray.attribs)
   attribs["filename"] = (InputAttribute,{ "dtype" : str, "default": "out",
                                           "help": "A string to specify the name of the file that is output. The file name is given by 'prefix'.'filename' + format_specifier. The format specifier may also include a number if multiple similar files are output."} )
   attribs["stride"] = (InputAttribute,{ "dtype" : int, "default": 1,
                                         "help": "The number of steps between successive writes." } )
   attribs["flush"] = (InputAttribute, {"dtype"    : int,    "default"  : 1,
                                   "help"     : "How often should streams be flushed. 1 means each time, zero means never." })

   def __init__(self, help=None,  default=None, dtype=None, dimension=None):
      """Initializes InputProperties.

      Just calls the parent initialization function with appropriate arguments.
      """

      super(InputProperties,self).__init__(help=help, default=default, dtype=str, dimension=dimension)

   def fetch(self):
      """Returns a PropertyOutput object."""

      return ipi.engine.outputs.PropertyOutput(filename=self.filename.fetch(),
        stride=self.stride.fetch(), flush=self.flush.fetch(), outlist=super(InputProperties,self).fetch())

   def store(self, prop):
      """Stores a PropertyOutput object."""

      super(InputProperties,self).store(prop.outlist)
      self.stride.store(prop.stride)
      self.flush.store(prop.flush)
      self.filename.store(prop.filename)

   def check(self):
      """Checks for optional parameters."""

      super(InputProperties,self).check()
      if self.stride.fetch() < 0:
         raise ValueError("The stride length for the properties file output must be positive.")


class InputTrajectory(InputValue):
   """Simple input class to describe output for trajectories.

   Storage class for TrajectoryOutput.

   Attributes:
      filename: The (base) name of the file to output to.
      stride: The number of steps that should be taken between outputting the
         data to file.
      format: The format of the trajectory output file.
      cell_units: The units that the cell parameters are given in.
      bead: If the trajectory is a per-bead property, this can be used to
         specify a single bead to output. If negative, it defaults to
         the centroid.
      flush: An integer describing how often the output streams are flushed,
         so that it doesn't wait for the buffer to fill before outputting to
         file.
   """

   default_help = """This class defines how one trajectory file should be output. Between each trajectory tag one string should be given, which specifies what data is to be output."""
   default_label = "TRAJECTORY"

   attribs = copy(InputValue.attribs)
   attribs["filename"] = (InputAttribute,{ "dtype" : str, "default": "traj",
                                           "help": "A string to specify the name of the file that is output. The file name is given by 'prefix'.'filename' + format_specifier. The format specifier may also include a number if multiple similar files are output."} )
   attribs["stride"] = (InputAttribute,{ "dtype" : int, "default": 1,
                                         "help": "The number of steps between successive writes." } )
   attribs["format"] = (InputAttribute,{ "dtype" : str, "default": "xyz",
                                       "help": "The output file format.",
                                       "options": ['xyz', 'pdb'] } )
   attribs["cell_units"] = (InputAttribute,{ "dtype" : str, "default": "",
                                       "help": "The units for the cell dimensions." } )
   attribs["bead"] = (InputAttribute,{ "dtype" : int, "default": -1,
                                         "help": "Print out only the specified bead. A negative value means print all." } )
   attribs["flush"] = (InputAttribute, {"dtype"    : int,    "default"  : 1,
                                   "help"     : "How often should streams be flushed. 1 means each time, zero means never." })

   def __init__(self, help=None,  default=None, dtype=None, dimension=None):
      """Initializes InputTrajectory.

      Just calls the parent initialization function with appropriate arguments.
      """

      super(InputTrajectory,self).__init__(help=help, default=default, dtype=str, dimension=dimension)

   def fetch(self):
      """Returns a TrajectoryOutput object."""

      return ipi.engine.outputs.TrajectoryOutput(filename=self.filename.fetch(), stride=self.stride.fetch(),
               flush=self.flush.fetch(), what=super(InputTrajectory,self).fetch(),
               format=self.format.fetch(), cell_units=self.cell_units.fetch(), ibead=self.bead.fetch())

   def store(self, traj):
      """Stores a PropertyOutput object."""

      super(InputTrajectory,self).store(traj.what)
      self.stride.store(traj.stride)
      self.flush.store(traj.flush)
      self.filename.store(traj.filename)
      self.format.store(traj.format)
      self.cell_units.store(traj.cell_units)
      self.bead.store(traj.ibead)

   def check(self):
      """Checks for optional parameters."""

      super(InputTrajectory,self).check()
      if self.stride.fetch() < 0:
         raise ValueError("The stride length for the trajectory file output must be positive.")


class InputCheckpoint(InputValue):
   """Simple input class to describe output for properties.

   Storage class for CheckpointOutput.

   Attributes:
      filename: The (base) name of the file to output to.
      stride: The number of steps that should be taken between outputting the
         data to file.
      overwrite: whether checkpoints should be overwritten, or multiple
         files output.
   """

   default_help = """This class defines how a checkpoint file should be output. Optionally, between the checkpoint tags, you can specify one integer giving the current step of the simulation. By default this integer will be zero."""
   default_label = "CHECKPOINT"

   attribs=copy(InputValue.attribs)
   attribs["filename"] = (InputAttribute,{ "dtype" : str, "default": "restart",
                                           "help": "A string to specify the name of the file that is output. The file name is given by 'prefix'.'filename' + format_specifier. The format specifier may also include a number if multiple similar files are output."} )
   attribs["stride"] = (InputAttribute,{ "dtype" : int, "default": 1,
                                         "help": "The number of steps between successive writes." } )
   attribs["overwrite"] = (InputAttribute,{ "dtype" : bool, "default": True,
                                            "help": "This specifies whether or not each consecutive checkpoint file will overwrite the old one."} )

   def __init__(self, help=None,  default=None, dtype=None, dimension=None):
      """Initializes InputCheckpoint.

      Just calls the parent initialization function with appropriate arguments.
      """

      super(InputCheckpoint,self).__init__(help=help, default=default, dtype=int, dimension=dimension)

   def fetch(self):
      """Returns a CheckpointOutput object."""

      step = super(InputCheckpoint,self).fetch()
      return ipi.engine.outputs.CheckpointOutput(self.filename.fetch(), self.stride.fetch(), self.overwrite.fetch(), step=step )

   def parse(self, xml=None, text=""):
      """Overwrites the standard parse function so that we can specify this tag
      in the input without any data.

      We can use the syntax <checkpoint /> to do this

      Args:
         xml: An xml node containing all the data for the parent tag.
         text: The data to read the data from. Will be None if we have not
            specified any data.
      """

      # just a quick hack to allow an empty element
      try:
         super(InputCheckpoint,self).parse(xml,text)
      except: #TODO make this except a specific exception, not every one
         self.value = 0  #This could hide actual errors, at least in theory.

   def store(self, chk):
      """Stores a CheckpointOutput object."""

      super(InputCheckpoint,self).store(chk.step)
      self.stride.store(chk.stride)
      self.filename.store(chk.filename)
      self.overwrite.store(chk.overwrite)

   def check(self):
      """Checks for optional parameters."""

      super(InputCheckpoint,self).check()
      if self.stride.fetch() < 0:
         raise ValueError("The stride length for the checkpoint file output must be positive.")


class InputOutputs(Input):
   """ List of outputs input class.

   An example of a dynamic input class: a variable number of tags might be
   present, corresponding to different output requests. This allows for
   instance to print multiple property outputs, with different content
   and/or output frequency.

   Attributes:
      prefix: A string that will be appended to all output files from this
         simulation.

   Dynamic fields:
      trajectory: Specifies a trajectory to be output
      properties: Specifies some properties to be output.
      checkpoint: Specifies a checkpoint file to be output.
   """

   attribs = { "prefix" : ( InputAttribute, { "dtype" : str,
                                          "default"  : "i-pi",
                                          "help"     : "A string that will be prepended to each output file name. The file name is given by 'prefix'.'filename' + format_specifier. The format specifier may also include a number if multiple similar files are output." })
             }

   dynamic = {  "properties" : (InputProperties, { "help" : "Each of the properties tags specify how to create a file in which one or more properties are written, one line per frame. " } ),
               "trajectory" : (InputTrajectory, { "help" : "Each of the trajectory tags specify how to create a trajectory file, containing a list of per-atom coordinate properties. " } ),
               "checkpoint" : (InputCheckpoint, { "help" : "Each of the checkpoint tags specify how to create a checkpoint file, which can be used to restart a simulation. " } ),
            }

   default_help = """This class defines how properties, trajectories and checkpoints should be output during the simulation. May contain zero, one or many instances of properties, trajectory or checkpoint tags, each giving instructions on how one output file should be created and managed."""
   default_label = "OUTPUTS"

   @classmethod
   def make_default(cls):
      """Used to make the default value of the outputs class for use when no
      output is specified.

      Needed since this is a fairly complicated default, with many mutable
      objects, and the default has to be generated by a function that does not
      use any mutable objects as arguments.
      """

      return [ ipi.engine.outputs.PropertyOutput(filename="i-pi.md", stride=10, outlist=[ "time", "step", "conserved", "temperature", "potential", "kinetic_cv" ] ),
               ipi.engine.outputs.TrajectoryOutput(filename="i-pi.pos", stride=100, what="positions", format="xyz"),
               ipi.engine.outputs.CheckpointOutput(filename="i-pi.checkpoint", stride=1000, overwrite=True)]

   def fetch(self):
      """Returns a list of the output objects included in this dynamic
      container.

      Returns:
         A list of tuples, with each tuple being of the form ('type', 'object')
         where 'type' is the type of output object and 'object' is a particular
         object of that type.
      """

      super(InputOutputs, self).fetch()
      outlist = [ p.fetch() for (n, p) in self.extra ]
      prefix = self.prefix.fetch()
      if not prefix == "":
         for p in outlist:
            p.filename = prefix + "." + p.filename

      return outlist

   def store(self, plist):
      """ Stores a list of the output objects, creating a sequence of
      dynamic containers.

      Args:
         plist: A list of tuples, with each tuple being of the form
            ('type', 'object') where 'type' is the type of forcefield and
            'object' is a particular object of that type.
      """

      super(InputOutputs, self).store()
      self.extra = []

      self.prefix.store("")
      for el in plist:
         if (isinstance(el, ipi.engine.outputs.PropertyOutput)):
            ip = InputProperties()
            ip.store(el)
            self.extra.append(("properties", ip))
         elif (isinstance(el, ipi.engine.outputs.TrajectoryOutput)):
            ip = InputTrajectory()
            ip.store(el)
            self.extra.append(("trajectory", ip))
         elif (isinstance(el, ipi.engine.outputs.CheckpointOutput)):
            ip = InputCheckpoint()
            ip.store(el)
            self.extra.append(("checkpoint", ip))
