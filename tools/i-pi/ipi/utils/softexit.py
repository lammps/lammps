"""Utility functions for killing the wrapper softly.

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
   Softexit: Concise class to manage cleaning up in case of an emergency exit.
"""

import traceback, sys
from ipi.utils.messages import verbosity, warning

__all__ = ['Softexit', 'softexit']


class Softexit(object):
   """Class to deal with stopping a simulation half way through.

   Holds the functions used to clean up a simulation that has been
   stopped early, either because of a SIGTERM signal or because the
   user has added an EXIT file to the directory in which it is 
   running. This will then properly shut down the socket interface,
   and print out a RESTART file for the appropriate time step.

   Attributes:
      flist: A list of functions used to close down the socket
         interface.
   """

   def __init__(self):
      """Initializes SoftExit."""

      self.flist = []

   def register(self, func):
      """Adds another function to flist.

      Args:
         func: The function to be added to flist.
      """

      self.flist.append(func)

   def trigger(self, message=""):
      """Halts the simulation.

      Prints out a warning message, then runs all the exit functions in flist
      before terminating the simulation.

      Args:
         message: The message to output to standard output.
      """

      if message != "":
         warning("Soft exit has been requested with message: '" + message + "'. Cleaning up.", verbosity.low)
      for f in self.flist:
         f()
      sys.exit()

softexit = Softexit()
