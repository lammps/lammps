"""Utility functions for outputting messages, diagnostics and errors'

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
   Verbosity: Concise class to check the selected level of output

Functions:
   banner:    Prints the program welcome "screen"
   help:      Prints the input syntax help
   info:      Prints some information to standard output, depending on the level of verbosity
   warning:   Same as info, but with a "!W!" prefix and optionally printing a stack trace
"""

import traceback, sys

__all__ = ['Verbosity', 'verbosity',' help', 'banner', 'info', 'warning']


VERB_QUIET  = 0
VERB_LOW    = 1
VERB_MEDIUM = 2
VERB_HIGH   = 3
VERB_DEBUG  = 4

class Verbosity(object):
   """Class used to determine what to print to standard output.

   Attributes:
      level: Determines what level of output to print.
   """

   level = "low"

   def __getattr__(self, name):
      """Determines whether a certain verbosity level is
      less than or greater than the stored value.

      Used to decide whether or not a certain info or warning string
      should be output.

      Args:
         name: The verbosity level at which the info/warning string
            will be output.
      """

      if name is "quiet":
         return self.level >= VERB_QUIET
      elif name is "low":
         return self.level >= VERB_LOW
      elif name is "medium":
         return self.level >= VERB_MEDIUM
      elif name is "high":
         return self.level >= VERB_HIGH
      elif name is "debug":
         return self.level >= VERB_DEBUG

   def __setattr__(self, name, value):
      """Sets the verbosity level

      Args:
         name: The name of what to set. Should always be 'level'.
         value: The value to set the verbosity to.

      Raises:
         ValueError: Raised if either the name or the level is not
            a valid option.
      """

      if name == "level":
         if value == "quiet":
            level = VERB_QUIET
         elif value == "low":
            level = VERB_LOW
         elif value == "medium":
            level = VERB_MEDIUM
         elif value == "high":
            level = VERB_HIGH
         elif value == "debug":
            level = VERB_DEBUG
         else: 
            raise ValueError("Invalid verbosity level " + str(value) + " specified.")
         super(Verbosity,self).__setattr__("level", level)


verbosity = Verbosity()

def help():
   """Prints out a help string."""

   print """usage:  %s input """%sys.argv[0]

def banner():
   """Prints out a banner."""

   print """
 ____       ____       ____       ____  
/    \     /    \     /    \     /    \  
|  #################################  | 
\__#_/     \____/     \____/     \_#__/     
   #    _        _______  _____    #                 
   #   (_)      |_   __ \|_   _|   #       v. 1.0                       
   #   __  ______ | |__) | | |     #                                         
   Y  [  ||______||  ___/  | |     #      A Python interface for (ab initio)  
  0 0  | |       _| |_    _| |_    #      (path integral) molecular dynamics. 
   #  [___]     |_____|  |_____|   #
 __#_       ____       ____       _#__                
/  # \     /    \     /    \     / #  \                                      
|  #################################  |                                      
\____/     \____/     \____/     \____/        

   """


def info(text="", show=True ):
   """Prints a warning message.

   Args:
      text: The text of the information message.
      show: A boolean describing whether or not the message should be
         printed.
   """

   if not show:
      return
   print text

def warning(text="", show=True):
   """Prints a warning message.

   Args:
      text: The text of the information message.
      show: A boolean describing whether or not the message should be
         printed.
   """

   if not show:
      return
   if verbosity.debug:
      traceback.print_stack(file=sys.stdout)
   print " !W! " + text
