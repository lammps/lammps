"""Deals with creating the normal mode representation arrays.

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
   InputNormalModes: Deals with creating the normal mode objects.
"""

import numpy as np
from copy import copy
from ipi.engine.normalmodes import *
from ipi.utils.inputvalue import *
from ipi.utils.units import *

__all__ = ['InputNormalModes']

class InputNormalModes(InputArray):
   """ Storage class for NormalModes engine.

   Describes how normal-modes transformation and integration should be
   performed.

   Attributes:
      mode: Specifies the method by which the dynamical masses are created.
      transform: Specifies whether the normal mode calculation will be
         done using a FFT transform or a matrix multiplication.
   """

   attribs = copy(InputArray.attribs)
   attribs["mode"] = (InputAttribute, {"dtype"   : str,
                                       "default" : "rpmd",
                                       "help"    : "Specifies the technique to be used to calculate the dynamical masses. 'rpmd' simply assigns the bead masses the physical mass. 'manual' sets all the normal mode frequencies except the centroid normal mode manually. 'pa-cmd' takes an argument giving the frequency to set all the non-centroid normal modes to. 'wmax-cmd' is similar to 'pa-cmd', except instead of taking one argument it takes two ([wmax,wtarget]). The lowest-lying normal mode will be set to wtarget for a free particle, and all the normal modes will coincide at frequency wmax. ",
                                       "options" : ['pa-cmd', 'wmax-cmd', 'manual', 'rpmd']})
   attribs["transform"] = (InputValue,{"dtype"   : str,
                                       "default" : "fft",
                                       "help"    : "Specifies whether to calculate the normal mode transform using a fast Fourier transform or a matrix multiplication. For small numbers of beads the matrix multiplication may be faster.",
                                       "options" : ['fft', 'matrix']})

   default_help = "Deals with the normal mode transformations, including the adjustment of bead masses to give the desired ring polymer normal mode frequencies if appropriate. Takes as arguments frequencies, of which different numbers must be specified and which are used to scale the normal mode frequencies in different ways depending on which 'mode' is specified."
   default_label = "NORMALMODES"

   def __init__(self, help=None, dimension=None, default=None, dtype=None):
      """ Initializes InputNormalModes.

      Just calls the parent initialization function with appropriate arguments.
      """

      super(InputNormalModes,self).__init__(help=help, default=default, dtype=float, dimension="frequency")

   def store(self, nm):
      """Takes a normal modes instance and stores a minimal representation
      of it.

      Args:
         nm: A normal modes object.
      """

      super(InputNormalModes,self).store(nm.nm_freqs)
      self.mode.store(nm.mode)
      self.transform.store(nm.transform_method)

   def fetch(self):
      """Creates a normal modes object.

      Returns:
         A normal modes object.
      """

      super(InputNormalModes,self).check()
      return NormalModes(self.mode.fetch(), self.transform.fetch(), super(InputNormalModes,self).fetch() )
