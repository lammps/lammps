"""Deals with creating the cell class.

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


Generates an cell class from a cell vector.

Classes:
   InputCell: Deals with creating the Cell object from a file, and
      writing the checkpoints.
"""

import numpy as np
from copy import copy
from ipi.engine.cell import *
from ipi.utils.inputvalue import *
from ipi.utils.units import UnitMap
from ipi.utils.messages import verbosity, warning

__all__ = [ 'InputCell' ]

class InputCell(InputArray):
   """Cell input class.

   Handles generating the appropriate cell class from the xml input file,
   and generating the xml checkpoint tags and data from an instance of the
   object.
   """

   attribs = copy(InputArray.attribs)

   default_help = "Deals with the cell parameters. Takes as array which can be used to initialize the cell vector matrix."
   default_label = "CELL"

   def __init__(self, help=None, dimension=None, units=None, default=None, dtype=None):
      """Initializes InputCell.

      Just calls the parent initialization function with appropriate arguments.
      """

      super(InputCell,self).__init__(dtype=float, dimension="length", default=default, help=help)

   def store(self, cell):
      """Takes a Cell instance and stores of minimal representation of it.

      Args:
         cell: A cell object.
      """

      super(InputCell,self).store(cell.h)
      self.shape.store((3,3))

   def fetch(self):
      """Creates a cell object.

      Returns:
         A cell object of the appropriate type and with the appropriate
         properties given the attributes of the InputCell object.
      """

      h = super(InputCell,self).fetch()
      h.shape = (3,3)

      return Cell(h=h)
