# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   https://www.lammps.org/ Sandia National Laboratories
#   LAMMPS Development team: developers@lammps.org
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
# -------------------------------------------------------------------------

################################################################################
# IPython/Jupyter Notebook additions
# Written by Richard Berger <richard.berger@outlook.com>
################################################################################

import io
import os
import sys
import tempfile
from IPython.core.magic import (Magics, magics_class, cell_magic)
import IPython.core.magic_arguments as magic_arguments

class OutputCapture(object):
  """ Utility class to capture LAMMPS library output """
  def __init__(self):
    self.stdout_fd = 1
    self.captured_output = ""

  def __enter__(self):
    self.tmpfile = tempfile.TemporaryFile(mode='w+b')

    sys.stdout.flush()

    # make copy of original stdout
    self.stdout_orig = os.dup(self.stdout_fd)

    # replace stdout and redirect to temp file
    os.dup2(self.tmpfile.fileno(), self.stdout_fd)
    return self

  def __exit__(self, exc_type, exc_value, traceback):
    os.dup2(self.stdout_orig, self.stdout_fd)
    os.close(self.stdout_orig)
    self.tmpfile.close()

  @property
  def output(self):
    sys.stdout.flush()
    self.tmpfile.flush()
    self.tmpfile.seek(0, io.SEEK_SET)
    self.captured_output = self.tmpfile.read().decode('utf-8')
    return self.captured_output

# -------------------------------------------------------------------------

@magics_class
class LammpsMagics(Magics):
    @magic_arguments.magic_arguments()
    @magic_arguments.argument('output', type=str, default='', nargs='?',
        help="""The name of the variable in which to store output.

        If unspecified, captured output is discarded.
        """
    )
    @cell_magic
    def capture_lammps_output(self, line, cell):
        """run the cell, capturing LAMMPS stdout and stderr."""
        args = magic_arguments.parse_argstring(self.capture_lammps_output, line)
        with OutputCapture() as capture:
            self.shell.run_cell(cell)
            if args.output:
                self.shell.user_ns[args.output] = capture.output
