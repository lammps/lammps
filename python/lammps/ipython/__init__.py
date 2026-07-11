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

from .wrapper import wrapper
from .magics import LammpsMagics

def load_ipython_extension(ipython):
    ipython.register_magics(LammpsMagics)
