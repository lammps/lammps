/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Xiaowng Zhou (Sandia)
------------------------------------------------------------------------- */

#include "pair_eam_he_opt.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   inherit optimized compute() from PairEAMOpt,
   select the eam/he file format as in PairEAMHE
------------------------------------------------------------------------- */

PairEAMHEOpt::PairEAMHEOpt(LAMMPS *lmp) : PairEAMOpt(lmp)
{
  fileformat = FS;
  one_coeff = 1;
  he_flag = 1;
}
