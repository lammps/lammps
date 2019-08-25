/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
     James Fischer, High Performance Technologies, Inc.
     Charles Cornwell, High Performance Technologies, Inc.
     David Richie, Stone Ridge Technology
     Vincent Natoli, Stone Ridge Technology
------------------------------------------------------------------------- */

#include "pair_eam_fs_opt.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   multiple inheritance from two parent classes
   invoke constructor of grandparent class, then of each parent
   inherit optimized compute() from PairEAMOpt
   inherit everything else from PairEAMFS
------------------------------------------------------------------------- */

PairEAMFSOpt::PairEAMFSOpt(LAMMPS *lmp) :
  PairEAM(lmp), PairEAMFS(lmp), PairEAMOpt(lmp) {}
