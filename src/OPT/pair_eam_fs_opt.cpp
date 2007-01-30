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

#include "pair_eam_fs_opt.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   multiple inheritance from two parent classes
   invoke constructor of grandparent class, then of each parent
------------------------------------------------------------------------- */

PairEAMFSOpt::PairEAMFSOpt(LAMMPS *lmp) :
  PairEAM(lmp), PairEAMFS(lmp), PairEAMOpt(lmp) {}
