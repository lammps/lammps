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

#include "pair_rebo.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairREBO::PairREBO(LAMMPS *lmp) : PairAIREBO(lmp) {}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairREBO::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

  cutlj = 0.0;
  ljflag = torflag = 0;

  // this one parameter for C-C interactions is different in REBO vs AIREBO
  // see Favata, Micheletti, Ryu, Pugno, Comp Phys Comm (2016)

  PCCf_2_0 = 0.0;
}
