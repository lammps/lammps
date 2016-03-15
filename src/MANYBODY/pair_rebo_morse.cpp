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

#include "pair_rebo_morse.h"
#include "force.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairREBOMorse::PairREBOMorse(LAMMPS *lmp) : PairAIREBO(lmp) {}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairREBOMorse::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

  cutlj = 0.0;
  ljflag = torflag = 0;
  morseflag = 1;
}
