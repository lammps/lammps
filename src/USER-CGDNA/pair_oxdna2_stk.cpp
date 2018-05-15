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
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_oxdna2_stk.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairOxdna2Stk::PairOxdna2Stk(LAMMPS *lmp) : PairOxdnaStk(lmp)
{

}

/* ---------------------------------------------------------------------- */

PairOxdna2Stk::~PairOxdna2Stk()
{

}

/* ----------------------------------------------------------------------
   return temperature dependent oxDNA2 stacking strength
------------------------------------------------------------------------- */

double PairOxdna2Stk::stacking_strength(double T)
{
  double eps;

  eps = 1.3523 + 2.6717 * T;

  return eps;
}
