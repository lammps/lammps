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
   Contributing author: Mike Parks (SNL)
------------------------------------------------------------------------- */

#include "fix_peri_neigh_omp.h"
#include "fix_omp.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_request.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

void FixPeriNeighOMP::init()
{
  if (!first) return;

  // determine status of neighbor flag of the omp package command
  int ifix = modify->find_fix("package_omp");
  int use_omp = 0;
  if (ifix >=0) {
     FixOMP * fix = static_cast<FixOMP *>(lmp->modify->fix[ifix]);
     if (fix->get_neighbor()) use_omp = 1;
  }

  // need a full neighbor list once

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->omp = use_omp;
  neighbor->requests[irequest]->occasional = 1;
}
