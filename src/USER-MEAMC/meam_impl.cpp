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
   Contributing author: Sebastian Hütter (OvGU)
------------------------------------------------------------------------- */

#include "meam.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

MEAM::MEAM(Memory *mem) : memory(mem)
{
  nmax = 0;
  rho = rho0 = rho1 = rho2 = rho3 = frhop = NULL;
  gamma = dgamma1 = dgamma2 = dgamma3 = arho2b = NULL;
  arho1 = arho2 = arho3 = arho3b = t_ave = tsq_ave = NULL;

  maxneigh = 0;
  scrfcn = dscrfcn = fcpair = NULL;
}

MEAM::~MEAM()
{
  meam_cleanup();

  memory->destroy(rho);
  memory->destroy(rho0);
  memory->destroy(rho1);
  memory->destroy(rho2);
  memory->destroy(rho3);
  memory->destroy(frhop);
  memory->destroy(gamma);
  memory->destroy(dgamma1);
  memory->destroy(dgamma2);
  memory->destroy(dgamma3);
  memory->destroy(arho2b);

  memory->destroy(arho1);
  memory->destroy(arho2);
  memory->destroy(arho3);
  memory->destroy(arho3b);
  memory->destroy(t_ave);
  memory->destroy(tsq_ave);

  memory->destroy(scrfcn);
  memory->destroy(dscrfcn);
  memory->destroy(fcpair);
}


