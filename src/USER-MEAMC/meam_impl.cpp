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

MEAM::MEAM(Memory* mem)
  : memory(mem)
{
  phir = phirar = phirar1 = phirar2 = phirar3 = phirar4 = phirar5 = phirar6 = NULL;

  nmax = 0;
  rho = rho0 = rho1 = rho2 = rho3 = frhop = NULL;
  gamma = dgamma1 = dgamma2 = dgamma3 = arho2b = NULL;
  arho1 = arho2 = arho3 = arho3b = t_ave = tsq_ave = NULL;

  maxneigh = 0;
  scrfcn = dscrfcn = fcpair = NULL;
}

MEAM::~MEAM()
{
  memory->destroy(this->phirar6);
  memory->destroy(this->phirar5);
  memory->destroy(this->phirar4);
  memory->destroy(this->phirar3);
  memory->destroy(this->phirar2);
  memory->destroy(this->phirar1);
  memory->destroy(this->phirar);
  memory->destroy(this->phir);

  memory->destroy(this->rho);
  memory->destroy(this->rho0);
  memory->destroy(this->rho1);
  memory->destroy(this->rho2);
  memory->destroy(this->rho3);
  memory->destroy(this->frhop);
  memory->destroy(this->gamma);
  memory->destroy(this->dgamma1);
  memory->destroy(this->dgamma2);
  memory->destroy(this->dgamma3);
  memory->destroy(this->arho2b);

  memory->destroy(this->arho1);
  memory->destroy(this->arho2);
  memory->destroy(this->arho3);
  memory->destroy(this->arho3b);
  memory->destroy(this->t_ave);
  memory->destroy(this->tsq_ave);

  memory->destroy(this->scrfcn);
  memory->destroy(this->dscrfcn);
  memory->destroy(this->fcpair);
}
