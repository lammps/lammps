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

#include <stdio.h>
#include <string.h>
#include "fix_dpd_energy_kokkos.h"
#include "atom_kokkos.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "modify.h"
#include "error.h"
#include "pair_dpd_fdt_energy.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDPDenergyKokkos::FixDPDenergyKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixDPDenergy(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixDPDenergyKokkos::initial_integrate(int vflag)
{
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  t_efloat_1d uCond = atomKK
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *duCond = pairDPDE->duCond;
  double *duMech = pairDPDE->duMech;

  for (int i = 0; i < nlocal; i++){
    uCond[i] += 0.5*update->dt*duCond[i];
    uMech[i] += 0.5*update->dt*duMech[i];
  }
}

/* ---------------------------------------------------------------------- */

void FixDPDenergyKokkos::final_integrate()
{
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *duCond = pairDPDE->duCond;
  double *duMech = pairDPDE->duMech;

  for (int i = 0; i < nlocal; i++){
    uCond[i] += 0.5*update->dt*duCond[i];
    uMech[i] += 0.5*update->dt*duMech[i];
  }
}
