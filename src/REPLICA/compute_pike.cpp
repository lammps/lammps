// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_pike.h"

#include "atom.h"
#include "update.h"
#include "force.h"
#include "error.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "universe.h"
#include "update.h"
#include "utils.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <mpi.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePIKE::ComputePIKE(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 7) error->all(FLERR,"Illegal compute pike command");

  int iarg = 3;

  for (int iarg = 3; iarg < narg - 1; iarg += 2) {
    if (strcmp(arg[iarg], "dp_pimd") == 0) {
      id_fix = arg[iarg+1];
    }

    if (strcmp(arg[iarg], "bead") == 0) {
      if (strcmp(arg[iarg+1], "yes") == 0) bead_flag = 1;
      else if (strcmp(arg[iarg+1], "no") == 0) bead_flag = 0;
    }
  }

  printf("bead_flag=%d\n", bead_flag);
  scalar_flag = 1;
  extscalar = 1;
}

/* ---------------------------------------------------------------------- */

void ComputePIKE::init()
{
  pfactor = 0.5 * force->mvv2e;
  ifix = modify->find_fix(id_fix);
  fix_dppimd = (FixDPPimd *) modify->fix[ifix];
}

/* ---------------------------------------------------------------------- */

double ComputePIKE::compute_scalar()
{
  printf("coming into compute_scalar\n");
  invoked_scalar = update->ntimestep;

  double **v = atom->v;
  double *rmass = atom->rmass;
  printf("before fix_dppimd->mass\n");
  double *mass = fix_dppimd->mass;
  //for(int i=0; i<atom->nlocal; i++) printf("%.2f ", mass[atom->type[i]]);
  printf("\n");
  printf("after fix_dppimd->mass\n");
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  printf("2 oming into compute_scalar\n");
  double ke = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        ke += rmass[i] * (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        ke += mass[type[i]] *
          (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
  }

//   MPI_Allreduce(&ke,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  double ke_total = 0.0;
  printf("bead_flag=%d\n", this->bead_flag);
  if (bead_flag == 1){
      printf("yes bead!\n");
      MPI_Allreduce(&ke,&ke_total,1,MPI_DOUBLE,MPI_SUM,world);
  }
  else if (bead_flag == 0){
      printf("no bead!\n");
    //   MPI_Allreduce(&ke,&ke_total,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
  }
  scalar *= pfactor;

  return scalar;
}
