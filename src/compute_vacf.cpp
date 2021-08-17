// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government ret
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_vacf.h"

#include "atom.h"
#include "update.h"
#include "group.h"
#include "modify.h"
#include "fix_store.h"
#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeVACF::ComputeVACF(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  id_fix(nullptr)
{
  if (narg < 3) error->all(FLERR,"Illegal compute vacf command");

  vector_flag = 1;
  size_vector = 4;
  extvector = 0;
  create_attribute = 1;

  // create a new fix STORE style
  // id = compute-ID + COMPUTE_STORE, fix group = compute group

  id_fix = utils::strdup(id + std::string("_COMPUTE_STORE"));
  fix = (FixStore *) modify->add_fix(fmt::format("{} {} STORE peratom 1 3",
                                                 id_fix, group->names[igroup]));

  // store current velocities in fix store array
  // skip if reset from restart file

  if (fix->restart_reset) fix->restart_reset = 0;
  else {
    double **voriginal = fix->astore;

    double **v = atom->v;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        voriginal[i][0] = v[i][0];
        voriginal[i][1] = v[i][1];
        voriginal[i][2] = v[i][2];
      } else voriginal[i][0] = voriginal[i][1] = voriginal[i][2] = 0.0;
  }

  // displacement vector

  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputeVACF::~ComputeVACF()
{
  // check nfix in case all fixes have already been deleted

  if (modify->nfix) modify->delete_fix(id_fix);

  delete [] id_fix;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeVACF::init()
{
  // set fix which stores original atom velocities

  int ifix = modify->find_fix(id_fix);
  if (ifix < 0) error->all(FLERR,"Could not find compute vacf fix ID");
  fix = (FixStore *) modify->fix[ifix];

  // nvacf = # of atoms in group

  nvacf = group->count(igroup);
}

/* ---------------------------------------------------------------------- */

void ComputeVACF::compute_vector()
{
  invoked_vector = update->ntimestep;

  double **voriginal = fix->astore;

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double vxsq,vysq,vzsq;
  double vacf[4];
  vacf[0] = vacf[1] = vacf[2] = vacf[3] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      vxsq = v[i][0] * voriginal[i][0];
      vysq = v[i][1] * voriginal[i][1];
      vzsq = v[i][2] * voriginal[i][2];
      vacf[0] += vxsq;
      vacf[1] += vysq;
      vacf[2] += vzsq;
      vacf[3] += vxsq + vysq + vzsq;
    }

  MPI_Allreduce(vacf,vector,4,MPI_DOUBLE,MPI_SUM,world);
  if (nvacf) {
    vector[0] /= nvacf;
    vector[1] /= nvacf;
    vector[2] /= nvacf;
    vector[3] /= nvacf;
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
------------------------------------------------------------------------- */

void ComputeVACF::set_arrays(int i)
{
  double **voriginal = fix->astore;
  double **v = atom->v;
  voriginal[i][0] = v[i][0];
  voriginal[i][1] = v[i][1];
  voriginal[i][2] = v[i][2];
}

