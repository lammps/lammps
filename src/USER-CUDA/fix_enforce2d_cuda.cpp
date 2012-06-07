/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

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

#include <cstring>
#include "fix_enforce2d_cuda.h"
#include "fix_enforce2d_cuda_cu.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "error.h"
#include "cuda.h"
#include "cuda_modify_flags.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace FixConstCuda;

/* ---------------------------------------------------------------------- */

FixEnforce2DCuda::FixEnforce2DCuda(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

  if (narg != 3) error->all(FLERR,"Illegal fix enforce2d command");
}

/* ---------------------------------------------------------------------- */

int FixEnforce2DCuda::setmask()
{
  int mask = 0;
  mask |= POST_FORCE_CUDA;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE_CUDA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEnforce2DCuda::init()
{
  if (domain->dimension == 3)
    error->all(FLERR,"Cannot use fix enforce2d/cuda with 3d simulation");
  if (atom->omega_flag)
    error->warning(FLERR,"Enforce2d/cuda does not support omega_flag on gpu yet. Will be handled on cpu.");

  if (atom->angmom_flag)
    error->warning(FLERR,"Enforce2d/cuda does not support angmom_flag (angular momentum) on gpu yet. Will be handled on cpu.");

  if (atom->torque_flag)
    error->warning(FLERR,"Enforce2d/cuda does not support torque_flag on gpu yet. Will be handled on cpu.");
}

/* ---------------------------------------------------------------------- */

void FixEnforce2DCuda::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
  {
    Cuda_FixEnforce2dCuda_Init(&cuda->shared_data);
    cuda->cu_f->upload();
    cuda->cu_v->upload();
    post_force(vflag);
    cuda->cu_f->download();
    cuda->cu_v->download();
  }
  else {
    int nlevels_respa = ((Respa *) update->integrate)->nlevels;
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(vflag,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixEnforce2DCuda::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixEnforce2DCuda::post_force(int vflag)
{
  Cuda_FixEnforce2dCuda_PostForce(&cuda->shared_data, groupbit);

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (atom->omega_flag) {
    double **omega = atom->omega;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        omega[i][0] = 0.0;
        omega[i][1] = 0.0;
      }
  }

  if (atom->angmom_flag) {
    double **angmom = atom->angmom;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        angmom[i][0] = 0.0;
        angmom[i][1] = 0.0;
      }
  }

  if (atom->torque_flag) {
    double **torque = atom->torque;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        torque[i][0] = 0.0;
        torque[i][1] = 0.0;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixEnforce2DCuda::post_force_respa(int vflag, int ilevel, int iloop)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixEnforce2DCuda::min_post_force(int vflag)
{
  post_force(vflag);
}
