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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "fix_viscous_cuda.h"
#include "fix_viscous_cuda_cu.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "cuda_modify_flags.h"
#include "user_cuda.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace FixConstCuda;

/* ---------------------------------------------------------------------- */

FixViscousCuda::FixViscousCuda(LAMMPS *lmp, int narg, char **arg) :
  FixViscous(lmp, narg, arg)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

        cu_gamma=NULL;
}

/* ---------------------------------------------------------------------- */

FixViscousCuda::~FixViscousCuda()
{
        delete cu_gamma;
}

/* ---------------------------------------------------------------------- */

int FixViscousCuda::setmask()
{
  int mask = 0;
  mask |= POST_FORCE_CUDA;
 // mask |= POST_FORCE_RESPA;
 // mask |= MIN_POST_FORCE;
  return mask;
}


/* ---------------------------------------------------------------------- */

void FixViscousCuda::setup(int vflag)
{
   if(not cu_gamma)
   cu_gamma = new cCudaData<double, F_CFLOAT, x> (gamma,atom->ntypes+1);
   Cuda_FixViscousCuda_Init(&cuda->shared_data);
   cu_gamma->upload();
 // if (strcmp(update->integrate_style,"verlet/cuda") == 0)
    post_force(vflag);
 /* else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }*/
}

/* ---------------------------------------------------------------------- */

void FixViscousCuda::min_setup(int vflag)
{
  Cuda_FixViscousCuda_Init(&cuda->shared_data);
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixViscousCuda::post_force(int vflag)
{
  // apply drag force to atoms in group
  // direction is opposed to velocity vector
  // magnitude depends on atom type

  Cuda_FixViscousCuda_PostForce(&cuda->shared_data, groupbit,cu_gamma->dev_data());
}
