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
#include <cstring>
#include <cstdlib>
#include "fix_set_force_cuda.h"
#include "fix_set_force_cuda_cu.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "force.h"
#include "user_cuda.h"
#include "memory.h"
#include "cuda_modify_flags.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace FixConstCuda;

/* ---------------------------------------------------------------------- */

FixSetForceCuda::FixSetForceCuda(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  cuda = lmp->cuda;
  if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

  if (narg != 6) error->all(FLERR,"Illegal fix setforce/cuda command");

  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;

  flagx = flagy = flagz = 1;
  if (strcmp(arg[3],"NULL") == 0) flagx = 0;
  else xvalue = force->numeric(FLERR,arg[3]);
  if (strcmp(arg[4],"NULL") == 0) flagy = 0;
  else yvalue = force->numeric(FLERR,arg[4]);
  if (strcmp(arg[5],"NULL") == 0) flagz = 0;
  else zvalue = force->numeric(FLERR,arg[5]);

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  cu_foriginal=NULL;
}

/* ---------------------------------------------------------------------- */

int FixSetForceCuda::setmask()
{
  int mask = 0;
  mask |= POST_FORCE_CUDA;
  mask |= THERMO_ENERGY_CUDA;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE_CUDA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSetForceCuda::init()
{
  if(not cu_foriginal)
  cu_foriginal = new cCudaData<double, F_FLOAT, x> (foriginal,3);
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixSetForceCuda::setup(int vflag)
{
  MYDBG( printf("# CUDA: FixSetForceCuda::setup\n"); )

  if (strstr(update->integrate_style,"verlet"))
  {
    Cuda_FixSetForceCuda_Init(&cuda->shared_data);
    cuda->cu_f->upload();
    post_force(vflag);
    cuda->cu_f->download();

  }
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    cuda->cu_f->download();
    post_force_respa(vflag,nlevels_respa-1,0);
    cuda->cu_f->upload();
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
  MYDBG( printf("# CUDA: FixSetForceCuda::setup done\n"); )
}

/* ---------------------------------------------------------------------- */

void FixSetForceCuda::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSetForceCuda::post_force(int vflag)
{
  MYDBG( printf("# CUDA: FixSetForceCuda::postforce start\n"); )
  force_flag = 0;
  cu_foriginal->memset_device(0);
  Cuda_FixSetForceCuda_PostForce(&cuda->shared_data, groupbit, xvalue, yvalue,zvalue,(F_FLOAT*) cu_foriginal->dev_data(),flagx,flagy,flagz);
  cu_foriginal->download();
}

/* ---------------------------------------------------------------------- */

void FixSetForceCuda::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
  else {
          cuda->cu_f->download();
          cuda->cu_mask->download();

    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
    force_flag = 0;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        foriginal[0] += f[i][0];
        foriginal[1] += f[i][1];
        foriginal[2] += f[i][2];
        if (flagx) f[i][0] = 0.0;
        if (flagy) f[i][1] = 0.0;
        if (flagz) f[i][2] = 0.0;
      }
          cuda->cu_f->upload();
  }
}

/* ---------------------------------------------------------------------- */

void FixSetForceCuda::min_post_force(int vflag)
{
  post_force(vflag);
}


/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixSetForceCuda::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,3,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n+1];
}
