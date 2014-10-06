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
#include "fix_freeze_cuda.h"
#include "fix_freeze_cuda_cu.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "user_cuda.h"
#include "memory.h"
#include "modify.h"
#include "cuda_modify_flags.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace FixConstCuda;

/* ---------------------------------------------------------------------- */

FixFreezeCuda::FixFreezeCuda(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  cuda = lmp->cuda;
  if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

  if (narg != 3) error->all(FLERR,"Illegal fix freeze command");

  if (!atom->torque_flag)
    error->all(FLERR,"Fix freeze requires atom attribute torque");

  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;



  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  cu_foriginal=NULL;
}

/* ---------------------------------------------------------------------- */

int FixFreezeCuda::setmask()
{
  int mask = 0;
  mask |= POST_FORCE_CUDA;
  mask |= THERMO_ENERGY_CUDA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFreezeCuda::init()
{
  if(not cu_foriginal)
  cu_foriginal = new cCudaData<double, F_FLOAT, x> (foriginal,3);
  int count = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"freeze") == 0) count++;
  if (count > 1) error->all(FLERR,"More than one fix freeze");
}

/* ---------------------------------------------------------------------- */

void FixFreezeCuda::setup(int vflag)
{
  MYDBG( printf("# CUDA: FixFreezeCuda::setup\n"); )

  if (strstr(update->integrate_style,"verlet"))
  {
    Cuda_FixFreezeCuda_Init(&cuda->shared_data);
    cuda->cu_f->upload();
    post_force(vflag);
    cuda->cu_f->download();

  }

  MYDBG( printf("# CUDA: FixFreezeCuda::setup done\n"); )
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void FixFreezeCuda::post_force(int vflag)
{
  MYDBG( printf("# CUDA: FixFreezeCuda::postforce start\n"); )
  force_flag = 0;
  cu_foriginal->memset_device(0);
  Cuda_FixFreezeCuda_PostForce(&cuda->shared_data, groupbit, (F_FLOAT*) cu_foriginal->dev_data());
  cu_foriginal->download();
}

/* ---------------------------------------------------------------------- */



/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixFreezeCuda::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,3,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n+1];
}
