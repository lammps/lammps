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
#include "fix_addforce_cuda.h"
#include "fix_addforce_cuda_cu.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "domain.h"
#include "cuda.h"
#include "memory.h"
#include "cuda_modify_flags.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace FixConstCuda;

/* ---------------------------------------------------------------------- */

FixAddForceCuda::FixAddForceCuda(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

  if (narg < 6) error->all(FLERR,"Illegal fix addforce/cuda command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  xvalue = atof(arg[3]);
  yvalue = atof(arg[4]);
  zvalue = atof(arg[5]);

  // optional args

  iregion = -1;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix addforce/cuda command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1) error->all(FLERR,"Fix addforce/cuda region ID does not exist");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix addforce/cuda command");
  }
  
  if(iregion!=-1) error->all(FLERR,"Error: fix addforce/cuda does not currently support 'region' option");
  
  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  cu_foriginal = NULL;
}

/* ---------------------------------------------------------------------- */

int FixAddForceCuda::setmask()
{
  int mask = 0;
  mask |= POST_FORCE_CUDA;
  mask |= THERMO_ENERGY_CUDA;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE_CUDA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddForceCuda::init()
{
  if(not cu_foriginal)
  cu_foriginal = new cCudaData<double, F_FLOAT, x> (foriginal,4);    
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixAddForceCuda::setup(int vflag)
{
  MYDBG( printf("# CUDA: FixAddForceCuda::setup\n"); )
	
  if (strstr(update->integrate_style,"verlet"))
  {
    Cuda_FixAddForceCuda_Init(&cuda->shared_data);
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
  MYDBG( printf("# CUDA: FixAddForceCuda::setup done\n"); )
}

/* ---------------------------------------------------------------------- */

void FixAddForceCuda::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddForceCuda::post_force(int vflag)
{
  MYDBG( printf("# CUDA: FixAddForceCuda::postforce start\n"); )
  force_flag = 0;
  cu_foriginal->memset_device(0);
  Cuda_FixAddForceCuda_PostForce(&cuda->shared_data, groupbit, xvalue, yvalue,zvalue,(F_FLOAT*) cu_foriginal->dev_data());
  cu_foriginal->download();
}

/* ---------------------------------------------------------------------- */

void FixAddForceCuda::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddForceCuda::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixAddForceCuda::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[0];
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixAddForceCuda::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n+1];
}
