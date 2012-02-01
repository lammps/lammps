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


#include "mpi.h"
#include <cstring>
#include <cstdlib>
#include "fix_aveforce_cuda.h"
#include "fix_aveforce_cuda_cu.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "domain.h"
#include "cuda.h"
#include "cuda_modify_flags.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace FixConstCuda;

/* ---------------------------------------------------------------------- */

FixAveForceCuda::FixAveForceCuda(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

  if (narg != 6) error->all(FLERR,"Illegal fix aveforce command");

  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;

  xflag = yflag = zflag = 1;
  if (strcmp(arg[3],"NULL") == 0) xflag = 0;
  else xvalue = atof(arg[3]);
  if (strcmp(arg[4],"NULL") == 0) yflag = 0;
  else yvalue = atof(arg[4]);
  if (strcmp(arg[5],"NULL") == 0) zflag = 0;
  else zvalue = atof(arg[5]);

  // optional args

  iregion = -1;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix aveforce command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1) error->all(FLERR,"Fix aveforce region ID does not exist");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix aveforce command");

  }
  
  if(iregion!=-1) error->all(FLERR,"Error: fix aveforce/cuda does not currently support 'region' option");

  foriginal_all[0] = foriginal_all[1] = foriginal_all[2] = foriginal_all[3] = 0.0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  cu_foriginal = NULL;
  
}

/* ---------------------------------------------------------------------- */

int FixAveForceCuda::setmask()
{
  int mask = 0;
  mask |= POST_FORCE_CUDA;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE_CUDA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveForceCuda::init()
{
  if(not cu_foriginal)
  cu_foriginal = new cCudaData<double, F_FLOAT, x> (foriginal,4);    
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // ncount = total # of atoms in group

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
}

/* ---------------------------------------------------------------------- */

void FixAveForceCuda::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
  {
    Cuda_FixAveForceCuda_Init(&cuda->shared_data);
    cuda->cu_f->upload();
    post_force(vflag);
    cuda->cu_f->download();
    
  }
  else
  {
    cuda->cu_f->download();
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(vflag,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
    cuda->cu_f->upload();
  }
}

/* ---------------------------------------------------------------------- */

void FixAveForceCuda::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAveForceCuda::post_force(int vflag)
{
  // sum forces on participating atoms

  cu_foriginal->memset_device(0);
  Cuda_FixAveForceCuda_PostForce_FOrg(&cuda->shared_data, groupbit,(F_FLOAT*) cu_foriginal->dev_data());
  cu_foriginal->download();

  // average the force on participating atoms
  // add in requested amount

  MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
  int ncount = static_cast<int> (foriginal_all[3]);
  if (ncount == 0) return;
  double fave[3];
  fave[0] = foriginal_all[0]/ncount + xvalue;
  fave[1] = foriginal_all[1]/ncount + yvalue;
  fave[2] = foriginal_all[2]/ncount + zvalue;

  // set force of all participating atoms to same value
  // only for active dimensions

  Cuda_FixAveForceCuda_PostForce_Set(&cuda->shared_data, groupbit,xflag,yflag,zflag,fave[0],fave[1],fave[2]);
}

/* ---------------------------------------------------------------------- */

void FixAveForceCuda::post_force_respa(int vflag, int ilevel, int iloop)
{
  // ave + extra force on outermost level
  // just ave on inner levels
  if (ilevel == nlevels_respa-1) post_force(vflag);
  else {
    cuda->cu_f->download();
    cuda->cu_mask->download();
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double foriginal[4];
    foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	foriginal[0] += f[i][0];
	foriginal[1] += f[i][1];
	foriginal[2] += f[i][2];
	foriginal[3] += 1;
	
      }

    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    int ncount = static_cast<int> (foriginal_all[3]);
    if (ncount == 0) return;
    double fave[3];
    fave[0] = foriginal_all[0]/ncount;
    fave[1] = foriginal_all[1]/ncount;
    fave[2] = foriginal_all[2]/ncount;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	if (xflag) f[i][0] = fave[0];
	if (yflag) f[i][1] = fave[1];
	if (zflag) f[i][2] = fave[2];
      }
    cuda->cu_f->upload();
  }
}

/* ---------------------------------------------------------------------- */

void FixAveForceCuda::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixAveForceCuda::compute_vector(int n)
{
  return foriginal_all[n];
}
