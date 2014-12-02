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
#include "domain.h"
#include "user_cuda.h"
#include "cuda_modify_flags.h"
#include "variable.h"
#include "input.h"
#include "modify.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace FixConstCuda;

enum{NONE,CONSTANT,EQUAL};

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

  xstr = ystr = zstr = NULL;
  xvalue = yvalue = zvalue = 0;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else if (strcmp(arg[3],"NULL") == 0) {
    xstyle = NONE;
  } else {
    xvalue = force->numeric(FLERR,arg[3]);
    xstyle = CONSTANT;
  }
  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[4][2]);
  } else if (strcmp(arg[4],"NULL") == 0) {
    ystyle = NONE;
  } else {
    yvalue = force->numeric(FLERR,arg[4]);
    ystyle = CONSTANT;
  }
  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[5][2]);
  } else if (strcmp(arg[5],"NULL") == 0) {
    zstyle = NONE;
  } else {
    zvalue = force->numeric(FLERR,arg[5]);
    zstyle = CONSTANT;
  }

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

FixAveForceCuda::~FixAveForceCuda()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
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
  cu_foriginal = new cCudaData<double, F_CFLOAT, x> (foriginal,4);

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix aveforce does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else error->all(FLERR,"Variable for fix aveforce is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix aveforce does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else error->all(FLERR,"Variable for fix aveforce is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix aveforce does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else error->all(FLERR,"Variable for fix aveforce is invalid style");
  }

  if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL) varflag = EQUAL;
  else varflag = CONSTANT;

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
  Cuda_FixAveForceCuda_PostForce_FOrg(&cuda->shared_data, groupbit,(F_CFLOAT*) cu_foriginal->dev_data());
  cu_foriginal->download();

  // average the force on participating atoms
  // add in requested amount

  MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
  int ncount = static_cast<int> (foriginal_all[3]);
  if (ncount == 0) return;

  if (varflag == EQUAL) {
    unsigned int datamask = EMPTY_MASK;
    if (xstyle == EQUAL) datamask &= input->variable->data_mask(xstr);
    if (ystyle == EQUAL) datamask &= input->variable->data_mask(ystr);
    if (zstyle == EQUAL) datamask &= input->variable->data_mask(zstr);

    cuda->download(datamask);
    modify->clearstep_compute();
    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    modify->addstep_compute(update->ntimestep + 1);
  }

  double fave[3];
  fave[0] = foriginal_all[0]/ncount + xvalue;
  fave[1] = foriginal_all[1]/ncount + yvalue;
  fave[2] = foriginal_all[2]/ncount + zvalue;

  // set force of all participating atoms to same value
  // only for active dimensions

  Cuda_FixAveForceCuda_PostForce_Set(&cuda->shared_data, groupbit,!(xstyle==NONE),!(ystyle==NONE),!(zstyle==NONE),fave[0],fave[1],fave[2]);
}

/* ---------------------------------------------------------------------- */

void FixAveForceCuda::post_force_respa(int vflag, int ilevel, int iloop)
{

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
