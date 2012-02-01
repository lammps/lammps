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
#include <cmath>
#include "fix_temp_rescale_cuda.h"
#include "fix_temp_rescale_cuda_cu.h"
#include "atom.h"
#include "force.h"
#include "group.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "comm.h"
#include "modify.h"
#include "compute.h"
#include "error.h"
#include "cuda.h"
#include "cuda_modify_flags.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace FixConstCuda;

enum{NOBIAS,BIAS};

/* ---------------------------------------------------------------------- */

FixTempRescaleCuda::FixTempRescaleCuda(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

  if (narg < 8) error->all(FLERR,"Illegal fix temp/rescale/cuda command");

  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix temp/rescale/cuda command");

  scalar_flag = 1;
  global_freq = nevery;
  extscalar = 1;

  t_start = atof(arg[4]);
  t_stop = atof(arg[5]);
  t_window = atof(arg[6]);
  fraction = atof(arg[7]);

  // create a new compute temp
  // id = fix-ID + temp, compute group = fix group

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");

  char **newarg = new char*[6];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "temp/cuda";
  modify->add_compute(3,newarg);
  delete [] newarg;
  tflag = 1;

  energy = 0.0;
}

/* ---------------------------------------------------------------------- */

FixTempRescaleCuda::~FixTempRescaleCuda()
{
  // delete temperature if fix created it

  if (tflag) modify->delete_compute(id_temp);
  delete [] id_temp;
}

/* ---------------------------------------------------------------------- */

int FixTempRescaleCuda::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP_CUDA;
  mask |= THERMO_ENERGY_CUDA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTempRescaleCuda::init()
{
  int icompute = modify->find_compute(id_temp);
  if (icompute < 0) 
    error->all(FLERR,"Temperature ID for fix temp/rescale/cuda does not exist");
  temperature = modify->compute[icompute];
  if(not temperature->cudable) 
	error->warning(FLERR,"Fix temp/rescale/cuda uses non cudable temperature compute");
  if (temperature->tempbias) which = BIAS;
  else which = NOBIAS;
}

/* ---------------------------------------------------------------------- */

void FixTempRescaleCuda::end_of_step()
{
  double t_current;
  if(not temperature->cudable) {cuda->cu_x->download();cuda->cu_v->download();}
  t_current = temperature->compute_scalar();
  if (t_current == 0.0)
    error->all(FLERR,"Computed temperature for fix temp/rescale/cuda cannot be 0.0");

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  double t_target = t_start + delta * (t_stop-t_start);

  // rescale velocity of appropriate atoms if outside window

  if (fabs(t_current-t_target) > t_window) {
    t_target = t_current - fraction*(t_current-t_target);
    double factor = sqrt(t_target/t_current);
    double efactor = 0.5 * force->boltz * temperature->dof;

    double **v = atom->v;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    if (which == NOBIAS) {
      energy += (t_current-t_target) * efactor;

	Cuda_FixTempRescaleCuda_EndOfStep(&cuda->shared_data, groupbit,factor);

    } else if (which == BIAS) {
      energy += (t_current-t_target) * efactor;
      if(not temperature->cudable)
      {
      	cuda->cu_x->download();cuda->cu_v->download();
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  temperature->remove_bias(i,v[i]);
	  v[i][0] *= factor;
	  v[i][1] *= factor;
	  v[i][2] *= factor;
	  temperature->restore_bias(i,v[i]);
	}
        }
	  cuda->cu_v->upload();
      }
      else
      {
	    temperature->remove_bias_all();
	    Cuda_FixTempRescaleCuda_EndOfStep(&cuda->shared_data, groupbit,factor);
	    temperature->restore_bias_all();
      }
    }

  }
}

/* ---------------------------------------------------------------------- */

int FixTempRescaleCuda::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete [] id_temp;
    int n = strlen(arg[1]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0) error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,"Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning(FLERR,"Group for fix_modify temp != fix group");
    if(not temperature->cudable) 
	  error->warning(FLERR,"Fix temp/rescale/cuda uses non cudable temperature compute");
    return 2;
  }
  return 0;
}


/* ---------------------------------------------------------------------- */

void FixTempRescaleCuda::reset_target(double t_new)
{
  t_start = t_stop = t_new;
}

/* ---------------------------------------------------------------------- */

double FixTempRescaleCuda::compute_scalar()
{
  return energy;
}
