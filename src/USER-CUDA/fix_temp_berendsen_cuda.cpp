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
#include <cstdlib>
#include <cmath>
#include "fix_temp_berendsen_cuda.h"
#include "fix_temp_berendsen_cuda_cu.h"
#include "atom.h"
#include "force.h"
#include "group.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "compute.h"
#include "error.h"
#include "user_cuda.h"
#include "cuda_modify_flags.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace FixConstCuda;

enum{NOBIAS,BIAS};

/* ---------------------------------------------------------------------- */

FixTempBerendsenCuda::FixTempBerendsenCuda(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

  if (narg != 6) error->all(FLERR,"Illegal fix temp/berendsen/cuda command");

  // Berendsen thermostat should be applied every step

  nevery = 1;

  t_start = force->numeric(FLERR,arg[3]);
  t_stop = force->numeric(FLERR,arg[4]);
  t_period = force->numeric(FLERR,arg[5]);

  // error checks

  if (t_period <= 0.0) error->all(FLERR,"Fix temp/berendsen/cuda period must be > 0.0");

  // create a new compute temp style
  // id = fix-ID + temp, compute group = fix group

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");

  char **newarg = new char*[3];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "temp/cuda";
  modify->add_compute(3,newarg);
  delete [] newarg;
  tflag = 1;
}

/* ---------------------------------------------------------------------- */

FixTempBerendsenCuda::~FixTempBerendsenCuda()
{
  // delete temperature if fix created it

  if (tflag) modify->delete_compute(id_temp);
  delete [] id_temp;
}

/* ---------------------------------------------------------------------- */

int FixTempBerendsenCuda::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP_CUDA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTempBerendsenCuda::init()
{
  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR,"Temperature ID for fix temp/berendsen/cuda does not exist");
  temperature = modify->compute[icompute];
  if(not temperature->cudable)
        error->warning(FLERR,"Fix temp/berendsen/cuda uses non cudable temperature compute");
  if (temperature->tempbias) which = BIAS;
  else which = NOBIAS;

  //temperature->init();        //not in original berendsen possible error?
}

/* ---------------------------------------------------------------------- */

void FixTempBerendsenCuda::end_of_step()
{
  double t_current;
  if(not temperature->cudable) {cuda->cu_x->download();cuda->cu_v->download();}
  t_current = temperature->compute_scalar();
  if (t_current == 0.0)
    error->all(FLERR,"Computed temperature for fix temp/berendsen/cuda cannot be 0.0");

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  t_target = t_start + delta * (t_stop-t_start);

  // rescale velocities by lamda

  double lamda = sqrt(1.0 + update->dt/t_period*(t_target/t_current - 1.0));

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (which == NOBIAS) {
        Cuda_FixTempBerendsenCuda_EndOfStep(&cuda->shared_data, groupbit,lamda);

    } else {
      if(not temperature->cudable)
      {
              cuda->cu_x->download();cuda->cu_v->download();
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          temperature->remove_bias(i,v[i]);
           v[i][0] *= lamda;
          v[i][1] *= lamda;
          v[i][2] *= lamda;
          temperature->restore_bias(i,v[i]);
        }
        }
          cuda->cu_v->upload();
      }
      else
          {
              temperature->remove_bias_all();
            Cuda_FixTempBerendsenCuda_EndOfStep(&cuda->shared_data, groupbit,lamda);
            temperature->restore_bias_all();
          }
    }


}

/* ---------------------------------------------------------------------- */

int FixTempBerendsenCuda::modify_param(int narg, char **arg)
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
    return 2;
  }
  return 0;
}


/* ---------------------------------------------------------------------- */

void FixTempBerendsenCuda::reset_target(double t_new)
{
  t_start = t_stop = t_new;
}
