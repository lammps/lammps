/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "fix_nvt_sllod_omp.h"
#include <cmath>
#include <cstring>
#include "math_extra.h"
#include "atom.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "fix_deform.h"
#include "compute.h"
#include "error.h"
#include "domain.h"
#include "timer.h"

using namespace LAMMPS_NS;
using namespace FixConst;

typedef struct { double x,y,z; } dbl3_t;

/* ---------------------------------------------------------------------- */

FixNVTSllodOMP::FixNVTSllodOMP(LAMMPS *lmp, int narg, char **arg) :
  FixNHOMP(lmp, narg, arg)
{
  if (!tstat_flag)
    error->all(FLERR,"Temperature control must be used with fix nvt/sllod");
  if (pstat_flag)
    error->all(FLERR,"Pressure control can not be used with fix nvt/sllod");

  // default values

  if (mtchain_default_flag) mtchain = 1;


  // create a new compute temp style
  // id = fix-ID + temp

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");

  char **newarg = new char*[3];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "temp/deform";

  modify->add_compute(3,newarg);
  delete [] newarg;
  tcomputeflag = 1;
}

/* ---------------------------------------------------------------------- */

void FixNVTSllodOMP::init()
{
  FixNHOMP::init();

  if (!temperature->tempbias)
    error->all(FLERR,"Temperature for fix nvt/sllod/omp does not have a bias");

  nondeformbias = 0;
  if (strcmp(temperature->style,"temp/deform") != 0) nondeformbias = 1;

  // check fix deform remap settings

  int i;
  for (i = 0; i < modify->nfix; i++)
    if (strncmp(modify->fix[i]->style,"deform",6) == 0) {
      if (((FixDeform *) modify->fix[i])->remapflag != Domain::V_REMAP)
        error->all(FLERR,"Using fix nvt/sllod/omp with inconsistent fix "
                   "deform remap option");
      break;
    }
  if (i == modify->nfix)
    error->all(FLERR,"Using fix nvt/sllod/omp with no fix deform defined");
}

/* ----------------------------------------------------------------------
   perform half-step scaling of velocities
-----------------------------------------------------------------------*/

void FixNVTSllodOMP::nh_v_temp()
{
  // remove and restore bias = streaming velocity = Hrate*lamda + Hratelo
  // thermostat thermal velocity only
  // vdelu = SLLOD correction = Hrate*Hinv*vthermal
  // for non temp/deform BIAS:
  //   calculate temperature since some computes require temp
  //   computed on current nlocal atoms to remove bias

  dbl3_t * _noalias const v = (dbl3_t *) atom->v[0];
  const int * _noalias const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;
  int i;

  if (nondeformbias) temperature->compute_scalar();

  double h_two[6];
  MathExtra::multiply_shape_shape(domain->h_rate,domain->h_inv,h_two);

#if defined(_OPENMP)
#pragma omp parallel for default(none) private(i) shared(h_two) schedule(static)
#endif
  for (i = 0; i < nlocal; i++) {
    double vdelu0,vdelu1,vdelu2,buf[3];
    if (mask[i] & groupbit) {
      vdelu0 = h_two[0]*v[i].x + h_two[5]*v[i].y + h_two[4]*v[i].z;
      vdelu1 = h_two[1]*v[i].y + h_two[3]*v[i].z;
      vdelu2 = h_two[2]*v[i].z;
      temperature->remove_bias_thr(i,&v[i].x,buf);
      v[i].x = v[i].x*factor_eta - dthalf*vdelu0;
      v[i].y = v[i].y*factor_eta - dthalf*vdelu1;
      v[i].z = v[i].z*factor_eta - dthalf*vdelu2;
      temperature->restore_bias_thr(i,&v[i].x,buf);
    }
  }
}
