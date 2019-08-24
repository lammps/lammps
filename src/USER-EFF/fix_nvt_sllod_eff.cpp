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

#include <cmath>
#include <cstring>
#include <cstdlib>
#include "fix_nvt_sllod_eff.h"
#include "math_extra.h"
#include "atom.h"
#include "domain.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "fix_deform.h"
#include "compute.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVTSllodEff::FixNVTSllodEff(LAMMPS *lmp, int narg, char **arg) :
  FixNHEff(lmp, narg, arg)
{
  if (!tstat_flag)
    error->all(FLERR,"Temperature control must be used with fix nvt/sllod/eff");
  if (pstat_flag)
    error->all(FLERR,"Pressure control can not be used with fix nvt/sllod/eff");

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
  newarg[2] = (char *) "temp/deform/eff";

  modify->add_compute(3,newarg);
  delete [] newarg;
  tcomputeflag = 1;
}

/* ---------------------------------------------------------------------- */

void FixNVTSllodEff::init()
{
  FixNHEff::init();

  if (!temperature->tempbias)
    error->all(FLERR,"Temperature for fix nvt/sllod/eff does not have a bias");

  nondeformbias = 0;
  if (strcmp(temperature->style,"temp/deform/eff") != 0) nondeformbias = 1;

  // check fix deform remap settings

  int i;
  for (i = 0; i < modify->nfix; i++)
    if (strncmp(modify->fix[i]->style,"deform",6) == 0) {
      if (((FixDeform *) modify->fix[i])->remapflag != Domain::V_REMAP)
        error->all(FLERR,"Using fix nvt/sllod/eff with inconsistent fix deform "
                   "remap option");
      break;
    }
  if (i == modify->nfix)
    error->all(FLERR,"Using fix nvt/sllod/eff with no fix deform defined");
}


/* ----------------------------------------------------------------------
   perform half-step scaling of velocities
-----------------------------------------------------------------------*/

void FixNVTSllodEff::nh_v_temp()
{
  // remove and restore bias = streaming velocity = Hrate*lamda + Hratelo
  // thermostat thermal velocity only
  // vdelu = SLLOD correction = Hrate*Hinv*vthermal
  // for non temp/deform BIAS:
  //   calculate temperature since some computes require temp
  //   computed on current nlocal atoms to remove bias

  if (nondeformbias) temperature->compute_scalar();

  double **v = atom->v;
  double *ervel = atom->ervel;
  int *spin = atom->spin;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double h_two[6],vdelu[3];
  MathExtra::multiply_shape_shape(domain->h_rate,domain->h_inv,h_two);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      temperature->remove_bias(i,v[i]);
      vdelu[0] = h_two[0]*v[i][0] + h_two[5]*v[i][1] + h_two[4]*v[i][2];
      vdelu[1] = h_two[1]*v[i][1] + h_two[3]*v[i][2];
      vdelu[2] = h_two[2]*v[i][2];
      v[i][0] = v[i][0]*factor_eta - dthalf*vdelu[0];
      v[i][1] = v[i][1]*factor_eta - dthalf*vdelu[1];
      v[i][2] = v[i][2]*factor_eta - dthalf*vdelu[2];
      temperature->restore_bias(i,v[i]);
      if (abs(spin[i])==1)
        ervel[i] = ervel[i]*factor_eta -
          dthalf*sqrt(vdelu[0]*vdelu[0]+vdelu[1]*vdelu[1]+vdelu[2]*vdelu[2]);
    }
  }
}
