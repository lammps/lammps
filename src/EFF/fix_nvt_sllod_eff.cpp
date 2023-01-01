// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_nvt_sllod_eff.h"

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "fix_deform.h"
#include "group.h"
#include "math_extra.h"
#include "modify.h"

#include <cmath>
#include <cstring>

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

  psllod_flag = 0;
  if (mtchain_default_flag) mtchain = 1;

  // select SLLOD/p-SLLOD/g-SLLOD variant

  int iarg = 3;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"psllod") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix nvt/sllod/eff psllod", error);
      psllod_flag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else iarg++;
  }

  // create a new compute temp style
  // id = fix-ID + temp

  id_temp = utils::strdup(std::string(id) + "_temp");
  modify->add_compute(fmt::format("{} {} tmp/deform/eff",
                                  id_temp,group->names[igroup]));
  tcomputeflag = 1;
  nondeformbias = 0;
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
      if ((dynamic_cast<FixDeform *>(modify->fix[i]))->remapflag != Domain::V_REMAP)
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
      if (!psllod_flag) temperature->remove_bias(i,v[i]);
      vdelu[0] = h_two[0]*v[i][0] + h_two[5]*v[i][1] + h_two[4]*v[i][2];
      vdelu[1] = h_two[1]*v[i][1] + h_two[3]*v[i][2];
      vdelu[2] = h_two[2]*v[i][2];
      if (psllod_flag) temperature->remove_bias(i,v[i]);
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
