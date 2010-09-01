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

/* ----------------------------------------------------------------------
   Contributing author: Andres Jaramillo-Botero
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "fix_temp_rescale_eff.h"
#include "atom.h"
#include "force.h"
#include "group.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "modify.h"
#include "compute.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NOBIAS,BIAS};

/* ---------------------------------------------------------------------- */

FixTempRescaleEff::FixTempRescaleEff(LAMMPS *lmp, int narg, char **arg) :
  FixTempRescale(lmp, narg, arg)
{
  // create a new compute temp/eff, wiping out one parent class just created
  // id = fix-ID + temp, compute group = fix group

  modify->delete_compute(id_temp);

  char **newarg = new char*[6];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "temp/eff";
  modify->add_compute(3,newarg);
  delete [] newarg;
}

/* ---------------------------------------------------------------------- */

void FixTempRescaleEff::end_of_step()
{
  double t_current = temperature->compute_scalar();
  if (t_current == 0.0)
    error->all("Computed temperature for fix temp/rescale/eff cannot be 0.0");

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  double t_target = t_start + delta * (t_stop-t_start);

  // rescale velocity of appropriate atoms if outside window
  // for BIAS:
  //   temperature is current, so do not need to re-compute
  //   OK to not test returned v = 0, since factor is multiplied by v

  if (fabs(t_current-t_target) > t_window) {
    t_target = t_current - fraction*(t_current-t_target);
    double factor = sqrt(t_target/t_current);
    double efactor = 0.5 * force->boltz * temperature->dof;

    double **v = atom->v;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int *spin = atom->spin;
    double *ervel = atom->ervel;

    if (which == NOBIAS) {
      energy += (t_current-t_target) * efactor;
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  v[i][0] *= factor;
	  v[i][1] *= factor;
	  v[i][2] *= factor;
          if (spin[i]) 
            ervel[i] *= factor;
	}
      }
    } else {
      energy += (t_current-t_target) * efactor;
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  temperature->remove_bias(i,v[i]);
	  v[i][0] *= factor;
	  v[i][1] *= factor;
	  v[i][2] *= factor;
          if (spin[i])
            ervel[i] *= factor;          
	  temperature->restore_bias(i,v[i]);
	}
      }
    }

  }
}
