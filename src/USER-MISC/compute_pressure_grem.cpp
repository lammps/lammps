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

#include "compute_pressure_grem.h"
#include <cstring>
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "force.h"
#include "kspace.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   Last argument is the id of the gREM fix
------------------------------------------------------------------------- */

ComputePressureGrem::ComputePressureGrem(LAMMPS *lmp, int narg, char **arg) :
  ComputePressure(lmp, narg-1, arg)
{
  int len = strlen(arg[narg-1])+1;
  fix_grem = new char[len];
  strcpy(fix_grem,arg[narg-1]);
}

/* ---------------------------------------------------------------------- */

ComputePressureGrem::~ComputePressureGrem()
{
  delete [] fix_grem;
}
/* ---------------------------------------------------------------------- */

void ComputePressureGrem::init()
{
  ComputePressure::init();

  // Initialize hook to gREM fix
  int ifix = modify->find_fix(fix_grem);
  if (ifix < 0)
    error->all(FLERR,"Fix grem ID for compute PRESSURE/GREM does not exist");

  int dim;
  scale_grem = (double *)modify->fix[ifix]->extract("scale_grem",dim);

  if (scale_grem == NULL || dim != 0)
    error->all(FLERR,"Cannot extract gREM scale factor from fix grem");
}

/* ----------------------------------------------------------------------
   compute total pressure, averaged over Pxx, Pyy, Pzz
------------------------------------------------------------------------- */

double ComputePressureGrem::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if (update->vflag_global != invoked_scalar)
    error->all(FLERR,"Virial was not tallied on needed timestep");

  // invoke temperature if it hasn't been already

  double t;
  if (keflag) {
    if (temperature->invoked_scalar != update->ntimestep)
      t = temperature->compute_scalar() / (*scale_grem);
    else t = temperature->scalar / (*scale_grem);
  }

  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    virial_compute(3,3);
    if (keflag)
      scalar = (temperature->dof * boltz * t +
                virial[0] + virial[1] + virial[2]) / 3.0 * inv_volume * nktv2p;
    else
      scalar = (virial[0] + virial[1] + virial[2]) / 3.0 * inv_volume * nktv2p;
  } else {
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
    virial_compute(2,2);
    if (keflag)
      scalar = (temperature->dof * boltz * t +
                virial[0] + virial[1]) / 2.0 * inv_volume * nktv2p;
    else
      scalar = (virial[0] + virial[1]) / 2.0 * inv_volume * nktv2p;
  }

  return scalar;
}

/* ----------------------------------------------------------------------
   compute pressure tensor
   assume KE tensor has already been computed
------------------------------------------------------------------------- */

void ComputePressureGrem::compute_vector()
{
  invoked_vector = update->ntimestep;
  if (update->vflag_global != invoked_vector)
    error->all(FLERR,"Virial was not tallied on needed timestep");

  if (force->kspace && kspace_virial && force->kspace->scalar_pressure_flag)
    error->all(FLERR,"Must use 'kspace_modify pressure/scalar no' for "
               "tensor components with kspace_style msm");

  // invoke temperature if it hasn't been already

  double ke_tensor[6];
  if (keflag) {
    if (temperature->invoked_vector != update->ntimestep)
      temperature->compute_vector();
    for (int i = 0; i < 6; ++i)
      ke_tensor[i] = temperature->vector[i] / (*scale_grem);
  }

  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    virial_compute(6,3);
    if (keflag) {
      for (int i = 0; i < 6; i++)
        vector[i] = (ke_tensor[i] + virial[i]) * inv_volume * nktv2p;
    } else
      for (int i = 0; i < 6; i++)
        vector[i] = virial[i] * inv_volume * nktv2p;
  } else {
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
    virial_compute(4,2);
    if (keflag) {
      vector[0] = (ke_tensor[0] + virial[0]) * inv_volume * nktv2p;
      vector[1] = (ke_tensor[1] + virial[1]) * inv_volume * nktv2p;
      vector[3] = (ke_tensor[3] + virial[3]) * inv_volume * nktv2p;
      vector[2] = vector[4] = vector[5] = 0.0;
    } else {
      vector[0] = virial[0] * inv_volume * nktv2p;
      vector[1] = virial[1] * inv_volume * nktv2p;
      vector[3] = virial[3] * inv_volume * nktv2p;
      vector[2] = vector[4] = vector[5] = 0.0;
    }
  }
}

