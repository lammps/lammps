// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include "fix_langevin_spin.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "modify.h"
#include "random_mars.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixLangevinSpin::FixLangevinSpin(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), random(nullptr)
{
  if (narg != 6) error->all(FLERR,"Illegal langevin/spin command");

  temp = utils::numeric(FLERR,arg[3],false,lmp);
  alpha_t = utils::numeric(FLERR,arg[4],false,lmp);
  seed = utils::inumeric(FLERR,arg[5],false,lmp);

  if (alpha_t < 0.0) {
    error->all(FLERR,"Illegal langevin/spin command");
  } else if (alpha_t == 0.0) {
    tdamp_flag = 0;
  } else {
    tdamp_flag = 1;
  }

  if (temp < 0.0) {
    error->all(FLERR,"Illegal langevin/spin command");
  } else if (temp == 0.0) {
    temp_flag = 0;
  } else {
    temp_flag = 1;
  }

  // initialize Marsaglia RNG with processor-unique seed

  // random = new RanPark(lmp,seed + comm->me);
  random = new RanMars(lmp,seed + comm->me);
}

/* ---------------------------------------------------------------------- */

FixLangevinSpin::~FixLangevinSpin()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

int FixLangevinSpin::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::init()
{
  // fix_langevin_spin has to be the last defined fix

  int flag_force = 0;
  int flag_lang = 0;
  for (int i = 0; i < modify->nfix; i++) {
     if (strcmp("precession/spin",modify->fix[i]->style)==0) flag_force = MAX(flag_force,i);
     if (strcmp("langevin/spin",modify->fix[i]->style)==0) flag_lang = i;
  }
  if (flag_force >= flag_lang) error->all(FLERR,"Fix langevin/spin has to come after all other spin fixes");

  gil_factor = 1.0/(1.0+(alpha_t)*(alpha_t));
  dts = 0.25 * update->dt;

  double hbar = force->hplanck/MY_2PI;  // eV/(rad.THz)
  double kb = force->boltz;             // eV/K

  D = (alpha_t*gil_factor*kb*temp);
  D /= (hbar*dts);
  sigma = sqrt(2.0*D);
}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^respa")) {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(nlevels_respa-1);
  } else post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::add_tdamping(double spi[3], double fmi[3])
{
  double cpx = fmi[1]*spi[2] - fmi[2]*spi[1];
  double cpy = fmi[2]*spi[0] - fmi[0]*spi[2];
  double cpz = fmi[0]*spi[1] - fmi[1]*spi[0];

  // adding the transverse damping

  fmi[0] -= alpha_t*cpx;
  fmi[1] -= alpha_t*cpy;
  fmi[2] -= alpha_t*cpz;
}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::add_temperature(double fmi[3])
{
  double rx = sigma*random->gaussian();
  double ry = sigma*random->gaussian();
  double rz = sigma*random->gaussian();

  // adding the random field

  fmi[0] += rx;
  fmi[1] += ry;
  fmi[2] += rz;

  // adding gilbert's prefactor

  fmi[0] *= gil_factor;
  fmi[1] *= gil_factor;
  fmi[2] *= gil_factor;
}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::compute_single_langevin(int i, double spi[3], double fmi[3])
{
  int *mask = atom->mask;
  if (mask[i] & groupbit) {
    if (tdamp_flag) add_tdamping(spi,fmi);
    if (temp_flag) add_temperature(fmi);
  }
}
