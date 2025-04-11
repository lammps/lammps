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

#include "fix_sph.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSPH::FixSPH(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if ((atom->esph_flag != 1) || (atom->rho_flag != 1) || (atom->vest_flag != 1))
    error->all(FLERR, "Fix sph requires atom attributes energy, density, and velocity estimates, e.g. in atom_style sph");

  if (narg != 3)
    error->all(FLERR,"Illegal number of arguments for fix sph command");

  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixSPH::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPH::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ---------------------------------------------------------------------- */

void FixSPH::setup_pre_force(int /*vflag*/)
{
  remap_v_flag = domain->deform_vremap;
  if (remap_v_flag && (!comm->ghost_velocity))
    error->all(FLERR, "Fix sph requires ghost atoms store velocity when deforming with remap v");

  // set vest equal to v
  double **v = atom->v;
  double **vest = atom->vest;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      vest[i][0] = v[i][0];
      vest[i][1] = v[i][1];
      vest[i][2] = v[i][2];
    }
  }
}

/* ----------------------------------------------------------------------
 allow for both per-type and per-atom mass
 ------------------------------------------------------------------------- */

void FixSPH::initial_integrate(int /*vflag*/)
{
  // update v and x and rho and e of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **vest = atom->vest;
  double *rho = atom->rho;
  double *drho = atom->drho;
  double *esph = atom->esph;
  double *desph = atom->desph;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int rmass_flag = atom->rmass_flag;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int i;
  double dtfm;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass_flag) {
        dtfm = dtf / rmass[i];
      } else {
        dtfm = dtf / mass[type[i]];
      }

      esph[i] += dtf * desph[i]; // half-step update of particle internal energy
      rho[i] += dtf * drho[i]; // ... and density

      // extrapolate velocity for use with velocity-dependent potentials, e.g. SPH
      // if velocities are remapped, perform this extrapolation after communication
      if (remap_v_flag) {
        vest[i][0] = dtfm * f[i][0];
        vest[i][1] = dtfm * f[i][1];
        vest[i][2] = dtfm * f[i][2];
      } else {
        vest[i][0] = v[i][0] + 2.0 * dtfm * f[i][0];
        vest[i][1] = v[i][1] + 2.0 * dtfm * f[i][1];
        vest[i][2] = v[i][2] + 2.0 * dtfm * f[i][2];
      }

      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSPH::pre_force(int /*vflag*/)
{
  // if velocities are remapped, calculate estimates here
  // note that vest currently stores dtfm * force
  if (!remap_v_flag) return;

  double **v = atom->v;
  double **vest = atom->vest;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  int nall = nlocal + atom->nghost;
  for (int i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      vest[i][0] += v[i][0];
      vest[i][1] += v[i][1];
      vest[i][2] += v[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSPH::final_integrate()
{
  // update v, rho, and e of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double *esph = atom->esph;
  double *desph = atom->desph;
  double *rho = atom->rho;
  double *drho = atom->drho;
  int *type = atom->type;
  int *mask = atom->mask;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;
  double dtfm;
  double *rmass = atom->rmass;
  int rmass_flag = atom->rmass_flag;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      if (rmass_flag) {
        dtfm = dtf / rmass[i];
      } else {
        dtfm = dtf / mass[type[i]];
      }
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      esph[i] += dtf * desph[i];
      rho[i] += dtf * drho[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSPH::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
