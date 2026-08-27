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

/* ------------------------------------------------------------------------
   Contributing author: Zhengtao Huang (The University of Hong Kong)
                        hzt990224@gmail.com

   Velocity-Verlet integration of inertial spin dynamics on atom_style spin.

   The spin of atom i is the full three-component vector

       S_i = sp[i][3] * (sp[i][0], sp[i][1], sp[i][2])          in muB

   whose modulus is a dynamical degree of freedom, and it obeys the Newtonian
   equation of motion

       m_s d^2 S_i / dt^2 = F_i,      F_i = -dU/dS_i            in eV/muB

   with a spin mass m_s and a spin velocity dS_i/dt.  This is a different
   physical model from fix nve/spin, which integrates the fixed-modulus
   Landau-Lifshitz precession of a unit vector.

   The magnetic interaction is read from the existing per-atom array fm in
   rad.THz.  A TSPIN-compatible interaction encodes the full magnetic force as

       fm_i = (|S_i|/hbar) F_i,       F_i = -dU/dS_i

   so the force needed by the equation of motion is recovered with

       F_i = (hbar / |S_i|) fm_i

   with all three components kept.  This is a stronger contract than the torque
   convention used by fixed-modulus Landau-Lifshitz styles, where a component
   parallel to the spin drops out of fm x S and need not represent a physical
   radial force.
------------------------------------------------------------------------- */

#include "fix_nve_tspin.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "modify.h"
#include "tspin.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_2PI;

/* ----------------------------------------------------------------------
   parse a degree-of-freedom keyword value
   accepts the (moving,frozen) wording of fix nve/spin as well as booleans
------------------------------------------------------------------------- */

static int dof_moving(int iarg, char **arg, LAMMPS *lmp)
{
  const std::string value = arg[iarg];
  if ((value == "moving") || (value == "mobile")) return 1;
  if (value == "frozen") return 0;
  return utils::logical(FLERR, arg[iarg], false, lmp) ? 1 : 0;
}

/* ---------------------------------------------------------------------- */

FixNVETSpin::FixNVETSpin(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  if (!atom->sp_flag) error->all(FLERR, "Fix {} requires atom style spin", style);

  time_integrate = 1;

  lattice_flag = 1;
  spin_flag = 1;
  spinmass_flag = 0;
  spinmass = 1.0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "lattice") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, std::string("fix ") + style + " lattice", error);
      lattice_flag = dof_moving(iarg + 1, arg, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "spin") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, std::string("fix ") + style + " spin", error);
      spin_flag = dof_moving(iarg + 1, arg, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "spinmass") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, std::string("fix ") + style + " spinmass", error);
      spinmass = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (spinmass <= 0.0)
        error->all(FLERR, iarg + 1, "Fix {} spinmass value must be > 0", style);
      spinmass_flag = 1;
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown fix {} keyword: {}", style, arg[iarg]);
    }
  }

  // spin velocity and spin mass live in custom per-atom properties, so that
  // they migrate and are stored in restart files without a new atom style

}

/* ----------------------------------------------------------------------
   the spin velocity and the spin mass live in custom per-atom properties, so
   that they migrate with the atoms and are written to restart files without a
   new atom style.  Modify::add_fix() may only be called from here, once this
   fix is fully registered.
------------------------------------------------------------------------- */

void FixNVETSpin::post_constructor()
{
  tspin_create_state(modify, error);
  tspin_bind_state(atom, error, std::string("Fix ") + style, index_vs, index_sm);

  // assign the spin masses right away, so that the velocity/tspin command can
  // be used before the first run

  if (spinmass_flag) set_spin_mass();
}

/* ---------------------------------------------------------------------- */

int FixNVETSpin::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVETSpin::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  hbar = force->hplanck / MY_2PI;

  // the state fix can be removed between runs, so bind to it again here

  tspin_bind_state(atom, error, std::string("Fix ") + style, index_vs, index_sm);

  // the magnetic force is not stored per rRESPA level, so the spin update
  // cannot be split over the levels

  if (utils::strmatch(update->integrate_style, "^respa"))
    error->all(FLERR, "Fix {} is not compatible with run_style respa", style);
}

/* ----------------------------------------------------------------------
   assign the spin mass as a multiple of the atomic mass of the atom type
------------------------------------------------------------------------- */

void FixNVETSpin::set_spin_mass()
{
  double **sp = atom->sp;
  double *s_mass = atom->dvector[index_sm];
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if ((sp[i][3] <= TSPIN_EPS) && (s_mass[i] <= 0.0)) continue;
    s_mass[i] = spinmass * (rmass ? rmass[i] : mass[type[i]]);
  }
}

/* ---------------------------------------------------------------------- */

void FixNVETSpin::setup(int vflag)
{
  if (spinmass_flag) set_spin_mass();
  check_spin_mass();
  Fix::setup(vflag);
}

/* ----------------------------------------------------------------------
   every magnetic atom in the group must have a spin mass by now
------------------------------------------------------------------------- */

void FixNVETSpin::check_spin_mass()
{
  if (!spin_flag) return;

  double **sp = atom->sp;
  double *s_mass = atom->dvector[index_sm];
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if ((mask[i] & groupbit) && (sp[i][3] > TSPIN_EPS) && (s_mass[i] <= 0.0)) flag = 1;

  int flagall = 0;
  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, world);
  if (flagall)
    error->all(FLERR,
               "Fix {} requires a spin mass for every magnetic atom in the group; use the "
               "spinmass keyword or the velocity/tspin command", style);
}

/* ----------------------------------------------------------------------
   half-step update of the spin velocities
------------------------------------------------------------------------- */

void FixNVETSpin::spin_kick()
{
  double **sp = atom->sp;
  double **fm = atom->fm;
  double **s_dot = atom->darray[index_vs];
  double *s_mass = atom->dvector[index_sm];
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (s_mass[i] <= 0.0) continue;

    // The fm encoding is singular only at exactly zero modulus.  Do not turn
    // a small modulus into an absorbing state: omit the kick at the isolated
    // zero crossing and let the drift carry the full spin vector through it.

    if (sp[i][3] == 0.0) continue;

    // A TSPIN-compatible fm encodes the full spin force in rad.THz;
    // hbar/|S| converts it back to eV/muB, radial component included.

    const double dtfsp = dtf * hbar / (s_mass[i] * sp[i][3]);

    s_dot[i][0] += dtfsp * fm[i][0];
    s_dot[i][1] += dtfsp * fm[i][1];
    s_dot[i][2] += dtfsp * fm[i][2];
  }
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVETSpin::initial_integrate(int /*vflag*/)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (lattice_flag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        const double dtfm = dtf / (rmass ? rmass[i] : mass[type[i]]);
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  }

  if (!spin_flag) return;

  spin_kick();

  // rebuild the un-normalized spin, advance it, then split it again into a
  // direction and a modulus, so that sp keeps its atom_style spin meaning

  double **sp = atom->sp;
  double **s_dot = atom->darray[index_vs];
  double *s_mass = atom->dvector[index_sm];

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (s_mass[i] <= 0.0) continue;

    double s[3];
    s[0] = sp[i][3] * sp[i][0] + dtv * s_dot[i][0];
    s[1] = sp[i][3] * sp[i][1] + dtv * s_dot[i][1];
    s[2] = sp[i][3] * sp[i][2] + dtv * s_dot[i][2];

    const double smag = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
    sp[i][3] = smag;
    if (smag > 0.0) {
      sp[i][0] = s[0] / smag;
      sp[i][1] = s[1] / smag;
      sp[i][2] = s[2] / smag;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVETSpin::final_integrate()
{
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (lattice_flag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        const double dtfm = dtf / (rmass ? rmass[i] : mass[type[i]]);
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }
  }

  if (spin_flag) spin_kick();
}

/* ---------------------------------------------------------------------- */

void FixNVETSpin::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  hbar = force->hplanck / MY_2PI;
}
