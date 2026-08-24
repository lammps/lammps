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
   Contributing author: AUTHOR_NAME_TBD (AFFILIATION_TBD)

   Velocity-Verlet integration of inertial spin dynamics.  Each spin obeys
   the Newtonian equation of motion

       s_mass d^2 s / dt^2 = fm

   for the un-normalized spin vector s = sp[3] * (sp[0], sp[1], sp[2]).

   The force is assembled from two contributions.  The SPIN pair styles and
   fix precession/spin express their interaction through the precession field
   fm, an angular frequency in rad.THz whose associated energy depends on the
   spin direction only.  The corresponding force on the full spin vector is
   therefore transverse and inversely proportional to the modulus,

       F_perp = (hbar/|S|) [ fm - (fm.s^) s^ ]

   Styles whose energy is a genuine function of the full spin vector, such as
   fix spring/tspin and fix langevin/tspin, instead accumulate their force
   directly into f_spin, in energy units, and it is used as it is.  After
   every update
   the spin direction sp[0..2] and the spin modulus sp[3] are recomputed, so
   the modulus is a dynamical degree of freedom.  This is in contrast to
   fix nve/spin, which integrates the fixed-modulus Landau-Lifshitz
   precession of a unit vector.
------------------------------------------------------------------------- */

#include "fix_nve_tspin.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "respa.h"
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
    Fix(lmp, narg, arg), step_respa(nullptr)
{
  if (!atom->tsp_flag)
    error->all(FLERR, "Fix {} requires atom style tspin", style);

  dynamic_group_allow = 1;
  time_integrate = 1;

  lattice_flag = 1;
  spin_flag = 1;
  spinmass_flag = 0;
  spinmass = 1.0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "lattice") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, std::string("fix ") + style + " lattice", error);
      lattice_flag = dof_moving(iarg + 1, arg, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "spin") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, std::string("fix ") + style + " spin", error);
      spin_flag = dof_moving(iarg + 1, arg, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "spinmass") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, std::string("fix ") + style + " spinmass", error);
      spinmass = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (spinmass <= 0.0)
        error->all(FLERR, iarg + 1, "Fix {} spinmass value must be > 0", style);
      spinmass_flag = 1;
      iarg += 2;

    } else {
      error->all(FLERR, iarg, "Unknown fix {} keyword: {}", style, arg[iarg]);
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixNVETSpin::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVETSpin::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  // fm is in rad.THz; hbar converts it into an energy gradient

  dtf_spin = dtf * force->hplanck / MY_2PI;

  if (utils::strmatch(update->integrate_style, "^respa"))
    step_respa = (dynamic_cast<Respa *>(update->integrate))->step;
}

/* ----------------------------------------------------------------------
   assign the spin mass as a multiple of the atomic mass of the atom type
------------------------------------------------------------------------- */

void FixNVETSpin::set_spin_mass()
{
  double **sp = atom->sp;
  double *s_mass = atom->s_mass;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (sp[i][3] <= 0.0) continue;
    s_mass[i] = spinmass * (rmass ? rmass[i] : mass[type[i]]);
  }
}

/* ---------------------------------------------------------------------- */

void FixNVETSpin::setup(int vflag)
{
  if (spinmass_flag) set_spin_mass();

  // every magnetic atom in the group must have a spin mass, otherwise its
  // spin would be silently left out of the integration

  double **sp = atom->sp;
  double *s_mass = atom->s_mass;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if ((mask[i] & groupbit) && (sp[i][3] > 0.0) && (s_mass[i] <= 0.0)) flag = 1;

  int flagall = 0;
  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, world);
  if (flagall && spin_flag)
    error->all(FLERR,
               "Fix {} requires a spin mass for every magnetic atom in the group; "
               "use the spinmass keyword or the velocity/tspin command", style);

  Fix::setup(vflag);
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVETSpin::initial_integrate(int /*vflag*/)
{
  double dtfm;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // update v and x of atoms in group

  if (lattice_flag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / (rmass ? rmass[i] : mass[type[i]]);
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  }

  // update v_s and the spin of atoms in group

  if (spin_flag) {
    double **sp = atom->sp;
    double **s_dot = atom->v_s;
    double **fm = atom->fm;
    double **f_spin = atom->f_spin;
    double *s_mass = atom->s_mass;

    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      if ((sp[i][3] <= 0.0) || (s_mass[i] <= 0.0)) continue;

      const double dtfsm = dtf / s_mass[i];
      const double dtfsp = dtf_spin / (s_mass[i] * sp[i][3]);
      const double fdots = fm[i][0] * sp[i][0] + fm[i][1] * sp[i][1] + fm[i][2] * sp[i][2];

      s_dot[i][0] += dtfsp * (fm[i][0] - fdots * sp[i][0]) + dtfsm * f_spin[i][0];
      s_dot[i][1] += dtfsp * (fm[i][1] - fdots * sp[i][1]) + dtfsm * f_spin[i][1];
      s_dot[i][2] += dtfsp * (fm[i][2] - fdots * sp[i][2]) + dtfsm * f_spin[i][2];

      // rebuild the un-normalized spin, advance it, then split it again
      // into a direction and a modulus

      double s[3];
      s[0] = sp[i][3] * sp[i][0] + dtv * s_dot[i][0];
      s[1] = sp[i][3] * sp[i][1] + dtv * s_dot[i][1];
      s[2] = sp[i][3] * sp[i][2] + dtv * s_dot[i][2];

      const double smag = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
      if (smag > 0.0) {
        sp[i][0] = s[0] / smag;
        sp[i][1] = s[1] / smag;
        sp[i][2] = s[2] / smag;
        sp[i][3] = smag;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVETSpin::final_integrate()
{
  double dtfm;

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // update v of atoms in group

  if (lattice_flag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / (rmass ? rmass[i] : mass[type[i]]);
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }
  }

  // update v_s of atoms in group

  if (spin_flag) {
    double **sp = atom->sp;
    double **s_dot = atom->v_s;
    double **fm = atom->fm;
    double **f_spin = atom->f_spin;
    double *s_mass = atom->s_mass;

    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      if ((sp[i][3] <= 0.0) || (s_mass[i] <= 0.0)) continue;

      const double dtfsm = dtf / s_mass[i];
      const double dtfsp = dtf_spin / (s_mass[i] * sp[i][3]);
      const double fdots = fm[i][0] * sp[i][0] + fm[i][1] * sp[i][1] + fm[i][2] * sp[i][2];

      s_dot[i][0] += dtfsp * (fm[i][0] - fdots * sp[i][0]) + dtfsm * f_spin[i][0];
      s_dot[i][1] += dtfsp * (fm[i][1] - fdots * sp[i][1]) + dtfsm * f_spin[i][1];
      s_dot[i][2] += dtfsp * (fm[i][2] - fdots * sp[i][2]) + dtfsm * f_spin[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVETSpin::initial_integrate_respa(int vflag, int ilevel, int /*iloop*/)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dtf_spin = dtf * force->hplanck / MY_2PI;

  // innermost level - update of v, x, v_s and the spins
  // all other levels - update of v and v_s

  if (ilevel == 0)
    initial_integrate(vflag);
  else
    final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVETSpin::final_integrate_respa(int ilevel, int /*iloop*/)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dtf_spin = dtf * force->hplanck / MY_2PI;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVETSpin::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dtf_spin = dtf * force->hplanck / MY_2PI;
}
