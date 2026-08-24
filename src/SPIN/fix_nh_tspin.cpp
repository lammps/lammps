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

   Nose-Hoover machinery for inertial spin dynamics.  This is a base class
   and registers no fix style of its own; fix nvt/tspin, fix npt/tspin and
   fix nph/tspin derive from it.

   In addition to the lattice degrees of freedom handled by the parent
   class, the spins are integrated with velocity-Verlet and their spin
   velocities are coupled to a second Nose-Hoover chain.  That chain is
   independent of the lattice chain but is driven by the same target
   temperature and the same damping parameter, so that a spin-lattice
   simulation is thermostatted as a whole.
------------------------------------------------------------------------- */

#include "fix_nh_tspin.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathConst::MY_2PI;

/* ---------------------------------------------------------------------- */

FixNHTSpin::FixNHTSpin(LAMMPS *lmp, int narg, char **arg) :
    FixNH(lmp, narg, arg), lattice_flag(1), spin_flag(1), spinmass_flag(0), spinmass(1.0),
    dtf_spin(0.0), t_current_spin(0.0), ke_target_spin(0.0), tdof_spin(0.0), etas(nullptr),
    etas_dot(nullptr), etas_dotdot(nullptr), etas_mass(nullptr)
{
  if (!atom->tsp_flag) error->all(FLERR, "Fix {} requires atom style tspin", style);

  // parse the keywords that the parent class skipped

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "spinmass") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, std::string("fix ") + style + " spinmass", error);
      spinmass = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (spinmass <= 0.0) error->all(FLERR, iarg + 1, "Fix {} spinmass must be > 0", style);
      spinmass_flag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "lattice") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, std::string("fix ") + style + " lattice", error);
      const std::string value = arg[iarg + 1];
      if (value == "frozen")
        lattice_flag = 0;
      else if ((value == "moving") || (value == "mobile"))
        lattice_flag = 1;
      else
        lattice_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp) ? 1 : 0;
      iarg += 2;
    } else if (strcmp(arg[iarg], "spin") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, std::string("fix ") + style + " spin", error);
      const std::string value = arg[iarg + 1];
      if (value == "frozen")
        spin_flag = 0;
      else if ((value == "moving") || (value == "mobile"))
        spin_flag = 1;
      else
        spin_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp) ? 1 : 0;
      iarg += 2;
    } else {
      iarg++;    // every other keyword belongs to the parent class
    }
  }

  if (!spin_flag) error->all(FLERR, "Fix {} cannot be used with spin frozen", style);
  if (pstat_flag && !lattice_flag)
    error->all(FLERR, "Fix {} cannot use a barostat with lattice frozen", style);

  // the spin chain has the same length as the lattice chain
  // etas_dot has one extra element so the update can read etas_dot[ich+1]

  memory->create(etas, mtchain, "nh/tspin:etas");
  memory->create(etas_dot, mtchain + 1, "nh/tspin:etas_dot");
  memory->create(etas_dotdot, mtchain, "nh/tspin:etas_dotdot");
  memory->create(etas_mass, mtchain, "nh/tspin:etas_mass");

  for (int ich = 0; ich < mtchain; ich++) etas[ich] = etas_dot[ich] = etas_dotdot[ich] = 0.0;
  etas_dot[mtchain] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixNHTSpin::~FixNHTSpin()
{
  memory->destroy(etas);
  memory->destroy(etas_dot);
  memory->destroy(etas_dotdot);
  memory->destroy(etas_mass);
}

/* ---------------------------------------------------------------------- */

void FixNHTSpin::init()
{
  FixNH::init();
  reset_dt();
}

/* ---------------------------------------------------------------------- */

void FixNHTSpin::reset_dt()
{
  FixNH::reset_dt();

  // fm is in rad.THz; hbar converts it into an energy gradient

  dtf_spin = dtf * force->hplanck / MY_2PI;
}

/* ----------------------------------------------------------------------
   assign the spin mass as a multiple of the atomic mass of the atom type
------------------------------------------------------------------------- */

void FixNHTSpin::set_spin_mass()
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

void FixNHTSpin::setup(int vflag)
{
  if (spinmass_flag) set_spin_mass();

  double **sp = atom->sp;
  double *s_mass = atom->s_mass;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  // count the thermostatted spin degrees of freedom and check the spin masses

  int ndof = 0, flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit) || (sp[i][3] <= 0.0)) continue;
    if (s_mass[i] > 0.0)
      ndof += 3;
    else
      flag = 1;
  }

  int ndof_all = 0, flagall = 0;
  MPI_Allreduce(&ndof, &ndof_all, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, world);
  tdof_spin = ndof_all;

  if (flagall)
    error->all(FLERR,
               "Fix {} requires a spin mass for every magnetic atom in the group; "
               "use the spinmass keyword or the velocity/tspin command",
               style);
  if (tdof_spin == 0.0)
    error->all(FLERR, "Fix {} has no spin degrees of freedom to thermostat", style);

  FixNH::setup(vflag);

  // initialize the spin chain masses from the target temperature that the
  // parent class has just computed

  if (tstat_flag) {
    etas_mass[0] = tdof_spin * boltz * t_target / (t_freq * t_freq);
    for (int ich = 1; ich < mtchain; ich++)
      etas_mass[ich] = boltz * t_target / (t_freq * t_freq);
    for (int ich = 1; ich < mtchain; ich++)
      etas_dotdot[ich] =
          (etas_mass[ich - 1] * etas_dot[ich - 1] * etas_dot[ich - 1] - boltz * t_target) /
          etas_mass[ich];
  }
}

/* ----------------------------------------------------------------------
   instantaneous spin temperature from the spin kinetic energy
------------------------------------------------------------------------- */

double FixNHTSpin::compute_spin_temp()
{
  double **sp = atom->sp;
  double **s_dot = atom->v_s;
  double *s_mass = atom->s_mass;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  double ke = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if ((sp[i][3] <= 0.0) || (s_mass[i] <= 0.0)) continue;
    ke += s_mass[i] *
        (s_dot[i][0] * s_dot[i][0] + s_dot[i][1] * s_dot[i][1] + s_dot[i][2] * s_dot[i][2]);
  }

  double keall = 0.0;
  MPI_Allreduce(&ke, &keall, 1, MPI_DOUBLE, MPI_SUM, world);

  if (tdof_spin == 0.0) return 0.0;
  return keall * force->mvv2e / (tdof_spin * boltz);
}

/* ----------------------------------------------------------------------
   scale all thermostatted spin velocities by factor
------------------------------------------------------------------------- */

void FixNHTSpin::nh_vs_scale(double factor)
{
  double **sp = atom->sp;
  double **s_dot = atom->v_s;
  double *s_mass = atom->s_mass;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if ((sp[i][3] <= 0.0) || (s_mass[i] <= 0.0)) continue;
    s_dot[i][0] *= factor;
    s_dot[i][1] *= factor;
    s_dot[i][2] *= factor;
  }
}

/* ----------------------------------------------------------------------
   propagate the spin Nose-Hoover chain by dthalf
   same scheme as FixNH::nhc_temp_integrate, acting on the spin velocities
------------------------------------------------------------------------- */

void FixNHTSpin::nhc_spin_integrate()
{
  int ich;
  double expfac;
  double kecurrent = tdof_spin * boltz * t_current_spin;

  etas_mass[0] = tdof_spin * boltz * t_target / (t_freq * t_freq);
  for (ich = 1; ich < mtchain; ich++) etas_mass[ich] = boltz * t_target / (t_freq * t_freq);

  if (etas_mass[0] > 0.0)
    etas_dotdot[0] = (kecurrent - ke_target_spin) / etas_mass[0];
  else
    etas_dotdot[0] = 0.0;

  const double ncfac = 1.0 / nc_tchain;
  for (int iloop = 0; iloop < nc_tchain; iloop++) {

    for (ich = mtchain - 1; ich > 0; ich--) {
      expfac = exp(-ncfac * dt8 * etas_dot[ich + 1]);
      etas_dot[ich] *= expfac;
      etas_dot[ich] += etas_dotdot[ich] * ncfac * dt4;
      etas_dot[ich] *= tdrag_factor;
      etas_dot[ich] *= expfac;
    }

    expfac = exp(-ncfac * dt8 * etas_dot[1]);
    etas_dot[0] *= expfac;
    etas_dot[0] += etas_dotdot[0] * ncfac * dt4;
    etas_dot[0] *= tdrag_factor;
    etas_dot[0] *= expfac;

    const double factor_etas = exp(-ncfac * dthalf * etas_dot[0]);
    nh_vs_scale(factor_etas);

    // rescale the temperature to account for the velocity scaling

    t_current_spin *= factor_etas * factor_etas;
    kecurrent = tdof_spin * boltz * t_current_spin;

    if (etas_mass[0] > 0.0)
      etas_dotdot[0] = (kecurrent - ke_target_spin) / etas_mass[0];
    else
      etas_dotdot[0] = 0.0;

    for (ich = 0; ich < mtchain; ich++) etas[ich] += ncfac * dthalf * etas_dot[ich];

    etas_dot[0] *= expfac;
    etas_dot[0] += etas_dotdot[0] * ncfac * dt4;
    etas_dot[0] *= expfac;

    for (ich = 1; ich < mtchain; ich++) {
      expfac = exp(-ncfac * dt8 * etas_dot[ich + 1]);
      etas_dot[ich] *= expfac;
      etas_dotdot[ich] =
          (etas_mass[ich - 1] * etas_dot[ich - 1] * etas_dot[ich - 1] - boltz * t_target) /
          etas_mass[ich];
      etas_dot[ich] += etas_dotdot[ich] * ncfac * dt4;
      etas_dot[ich] *= expfac;
    }
  }
}

/* ----------------------------------------------------------------------
   the parent class calls this to propagate the lattice chain
   also propagate the spin chain, using the same target temperature
   with a frozen lattice the lattice chain is skipped entirely, so that it
   does not drift and pollute the conserved quantity
------------------------------------------------------------------------- */

void FixNHTSpin::nhc_temp_integrate()
{
  if (lattice_flag) FixNH::nhc_temp_integrate();

  ke_target_spin = tdof_spin * boltz * t_target;
  t_current_spin = compute_spin_temp();
  nhc_spin_integrate();
}

/* ----------------------------------------------------------------------
   half-step update of the spin velocities
------------------------------------------------------------------------- */

void FixNHTSpin::spin_kick()
{
  double **sp = atom->sp;
  double **s_dot = atom->v_s;
  double **fm = atom->fm;
  double **f_spin = atom->f_spin;
  double *s_mass = atom->s_mass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

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

/* ----------------------------------------------------------------------
   half-step velocity update, for the lattice and for the spins
------------------------------------------------------------------------- */

void FixNHTSpin::nve_v()
{
  if (lattice_flag) FixNH::nve_v();
  spin_kick();
}

/* ----------------------------------------------------------------------
   full-step position update, for the lattice and for the spins
------------------------------------------------------------------------- */

void FixNHTSpin::nve_x()
{
  if (lattice_flag) FixNH::nve_x();

  double **sp = atom->sp;
  double **s_dot = atom->v_s;
  double *s_mass = atom->s_mass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if ((sp[i][3] <= 0.0) || (s_mass[i] <= 0.0)) continue;

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

/* ----------------------------------------------------------------------
   the lattice thermostat must not rescale the velocities of a frozen
   lattice; the spin velocities are rescaled by the spin chain instead
------------------------------------------------------------------------- */

void FixNHTSpin::nh_v_temp()
{
  if (lattice_flag) FixNH::nh_v_temp();
}

/* ----------------------------------------------------------------------
   add the energy of the spin chain to the conserved quantity
------------------------------------------------------------------------- */

double FixNHTSpin::compute_scalar()
{
  double energy = lattice_flag ? FixNH::compute_scalar() : 0.0;

  if (!tstat_flag) return energy;

  energy += tdof_spin * boltz * t_target * etas[0];
  for (int ich = 1; ich < mtchain; ich++) energy += boltz * t_target * etas[ich];
  for (int ich = 0; ich < mtchain; ich++)
    energy += 0.5 * etas_mass[ich] * etas_dot[ich] * etas_dot[ich] * force->mvv2e;

  return energy;
}

/* ---------------------------------------------------------------------- */

int FixNHTSpin::size_restart_global()
{
  return FixNH::size_restart_global() + 2 * mtchain;
}

/* ---------------------------------------------------------------------- */

void FixNHTSpin::write_restart(FILE *fp)
{
  const int nsize = size_restart_global();
  auto list = new double[nsize];

  // the parent class writes its own state first, in its own layout

  int n = pack_restart_data(list);
  for (int ich = 0; ich < mtchain; ich++) list[n++] = etas[ich];
  for (int ich = 0; ich < mtchain; ich++) list[n++] = etas_dot[ich];

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(list, sizeof(double), n, fp);
  }
  delete[] list;
}

/* ---------------------------------------------------------------------- */

void FixNHTSpin::restart(char *buf)
{
  FixNH::restart(buf);

  // the spin chain state follows whatever the parent class consumed

  auto list = (double *) buf;
  int n = FixNH::size_restart_global();
  for (int ich = 0; ich < mtchain; ich++) etas[ich] = list[n++];
  for (int ich = 0; ich < mtchain; ich++) etas_dot[ich] = list[n++];
}
