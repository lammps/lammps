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

#include "fix_pimd_nvt.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathConst::THIRD;

enum { PHYSICAL, NORMAL };
enum { BAOAB, OBABO };
enum { SINGLE_PROC, MULTI_PROC };

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::init_nvt_defaults()
{
  pilescale = 1.0;
  tstat_flag = 1;

  mtchain = 3;
  nc_tchain = 1;
  t_period = 0.0;
  drag = 0.0;

  factor_eta = 1.0;
  tdrag_factor = 1.0;
  t_freq = 0.0;
  tdof = 0.0;
  tdof_override_flag = 0;
  tdof_override = 0.0;
  ke_target = 0.0;
  ecouple_work = 0.0;
  dthalf = dt4 = dt8 = 0.0;

  fixedpoint[0] = fixedpoint[1] = fixedpoint[2] = 0.0;

  scalar_flag = 1;
  extscalar = 1;
  ecouple_flag = 1;
  thermo_modify_colname = 1;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::parse_nvt_arguments(int narg, char **arg, const KeywordParser &subclass_parser)
{
  for (int i = 3; i < narg;) {
    if (subclass_parser && subclass_parser(narg, arg, i)) continue;
    if (parse_nvt_keyword(narg, arg, i)) continue;
    if (FixPIMDNVE::parse_common_keyword(narg, arg, i)) continue;
    error->all(FLERR, "Unknown keyword {} for fix {}", arg[i], style);
  }

  if (t_period <= 0.0) error->all(FLERR, "Temperature damping for fix {} must be > 0.0", style);
}

/* ---------------------------------------------------------------------- */

bool FixPIMDNVT::parse_nvt_keyword(int narg, char **arg, int &i)
{
  if (strcmp(arg[i], "method") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} method", style), error);
    if (strcmp(arg[i + 1], "nmpimd") != 0) error->all(FLERR, "Fix {} only supports method nmpimd", style);
    i += 2;
    return true;
  }
  if (strcmp(arg[i], "ensemble") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} ensemble", style), error);
    if (strcmp(arg[i + 1], "nvt") == 0)
      error->all(FLERR, "Fix {} is already NVT; remove the ensemble keyword", style);
    if (strcmp(arg[i + 1], "uvt") == 0)
      error->all(FLERR, "Fix {} does not support ensemble uvt; use fix pimd/uvt instead", style);
    error->all(FLERR, "Fix {} only supports the NVT ensemble", style);
  }
  if (strcmp(arg[i], "thermostat") == 0) {
    if (i + 2 > narg)
      utils::missing_cmd_args(FLERR, fmt::format("fix {} thermostat", style), error);
    if (strcmp(arg[i + 1], "NHC") != 0)
      error->all(FLERR, "Fix {} only supports thermostat NHC", style);
    i += 2;
    return true;
  }
  if (strcmp(arg[i], "Tdamp") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} Tdamp", style), error);
    t_period = utils::numeric(FLERR, arg[i + 1], false, lmp);
    i += 2;
    return true;
  }
  if (strcmp(arg[i], "tchain") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} tchain", style), error);
    mtchain = utils::inumeric(FLERR, arg[i + 1], false, lmp);
    i += 2;
    return true;
  }
  if (strcmp(arg[i], "tloop") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} tloop", style), error);
    nc_tchain = utils::inumeric(FLERR, arg[i + 1], false, lmp);
    i += 2;
    return true;
  }
  if (strcmp(arg[i], "drag") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} drag", style), error);
    drag = utils::numeric(FLERR, arg[i + 1], false, lmp);
    i += 2;
    return true;
  }
  if (strcmp(arg[i], "tdof") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} tdof", style), error);
    tdof_override = utils::numeric(FLERR, arg[i + 1], false, lmp);
    if (tdof_override <= 0.0) error->all(FLERR, "Temperature DOF override for fix {} must be > 0.0", style);
    tdof_override_flag = 1;
    i += 2;
    return true;
  }
  if ((strcmp(arg[i], "barostat") == 0) || (strcmp(arg[i], "iso") == 0) ||
      (strcmp(arg[i], "aniso") == 0) || (strcmp(arg[i], "taup") == 0) ||
      (strcmp(arg[i], "fixedpoint") == 0)) {
    error->all(FLERR, "Pressure control is not supported by fix {}", style);
  }
  if ((strcmp(arg[i], "seed") == 0) || (strcmp(arg[i], "PILE_L_temp") == 0)) {
    error->all(FLERR, "Legacy thermostat options are not supported by fix {}", style);
  }
  return false;
}

/* ---------------------------------------------------------------------- */

FixPIMDNVT::FixPIMDNVT(LAMMPS *lmp, int narg, char **arg, bool) :
    FixPIMDNVE(lmp, narg, arg, true), eta(nullptr), eta_dot(nullptr), eta_dotdot(nullptr),
    eta_mass(nullptr), tau_k(nullptr)
{
  init_nvt_defaults();

  if (narg < 4) utils::missing_cmd_args(FLERR, std::string("fix ") + style, error);
}

/* ---------------------------------------------------------------------- */

FixPIMDNVT::FixPIMDNVT(LAMMPS *lmp, int narg, char **arg) :
    FixPIMDNVT(lmp, narg, arg, true)
{
  parse_nvt_arguments(narg, arg, {});
  finish_nuclear_constructor_setup();
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::finish_nuclear_constructor_setup()
{
  if (tstat_flag) {
    eta = new double[mtchain];
    eta_dot = new double[mtchain + 1];
    eta_dot[mtchain] = 0.0;
    eta_dotdot = new double[mtchain];
    for (int ich = 0; ich < mtchain; ich++) eta[ich] = eta_dot[ich] = eta_dotdot[ich] = 0.0;
    eta_mass = new double[mtchain];
    size_vector += 4 * mtchain;
  }

  finish_constructor_setup();

  fixedpoint[0] = 0.5 * (domain->boxlo[0] + domain->boxhi[0]);
  fixedpoint[1] = 0.5 * (domain->boxlo[1] + domain->boxhi[1]);
  fixedpoint[2] = 0.5 * (domain->boxlo[2] + domain->boxhi[2]);
}

/* ---------------------------------------------------------------------- */

FixPIMDNVT::~FixPIMDNVT()
{
  delete[] tau_k;
  delete[] eta;
  delete[] eta_dot;
  delete[] eta_dotdot;
  delete[] eta_mass;
}

/* ---------------------------------------------------------------------- */

double FixPIMDNVT::chain0_target_energy() const
{
  return static_cast<double>(np) * ke_target;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::setup_subclass_state()
{
  dthalf = 0.5 * update->dt;
  dt4 = 0.25 * update->dt;
  dt8 = 0.125 * update->dt;
  nhc_init();
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::initial_integrate(int /*vflag*/)
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  imageint *image = atom->image;
  if (mapflag) {
    for (int i = 0; i < nlocal; i++) domain->unmap(x[i], image[i]);
  }
  if (integrator == OBABO) {
    thermostat_step();
    force_half_step();
    inter_replica_comm(x);
    if (cmode == SINGLE_PROC)
      nmpimd_transform(bufsortedall, x, M_x2xp[universe->iworld]);
    else if (cmode == MULTI_PROC)
      nmpimd_transform(bufbeads, x, M_x2xp[universe->iworld]);
    centroid_position_half_step();
    a_step();
    centroid_position_half_step();
    a_step();
  } else if (integrator == BAOAB) {
    force_half_step();
    inter_replica_comm(x);
    if (cmode == SINGLE_PROC)
      nmpimd_transform(bufsortedall, x, M_x2xp[universe->iworld]);
    else if (cmode == MULTI_PROC)
      nmpimd_transform(bufbeads, x, M_x2xp[universe->iworld]);
    centroid_position_half_step();
    a_step();
    thermostat_step();
    centroid_position_half_step();
    a_step();
  } else {
    error->universe_all(FLERR, fmt::format("Unknown integrator parameter for fix {}. Only obabo "
                                           "and baoab integrators are supported!",
                                           style));
  }
  collect_xc();

  compute_spring_energy();
  compute_t_prim();
  compute_p_prim();
  inter_replica_comm(x);
  if (cmode == SINGLE_PROC)
    nmpimd_transform(bufsortedall, x, M_xp2x[universe->iworld]);
  else if (cmode == MULTI_PROC)
    nmpimd_transform(bufbeads, x, M_xp2x[universe->iworld]);

  if (mapflag) {
    for (int i = 0; i < nlocal; i++) domain->unmap_inv(x[i], image[i]);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::final_integrate()
{
  force_half_step();

  if (integrator == OBABO) {
    thermostat_step();
  } else if (integrator == BAOAB) {

  } else {
    error->universe_all(FLERR, fmt::format("Unknown integrator parameter for fix {}", style));
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::nhc_init()
{
  if (kt <= 0.0 || hbar <= 0.0)
    error->universe_all(FLERR, fmt::format("Fix {} requires positive kt and hbar in nhc_init",
                                           style));

  const double beta_local = 1.0 / kt;
  const double omega_np_local = np / beta_local / hbar;
  const double omega_np_dt_half = omega_np_local * update->dt * 0.5;

  if (fmmode == PHYSICAL) {
    for (int i = 0; i < np; i++) {
      _omega_k[i] = omega_np_local * sqrt(lam[i]) / sqrt(fmass);
      Lan_c[i] = cos(sqrt(lam[i]) * omega_np_dt_half);
      Lan_s[i] = sin(sqrt(lam[i]) * omega_np_dt_half);
    }
  } else {
    for (int i = 0; i < np; i++) {
      _omega_k[i] = omega_np_local / sqrt(fmass);
      Lan_c[i] = cos(omega_np_dt_half);
      Lan_s[i] = sin(omega_np_dt_half);
    }
  }

  if (tstat_flag) {
    t_freq = 1.0 / t_period;
    tdrag_factor = 1.0 - (update->dt * t_freq * drag / nc_tchain);
  }

  int fix_dof = 0;
  for (auto &ifix : modify->get_fix_list())
    if (ifix->dof_flag) fix_dof += ifix->dof(igroup);
  int extra_dof = (removecomflag && universe->iworld == 0) ? domain->dimension : 0;
  tdof = domain->dimension * group->count(igroup);
  tdof -= extra_dof + fix_dof;
  if (tdof_override_flag) tdof = tdof_override;
  if (tdof <= 0.0) error->all(FLERR, "Temperature DOF for fix {} must be > 0.0", style);

  ke_target = tdof * force->boltz * temp;

  delete[] tau_k;
  tau_k = new double[np];
  tau_k[0] = t_period;
  for (int i = 1; i < np; i++) tau_k[i] = 0.5 / pilescale / _omega_k[i];

  if (tstat_flag && nc_tchain == 1 && np > 1) {
    double tau_min = tau_k[1];
    for (int i = 2; i < np; i++) tau_min = MIN(tau_min, tau_k[i]);
    if (tau_min > 0.0) {
      const int required_tloop = MAX(1, static_cast<int>(ceil(update->dt / tau_min)));
      if (required_tloop > nc_tchain) {
        nc_tchain = required_tloop;
        if (universe->me == 0) {
          utils::logmesg(
              lmp,
              fmt::format("  Auto-increased NHC tloop to {:d} so dt/tau_min = {:.6f} does not "
                          "under-resolve the fastest internal thermostat mode.\n",
                          nc_tchain, update->dt / tau_min));
        }
      }
    }
  }

  if (tstat_flag) {
    const double omega_np_local_unscaled = np / beta_local / hbar;
    const double chain0_target = chain0_target_energy();
    const double chain_target = chain_target_energy();

    // Preserve the old NHC parameterization from the monolithic pimd/langevin
    // implementation: the first chain mass is controlled by Tdamp, while
    // higher chain masses use Tdamp for classical MD and omega_np for NMPIMD.
    eta_mass[0] = chain0_target / (t_freq * t_freq);
    for (int ich = 1; ich < mtchain; ich++) {
      if (np == 1)
        eta_mass[ich] = chain_target / (t_freq * t_freq);
      else
        eta_mass[ich] = chain_target / (omega_np_local_unscaled * omega_np_local_unscaled);
      if (eta_mass[ich] > 0.0)
        eta_dotdot[ich] =
            (eta_mass[ich - 1] * eta_dot[ich - 1] * eta_dot[ich - 1] - chain_target) /
            eta_mass[ich];
      else
        eta_dotdot[ich] = 0.0;
    }
    if (!thermostat_chain_active()) {
      for (int ich = 1; ich < mtchain; ich++) eta_dotdot[ich] = 0.0;
    }
  }

  std::string out = "Initializing path-integral Nose-Hoover thermostat chain...\n";
  out += "  Bead ID    |    omega    |    timescale\n";
  for (int i = 0; i < np; i++) {
    out += fmt::format("      {:d}     {:.8e} {:.8e}\n", i, _omega_k[i], tau_k[i]);
  }
  out += "  NHC thermostat successfully initialized!\n\n";
  if (universe->me == 0) utils::logmesg(lmp, out);
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::o_step()
{
  if (tstat_flag) nhc_temp_integrate();
}

/* ---------------------------------------------------------------------- */

double FixPIMDNVT::compute_nuclear_kinetic_energy() const
{
  int *mask = atom->mask;
  int *type = atom->type;
  double **v = atom->v;
  int nlocal = atom->nlocal;
  double kecurrent = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      kecurrent += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * mass[type[i]];
  }
  double ketotal = 0.0;
  kecurrent *= force->mvv2e;
  MPI_Allreduce(&kecurrent, &ketotal, 1, MPI_DOUBLE, MPI_SUM, world);
  return ketotal;
}

/* ---------------------------------------------------------------------- */

bool FixPIMDNVT::thermostat_chain_active() const
{
  return true;
}

/* ---------------------------------------------------------------------- */

double FixPIMDNVT::chain_target_energy() const
{
  return thermostat_chain_active() ? static_cast<double>(np) * force->boltz * temp : 0.0;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::update_chain0_acceleration(double extra_ke)
{
  if (!thermostat_chain_active()) return;

  const double chain0_target = chain0_target_energy();
  if (eta_mass[0] > 0.0)
    eta_dotdot[0] = (compute_nuclear_kinetic_energy() + extra_ke - chain0_target) / eta_mass[0];
  else
    eta_dotdot[0] = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::propagate_chain_tail_halfstep(double ncfac)
{
  for (int ich = mtchain - 1; ich > 0; ich--) {
    double expfac = exp(-ncfac * dt8 * eta_dot[ich + 1]);
    eta_dot[ich] *= expfac;
    eta_dot[ich] += eta_dotdot[ich] * ncfac * dt4;
    eta_dot[ich] *= tdrag_factor;
    eta_dot[ich] *= expfac;
  }
}

/* ---------------------------------------------------------------------- */

double FixPIMDNVT::propagate_chain0_halfstep(double ncfac, bool apply_velocity_scaling)
{
  double expfac = exp(-ncfac * dt8 * eta_dot[1]);
  eta_dot[0] *= expfac;
  eta_dot[0] += eta_dotdot[0] * ncfac * dt4;
  eta_dot[0] *= tdrag_factor;
  eta_dot[0] *= expfac;

  factor_eta = exp(-ncfac * dthalf * eta_dot[0]);
  if (apply_velocity_scaling) nh_v_temp();
  return expfac;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::update_scaled_nuclear_kinetic(double &t_current, double &kecurrent) const
{
  t_current *= factor_eta * factor_eta;
  kecurrent = tdof * force->boltz * t_current;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::advance_chain_positions(double ncfac)
{
  for (int ich = 0; ich < mtchain; ich++) eta[ich] += ncfac * dthalf * eta_dot[ich];
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::complete_chain0_halfstep(double ncfac, double expfac)
{
  eta_dot[0] *= expfac;
  eta_dot[0] += eta_dotdot[0] * ncfac * dt4;
  eta_dot[0] *= expfac;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::update_outer_chain_accelerations(double chain_target)
{
  for (int ich = 1; ich < mtchain; ich++) {
    if (eta_mass[ich] > 0.0)
      eta_dotdot[ich] =
          (eta_mass[ich - 1] * eta_dot[ich - 1] * eta_dot[ich - 1] - chain_target) /
          eta_mass[ich];
    else
      eta_dotdot[ich] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::complete_chain_tail_halfstep(double ncfac, double chain_target)
{
  for (int ich = 1; ich < mtchain; ich++) {
    double expfac = exp(-ncfac * dt8 * eta_dot[ich + 1]);
    eta_dot[ich] *= expfac;
    if (eta_mass[ich] > 0.0)
      eta_dotdot[ich] =
          (eta_mass[ich - 1] * eta_dot[ich - 1] * eta_dot[ich - 1] - chain_target) /
          eta_mass[ich];
    else
      eta_dotdot[ich] = 0.0;
    eta_dot[ich] += eta_dotdot[ich] * ncfac * dt4;
    eta_dot[ich] *= expfac;
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::nhc_temp_integrate()
{
  double kecurrent = compute_nuclear_kinetic_energy();
  double t_current = kecurrent / force->boltz / tdof;
  if (!thermostat_chain_active()) return;

  update_chain0_acceleration(0.0);

  double ncfac = 1.0 / nc_tchain;
  const double chain_target = chain_target_energy();
  for (int iloop = 0; iloop < nc_tchain; iloop++) {
    propagate_chain_tail_halfstep(ncfac);
    double expfac = propagate_chain0_halfstep(ncfac, true);

    update_scaled_nuclear_kinetic(t_current, kecurrent);
    if (eta_mass[0] > 0.0)
      eta_dotdot[0] = (kecurrent - chain0_target_energy()) / eta_mass[0];
    else
      eta_dotdot[0] = 0.0;

    advance_chain_positions(ncfac);
    complete_chain0_halfstep(ncfac, expfac);
    complete_chain_tail_halfstep(ncfac, chain_target);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::thermostat_step()
{
  if (tstat_flag) {
    o_step();
    if (removecomflag) remove_com_motion();
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::force_half_step()
{
  b_step();
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::centroid_position_half_step()
{
  qc_step();
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::nh_v_temp()
{
  const double work_delta = thermostat_work_delta(factor_eta);

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      v[i][0] *= factor_eta;
      v[i][1] *= factor_eta;
      v[i][2] *= factor_eta;
    }
  }

  ecouple_work += work_delta;
}

/* ---------------------------------------------------------------------- */

double FixPIMDNVT::thermostat_work_delta(double scale_factor) const
{
  const double kinetic_before_local = local_kinetic_energy_sum(true);
  const double kinetic_after_local = kinetic_before_local * scale_factor * scale_factor;
  double work_delta_local = kinetic_before_local - kinetic_after_local;
  double work_delta = 0.0;

  MPI_Allreduce(&work_delta_local, &work_delta, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  return work_delta;
}

/* ---------------------------------------------------------------------- */

double FixPIMDNVT::compute_scalar()
{
  return tstat_flag ? ecouple_work : 0.0;
}

/* ---------------------------------------------------------------------- */

std::string FixPIMDNVT::get_thermo_colname(int n)
{
  if (n == -1) return fmt::format("f_{}:ecouple", id);
  return Fix::get_thermo_colname(n);
}

/* ---------------------------------------------------------------------- */

int FixPIMDNVT::base_restart_size() const
{
  int nsize = 2;
  if (tstat_flag) nsize += 1 + 2 * mtchain;
  return nsize;
}

/* ---------------------------------------------------------------------- */

int FixPIMDNVT::pack_base_restart(double *list) const
{
  int n = 0;
  list[n++] = tstat_flag;
  list[n++] = ecouple_work;
  if (tstat_flag) {
    list[n++] = mtchain;
    for (int ich = 0; ich < mtchain; ich++) list[n++] = eta[ich];
    for (int ich = 0; ich < mtchain; ich++) list[n++] = eta_dot[ich];
  }
  return n;
}

/* ---------------------------------------------------------------------- */

int FixPIMDNVT::unpack_base_restart(const double *list)
{
  int n = 0;
  int flag = static_cast<int>(list[n++]);
  ecouple_work = list[n++];
  if (flag) {
    int m = static_cast<int>(list[n++]);
    if (tstat_flag && m == mtchain) {
      for (int ich = 0; ich < mtchain; ich++) eta[ich] = list[n++];
      for (int ich = 0; ich < mtchain; ich++) eta_dot[ich] = list[n++];
    } else {
      n += 2 * m;
    }
  }
  return n;
}

/* ---------------------------------------------------------------------- */

int FixPIMDNVT::nuclear_vector_size() const
{
  int nsize = FixPIMDNVE::nuclear_vector_size();
  if (tstat_flag) nsize += 4 * mtchain;
  return nsize;
}

/* ---------------------------------------------------------------------- */

double FixPIMDNVT::compute_nuclear_vector(int n) const
{
  const int prefix = FixPIMDNVE::nuclear_vector_size();
  if (n < prefix) return FixPIMDNVE::compute_nuclear_vector(n);
  n -= prefix;

  int ilen;
  if (tstat_flag) {
    ilen = mtchain;
    if (n < ilen) return eta[n];
    n -= ilen;
    ilen = mtchain;
    if (n < ilen) return eta_dot[n];
    n -= ilen;
  }

  const double chain0_target = chain0_target_energy();
  const double chain_target = chain_target_energy();
  if (tstat_flag) {
    ilen = mtchain;
    if (n < ilen) {
      if (n == 0) return chain0_target * eta[0];
      return chain_target * eta[n];
    }
    n -= ilen;
    ilen = mtchain;
    if (n < ilen) return 0.5 * eta_mass[n] * eta_dot[n] * eta_dot[n];
  }

  return 0.0;
}
