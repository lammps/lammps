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

#include "fix_pimd_uvt.h"

#include "atom.h"
#include "compute.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "input.h"
#include "modify.h"
#include "universe.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

enum { TPRPMD = 1 };
enum { UVT = 1 };

/* ---------------------------------------------------------------------- */

FixPIMDUVT::FixPIMDUVT(LAMMPS *lmp, int narg, char **arg) :
    FixPIMDNVT(lmp, narg, arg, true), ustat_flag(1), mu(-3.5), Ne(nullptr),
    Ne_dot(nullptr), Ne_mass(nullptr), u_start(0.0), u_stop(0.0), u_current(0.0),
    u_target(0.0), u_freq(0.0), u_period(0.0), dedn_name(nullptr), dedn_which(ArgInfo::NONE),
    dedn_index(0), dedn_var(-1), dedn_compute(nullptr), dedn_fix(nullptr), dedn_current(0.0)
{
  parse_nvt_arguments(narg, arg, [this](int parse_narg, char **parse_arg, int &i) {
    return parse_uvt_keyword(parse_narg, parse_arg, i);
  });
  finish_nuclear_constructor_setup();
  finish_uvt_constructor_setup();
}

/* ---------------------------------------------------------------------- */

FixPIMDUVT::~FixPIMDUVT()
{
  delete[] Ne;
  delete[] Ne_dot;
  delete[] Ne_mass;
  delete[] dedn_name;
}

/* ---------------------------------------------------------------------- */

bool FixPIMDUVT::parse_uvt_keyword(int narg, char **arg, int &i)
{
  if (strcmp(arg[i], "mu") == 0) {
    if (i + 4 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} mu", style), error);
    u_start = utils::numeric(FLERR, arg[i + 1], false, lmp);
    u_stop = utils::numeric(FLERR, arg[i + 2], false, lmp);
    u_target = u_start;
    mu = u_target;
    u_period = utils::numeric(FLERR, arg[i + 3], false, lmp);
    i += 4;
    return true;
  }
  if (strcmp(arg[i], "ne") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} ne", style), error);
    if (!Ne) Ne = new double[1];
    *Ne = utils::numeric(FLERR, arg[i + 1], false, lmp);
    i += 2;
    return true;
  }
  if (strcmp(arg[i], "ne_velocity") == 0) {
    if (i + 2 > narg)
      utils::missing_cmd_args(FLERR, fmt::format("fix {} ne_velocity", style), error);
    if (!Ne_dot) Ne_dot = new double[1];
    *Ne_dot = utils::numeric(FLERR, arg[i + 1], false, lmp);
    i += 2;
    return true;
  }
  if (strcmp(arg[i], "dedn") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} dedn", style), error);
    parse_dedn_source(arg[i + 1]);
    i += 2;
    return true;
  }
  return false;
}

/* ---------------------------------------------------------------------- */

void FixPIMDUVT::finish_uvt_constructor_setup()
{
  if (ensemble != UVT) error->all(FLERR, "Fix {} only supports ensemble uvt", style);
  if (!Ne) error->all(FLERR, "Missing ne keyword for fix {}", style);
  if (!dedn_name) error->all(FLERR, "Missing dedn keyword for fix {}", style);
  if (u_period <= 0.0)
    error->all(FLERR, "Chemical-potential damping for fix {} must be > 0.0", style);

  if (!Ne_dot) {
    Ne_dot = new double[1];
    *Ne_dot = 0.0;
  }
  if (!Ne_mass) {
    Ne_mass = new double[1];
    *Ne_mass = 0.0;
  }

  const int old_size = size_vector;
  size_vector += 6;
  delete[] extlist;
  extlist = new int[size_vector];
  for (int i = 0; i < old_size; i++) extlist[i] = 1;
  for (int i = old_size; i < size_vector; i++) extlist[i] = 0;

  u_freq = 1.0 / u_period;
  u_current = u_start;
}

/* ---------------------------------------------------------------------- */

void FixPIMDUVT::setup_subclass_state()
{
  FixPIMDNVT::setup_subclass_state();
  resolve_dedn_source();
  u_freq = 1.0 / u_period;
  *Ne_mass = tdof * force->boltz * temp / (u_freq * u_freq);
}

/* ---------------------------------------------------------------------- */

void FixPIMDUVT::after_force_transform_hook()
{
  if (ustat_flag) refresh_dedn_cache();
}

/* ---------------------------------------------------------------------- */

void FixPIMDUVT::thermostat_step()
{
  if (!tstat_flag) return;

  compute_mu_target();
  nhc_mu_integrate();
  if (removecomflag) remove_com_motion();
}

/* ---------------------------------------------------------------------- */

void FixPIMDUVT::force_half_step()
{
  b_step();

  if (ustat_flag) {
    double dtfm = dthalf / *Ne_mass;
    if (universe->iworld == 0) *Ne_dot += dtfm * (-dedn_current + u_target);
    MPI_Bcast(Ne_dot, 1, MPI_DOUBLE, 0, universe->uworld);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDUVT::centroid_position_half_step()
{
  qc_step();

  if (ustat_flag) {
    if (universe->iworld == 0) *Ne += dtv * (*Ne_dot);
    MPI_Bcast(Ne, 1, MPI_DOUBLE, 0, universe->uworld);
  }
}

/* ---------------------------------------------------------------------- */

bool FixPIMDUVT::thermostat_chain_active() const
{
  return !nuclear_thermostat_off() || ne_thermostat_participates();
}

/* ---------------------------------------------------------------------- */

bool FixPIMDUVT::ne_thermostat_participates() const
{
  if (!ustat_flag) return false;
  if (method != TPRPMD) return true;
  return (np == 1) ? true : (universe->iworld != 0);
}

/* ---------------------------------------------------------------------- */

double FixPIMDUVT::ne_kinetic_current_share() const
{
  if (!ne_thermostat_participates()) return 0.0;
  return 0.5 * inverse_np * (*Ne_mass) * (*Ne_dot) * (*Ne_dot);
}

/* ---------------------------------------------------------------------- */

int FixPIMDUVT::subclass_restart_size() const
{
  return 1 + (ustat_flag ? 2 : 0);
}

/* ---------------------------------------------------------------------- */

int FixPIMDUVT::pack_subclass_restart(double *list, int n) const
{
  list[n++] = ustat_flag;
  if (ustat_flag) {
    list[n++] = *Ne;
    list[n++] = *Ne_dot;
  }
  return n;
}

/* ---------------------------------------------------------------------- */

int FixPIMDUVT::unpack_subclass_restart(const double *list, int n)
{
  int flag = static_cast<int>(list[n++]);
  if (flag && ustat_flag) {
    *Ne = list[n++];
    *Ne_dot = list[n++];
  }
  return n;
}

/* ---------------------------------------------------------------------- */

int FixPIMDUVT::subclass_vector_size() const
{
  return 6;
}

/* ---------------------------------------------------------------------- */

double FixPIMDUVT::compute_subclass_vector(int n) const
{
  const double kt_local = force->boltz * temp;
  if (n == 0) return *Ne;
  if (n == 1) return *Ne_dot;
  if (n == 2) return dedn_current;
  if (n == 3) return u_target;
  if (n == 4) return 0.5 * (*Ne_mass) * (*Ne_dot) * (*Ne_dot);
  if (n == 5) return kt_local * eta[0] - u_target * (*Ne);
  return 0.0;
}

/* ---------------------------------------------------------------------- */

void FixPIMDUVT::compute_mu_target()
{
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  u_target = u_start + delta * (u_stop - u_start);
  mu = u_target;
}

/* ---------------------------------------------------------------------- */

void FixPIMDUVT::nhc_mu_integrate()
{
  double kecurrent = compute_nuclear_kinetic_energy();
  double t_current = kecurrent / force->boltz / tdof;

  if (thermostat_chain_active()) update_chain0_acceleration(ne_kinetic_current_share());

  double ncfac = 1.0 / nc_tchain;
  const double chain_target = chain_target_energy();
  for (int iloop = 0; iloop < nc_tchain; iloop++) {
    bool active = thermostat_chain_active();
    double expfac = 1.0;
    if (active) {
      propagate_chain_tail_halfstep(ncfac);
      expfac = propagate_chain0_halfstep(ncfac, true);
    }

    double eta_dot_k = eta_dot[0];
    double eta_dot_ave = 0.0;
    MPI_Allreduce(&eta_dot_k, &eta_dot_ave, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
    if (method == TPRPMD && np != 1)
      eta_dot_ave /= (static_cast<double>(np) - 1);
    else
      eta_dot_ave *= inverse_np;

    *Ne_dot *= exp(-ncfac * dthalf * eta_dot_ave);

    if (active) {
      update_scaled_nuclear_kinetic(t_current, kecurrent);
      if (eta_mass[0] > 0.0)
        eta_dotdot[0] = (kecurrent + ne_kinetic_current_share() - chain0_target_energy()) /
            eta_mass[0];
      else
        eta_dotdot[0] = 0.0;

      advance_chain_positions(ncfac);
      complete_chain0_halfstep(ncfac, expfac);
      update_outer_chain_accelerations(chain_target);
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixPIMDUVT::evaluate_dedn()
{
  if (dedn_which == ArgInfo::VARIABLE) {
    if (dedn_index == 0) return input->variable->compute_equal(dedn_var);
    double *varvec;
    int nvec = input->variable->compute_vector(dedn_var, &varvec);
    if (dedn_index > nvec)
      error->all(FLERR, "Variable {} for fix {} vector is accessed out-of-range{}", dedn_name,
                 style, utils::errorurl(20));
    return varvec[dedn_index - 1];
  }

  if (dedn_which == ArgInfo::COMPUTE) {
    if (dedn_index == 0) {
      dedn_compute->compute_scalar();
      dedn_compute->invoked_scalar = update->ntimestep;
      dedn_compute->invoked_flag |= Compute::INVOKED_SCALAR;
      return dedn_compute->scalar;
    }
    dedn_compute->compute_vector();
    dedn_compute->invoked_vector = update->ntimestep;
    dedn_compute->invoked_flag |= Compute::INVOKED_VECTOR;
    return dedn_compute->vector[dedn_index - 1];
  }

  if (dedn_index == 0) return dedn_fix->compute_scalar();
  return dedn_fix->compute_vector(dedn_index - 1);
}

/* ---------------------------------------------------------------------- */

void FixPIMDUVT::refresh_dedn_cache()
{
  double dedn_local = evaluate_dedn();
  double dedn_avg = 0.0;
  MPI_Allreduce(&dedn_local, &dedn_avg, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  dedn_current = dedn_avg * inverse_np;
  u_current = dedn_current;
}

/* ---------------------------------------------------------------------- */

void FixPIMDUVT::parse_dedn_source(const char *arg)
{
  ArgInfo argi(arg);
  if ((argi.get_type() == ArgInfo::UNKNOWN) || (argi.get_type() == ArgInfo::NONE) ||
      (argi.get_dim() > 1))
    error->all(FLERR, "Illegal dedn reference {} for fix {}", arg, style);

  delete[] dedn_name;
  dedn_name = argi.copy_name();
  dedn_which = argi.get_type();
  dedn_index = argi.get_index1();
}

/* ---------------------------------------------------------------------- */

void FixPIMDUVT::resolve_dedn_source()
{
  dedn_var = -1;
  dedn_compute = nullptr;
  dedn_fix = nullptr;

  if (dedn_index < 0)
    error->all(FLERR, "Illegal dedn index {} for fix {}", dedn_index, style);

  if (dedn_which == ArgInfo::VARIABLE) {
    dedn_var = input->variable->find(dedn_name);
    if (dedn_var < 0) error->all(FLERR, "Variable {} for fix {} does not exist", dedn_name, style);
    if (dedn_index == 0) {
      if (!input->variable->equalstyle(dedn_var))
        error->all(FLERR, "Variable {} for fix {} is not equal-style", dedn_name, style);
    } else if (!input->variable->vectorstyle(dedn_var)) {
      error->all(FLERR, "Variable {} for fix {} is not vector-style", dedn_name, style);
    }
  } else if (dedn_which == ArgInfo::COMPUTE) {
    dedn_compute = modify->get_compute_by_id(dedn_name);
    if (!dedn_compute) error->all(FLERR, "Compute {} for fix {} does not exist", dedn_name, style);
    if (dedn_index == 0) {
      if (!dedn_compute->scalar_flag)
        error->all(FLERR, "Compute {} for fix {} does not compute a scalar", dedn_name, style);
    } else {
      if (!dedn_compute->vector_flag)
        error->all(FLERR, "Compute {} for fix {} does not compute a vector", dedn_name, style);
      if (dedn_compute->size_vector > 0 && dedn_index > dedn_compute->size_vector)
        error->all(FLERR, "Compute {} for fix {} vector index {} is out of range", dedn_name,
                   style, dedn_index);
    }
  } else if (dedn_which == ArgInfo::FIX) {
    dedn_fix = modify->get_fix_by_id(dedn_name);
    if (!dedn_fix) error->all(FLERR, "Fix {} for fix {} does not exist", dedn_name, style);
    if (dedn_index == 0) {
      if (!dedn_fix->scalar_flag)
        error->all(FLERR, "Fix {} for fix {} does not compute a scalar", dedn_name, style);
    } else {
      if (!dedn_fix->vector_flag)
        error->all(FLERR, "Fix {} for fix {} does not compute a vector", dedn_name, style);
      if (dedn_fix->size_vector > 0 && dedn_index > dedn_fix->size_vector)
        error->all(FLERR, "Fix {} for fix {} vector index {} is out of range", dedn_name, style,
                   dedn_index);
    }
  } else {
    error->all(FLERR, "dedn for fix {} must be a variable, compute, or fix reference", style);
  }
}
