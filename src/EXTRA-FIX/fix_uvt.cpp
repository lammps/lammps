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

#include "fix_uvt.h"

#include "arg_info.h"
#include "error.h"
#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "group.h"
#include "input.h"
#include "kspace.h"
#include "neighbor.h"
#include "modify.h"
#include "force.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>
#include <string>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;

namespace {
enum { ISO, ANISO, TRICLINIC };
enum { NOBIAS, BIAS };

struct UvtArgCache {
  std::vector<std::string> storage;
  std::vector<char *> argv;
  int argc = 0;
  int source_narg = 0;
  char **source_arg = nullptr;
};

thread_local UvtArgCache uvt_cache;
thread_local bool uvt_cache_ready = false;

int uvt_skip_count(const char *keyword)
{
  if (strcmp(keyword, "mu") == 0) return 4;
  if (strcmp(keyword, "ne") == 0) return 2;
  if (strcmp(keyword, "ne_velocity") == 0) return 2;
  if (strcmp(keyword, "dedn") == 0) return 2;
  if (strcmp(keyword, "dedn_defer") == 0) return 2;
  return 0;
}

void prepare_uvt_args(int narg, char **arg)
{
  if (uvt_cache_ready && uvt_cache.source_narg == narg && uvt_cache.source_arg == arg) return;

  uvt_cache.storage.clear();
  uvt_cache.argv.clear();
  uvt_cache.source_narg = narg;
  uvt_cache.source_arg = arg;

  for (int i = 0; i < narg;) {
    if (i >= 3) {
      const int skip = uvt_skip_count(arg[i]);
      if (skip > 0) {
        i += skip;
        continue;
      }
    }
    uvt_cache.storage.emplace_back(arg[i]);
    ++i;
  }

  uvt_cache.argv.reserve(uvt_cache.storage.size());
  for (auto &s : uvt_cache.storage) uvt_cache.argv.push_back(const_cast<char *>(s.c_str()));
  uvt_cache.argc = static_cast<int>(uvt_cache.argv.size());
  uvt_cache_ready = true;
}

int uvt_argc(int narg, char **arg)
{
  prepare_uvt_args(narg, arg);
  return uvt_cache.argc;
}

char **uvt_argv(int narg, char **arg)
{
  prepare_uvt_args(narg, arg);
  return uvt_cache.argv.data();
}

void release_uvt_args()
{
  uvt_cache_ready = false;
  uvt_cache.storage.clear();
  uvt_cache.argv.clear();
  uvt_cache.argc = 0;
  uvt_cache.source_narg = 0;
  uvt_cache.source_arg = nullptr;
}
}    // namespace

/* ---------------------------------------------------------------------- */

FixUVT::FixUVT(LAMMPS *lmp, int narg, char **arg) :
    FixNH(lmp, uvt_argc(narg, arg), uvt_argv(narg, arg)),
    u_start(0.0), u_stop(0.0), u_current(0.0), u_target(0.0), u_freq(0.0), ustat_flag(0),
    Ne(nullptr), Ne_dot(nullptr), Ne_mass(nullptr), dedn_name(nullptr), dedn_which(ArgInfo::NONE),
    dedn_index(0), dedn_var(-1), dedn_compute(nullptr), dedn_fix(nullptr), dedn_current(0.0),
    dedn_defer(0)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, std::string("fix ") + style, error);

  double u_period = 0.0;
  bool mu_seen = false;
  bool ne_seen = false;
  bool dedn_seen = false;

  for (int i = 3; i < narg;) {
    if (strcmp(arg[i], "mu") == 0) {
      if (i + 4 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} mu", style), error);
      ustat_flag = 1;
      mu_seen = true;
      u_start = utils::numeric(FLERR, arg[i + 1], false, lmp);
      u_target = u_start;
      u_stop = utils::numeric(FLERR, arg[i + 2], false, lmp);
      u_period = utils::numeric(FLERR, arg[i + 3], false, lmp);
      i += 4;
    } else if (strcmp(arg[i], "ne") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} ne", style), error);
      ne_seen = true;
      if (!Ne) Ne = new double[1];
      *Ne = utils::numeric(FLERR, arg[i + 1], false, lmp);
      i += 2;
    } else if (strcmp(arg[i], "ne_velocity") == 0) {
      if (i + 2 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} ne_velocity", style), error);
      if (!Ne_dot) Ne_dot = new double[1];
      *Ne_dot = utils::numeric(FLERR, arg[i + 1], false, lmp);
      i += 2;
    } else if (strcmp(arg[i], "dedn") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} dedn", style), error);
      dedn_seen = true;
      parse_dedn_source(arg[i + 1]);
      i += 2;
    } else if (strcmp(arg[i], "dedn_defer") == 0) {
      if (i + 2 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} dedn_defer", style), error);
      if (strcmp(arg[i + 1], "yes") == 0)
        dedn_defer = 1;
      else if (strcmp(arg[i + 1], "no") == 0)
        dedn_defer = 0;
      else
        error->all(FLERR, "Fix {} dedn_defer must be yes or no", style);
      i += 2;
    } else {
      ++i;
    }
  }

  if (!tstat_flag) error->all(FLERR, "Temperature control must be used with fix uvt");
  if (pstat_flag) error->all(FLERR, "Pressure control can not be used with fix uvt");
  if (!mu_seen) error->all(FLERR, "Missing mu keyword for fix {}", style);
  if (!ne_seen) error->all(FLERR, "Missing ne keyword for fix {}", style);
  if (!dedn_seen) error->all(FLERR, "Missing dedn keyword for fix {}", style);
  if (u_period <= 0.0) error->all(FLERR, "Chemical-potential damping for fix {} must be > 0.0", style);
  if (!Ne_dot) {
    Ne_dot = new double[1];
    *Ne_dot = 0.0;
  }

  u_freq = 1.0 / u_period;
  u_current = u_start;

  if (!Ne_mass) Ne_mass = new double[1];
  *Ne_mass = 0.0;
  size_vector += 6;

  id_temp = utils::strdup(std::string(id) + "_temp");
  modify->add_compute(fmt::format("{} {} temp", id_temp, group->names[igroup]));
  tcomputeflag = 1;

  release_uvt_args();
}

/* ---------------------------------------------------------------------- */

FixUVT::~FixUVT()
{
  if (copymode) return;

  delete[] Ne;
  delete[] Ne_dot;
  delete[] Ne_mass;
  delete[] dedn_name;
}

/* ---------------------------------------------------------------------- */

void FixUVT::init()
{
  FixNH::init();

  dedn_var = -1;
  dedn_compute = nullptr;
  dedn_fix = nullptr;

  if (dedn_which == ArgInfo::VARIABLE) {
    dedn_var = input->variable->find(dedn_name);
    if (dedn_var < 0)
      error->all(FLERR, "Variable {} for fix {} does not exist", dedn_name, style);
    if (dedn_index == 0) {
      if (!input->variable->equalstyle(dedn_var))
        error->all(FLERR, "Variable {} for fix {} is not equal-style", dedn_name, style);
    } else {
      if (!input->variable->vectorstyle(dedn_var))
        error->all(FLERR, "Variable {} for fix {} is not vector-style", dedn_name, style);
    }
  } else if (dedn_which == ArgInfo::COMPUTE) {
    dedn_compute = modify->get_compute_by_id(dedn_name);
    if (!dedn_compute)
      error->all(FLERR, "Compute {} for fix {} does not exist", dedn_name, style);
    if (dedn_index == 0 && !dedn_compute->scalar_flag)
      error->all(FLERR, "Compute {} for fix {} does not calculate a global scalar", dedn_name,
                 style);
    if (dedn_index > 0) {
      if (!dedn_compute->vector_flag)
        error->all(FLERR, "Compute {} for fix {} does not calculate a global vector", dedn_name,
                   style);
      if (dedn_index > dedn_compute->size_vector)
        error->all(FLERR, "Compute {} for fix {} vector is accessed out-of-range{}", dedn_name,
                   style, utils::errorurl(20));
    }
  } else if (dedn_which == ArgInfo::FIX) {
    dedn_fix = modify->get_fix_by_id(dedn_name);
    if (!dedn_fix)
      error->all(FLERR, "Fix {} for fix {} does not exist", dedn_name, style);
    if (dedn_index == 0 && !dedn_fix->scalar_flag)
      error->all(FLERR, "Fix {} for fix {} does not calculate a global scalar", dedn_name, style);
    if (dedn_index > 0) {
      if (!dedn_fix->vector_flag)
        error->all(FLERR, "Fix {} for fix {} does not calculate a global vector", dedn_name,
                   style);
      if (dedn_index > dedn_fix->size_vector)
        error->all(FLERR, "Fix {} for fix {} vector is accessed out-of-range{}", dedn_name, style,
                   utils::errorurl(20));
    }
  } else {
    error->all(FLERR, "dedn for fix {} must be a variable, compute, or fix reference", style);
  }

  compute_mu_target();
  if (dedn_defer)
    dedn_current = u_target;
  else
    dedn_current = evaluate_dedn();
  u_current = dedn_current;
}

/* ---------------------------------------------------------------------- */

void FixUVT::setup(int vflag)
{
  FixNH::setup(vflag);
  compute_mu_target();
  if (dedn_defer)
    dedn_current = u_target;
  else
    dedn_current = evaluate_dedn();
  u_current = dedn_current;
  *Ne_mass = tdof * boltz * t_target / (u_freq*u_freq);
}

/* ---------------------------------------------------------------------- */

void FixUVT::initial_integrate(int /*vflag*/)
{
  if (pstat_flag && mpchain) nhc_press_integrate();

  if (tstat_flag) {
    compute_temp_target();
    compute_mu_target();
    nhc_mu_integrate();
  }

  if (pstat_flag) {
    if (pstyle == ISO) {
      temperature->compute_scalar();
      pressure->compute_scalar();
    } else {
      temperature->compute_vector();
      pressure->compute_vector();
    }
    couple();
    pressure->addstep(update->ntimestep+1);
  }

  if (pstat_flag) {
    compute_press_target();
    nh_omega_dot();
    nh_v_press();
  }

  nve_v();
  if (pstat_flag) remap();
  nve_x();
  if (pstat_flag) {
    remap();
    if (kspace_flag) force->kspace->setup();
  }
}

/* ---------------------------------------------------------------------- */

void FixUVT::final_integrate()
{
  if (dedn_defer) {
    dedn_current = evaluate_dedn();
    u_current = dedn_current;
  }
  nve_v();

  if (which == BIAS && neighbor->ago == 0)
    t_current = temperature->compute_scalar();

  if (pstat_flag) nh_v_press();

  t_current = temperature->compute_scalar();
  tdof = temperature->dof;

  if (pstat_flag) {
    if (pstyle == ISO) pressure->compute_scalar();
    else {
      temperature->compute_vector();
      pressure->compute_vector();
    }
    couple();
    pressure->addstep(update->ntimestep+1);
  }

  if (pstat_flag) nh_omega_dot();

  if (tstat_flag) nhc_mu_integrate();
  if (pstat_flag && mpchain) nhc_press_integrate();
}

/* ---------------------------------------------------------------------- */

void FixUVT::initial_integrate_respa(int /*vflag*/, int ilevel, int /*iloop*/)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dthalf = 0.5 * step_respa[ilevel];

  if (ilevel == nlevels_respa-1) {
    if (pstat_flag && mpchain) nhc_press_integrate();
    if (tstat_flag) {
      compute_temp_target();
      compute_mu_target();
      nhc_mu_integrate();
    }

    if (pstat_flag) {
      if (pstyle == ISO) {
        temperature->compute_scalar();
        pressure->compute_scalar();
      } else {
        temperature->compute_vector();
        pressure->compute_vector();
      }
      couple();
      pressure->addstep(update->ntimestep+1);
    }

    if (pstat_flag) {
      compute_press_target();
      nh_omega_dot();
      nh_v_press();
    }

    nve_v();
  } else nve_v();

  if (ilevel == 0) {
    if (pstat_flag) remap();
    nve_x();
    if (pstat_flag) remap();
  }
}

/* ---------------------------------------------------------------------- */

void FixUVT::final_integrate_respa(int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) final_integrate();
  else nve_v();
}

/* ---------------------------------------------------------------------- */

double FixUVT::compute_scalar()
{
  double energy = FixNH::compute_scalar();
  double kt = boltz * t_target;
  energy += 0.5*(*Ne_mass)*(*Ne_dot)*(*Ne_dot) + kt * eta[0] - u_target*(*Ne);
  return energy;
}

/* ---------------------------------------------------------------------- */

double FixUVT::compute_vector(int n)
{
  const int base_n = size_vector - 6;
  if (n < base_n) return FixNH::compute_vector(n);
  n -= base_n;

  if (n == 0) return *Ne;
  if (n == 1) return *Ne_dot;
  if (n == 2) return dedn_current;
  if (n == 3) return u_target;
  if (n == 4) return 0.5*(*Ne_mass)*(*Ne_dot)*(*Ne_dot);
  if (n == 5) return -u_target*(*Ne);
  return 0.0;
}

/* ---------------------------------------------------------------------- */

std::string FixUVT::get_thermo_colname(int n)
{
  if (n == -1) return FixNH::get_thermo_colname(n);

  const int base_n = size_vector - 6;
  if (n < base_n) return FixNH::get_thermo_colname(n);
  n -= base_n;

  if (n == 0) return fmt::format("f_{}:Ne", id);
  if (n == 1) return fmt::format("f_{}:Ne_dot", id);
  if (n == 2) return fmt::format("f_{}:dEdN", id);
  if (n == 3) return fmt::format("f_{}:mu", id);
  if (n == 4) return fmt::format("f_{}:ke_Ne", id);
  if (n == 5) return fmt::format("f_{}:pe_mu", id);
  return "none";
}

/* ---------------------------------------------------------------------- */

int FixUVT::size_restart_global()
{
  return FixNH::size_restart_global() + 3;
}

/* ---------------------------------------------------------------------- */

int FixUVT::pack_restart_data(double *list)
{
  int n = FixNH::pack_restart_data(list);
  list[n++] = ustat_flag;
  list[n++] = *Ne;
  list[n++] = *Ne_dot;
  return n;
}

/* ---------------------------------------------------------------------- */

void FixUVT::restart(char *buf)
{
  int n = 0;
  auto *list = (double *) buf;
  int flag = static_cast<int>(list[n++]);
  if (flag) {
    int m = static_cast<int>(list[n++]);
    if (tstat_flag && m == mtchain) {
      for (int ich = 0; ich < mtchain; ich++)
        eta[ich] = list[n++];
      for (int ich = 0; ich < mtchain; ich++)
        eta_dot[ich] = list[n++];
    } else n += 2*m;
  }

  flag = static_cast<int>(list[n++]);
  if (flag) {
    omega[0] = list[n++];
    omega[1] = list[n++];
    omega[2] = list[n++];
    omega[3] = list[n++];
    omega[4] = list[n++];
    omega[5] = list[n++];
    omega_dot[0] = list[n++];
    omega_dot[1] = list[n++];
    omega_dot[2] = list[n++];
    omega_dot[3] = list[n++];
    omega_dot[4] = list[n++];
    omega_dot[5] = list[n++];
    vol0 = list[n++];
    t0 = list[n++];
    int m = static_cast<int>(list[n++]);
    if (pstat_flag && m == mpchain) {
      for (int ich = 0; ich < mpchain; ich++)
        etap[ich] = list[n++];
      for (int ich = 0; ich < mpchain; ich++)
        etap_dot[ich] = list[n++];
    } else n += 2*m;
    flag = static_cast<int>(list[n++]);
    if (flag) {
      h0_inv[0] = list[n++];
      h0_inv[1] = list[n++];
      h0_inv[2] = list[n++];
      h0_inv[3] = list[n++];
      h0_inv[4] = list[n++];
      h0_inv[5] = list[n++];
    }
    flag = static_cast<int>(list[n++]);
    if (flag) {
      p_isoch[0] = static_cast<int>(list[n++]);
      p_isoch[1] = static_cast<int>(list[n++]);
      p_isoch[2] = static_cast<int>(list[n++]);
      vol_start = list[n++];
    }
  }

  flag = static_cast<int>(list[n++]);
  if (flag) {
    *Ne = list[n++];
    *Ne_dot = list[n++];
  } else n += 2;
}

/* ---------------------------------------------------------------------- */

void *FixUVT::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "u_start") == 0) return &u_start;
  if (strcmp(str, "u_stop") == 0) return &u_stop;
  if (strcmp(str, "u_target") == 0) return &u_target;
  if (strcmp(str, "u_current") == 0) return &u_current;
  if (strcmp(str, "u_freq") == 0) return &u_freq;

  dim = 1;
  if (strcmp(str, "ne") == 0) return Ne;
  if (strcmp(str, "ne_dot") == 0) return Ne_dot;
  if (strcmp(str, "ne_mass") == 0) return Ne_mass;
  if (strcmp(str, "dedn") == 0) return &dedn_current;
  return FixNH::extract(str, dim);
}

/* ---------------------------------------------------------------------- */

double FixUVT::memory_usage()
{
  return FixNH::memory_usage();
}

/* ---------------------------------------------------------------------- */

void FixUVT::nve_v()
{
  FixNH::nve_v();
  if (!dedn_defer) {
    dedn_current = evaluate_dedn();
    u_current = dedn_current;
  }
  *Ne_dot += (dthalf / *Ne_mass) * (-dedn_current + u_target);
}

/* ---------------------------------------------------------------------- */

void FixUVT::nve_x()
{
  FixNH::nve_x();
  *Ne += dtv * (*Ne_dot);
}

/* ---------------------------------------------------------------------- */

void FixUVT::nhc_mu_integrate()
{
  int ich;
  double expfac;
  double kecurrent = tdof * boltz * t_current;
  const double ext_ke_target = ke_target + boltz * t_target;

  if (eta_mass_flag) {
    // The first thermostat variable controls the atomistic kinetic energy
    // plus the single quadratic Ne degree of freedom.
    eta_mass[0] = ext_ke_target / (t_freq*t_freq);
    *Ne_mass = tdof * boltz * t_target / (u_freq*u_freq);
    for (ich = 1; ich < mtchain; ich++)
      eta_mass[ich] = boltz * t_target / (t_freq*t_freq);
  }

  if (eta_mass[0] > 0.0)
    eta_dotdot[0] =
      (kecurrent + (*Ne_mass)*(*Ne_dot)*(*Ne_dot) - ext_ke_target) / eta_mass[0];
  else eta_dotdot[0] = 0.0;

  double ncfac = 1.0/nc_tchain;
  for (int iloop = 0; iloop < nc_tchain; iloop++) {
    for (ich = mtchain-1; ich > 0; ich--) {
      expfac = exp(-ncfac*dt8*eta_dot[ich+1]);
      eta_dot[ich] *= expfac;
      eta_dot[ich] += eta_dotdot[ich] * ncfac*dt4;
      eta_dot[ich] *= tdrag_factor;
      eta_dot[ich] *= expfac;
    }

    expfac = exp(-ncfac*dt8*eta_dot[1]);
    eta_dot[0] *= expfac;
    eta_dot[0] += eta_dotdot[0] * ncfac*dt4;
    eta_dot[0] *= tdrag_factor;
    eta_dot[0] *= expfac;

    factor_eta = exp(-ncfac*dthalf*eta_dot[0]);
    nh_v_temp();
    *Ne_dot *= factor_eta;

    t_current *= factor_eta*factor_eta;
    kecurrent = tdof * boltz * t_current;

    if (eta_mass[0] > 0.0)
      eta_dotdot[0] =
        (kecurrent + (*Ne_mass)*(*Ne_dot)*(*Ne_dot) - ext_ke_target) / eta_mass[0];
    else eta_dotdot[0] = 0.0;

    for (ich = 0; ich < mtchain; ich++)
      eta[ich] += ncfac*dthalf*eta_dot[ich];

    eta_dot[0] *= expfac;
    eta_dot[0] += eta_dotdot[0] * ncfac*dt4;
    eta_dot[0] *= expfac;

    for (ich = 1; ich < mtchain; ich++) {
      expfac = exp(-ncfac*dt8*eta_dot[ich+1]);
      eta_dot[ich] *= expfac;
      eta_dotdot[ich] = (eta_mass[ich-1]*eta_dot[ich-1]*eta_dot[ich-1]
                         - boltz * t_target) / eta_mass[ich];
      eta_dot[ich] += eta_dotdot[ich] * ncfac*dt4;
      eta_dot[ich] *= expfac;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixUVT::compute_mu_target()
{
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  u_target = u_start + delta * (u_stop-u_start);
  u_current = u_target;
}

/* ---------------------------------------------------------------------- */

double FixUVT::evaluate_dedn()
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
      if (!(dedn_compute->invoked_flag & Compute::INVOKED_SCALAR)) {
        dedn_compute->compute_scalar();
        dedn_compute->invoked_scalar = update->ntimestep;
        dedn_compute->invoked_flag |= Compute::INVOKED_SCALAR;
      }
      return dedn_compute->scalar;
    }

    if (!(dedn_compute->invoked_flag & Compute::INVOKED_VECTOR)) {
      dedn_compute->compute_vector();
      dedn_compute->invoked_vector = update->ntimestep;
      dedn_compute->invoked_flag |= Compute::INVOKED_VECTOR;
    }
    return dedn_compute->vector[dedn_index - 1];
  }

  if (dedn_index == 0) return dedn_fix->compute_scalar();
  return dedn_fix->compute_vector(dedn_index - 1);
}

/* ---------------------------------------------------------------------- */

void FixUVT::parse_dedn_source(const char *arg)
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
