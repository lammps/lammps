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

/* ----------------------------------------------------------------------
   Contributing authors: Yifan Li @ Princeton University (yifanl0716@gmail.com)
                         Li Fu @ Peking University
------------------------------------------------------------------------- */

#include "fix_pimd_uvt.h"

#include "atom.h"
#include "comm.h"
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

/* ---------------------------------------------------------------------- */

FixPIMDUVT::FixPIMDUVT(LAMMPS *lmp, int narg, char **arg) :
    FixPIMDNVT(lmp, narg, arg, true), mu(0.0), Ne(nullptr), Ne_dot(nullptr), Ne_mass(nullptr),
    u_freq(0.0), u_period(0.0), ne_ecouple_work(0.0), dedn_name(nullptr),
    dedn_which(ArgInfo::NONE), dedn_index(0), dedn_var(-1), dedn_compute(nullptr),
    dedn_fix(nullptr), dedn_current(0.0), dedn_bead(0.0)
{
  // process keywords

  bool mu_seen = false;
  for (int i = 3; i < narg;) {
    if (strcmp(arg[i], "mu") == 0) mu_seen = true;
    if (!parse_keyword(narg, arg, i))
      error->all(FLERR, "Unknown keyword {} for fix {}", arg[i], style);
  }

  finish_nuclear_constructor_setup();
  if (method != NMPIMD) error->all(FLERR, "Fix {} only supports method nmpimd", style);
  if (!mu_seen) error->all(FLERR, "Missing mu keyword for fix {}", style);
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

  const int nuclear_size = nuclear_vector_size();
  size_vector += 6;
  delete[] extlist;
  extlist = new int[size_vector];
  // Keep nuclear estimators first, then six physical outputs, then the NHC.
  for (int i = 0; i < size_vector; i++)
    extlist[i] = (i < nuclear_size || i >= nuclear_size + 6) ? 1 : 0;

  u_freq = 1.0 / u_period;
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

bool FixPIMDUVT::parse_keyword(int narg, char **arg, int &i)
{
  if (strcmp(arg[i], "mu") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} mu", style), error);
    mu = utils::numeric(FLERR, arg[i + 1], false, lmp);
    i += 2;
    return true;
  }
  if (strcmp(arg[i], "Udamp") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} Udamp", style), error);
    u_period = utils::numeric(FLERR, arg[i + 1], false, lmp);
    i += 2;
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
  return FixPIMDNVT::parse_keyword(narg, arg, i);
}

/* ---------------------------------------------------------------------- */

void FixPIMDUVT::setup_subclass_state()
{
  FixPIMDNVT::setup_subclass_state();
  resolve_dedn_source();
  modify->addstep_compute_all(update->ntimestep);
  u_freq = 1.0 / u_period;
  // The shared electronic coordinate uses the centroid's nuclear DOF.
  // Internal modes can have different DOF when COM motion is removed.
  if (universe->me == 0)
    *Ne_mass = tdof * force->boltz * temp / (u_freq * u_freq);
  MPI_Bcast(Ne_mass, 1, MPI_DOUBLE, 0, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMDUVT::b_step()
{
  FixPIMDNVE::b_step();

  // Force kicks are half-steps for both OBABO and BAOAB.
  double dtfm = dtv / *Ne_mass;
  if (universe->iworld == 0) *Ne_dot += dtfm * (-dedn_current + mu);
  MPI_Bcast(Ne_dot, 1, MPI_DOUBLE, 0, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMDUVT::qc_step()
{
  FixPIMDNVE::qc_step();

  if (universe->iworld == 0) *Ne += dtv * (*Ne_dot);
  MPI_Bcast(Ne, 1, MPI_DOUBLE, 0, universe->uworld);
}

/* ---------------------------------------------------------------------- */

bool FixPIMDUVT::thermostat_chain_active() const
{
  return true;
}

/* ---------------------------------------------------------------------- */

double FixPIMDUVT::thermostat_extra_kinetic_energy() const
{
  return (*Ne_mass) * (*Ne_dot) * (*Ne_dot);
}

/* ---------------------------------------------------------------------- */

void FixPIMDUVT::scale_ne_velocity(double scale_factor)
{
  const double kinetic_before = 0.5 * (*Ne_mass) * (*Ne_dot) * (*Ne_dot);
  *Ne_dot *= scale_factor;
  const double kinetic_after = 0.5 * (*Ne_mass) * (*Ne_dot) * (*Ne_dot);
  ne_ecouple_work += static_cast<double>(np) * (kinetic_before - kinetic_after);
}

/* ---------------------------------------------------------------------- */

double FixPIMDUVT::chain0_target_energy() const
{
  return FixPIMDNVT::chain0_target_energy() + force->boltz * temp;
}

/* ---------------------------------------------------------------------- */

int FixPIMDUVT::subclass_restart_size() const
{
  return 4;
}

/* ---------------------------------------------------------------------- */

int FixPIMDUVT::pack_subclass_restart(double *list, int n) const
{
  // Preserve the enabled flag in the existing restart format.
  list[n++] = 1.0;
  list[n++] = *Ne;
  list[n++] = *Ne_dot;
  list[n++] = ne_ecouple_work;
  return n;
}

/* ---------------------------------------------------------------------- */

int FixPIMDUVT::unpack_subclass_restart(const double *list, int n)
{
  int flag = static_cast<int>(list[n++]);
  if (flag) {
    *Ne = list[n++];
    *Ne_dot = list[n++];
    ne_ecouple_work = list[n++];
  }
  return n;
}

/* ---------------------------------------------------------------------- */

double FixPIMDUVT::compute_subclass_vector(int n) const
{
  if (n == 0) return *Ne;
  if (n == 1) {
    // Convert the nuclear kinetic energy from PT to physical temperature.
    // Count the shared electronic degree of freedom only once.
    double local[3] = {2.0 * local_kinetic_energy_sum(), comm->me == 0 ? tdof : 0.0,
                       universe->me == 0 ? (*Ne_mass) * (*Ne_dot) * (*Ne_dot) : 0.0};
    double total[3];
    MPI_Allreduce(local, total, 3, MPI_DOUBLE, MPI_SUM, universe->uworld);
    return (total[0] * inverse_np + total[2]) / ((total[1] + 1.0) * force->boltz);
  }
  if (n == 2) return *Ne_dot;
  if (n == 3) return dedn_bead;
  if (n == 4) return 0.5 * (*Ne_mass) * (*Ne_dot) * (*Ne_dot);
  if (n == 5) return -mu * (*Ne);
  return FixPIMDNVT::compute_subclass_vector(n - 6);
}

/* ---------------------------------------------------------------------- */

std::string FixPIMDUVT::get_thermo_colname(int n)
{
  if (n == -1) return FixPIMDNVT::get_thermo_colname(n);
  const int physical_index = n - nuclear_vector_size();
  if (physical_index == 0) return fmt::format("f_{}:Ne", id);
  if (physical_index == 1) return fmt::format("f_{}:T_ins", id);
  if (physical_index == 2) return fmt::format("f_{}:Ne_dot", id);
  if (physical_index == 3) return fmt::format("f_{}:dEdN", id);
  if (physical_index == 4) return fmt::format("f_{}:ke_Ne", id);
  if (physical_index == 5) return fmt::format("f_{}:pe_mu", id);
  return Fix::get_thermo_colname(n);
}

/* ---------------------------------------------------------------------- */

double FixPIMDUVT::compute_scalar()
{
  return ecouple_work + ne_ecouple_work;
}

/* ---------------------------------------------------------------------- */

void *FixPIMDUVT::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str, "ne") == 0) return Ne;
  if (strcmp(str, "ne_dot") == 0) return Ne_dot;
  if (strcmp(str, "ne_mass") == 0) return Ne_mass;
  if (strcmp(str, "dedn") == 0) return &dedn_current;
  return FixPIMDNVT::extract(str, dim);
}

/* ---------------------------------------------------------------------- */

void FixPIMDUVT::thermostat_extra_velocity_step()
{
  // The electronic coordinate is shared by all beads and couples to their
  // mean chain friction, while nuclear velocities use each bead's own chain.
  double eta_dot_k = eta_dot[0];
  double eta_dot_world = 0.0;
  MPI_Allreduce(&eta_dot_k, &eta_dot_world, 1, MPI_DOUBLE, MPI_SUM, world);
  eta_dot_k = (comm->me == 0) ? eta_dot_world / comm->nprocs : 0.0;
  double eta_dot_ave = 0.0;
  MPI_Allreduce(&eta_dot_k, &eta_dot_ave, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  eta_dot_ave *= inverse_np;

  const double ncfac = 1.0 / nc_tchain;
  scale_ne_velocity(exp(-ncfac * dthalf * eta_dot_ave));
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

void FixPIMDUVT::post_force(int vflag)
{
  FixPIMDNVE::post_force(vflag);

  // Evaluate on all ranks: the source may perform collective communication.
  // The source returns a global bead value; contribute it once per partition.
  modify->clearstep_compute();
  dedn_bead = evaluate_dedn();
  double dedn_local = dedn_bead;
  modify->addstep_compute(update->ntimestep + 1);
  if (comm->me != 0) dedn_local = 0.0;
  double dedn_avg = 0.0;
  MPI_Allreduce(&dedn_local, &dedn_avg, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  dedn_current = dedn_avg * inverse_np;
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
