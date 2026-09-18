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

#include "domain.h"
#include "error.h"
#include "fix_nh.h"
#include "force.h"
#include "group.h"
#include "modify.h"
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

// Thermostat-only composition: this object is never registered with Modify.
// FixPIMDNVT supplies the kinetic energy, chain masses, and timestep. We do not
// call FixNH::init/setup or its particle/barostat integration callbacks.
namespace LAMMPS_NS {
class PIMDNoseHoover : public FixNH {
 public:
  PIMDNoseHoover(LAMMPS *lmp, int narg, char **arg, FixPIMDNVT *owner) :
      FixNH(lmp, narg, arg), owner(owner)
  {
    owner->eta = eta;
    owner->eta_dot = eta_dot;
    owner->eta_dotdot = eta_dotdot;
    owner->eta_mass = eta_mass;
    eta_mass_flag = 0;
    // FixNH::init normally selects this mode. PIMD scales velocities without
    // a temperature-compute bias (NOBIAS = 0 in fix_nh.cpp).
    which = 0;
  }

  void integrate()
  {
    if (!owner->thermostat_chain_active()) {
      // CMD leaves the centroid unthermostatted, but every partition must
      // participate in the work-accounting collective used by nh_v_temp().
      for (int iloop = 0; iloop < owner->nc_tchain; iloop++)
        owner->ecouple_work += owner->thermostat_work_delta(1.0);
      return;
    }

    boltz = force->boltz;
    tdof = owner->tdof;
    nuclear_temperature = owner->compute_nuclear_kinetic_energy() / boltz / tdof;
    update_temperature();
    t_target = owner->np * owner->temp;
    ke_target = owner->chain0_target_energy();
    nc_tchain = owner->nc_tchain;
    tdrag_factor = owner->tdrag_factor;
    dthalf = owner->dthalf;
    dt4 = owner->dt4;
    dt8 = owner->dt8;
    FixNH::nhc_temp_integrate();
  }

 protected:
  void nh_v_temp() override
  {
    const double work_delta = owner->thermostat_work_delta(factor_eta);
    FixNH::nh_v_temp();
    owner->ecouple_work += work_delta;
    owner->thermostat_extra_velocity_step();
    nuclear_temperature *= factor_eta * factor_eta;
    update_temperature();

    // UVT scales the shared electronic velocity using the mean chain friction,
    // which can differ from this bead's nuclear scaling. Supply the updated
    // combined kinetic energy and suppress FixNH's subsequent uniform rescaling.
    factor_eta = 1.0;
  }

 private:
  void update_temperature()
  {
    t_current = (tdof * boltz * nuclear_temperature +
                 owner->thermostat_extra_kinetic_energy()) / (tdof * boltz);
  }

  FixPIMDNVT *owner;
  double nuclear_temperature;
};
}    // namespace LAMMPS_NS

/* ---------------------------------------------------------------------- */

FixPIMDNVT::FixPIMDNVT(LAMMPS *lmp, int narg, char **arg, bool defer_setup) :
    FixPIMDNVE(lmp, narg, arg, true), nhc(nullptr), eta(nullptr), eta_dot(nullptr), eta_dotdot(nullptr),
    eta_mass(nullptr), tau_k(nullptr)
{
  pilescale = 1.0;
  tstat_flag = 1;

  mtchain = 3;
  nc_tchain = 1;
  t_period = 0.0;
  drag = 0.0;

  tdrag_factor = 1.0;
  t_freq = 0.0;
  tdof = 0.0;
  tdof_override_flag = 0;
  tdof_override = 0.0;
  ke_target = 0.0;
  ecouple_work = 0.0;
  dthalf = dt4 = dt8 = 0.0;

  scalar_flag = 1;
  extscalar = 1;
  ecouple_flag = 1;
  thermo_modify_colname = 1;

  if (narg < 4) utils::missing_cmd_args(FLERR, std::string("fix ") + style, error);
  if (defer_setup) return;

  // process keywords

  for (int i = 3; i < narg;) {
    if (!parse_keyword(narg, arg, i))
      error->all(FLERR, "Unknown keyword {} for fix {}", arg[i], style);
  }

  finish_nuclear_constructor_setup();
}

/* ---------------------------------------------------------------------- */

bool FixPIMDNVT::parse_keyword(int narg, char **arg, int &i)
{
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
    if (mtchain < 1)
      error->all(FLERR, i + 1, "Invalid fix {} tchain argument: {}", style, mtchain);
    i += 2;
    return true;
  }
  if (strcmp(arg[i], "tloop") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} tloop", style), error);
    nc_tchain = utils::inumeric(FLERR, arg[i + 1], false, lmp);
    if (nc_tchain < 0)
      error->all(FLERR, i + 1, "Invalid fix {} tloop argument: {}", style, nc_tchain);
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
  return FixPIMDNVE::parse_keyword(narg, arg, i);
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::finish_nuclear_constructor_setup()
{
  if (t_period <= 0.0) error->all(FLERR, "Temperature damping for fix {} must be > 0.0", style);

  if (tstat_flag) {
    // Only ask the FixNH constructor to allocate a chain. Physical parameters
    // are initialized by nhc_init() and supplied to the adapter at each call.
    std::string chain_length = std::to_string(mtchain);
    std::string nh_id = std::string(id) + "_nhc";
    const char *args[] = {nh_id.c_str(), group->names[igroup], "nvt", "temp",
                          "1.0", "1.0", "1.0", "tchain", chain_length.c_str()};
    nhc = new PIMDNoseHoover(lmp, 9, const_cast<char **>(args), this);
    size_vector += 4 * mtchain;
  }

  finish_constructor_setup();
}

/* ---------------------------------------------------------------------- */

FixPIMDNVT::~FixPIMDNVT()
{
  delete[] tau_k;
  delete nhc;
}

/* ---------------------------------------------------------------------- */

double FixPIMDNVT::chain0_target_energy() const
{
  if (!thermostat_chain_active()) return 0.0;
  return static_cast<double>(np) * ke_target;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::setup_subclass_state()
{
  // FixNH propagates for dthalf per call: OBABO has two thermostat
  // half-steps, whereas BAOAB has one full thermostat step.
  dthalf = (integrator == BAOAB ? 1.0 : 0.5) * update->dt;
  dt4 = 0.5 * dthalf;
  dt8 = 0.25 * dthalf;
  nhc_init();
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::o_step()
{
  thermostat_step();
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::nhc_init()
{
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

  if (tstat_flag) {
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
        eta_mass[ich] = chain_target / (omega_np * omega_np);
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

double FixPIMDNVT::compute_nuclear_kinetic_energy() const
{
  // The thermostat equations use twice the kinetic energy.
  double kecurrent = 2.0 * local_kinetic_energy_sum();
  double ketotal = 0.0;
  MPI_Allreduce(&kecurrent, &ketotal, 1, MPI_DOUBLE, MPI_SUM, world);
  return ketotal;
}

/* ---------------------------------------------------------------------- */

bool FixPIMDNVT::thermostat_chain_active() const
{
  if (method == CMD && universe->iworld == 0) return false;
  return true;
}

/* ---------------------------------------------------------------------- */

double FixPIMDNVT::chain_target_energy() const
{
  return thermostat_chain_active() ? static_cast<double>(np) * force->boltz * temp : 0.0;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVT::thermostat_step()
{
  if (tstat_flag) {
    nhc->integrate();
    if (removecomflag) remove_com_motion();
  }
}

/* ---------------------------------------------------------------------- */

double FixPIMDNVT::thermostat_work_delta(double scale_factor) const
{
  const double kinetic_before_local = local_kinetic_energy_sum();
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

double FixPIMDNVT::compute_subclass_vector(int n) const
{
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
