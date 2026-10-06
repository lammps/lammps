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
   Contributing authors: Gareth Tribello (Queens U, Belfast)
                         Pablo Piaggi (EPFL)
------------------------------------------------------------------------- */

#include "fix_plumed.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "pair.h"
#include "respa.h"
#include "timer.h"
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <string>
#include <vector>

#include "plumed/wrapper/Plumed.h"

#if defined(__PLUMED_DEFAULT_KERNEL)
#define PLUMED_QUOTE_DIRECT(name) #name
#define PLUMED_QUOTE(macro) PLUMED_QUOTE_DIRECT(macro)
static const char plumed_default_kernel[] = "PLUMED_KERNEL=" PLUMED_QUOTE(__PLUMED_DEFAULT_KERNEL);
#endif

/* -------------------------------------------------------------------- */

using namespace LAMMPS_NS;
using namespace FixConst;

namespace {
enum PlumedTraceField {
  PLUMED_TRACE_NONE,
  PLUMED_TRACE_BOX,
  PLUMED_TRACE_POSITION,
  PLUMED_TRACE_VELOCITY,
  PLUMED_TRACE_PHYSICAL_FORCE,
  PLUMED_TRACE_FORCE_DELTA,
  PLUMED_TRACE_TOTAL_FORCE,
  PLUMED_TRACE_CENTROID_POSITION,
  PLUMED_TRACE_CENTROID_FORCE
};

struct PlumedTraceRecord {
  int found;
  int field;
  int component;
  int universe_rank;
  int bead_world;
  int world_rank;
  tagint tag;
  double value;
  double position[3];
  double velocity[3];
  double physical_force[3];
  double force_delta[3];
  double total_force[3];
  double box[9];
};

const char *plumed_trace_field_name(int field)
{
  if (field == PLUMED_TRACE_BOX) return "box";
  if (field == PLUMED_TRACE_POSITION) return "position";
  if (field == PLUMED_TRACE_VELOCITY) return "velocity";
  if (field == PLUMED_TRACE_PHYSICAL_FORCE) return "physical-force";
  if (field == PLUMED_TRACE_FORCE_DELTA) return "plumed-force-delta";
  if (field == PLUMED_TRACE_TOTAL_FORCE) return "total-force";
  if (field == PLUMED_TRACE_CENTROID_POSITION) return "centroid-position";
  if (field == PLUMED_TRACE_CENTROID_FORCE) return "centroid-force";
  return "unknown";
}
}    // namespace

FixPlumed::FixPlumed(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), p(nullptr), pimd_fix(nullptr), nlocal(-1), natoms(0),
    path_integral_mode(PATH_INTEGRAL_OFF), plumed_active(1), centroid_force_scale(0.0),
    bead_density_force_scale(0.0), gatindex(nullptr), masses(nullptr), charges(nullptr),
    centroid_coordinates(nullptr), centroid_positions(nullptr), centroid_forces(nullptr),
    centroid_forces_all(nullptr), centroid_virial_pending(nullptr),
    bead_bias_virial_pending(nullptr), forces_before_plumed(nullptr), nlevels_respa(0), bias(0.0),
    c_pe(nullptr), c_press(nullptr), plumedNeedsEnergy(0), id_pe(nullptr), id_press(nullptr),
    nonfinite_trace_prefix(nullptr), id_pimd(nullptr)
{

  if (!atom->tag_enable) error->all(FLERR, "Fix plumed requires atom tags");

  if (atom->tag_consecutive() == 0) error->all(FLERR, "Fix plumed requires consecutive atom IDs");

  if (igroup != 0 && comm->me == 0)
    error->warning(FLERR, "Fix group for fix plumed is not 'all'. Group will be ignored.");

  const char *plumedfile = nullptr;
  const char *outfile = nullptr;
  for (int i = 3; i < narg; i += 2) {
    if (i + 1 >= narg) utils::missing_cmd_args(FLERR, "fix plumed", error);
    if (strcmp(arg[i], "plumedfile") == 0) {
      plumedfile = arg[i + 1];
    } else if (strcmp(arg[i], "outfile") == 0) {
      outfile = arg[i + 1];
    } else if (strcmp(arg[i], "path_integral") == 0) {
      if (strcmp(arg[i + 1], "off") == 0)
        path_integral_mode = PATH_INTEGRAL_OFF;
      else if (strcmp(arg[i + 1], "centroid") == 0)
        path_integral_mode = PATH_INTEGRAL_CENTROID;
      else if (strcmp(arg[i + 1], "bead_mean") == 0)
        path_integral_mode = PATH_INTEGRAL_BEAD_MEAN;
      else if (strcmp(arg[i + 1], "bead_density") == 0)
        path_integral_mode = PATH_INTEGRAL_BEAD_DENSITY;
      else
        error->all(FLERR, "Unknown fix plumed path_integral value: {}", arg[i + 1]);
    } else if (strcmp(arg[i], "pimd_fix") == 0) {
      delete[] id_pimd;
      id_pimd = utils::strdup(arg[i + 1]);
    } else if (strcmp(arg[i], "nonfinite_trace") == 0) {
      delete[] nonfinite_trace_prefix;
      nonfinite_trace_prefix = utils::strdup(arg[i + 1]);
    } else {
      error->all(FLERR, "Unknown fix plumed keyword: {}", arg[i]);
    }
  }
  if (path_integral_mode != PATH_INTEGRAL_OFF && id_pimd == nullptr)
    error->all(FLERR, "Fix plumed path_integral mode requires the pimd_fix keyword");
  if (path_integral_mode == PATH_INTEGRAL_OFF && id_pimd != nullptr)
    error->all(FLERR, "Fix plumed pimd_fix requires a path_integral mode");

  plumed_active = (path_integral_mode != PATH_INTEGRAL_CENTROID) || (universe->iworld == 0);

#if defined(__PLUMED_DEFAULT_KERNEL)
  if (getenv("PLUMED_KERNEL") == nullptr) platform::putenv(plumed_default_kernel);
#endif

  // Check API version

  int api_version = 0;
  if (plumed_active) {
    try {
      p = new PLMD::Plumed;
      p->cmd("getApiVersion", &api_version);
    } catch (const std::exception &exception) {
      error->universe_one(FLERR, fmt::format("Could not initialize PLUMED: {}", exception.what()));
    }
  }
  if (path_integral_mode == PATH_INTEGRAL_CENTROID)
    MPI_Bcast(&api_version, 1, MPI_INT, 0, universe->uworld);
  if ((api_version < 5) || (api_version > 11))
    error->all(
        FLERR,
        "Incompatible API version for PLUMED in fix plumed. "
        "Only Plumed 2.4.x, 2.5.x, 2.6.x, 2.7.x, 2.8.x, 2.9.x, 2.10.x are tested and supported.");

#if !defined(MPI_STUBS)
  // If the -partition option is activated then enable
  // inter-partition communication

  try {
    if ((plumed_active) && (path_integral_mode != PATH_INTEGRAL_CENTROID) &&
        (universe->existflag == 1)) {
      MPI_Comm inter_comm;

      // Change MPI_COMM_WORLD to universe->uworld which seems more appropriate

      MPI_Comm_split(universe->uworld, comm->me, 0, &inter_comm);
      p->cmd("GREX setMPIIntracomm", &world);
      if (comm->me == 0) {
        // The inter-partition communicator is only defined for the root in
        //    each partition (a.k.a. world). This is due to the way in which
        //    it is defined inside plumed.
        p->cmd("GREX setMPIIntercomm", &inter_comm);
      }
      p->cmd("GREX init", nullptr);
      // PLUMED duplicates the communicator passed to setMPIIntercomm.
      MPI_Comm_free(&inter_comm);
    }

    // The general communicator is independent of the existence of partitions,
    // if there are partitions, world is defined within each partition,
    // whereas if partitions are not defined then world is equal to
    // MPI_COMM_WORLD.

    // plumed does not know about LAMMPS using the MPI STUBS library and will
    // fail if this is called under these circumstances
    if (plumed_active) p->cmd("setMPIComm", &world);
  } catch (const std::exception &exception) {
    error->universe_one(FLERR, fmt::format("Could not configure PLUMED MPI: {}", exception.what()));
  }
#endif

  // Set up units
  // LAMMPS units wrt kj/mol - nm - ps
  // Set up units

  if (plumed_active && strcmp(update->unit_style, "lj") == 0) {
    // LAMMPS units lj
    try {
      p->cmd("setNaturalUnits");
    } catch (const std::exception &exception) {
      error->universe_one(FLERR,
                          fmt::format("Could not configure PLUMED units: {}", exception.what()));
    }
  } else if (plumed_active) {

    // Conversion factor from LAMMPS energy units to kJ/mol (units of PLUMED)

    double energyUnits = 1.0;

    // LAMMPS units real :: kcal/mol;

    if (strcmp(update->unit_style, "real") == 0) {
      energyUnits = 4.184;

      // LAMMPS units metal :: eV;

    } else if (strcmp(update->unit_style, "metal") == 0) {
      energyUnits = 96.48530749925792;

      // LAMMPS units si :: Joule;

    } else if (strcmp(update->unit_style, "si") == 0) {
      energyUnits = 0.001;

      // LAMMPS units cgs :: erg;

    } else if (strcmp(update->unit_style, "cgs") == 0) {
      energyUnits = 6.0221418e13;

      // LAMMPS units electron :: Hartree;

    } else if (strcmp(update->unit_style, "electron") == 0) {
      energyUnits = 2625.5257;

    } else
      error->all(FLERR, "Fix plumed cannot handle {} units", update->unit_style);

    // Conversion factor from LAMMPS length units to nm (units of PLUMED)

    double lengthUnits = 0.1 / force->angstrom;

    // Conversion factor from LAMMPS time unit to ps (units of PLUMED)

    double timeUnits = 0.001 / force->femtosecond;

    try {
      p->cmd("setMDEnergyUnits", &energyUnits);
      p->cmd("setMDLengthUnits", &lengthUnits);
      p->cmd("setMDTimeUnits", &timeUnits);
    } catch (const std::exception &exception) {
      error->universe_one(FLERR,
                          fmt::format("Could not configure PLUMED units: {}", exception.what()));
    }
  }

  if (atom->natoms > MAXSMALLINT)
    error->all(FLERR, "Fix plumed can only handle up to 2.1 billion atoms");

  natoms = int(atom->natoms);
  double dt = update->dt;

  if (plumed_active) {
    try {
      if (outfile) {
        if ((path_integral_mode != PATH_INTEGRAL_CENTROID) && (universe->existflag == 1))
          p->cmd("setLogFile", fmt::format("{}.{}", outfile, universe->iworld).c_str());
        else
          p->cmd("setLogFile", outfile);
      }
      if (plumedfile) p->cmd("setPlumedDat", plumedfile);
      p->cmd("setMDEngine", "LAMMPS");
      p->cmd("setNatoms", &natoms);
      p->cmd("setTimestep", &dt);
      p->cmd("init");
    } catch (const std::exception &exception) {
      error->universe_one(FLERR, fmt::format("Could not configure PLUMED: {}", exception.what()));
    }
  }

  extscalar = 1;
  scalar_flag = 1;
  energy_global_flag = virial_global_flag = 1;
  thermo_energy = thermo_virial = 1;

  // Define compute to calculate potential energy

  if (path_integral_mode != PATH_INTEGRAL_CENTROID) {
    delete[] id_pe;
    id_pe = utils::strdup("plmd_pe");
    c_pe = modify->add_compute(std::string(id_pe) + " all pe");

    // Define compute to calculate pressure tensor

    delete[] id_press;
    id_press = utils::strdup("plmd_press");
    c_press = modify->add_compute(std::string(id_press) + " all pressure NULL virial");
  }

  for (const auto &fix : modify->get_fix_list()) {
    const char *const check_style = fix->style;

    // There must be only one

    if (strcmp(check_style, "plumed") == 0)
      error->all(FLERR, "There must be only one instance of fix plumed");

    // Avoid conflict with fixes that define internal pressure computes.
    // See comment in the setup method

    if (utils::strmatch(check_style, "^nph") || utils::strmatch(check_style, "^press/berendsen") ||
        utils::strmatch(check_style, "^npt") || utils::strmatch(check_style, "^deform/pressure") ||
        utils::strmatch(check_style, "^msst") || utils::strmatch(check_style, "^press/langevin") ||
        utils::strmatch(check_style, "^nphug") || utils::strmatch(check_style, "^tgnpt/drude") ||
        utils::strmatch(check_style, "^qbmsst") || utils::strmatch(check_style, "^rigid/nph") ||
        utils::strmatch(check_style, "^rigid/npt") || utils::strmatch(check_style, "^box/relax"))
      error->all(FLERR,
                 "Fix plumed must be defined before any other fixes like fix {} that compute "
                 "pressure internally",
                 check_style);
  }
  check_path_integral_compatibility();
}

/* ---------------------------------------------------------------------- */

void FixPlumed::check_path_integral_compatibility()
{
  int has_path_integral_fix = 0;
  for (const auto &fix : modify->get_fix_list())
    if (utils::strmatch(fix->style, "^pimd") || utils::strmatch(fix->style, "^ipi"))
      has_path_integral_fix = 1;

  if (path_integral_mode == PATH_INTEGRAL_OFF && has_path_integral_fix)
    error->all(FLERR, "Fix plumed is incompatible with path-integrals");
  if (path_integral_mode != PATH_INTEGRAL_OFF) {
    pimd_fix = modify->get_fix_by_id(id_pimd);
    if (!pimd_fix || strcmp(pimd_fix->style, "pimd/langevin") != 0)
      error->all(FLERR, "Fix plumed pimd_fix {} must use style pimd/langevin", id_pimd);
  }
}

/* ---------------------------------------------------------------------- */

void FixPlumed::check_normal_mode_post_force_order()
{
  if (!centroid_virial_pending && !bead_bias_virial_pending) return;

  int plumed_index = -1;
  for (int i = 0; i < modify->nfix; ++i) {
    if (modify->fix[i] == this) {
      plumed_index = i;
      break;
    }
  }

  for (int i = plumed_index + 1; i < modify->nfix; ++i) {
    if ((modify->fmask[i] & POST_FORCE) && strcmp(modify->fix[i]->style, "GROUP") != 0)
      error->all(FLERR,
                 "Fix plumed with NMPIMD path_integral modes must be defined after fix {} "
                 "because it has a post-force callback",
                 modify->fix[i]->style);
  }
}

FixPlumed::~FixPlumed()
{
  delete p;
  if (id_pe) modify->delete_compute(id_pe);
  if (id_press) modify->delete_compute(id_press);
  delete[] id_pe;
  delete[] id_press;
  delete[] id_pimd;
  delete[] nonfinite_trace_prefix;
  delete[] masses;
  delete[] charges;
  delete[] gatindex;
  delete[] centroid_positions;
  delete[] centroid_forces;
  delete[] centroid_forces_all;
  delete[] forces_before_plumed;
}

int FixPlumed::setmask()
{
  // set with a bitmask how and when apply the force from plumed
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

void FixPlumed::init()
{
  check_path_integral_compatibility();

  if (path_integral_mode != PATH_INTEGRAL_OFF) {
    if (utils::strmatch(update->integrate_style, "^respa"))
      error->all(FLERR, "Fix plumed path_integral modes do not support r-RESPA");

    int dim = -1;
    double *force_scale = nullptr;
    if (path_integral_mode == PATH_INTEGRAL_CENTROID) {
      force_scale = static_cast<double *>(pimd_fix->extract("centroid_bias_force_scale", dim));
      if (!force_scale) {
        force_scale =
            static_cast<double *>(pimd_fix->extract("normal_mode_centroid_force_scale", dim));
        centroid_virial_pending =
            static_cast<int *>(pimd_fix->extract("centroid_bias_virial_pending", dim));
        if (force_scale && !centroid_virial_pending)
          error->all(FLERR, "Fix plumed could not access the NMPIMD centroid virial state");
      }
    } else {
      force_scale = static_cast<double *>(pimd_fix->extract("bead_bias_force_scale", dim));
      int deferred_dim = -1;
      auto *defer_normal_mode_force =
          static_cast<int *>(pimd_fix->extract("defer_normal_mode_force", deferred_dim));
      if (defer_normal_mode_force) {
        int pending_dim = -1;
        bead_bias_virial_pending =
            static_cast<int *>(pimd_fix->extract("bead_bias_virial_pending", pending_dim));
        if (deferred_dim != 0 || !bead_bias_virial_pending || pending_dim != 0)
          error->all(FLERR, "Fix plumed could not access the NMPIMD bead-force state");
        *defer_normal_mode_force = 1;
      }
    }
    if (!force_scale || dim != 0) {
      if (path_integral_mode == PATH_INTEGRAL_CENTROID)
        error->all(FLERR,
                   "Fix plumed path_integral centroid requires method pimd with ensemble nvt or "
                   "method nmpimd with ensemble nvt, nph, or npt");
      error->all(FLERR,
                 "Fix plumed path_integral bead modes require method pimd with ensemble nvt or "
                 "method nmpimd with ensemble nvt, nph, or npt");
    }
    auto *beads = static_cast<int *>(pimd_fix->extract("nbeads", dim));
    if (!beads || dim != 0 || *beads != universe->nworlds)
      error->all(FLERR, "Fix plumed could not determine the PIMD bead count");
    check_normal_mode_post_force_order();

    if (path_integral_mode == PATH_INTEGRAL_BEAD_MEAN ||
        path_integral_mode == PATH_INTEGRAL_BEAD_DENSITY) {
      if (universe->existflag == 0 || universe->nworlds < 2)
        error->all(FLERR, "Fix plumed path_integral bead modes require multiple partitions");
      if (path_integral_mode == PATH_INTEGRAL_BEAD_DENSITY) bead_density_force_scale = *force_scale;
    } else {
      // The integrator samples H_ring at beta/P. A physical bias B therefore
      // enters the dynamical Hamiltonian as P*B; retain B as the reported bias.
      centroid_force_scale = *force_scale * (*beads);
      centroid_coordinates = static_cast<double *>(pimd_fix->extract("centroid_coordinates", dim));
      if (!centroid_coordinates || dim != 1)
        error->all(FLERR, "Fix plumed could not access the PIMD centroid coordinates");

      delete[] centroid_forces_all;
      centroid_forces_all = new double[3 * natoms]();
      if (plumed_active) {
        delete[] centroid_positions;
        delete[] centroid_forces;
        centroid_positions = new double[3 * natoms];
        centroid_forces = new double[3 * natoms];
      }
      for (int i = 0; i < 6; i++) virial[i] = 0.0;
      return;
    }
  }

  if (utils::strmatch(update->integrate_style, "^respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // This avoids nan pressure if compute_pressure is called
  // in a setup method

  for (int i = 0; i < 6; i++) virial[i] = 0.;

  c_pe = modify->get_compute_by_id(id_pe);
  if (!c_pe) {
    error->all(FLERR, "Potential energy compute ID {} for fix plumed does not exist", id_pe);
  } else {
    if (c_pe->peflag == 0)
      error->all(FLERR, "Compute ID {} does not compute potential energy", id_pe);
  }

  c_press = modify->get_compute_by_id(id_press);
  if (!c_press) {
    error->all(FLERR, "Pressure compute ID {} for fix plumed does not exist", id_press);
  } else {
    if (c_press->pressflag == 0)
      error->all(FLERR, "Compute ID {} does not compute pressure", id_press);
  }
}

void FixPlumed::setup(int vflag)
{
  // Here there is a crucial issue connected to constant pressure
  // simulations. The fix_nh will call the compute_pressure inside
  // the setup method, that is executed once and for all at the
  // beginning of the simulation. Since our fix has a contribution
  // to the virial, when this happens the variable virial must have
  // been calculated. In other words, the setup method of fix_plumed
  // has to be executed first. This creates a race condition with the
  // setup method of fix_nh. This is why in the constructor I check if
  // nh fixes have already been called.
  if (utils::strmatch(update->integrate_style, "^respa")) {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa - 1);
    post_force_respa(vflag, nlevels_respa - 1, 0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa - 1);
  } else {
    post_force(vflag);
  }
  if (centroid_virial_pending || bead_bias_virial_pending) pimd_fix->end_of_step();
}

void FixPlumed::min_setup(int vflag)
{
  if (path_integral_mode != PATH_INTEGRAL_OFF)
    error->all(FLERR, "Fix plumed path_integral modes do not support minimization");
  // This has to be checked.
  // For instance it might have problems with fix_box_relax
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPlumed::update_atom_data()
{
  int update_gatindex = 0;

  if (natoms != int(atom->natoms))
    error->all(FLERR, "Fix plumed does not support simulations with varying numbers of atoms");

  if (nlocal != atom->nlocal) {
    delete[] charges;
    delete[] masses;
    delete[] gatindex;
    delete[] forces_before_plumed;

    nlocal = atom->nlocal;
    const int local_capacity = nlocal > 0 ? nlocal : 1;
    gatindex = new int[local_capacity];
    masses = new double[local_capacity];
    charges = new double[local_capacity];
    forces_before_plumed = (path_integral_mode == PATH_INTEGRAL_BEAD_DENSITY ||
                            path_integral_mode == PATH_INTEGRAL_BEAD_MEAN)
        ? new double[nlocal > 0 ? 3 * nlocal : 1]
        : nullptr;
    update_gatindex = 1;
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (gatindex[i] != atom->tag[i] - 1) {
        update_gatindex = 1;
        break;
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &update_gatindex, 1, MPI_INT, MPI_SUM, world);

  if (update_gatindex) {
    for (int i = 0; i < nlocal; i++) gatindex[i] = atom->tag[i] - 1;
    if (atom->rmass_flag) {
      for (int i = 0; i < nlocal; i++) masses[i] = atom->rmass[i];
    } else {
      for (int i = 0; i < nlocal; i++) masses[i] = atom->mass[atom->type[i]];
    }
    if (atom->q_flag) {
      for (int i = 0; i < nlocal; i++) charges[i] = atom->q[i];
    } else {
      for (int i = 0; i < nlocal; i++) charges[i] = 0.0;
    }
    try {
      p->cmd("setAtomsNlocal", &nlocal);
      p->cmd("setAtomsGatindex", gatindex);
    } catch (const std::exception &exception) {
      error->universe_one(FLERR,
                          fmt::format("Could not update PLUMED atom data: {}", exception.what()));
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPlumed::trace_nonfinite_state(const char *stage, const double *force_before,
                                      bool current_force_is_physical)
{
  if (!nonfinite_trace_prefix) return;

  PlumedTraceRecord local{};
  local.universe_rank = universe->me;
  local.bead_world = universe->iworld;
  local.world_rank = comm->me;
  local.tag = -1;

  const double box[9] = {domain->boxlo[0], domain->boxlo[1], domain->boxlo[2],
                         domain->boxhi[0], domain->boxhi[1], domain->boxhi[2],
                         domain->xy,       domain->xz,       domain->yz};
  for (int i = 0; i < 9; ++i) {
    local.box[i] = box[i];
    if (!local.found && !std::isfinite(box[i])) {
      local.found = 1;
      local.field = PLUMED_TRACE_BOX;
      local.component = i;
      local.value = box[i];
    }
  }
  if (!local.found && (!(domain->xprd > 0.0) || !(domain->yprd > 0.0) || !(domain->zprd > 0.0))) {
    local.found = 1;
    local.field = PLUMED_TRACE_BOX;
    if (!(domain->xprd > 0.0)) {
      local.component = 0;
      local.value = domain->xprd;
    } else if (!(domain->yprd > 0.0)) {
      local.component = 1;
      local.value = domain->yprd;
    } else {
      local.component = 2;
      local.value = domain->zprd;
    }
  }

  int local_index = -1;
  auto select_atom_field = [&](int index, int field, int component, double value) {
    const tagint tag = atom->tag[index];
    if (!local.found || local.field == PLUMED_TRACE_BOX || tag < local.tag ||
        (tag == local.tag &&
         (field < local.field || (field == local.field && component < local.component)))) {
      if (local.field == PLUMED_TRACE_BOX) return;
      local.found = 1;
      local.field = field;
      local.component = component;
      local.tag = tag;
      local.value = value;
      local_index = index;
    }
  };

  if (!local.found) {
    for (int i = 0; i < atom->nlocal; ++i) {
      for (int d = 0; d < 3; ++d) {
        const double physical = force_before ? force_before[3 * i + d] : atom->f[i][d];
        const double delta = force_before ? atom->f[i][d] - physical : 0.0;
        if (!std::isfinite(atom->x[i][d]))
          select_atom_field(i, PLUMED_TRACE_POSITION, d, atom->x[i][d]);
        if (!std::isfinite(atom->v[i][d]))
          select_atom_field(i, PLUMED_TRACE_VELOCITY, d, atom->v[i][d]);
        if (force_before && !std::isfinite(physical))
          select_atom_field(i, PLUMED_TRACE_PHYSICAL_FORCE, d, physical);
        if (force_before && !std::isfinite(delta))
          select_atom_field(i, PLUMED_TRACE_FORCE_DELTA, d, delta);
        if (!std::isfinite(atom->f[i][d]))
          select_atom_field(
              i, current_force_is_physical ? PLUMED_TRACE_PHYSICAL_FORCE : PLUMED_TRACE_TOTAL_FORCE,
              d, atom->f[i][d]);
      }
    }
  }

  if (local_index >= 0) {
    for (int d = 0; d < 3; ++d) {
      local.position[d] = atom->x[local_index][d];
      local.velocity[d] = atom->v[local_index][d];
      local.physical_force[d] = force_before ? force_before[3 * local_index + d]
                                             : atom->f[local_index][d];
      local.force_delta[d] =
          force_before ? atom->f[local_index][d] - force_before[3 * local_index + d] : 0.0;
      local.total_force[d] = atom->f[local_index][d];
    }
  }

  int any_nonfinite = local.found;
  MPI_Allreduce(MPI_IN_PLACE, &any_nonfinite, 1, MPI_INT, MPI_MAX, universe->uworld);
  if (!any_nonfinite) return;

  std::vector<PlumedTraceRecord> records(universe->nprocs);
  MPI_Allgather(&local, sizeof(PlumedTraceRecord), MPI_BYTE, records.data(),
                sizeof(PlumedTraceRecord), MPI_BYTE, universe->uworld);
  const PlumedTraceRecord *winner = nullptr;
  for (const auto &record : records) {
    if (!record.found) continue;
    if (!winner) {
      winner = &record;
      continue;
    }
    const int record_box = record.field == PLUMED_TRACE_BOX;
    const int winner_box = winner->field == PLUMED_TRACE_BOX;
    if (record_box != winner_box) {
      if (record_box) winner = &record;
      continue;
    }
    if (record.bead_world != winner->bead_world) {
      if (record.bead_world < winner->bead_world) winner = &record;
      continue;
    }
    if (record.tag != winner->tag) {
      if (record.tag < winner->tag) winner = &record;
      continue;
    }
    if (record.field != winner->field) {
      if (record.field < winner->field) winner = &record;
      continue;
    }
    if (record.component != winner->component) {
      if (record.component < winner->component) winner = &record;
      continue;
    }
    if (record.universe_rank < winner->universe_rank) winner = &record;
  }

  if (!winner) return;
  const std::string path =
      fmt::format("{}.plumed.step{}.u{}.w{}.r{}.txt", nonfinite_trace_prefix, update->ntimestep,
                  winner->universe_rank, winner->bead_world, winner->world_rank);
  if (universe->me == winner->universe_rank) {
    const std::string report = fmt::format(
        "schema=plumed-nonfinite-trace-v1\nstep={}\nstage={}\nuniverse_rank={}\n"
        "bead_world={}\nworld_rank={}\natom_tag={}\nfield={}\ncomponent={}\n"
        "value={:.17g}\nposition={:.17g} {:.17g} {:.17g}\n"
        "velocity={:.17g} {:.17g} {:.17g}\nphysical_force={:.17g} {:.17g} {:.17g}\n"
        "plumed_force_delta={:.17g} {:.17g} {:.17g}\ntotal_force={:.17g} {:.17g} {:.17g}\n"
        "boxlo={:.17g} {:.17g} {:.17g}\nboxhi={:.17g} {:.17g} {:.17g}\n"
        "tilt_xy_xz_yz={:.17g} {:.17g} {:.17g}\n",
        update->ntimestep, stage, winner->universe_rank, winner->bead_world, winner->world_rank,
        winner->tag, plumed_trace_field_name(winner->field), winner->component, winner->value,
        winner->position[0], winner->position[1], winner->position[2], winner->velocity[0],
        winner->velocity[1], winner->velocity[2], winner->physical_force[0],
        winner->physical_force[1], winner->physical_force[2], winner->force_delta[0],
        winner->force_delta[1], winner->force_delta[2], winner->total_force[0],
        winner->total_force[1], winner->total_force[2], winner->box[0], winner->box[1],
        winner->box[2], winner->box[3], winner->box[4], winner->box[5], winner->box[6],
        winner->box[7], winner->box[8]);
    if (FILE *file = std::fopen(path.c_str(), "w")) {
      std::fwrite(report.data(), 1, report.size(), file);
      std::fclose(file);
    }
  }
  error->all(FLERR,
             fmt::format("Fix plumed nonfinite trace detected a {} at step {} stage {} on bead {} "
                         "atom {} component {}; diagnostic {}",
                         plumed_trace_field_name(winner->field), update->ntimestep, stage,
                         winner->bead_world, winner->tag, winner->component, path));
}

/* ---------------------------------------------------------------------- */

void FixPlumed::trace_centroid_nonfinite(const char *stage)
{
  if (!nonfinite_trace_prefix) return;

  PlumedTraceRecord local{};
  local.universe_rank = universe->me;
  local.bead_world = universe->iworld;
  local.world_rank = comm->me;
  local.tag = -1;
  int local_index = -1;
  if (plumed_active) {
    for (int i = 0; i < nlocal; ++i) {
      const tagint tag = atom->tag[i];
      for (int d = 0; d < 3; ++d) {
        const double position = centroid_positions[3 * i + d];
        const double force = centroid_forces[3 * i + d];
        const int field = !std::isfinite(position) ? PLUMED_TRACE_CENTROID_POSITION
                                                   : PLUMED_TRACE_CENTROID_FORCE;
        const double value = field == PLUMED_TRACE_CENTROID_POSITION ? position : force;
        if ((std::isfinite(position) && std::isfinite(force)) ||
            (local.found &&
             (tag > local.tag ||
              (tag == local.tag &&
               (field > local.field || (field == local.field && d >= local.component))))))
          continue;
        local.found = 1;
        local.field = field;
        local.component = d;
        local.tag = tag;
        local.value = value;
        local_index = i;
      }
    }
  }
  if (local_index >= 0) {
    for (int d = 0; d < 3; ++d) {
      local.position[d] = centroid_positions[3 * local_index + d];
      local.total_force[d] = centroid_forces[3 * local_index + d];
    }
  }

  int any_nonfinite = local.found;
  MPI_Allreduce(MPI_IN_PLACE, &any_nonfinite, 1, MPI_INT, MPI_MAX, universe->uworld);
  if (!any_nonfinite) return;

  std::vector<PlumedTraceRecord> records(universe->nprocs);
  MPI_Allgather(&local, sizeof(PlumedTraceRecord), MPI_BYTE, records.data(),
                sizeof(PlumedTraceRecord), MPI_BYTE, universe->uworld);
  const PlumedTraceRecord *winner = nullptr;
  for (const auto &record : records) {
    if (!record.found) continue;
    if (!winner || record.tag < winner->tag ||
        (record.tag == winner->tag &&
         (record.field < winner->field ||
          (record.field == winner->field && record.component < winner->component))))
      winner = &record;
  }
  if (!winner) return;

  const std::string path =
      fmt::format("{}.plumed.step{}.u{}.w{}.r{}.txt", nonfinite_trace_prefix, update->ntimestep,
                  winner->universe_rank, winner->bead_world, winner->world_rank);
  if (universe->me == winner->universe_rank) {
    const std::string report = fmt::format(
        "schema=plumed-nonfinite-trace-v1\nstep={}\nstage={}\nuniverse_rank={}\n"
        "bead_world={}\nworld_rank={}\natom_tag={}\nfield={}\ncomponent={}\n"
        "value={:.17g}\ncentroid_position={:.17g} {:.17g} {:.17g}\n"
        "centroid_force={:.17g} {:.17g} {:.17g}\n",
        update->ntimestep, stage, winner->universe_rank, winner->bead_world, winner->world_rank,
        winner->tag, plumed_trace_field_name(winner->field), winner->component, winner->value,
        winner->position[0], winner->position[1], winner->position[2], winner->total_force[0],
        winner->total_force[1], winner->total_force[2]);
    if (FILE *file = std::fopen(path.c_str(), "w")) {
      std::fwrite(report.data(), 1, report.size(), file);
      std::fclose(file);
    }
  }
  error->all(FLERR,
             fmt::format("Fix plumed nonfinite trace detected a {} at step {} stage {} on atom {} "
                         "component {}; diagnostic {}",
                         plumed_trace_field_name(winner->field), update->ntimestep, stage,
                         winner->tag, winner->component, path));
}

/* ---------------------------------------------------------------------- */

void FixPlumed::post_force(int /* vflag */)
{
  trace_nonfinite_state("pre-plumed", nullptr, true);

  if (path_integral_mode == PATH_INTEGRAL_CENTROID) {
    post_force_centroid();
    return;
  }

  update_atom_data();

  if (path_integral_mode == PATH_INTEGRAL_BEAD_DENSITY ||
      path_integral_mode == PATH_INTEGRAL_BEAD_MEAN)
    for (int i = 0; i < nlocal; i++)
      for (int d = 0; d < 3; d++) forces_before_plumed[3 * i + d] = atom->f[i][d];

  // set up local virial/box. plumed uses full 3x3 matrices
  double plmd_virial[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) plmd_virial[i][j] = 0.0;
  double box[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) box[i][j] = 0.0;
  box[0][0] = domain->h[0];
  box[1][1] = domain->h[1];
  box[2][2] = domain->h[2];
  box[2][1] = domain->h[3];
  box[2][0] = domain->h[4];
  box[1][0] = domain->h[5];

  // Make initial of virial of this fix zero
  // The following line is very important, otherwise
  // the compute pressure will include
  for (int i = 0; i < 6; ++i) virial[i] = 0.;

  // local variable with timestep:
  if (update->ntimestep > MAXSMALLINT)
    error->all(FLERR, "Fix plumed can only handle up to 2.1 billion timesteps");
  int step = int(update->ntimestep);

  // pass all pointers to plumed:
  p->cmd("setStep", &step);
  int plumedStopCondition = 0;
  p->cmd("setStopFlag", &plumedStopCondition);
  p->cmd("setPositions", &atom->x[0][0]);
  p->cmd("setBox", &box[0][0]);
  p->cmd("setForces", &atom->f[0][0]);
  p->cmd("setMasses", &masses[0]);
  p->cmd("setCharges", &charges[0]);
  p->cmd("getBias", &bias);

  // Pass virial to plumed
  // If energy is needed plmd_virial is equal to LAMMPS' virial
  // If energy is not needed plmd_virial is initialized to zero
  // In the first case the virial will be rescaled and an extra term will be added
  // In the latter case only an extra term will be added
  p->cmd("setVirial", &plmd_virial[0][0]);
  p->cmd("prepareCalc");

  plumedNeedsEnergy = 0;
  p->cmd("isEnergyNeeded", &plumedNeedsEnergy);
  if (path_integral_mode == PATH_INTEGRAL_BEAD_MEAN ||
      path_integral_mode == PATH_INTEGRAL_BEAD_DENSITY) {
    int energy_requests = plumedNeedsEnergy;
    MPI_Allreduce(MPI_IN_PLACE, &energy_requests, 1, MPI_INT, MPI_SUM, universe->uworld);
    if (energy_requests)
      error->universe_all(
          FLERR, "Fix plumed path_integral bead modes do not support energy-dependent actions");
  }

  // Pass potential energy and virial if needed
  double *virial_lmp;
  if (plumedNeedsEnergy) {
    // Error if tail corrections are included
    if (force->pair && force->pair->tail_flag && comm->me == 0)
      error->warning(FLERR,
                     "Tail corrections to the pair potential included."
                     " The energy cannot be biased correctly in this case."
                     " Remove the tail corrections by removing the"
                     " command: pair_modify tail yes");

    // compute the potential energy
    double pot_energy = 0.;
    c_pe->compute_scalar();
    pot_energy = c_pe->scalar;

    // Divide energy by number of processes
    // Plumed wants it this way
    pot_energy /= comm->nprocs;
    p->cmd("setEnergy", &pot_energy);

    // Compute pressure due to the virial (no kinetic energy term!)
    c_press->compute_vector();
    virial_lmp = c_press->vector;

    // Check if pressure is finite
    if (!std::isfinite(virial_lmp[0]) || !std::isfinite(virial_lmp[1]) ||
        !std::isfinite(virial_lmp[2]) || !std::isfinite(virial_lmp[3]) ||
        !std::isfinite(virial_lmp[4]) || !std::isfinite(virial_lmp[5]))
      error->all(FLERR, "Non-numeric virial - Plumed cannot work with that");

    // Convert pressure to virial per number of MPI processes
    // From now on all virials are divided by the number of MPI processes

    double nktv2p = force->nktv2p;
    double inv_volume;
    if (domain->dimension == 3) {
      inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    } else {
      inv_volume = 1.0 / (domain->xprd * domain->yprd);
    }
    for (int i = 0; i < 6; i++) virial_lmp[i] /= (inv_volume * nktv2p * comm->nprocs);
    // Convert virial from lammps to plumed representation
    plmd_virial[0][0] = -virial_lmp[0];
    plmd_virial[1][1] = -virial_lmp[1];
    plmd_virial[2][2] = -virial_lmp[2];
    plmd_virial[0][1] = -virial_lmp[3];
    plmd_virial[0][2] = -virial_lmp[4];
    plmd_virial[1][2] = -virial_lmp[5];
  }
  // do the real calculation:
  p->cmd("performCalc");
  trace_nonfinite_state("post-perform-pre-scale", forces_before_plumed, false);

  if (path_integral_mode == PATH_INTEGRAL_BEAD_MEAN ||
      path_integral_mode == PATH_INTEGRAL_BEAD_DENSITY) {
    p->cmd("getBias", &bias);
    MPI_Allreduce(MPI_IN_PLACE, &plumedStopCondition, 1, MPI_INT, MPI_MAX, universe->uworld);
    const double dynamical_scale = universe->nworlds *
        (path_integral_mode == PATH_INTEGRAL_BEAD_DENSITY ? bead_density_force_scale : 1.0);
    for (int i = 0; i < nlocal; i++)
      for (int d = 0; d < 3; d++)
        atom->f[i][d] = forces_before_plumed[3 * i + d] +
            dynamical_scale * (atom->f[i][d] - forces_before_plumed[3 * i + d]);
    if (path_integral_mode == PATH_INTEGRAL_BEAD_DENSITY) {
      double bead_density_bias = comm->me == 0 ? bias : 0.0;
      MPI_Allreduce(MPI_IN_PLACE, &bead_density_bias, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
      bias = universe->iworld == 0 ? bead_density_force_scale * bead_density_bias : 0.0;
    } else if (universe->iworld != 0) {
      bias = 0.0;
    }
  }
  trace_nonfinite_state("post-scale", forces_before_plumed, false);
  if (plumedStopCondition) timer->force_timeout();

  // retransform virial to lammps representation and assign it to this
  // fix's virial. If the energy is biased, Plumed is giving back the full
  // virial and therefore we have to subtract the initial virial i.e. virial_lmp.
  // The vector virial contains only the contribution added by plumed.
  // The calculation of the pressure will be done by a compute pressure
  // and will include this contribution.
  if (plumedNeedsEnergy) {
    virial[0] = -plmd_virial[0][0] - virial_lmp[0];
    virial[1] = -plmd_virial[1][1] - virial_lmp[1];
    virial[2] = -plmd_virial[2][2] - virial_lmp[2];
    virial[3] = -plmd_virial[0][1] - virial_lmp[3];
    virial[4] = -plmd_virial[0][2] - virial_lmp[4];
    virial[5] = -plmd_virial[1][2] - virial_lmp[5];
  } else {
    virial[0] = -plmd_virial[0][0];
    virial[1] = -plmd_virial[1][1];
    virial[2] = -plmd_virial[2][2];
    virial[3] = -plmd_virial[0][1];
    virial[4] = -plmd_virial[0][2];
    virial[5] = -plmd_virial[1][2];
  }
  if (path_integral_mode == PATH_INTEGRAL_BEAD_DENSITY)
    for (int i = 0; i < 6; i++) virial[i] *= bead_density_force_scale * universe->nworlds;
  else if (path_integral_mode == PATH_INTEGRAL_BEAD_MEAN)
    for (int i = 0; i < 6; i++) virial[i] *= universe->nworlds;
  if (bead_bias_virial_pending) *bead_bias_virial_pending = 1;

  // Ask for the computes in the next time step
  // such that the virial and energy are tallied.
  // This should be changed to something that triggers the
  // calculation only if plumed needs it.
  c_pe->addstep(update->ntimestep + 1);
  c_press->addstep(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixPlumed::post_force_centroid()
{
  int invalid_atom_count = natoms != int(atom->natoms);
  MPI_Allreduce(MPI_IN_PLACE, &invalid_atom_count, 1, MPI_INT, MPI_SUM, universe->uworld);
  if (invalid_atom_count)
    error->universe_all(
        FLERR, "Fix plumed path_integral centroid requires a fixed atom count on each bead");

  for (int i = 0; i < 6; i++) virial[i] = 0.0;
  bias = 0.0;
  int plumed_stop_condition = 0;
  int needs_energy = 0;
  double plmd_virial[3][3] = {};

  if (plumed_active) {
    update_atom_data();
    for (int i = 0; i < nlocal; i++) {
      const int index = atom->tag[i] - 1;
      for (int d = 0; d < 3; d++) {
        centroid_positions[3 * i + d] = centroid_coordinates[3 * index + d];
        centroid_forces[3 * i + d] = 0.0;
      }
    }

    double box[3][3] = {};
    box[0][0] = domain->h[0];
    box[1][1] = domain->h[1];
    box[2][2] = domain->h[2];
    box[2][1] = domain->h[3];
    box[2][0] = domain->h[4];
    box[1][0] = domain->h[5];

    if (update->ntimestep > MAXSMALLINT)
      error->all(FLERR, "Fix plumed can only handle up to 2.1 billion timesteps");
    int step = int(update->ntimestep);
    try {
      p->cmd("setStep", &step);
      p->cmd("setStopFlag", &plumed_stop_condition);
      p->cmd("setPositions", centroid_positions);
      p->cmd("setBox", &box[0][0]);
      p->cmd("setForces", centroid_forces);
      p->cmd("setMasses", masses);
      p->cmd("setCharges", charges);
      p->cmd("setVirial", &plmd_virial[0][0]);
      p->cmd("prepareCalc");
      p->cmd("isEnergyNeeded", &needs_energy);
    } catch (const std::exception &exception) {
      error->universe_one(FLERR, fmt::format("Could not prepare PLUMED: {}", exception.what()));
    }
  }
  trace_centroid_nonfinite("pre-perform-centroid");

  MPI_Bcast(&needs_energy, 1, MPI_INT, 0, universe->uworld);
  if (needs_energy)
    error->universe_all(
        FLERR, "Fix plumed path_integral centroid does not support energy-dependent actions");

  if (plumed_active) {
    try {
      p->cmd("performCalc");
      p->cmd("getBias", &bias);
    } catch (const std::exception &exception) {
      error->universe_one(FLERR, fmt::format("Could not execute PLUMED: {}", exception.what()));
    }
    for (int i = 0; i < 3 * natoms; i++) centroid_forces_all[i] = 0.0;
    for (int i = 0; i < nlocal; i++) {
      const int index = atom->tag[i] - 1;
      for (int d = 0; d < 3; d++) centroid_forces_all[3 * index + d] = centroid_forces[3 * i + d];
    }
    MPI_Allreduce(MPI_IN_PLACE, centroid_forces_all, 3 * natoms, MPI_DOUBLE, MPI_SUM, world);
    virial[0] = -plmd_virial[0][0];
    virial[1] = -plmd_virial[1][1];
    virial[2] = -plmd_virial[2][2];
    virial[3] = -plmd_virial[0][1];
    virial[4] = -plmd_virial[0][2];
    virial[5] = -plmd_virial[1][2];
  }
  trace_centroid_nonfinite("post-perform-centroid");

  for (int i = 0; i < 6; i++) virial[i] *= universe->nworlds;
  if (centroid_virial_pending) *centroid_virial_pending = 1;

  MPI_Bcast(centroid_forces_all, 3 * natoms, MPI_DOUBLE, 0, universe->uworld);
  for (int i = 0; i < atom->nlocal; i++) {
    const int index = atom->tag[i] - 1;
    for (int d = 0; d < 3; d++)
      atom->f[i][d] += centroid_forces_all[3 * index + d] * centroid_force_scale;
  }
  trace_nonfinite_state("post-scale-centroid", nullptr, false);

  MPI_Bcast(&plumed_stop_condition, 1, MPI_INT, 0, universe->uworld);
  if (plumed_stop_condition) timer->force_timeout();
}

void FixPlumed::post_force_respa(int vflag, int ilevel, int /* iloop */)
{
  if (ilevel == nlevels_respa - 1) post_force(vflag);
}

void FixPlumed::min_post_force(int vflag)
{ post_force(vflag); }

void FixPlumed::reset_dt()
{ error->all(FLERR, "Fix plumed is not compatible with changing the timestep"); }

double FixPlumed::compute_scalar()
{ return bias; }

void *FixPlumed::extract(const char *name, int &dim)
{
  dim = 0;
  if (strcmp(name, "pimd_physical_bias_energy") == 0 && path_integral_mode != PATH_INTEGRAL_OFF)
    return &bias;
  return nullptr;
}

int FixPlumed::modify_param(int narg, char **arg)
{
  if (path_integral_mode == PATH_INTEGRAL_CENTROID &&
      (strcmp(arg[0], "pe") == 0 || strcmp(arg[0], "press") == 0))
    error->all(FLERR,
               "Fix_modify pe and press are not supported with fix plumed path_integral centroid");

  if (strcmp(arg[0], "pe") == 0) {
    if (narg < 2) error->all(FLERR, "Fix_modify pe requires an argument");
    modify->delete_compute(id_pe);
    delete[] id_pe;
    id_pe = utils::strdup(arg[1]);

    c_pe = modify->get_compute_by_id(id_pe);
    if (!c_pe) error->all(FLERR, "Could not find fix_modify potential energy ID {}", id_pe);

    if (c_pe->peflag == 0)
      error->all(FLERR, "Fix_modify compute pe ID {} does not compute potential energy", id_pe);
    if (c_pe->igroup != 0 && comm->me == 0)
      error->warning(FLERR, "Potential energy compute {} for fix PLUMED is not for group all",
                     id_pe);

    return 2;

  } else if (strcmp(arg[0], "press") == 0) {
    if (narg < 2) error->all(FLERR, "Fix_modify press requires an argument");
    modify->delete_compute(id_press);
    delete[] id_press;
    id_press = utils::strdup(arg[1]);

    c_press = modify->get_compute_by_id(id_press);
    if (!c_press) error->all(FLERR, "Could not find fix_modify compute pressure ID {}", id_press);

    if (c_press->pressflag == 0)
      error->all(FLERR, "Fix_modify compute pressure ID {} does not compute pressure", id_press);
    if (c_press->igroup != 0 && comm->me == 0)
      error->warning(FLERR, "Virial for fix PLUMED is not for group all");

    return 2;
  }
  return 0;
}

double FixPlumed::memory_usage()
{
  double bytes = double((8 + 8 + 4) * atom->nlocal);
  if (path_integral_mode == PATH_INTEGRAL_BEAD_DENSITY ||
      path_integral_mode == PATH_INTEGRAL_BEAD_MEAN)
    bytes += double(3 * sizeof(double) * atom->nlocal);
  if (path_integral_mode == PATH_INTEGRAL_CENTROID) {
    bytes += double(3 * sizeof(double) * natoms);
    if (plumed_active) bytes += double(6 * sizeof(double) * natoms);
  }
  return bytes;
}
