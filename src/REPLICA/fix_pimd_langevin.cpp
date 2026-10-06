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
   Package      FixPIMDLangevin
   Purpose      Path Integral Molecular Dynamics with Langevin Thermostat

   Yifan Li @ Princeton University (yifanl0716@gmail.com)
   Current Features:
   - Multi-processor parallelism for each bead
   - White-noise Langevin thermostat
   - Bussi-Zykova-Parrinello barostat (isotropic and anisotropic)
   - Several quantum estimators
   Futher plans:
   - Triclinic barostat
------------------------------------------------------------------------- */

#include "fix_pimd_langevin.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "pimd_partition.h"
#include "random_mars.h"
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <map>
#include <string>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_2PI;
using MathConst::MY_PI;
using MathConst::MY_SQRT2;
using MathConst::THIRD;
using MathSpecial::powint;

namespace {
std::map<int, std::string> Barostats{{FixPIMDLangevin::MTTK, "MTTK"},
                                     {FixPIMDLangevin::BZP, "BZP"}};
std::map<int, std::string> Ensembles{{FixPIMDLangevin::NVE, "NVE"},
                                     {FixPIMDLangevin::NVT, "NVT"},
                                     {FixPIMDLangevin::NPH, "NPH"},
                                     {FixPIMDLangevin::NPT, "NPT"}};
}    // namespace

namespace {
enum NonfiniteTraceField { TRACE_NONE, TRACE_BOX, TRACE_POSITION, TRACE_VELOCITY, TRACE_FORCE };

struct NonfiniteTraceRecord {
  int found;
  int reason;
  int field;
  int component;
  int universe_rank;
  int bead_world;
  int world_rank;
  int previous_available;
  tagint tag;
  double value;
  double current[9];
  double previous[9];
  double box[9];
};

const char *trace_field_name(int field)
{
  if (field == TRACE_BOX) return "box";
  if (field == TRACE_POSITION) return "position";
  if (field == TRACE_VELOCITY) return "velocity";
  if (field == TRACE_FORCE) return "force";
  return "unknown";
}
}    // namespace

namespace {
constexpr int TAG_INTER_REPLICA_COUNT = 10;
constexpr int TAG_INTER_REPLICA_TAGS = 11;
constexpr int TAG_INTER_REPLICA_VALS = 12;
}    // namespace

/* ---------------------------------------------------------------------- */

FixPIMDLangevin::FixPIMDLangevin(LAMMPS *lmp, int narg, char **arg) :
    FixPIMDLangevin(lmp, narg, arg, false)
{
}

/* ---------------------------------------------------------------------- */

FixPIMDLangevin::FixPIMDLangevin(LAMMPS *lmp, int narg, char **arg, bool allow_esynch) :
    Fix(lmp, narg, arg), mass(nullptr), plansend(nullptr), planrecv(nullptr), tagsend(nullptr),
    tagrecv(nullptr), bufsend(nullptr), bufrecv(nullptr), bufbeads(nullptr), bufsorted(nullptr),
    bufsortedall(nullptr), counts(nullptr), displacements(nullptr), lam(nullptr), M_x2xp(nullptr),
    M_xp2x(nullptr), M_f2fp(nullptr), M_fp2f(nullptr), modeindex(nullptr), tau_k(nullptr),
    c1_k(nullptr), c2_k(nullptr), _omega_k(nullptr), Lan_s(nullptr), Lan_c(nullptr),
    random(nullptr), xc(nullptr), xcall(nullptr), x_unwrap(nullptr), id_pe(nullptr),
    id_press(nullptr), c_pe(nullptr), c_press(nullptr), nonfinite_trace_prefix(nullptr),
    nonfinite_trace_nmax(0), nonfinite_trace_last_x(nullptr), nonfinite_trace_last_v(nullptr),
    nonfinite_trace_last_f(nullptr), nonfinite_trace_last_tag(nullptr)
{
  nonfinite_trace_last_stage[0] = '\0';
  restart_global = 1;
  time_integrate = 1;

  ntotal = 0;
  maxlocal = maxsend = maxunwrap = maxxc = 0;

  sizeplan = 0;

  method = NMPIMD;
  ensemble = NVT;
  integrator = OBABO;
  thermostat = PILE_L;
  barostat = BZP;
  fmass = 1.0;
  np = universe->nworlds;
  inverse_np = 1.0 / np;
  normal_mode_centroid_force_scale = universe->iworld == 0 ? 1.0 / sqrt((double) np) : 0.0;
  sp = 1.0;
  temp = 298.15;
  Lan_temp = 298.15;
  tau = 1.0;
  tau_p = 1.0;
  Pext = 1.0;
  pdim = 0;
  pilescale = 1.0;
  tstat_flag = 1;
  pstat_flag = 0;
  centroid_bias_virial_pending = 0;
  defer_normal_mode_force = 0;
  normal_mode_force_pending = 0;
  bead_bias_virial_pending = 0;
  mapflag = 1;
  removecomflag = 1;
  fmmode = PHYSICAL;
  pstyle = ISO;
  pote = tote = totke = totenthalpy = total_spring_energy = 0.0;
  centroid_vir = vir = vir_ = 0.0;
  ke_bead = se_bead = pe_bead = tote = t_prim = t_vir = t_cv = p_prim = p_md = p_cv = 0.0;

  int seed = -1;
  int scale_flag = 0;

  if (domain->dimension != 3) error->all(FLERR, fmt::format("Fix {} requires a 3d system", style));
  if (!atom->mass || atom->rmass_flag)
    error->all(FLERR, fmt::format("Fix {} requires per-type atom masses", style));

  for (int i = 1; i < universe->nworlds; i++)
    if (universe->procs_per_world[i] != universe->procs_per_world[0])
      error->all(
          FLERR,
          fmt::format("Fix {} requires the same number of processors in every partition", style));

  for (int i = 0; i < 6; i++) {
    vw[i] = 0.0;
    p_flag[i] = 0;
    p_target[i] = 0.0;
  }

  for (int i = 3; i < narg; i += 2) {
    if (strcmp(arg[i], "method") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} method", style), error);
      if (strcmp(arg[i + 1], "nmpimd") == 0)
        method = NMPIMD;
      else if (strcmp(arg[i + 1], "pimd") == 0) {
        method = PIMD;
        if (scale_flag)
          error->all(FLERR,
                     "Scale parameter of the PILE_L thermostat is not supported with method pimd");
      } else
        error->universe_all(FLERR, fmt::format("Unknown method parameter for fix {}", style));
    } else if (strcmp(arg[i], "integrator") == 0) {
      if (i + 2 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} integrator", style), error);
      if (strcmp(arg[i + 1], "obabo") == 0)
        integrator = OBABO;
      else if (strcmp(arg[i + 1], "baoab") == 0)
        integrator = BAOAB;
      else
        error->universe_all(FLERR,
                            fmt::format("Unknown integrator parameter for fix {}. Only obabo and "
                                        "baoab integrators are supported!",
                                        style));
    } else if (strcmp(arg[i], "ensemble") == 0) {
      if (i + 2 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} ensemble", style), error);
      if (strcmp(arg[i + 1], "nve") == 0) {
        ensemble = NVE;
        tstat_flag = 0;
        pstat_flag = 0;
      } else if (strcmp(arg[i + 1], "nvt") == 0) {
        ensemble = NVT;
        tstat_flag = 1;
        pstat_flag = 0;
      } else if (strcmp(arg[i + 1], "nph") == 0) {
        ensemble = NPH;
        tstat_flag = 0;
        pstat_flag = 1;
      } else if (strcmp(arg[i + 1], "npt") == 0) {
        ensemble = NPT;
        tstat_flag = 1;
        pstat_flag = 1;
      } else
        error->universe_all(FLERR,
                            fmt::format("Unknown ensemble parameter for fix {}. Only nve, nvt, "
                                        "nph, and npt ensembles are supported!",
                                        style));
    } else if (strcmp(arg[i], "fmass") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} fmass", style), error);
      fmass = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (fmass <= 0.0 || fmass > np)
        error->all(FLERR, fmt::format("Invalid fmass value for fix {}", style));
    } else if (strcmp(arg[i], "sp") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} sp", style), error);
      sp = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (sp <= 0.0) error->all(FLERR, fmt::format("Invalid sp value for fix {}", style));
    } else if (strcmp(arg[i], "fmmode") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} fmmode", style), error);
      if (strcmp(arg[i + 1], "physical") == 0)
        fmmode = PHYSICAL;
      else if (strcmp(arg[i + 1], "normal") == 0)
        fmmode = NORMAL;
      else
        error->universe_all(FLERR,
                            fmt::format("Unknown fictitious mass mode for fix {}. Only physical "
                                        "mass and normal mode mass are supported!",
                                        style));
    } else if (strcmp(arg[i], "scale") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} scale", style), error);
      if (method == PIMD)
        error->all(FLERR,
                   "Scale parameter of the PILE_L thermostat is not supported with method pimd");
      scale_flag = 1;
      pilescale = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (pilescale < 0.0)
        error->universe_all(FLERR, fmt::format("Invalid PILE_L scale value for fix {}", style));
    } else if (strcmp(arg[i], "temp") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} temp", style), error);
      temp = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (temp <= 0.0) error->all(FLERR, fmt::format("Invalid temp value for fix {}", style));
    } else if (strcmp(arg[i], "thermostat") == 0) {
      if (i + 3 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} thermostat", style), error);
      if (strcmp(arg[i + 1], "PILE_L") == 0) {
        thermostat = PILE_L;
        seed = utils::inumeric(FLERR, arg[i + 2], false, lmp);
        if (seed <= 0)
          error->all(FLERR, fmt::format("Invalid thermostat seed value for fix {}", style));
        i++;
      } else
        error->all(FLERR, fmt::format("Unknown thermostat parameter for fix {}", style));
    } else if (strcmp(arg[i], "tau") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} tau", style), error);
      tau = utils::numeric(FLERR, arg[i + 1], false, lmp);
    } else if (strcmp(arg[i], "barostat") == 0) {
      if (i + 2 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} barostat", style), error);
      if (strcmp(arg[i + 1], "MTTK") == 0) {
        error->all(FLERR, fmt::format("MTTK barostat is not supported by fix {}", style));
      } else if (strcmp(arg[i + 1], "BZP") == 0) {
        barostat = BZP;
      } else
        error->universe_all(FLERR, fmt::format("Unknown barostat parameter for fix {}", style));
    } else if (strcmp(arg[i], "iso") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} iso", style), error);
      pstyle = ISO;
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      Pext = utils::numeric(FLERR, arg[i + 1], false, lmp);
      p_target[0] = p_target[1] = p_target[2] = Pext;
    } else if (strcmp(arg[i], "aniso") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} aniso", style), error);
      pstyle = ANISO;
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      Pext = utils::numeric(FLERR, arg[i + 1], false, lmp);
      p_target[0] = p_target[1] = p_target[2] = Pext;
    } else if (strcmp(arg[i], "x") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} x", style), error);
      pstyle = ANISO;
      p_flag[0] = 1;
      p_target[0] = utils::numeric(FLERR, arg[i + 1], false, lmp);
    } else if (strcmp(arg[i], "y") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} y", style), error);
      pstyle = ANISO;
      p_flag[1] = 1;
      p_target[1] = utils::numeric(FLERR, arg[i + 1], false, lmp);
    } else if (strcmp(arg[i], "z") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} z", style), error);
      pstyle = ANISO;
      p_flag[2] = 1;
      p_target[2] = utils::numeric(FLERR, arg[i + 1], false, lmp);
    } else if (strcmp(arg[i], "taup") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} taup", style), error);
      tau_p = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (tau_p <= 0.0)
        error->universe_all(FLERR, fmt::format("Invalid tau_p value for fix {}", style));
    } else if (strcmp(arg[i], "fixcom") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} fixcom", style), error);
      if (strcmp(arg[i + 1], "yes") == 0)
        removecomflag = 1;
      else if (strcmp(arg[i + 1], "no") == 0)
        removecomflag = 0;
      else
        error->all(FLERR, fmt::format("Unknown fixcom value {} for fix {}", arg[i + 1], style));
    } else if (strcmp(arg[i], "nonfinite_trace") == 0) {
      if (i + 2 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} nonfinite_trace", style), error);
      delete[] nonfinite_trace_prefix;
      nonfinite_trace_prefix = utils::strdup(arg[i + 1]);
    } else if (strcmp(arg[i], "esynch") == 0) {
      if (!allow_esynch)
        error->all(FLERR, fmt::format("Unknown keyword {} for fix {}", arg[i], style));
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} esynch", style), error);
    } else if (strcmp(arg[i], "") != 0) {
      error->universe_all(FLERR, fmt::format("Unknown keyword {} for fix {}", arg[i], style));
    }
  }

  if (tstat_flag && seed < 0)
    error->all(FLERR,
               fmt::format("Thermostat seed must be specified for fix {} with {} ensemble", style,
                           Ensembles[ensemble]));

  pdim = p_flag[0] + p_flag[1] + p_flag[2];

  if (pstat_flag && !pdim)
    error->universe_all(
        FLERR, fmt::format("Must use pressure coupling with {} ensemble", Ensembles[ensemble]));
  if (!pstat_flag && pdim)
    error->universe_all(
        FLERR, fmt::format("Must not use pressure coupling with {} ensemble", Ensembles[ensemble]));

  if (method == PIMD && pstat_flag)
    error->universe_all(FLERR,
                        "Pressure control has not been supported for method pimd yet. Please set "
                        "method to nmpimd.");

  if (method == PIMD && fmmode == NORMAL)
    error->universe_all(
        FLERR, "Normal mode mass is not supported for method pimd. Please set method to nmpimd.");

  /* Initiation */

  global_freq = 1;
  vector_flag = 1;
  if (!pstat_flag) {
    size_vector = 10;
  } else if (pstat_flag) {
    if (pstyle == ISO) {
      size_vector = 15;
    } else if (pstyle == ANISO) {
      size_vector = 17;
    }
  }
  extvector = -1;
  extlist = new int[size_vector];
  for (int i = 0; i < size_vector; i++) extlist[i] = 1;
  for (int i = 7; i < 10; i++) extlist[i] = 0;
  if (pstat_flag) {
    extlist[10] = 0;
    if (pstyle == ANISO) extlist[11] = extlist[12] = 0;
  }
  kt = force->boltz * temp;
  if (pstat_flag) FixPIMDLangevin::baro_init();

  // some initilizations

  id_pe = utils::strdup(std::string(id) + "_pimd_pe");
  modify->add_compute(std::string(id_pe) + " all pe");

  id_press = utils::strdup(std::string(id) + "_pimd_press");
  modify->add_compute(std::string(id_press) + " all pressure thermo_temp virial");

  vol0 = domain->xprd * domain->yprd * domain->zprd;

  fixedpoint[0] = 0.5 * (domain->boxlo[0] + domain->boxhi[0]);
  fixedpoint[1] = 0.5 * (domain->boxlo[1] + domain->boxhi[1]);
  fixedpoint[2] = 0.5 * (domain->boxlo[2] + domain->boxhi[2]);
  if (pstat_flag) p_hydro = (p_target[0] + p_target[1] + p_target[2]) / pdim;

  // initialize Marsaglia RNG with processor-unique seed

  if (tstat_flag) {
    if (integrator == BAOAB || integrator == OBABO) {
      Lan_temp = temp;
      random = new RanMars(lmp, seed + universe->me);
    }
  }

  me = comm->me;
  nprocs = comm->nprocs;
  if (nprocs == 1)
    cmode = SINGLE_PROC;
  else
    cmode = MULTI_PROC;

  nprocs_universe = universe->nprocs;
  nreplica = universe->nworlds;
  ireplica = universe->iworld;

  if (nreplica == 1)
    mapflag = 0;
  else
    mapflag = 1;

  int *iroots = new int[nreplica];
  MPI_Group uworldgroup, rootgroup;

  for (int i = 0; i < nreplica; i++) iroots[i] = universe->root_proc[i];
  MPI_Comm_group(universe->uworld, &uworldgroup);
  MPI_Group_incl(uworldgroup, nreplica, iroots, &rootgroup);
  MPI_Comm_create(universe->uworld, rootgroup, &rootworld);
  if (rootgroup != MPI_GROUP_NULL) MPI_Group_free(&rootgroup);
  if (uworldgroup != MPI_GROUP_NULL) MPI_Group_free(&uworldgroup);
  delete[] iroots;

  ntotal = atom->natoms;
  if (atom->nmax > maxlocal) reallocate();
  if (atom->nmax > maxunwrap) reallocate_x_unwrap();
  if (atom->nmax > maxxc) reallocate_xc();
  memory->create(xcall, ntotal * 3, "FixPIMDLangevin:xcall");
}

/* ---------------------------------------------------------------------- */

FixPIMDLangevin::~FixPIMDLangevin()
{
  modify->delete_compute(id_pe);
  modify->delete_compute(id_press);
  delete[] id_pe;
  delete[] nonfinite_trace_prefix;
  memory->destroy(nonfinite_trace_last_x);
  memory->destroy(nonfinite_trace_last_v);
  memory->destroy(nonfinite_trace_last_f);
  memory->destroy(nonfinite_trace_last_tag);
  delete[] id_press;
  delete[] extlist;
  delete random;
  delete[] mass;
  delete[] _omega_k;
  delete[] Lan_c;
  delete[] Lan_s;
  delete[] tau_k;
  delete[] c1_k;
  delete[] c2_k;
  delete[] plansend;
  delete[] planrecv;
  delete[] modeindex;
  memory->destroy(xcall);
  if (cmode == SINGLE_PROC) {
    memory->destroy(bufsorted);
    memory->destroy(bufsortedall);
    memory->destroy(counts);
    memory->destroy(displacements);
  }

  memory->sfree(lam);
  memory->destroy(M_x2xp);
  memory->destroy(M_xp2x);
  memory->destroy(xc);
  memory->destroy(x_unwrap);
  memory->destroy(bufsend);
  memory->destroy(bufrecv);
  memory->destroy(tagsend);
  memory->destroy(tagrecv);
  memory->destroy(bufbeads);
  if (rootworld != MPI_COMM_NULL) MPI_Comm_free(&rootworld);
}

/* ---------------------------------------------------------------------- */

int FixPIMDLangevin::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::trace_nonfinite_state(const char *stage, const char *basis)
{
  if (!nonfinite_trace_prefix) return;

  NonfiniteTraceRecord local{};
  local.universe_rank = universe->me;
  local.bead_world = universe->iworld;
  local.world_rank = comm->me;
  local.tag = -1;

  const double box[9] = {domain->boxlo[0], domain->boxlo[1], domain->boxlo[2],
                         domain->boxhi[0], domain->boxhi[1], domain->boxhi[2],
                         domain->xy,       domain->xz,       domain->yz};
  for (int i = 0; i < 9; ++i) local.box[i] = box[i];

  for (int i = 0; i < 9; ++i) {
    if (!std::isfinite(box[i])) {
      local.found = 1;
      local.reason = 1;
      local.field = TRACE_BOX;
      local.component = i;
      local.value = box[i];
      break;
    }
  }
  if (!local.found && (!(domain->xprd > 0.0) || !(domain->yprd > 0.0) || !(domain->zprd > 0.0))) {
    local.found = 1;
    local.reason = 2;
    local.field = TRACE_BOX;
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
    if (!local.found || local.field == TRACE_BOX || tag < local.tag ||
        (tag == local.tag &&
         (field < local.field || (field == local.field && component < local.component)))) {
      if (local.field == TRACE_BOX) return;
      local.found = 1;
      local.reason = 1;
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
        if (!std::isfinite(atom->x[i][d])) select_atom_field(i, TRACE_POSITION, d, atom->x[i][d]);
        if (!std::isfinite(atom->v[i][d])) select_atom_field(i, TRACE_VELOCITY, d, atom->v[i][d]);
        if (!std::isfinite(atom->f[i][d])) select_atom_field(i, TRACE_FORCE, d, atom->f[i][d]);
      }
    }
  }

  if (local_index >= 0) {
    for (int d = 0; d < 3; ++d) {
      local.current[d] = atom->x[local_index][d];
      local.current[3 + d] = atom->v[local_index][d];
      local.current[6 + d] = atom->f[local_index][d];
    }
    if (local_index < nonfinite_trace_nmax && nonfinite_trace_last_tag[local_index] == local.tag) {
      local.previous_available = 1;
      for (int d = 0; d < 3; ++d) {
        local.previous[d] = nonfinite_trace_last_x[local_index][d];
        local.previous[3 + d] = nonfinite_trace_last_v[local_index][d];
        local.previous[6 + d] = nonfinite_trace_last_f[local_index][d];
      }
    }
  }

  int any_nonfinite = local.found;
  MPI_Allreduce(MPI_IN_PLACE, &any_nonfinite, 1, MPI_INT, MPI_MAX, universe->uworld);
  if (any_nonfinite) {

    std::vector<NonfiniteTraceRecord> records(universe->nprocs);
    MPI_Allgather(&local, sizeof(NonfiniteTraceRecord), MPI_BYTE, records.data(),
                  sizeof(NonfiniteTraceRecord), MPI_BYTE, universe->uworld);

    const NonfiniteTraceRecord *winner = nullptr;
    for (const auto &record : records) {
      if (!record.found) continue;
      if (!winner) {
        winner = &record;
        continue;
      }
      const int record_box = record.field == TRACE_BOX;
      const int winner_box = winner->field == TRACE_BOX;
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

    if (winner) {
      const std::string path =
          fmt::format("{}.pimd.step{}.u{}.w{}.r{}.txt", nonfinite_trace_prefix, update->ntimestep,
                      winner->universe_rank, winner->bead_world, winner->world_rank);
      if (universe->me == winner->universe_rank) {
        const std::string report = fmt::format(
            "schema=pimd-nonfinite-trace-v1\nstep={}\nstage={}\nbasis={}\n"
            "universe_rank={}\nbead_world={}\nworld_rank={}\natom_tag={}\nfield={}\n"
            "component={}\nreason={}\nvalue={:.17g}\n"
            "current_x={:.17g} {:.17g} {:.17g}\ncurrent_v={:.17g} {:.17g} {:.17g}\n"
            "current_f={:.17g} {:.17g} {:.17g}\nprevious_available={}\n"
            "previous_finite_stage={}\nprevious_x={:.17g} {:.17g} {:.17g}\n"
            "previous_v={:.17g} {:.17g} {:.17g}\nprevious_f={:.17g} {:.17g} {:.17g}\n"
            "boxlo={:.17g} {:.17g} {:.17g}\nboxhi={:.17g} {:.17g} {:.17g}\n"
            "tilt_xy_xz_yz={:.17g} {:.17g} {:.17g}\n",
            update->ntimestep, stage, basis, winner->universe_rank, winner->bead_world,
            winner->world_rank, winner->tag, trace_field_name(winner->field), winner->component,
            winner->reason == 2 ? "degenerate-box" : "non-finite", winner->value,
            winner->current[0], winner->current[1], winner->current[2], winner->current[3],
            winner->current[4], winner->current[5], winner->current[6], winner->current[7],
            winner->current[8], winner->previous_available, nonfinite_trace_last_stage,
            winner->previous[0], winner->previous[1], winner->previous[2], winner->previous[3],
            winner->previous[4], winner->previous[5], winner->previous[6], winner->previous[7],
            winner->previous[8], winner->box[0], winner->box[1], winner->box[2], winner->box[3],
            winner->box[4], winner->box[5], winner->box[6], winner->box[7], winner->box[8]);
        if (FILE *file = std::fopen(path.c_str(), "w")) {
          std::fwrite(report.data(), 1, report.size(), file);
          std::fclose(file);
        }
      }
      error->all(FLERR,
                 fmt::format("Fix pimd/langevin nonfinite trace detected a {} at step {} stage {} "
                             "on bead {} atom {} component {}; diagnostic {}",
                             trace_field_name(winner->field), update->ntimestep, stage,
                             winner->bead_world, winner->tag, winner->component, path));
    }
  }
  if (atom->nmax > nonfinite_trace_nmax) {
    memory->destroy(nonfinite_trace_last_x);
    memory->destroy(nonfinite_trace_last_v);
    memory->destroy(nonfinite_trace_last_f);
    memory->destroy(nonfinite_trace_last_tag);
    nonfinite_trace_nmax = atom->nmax;
    memory->create(nonfinite_trace_last_x, nonfinite_trace_nmax, 3,
                   "FixPIMDLangevin:nonfinite_trace_last_x");
    memory->create(nonfinite_trace_last_v, nonfinite_trace_nmax, 3,
                   "FixPIMDLangevin:nonfinite_trace_last_v");
    memory->create(nonfinite_trace_last_f, nonfinite_trace_nmax, 3,
                   "FixPIMDLangevin:nonfinite_trace_last_f");
    memory->create(nonfinite_trace_last_tag, nonfinite_trace_nmax,
                   "FixPIMDLangevin:nonfinite_trace_last_tag");
    for (int i = 0; i < nonfinite_trace_nmax; ++i) nonfinite_trace_last_tag[i] = -1;
  }
  for (int i = 0; i < atom->nlocal; ++i) {
    nonfinite_trace_last_tag[i] = atom->tag[i];
    for (int d = 0; d < 3; ++d) {
      nonfinite_trace_last_x[i][d] = atom->x[i][d];
      nonfinite_trace_last_v[i][d] = atom->v[i][d];
      nonfinite_trace_last_f[i][d] = atom->f[i][d];
    }
  }
  std::snprintf(nonfinite_trace_last_stage, sizeof(nonfinite_trace_last_stage), "%s", stage);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::init()
{
  defer_normal_mode_force = 0;
  normal_mode_force_pending = 0;
  bead_bias_virial_pending = 0;

  bigint min_natoms, max_natoms;
  MPI_Allreduce(&atom->natoms, &min_natoms, 1, MPI_LMP_BIGINT, MPI_MIN, universe->uworld);
  MPI_Allreduce(&atom->natoms, &max_natoms, 1, MPI_LMP_BIGINT, MPI_MAX, universe->uworld);
  if (min_natoms != max_natoms)
    error->all(FLERR,
               fmt::format("Fix {} requires the same number of atoms in every partition", style));

  double cell[9] = {domain->boxlo[0], domain->boxlo[1], domain->boxlo[2],
                    domain->boxhi[0], domain->boxhi[1], domain->boxhi[2],
                    domain->xy,       domain->xz,       domain->yz};
  double min_cell[9], max_cell[9];
  MPI_Allreduce(cell, min_cell, 9, MPI_DOUBLE, MPI_MIN, universe->uworld);
  MPI_Allreduce(cell, max_cell, 9, MPI_DOUBLE, MPI_MAX, universe->uworld);
  int cell_mismatch = 0;
  for (int i = 0; i < 9; i++)
    if (min_cell[i] != max_cell[i]) cell_mismatch = 1;

  int cell_flags[8] = {domain->triclinic,      domain->triclinic_general, domain->boundary[0][0],
                       domain->boundary[0][1], domain->boundary[1][0],    domain->boundary[1][1],
                       domain->boundary[2][0], domain->boundary[2][1]};
  int min_flags[8], max_flags[8];
  MPI_Allreduce(cell_flags, min_flags, 8, MPI_INT, MPI_MIN, universe->uworld);
  MPI_Allreduce(cell_flags, max_flags, 8, MPI_INT, MPI_MAX, universe->uworld);
  for (int i = 0; i < 8; i++)
    if (min_flags[i] != max_flags[i]) cell_mismatch = 1;
  if (cell_mismatch)
    error->all(FLERR,
               fmt::format("Fix {} requires the same simulation cell in every partition", style));

  int min_ntypes, max_ntypes;
  MPI_Allreduce(&atom->ntypes, &min_ntypes, 1, MPI_INT, MPI_MIN, universe->uworld);
  MPI_Allreduce(&atom->ntypes, &max_ntypes, 1, MPI_INT, MPI_MAX, universe->uworld);
  if (min_ntypes != max_ntypes)
    error->all(FLERR,
               fmt::format("Fix {} requires the same atom masses in every partition", style));

  std::vector<double> min_mass(atom->ntypes), max_mass(atom->ntypes);
  MPI_Allreduce(atom->mass + 1, min_mass.data(), atom->ntypes, MPI_DOUBLE, MPI_MIN,
                universe->uworld);
  MPI_Allreduce(atom->mass + 1, max_mass.data(), atom->ntypes, MPI_DOUBLE, MPI_MAX,
                universe->uworld);
  for (int i = 0; i < atom->ntypes; i++)
    if (min_mass[i] != max_mass[i])
      error->all(FLERR,
                 fmt::format("Fix {} requires the same atom masses in every partition", style));

  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, fmt::format("Fix {} requires an atom map, see atom_modify", style));
  if (atom->tag_consecutive() == 0)
    error->all(FLERR, "Atom IDs must be consecutive for fix {}", style);
  PIMDUtils::check_atom_consistency(lmp, style, groupbit);

  if (universe->me == 0 && universe->uscreen)
    utils::print(universe->uscreen, "Fix {}: initializing Path-Integral ...\n", style);

  // prepare the constants

  masstotal = group->mass(igroup);

  double planck = sp * force->hplanck;
  hbar = planck / MY_2PI;
  beta = 1.0 / (force->boltz * temp);
  double _fbond = 1.0 * np * np / (beta * beta * hbar * hbar);

  omega_np = np / (hbar * beta) * sqrt(force->mvv2e);
  beta_np = 1.0 / force->boltz / temp * inverse_np;
  fbond = _fbond * force->mvv2e;

  if ((universe->me == 0) && (universe->uscreen))
    utils::print(universe->uscreen, "Fix {}: -P/(beta^2 * hbar^2) = {:20.7e} (kcal/mol/A^2)\n\n",
                 style, fbond);

  if (integrator == OBABO) {
    dtf = 0.5 * update->dt * force->ftm2v;
    dtv = 0.5 * update->dt;
    dtv2 = dtv * dtv;
    dtv3 = THIRD * dtv2 * dtv * force->ftm2v;
  } else if (integrator == BAOAB) {
    dtf = 0.5 * update->dt * force->ftm2v;
    dtv = 0.5 * update->dt;
    dtv2 = dtv * dtv;
    dtv3 = THIRD * dtv2 * dtv * force->ftm2v;
  } else {
    error->universe_all(FLERR, fmt::format("Unknown integrator parameter for fix {}", style));
  }

  comm_init();

  // init() runs once per run command, so release the array of the previous one

  delete[] mass;
  mass = new double[atom->ntypes + 1];

  nmpimd_init();

  langevin_init();

  c_pe = modify->get_compute_by_id(id_pe);
  if (!c_pe) {
    error->universe_all(
        FLERR,
        fmt::format("Potential energy compute ID {} for fix {} does not exist", id_pe, style));
  } else {
    if (c_pe->peflag == 0)
      error->universe_all(
          FLERR,
          fmt::format("Compute ID {} for fix {} does not compute potential energy", id_pe, style));
  }

  c_press = modify->get_compute_by_id(id_press);
  if (!c_press) {
    error->universe_all(
        FLERR, fmt::format("Could not find fix {} pressure compute ID {}", style, id_press));
  } else {
    if (c_press->pressflag == 0)
      error->universe_all(
          FLERR,
          fmt::format("Compute ID {} for fix {} does not compute pressure", id_press, style));
  }

  t_prim = t_vir = t_cv = p_prim = p_vir = p_cv = p_md = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::setup(int vflag)
{
  trace_nonfinite_state("setup-entry", "bead-x-v-physical-f");

  int nlocal = atom->nlocal;
  double **x = atom->x;
  imageint *image = atom->image;
  if (mapflag) {
    for (int i = 0; i < nlocal; i++) domain->unmap(x[i], image[i]);
  }

  if (method == NMPIMD) {
    inter_replica_comm(x);
    if (cmode == SINGLE_PROC)
      nmpimd_transform(bufsortedall, x, M_x2xp[universe->iworld]);
    else if (cmode == MULTI_PROC)
      nmpimd_transform(bufbeads, x, M_x2xp[universe->iworld]);
    trace_nonfinite_state("setup-bead-to-normal-post", "normal-mode");
  } else if (method == PIMD) {
    prepare_coordinates();
  } else {
    error->universe_all(
        FLERR,
        fmt::format("Unknown method parameter for fix {}. Only nmpimd and pimd are supported!",
                    style));
  }
  collect_xc();
  compute_spring_energy();
  compute_t_prim();
  compute_p_prim();
  if (method == NMPIMD) {
    inter_replica_comm(x);
    if (cmode == SINGLE_PROC)
      nmpimd_transform(bufsortedall, x, M_xp2x[universe->iworld]);
    else if (cmode == MULTI_PROC)
      nmpimd_transform(bufbeads, x, M_xp2x[universe->iworld]);
    trace_nonfinite_state("setup-normal-to-bead-post", "bead-x-normal-vf");
  }
  if (mapflag) {
    for (int i = 0; i < nlocal; i++) domain->unmap_inv(x[i], image[i]);
  }

  post_force(vflag);
  compute_totke();
  end_of_step();
  c_pe->addstep(update->ntimestep + 1);
  c_press->addstep(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::initial_integrate(int /*vflag*/)
{
  trace_nonfinite_state("initial-entry", "bead-x-normal-vf");
  prepare_normal_mode_forces();
  trace_nonfinite_state("initial-force-ready", "bead-x-normal-vf");

  int nlocal = atom->nlocal;
  double **x = atom->x;
  imageint *image = atom->image;
  if (mapflag) {
    for (int i = 0; i < nlocal; i++) domain->unmap(x[i], image[i]);
  }
  if (integrator == OBABO) {
    if (tstat_flag) {
      o_step();
      if (removecomflag) remove_com_motion();
      if (pstat_flag) press_o_step();
    }
    trace_nonfinite_state("initial-obabo-o1-post", "bead-x-normal-vf");
    if (pstat_flag) {
      compute_totke();
      compute_p_cv();
      press_v_step();
    }
    b_step();
    trace_nonfinite_state("initial-obabo-b-post", "bead-x-normal-vf");
    if (method == NMPIMD) {
      inter_replica_comm(x);
      if (cmode == SINGLE_PROC)
        nmpimd_transform(bufsortedall, x, M_x2xp[universe->iworld]);
      else if (cmode == MULTI_PROC)
        nmpimd_transform(bufbeads, x, M_x2xp[universe->iworld]);
      trace_nonfinite_state("initial-obabo-bead-to-normal-post", "normal-mode");
      qc_step();
      a_step();
      qc_step();
      a_step();
      trace_nonfinite_state("initial-obabo-a-post", "normal-mode");
    } else if (method == PIMD) {
      q_step();
      q_step();
    } else {
      error->universe_all(
          FLERR,
          fmt::format("Unknown method parameter for fix {}. Only nmpimd and pimd are supported!",
                      style));
    }
  } else if (integrator == BAOAB) {
    if (pstat_flag) {
      compute_totke();
      compute_p_cv();
      press_v_step();
    }
    b_step();
    trace_nonfinite_state("initial-baoab-b-post", "bead-x-normal-vf");
    if (method == NMPIMD) {
      inter_replica_comm(x);
      if (cmode == SINGLE_PROC)
        nmpimd_transform(bufsortedall, x, M_x2xp[universe->iworld]);
      else if (cmode == MULTI_PROC)
        nmpimd_transform(bufbeads, x, M_x2xp[universe->iworld]);
      trace_nonfinite_state("initial-baoab-bead-to-normal-post", "normal-mode");
      qc_step();
      a_step();
      trace_nonfinite_state("initial-baoab-a1-post", "normal-mode");
    } else if (method == PIMD) {
      q_step();
    } else {
      error->universe_all(
          FLERR,
          fmt::format("Unknown method parameter for fix {}. Only nmpimd and pimd are supported!",
                      style));
    }
    if (tstat_flag) {
      o_step();
      if (removecomflag) remove_com_motion();
      if (pstat_flag) press_o_step();
    }
    trace_nonfinite_state("initial-baoab-o-post", "normal-mode");
    if (method == NMPIMD) {
      qc_step();
      a_step();
      trace_nonfinite_state("initial-baoab-a2-post", "normal-mode");
    } else if (method == PIMD) {
      q_step();
      trace_nonfinite_state("initial-baoab-q2-post", "bead");
    } else {
      error->universe_all(
          FLERR,
          fmt::format("Unknown method parameter for fix {}. Only nmpimd and pimd are supported!",
                      style));
    }
  } else {
    error->universe_all(FLERR,
                        fmt::format("Unknown integrator parameter for fix {}. Only obabo and baoab "
                                    "integrators are supported!",
                                    style));
  }
  collect_xc();

  if (method == NMPIMD) {
    compute_spring_energy();
    compute_t_prim();
    compute_p_prim();
  }

  if (method == NMPIMD) {
    inter_replica_comm(x);
    if (cmode == SINGLE_PROC)
      nmpimd_transform(bufsortedall, x, M_xp2x[universe->iworld]);
    else if (cmode == MULTI_PROC)
      nmpimd_transform(bufbeads, x, M_xp2x[universe->iworld]);
    trace_nonfinite_state("initial-normal-to-bead-post", "bead-x-normal-vf");
  }

  if (mapflag) {
    for (int i = 0; i < nlocal; i++) { domain->unmap_inv(x[i], image[i]); }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::final_integrate()
{
  trace_nonfinite_state("final-entry", "bead-x-normal-v-physical-f");
  prepare_normal_mode_forces();
  trace_nonfinite_state("final-force-ready", "bead-x-normal-vf");

  if (pstat_flag) {
    compute_totke();
    compute_p_cv();
    press_v_step();
  }
  b_step();
  trace_nonfinite_state("final-b-post", "bead-x-normal-vf");
  if (integrator == OBABO) {
    if (tstat_flag) {
      o_step();
      if (removecomflag) remove_com_motion();
      if (pstat_flag) press_o_step();
    }
    trace_nonfinite_state("final-obabo-o-post", "bead-x-normal-vf");
  } else if (integrator == BAOAB) {

  } else {
    error->universe_all(FLERR, fmt::format("Unknown integrator parameter for fix {}", style));
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::prepare_coordinates()
{ inter_replica_comm(atom->x); }

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::post_force(int /*flag*/)
{
  trace_nonfinite_state("post-physical-force", "bead-x-normal-v-physical-f");

  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **f = atom->f;
  imageint *image = atom->image;
  tagint *tag = atom->tag;

  if (atom->nmax > maxunwrap) reallocate_x_unwrap();
  if (atom->nmax > maxxc) reallocate_xc();
  for (int i = 0; i < nlocal; i++) {
    x_unwrap[i][0] = x[i][0];
    x_unwrap[i][1] = x[i][1];
    x_unwrap[i][2] = x[i][2];
  }
  if (mapflag) {
    for (int i = 0; i < nlocal; i++) { domain->unmap(x_unwrap[i], image[i]); }
  }
  for (int i = 0; i < nlocal; i++) {
    xc[i][0] = xcall[3 * (tag[i] - 1) + 0];
    xc[i][1] = xcall[3 * (tag[i] - 1) + 1];
    xc[i][2] = xcall[3 * (tag[i] - 1) + 2];
  }

  compute_vir();
  compute_xf_vir();
  compute_cvir();
  compute_t_vir();

  if (method == PIMD) {
    if (mapflag) {
      for (int i = 0; i < nlocal; i++) { domain->unmap(x[i], image[i]); }
    }
    prepare_coordinates();
    spring_force();
    compute_spring_energy();
    compute_t_prim();
    compute_p_prim();
    if (mapflag) {
      for (int i = 0; i < nlocal; i++) { domain->unmap_inv(x[i], image[i]); }
    }
  }
  if (method == NMPIMD) {
    if (defer_normal_mode_force) {
      normal_mode_force_pending = 1;
    } else {
      inter_replica_comm(f);
      if (cmode == SINGLE_PROC)
        nmpimd_transform(bufsortedall, f, M_x2xp[universe->iworld]);
      else if (cmode == MULTI_PROC)
        nmpimd_transform(bufbeads, f, M_x2xp[universe->iworld]);
      trace_nonfinite_state("post-force-transform-post", "bead-x-normal-vf");
    }
  }

  c_pe->addstep(update->ntimestep + 1);
  c_press->addstep(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::end_of_step()
{
  if (bead_bias_virial_pending) prepare_normal_mode_forces();
  compute_pote();
  compute_totke();
  compute_p_cv();
  compute_tote();
  if (pstat_flag) compute_totenthalpy();
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::collect_xc()
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  tagint *tag = atom->tag;

  if (method == PIMD) {
    inter_replica_comm(x);
    if (ireplica == 0) {
      for (int i = 0; i < ntotal * 3; i++) xcall[i] = 0.0;
      if (cmode == SINGLE_PROC) {
        for (int i = 0; i < ntotal; i++) {
          for (int d = 0; d < 3; d++) {
            for (int bead = 0; bead < np; bead++)
              xcall[3 * i + d] += bufsortedall[bead * ntotal + i][d] * inverse_np;
          }
        }
      } else {
        for (int i = 0; i < nlocal; i++) {
          const int index = tag[i] - 1;
          for (int d = 0; d < 3; d++) {
            for (int bead = 0; bead < np; bead++)
              xcall[3 * index + d] += bufbeads[bead][3 * i + d] * inverse_np;
          }
        }
        MPI_Allreduce(MPI_IN_PLACE, xcall, ntotal * 3, MPI_DOUBLE, MPI_SUM, world);
      }
    }
  } else {
    if (ireplica == 0) {
      for (int i = 0; i < ntotal; i++) xcall[3 * i + 0] = xcall[3 * i + 1] = xcall[3 * i + 2] = 0.0;

      const double sqrtnp = sqrt((double) np);
      for (int i = 0; i < nlocal; i++) {
        xcall[3 * (tag[i] - 1) + 0] = x[i][0] / sqrtnp;
        xcall[3 * (tag[i] - 1) + 1] = x[i][1] / sqrtnp;
        xcall[3 * (tag[i] - 1) + 2] = x[i][2] / sqrtnp;
      }

      if (cmode == MULTI_PROC)
        MPI_Allreduce(MPI_IN_PLACE, xcall, ntotal * 3, MPI_DOUBLE, MPI_SUM, world);
    }
  }
  MPI_Bcast(xcall, ntotal * 3, MPI_DOUBLE, 0, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void *FixPIMDLangevin::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "centroid_coordinates") == 0) {
    dim = 1;
    return xcall;
  }
  if (strcmp(str, "nvt_unwrapped_coordinates") == 0 && ensemble == NVT) {
    // Refreshed in post_force; callers must not retain the reallocatable pointer.
    dim = 2;
    return x_unwrap;
  }
  if (strcmp(str, "nbeads") == 0) return &np;
  if (strcmp(str, "centroid_bias_force_scale") == 0 && method == PIMD && ensemble == NVT)
    return &inverse_np;
  if (strcmp(str, "bead_bias_force_scale") == 0 &&
      ((method == PIMD && ensemble == NVT) ||
       (method == NMPIMD && (ensemble == NVT || ensemble == NPH || ensemble == NPT))))
    return &inverse_np;
  if (strcmp(str, "normal_mode_centroid_force_scale") == 0 && method == NMPIMD &&
      (ensemble == NVT || ensemble == NPH || ensemble == NPT))
    return &normal_mode_centroid_force_scale;
  if (strcmp(str, "centroid_bias_virial_pending") == 0 && method == NMPIMD &&
      (ensemble == NVT || ensemble == NPH || ensemble == NPT))
    return &centroid_bias_virial_pending;
  if (strcmp(str, "defer_normal_mode_force") == 0 && method == NMPIMD &&
      (ensemble == NVT || ensemble == NPH || ensemble == NPT))
    return &defer_normal_mode_force;
  if (strcmp(str, "bead_bias_virial_pending") == 0 && method == NMPIMD &&
      (ensemble == NVT || ensemble == NPH || ensemble == NPT))
    return &bead_bias_virial_pending;
  return nullptr;
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::b_step()
{
  // used for both NMPIMD and PIMD
  // For NMPIMD, force only includes the contribution of external potential.
  // For PIMD, force includes the contributions of external potential and spring force.
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  double **v = atom->v;
  double **f = atom->f;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double dtfm = dtf / mass[type[i]];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::qc_step()
{
  // used for NMPIMD
  // evolve the centroid mode
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double oldlo, oldhi;

  if (!pstat_flag) {
    if (universe->iworld == 0) {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
        }
      }
    }
  } else {
    if (universe->iworld == 0) {
      double expp[3], expq[3], velocity_factor[3];
      if (pstyle == ISO) {
        vw[1] = vw[0];
        vw[2] = vw[0];
      }
      for (int j = 0; j < 3; j++) {
        const double eta = dtv * vw[j];
        expq[j] = exp(eta);
        expp[j] = exp(-eta);
        velocity_factor[j] = dtv;
        if (eta != 0.0) velocity_factor[j] *= sinh(eta) / eta;
      }
      if (barostat == BZP) {
        for (int i = 0; i < nlocal; i++) {
          if (mask[i] & groupbit) {
            for (int j = 0; j < 3; j++) {
              if (p_flag[j]) {
                x[i][j] = expq[j] * x[i][j] + velocity_factor[j] * v[i][j];
                v[i][j] = expp[j] * v[i][j];
              } else {
                x[i][j] += dtv * v[i][j];
              }
            }
          }
        }
        oldlo = domain->boxlo[0];
        oldhi = domain->boxhi[0];

        domain->boxlo[0] = (oldlo - fixedpoint[0]) * expq[0] + fixedpoint[0];
        domain->boxhi[0] = (oldhi - fixedpoint[0]) * expq[0] + fixedpoint[0];

        oldlo = domain->boxlo[1];
        oldhi = domain->boxhi[1];
        domain->boxlo[1] = (oldlo - fixedpoint[1]) * expq[1] + fixedpoint[1];
        domain->boxhi[1] = (oldhi - fixedpoint[1]) * expq[1] + fixedpoint[1];

        oldlo = domain->boxlo[2];
        oldhi = domain->boxhi[2];
        domain->boxlo[2] = (oldlo - fixedpoint[2]) * expq[2] + fixedpoint[2];
        domain->boxhi[2] = (oldhi - fixedpoint[2]) * expq[2] + fixedpoint[2];
      }
    }
    MPI_Barrier(universe->uworld);
    MPI_Bcast(&domain->boxlo[0], 3, MPI_DOUBLE, 0, universe->uworld);
    MPI_Bcast(&domain->boxhi[0], 3, MPI_DOUBLE, 0, universe->uworld);
    domain->set_global_box();
    domain->set_local_box();
    if (force->kspace) force->kspace->setup();
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::a_step()
{
  // used for NMPIMD
  // use analytical solution of harmonic oscillator to evolve the non-centroid modes
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double x0, x1, x2, v0, v1, v2;    // three components of x[i] and v[i]

  if (universe->iworld != 0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        x0 = x[i][0];
        x1 = x[i][1];
        x2 = x[i][2];
        v0 = v[i][0];
        v1 = v[i][1];
        v2 = v[i][2];
        x[i][0] = Lan_c[universe->iworld] * x0 +
            1.0 / _omega_k[universe->iworld] * Lan_s[universe->iworld] * v0;
        x[i][1] = Lan_c[universe->iworld] * x1 +
            1.0 / _omega_k[universe->iworld] * Lan_s[universe->iworld] * v1;
        x[i][2] = Lan_c[universe->iworld] * x2 +
            1.0 / _omega_k[universe->iworld] * Lan_s[universe->iworld] * v2;
        v[i][0] = -1.0 * _omega_k[universe->iworld] * Lan_s[universe->iworld] * x0 +
            Lan_c[universe->iworld] * v0;
        v[i][1] = -1.0 * _omega_k[universe->iworld] * Lan_s[universe->iworld] * x1 +
            Lan_c[universe->iworld] * v1;
        v[i][2] = -1.0 * _omega_k[universe->iworld] * Lan_s[universe->iworld] * x2 +
            Lan_c[universe->iworld] * v2;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::q_step()
{
  // used for PIMD
  // evolve all beads
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;

  if (!pstat_flag) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::baro_init()
{
  vw[0] = vw[1] = vw[2] = vw[3] = vw[4] = vw[5] = 0.0;
  if (pstyle == ISO) {
    W = 3 * (group->count(igroup)) * tau_p * tau_p * np * kt;
  }    // consistent with the definition in i-Pi
  else if (pstyle == ANISO) {
    W = group->count(igroup) * tau_p * tau_p * np * kt;
  }
  Vcoeff = 1.0;
  std::string out = fmt::format("\nInitializing PIMD {:s} barostat...\n", Barostats[barostat]);
  out += fmt::format("  The barostat mass is W = {:.16e}\n", W);
  if (universe->me == 0) utils::logmesg(lmp, out);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::press_v_step()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **f = atom->f;
  double **v = atom->v;
  int *type = atom->type;
  double volume = domain->xprd * domain->yprd * domain->zprd;

  if (pstyle == ISO) {
    if (barostat == BZP) {
      vw[0] += dtv * 3 * (volume * np * (p_cv - p_hydro) / force->nktv2p + Vcoeff / beta_np) / W;
      if (universe->iworld == 0) {
        double dvw_proc = 0.0, dvw = 0.0;
        for (int i = 0; i < nlocal; i++) {
          if (mask[i] & groupbit) {
            for (int j = 0; j < 3; j++) {
              dvw_proc +=
                  dtv2 * f[i][j] * v[i][j] / W + dtv3 * f[i][j] * f[i][j] / mass[type[i]] / W;
            }
          }
        }
        MPI_Allreduce(&dvw_proc, &dvw, 1, MPI_DOUBLE, MPI_SUM, world);
        vw[0] += dvw;
      }
      MPI_Barrier(universe->uworld);
      MPI_Bcast(&vw[0], 1, MPI_DOUBLE, 0, universe->uworld);
    } else if (barostat == MTTK) {
      double mtk_term1 = 2.0 / group->count(igroup) * totke / 3.0;
      vw[0] += 0.5 * dtv * (volume * np * (p_md - p_hydro) + mtk_term1) / W;
    }
  } else if (pstyle == ANISO) {
    compute_stress_tensor();
    for (int ii = 0; ii < 3; ii++) {
      if (p_flag[ii]) {
        vw[ii] += dtv *
            (volume * np * (stress_tensor[ii] - p_target[ii]) / force->nktv2p + Vcoeff / beta_np) /
            W;
        if (universe->iworld == 0) {
          double dvw_proc = 0.0, dvw = 0.0;
          for (int i = 0; i < nlocal; i++) {
            if (mask[i] & groupbit) {
              dvw_proc +=
                  dtv2 * f[i][ii] * v[i][ii] / W + dtv3 * f[i][ii] * f[i][ii] / mass[type[i]] / W;
            }
          }
          MPI_Allreduce(&dvw_proc, &dvw, 1, MPI_DOUBLE, MPI_SUM, world);
          vw[ii] += dvw;
        }
      }
    }
    MPI_Barrier(universe->uworld);
    MPI_Bcast(vw, 3, MPI_DOUBLE, 0, universe->uworld);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::press_o_step()
{
  if (pstyle == ISO) {
    if (universe->me == 0) vw[0] = c1 * vw[0] + c2 * sqrt(1.0 / W / beta_np) * random->gaussian();
    MPI_Barrier(universe->uworld);
    MPI_Bcast(&vw[0], 1, MPI_DOUBLE, 0, universe->uworld);
  } else if (pstyle == ANISO) {
    if (universe->me == 0) {
      for (int ii = 0; ii < 3; ii++) {
        if (p_flag[ii]) vw[ii] = c1 * vw[ii] + c2 * sqrt(1.0 / W / beta_np) * random->gaussian();
      }
    }
    MPI_Barrier(universe->uworld);
    MPI_Bcast(&vw, 3, MPI_DOUBLE, 0, universe->uworld);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::langevin_init()
{
  double beta = 1.0 / kt;
  const double _omega_np = np / beta / hbar;
  double _omega_np_dt_half = _omega_np * update->dt * 0.5;

  _omega_k = new double[np];
  Lan_c = new double[np];
  Lan_s = new double[np];
  if (method == NMPIMD) {
    if (fmmode == PHYSICAL) {
      for (int i = 0; i < np; i++) {
        _omega_k[i] = _omega_np * sqrt(lam[i]) / sqrt(fmass);
        Lan_c[i] = cos(sqrt(lam[i]) * _omega_np_dt_half);
        Lan_s[i] = sin(sqrt(lam[i]) * _omega_np_dt_half);
      }
    } else if (fmmode == NORMAL) {
      for (int i = 0; i < np; i++) {
        _omega_k[i] = _omega_np / sqrt(fmass);
        Lan_c[i] = cos(_omega_np_dt_half);
        Lan_s[i] = sin(_omega_np_dt_half);
      }
    } else {
      error->universe_all(FLERR, "Unknown fmmode setting; only physical and normal are supported!");
    }
  }

  if (tau > 0)
    gamma = 1.0 / tau;
  else
    gamma = np / beta / hbar;

  if (integrator == OBABO)
    c1 = exp(-gamma * 0.5 * update->dt);    // tau is the damping time of the centroid mode.
  else if (integrator == BAOAB)
    c1 = exp(-gamma * update->dt);
  else
    error->universe_all(FLERR,
                        fmt::format("Unknown integrator parameter for fix {}. Only obabo and baoab "
                                    "integrators are supported!",
                                    style));

  c2 = sqrt(1.0 - c1 * c1);    // note that c1 and c2 here only works for the centroid mode.

  if (thermostat == PILE_L) {
    std::string out = "Initializing PI Langevin equation thermostat...\n";
    out += "  Bead ID    |    omega    |    tau    |    c1    |    c2\n";
    if (method == NMPIMD) {
      tau_k = new double[np];
      c1_k = new double[np];
      c2_k = new double[np];
      tau_k[0] = 1.0 / gamma;
      c1_k[0] = c1;
      c2_k[0] = c2;
      for (int i = 1; i < np; i++) {
        if (pilescale == 0.0) {
          tau_k[i] = std::numeric_limits<double>::infinity();
          c1_k[i] = 1.0;
          c2_k[i] = 0.0;
          continue;
        }
        tau_k[i] = 0.5 / pilescale / _omega_k[i];
        if (integrator == OBABO)
          c1_k[i] = exp(-0.5 * update->dt / tau_k[i]);
        else if (integrator == BAOAB)
          c1_k[i] = exp(-1.0 * update->dt / tau_k[i]);
        else
          error->universe_all(FLERR,
                              fmt::format("Unknown integrator parameter for fix {}. Only obabo and "
                                          "baoab integrators are supported!",
                                          style));
        c2_k[i] = sqrt(1.0 - c1_k[i] * c1_k[i]);
      }
      for (int i = 0; i < np; i++) {
        out += fmt::format("      {:d}     {:.8e} {:.8e} {:.8e} {:.8e}\n", i, _omega_k[i], tau_k[i],
                           c1_k[i], c2_k[i]);
      }
    } else if (method == PIMD) {
      for (int i = 0; i < np; i++) {
        out += fmt::format("      {:d}     {:.8e} {:.8e} {:.8e} {:.8e}\n", i,
                           _omega_np / sqrt(fmass), tau, c1, c2);
      }
    }
    if (thermostat == PILE_L) out += "  PILE_L thermostat successfully initialized!\n";
    out += "\n";
    if (universe->me == 0) utils::logmesg(lmp, out);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::o_step()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  double beta_np = 1.0 / force->boltz / Lan_temp * inverse_np * force->mvv2e;
  if (thermostat == PILE_L) {
    if (method == NMPIMD) {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          atom->v[i][0] = c1_k[universe->iworld] * atom->v[i][0] +
              c2_k[universe->iworld] * sqrt(1.0 / mass[type[i]] / beta_np) * random->gaussian();
          atom->v[i][1] = c1_k[universe->iworld] * atom->v[i][1] +
              c2_k[universe->iworld] * sqrt(1.0 / mass[type[i]] / beta_np) * random->gaussian();
          atom->v[i][2] = c1_k[universe->iworld] * atom->v[i][2] +
              c2_k[universe->iworld] * sqrt(1.0 / mass[type[i]] / beta_np) * random->gaussian();
        }
      }
    } else if (method == PIMD) {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          atom->v[i][0] =
              c1 * atom->v[i][0] + c2 * sqrt(1.0 / mass[type[i]] / beta_np) * random->gaussian();
          atom->v[i][1] =
              c1 * atom->v[i][1] + c2 * sqrt(1.0 / mass[type[i]] / beta_np) * random->gaussian();
          atom->v[i][2] =
              c1 * atom->v[i][2] + c2 * sqrt(1.0 / mass[type[i]] / beta_np) * random->gaussian();
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Normal Mode PIMD
   ------------------------------------------------------------------------- */

void FixPIMDLangevin::nmpimd_init()
{
  // init() calls this once per run, so release what an earlier run allocated

  memory->destroy(M_x2xp);
  memory->destroy(M_xp2x);
  memory->sfree(lam);
  lam = nullptr;

  memory->create(M_x2xp, np, np, "fix_feynman:M_x2xp");
  memory->create(M_xp2x, np, np, "fix_feynman:M_xp2x");

  lam = (double *) memory->smalloc(sizeof(double) * np, "FixPIMDLangevin::lam");

  // Set up  eigenvalues
  for (int i = 0; i < np; i++) {
    double sin_tmp = sin(i * MY_PI / np);
    lam[i] = 4 * sin_tmp * sin_tmp;
  }

  // Set up eigenvectors for degenerated modes
  const double sqrtnp = sqrt((double) np);
  for (int j = 0; j < np; j++) {
    for (int i = 1; i < int(np / 2) + 1; i++) {
      M_x2xp[i][j] = MY_SQRT2 * cos(MY_2PI * double(i) * double(j) / double(np)) / sqrtnp;
    }
    for (int i = int(np / 2) + 1; i < np; i++) {
      M_x2xp[i][j] = MY_SQRT2 * sin(MY_2PI * double(i) * double(j) / double(np)) / sqrtnp;
    }
  }

  // Set up eigenvectors for non-degenerated modes
  for (int i = 0; i < np; i++) {
    M_x2xp[0][i] = 1.0 / sqrtnp;
    if (np % 2 == 0) M_x2xp[np / 2][i] = 1.0 / sqrtnp * powint(-1.0, i);
  }

  // Set up Ut
  for (int i = 0; i < np; i++)
    for (int j = 0; j < np; j++) { M_xp2x[i][j] = M_x2xp[j][i]; }

  // Set up fictitious masses
  int iworld = universe->iworld;
  for (int i = 1; i <= atom->ntypes; i++) {
    mass[i] = atom->mass[i];
    mass[i] *= fmass;
    if (iworld) {
      if (fmmode == PHYSICAL) {
        mass[i] *= 1.0;
      } else if (fmmode == NORMAL) {
        mass[i] *= lam[iworld];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::nmpimd_transform(double **src, double **des, double *vector)
{
  if (cmode == SINGLE_PROC) {
    for (int i = 0; i < ntotal; i++) {
      for (int d = 0; d < 3; d++) {
        bufsorted[i][d] = 0.0;
        for (int j = 0; j < nreplica; j++) {
          bufsorted[i][d] += src[j * ntotal + i][d] * vector[j];
        }
      }
    }
    for (int i = 0; i < ntotal; i++) {
      tagint tagtmp = atom->tag[i];
      for (int d = 0; d < 3; d++) { des[i][d] = bufsorted[tagtmp - 1][d]; }
    }
  } else if (cmode == MULTI_PROC) {
    int n = atom->nlocal;
    int m = 0;

    for (int i = 0; i < n; i++) {
      for (int d = 0; d < 3; d++) {
        des[i][d] = 0.0;
        for (int j = 0; j < np; j++) { des[i][d] += (src[j][m] * vector[j]); }
        m++;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::prepare_normal_mode_forces()
{
  if (method != NMPIMD || !normal_mode_force_pending) return;
  trace_nonfinite_state("deferred-force-transform-pre", "bead-x-normal-v-physical-f");

  if (bead_bias_virial_pending) {
    compute_vir();
    compute_xf_vir();
    compute_cvir();
    compute_t_vir();
    bead_bias_virial_pending = 0;
  }

  double **f = atom->f;
  inter_replica_comm(f);
  if (cmode == SINGLE_PROC)
    nmpimd_transform(bufsortedall, f, M_x2xp[universe->iworld]);
  else if (cmode == MULTI_PROC)
    nmpimd_transform(bufbeads, f, M_x2xp[universe->iworld]);
  trace_nonfinite_state("deferred-force-transform-post", "bead-x-normal-vf");
  normal_mode_force_pending = 0;
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::spring_force()
{
  spring_energy = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double *_mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  tagint *tagtmp = atom->tag;

  int *mask = atom->mask;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      const double *xlast;
      const double *xnext;
      if (cmode == SINGLE_PROC) {
        xlast = bufsortedall[x_last * ntotal + tagtmp[i] - 1];
        xnext = bufsortedall[x_next * ntotal + tagtmp[i] - 1];
      } else {
        xlast = &bufbeads[x_last][3 * i];
        xnext = &bufbeads[x_next][3 * i];
      }

      double delx1 = xlast[0] - x[i][0];
      double dely1 = xlast[1] - x[i][1];
      double delz1 = xlast[2] - x[i][2];

      double delx2 = xnext[0] - x[i][0];
      double dely2 = xnext[1] - x[i][1];
      double delz2 = xnext[2] - x[i][2];

      double ff = fbond * _mass[type[i]];
      // double ff = 0;

      double dx = delx1 + delx2;
      double dy = dely1 + dely2;
      double dz = delz1 + delz2;

      f[i][0] += dx * ff;
      f[i][1] += dy * ff;
      f[i][2] += dz * ff;

      spring_energy += 0.5 * ff * (delx2 * delx2 + dely2 * dely2 + delz2 * delz2);
    }
  }
}

/* ----------------------------------------------------------------------
   Comm operations
   ------------------------------------------------------------------------- */

void FixPIMDLangevin::comm_init()
{
  if (np != universe->nworlds)
    error->all(FLERR, "Fix pimd/langevin: np must equal universe->nworlds");

  int nlocal = atom->nlocal;
  if (cmode == SINGLE_PROC) {
    // init() calls this once per run, so release what an earlier run allocated
    memory->destroy(counts);
    memory->destroy(displacements);
    memory->create(counts, nreplica, "FixPIMDLangevin:counts");
    memory->create(displacements, nreplica, "FixPIMDLangevin:displacements");
    for (int i = 0; i < nreplica; i++) counts[i] = 3 * nlocal;
    displacements[0] = 0;
    for (int i = 0; i < nreplica - 1; i++) displacements[i + 1] = displacements[i] + counts[i];
  }
  if (sizeplan) {
    delete[] plansend;
    delete[] planrecv;
    delete[] modeindex;
  }

  sizeplan = np - 1;
  plansend = new int[sizeplan];
  planrecv = new int[sizeplan];
  modeindex = new int[sizeplan];
  for (int i = 0; i < sizeplan; i++) {

    // send to the (i+1)-th "next" replica, same local rank within that replica
    plansend[i] = universe->me + comm->nprocs * (i + 1);
    if (plansend[i] >= universe->nprocs) plansend[i] -= universe->nprocs;

    // receive from the (i+1)-th "previous" replica, same local rank within that replica
    planrecv[i] = universe->me - comm->nprocs * (i + 1);
    if (planrecv[i] < 0) planrecv[i] += universe->nprocs;

    // where to store what we receive this round:
    // this is the replica index you are pulling from in this step
    modeindex[i] = (universe->iworld + i + 1) % universe->nworlds;
  }

  x_next = (universe->iworld + 1 + universe->nworlds) % (universe->nworlds);
  x_last = (universe->iworld - 1 + universe->nworlds) % (universe->nworlds);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::reallocate_xc()
{
  maxxc = atom->nmax;
  memory->destroy(xc);
  memory->create(xc, maxxc, 3, "FixPIMDLangevin:xc");
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::reallocate_x_unwrap()
{
  maxunwrap = atom->nmax;
  memory->destroy(x_unwrap);
  memory->create(x_unwrap, maxunwrap, 3, "FixPIMDLangevin:x_unwrap");
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::reallocate()
{
  maxlocal = atom->nmax;
  ntotal = atom->natoms;
  if (cmode == SINGLE_PROC) {
    memory->destroy(bufsorted);
    memory->destroy(bufsortedall);
    memory->create(bufsorted, ntotal, 3, "FixPIMDLangevin:bufsorted");
    memory->create(bufsortedall, nreplica * ntotal, 3, "FixPIMDLangevin:bufsortedall");
  } else if (cmode == MULTI_PROC) {
    memory->destroy(bufsend);
    memory->destroy(bufrecv);
    memory->destroy(tagsend);
    memory->destroy(tagrecv);
    memory->destroy(bufbeads);
    memory->create(bufsend, maxlocal * 3, "FixPIMDLangevin:bufsend");
    memory->create(bufrecv, maxlocal * 3, "FixPIMDLangevin:bufrecv");
    memory->create(tagsend, maxlocal, "FixPIMDLangevin:tagsend");
    memory->create(tagrecv, maxlocal, "FixPIMDLangevin:tagrecv");
    memory->create(bufbeads, nreplica, maxlocal * 3, "FixPIMDLangevin:bufrecv");
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::inter_replica_comm(double **ptr)
{
  if (atom->nmax > maxlocal) reallocate();
  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  int i, m;

  // communicate values from the other beads
  if (cmode == SINGLE_PROC) {
    m = 0;
    for (i = 0; i < nlocal; i++) {
      tagint tagtmp = tag[i];
      bufsorted[tagtmp - 1][0] = ptr[i][0];
      bufsorted[tagtmp - 1][1] = ptr[i][1];
      bufsorted[tagtmp - 1][2] = ptr[i][2];
      m++;
    }
    MPI_Allgatherv(bufsorted[0], 3 * m, MPI_DOUBLE, bufsortedall[0], counts, displacements,
                   MPI_DOUBLE, universe->uworld);
  } else if (cmode == MULTI_PROC) {
    // buffers are (re)allocated as needed in reallocate()
    // copy local values
    for (i = 0; i < nlocal; i++) {
      bufbeads[ireplica][3 * i + 0] = ptr[i][0];
      bufbeads[ireplica][3 * i + 1] = ptr[i][1];
      bufbeads[ireplica][3 * i + 2] = ptr[i][2];
    }

    // Loop over replica comm plans
    for (int iplan = 0; iplan < sizeplan; iplan++) {

      // 1) exchange local counts between the paired ranks in universe->uworld
      int nsend = 0;
      MPI_Sendrecv((void *) &nlocal, 1, MPI_INT, plansend[iplan], TAG_INTER_REPLICA_COUNT,
                   (void *) &nsend, 1, MPI_INT, planrecv[iplan], TAG_INTER_REPLICA_COUNT,
                   universe->uworld, MPI_STATUS_IGNORE);

      // 2) ensure buffers sized for nsend
      if (nsend > maxsend) {
        maxsend = nsend + 200;
        tagsend = (tagint *) memory->srealloc(tagsend, sizeof(tagint) * maxsend,
                                              "FixPIMDLangevin:tagsend");
        bufsend = (double *) memory->srealloc(bufsend, sizeof(double) * 3 * maxsend,
                                              "FixPIMDLangevin:bufsend");
      }

      // 3) exchange tags:
      //    send my local tags (atom->tag[0..nlocal-1])
      //    receive remote rank's local tags into tagsend[0..nsend-1]
      MPI_Sendrecv((void *) atom->tag, nlocal, MPI_LMP_TAGINT, plansend[iplan],
                   TAG_INTER_REPLICA_TAGS, (void *) tagsend, nsend, MPI_LMP_TAGINT, planrecv[iplan],
                   TAG_INTER_REPLICA_TAGS, universe->uworld, MPI_STATUS_IGNORE);

      // 4) pack ptr for the tags the remote rank needs from me
      //    For each received tag, find my local index and copy ptr[index][0..2]
      std::vector<int> miss_idx;
      std::vector<tagint> miss_tag;
      miss_idx.reserve(nsend);
      miss_tag.reserve(nsend);

      for (int i = 0; i < nsend; i++) {
        const int idx = atom->map(tagsend[i]);
        if (idx >= 0 && idx < nlocal) {
          bufsend[3 * i + 0] = ptr[idx][0];
          bufsend[3 * i + 1] = ptr[idx][1];
          bufsend[3 * i + 2] = ptr[idx][2];
        } else {
          miss_idx.push_back(i);             // remember which slot in bufsend needs collect
          miss_tag.push_back(tagsend[i]);    // remember which tag that slot corresponds to
        }
      }

      // 5) collect missing tags within this world (local-only claiming)
      // The shared collector uses collective point-to-point exchanges, so all ranks
      // participate whenever any rank has missing tags.
      int has_missing_tags = miss_tag.empty() ? 0 : 1;
      MPI_Allreduce(MPI_IN_PLACE, &has_missing_tags, 1, MPI_INT, MPI_MAX, world);

      std::vector<double> collected_values;
      if (has_missing_tags)
        PIMDUtils::collect_atom_vectors(lmp, style, miss_tag, ptr, collected_values);

      if (!miss_tag.empty()) {
        for (int k = 0; k < static_cast<int>(miss_tag.size()); k++) {
          const int i = miss_idx[k];
          bufsend[3 * i + 0] = collected_values[3 * k + 0];
          bufsend[3 * i + 1] = collected_values[3 * k + 1];
          bufsend[3 * i + 2] = collected_values[3 * k + 2];
        }
      }

      // 6) exchange packed x/f buffers:
      //    - send bufsend (3*nsend) to planrecv[iplan]
      //    - receive bufrecv (3*nlocal) from plansend[iplan]
      //
      // Exchange values with the corresponding rank in the target bead.
      MPI_Sendrecv((void *) bufsend, 3 * nsend, MPI_DOUBLE, planrecv[iplan], TAG_INTER_REPLICA_VALS,
                   (void *) bufrecv, 3 * nlocal, MPI_DOUBLE, plansend[iplan],
                   TAG_INTER_REPLICA_VALS, universe->uworld, MPI_STATUS_IGNORE);

      // 7) store received x/f for this plan into bufbeads[modeindex[iplan]]
      memcpy(bufbeads[modeindex[iplan]], bufrecv, sizeof(double) * 3 * nlocal);
    }
  }
}

void FixPIMDLangevin::remove_com_motion()
{
  if (method == NMPIMD) {
    if (universe->iworld == 0) {
      double **v = atom->v;
      int *mask = atom->mask;
      int nlocal = atom->nlocal;
      if (dynamic) masstotal = group->mass(igroup);
      double vcm[3];
      group->vcm(igroup, masstotal, vcm);
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          v[i][0] -= vcm[0];
          v[i][1] -= vcm[1];
          v[i][2] -= vcm[2];
        }
      }
    }
  } else if (method == PIMD) {
    double **v = atom->v;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    if (dynamic) masstotal = group->mass(igroup);
    double vcm[3];
    group->vcm(igroup, masstotal, vcm);
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        v[i][0] -= vcm[0];
        v[i][1] -= vcm[1];
        v[i][2] -= vcm[2];
      }
    }
  } else {
    error->all(
        FLERR,
        fmt::format("Unknown method for fix {}. Only nmpimd and pimd are supported!", style));
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_xf_vir()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double xf = 0.0;
  vir_ = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      for (int j = 0; j < 3; j++) { xf += x_unwrap[i][j] * atom->f[i][j]; }
    }
  }
  MPI_Allreduce(&xf, &vir_, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_cvir()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double xcf = 0.0;
  centroid_vir = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      for (int j = 0; j < 3; j++) { xcf += (x_unwrap[i][j] - xc[i][j]) * atom->f[i][j]; }
    }
  }
  MPI_Allreduce(&xcf, &centroid_vir, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  if (pstyle == ANISO) {
    for (int i = 0; i < 6; i++) c_vir_tensor[i] = 0.0;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        c_vir_tensor[0] += (x_unwrap[i][0] - xc[i][0]) * atom->f[i][0];
        c_vir_tensor[1] += (x_unwrap[i][1] - xc[i][1]) * atom->f[i][1];
        c_vir_tensor[2] += (x_unwrap[i][2] - xc[i][2]) * atom->f[i][2];
        c_vir_tensor[3] += (x_unwrap[i][0] - xc[i][0]) * atom->f[i][1];
        c_vir_tensor[4] += (x_unwrap[i][0] - xc[i][0]) * atom->f[i][2];
        c_vir_tensor[5] += (x_unwrap[i][1] - xc[i][1]) * atom->f[i][2];
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &c_vir_tensor, 6, MPI_DOUBLE, MPI_SUM, universe->uworld);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_vir()
{
  double volume = domain->xprd * domain->yprd * domain->zprd;
  c_press->compute_vector();
  virial[0] = c_press->vector[0] * volume;
  virial[1] = c_press->vector[1] * volume;
  virial[2] = c_press->vector[2] * volume;
  virial[3] = c_press->vector[3] * volume;
  virial[4] = c_press->vector[4] * volume;
  virial[5] = c_press->vector[5] * volume;
  for (int i = 0; i < 6; i++) virial[i] /= universe->procs_per_world[universe->iworld];
  double vir_bead = (virial[0] + virial[1] + virial[2]);
  MPI_Allreduce(&vir_bead, &vir, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  MPI_Allreduce(MPI_IN_PLACE, &virial[0], 6, MPI_DOUBLE, MPI_SUM, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_stress_tensor()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  if (universe->iworld == 0) {
    double inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    for (int i = 0; i < 6; i++) ke_tensor[i] = 0.0;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        ke_tensor[0] += 0.5 * mass[type[i]] * atom->v[i][0] * atom->v[i][0] * force->mvv2e;
        ke_tensor[1] += 0.5 * mass[type[i]] * atom->v[i][1] * atom->v[i][1] * force->mvv2e;
        ke_tensor[2] += 0.5 * mass[type[i]] * atom->v[i][2] * atom->v[i][2] * force->mvv2e;
        ke_tensor[3] += 0.5 * mass[type[i]] * atom->v[i][0] * atom->v[i][1] * force->mvv2e;
        ke_tensor[4] += 0.5 * mass[type[i]] * atom->v[i][0] * atom->v[i][2] * force->mvv2e;
        ke_tensor[5] += 0.5 * mass[type[i]] * atom->v[i][1] * atom->v[i][2] * force->mvv2e;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &ke_tensor, 6, MPI_DOUBLE, MPI_SUM, world);
    for (int i = 0; i < 6; i++) {
      stress_tensor[i] =
          inv_volume * ((2 * ke_tensor[i] - c_vir_tensor[i]) * force->nktv2p + virial[i]) / np;
    }
  }
  MPI_Bcast(&stress_tensor, 6, MPI_DOUBLE, 0, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_totke()
{
  double kine = 0.0;
  totke = ke_bead = 0.0;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      for (int j = 0; j < 3; j++) { kine += 0.5 * mass[type[i]] * atom->v[i][j] * atom->v[i][j]; }
    }
  }
  kine *= force->mvv2e;
  MPI_Allreduce(&kine, &ke_bead, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&ke_bead, &totke, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  totke /= universe->procs_per_world[universe->iworld];
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_spring_energy()
{
  if (method == NMPIMD) {
    spring_energy = 0.0;
    total_spring_energy = se_bead = 0.0;

    double **x = atom->x;
    double *_mass = atom->mass;
    int *type = atom->type;
    int nlocal = atom->nlocal;
    int *mask = atom->mask;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        spring_energy += 0.5 * _mass[type[i]] * fbond * lam[universe->iworld] *
            (x[i][0] * x[i][0] + x[i][1] * x[i][1] + x[i][2] * x[i][2]);
      }
    }
    MPI_Allreduce(&spring_energy, &se_bead, 1, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(&se_bead, &total_spring_energy, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
    total_spring_energy /= universe->procs_per_world[universe->iworld];
  } else if (method == PIMD) {
    total_spring_energy = se_bead = 0.0;
    MPI_Allreduce(&spring_energy, &se_bead, 1, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(&se_bead, &total_spring_energy, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
    total_spring_energy /= universe->procs_per_world[universe->iworld];
  } else {
    error->universe_all(
        FLERR,
        fmt::format("Unknown method parameter for fix {}. Only nmpimd and pimd are supported!",
                    style));
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_pote()
{
  pe_bead = 0.0;
  pote = 0.0;
  c_pe->compute_scalar();
  pe_bead = c_pe->scalar;
  // PLUMED reports the physical bias B, while this beta/P integrator evolves
  // H_ring + P*B. Complete the energy tally without changing f_plumed or COLVAR.
  for (const auto &fix : modify->get_fix_list()) {
    int dim = -1;
    auto *physical_bias = static_cast<double *>(fix->extract("pimd_physical_bias_energy", dim));
    if (physical_bias && dim == 0)
      pe_bead += (np - (fix->thermo_energy ? 1 : 0)) * (*physical_bias);
  }
  double pot_energy_partition = pe_bead / universe->procs_per_world[universe->iworld];
  MPI_Allreduce(&pot_energy_partition, &pote, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_tote()
{ tote = totke + pote + total_spring_energy; }

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_t_prim()
{
  t_prim = 1.5 * group->count(igroup) * np * force->boltz * temp - total_spring_energy * inverse_np;
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_t_vir()
{
  t_vir = -0.5 * inverse_np * vir_;
  t_cv = 1.5 * group->count(igroup) * force->boltz * temp - 0.5 * inverse_np * centroid_vir;
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_p_prim()
{
  double inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
  p_prim = group->count(igroup) * np * force->boltz * temp * inv_volume -
      1.0 / 1.5 * inv_volume * total_spring_energy;
  p_prim *= force->nktv2p;
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_p_cv()
{
  if (centroid_bias_virial_pending) {
    compute_vir();
    centroid_bias_virial_pending = 0;
  }
  double inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
  p_md = THIRD * inv_volume * (totke + vir);
  if (method == NMPIMD) {
    if (universe->iworld == 0) {
      p_cv = THIRD * inv_volume * ((2.0 * ke_bead - centroid_vir) * force->nktv2p + vir) / np;
    }
    MPI_Bcast(&p_cv, 1, MPI_DOUBLE, 0, universe->uworld);
  } else if (method == PIMD) {
    p_cv = THIRD * inv_volume * ((2.0 * totke / np - centroid_vir) * force->nktv2p + vir) / np;
  } else {
    error->universe_all(
        FLERR,
        fmt::format("Unknown method parameter for fix {}. Only nmpimd and pimd are supported!",
                    style));
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::compute_totenthalpy()
{
  double volume = domain->xprd * domain->yprd * domain->zprd;
  if (barostat == BZP) {
    if (pstyle == ISO) {
      totenthalpy = tote + 0.5 * W * vw[0] * vw[0] * inverse_np + p_hydro * volume / force->nktv2p -
          Vcoeff * kt * log(volume);
    } else if (pstyle == ANISO) {
      totenthalpy = tote + 0.5 * W * vw[0] * vw[0] * inverse_np +
          0.5 * W * vw[1] * vw[1] * inverse_np + 0.5 * W * vw[2] * vw[2] * inverse_np +
          p_hydro * volume / force->nktv2p - Vcoeff * kt * log(volume);
    }
  } else if (barostat == MTTK)
    totenthalpy = tote + 1.5 * W * vw[0] * vw[0] * inverse_np + p_hydro * (volume - vol0);
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixPIMDLangevin::write_restart(FILE *fp)
{
  int nsize = size_restart_global();

  double *list;
  memory->create(list, nsize, "FixPIMDLangevin:list");

  int n = pack_restart_data(list);

  if (tstat_flag) {
    if (comm->me == 0) list[n] = comm->nprocs;

    double state[RanMars::STATE_SIZE];
    random->get_state(state);
    MPI_Gather(state, RanMars::STATE_SIZE, MPI_DOUBLE, list + n + 1, RanMars::STATE_SIZE,
               MPI_DOUBLE, 0, world);
  }

  if (comm->me == 0) {
    int size = nsize * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(list, sizeof(double), nsize, fp);
  }

  memory->destroy(list);
}
/* ---------------------------------------------------------------------- */

int FixPIMDLangevin::size_restart_global()
{
  int nsize = 6;
  if (tstat_flag) nsize += 1 + comm->nprocs * RanMars::STATE_SIZE;

  return nsize;
}

/* ---------------------------------------------------------------------- */

int FixPIMDLangevin::pack_restart_data(double *list)
{
  int n = 0;
  for (int i = 0; i < 6; i++) list[n++] = vw[i];
  return n;
}

/* ---------------------------------------------------------------------- */

void FixPIMDLangevin::restart(char *buf, int nbytes)
{
  constexpr int barostat_size = 6;
  if (nbytes < barostat_size * static_cast<int>(sizeof(double)) || nbytes % sizeof(double) != 0)
    error->all(FLERR, "Invalid fix pimd/langevin restart data size");

  auto *list = reinterpret_cast<double *>(buf);
  for (int i = 0; i < barostat_size; i++) vw[i] = list[i];

  // Older LAMMPS versions stored only the barostat state.
  const int nvalues = nbytes / sizeof(double);
  if (nvalues == barostat_size) {
    if (tstat_flag && comm->me == 0)
      error->warning(FLERR, "Legacy fix pimd/langevin restart has no PILE_L RNG state");
    return;
  }

  const double stored_nprocs = list[barostat_size];
  const int available = nvalues - barostat_size - 1;
  if (!std::isfinite(stored_nprocs) || stored_nprocs < 1 ||
      stored_nprocs > available / RanMars::STATE_SIZE || stored_nprocs != std::floor(stored_nprocs))
    error->all(FLERR, "Invalid fix pimd/langevin restart RNG count");
  const int restart_nprocs = static_cast<int>(stored_nprocs);
  if (available != restart_nprocs * RanMars::STATE_SIZE)
    error->all(FLERR, "Invalid fix pimd/langevin restart RNG data size");

  if (tstat_flag) {
    if (restart_nprocs != comm->nprocs) {
      if (comm->me == 0)
        error->warning(FLERR, "Different number of procs. Cannot restore PILE_L RNG state.");
    } else {
      double state[RanMars::STATE_SIZE];
      const double *saved = list + barostat_size + 1 + comm->me * RanMars::STATE_SIZE;
      std::copy(saved, saved + RanMars::STATE_SIZE, state);
      // The development implementation used the same 105-value layout with
      // an unused zero in the marker slot. It also stored the Gaussian cache.
      state[0] = -1.0;
      random->set_state(state);
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixPIMDLangevin::compute_vector(int n)
{
  if (n == 0) return ke_bead;
  if (n == 1) return se_bead;
  if (n == 2) return pe_bead;
  if (n == 3) return tote;
  if (n == 4) return t_prim;
  if (n == 5) return t_vir;
  if (n == 6) return t_cv;
  if (n == 7) return p_prim;
  if (n == 8) return p_md;
  if (n == 9) return p_cv;

  if (pstat_flag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    if (pstyle == ISO) {
      if (n == 10) return vw[0];
      if (barostat == BZP) {
        if (n == 11) return 0.5 * W * vw[0] * vw[0];
      } else if (barostat == MTTK) {
        if (n == 11) return 1.5 * W * vw[0] * vw[0];
      }
      if (n == 12) { return np * Pext * volume / force->nktv2p; }
      if (n == 13) { return -Vcoeff * np * kt * log(volume); }
      if (n == 14) return totenthalpy;
    } else if (pstyle == ANISO) {
      if (n == 10) return vw[0];
      if (n == 11) return vw[1];
      if (n == 12) return vw[2];
      if (n == 13) return 0.5 * W * (vw[0] * vw[0] + vw[1] * vw[1] + vw[2] * vw[2]);
      if (n == 14) { return np * Pext * volume / force->nktv2p; }
      if (n == 15) {
        double volume = domain->xprd * domain->yprd * domain->zprd;
        return -Vcoeff * np * kt * log(volume);
      }
      if (n == 16) return totenthalpy;
    }
  }
  return 0.0;
}
