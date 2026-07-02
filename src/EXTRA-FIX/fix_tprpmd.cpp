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

#include "fix_tprpmd.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "universe.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_2PI;
using MathConst::MY_PI;
using MathConst::MY_SQRT2;
using MathConst::THIRD;
using MathSpecial::powint;

enum { TPPIMD, TPRPMD };
enum { PHYSICAL, NORMAL };
enum { BAOAB, OBABO };
enum { ISO, ANISO, TRICLINIC };
enum { NHC };
enum { MTTK, BZP };
enum { NVE, NVT, NPH, NPT, UVT };
enum { SINGLE_PROC, MULTI_PROC };

/* ---------------------------------------------------------------------- */

FixTPRPMD::FixTPRPMD(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), mass(nullptr), plansend(nullptr), planrecv(nullptr), tagsend(nullptr),
    tagrecv(nullptr), bufsend(nullptr), bufrecv(nullptr), bufbeads(nullptr), bufsorted(nullptr),
    bufsortedall(nullptr), tagsendall(nullptr),
    tagrecvall(nullptr), bufsendall(nullptr), bufrecvall(nullptr), counts(nullptr),
    displacements(nullptr), rootworld(MPI_COMM_NULL), lam(nullptr), M_x2xp(nullptr), M_xp2x(nullptr), modeindex(nullptr), tau_k(nullptr), _omega_k(nullptr), Lan_s(nullptr),
    Lan_c(nullptr), xc(nullptr), xcall(nullptr),
    x_unwrap(nullptr), id_pe(nullptr), id_press(nullptr), c_pe(nullptr), c_press(nullptr),
    eta(nullptr), eta_dot(nullptr), eta_dotdot(nullptr), eta_mass(nullptr), Ne_dot(nullptr),
    Ne_mass(nullptr), Ne(nullptr), u_start(0.0), u_stop(0.0), u_current(0.0), u_target(0.0),
    dedn_name(nullptr), dedn_which(ArgInfo::NONE), dedn_index(0), dedn_var(-1),
    dedn_compute(nullptr), dedn_fix(nullptr), dedn_current(0.0)
{
  restart_global = 1;
  time_integrate = 1;
  global_freq = 1;
  vector_flag = 1;
  extvector = -1;
  size_vector = 10;

  ntotal = 0;
  maxlocal = maxunwrap = maxxc = 0;
  sizeplan = 0;

  method = TPRPMD;
  ensemble = UVT;
  integrator = OBABO;
  thermostat = NHC;
  lj_epsilon = 1;
  lj_sigma = 1;
  lj_mass = 1;
  other_planck = 1;
  other_mvv2e = 1;
  fmass = 1.0;
  np = universe->nworlds;
  inverse_np = 1.0 / np;
  sp = 1.0;
  temp = 298.15;
  tau = 1.0;
  pilescale = 1.0;
  tstat_flag = 1;
  ustat_flag = 1;
  mapflag = 1;
  removecomflag = 1;
  fmmode = PHYSICAL;
  pote = tote = totke = total_spring_energy = 0.0;
  centroid_vir = vir = vir_ = 0.0;
  ke_bead = se_bead = pe_bead = tote = t_prim = t_vir = t_cv = p_prim = p_md = p_cv = 0.0;

  mtchain = 3;
  nc_tchain = 1;
  t_period = 0.0;
  u_period = 0.0;
  drag = 0.0;
  mu = -3.5;

  if (domain->dimension != 3) error->universe_all(FLERR, "Fix tprpmd requires a 3d system");
  if (narg < 4) utils::missing_cmd_args(FLERR, std::string("fix ") + style, error);

  bool mu_seen = false;
  bool ne_seen = false;
  bool dedn_seen = false;

  for (int i = 3; i < narg;) {
    if (strcmp(arg[i], "method") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} method", style), error);
      if (strcmp(arg[i + 1], "tp-pimd") == 0)
        method = TPPIMD;
      else if (strcmp(arg[i + 1], "tp-rpmd") == 0)
        method = TPRPMD;
      else if (strcmp(arg[i + 1], "tprpmd") == 0)
        method = TPRPMD;
      else
        error->all(FLERR, "Fix {} only supports method tp-pimd or tp-rpmd", style);
      i += 2;
    } else if (strcmp(arg[i], "ensemble") == 0) {
      if (i + 2 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} ensemble", style), error);
      if (strcmp(arg[i + 1], "uvt") != 0)
        error->all(FLERR, "Fix {} only supports ensemble uvt", style);
      i += 2;
    } else if (strcmp(arg[i], "thermostat") == 0) {
      if (i + 2 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} thermostat", style), error);
      if (strcmp(arg[i + 1], "NHC") != 0)
        error->all(FLERR, "Fix {} only supports thermostat NHC", style);
      i += 2;
    } else if (strcmp(arg[i], "integrator") == 0) {
      if (i + 2 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} integrator", style), error);
      if (strcmp(arg[i + 1], "obabo") == 0)
        integrator = OBABO;
      else if (strcmp(arg[i + 1], "baoab") == 0)
        integrator = BAOAB;
      else
        error->all(FLERR, "Unknown integrator parameter for fix {}", style);
      i += 2;
    } else if (strcmp(arg[i], "temp") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} temp", style), error);
      temp = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (temp < 0.0) error->all(FLERR, "Invalid temp value for fix {}", style);
      i += 2;
    } else if (strcmp(arg[i], "Tdamp") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} Tdamp", style), error);
      t_period = utils::numeric(FLERR, arg[i + 1], false, lmp);
      i += 2;
    } else if (strcmp(arg[i], "tchain") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} tchain", style), error);
      mtchain = utils::inumeric(FLERR, arg[i + 1], false, lmp);
      i += 2;
    } else if (strcmp(arg[i], "tloop") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} tloop", style), error);
      nc_tchain = utils::inumeric(FLERR, arg[i + 1], false, lmp);
      i += 2;
    } else if (strcmp(arg[i], "drag") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} drag", style), error);
      drag = utils::numeric(FLERR, arg[i + 1], false, lmp);
      i += 2;
    } else if (strcmp(arg[i], "fmass") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} fmass", style), error);
      fmass = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (fmass < 0.0 || fmass > np) error->all(FLERR, "Invalid fmass value for fix {}", style);
      i += 2;
    } else if (strcmp(arg[i], "fmmode") == 0) {
      if (i + 2 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} fmmode", style), error);
      if (strcmp(arg[i + 1], "physical") == 0)
        fmmode = PHYSICAL;
      else if (strcmp(arg[i + 1], "normal") == 0)
        fmmode = NORMAL;
      else
        error->all(FLERR, "Unknown fictitious mass mode for fix {}", style);
      i += 2;
    } else if (strcmp(arg[i], "sp") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} sp", style), error);
      sp = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (sp < 0.0) error->all(FLERR, "Invalid sp value for fix {}", style);
      i += 2;
    } else if (strcmp(arg[i], "lj") == 0) {
      if (i + 6 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} lj", style), error);
      lj_epsilon = utils::numeric(FLERR, arg[i + 1], false, lmp);
      lj_sigma = utils::numeric(FLERR, arg[i + 2], false, lmp);
      lj_mass = utils::numeric(FLERR, arg[i + 3], false, lmp);
      other_planck = utils::numeric(FLERR, arg[i + 4], false, lmp);
      other_mvv2e = utils::numeric(FLERR, arg[i + 5], false, lmp);
      i += 6;
    } else if (strcmp(arg[i], "tau") == 0) {
      if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} tau", style), error);
      tau = utils::numeric(FLERR, arg[i + 1], false, lmp);
      i += 2;
    } else if (strcmp(arg[i], "mu") == 0) {
      if (i + 4 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} mu", style), error);
      mu_seen = true;
      u_start = utils::numeric(FLERR, arg[i + 1], false, lmp);
      u_stop = utils::numeric(FLERR, arg[i + 2], false, lmp);
      u_target = u_start;
      mu = u_target;
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
    } else if (strcmp(arg[i], "removecom") == 0) {
      if (i + 2 > narg)
        utils::missing_cmd_args(FLERR, fmt::format("fix {} removecom", style), error);
      removecomflag = utils::logical(FLERR, arg[i + 1], false, lmp);
      i += 2;
    } else if ((strcmp(arg[i], "barostat") == 0) || (strcmp(arg[i], "iso") == 0) ||
               (strcmp(arg[i], "aniso") == 0) || (strcmp(arg[i], "taup") == 0) ||
               (strcmp(arg[i], "fixedpoint") == 0)) {
      error->all(FLERR, "Pressure control is not supported by fix {}", style);
    } else if ((strcmp(arg[i], "seed") == 0) || (strcmp(arg[i], "PILE_L_temp") == 0)) {
      error->all(FLERR, "Legacy thermostat options are not supported by fix {}", style);
    } else {
      error->all(FLERR, "Unknown keyword {} for fix {}", arg[i], style);
    }
  }

  if (!mu_seen) error->all(FLERR, "Missing mu keyword for fix {}", style);
  if (!ne_seen) error->all(FLERR, "Missing ne keyword for fix {}", style);
  if (!dedn_seen) error->all(FLERR, "Missing dedn keyword for fix {}", style);
  if (t_period <= 0.0) error->all(FLERR, "Temperature damping for fix {} must be > 0.0", style);
  if (u_period <= 0.0) error->all(FLERR, "Chemical-potential damping for fix {} must be > 0.0", style);

  if (tstat_flag) {
    int ich;
    eta = new double[mtchain];
    eta_dot = new double[mtchain + 1];
    eta_dot[mtchain] = 0.0;
    eta_dotdot = new double[mtchain];
    for (ich = 0; ich < mtchain; ich++) eta[ich] = eta_dot[ich] = eta_dotdot[ich] = 0.0;
    eta_mass = new double[mtchain];
    size_vector += 4 * mtchain;
  }

  if (!Ne_dot) {
    Ne_dot = new double[1];
    *Ne_dot = 0.0;
  }
  Ne_mass = new double[1];
  *Ne_mass = 0.0;
  size_vector += 6;
  extlist = new int[size_vector];
  for (int i = 0; i < size_vector; i++) extlist[i] = 1;
  for (int i = size_vector - 6; i < size_vector; i++) extlist[i] = 0;
  u_freq = 1.0 / u_period;
  u_current = u_start;

  id_pe = utils::strdup(std::string(id) + "_pimd_pe");
  modify->add_compute(fmt::format("{} all pe", id_pe));

  id_press = utils::strdup(std::string(id) + "_pimd_press");
  modify->add_compute(fmt::format("{} all pressure NULL virial", id_press));

  ntotal = atom->natoms;
  nreplica = np;

  if (mass == nullptr) mass = new double[atom->ntypes + 1];
  for (int i = 1; i <= atom->ntypes; i++) mass[i] = atom->mass[i] * fmass;

  fixedpoint[0] = 0.5 * (domain->boxlo[0] + domain->boxhi[0]);
  fixedpoint[1] = 0.5 * (domain->boxlo[1] + domain->boxhi[1]);
  fixedpoint[2] = 0.5 * (domain->boxlo[2] + domain->boxhi[2]);
}

/* ---------------------------------------------------------------------- */

FixTPRPMD::~FixTPRPMD()
{
  modify->delete_compute(id_pe);
  modify->delete_compute(id_press);
  delete[] id_pe;
  delete[] id_press;
  delete[] extlist;
  delete[] mass;
  delete[] _omega_k;
  delete[] Lan_c;
  delete[] Lan_s;
  delete[] tau_k;
  delete[] plansend;
  delete[] planrecv;
  delete[] modeindex;
  memory->sfree(lam);
  memory->destroy(xcall);
  if (cmode == SINGLE_PROC) {
    memory->destroy(bufsorted);
    memory->destroy(bufsortedall);
    memory->destroy(counts);
    memory->destroy(displacements);
  }

  if (cmode == MULTI_PROC) {
    memory->destroy(bufsendall);
    memory->destroy(bufrecvall);
    memory->destroy(tagsendall);
    memory->destroy(tagrecvall);
    memory->destroy(counts);
    memory->destroy(displacements);
  }
  memory->destroy(M_x2xp);
  memory->destroy(M_xp2x);
  memory->destroy(xc);
  memory->destroy(x_unwrap);
  memory->destroy(bufsend);
  memory->destroy(bufrecv);
  memory->destroy(tagsend);
  memory->destroy(tagrecv);
  memory->destroy(bufbeads);

  if (thermostat == NHC) {
    if (tstat_flag) {
      delete [] eta;
      delete [] eta_dot;
      delete [] eta_dotdot;
      delete [] eta_mass;
    }
    if (ustat_flag) {
      delete [] Ne;
      delete [] Ne_dot;
      delete [] Ne_mass;
    }
  }
  if (rootworld != MPI_COMM_NULL) MPI_Comm_free(&rootworld);
  delete[] dedn_name;
}

/* ---------------------------------------------------------------------- */

bool FixTPRPMD::nuclear_thermostat_off() const
{
  return method == TPRPMD && np != 1 && universe->iworld == 0;
}

/* ---------------------------------------------------------------------- */

bool FixTPRPMD::ne_thermostat_participates() const
{
  if (!ustat_flag) return false;
  if (method != TPRPMD) return true;
  return (np == 1) ? true : (universe->iworld != 0);
}

/* ---------------------------------------------------------------------- */

bool FixTPRPMD::nhc_chain_active() const
{
  return !nuclear_thermostat_off() || ne_thermostat_participates();
}

/* ---------------------------------------------------------------------- */

double FixTPRPMD::chain0_target_energy() const
{
  return nuclear_thermostat_off() ? 0.0 : static_cast<double>(np) * ke_target;
}

/* ---------------------------------------------------------------------- */

double FixTPRPMD::ne_kinetic_current_share() const
{
  if (!ne_thermostat_participates()) return 0.0;
  return 0.5 * inverse_np * (*Ne_mass) * (*Ne_dot) * (*Ne_dot);
}

/* ---------------------------------------------------------------------- */

int FixTPRPMD::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::init()
{
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, "Fix tprpmd requires an atom map, see atom_modify");

  if (comm->nprocs != 1)
    error->all(FLERR, "Fix tprpmd currently requires one MPI rank per bead");

  if (universe->me == 0 && universe->uscreen)
    fprintf(universe->uscreen, "Fix tprpmd: initializing Path-Integral ...\n");

  // prepare the constants

  masstotal = group->mass(igroup);

  double planck;
  if (strcmp(update->unit_style, "lj") == 0) {
    double planck_star = sqrt(lj_epsilon) * sqrt(lj_mass) * lj_sigma * sqrt(other_mvv2e);
    planck = other_planck / planck_star;
  } else {
    planck = force->hplanck;
  }
  planck *= sp;
  hbar = planck / (MY_2PI);
  kt = force->boltz * temp;
  beta = 1.0 / kt;
  double _fbond = 1.0 * np * np / (beta * beta * hbar * hbar);

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
        error->all(FLERR, "Compute {} for fix {} vector index {} is out of range", dedn_name, style,
                   dedn_index);
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

  omega_np = np / (hbar * beta) * sqrt(force->mvv2e);
  beta_np = 1.0 / force->boltz / temp * inverse_np;
  fbond = _fbond * force->mvv2e;

  if ((universe->me == 0) && (universe->uscreen))
    fprintf(universe->uscreen,
            "Fix tprpmd: -P/(beta^2 * hbar^2) = %20.7lE (kcal/mol/A^2)\n\n", fbond);

  me = comm->me;
  nprocs = comm->nprocs;
  cmode = (nprocs == 1) ? SINGLE_PROC : MULTI_PROC;

  nprocs_universe = universe->nprocs;
  nreplica = universe->nworlds;
  ireplica = universe->iworld;
  mapflag = (nreplica == 1) ? 0 : 1;

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

  if (integrator == OBABO) {
    dtf = 0.5 * update->dt * force->ftm2v;
    dtv = 0.5 * update->dt;
    dtv2 = dtv * dtv;
    dtv3 = THIRD * dtv2 * dtv * force->ftm2v;
    dthalf = 0.5 * update->dt;
    dt4 = 0.25 * update->dt;
    dt8 = 0.125 * update->dt;
  } else if (integrator == BAOAB) {
    dtf = 0.5 * update->dt * force->ftm2v;
    dtv = 0.5 * update->dt;
    dtv2 = dtv * dtv;
    dtv3 = THIRD * dtv2 * dtv * force->ftm2v;
    dthalf = 0.5 * update->dt;
    dt4 = 0.25 * update->dt;
    dt8 = 0.125 * update->dt;
  } else {
    error->universe_all(FLERR, "Unknown integrator parameter for fix tprpmd");
  }

  comm_init();
  if (mass == nullptr) mass = new double[atom->ntypes + 1];
  if (xcall == nullptr) memory->create(xcall, ntotal * 3, "FixTPRPMD:xcall");
  nmpimd_init();
  nhc_init();

  c_pe = modify->get_compute_by_id(id_pe);
  if (!c_pe)
    error->universe_all(
        FLERR, fmt::format("Could not find fix {} potential energy compute ID {}", style, id_pe));

  c_press = modify->get_compute_by_id(id_press);
  if (!c_press)
    error->universe_all(
        FLERR, fmt::format("Could not find fix {} pressure compute ID {}", style, id_press));

  t_prim = t_vir = t_cv = p_cv = p_md = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::setup(int vflag)
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  imageint *image = atom->image;

  compute_mu_target();
  *Ne_mass = tdof * force->boltz * temp / (u_freq * u_freq);

  if (mapflag) {
    for (int i = 0; i < nlocal; i++) domain->unmap(x[i], image[i]);
  }

  inter_replica_comm(x);
  if (cmode == SINGLE_PROC)
    nmpimd_transform(bufsortedall, x, M_x2xp[universe->iworld]);
  else if (cmode == MULTI_PROC)
    nmpimd_transform(bufbeads, x, M_x2xp[universe->iworld]);
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

  post_force(vflag);
  compute_totke();
  end_of_step();
  c_pe->addstep(update->ntimestep + 1);
  c_press->addstep(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::initial_integrate(int /*vflag*/)
{
  double *Ne = this->Ne;
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
    }
    b_step();
    inter_replica_comm(x);
    if (cmode == SINGLE_PROC)
      nmpimd_transform(bufsortedall, x, M_x2xp[universe->iworld]);
    else if (cmode == MULTI_PROC)
      nmpimd_transform(bufbeads, x, M_x2xp[universe->iworld]);
    qc_step();
    a_step();
    qc_step();
    a_step();
  } else if (integrator == BAOAB) {
    b_step();
    inter_replica_comm(x);
    if (cmode == SINGLE_PROC)
      nmpimd_transform(bufsortedall, x, M_x2xp[universe->iworld]);
    else if (cmode == MULTI_PROC)
      nmpimd_transform(bufbeads, x, M_x2xp[universe->iworld]);
    qc_step();
    a_step();
    if (tstat_flag) {
      o_step();
      if (removecomflag) remove_com_motion();
    }
    qc_step();
    a_step();
  } else {
    error->universe_all(FLERR,
                        "Unknown integrator parameter for fix tprpmd. Only obabo and baoab "
                        "integrators are supported!");
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
    for (int i = 0; i < nlocal; i++) { domain->unmap_inv(x[i], image[i]); }
  }
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::final_integrate()
{
  b_step();

  if (integrator == OBABO) {
    if (tstat_flag) {
      o_step();
      if (removecomflag) remove_com_motion();
    }
  } else if (integrator == BAOAB) {

  } else {
    error->universe_all(FLERR, "Unknown integrator parameter for fix tprpmd");
  }
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::post_force(int /*flag*/)
{
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

  compute_pote();
  inter_replica_comm(f);
  if (cmode == SINGLE_PROC)
    nmpimd_transform(bufsortedall, f, M_x2xp[universe->iworld]);
  else if (cmode == MULTI_PROC)
    nmpimd_transform(bufbeads, f, M_x2xp[universe->iworld]);

  // Sample the FD dE/dN signal once per force evaluation and reuse it
  // throughout the rest of this MD step.
  if (ustat_flag) refresh_dedn_cache();

  c_pe->addstep(update->ntimestep + 1);
  c_press->addstep(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::end_of_step()
{
  compute_totke();
  compute_p_cv();
  compute_tote();
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::collect_xc()
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  tagint *tag = atom->tag;
  if (ireplica == 0) {
    if (cmode == SINGLE_PROC) {
      for (int i = 0; i < nlocal; i++) {
        xcall[3 * i + 0] = xcall[3 * i + 1] = xcall[3 * i + 2] = 0.0;
      }
    } else if (cmode == MULTI_PROC) {
      for (int i = 0; i < ntotal; i++) {
        xcall[3 * i + 0] = xcall[3 * i + 1] = xcall[3 * i + 2] = 0.0;
      }
    }

    const double sqrtnp = sqrt((double) np);
    for (int i = 0; i < nlocal; i++) {
      xcall[3 * (tag[i] - 1) + 0] = x[i][0] / sqrtnp;
      xcall[3 * (tag[i] - 1) + 1] = x[i][1] / sqrtnp;
      xcall[3 * (tag[i] - 1) + 2] = x[i][2] / sqrtnp;
    }

    if (cmode == MULTI_PROC) {
      MPI_Allreduce(MPI_IN_PLACE, xcall, ntotal * 3, MPI_DOUBLE, MPI_SUM, world);
    }
  }
  MPI_Bcast(xcall, ntotal * 3, MPI_DOUBLE, 0, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::b_step()
{
  // The force here is the external potential contribution in normal-mode space;
  // the harmonic ring-polymer propagation is handled separately in the A/QC steps.
  int n = atom->nlocal;
  int *type = atom->type;
  double **v = atom->v;
  double **f = atom->f;

  for (int i = 0; i < n; i++) {
    double dtfm = dtf / mass[type[i]];
    v[i][0] += dtfm * f[i][0];
    v[i][1] += dtfm * f[i][1];
    v[i][2] += dtfm * f[i][2];
  }

  if (ustat_flag) {
    double dtfm = dthalf / *Ne_mass;
    if (universe->iworld == 0) *Ne_dot += dtfm * (-dedn_current + u_target);
    MPI_Bcast(Ne_dot, 1, MPI_DOUBLE, 0, universe->uworld);
  }
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::qc_step()
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;
  double *Ne = this->Ne;
  if (universe->iworld == 0) {
    for (int i = 0; i < nlocal; i++) {
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
    if (ustat_flag) *Ne += dtv * (*Ne_dot);
  }
  if (ustat_flag) MPI_Bcast(Ne, 1, MPI_DOUBLE, 0, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::a_step()
{
  // Use the analytical harmonic-oscillator update for the non-centroid modes.
  int n = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;
  double x0, x1, x2, v0, v1, v2;    // three components of x[i] and v[i]

  if (universe->iworld != 0) {
    for (int i = 0; i < n; i++) {
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

/* ---------------------------------------------------------------------- */

void FixTPRPMD::nhc_init()
{
  if (kt <= 0.0 || hbar <= 0.0)
    error->universe_all(FLERR, "Fix tprpmd requires positive kt and hbar in nhc_init");

  const double beta_local = 1.0 / kt;
  const double _omega_np = np / beta_local / hbar;
  double _omega_np_dt_half = _omega_np * update->dt * 0.5;

  _omega_k = new double[np];
  Lan_c = new double[np];
  Lan_s = new double[np];
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

  if (tstat_flag) {
    t_freq = 1.0 / t_period;
    tdrag_factor = 1.0 - (update->dt * t_freq * drag / nc_tchain);
  }

  if (ustat_flag) {
    u_freq = 1.0 / u_period;
  }

  int fix_dof = 0;
  for (auto &ifix : modify->get_fix_list())
    if (ifix->dof_flag){
      fix_dof += ifix->dof(igroup);
    }
  int extra_dof = domain->dimension;
  tdof = domain->dimension * group->count(igroup);
  tdof -= extra_dof + fix_dof;

  ke_target = tdof * force->boltz * temp;

  if (tstat_flag) {
    const double chain0_target = chain0_target_energy();
    const double chain_target = nhc_chain_active() ? static_cast<double>(np) * force->boltz * temp : 0.0;
    eta_mass[0] = chain0_target / (t_freq * t_freq);
    for (int ich = 1; ich < mtchain; ich++) {
      eta_mass[ich] = (np == 1) ? chain_target / (t_freq * t_freq)
                                : chain_target / (_omega_np * _omega_np);
      if (eta_mass[ich] > 0.0)
        eta_dotdot[ich] =
            (eta_mass[ich - 1] * eta_dot[ich - 1] * eta_dot[ich - 1] - chain_target) /
            eta_mass[ich];
      else
        eta_dotdot[ich] = 0.0;
    }
    if (!nhc_chain_active()) {
      for (int ich = 1; ich < mtchain; ich++) eta_dotdot[ich] = 0;
    }
  }

  if (ustat_flag) *Ne_mass = tdof * force->boltz * temp / (u_freq * u_freq);

  std::string out = "Initializing TP path-integral Nose-Hoover thermostat chain...\n";
  out += "  Bead ID    |    omega    |    tau\n";
  tau_k = new double[np];
  tau_k[0] = tau;
  for (int i = 1; i < np; i++) tau_k[i] = 0.5 / pilescale / _omega_k[i];
  for (int i = 0; i < np; i++) {
    out += fmt::format("      {:d}     {:.8e} {:.8e}\n", i, _omega_k[i], tau_k[i]);
  }
  out += "  NHC thermostat successfully initialized!\n\n";
  if (universe->me == 0) utils::logmesg(lmp, out);
}


/* ---------------------------------------------------------------------- */

void FixTPRPMD::o_step()
{
  if (tstat_flag) {
    if (ustat_flag) {
      compute_mu_target();
      nhc_mu_integrate();
    }
    else nhc_temp_integrate();
  }
}

void FixTPRPMD::nhc_temp_integrate()
{
  int ich;
  double expfac;
  int *type = atom->type;
  double **v = atom->v;
  int nlocal = atom->nlocal;
  double kecurrent = 0, t_current;

  for (int i = 0; i < nlocal; i++) {
    kecurrent += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * mass[type[i]];
  }
  kecurrent *= force->mvv2e;

  t_current = kecurrent / force->boltz / tdof;

  if (nuclear_thermostat_off()) return;

  const double chain0_target = chain0_target_energy();
  const double chain_target = nhc_chain_active() ? static_cast<double>(np) * force->boltz * temp : 0.0;
  if (eta_mass[0] > 0.0)
    eta_dotdot[0] = (kecurrent - chain0_target) / eta_mass[0];
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

    // rescale temperature due to velocity scaling
    // should not be necessary to explicitly recompute the temperature

    t_current *= factor_eta*factor_eta;
    kecurrent = tdof * force->boltz * t_current;

    if (eta_mass[0] > 0.0)
      eta_dotdot[0] = (kecurrent - chain0_target) / eta_mass[0];
    else eta_dotdot[0] = 0.0;

    for (ich = 0; ich < mtchain; ich++)
      eta[ich] += ncfac*dthalf*eta_dot[ich];

    eta_dot[0] *= expfac;
    eta_dot[0] += eta_dotdot[0] * ncfac*dt4;
    eta_dot[0] *= expfac;

    for (ich = 1; ich < mtchain; ich++) {
      expfac = exp(-ncfac*dt8*eta_dot[ich+1]);
      eta_dot[ich] *= expfac;
      eta_dotdot[ich] = (eta_mass[ich-1]*eta_dot[ich-1]*eta_dot[ich-1] - chain_target) /
          eta_mass[ich];
      eta_dot[ich] += eta_dotdot[ich] * ncfac*dt4;
      eta_dot[ich] *= expfac;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::nhc_mu_integrate()
{
  int ich;
  double expfac;
  int *type = atom->type;
  double **v = atom->v;
  int nlocal = atom->nlocal;
  double kecurrent = 0, t_current;

  for (int i = 0; i < nlocal; i++) {
    kecurrent += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * mass[type[i]];
  }
  kecurrent *= force->mvv2e;

  t_current = kecurrent / force->boltz / tdof;

  const double chain0_target = chain0_target_energy();
  const double chain_target = nhc_chain_active() ? static_cast<double>(np) * force->boltz * temp : 0.0;
  if (method == TPRPMD && np != 1 && universe->iworld == 0)
    ;
  else {
    if (eta_mass[0] > 0.0)
      eta_dotdot[0] = (kecurrent + ne_kinetic_current_share() - chain0_target) / eta_mass[0];
    else eta_dotdot[0] = 0.0;
  }

  double ncfac = 1.0 / nc_tchain;
  for (int iloop = 0; iloop < nc_tchain; iloop++) {

    if (method == TPRPMD && np != 1 && universe->iworld == 0)
      ;
    else {

      for (ich = mtchain - 1; ich > 0; ich--) {
        expfac = exp(-ncfac * dt8 * eta_dot[ich + 1]);
        eta_dot[ich] *= expfac;
        eta_dot[ich] += eta_dotdot[ich] * ncfac * dt4;
        eta_dot[ich] *= tdrag_factor;
        eta_dot[ich] *= expfac;
      }

      expfac = exp(-ncfac * dt8 * eta_dot[1]);
      eta_dot[0] *= expfac;
      eta_dot[0] += eta_dotdot[0] * ncfac * dt4;
      eta_dot[0] *= tdrag_factor;
      eta_dot[0] *= expfac;

      factor_eta = exp(-ncfac * dthalf * eta_dot[0]);
      nh_v_temp();
    }

    double eta_dot_k = eta_dot[0], eta_dot_ave = 0;
    MPI_Allreduce(&eta_dot_k, &eta_dot_ave, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);

    if (method == TPRPMD && np != 1) {
      eta_dot_ave = eta_dot_ave / (static_cast<double>(np) - 1);
    } else {
      eta_dot_ave *= inverse_np;
    }

    *Ne_dot *= exp(-ncfac * dthalf * eta_dot_ave);

    // rescale temperature due to velocity scaling
    // should not be necessary to explicitly recompute the temperature
    if (method == TPRPMD && np != 1 && universe->iworld == 0)
      ;
    else {
      t_current *= factor_eta * factor_eta;
      kecurrent = tdof * force->boltz * t_current;

      if (eta_mass[0] > 0.0)
        eta_dotdot[0] =
            (kecurrent + ne_kinetic_current_share() - chain0_target) / eta_mass[0];
      else eta_dotdot[0] = 0.0;

      for (ich = 0; ich < mtchain; ich++)
        eta[ich] += ncfac * dthalf * eta_dot[ich];

      eta_dot[0] *= expfac;
      eta_dot[0] += eta_dotdot[0] * ncfac * dt4;
      eta_dot[0] *= expfac;

      for (ich = 1; ich < mtchain; ich++) {
        expfac = exp(-ncfac * dt8 * eta_dot[ich + 1]);
        eta_dot[ich] *= expfac;
        eta_dotdot[ich] =
            (eta_mass[ich - 1] * eta_dot[ich - 1] * eta_dot[ich - 1] - chain_target) /
            eta_mass[ich];
        eta_dot[ich] += eta_dotdot[ich] * ncfac * dt4;
        eta_dot[ich] *= expfac;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::nh_v_temp()
{
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
}

/* ----------------------------------------------------------------------
   Normal-mode ring-polymer operations
   ------------------------------------------------------------------------- */

void FixTPRPMD::nmpimd_init()
{
  memory->create(M_x2xp, np, np, "fix_feynman:M_x2xp");
  memory->create(M_xp2x, np, np, "fix_feynman:M_xp2x");

  lam = (double *) memory->smalloc(sizeof(double) * np, "FixTPRPMD::lam");

  // Set up  eigenvalues
  for (int i = 0; i < np; i++) {
    double sin_tmp = sin(i * MY_PI / np);
    lam[i] = 4 * sin_tmp * sin_tmp;
  }

  for (int i = 0; i < np; i++) {
    if (!std::isfinite(lam[i])) error->universe_all(FLERR, "Fix tprpmd encountered invalid lambda");
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

void FixTPRPMD::nmpimd_transform(double **src, double **des, double *vector)
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

    for (int i = 0; i < n; i++)
      for (int d = 0; d < 3; d++) {
        des[i][d] = 0.0;
        for (int j = 0; j < np; j++) { des[i][d] += (src[j][m] * vector[j]); }
        m++;
      }
  }
}

/* ----------------------------------------------------------------------
   Comm operations
   ------------------------------------------------------------------------- */

void FixTPRPMD::comm_init()
{
  if (np != universe->nworlds) error->all(FLERR, "Fix tprpmd: np must equal universe->nworlds");

  int nlocal = atom->nlocal;
  if (cmode == SINGLE_PROC) {
    memory->destroy(counts);
    memory->destroy(displacements);
    memory->create(counts, nreplica, "FixTPRPMD:counts");
    memory->create(displacements, nreplica, "FixTPRPMD:displacements");
    for (int i = 0; i < nreplica; i++) counts[i] = 3 * nlocal;
    displacements[0] = 0;
    for (int i = 0; i < nreplica - 1; i++) displacements[i + 1] = displacements[i] + counts[i];
  }

  if (sizeplan) {
    delete[] plansend;
    delete[] planrecv;
  }

  sizeplan = np - 1;
  plansend = new int[sizeplan];
  planrecv = new int[sizeplan];
  modeindex = new int[sizeplan];
  for (int i = 0; i < sizeplan; i++) {
    int isend = ireplica + i + 1;
    if (isend >= nreplica) isend -= nreplica;

    int irecv = ireplica - (i + 1);
    if (irecv < 0) irecv += nreplica;

    plansend[i] = universe->root_proc[isend];
    planrecv[i] = universe->root_proc[irecv];
    modeindex[i] = irecv;
  }

  x_next = (universe->iworld + 1 + universe->nworlds) % (universe->nworlds);
  x_last = (universe->iworld - 1 + universe->nworlds) % (universe->nworlds);
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::reallocate_xc()
{
  maxxc = atom->nmax;
  memory->destroy(xc);
  memory->create(xc, maxxc, 3, "FixTPRPMD:xc");
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::reallocate_x_unwrap()
{
  maxunwrap = atom->nmax;
  memory->destroy(x_unwrap);
  memory->create(x_unwrap, maxunwrap, 3, "FixTPRPMD:x_unwrap");
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::reallocate()
{
  maxlocal = atom->nmax;
  ntotal = atom->natoms;
  if (cmode == SINGLE_PROC) {
    memory->destroy(bufsorted);
    memory->destroy(bufsortedall);
    memory->create(bufsorted, ntotal, 3, "FixTPRPMD:bufsorted");
    memory->create(bufsortedall, nreplica * ntotal, 3, "FixTPRPMD:bufsortedall");
  } else if (cmode == MULTI_PROC) {
    memory->destroy(bufsend);
    memory->destroy(bufrecv);
    memory->destroy(tagsend);
    memory->destroy(tagrecv);
    memory->destroy(bufbeads);
    memory->create(bufsend, maxlocal, 3, "FixTPRPMD:bufsend");
    memory->create(bufrecv, maxlocal, 3, "FixTPRPMD:bufrecv");
    memory->create(tagsend, maxlocal, "FixTPRPMD:tagsend");
    memory->create(tagrecv, maxlocal, "FixTPRPMD:tagrecv");
    memory->create(bufbeads, nreplica, maxlocal * 3, "FixTPRPMD:bufbeads");
  }
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::inter_replica_comm(double **ptr)
{
  MPI_Request requests[2];
  MPI_Status statuses[2];
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
    for (i = 0; i < nlocal; i++) {
      bufbeads[ireplica][3 * i + 0] = ptr[i][0];
      bufbeads[ireplica][3 * i + 1] = ptr[i][1];
      bufbeads[ireplica][3 * i + 2] = ptr[i][2];
    }
    m = 0;
    for (i = 0; i < nlocal; i++) {
      tagsend[m] = tag[i];
      bufsend[m][0] = ptr[i][0];
      bufsend[m][1] = ptr[i][1];
      bufsend[m][2] = ptr[i][2];
      m++;
    }
    MPI_Gather(&m, 1, MPI_INT, counts, 1, MPI_INT, 0, world);
    displacements[0] = 0;
    for (i = 0; i < nprocs - 1; i++) displacements[i + 1] = displacements[i] + counts[i];
    MPI_Gatherv(tagsend, m, MPI_LMP_TAGINT, tagsendall, counts, displacements, MPI_LMP_TAGINT, 0,
                world);
    for (i = 0; i < nprocs; i++) counts[i] *= 3;
    for (i = 0; i < nprocs - 1; i++) displacements[i + 1] = displacements[i] + counts[i];
    MPI_Gatherv(bufsend[0], 3 * m, MPI_DOUBLE, bufsendall[0], counts, displacements, MPI_DOUBLE, 0,
                world);
    for (int iplan = 0; iplan < sizeplan; iplan++) {
      if (me == 0) {
        MPI_Irecv(bufrecvall[0], 3 * ntotal, MPI_DOUBLE, planrecv[iplan], 0, universe->uworld,
                  &requests[0]);
        MPI_Irecv(tagrecvall, ntotal, MPI_LMP_TAGINT, planrecv[iplan], 0, universe->uworld,
                  &requests[1]);
        MPI_Send(bufsendall[0], 3 * ntotal, MPI_DOUBLE, plansend[iplan], 0, universe->uworld);
        MPI_Send(tagsendall, ntotal, MPI_LMP_TAGINT, plansend[iplan], 0, universe->uworld);
        MPI_Waitall(2, requests, statuses);
      }
      MPI_Bcast(tagrecvall, ntotal, MPI_LMP_TAGINT, 0, world);
      MPI_Bcast(bufrecvall[0], 3 * ntotal, MPI_DOUBLE, 0, world);
      for (i = 0; i < ntotal; i++) {
        m = atom->map(tagrecvall[i]);
        if (m < 0 || m >= nlocal) continue;
        bufbeads[modeindex[iplan]][3 * m + 0] = bufrecvall[i][0];
        bufbeads[modeindex[iplan]][3 * m + 1] = bufrecvall[i][1];
        bufbeads[modeindex[iplan]][3 * m + 2] = bufrecvall[i][2];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::remove_com_motion()
{
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
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::compute_xf_vir()
{
  int nlocal = atom->nlocal;
  double xf = 0.0;
  vir_ = 0.0;
  for (int i = 0; i < nlocal; i++) {
    for (int j = 0; j < 3; j++) { xf += x_unwrap[i][j] * atom->f[i][j]; }
  }
  MPI_Allreduce(&xf, &vir_, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::compute_cvir()
{
  int nlocal = atom->nlocal;
  double xcf = 0.0;
  centroid_vir = 0.0;
  for (int i = 0; i < nlocal; i++) {
    for (int j = 0; j < 3; j++) { xcf += (x_unwrap[i][j] - xc[i][j]) * atom->f[i][j]; }
  }
  MPI_Allreduce(&xcf, &centroid_vir, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::compute_vir()
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

void FixTPRPMD::compute_totke()
{
  double kine = 0.0;
  totke = ke_bead = 0.0;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  for (int i = 0; i < nlocal; i++) {
    for (int j = 0; j < 3; j++) { kine += 0.5 * mass[type[i]] * atom->v[i][j] * atom->v[i][j]; }
  }
  kine *= force->mvv2e;
  MPI_Allreduce(&kine, &ke_bead, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&ke_bead, &totke, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  totke /= universe->procs_per_world[universe->iworld];
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::compute_spring_energy()
{
  spring_energy = 0.0;
  total_spring_energy = se_bead = 0.0;

  double **x = atom->x;
  double *_mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    spring_energy += 0.5 * _mass[type[i]] * fbond * lam[universe->iworld] *
        (x[i][0] * x[i][0] + x[i][1] * x[i][1] + x[i][2] * x[i][2]);
  }
  MPI_Allreduce(&spring_energy, &se_bead, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&se_bead, &total_spring_energy, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  total_spring_energy /= universe->procs_per_world[universe->iworld];
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::compute_pote()
{
  pe_bead = 0.0;
  pote = 0.0;
  c_pe->compute_scalar();
  pe_bead = c_pe->scalar;
  double pot_energy_partition = pe_bead / universe->procs_per_world[universe->iworld];
  MPI_Allreduce(&pot_energy_partition, &pote, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::compute_tote()
{
  tote = totke + pote + total_spring_energy;
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::compute_t_prim()
{
  t_prim = 1.5 * atom->natoms * np * force->boltz * temp - total_spring_energy * inverse_np;
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::compute_t_vir()
{
  t_vir = -0.5 * inverse_np * vir_;
  t_cv = 1.5 * atom->natoms * force->boltz * temp - 0.5 * inverse_np * centroid_vir;
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::compute_p_prim()
{
  double inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
  p_prim = atom->natoms * np * force->boltz * temp * inv_volume -
      1.0 / 1.5 * inv_volume * total_spring_energy;
  p_prim *= force->nktv2p;
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::compute_p_cv()
{
  double inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
  p_md = THIRD * inv_volume * (totke + vir);
  if (universe->iworld == 0) {
    p_cv = THIRD * inv_volume * ((2.0 * ke_bead - centroid_vir) * force->nktv2p + vir) / np;
  }
  MPI_Bcast(&p_cv, 1, MPI_DOUBLE, 0, universe->uworld);
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixTPRPMD::write_restart(FILE *fp)
{
  int nsize = size_restart_global();

  double *list;
  memory->create(list, nsize, "FixTPRPMD:list");

  pack_restart_data(list);

  if (comm->me == 0) {
    int size = nsize * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(list, sizeof(double), nsize, fp);
  }

  memory->destroy(list);
}
/* ---------------------------------------------------------------------- */

int FixTPRPMD::size_restart_global()
{
  int nsize = 2;
  if (tstat_flag) nsize += 1 + 2 * mtchain;
  if (ustat_flag) nsize += 2;
  return nsize;
}

/* ---------------------------------------------------------------------- */

int FixTPRPMD::pack_restart_data(double *list)
{
  int n = 0;
  list[n++] = tstat_flag;
  if (tstat_flag) {
    list[n++] = mtchain;
    for (int ich = 0; ich < mtchain; ich++) list[n++] = eta[ich];
    for (int ich = 0; ich < mtchain; ich++) list[n++] = eta_dot[ich];
  }
  list[n++] = ustat_flag;
  if (ustat_flag) {
    list[n++] = *Ne;
    list[n++] = *Ne_dot;
  }
  return n;
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::restart(char *buf)
{
  int n = 0;
  auto list = (double *) buf;
  int flag = static_cast<int>(list[n++]);
  if (flag) {
    int m = static_cast<int>(list[n++]);
    if (tstat_flag && m == mtchain) {
      if (!nhc_chain_active()) {
        for (int ich = 0; ich < mtchain; ich++) eta[ich] = eta_dot[ich] = 0.0;
        n += 2 * m;
      } else {
        for (int ich = 0; ich < mtchain; ich++) eta[ich] = list[n++];
        for (int ich = 0; ich < mtchain; ich++) eta_dot[ich] = list[n++];
      }
    } else {
      n += 2 * m;
    }
  }
  flag = static_cast<int>(list[n++]);
  if (flag && ustat_flag) {
    *Ne = list[n++];
    *Ne_dot = list[n++];
  }
}

/* ----------------------------------------------------------------------
   if thermostat == NHC, return a single element of the following vectors, in this order:
      eta[tchain], eta_dot[tchain], PE_eta[tchain], KE_eta_dot[tchain]
------------------------------------------------------------------------- */

double FixTPRPMD::compute_vector(int n)
{
  if (n == 0) return ke_bead;
  n -= 1;
  if (n == 0) return se_bead;
  n -= 1;
  if (n == 0) return pe_bead;
  n -= 1;
  if (n == 0) return tote;
  n -= 1;
  if (n == 0) return t_prim;
  n -= 1;
  if (n == 0) return t_vir;
  n -= 1;
  if (n == 0) return t_cv;
  n -= 1;
  if (n == 0) return p_prim;
  n -= 1;
  if (n == 0) return p_md;
  n -= 1;
  if (n == 0) return p_cv;
  n -= 1;

  int ilen;
  double *Ne = this->Ne;

  if (tstat_flag) {
    ilen = mtchain;
    if (n < ilen) return eta[n];
    n -= ilen;
    ilen = mtchain;
    if (n < ilen) return eta_dot[n];
    n -= ilen;
  }

  double kt = force->boltz * temp;
  int ich;

  if (tstat_flag) {
    ilen = mtchain;
    if (n < ilen) {
      ich = n;
      if (ich == 0) return chain0_target_energy() * eta[0];
      return (nhc_chain_active() ? kt : 0.0) * eta[ich];
    }
    n -= ilen;
    ilen = mtchain;
    if (n < ilen) {
      ich = n;
      if (ich == 0) return 0.5 * eta_mass[0] * eta_dot[0] * eta_dot[0];
      return 0.5 * eta_mass[ich] * eta_dot[ich] * eta_dot[ich];
    }
    n -= ilen;
  }
  if (ustat_flag) {
    ilen = 1;
    if (n < ilen) return *Ne;
    n -= ilen;
    ilen = 1;
    if (n < ilen) return *Ne_dot;
    n -= ilen;
    ilen = 1;
    if (n < ilen) return dedn_current;
    n -= ilen;
    ilen = 1;
    if (n < ilen) return u_target;
    n -= ilen;
    ilen = 1;
    if (n < ilen) return 0.5 * (*Ne_mass) * (*Ne_dot) * (*Ne_dot);
    n -= ilen;
    ilen = 1;
    if (n < ilen) return kt * eta[0] - u_target * (*Ne);
  }

  return 0.0;
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::compute_mu_target()
{
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  u_target = u_start + delta * (u_stop - u_start);
  mu = u_target;
}

/* ---------------------------------------------------------------------- */

double FixTPRPMD::evaluate_dedn()
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
      // Force a fresh compute evaluation when the fix explicitly refreshes
      // its cached dE/dN value for the current force state.
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

void FixTPRPMD::refresh_dedn_cache()
{
  double dedn_local = evaluate_dedn();
  double dedn_avg = 0.0;
  MPI_Allreduce(&dedn_local, &dedn_avg, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  dedn_current = dedn_avg * inverse_np;
  u_current = dedn_current;
}

/* ---------------------------------------------------------------------- */

void FixTPRPMD::parse_dedn_source(const char *arg)
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
