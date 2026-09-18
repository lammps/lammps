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

#include "fix_pimd_nve.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_2PI;
using MathConst::MY_PI;
using MathConst::MY_SQRT2;
using MathSpecial::powint;

namespace {
constexpr int TAG_INTER_REPLICA_COUNT = 10;
constexpr int TAG_INTER_REPLICA_TAGS  = 11;
constexpr int TAG_INTER_REPLICA_VALS  = 12;

constexpr int TAG_RING_MISS_COUNT = 400;
constexpr int TAG_RING_MISS_TAGS  = 401;
constexpr int TAG_RING_REP_COUNT  = 402;
constexpr int TAG_RING_REP_TAGS   = 403;
constexpr int TAG_RING_REP_VALS   = 404;
} // namespace

enum { PHYSICAL, NORMAL };
enum { SINGLE_PROC, MULTI_PROC };

/* ---------------------------------------------------------------------- */

FixPIMDNVE::FixPIMDNVE(LAMMPS *lmp, int narg, char **arg, bool defer_setup) :
    Fix(lmp, narg, arg), mass(nullptr), rootworld(MPI_COMM_NULL), plansend(nullptr),
    planrecv(nullptr), tagsend(nullptr), tagrecv(nullptr), bufsend(nullptr), bufrecv(nullptr),
    bufbeads(nullptr), bufsorted(nullptr), bufsortedall(nullptr), counts(nullptr),
    displacements(nullptr), lam(nullptr), M_x2xp(nullptr), M_xp2x(nullptr), modeindex(nullptr),
    _omega_k(nullptr), Lan_s(nullptr), Lan_c(nullptr), xc(nullptr), xcall(nullptr),
    x_unwrap(nullptr), id_pe(nullptr), id_press(nullptr), c_pe(nullptr), c_press(nullptr)
{
  restart_global = 1;
  time_integrate = 1;
  global_freq = 1;
  vector_flag = 1;
  extvector = -1;
  size_vector = 10;

  ntotal = 0;
  maxlocal = maxunwrap = maxxc = 0;
  sizeplan = maxsend = 0;
  me = nprocs = ireplica = nreplica = nprocs_universe = 0;
  x_last = x_next = 0;
  cmode = -1;

  method = NMPIMD;
  integrator = OBABO;
  lj_epsilon = 1.0;
  lj_sigma = 1.0;
  lj_mass = 1.0;
  other_planck = 1.0;
  other_mvv2e = 1.0;
  fmass = 1.0;
  np = universe->nworlds;
  inverse_np = 1.0 / np;
  sp = 1.0;
  temp = 298.15;
  mapflag = 1;
  removecomflag = 1;
  fmmode = PHYSICAL;

  pote = tote = totke = total_spring_energy = 0.0;
  centroid_vir = vir = vir_ = 0.0;
  ke_bead = se_bead = pe_bead = t_prim = t_vir = t_cv = p_prim = p_md = p_cv = 0.0;
  kt = 0.0;
  beta = 0.0;
  beta_np = 0.0;
  hbar = 0.0;
  omega_np = 0.0;
  fbond = 0.0;

  if (domain->dimension != 3)
    error->universe_all(FLERR, fmt::format("Fix {} requires a 3d system", style));
  if (narg < 3) utils::missing_cmd_args(FLERR, std::string("fix ") + style, error);
  // Derived styles initialize their own options before parsing the full command.
  if (defer_setup) return;

  // process keywords

  for (int i = 3; i < narg;) {
    if (!parse_keyword(narg, arg, i))
      error->all(FLERR, "Unknown keyword {} for fix {}", arg[i], style);
  }
  if (method == CMD) error->all(FLERR, "Fix pimd/nve does not support method cmd");
  finish_constructor_setup();
}

/* ---------------------------------------------------------------------- */

bool FixPIMDNVE::parse_keyword(int narg, char **arg, int &i)
{
  if (strcmp(arg[i], "method") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} method", style), error);
    if (strcmp(arg[i + 1], "nmpimd") == 0)
      method = NMPIMD;
    else if (strcmp(arg[i + 1], "pimd") == 0)
      method = PIMD;
    else if (strcmp(arg[i + 1], "cmd") == 0)
      method = CMD;
    else
      error->all(FLERR, "Unknown method parameter for fix {}", style);
    i += 2;
    return true;
  }
  if (strcmp(arg[i], "integrator") == 0) {
    if (i + 2 > narg)
      utils::missing_cmd_args(FLERR, fmt::format("fix {} integrator", style), error);
    if (strcmp(arg[i + 1], "obabo") == 0)
      integrator = OBABO;
    else if (strcmp(arg[i + 1], "baoab") == 0)
      integrator = BAOAB;
    else
      error->all(FLERR, "Unknown integrator parameter for fix {}", style);
    i += 2;
    return true;
  }
  if (strcmp(arg[i], "temp") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} temp", style), error);
    temp = utils::numeric(FLERR, arg[i + 1], false, lmp);
    if (temp < 0.0) error->all(FLERR, "Invalid temp value for fix {}", style);
    i += 2;
    return true;
  }
  if (strcmp(arg[i], "fmass") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} fmass", style), error);
    fmass = utils::numeric(FLERR, arg[i + 1], false, lmp);
    if (fmass <= 0.0 || fmass > np) error->all(FLERR, "Invalid fmass value for fix {}", style);
    i += 2;
    return true;
  }
  if (strcmp(arg[i], "fmmode") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} fmmode", style), error);
    if (strcmp(arg[i + 1], "physical") == 0)
      fmmode = PHYSICAL;
    else if (strcmp(arg[i + 1], "normal") == 0)
      fmmode = NORMAL;
    else
      error->all(FLERR, "Unknown fictitious mass mode for fix {}", style);
    i += 2;
    return true;
  }
  if (strcmp(arg[i], "sp") == 0) {
    if (i + 2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} sp", style), error);
    sp = utils::numeric(FLERR, arg[i + 1], false, lmp);
    if (sp < 0.0) error->all(FLERR, "Invalid sp value for fix {}", style);
    i += 2;
    return true;
  }
  if (strcmp(arg[i], "lj") == 0) {
    if (i + 6 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix {} lj", style), error);
    lj_epsilon = utils::numeric(FLERR, arg[i + 1], false, lmp);
    lj_sigma = utils::numeric(FLERR, arg[i + 2], false, lmp);
    lj_mass = utils::numeric(FLERR, arg[i + 3], false, lmp);
    other_planck = utils::numeric(FLERR, arg[i + 4], false, lmp);
    other_mvv2e = utils::numeric(FLERR, arg[i + 5], false, lmp);
    i += 6;
    return true;
  }
  if (strcmp(arg[i], "removecom") == 0) {
    if (i + 2 > narg)
      utils::missing_cmd_args(FLERR, fmt::format("fix {} removecom", style), error);
    removecomflag = utils::logical(FLERR, arg[i + 1], false, lmp);
    i += 2;
    return true;
  }
  return false;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::finish_constructor_setup()
{
  extlist = new int[size_vector];
  for (int i = 0; i < size_vector; i++) extlist[i] = 1;

  id_pe = utils::strdup(std::string(id) + "_pimd_pe");
  modify->add_compute(fmt::format("{} all pe", id_pe));

  id_press = utils::strdup(std::string(id) + "_pimd_press");
  modify->add_compute(fmt::format("{} all pressure NULL virial", id_press));

  ntotal = atom->natoms;
  nreplica = np;

  if (mass == nullptr) mass = new double[atom->ntypes + 1];
  for (int i = 1; i <= atom->ntypes; i++) mass[i] = atom->mass[i] * fmass;
}

/* ---------------------------------------------------------------------- */

FixPIMDNVE::~FixPIMDNVE()
{
  if (id_pe) modify->delete_compute(id_pe);
  if (id_press) modify->delete_compute(id_press);
  delete[] id_pe;
  delete[] id_press;
  delete[] extlist;
  delete[] mass;
  delete[] _omega_k;
  delete[] Lan_c;
  delete[] Lan_s;
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

int FixPIMDNVE::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::init()
{
  if (atom->tag_consecutive() == 0)
    error->all(FLERR, "Atom IDs must be consecutive for fix {}", style);
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, "Fix {} requires an atom map, see atom_modify", style);

  if (universe->me == 0 && universe->uscreen)
    fprintf(universe->uscreen, "Fix %s: initializing Path-Integral ...\n", style);

  masstotal = group->mass(igroup);

  double planck;
  if (strcmp(update->unit_style, "lj") == 0) {
    double planck_star = sqrt(lj_epsilon) * sqrt(lj_mass) * lj_sigma * sqrt(other_mvv2e);
    planck = other_planck / planck_star;
  } else {
    planck = force->hplanck;
  }
  planck *= sp;
  hbar = planck / MY_2PI;
  kt = force->boltz * temp;
  beta = 1.0 / kt;
  double bond_prefactor = static_cast<double>(np) * static_cast<double>(np) / (beta * beta * hbar * hbar);

  // hbar and kBT use the same energy unit, so their ratio is a frequency.
  // mvv2e belongs in the spring energy/force coefficient, not this frequency.
  omega_np = np / (hbar * beta);
  beta_np = 1.0 / force->boltz / temp * inverse_np;
  fbond = bond_prefactor * force->mvv2e;

  if ((universe->me == 0) && (universe->uscreen))
    fprintf(universe->uscreen, "Fix %s: -P/(beta^2 * hbar^2) = %20.7lE (kcal/mol/A^2)\n\n", style,
            fbond);

  me = comm->me;
  nprocs = comm->nprocs;
  cmode = (nprocs == 1) ? SINGLE_PROC : MULTI_PROC;
  if (method == PIMD && cmode == MULTI_PROC)
    error->universe_all(FLERR, "Method pimd only supports a single processor per bead");
  if (method == PIMD && fmmode == NORMAL)
    error->universe_all(FLERR, "Normal mode mass is not supported for method pimd");
  nprocs_universe = universe->nprocs;
  nreplica = universe->nworlds;
  ireplica = universe->iworld;
  mapflag = (nreplica == 1) ? 0 : 1;

  int *iroots = new int[nreplica];
  MPI_Group uworldgroup, rootgroup;
  for (int i = 0; i < nreplica; i++) iroots[i] = universe->root_proc[i];
  MPI_Comm_group(universe->uworld, &uworldgroup);
  MPI_Group_incl(uworldgroup, nreplica, iroots, &rootgroup);
  if (rootworld != MPI_COMM_NULL) MPI_Comm_free(&rootworld);
  MPI_Comm_create(universe->uworld, rootgroup, &rootworld);
  if (rootgroup != MPI_GROUP_NULL) MPI_Group_free(&rootgroup);
  if (uworldgroup != MPI_GROUP_NULL) MPI_Group_free(&uworldgroup);
  delete[] iroots;

  ntotal = atom->natoms;
  if (atom->nmax > maxlocal) reallocate();
  if (atom->nmax > maxunwrap) reallocate_x_unwrap();
  if (atom->nmax > maxxc) reallocate_xc();

  dtf = 0.5 * update->dt * force->ftm2v;
  dtv = 0.5 * update->dt;
  dtv2 = dtv * dtv;
  dtv3 = (1.0 / 3.0) * dtv2 * dtv * force->ftm2v;

  comm_init();
  if (mass == nullptr) mass = new double[atom->ntypes + 1];
  if (xcall == nullptr) memory->create(xcall, ntotal * 3, "FixPIMDNVE:xcall");
  nmpimd_init();

  c_pe = modify->get_compute_by_id(id_pe);
  if (!c_pe)
    error->universe_all(
        FLERR, fmt::format("Could not find fix {} potential energy compute ID {}", style, id_pe));

  c_press = modify->get_compute_by_id(id_press);
  if (!c_press)
    error->universe_all(
        FLERR, fmt::format("Could not find fix {} pressure compute ID {}", style, id_press));

  if (!c_pe->peflag)
    error->all(FLERR, "Compute ID {} for fix {} does not compute potential energy", id_pe, style);
  if (!c_press->pressflag)
    error->all(FLERR, "Compute ID {} for fix {} does not compute pressure", id_press, style);

  setup_subclass_state();
  t_prim = t_vir = t_cv = p_prim = p_cv = p_md = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::setup(int vflag)
{
  if (method == NMPIMD || method == CMD) {
    unmap_coordinates(atom->x, atom->image);
    // Forward: bead coordinates to normal modes.
    inter_replica_comm(atom->x);
    nmpimd_transform(normal_mode_transform_buffer(), atom->x, M_x2xp[universe->iworld]);
  } else if (method == PIMD) {
    unmap_coordinates(atom->x, atom->image);
    prepare_coordinates();
    spring_force();
  } else {
    error->universe_all(FLERR, fmt::format("Unknown method parameter for fix {}", style));
  }
  after_force_transform_hook();
  collect_xc();
  compute_spring_energy();
  compute_t_prim();
  compute_p_prim();
  if (method == NMPIMD || method == CMD) {
    // Backward: normal modes to bead coordinates.
    inter_replica_comm(atom->x);
    nmpimd_transform(normal_mode_transform_buffer(), atom->x, M_xp2x[universe->iworld]);
  }
  remap_coordinates(atom->x, atom->image);

  post_force(vflag);
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::initial_integrate(int /*vflag*/)
{
  if (integrator != OBABO && integrator != BAOAB)
    error->universe_all(FLERR, fmt::format("Unknown integrator parameter for fix {}", style));

  if (integrator == OBABO) o_step();
  b_step();
  unmap_coordinates(atom->x, atom->image);
  if (method == NMPIMD || method == CMD) {
    // Forward: bead coordinates to normal modes.
    inter_replica_comm(atom->x);
    nmpimd_transform(normal_mode_transform_buffer(), atom->x, M_x2xp[universe->iworld]);
    qc_step();
    a_step();
    if (integrator == BAOAB) o_step();
    qc_step();
    a_step();
    collect_xc();
    compute_spring_energy();
    compute_t_prim();
    compute_p_prim();
    // Backward: normal modes to bead coordinates.
    inter_replica_comm(atom->x);
    nmpimd_transform(normal_mode_transform_buffer(), atom->x, M_xp2x[universe->iworld]);
  } else if (method == PIMD) {
    q_step();
    if (integrator == BAOAB) o_step();
    q_step();
    collect_xc();
  } else {
    error->universe_all(FLERR, fmt::format("Unknown method parameter for fix {}", style));
  }
  remap_coordinates(atom->x, atom->image);
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::final_integrate()
{
  b_step();
  if (integrator == OBABO) o_step();
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::o_step() {}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::prepare_coordinates()
{
  inter_replica_comm(atom->x);
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::post_force(int /*flag*/)
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  imageint *image = atom->image;
  tagint *tag = atom->tag;

  if (atom->nmax > maxunwrap) reallocate_x_unwrap();
  if (atom->nmax > maxxc) reallocate_xc();

  for (int i = 0; i < nlocal; i++) {
    x_unwrap[i][0] = x[i][0];
    x_unwrap[i][1] = x[i][1];
    x_unwrap[i][2] = x[i][2];
  }
  unmap_coordinates(x_unwrap, image);
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
    unmap_coordinates(atom->x, atom->image);
    prepare_coordinates();
    spring_force();
    compute_spring_energy();
    compute_t_prim();
    remap_coordinates(atom->x, atom->image);
  }

  compute_pote();
  if (method == NMPIMD || method == CMD) {
    // Forward: bead forces to normal-mode forces.
    inter_replica_comm(atom->f);
    nmpimd_transform(normal_mode_transform_buffer(), atom->f, M_x2xp[universe->iworld]);
  }
  after_force_transform_hook();

  c_pe->addstep(update->ntimestep + 1);
  c_press->addstep(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::end_of_step()
{
  compute_totke();
  compute_p_cv();
  compute_tote();
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::setup_subclass_state() {}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::after_force_transform_hook() {}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::unmap_coordinates(double **coords, imageint *image)
{
  if (!mapflag) return;

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) domain->unmap(coords[i], image[i]);
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::remap_coordinates(double **coords, imageint *image)
{
  if (!mapflag) return;

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) domain->unmap_inv(coords[i], image[i]);
}

/* ---------------------------------------------------------------------- */

double **FixPIMDNVE::normal_mode_transform_buffer()
{
  if (cmode == SINGLE_PROC) return bufsortedall;
  return bufbeads;
}

/* ---------------------------------------------------------------------- */

int FixPIMDNVE::subclass_restart_size() const
{
  return 0;
}

/* ---------------------------------------------------------------------- */

int FixPIMDNVE::pack_subclass_restart(double *, int n) const
{
  return n;
}

/* ---------------------------------------------------------------------- */

int FixPIMDNVE::unpack_subclass_restart(const double *, int n)
{
  return n;
}

/* ---------------------------------------------------------------------- */

double FixPIMDNVE::compute_subclass_vector(int) const
{
  return 0.0;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::collect_xc()
{
  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  double **x = atom->x;
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

void FixPIMDNVE::b_step()
{
  // For NMPIMD, force only includes the contribution of external potential.
  // For PIMD, force includes the contributions of external potential and spring force.
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  double **v = atom->v;
  double **f = atom->f;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    double dtfm = dtf / mass[type[i]];
    v[i][0] += dtfm * f[i][0];
    v[i][1] += dtfm * f[i][1];
    v[i][2] += dtfm * f[i][2];
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::q_step()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    x[i][0] += dtv * v[i][0];
    x[i][1] += dtv * v[i][1];
    x[i][2] += dtv * v[i][2];
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::qc_step()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  if (universe->iworld == 0) {
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::a_step()
{
  int n = atom->nlocal;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;

  if (universe->iworld != 0) {
    for (int i = 0; i < n; i++) {
      if (!(mask[i] & groupbit)) continue;
      double x0 = x[i][0];
      double x1 = x[i][1];
      double x2 = x[i][2];
      double v0 = v[i][0];
      double v1 = v[i][1];
      double v2 = v[i][2];
      x[i][0] = Lan_c[universe->iworld] * x0 +
          1.0 / _omega_k[universe->iworld] * Lan_s[universe->iworld] * v0;
      x[i][1] = Lan_c[universe->iworld] * x1 +
          1.0 / _omega_k[universe->iworld] * Lan_s[universe->iworld] * v1;
      x[i][2] = Lan_c[universe->iworld] * x2 +
          1.0 / _omega_k[universe->iworld] * Lan_s[universe->iworld] * v2;
      v[i][0] = -_omega_k[universe->iworld] * Lan_s[universe->iworld] * x0 +
          Lan_c[universe->iworld] * v0;
      v[i][1] = -_omega_k[universe->iworld] * Lan_s[universe->iworld] * x1 +
          Lan_c[universe->iworld] * v1;
      v[i][2] = -_omega_k[universe->iworld] * Lan_s[universe->iworld] * x2 +
          Lan_c[universe->iworld] * v2;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::spring_force()
{
  spring_energy = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double *_mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    double delx1 = bufsortedall[x_last * nlocal + tag[i] - 1][0] - x[i][0];
    double dely1 = bufsortedall[x_last * nlocal + tag[i] - 1][1] - x[i][1];
    double delz1 = bufsortedall[x_last * nlocal + tag[i] - 1][2] - x[i][2];

    double delx2 = bufsortedall[x_next * nlocal + tag[i] - 1][0] - x[i][0];
    double dely2 = bufsortedall[x_next * nlocal + tag[i] - 1][1] - x[i][1];
    double delz2 = bufsortedall[x_next * nlocal + tag[i] - 1][2] - x[i][2];

    double ff = fbond * _mass[type[i]];
    f[i][0] += (delx1 + delx2) * ff;
    f[i][1] += (dely1 + dely2) * ff;
    f[i][2] += (delz1 + delz2) * ff;

    spring_energy += 0.5 * ff * (delx2 * delx2 + dely2 * dely2 + delz2 * delz2);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::nmpimd_init()
{
  if (kt <= 0.0 || hbar <= 0.0)
    error->universe_all(FLERR, fmt::format("Fix {} requires positive kt and hbar", style));

  memory->destroy(M_x2xp);
  memory->destroy(M_xp2x);
  memory->sfree(lam);
  delete[] _omega_k;
  delete[] Lan_c;
  delete[] Lan_s;

  memory->create(M_x2xp, np, np, "fix_pimd_nve:M_x2xp");
  memory->create(M_xp2x, np, np, "fix_pimd_nve:M_xp2x");
  lam = (double *) memory->smalloc(sizeof(double) * np, "FixPIMDNVE::lam");
  _omega_k = new double[np];
  Lan_c = new double[np];
  Lan_s = new double[np];

  for (int i = 0; i < np; i++) {
    double sin_tmp = sin(i * MY_PI / np);
    lam[i] = 4.0 * sin_tmp * sin_tmp;
  }

  const double sqrtnp = sqrt((double) np);
  for (int j = 0; j < np; j++) {
    for (int i = 1; i < int(np / 2) + 1; i++) {
      M_x2xp[i][j] = MY_SQRT2 * cos(MY_2PI * double(i) * double(j) / double(np)) / sqrtnp;
    }
    for (int i = int(np / 2) + 1; i < np; i++) {
      M_x2xp[i][j] = MY_SQRT2 * sin(MY_2PI * double(i) * double(j) / double(np)) / sqrtnp;
    }
  }
  for (int i = 0; i < np; i++) {
    M_x2xp[0][i] = 1.0 / sqrtnp;
    if (np % 2 == 0) M_x2xp[np / 2][i] = 1.0 / sqrtnp * powint(-1.0, i);
  }
  for (int i = 0; i < np; i++)
    for (int j = 0; j < np; j++) M_xp2x[i][j] = M_x2xp[j][i];

  init_normal_mode_coefficients();

  int iworld = universe->iworld;
  for (int i = 1; i <= atom->ntypes; i++) {
    mass[i] = atom->mass[i] * fmass;
    if (iworld != 0 && fmmode == NORMAL) mass[i] *= lam[iworld];
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::init_normal_mode_coefficients()
{
  if (fmmode != PHYSICAL && fmmode != NORMAL)
    error->universe_all(FLERR, "Unknown fmmode setting; only physical and normal are supported!");

  // The centroid is free; only internal modes have a spring frequency.
  _omega_k[0] = 0.0;
  Lan_c[0] = 1.0;
  Lan_s[0] = 0.0;
  for (int i = 1; i < np; i++) {
    _omega_k[i] = omega_np / sqrt(fmass);
    if (fmmode == PHYSICAL) _omega_k[i] *= sqrt(lam[i]);
    const double angle = _omega_k[i] * update->dt * 0.5;
    Lan_c[i] = cos(angle);
    Lan_s[i] = sin(angle);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::nmpimd_transform(double **src, double **des, double *vector)
{
  if (cmode == SINGLE_PROC) {
    for (int i = 0; i < ntotal; i++) {
      for (int d = 0; d < 3; d++) {
        bufsorted[i][d] = 0.0;
        for (int j = 0; j < nreplica; j++) bufsorted[i][d] += src[j * ntotal + i][d] * vector[j];
      }
    }
    for (int i = 0; i < ntotal; i++) {
      tagint tagtmp = atom->tag[i];
      for (int d = 0; d < 3; d++) des[i][d] = bufsorted[tagtmp - 1][d];
    }
  } else {
    int n = atom->nlocal;
    int m = 0;
    for (int i = 0; i < n; i++)
      for (int d = 0; d < 3; d++) {
        des[i][d] = 0.0;
        for (int j = 0; j < np; j++) des[i][d] += src[j][m] * vector[j];
        m++;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::comm_init()
{
  if (np != universe->nworlds)
    error->all(FLERR, "Fix {}: np must equal universe->nworlds", style);

  // Paired ranks require the same number of processes in every bead partition.
  int minprocs, maxprocs;
  MPI_Allreduce(&nprocs, &minprocs, 1, MPI_INT, MPI_MIN, universe->uworld);
  MPI_Allreduce(&nprocs, &maxprocs, 1, MPI_INT, MPI_MAX, universe->uworld);
  if (minprocs != maxprocs)
    error->universe_all(FLERR, "Fix PIMD requires the same number of processors per bead");

  int nlocal = atom->nlocal;
  if (cmode == SINGLE_PROC) {
    memory->destroy(counts);
    memory->destroy(displacements);
    memory->create(counts, nreplica, "FixPIMDNVE:counts");
    memory->create(displacements, nreplica, "FixPIMDNVE:displacements");
    for (int i = 0; i < nreplica; i++) counts[i] = 3*nlocal;
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

void FixPIMDNVE::reallocate_xc()
{
  maxxc = atom->nmax;
  memory->destroy(xc);
  memory->create(xc, maxxc, 3, "FixPIMDNVE:xc");
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::reallocate_x_unwrap()
{
  maxunwrap = atom->nmax;
  memory->destroy(x_unwrap);
  memory->create(x_unwrap, maxunwrap, 3, "FixPIMDNVE:x_unwrap");
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::reallocate()
{
  maxlocal = atom->nmax;
  ntotal = atom->natoms;
  maxsend = maxlocal;
  if (cmode == SINGLE_PROC) {
    memory->destroy(bufsorted);
    memory->destroy(bufsortedall);
    memory->create(bufsorted, ntotal, 3, "FixPIMDNVE:bufsorted");
    memory->create(bufsortedall, nreplica * ntotal, 3, "FixPIMDNVE:bufsortedall");
  } else if (cmode == MULTI_PROC) {
    memory->destroy(bufsend);
    memory->destroy(bufrecv);
    memory->destroy(tagsend);
    memory->destroy(tagrecv);
    memory->destroy(bufbeads);
    memory->create(bufsend, maxlocal*3, "FixPIMDNVE:bufsend");
    memory->create(bufrecv, maxlocal*3, "FixPIMDNVE:bufrecv");
    memory->create(tagsend, maxlocal, "FixPIMDNVE:tagsend");
    memory->create(tagrecv, maxlocal, "FixPIMDNVE:tagrecv");
    memory->create(bufbeads, nreplica, maxlocal * 3, "FixPIMDNVE:bufbeads");
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::inter_replica_comm(double **ptr)
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
      MPI_Sendrecv((void*)&nlocal, 1, MPI_INT,
                  plansend[iplan], TAG_INTER_REPLICA_COUNT,
                  (void*)&nsend, 1, MPI_INT,
                  planrecv[iplan], TAG_INTER_REPLICA_COUNT,
                  universe->uworld, MPI_STATUS_IGNORE);

      // 2) ensure buffers sized for nsend
      if (nsend > maxsend) {
        maxsend = nsend + 200;
        tagsend = (tagint *) memory->srealloc(tagsend, sizeof(tagint) * maxsend,
                                            "FixPIMDNVE:tagsend");
        bufsend = (double *) memory->srealloc(bufsend, sizeof(double) * 3 * maxsend,
                                            "FixPIMDNVE:bufsend");
      }

      // 3) exchange tags:
      //    send my local tags (atom->tag[0..nlocal-1])
      //    receive remote rank's local tags into tagsend[0..nsend-1]
      MPI_Sendrecv((void*)atom->tag, nlocal, MPI_LMP_TAGINT,
                  plansend[iplan], TAG_INTER_REPLICA_TAGS,
                  (void*)tagsend, nsend, MPI_LMP_TAGINT,
                  planrecv[iplan], TAG_INTER_REPLICA_TAGS,
                  universe->uworld, MPI_STATUS_IGNORE);

      // 4) pack ptr for the tags the remote rank needs from me
      //    For each received tag, find my local index and copy ptr[index][0..2]
      std::vector<int> miss_idx;
      std::vector<tagint> miss_tag;
      miss_idx.reserve(nsend);
      miss_tag.reserve(nsend);

      for (int i = 0; i < nsend; i++) {
        const int idx = atom->map(tagsend[i]);
        if (idx >= 0 && idx < nlocal) {
          bufsend[3*i + 0] = ptr[idx][0];
          bufsend[3*i + 1] = ptr[idx][1];
          bufsend[3*i + 2] = ptr[idx][2];
        } else {
          miss_idx.push_back(i);   // remember which slot in bufsend needs collect
          miss_tag.push_back(tagsend[i]);   // remember which tag that slot corresponds to
        }
      }

      // 5) collect missing tags within this world (local-only claiming)
      int missing = !miss_tag.empty();
      int any_missing;
      MPI_Allreduce(&missing, &any_missing, 1, MPI_INT, MPI_MAX, world);
      // Every rank must participate in the ring, including ranks with no requests.
      if (any_missing) {
        std::vector<tagint> rep_tag;
        std::vector<double> rep_val;
        ring_collect(miss_tag, ptr, rep_tag, rep_val);

        // fill missing slots in bufsend by tag lookup (missing is small)
        // Use a simple O(N^2) search since missing tags expected to be few
        for (int k = 0; k < (int)miss_tag.size(); k++) {
          const tagint t = miss_tag[k];
          int pos = -1;
          for (int j = 0; j < (int)rep_tag.size(); j++) {
            if (rep_tag[j] == t) { pos = j; break; }
          }
          if (pos < 0) {
            auto mesg = fmt::format("collect failed: tag {} not returned on world [{}] rank [{}]\n",
                                  (int)t, universe->iworld, comm->me);
            error->universe_one(FLERR, mesg);
          }

          const int i = miss_idx[k];
          bufsend[3*i + 0] = rep_val[3*pos + 0];
          bufsend[3*i + 1] = rep_val[3*pos + 1];
          bufsend[3*i + 2] = rep_val[3*pos + 2];
        }
      }

      // 6) exchange packed x/f buffers:
      //    - send bufsend (3*nsend) to planrecv[iplan]
      //    - receive bufrecv (3*nlocal) from plansend[iplan]
      //
      // This mirrors your reference's direction choices.
      MPI_Sendrecv((void*)bufsend, 3*nsend, MPI_DOUBLE,
                  planrecv[iplan], TAG_INTER_REPLICA_VALS,
                  (void*)bufrecv, 3*nlocal, MPI_DOUBLE,
                  plansend[iplan], TAG_INTER_REPLICA_VALS,
                  universe->uworld, MPI_STATUS_IGNORE);

      // 6) store received x/f for this plan into bufbeads[modeindex[iplan]]
      memcpy(bufbeads[modeindex[iplan]], bufrecv, sizeof(double) * 3 * nlocal);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::ring_collect(const std::vector<tagint> &miss_tag, double **ptr,
                              std::vector<tagint> &rep_tag, std::vector<double> &rep_val)
{
  // ring-collection: collect missing atoms from other ranks in this world
  // by passing missing tag lists and found values in a ring
  const int me = comm->me;
  const int P  = comm->nprocs;
  const int next = (me + 1) % P;
  const int prev = (me - 1 + P) % P;
  const int nlocal = atom->nlocal;

  // Token state for this rank's missing tags
  std::vector<tagint> tok_missing = miss_tag;
  std::vector<tagint> tok_found_tags;
  std::vector<double> tok_found_vals;   // 3 per found tag

  // If missing is tiny as expected, reserve to reduce realloc
  tok_found_tags.reserve(tok_missing.size());
  tok_found_vals.reserve(3 * tok_missing.size());

  // Move one hop per iteration; after P hops, your token returns to you.
  for (int hop = 0; hop < P; hop++) {

    // 1) exchange sizes
    int sm = (int) tok_missing.size();
    int sf = (int) tok_found_tags.size();
    int rm = 0, rf = 0;

    MPI_Sendrecv(&sm, 1, MPI_INT, next, TAG_RING_MISS_COUNT,
                 &rm, 1, MPI_INT, prev, TAG_RING_MISS_COUNT,
                 world, MPI_STATUS_IGNORE);

    MPI_Sendrecv(&sf, 1, MPI_INT, next, TAG_RING_MISS_TAGS,
                 &rf, 1, MPI_INT, prev, TAG_RING_MISS_TAGS,
                 world, MPI_STATUS_IGNORE);

    // 2) prepare recv buffers
    std::vector<tagint> in_missing(rm);
    std::vector<tagint> in_found_tags(rf);
    std::vector<double> in_found_vals(3 * (size_t)rf);

    // 3) exchange payloads
    MPI_Sendrecv(tok_missing.data(), sm, MPI_LMP_TAGINT, next, TAG_RING_REP_COUNT,
                 in_missing.data(), rm, MPI_LMP_TAGINT, prev, TAG_RING_REP_COUNT,
                 world, MPI_STATUS_IGNORE);

    MPI_Sendrecv(tok_found_tags.data(), sf, MPI_LMP_TAGINT, next, TAG_RING_REP_TAGS,
                 in_found_tags.data(), rf, MPI_LMP_TAGINT, prev, TAG_RING_REP_TAGS,
                 world, MPI_STATUS_IGNORE);

    MPI_Sendrecv(tok_found_vals.data(), 3*sf, MPI_DOUBLE, next, TAG_RING_REP_VALS,
                 in_found_vals.data(), 3*rf, MPI_DOUBLE, prev, TAG_RING_REP_VALS,
                 world, MPI_STATUS_IGNORE);

    // 4) process received token: claim only if local owner
    std::vector<tagint> out_missing;
    out_missing.reserve(in_missing.size());

    for (tagint t : in_missing) {
      const int idx = atom->map(t);

      // local-only claim (ignore ghosts)
      // When excecuting this function at the end of initial_integrate,
      // where the coordinates of local atoms are updated while those of ghost atoms are not,
      // considering ghost atoms lead to incorrect coordinates.
      if (idx >= 0 && idx < nlocal) {
        in_found_tags.push_back(t);
        in_found_vals.push_back(ptr[idx][0]);
        in_found_vals.push_back(ptr[idx][1]);
        in_found_vals.push_back(ptr[idx][2]);
      } else {
        out_missing.push_back(t);
      }
    }

    // 5) forward updated token
    tok_missing.swap(out_missing);
    tok_found_tags.swap(in_found_tags);
    tok_found_vals.swap(in_found_vals);
  }

  // After full ring, the token for this rank should be back here.
  // Now we have the resolved list for this rank.
  rep_tag.swap(tok_found_tags);
  rep_val.swap(tok_found_vals);

  // If anything still missing, it's a real error (tag not present in this world).
  if (!tok_missing.empty()) {
    // Print a small sample to help debug
    const tagint t0 = tok_missing[0];
    auto mesg = fmt::format(
        "ring_collect: unresolved {} tags after {} hops on world [{}] rank [{}]. "
        "Example tag = {}.\n",
        (int)tok_missing.size(), P, universe->iworld, me, (int)t0);
    error->universe_one(FLERR, mesg);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::remove_com_motion()
{
  // Do not remove center-of-mass motion for CMD.
  if (method == CMD) return;

  // Cartesian PIMD: every bead; NMPIMD: only the centroid mode.
  if (method == PIMD || universe->iworld == 0) {
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

double FixPIMDNVE::local_kinetic_energy_sum() const
{
  double kine = 0.0;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    for (int j = 0; j < 3; j++) kine += 0.5 * mass[type[i]] * atom->v[i][j] * atom->v[i][j];
  }
  return kine * force->mvv2e;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::reduce_bead_and_total(double local_value, double &bead_value, double &total_value) const
{
  MPI_Allreduce(&local_value, &bead_value, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&bead_value, &total_value, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  total_value /= universe->procs_per_world[universe->iworld];
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::compute_xf_vir()
{
  vir_ = 0.0;
  double xf = 0.0;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    for (int j = 0; j < 3; j++) xf += x_unwrap[i][j] * atom->f[i][j];
  }
  MPI_Allreduce(&xf, &vir_, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::compute_cvir()
{
  centroid_vir = 0.0;
  double xcf = 0.0;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    for (int j = 0; j < 3; j++) xcf += (x_unwrap[i][j] - xc[i][j]) * atom->f[i][j];
  }
  MPI_Allreduce(&xcf, &centroid_vir, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::compute_vir()
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
  double vir_bead = virial[0] + virial[1] + virial[2];
  MPI_Allreduce(&vir_bead, &vir, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  MPI_Allreduce(MPI_IN_PLACE, &virial[0], 6, MPI_DOUBLE, MPI_SUM, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::compute_totke()
{
  totke = ke_bead = 0.0;
  double kine = local_kinetic_energy_sum();
  reduce_bead_and_total(kine, ke_bead, totke);
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::compute_spring_energy()
{
  total_spring_energy = se_bead = 0.0;
  if (method == NMPIMD || method == CMD) {
    spring_energy = 0.0;
    double **x = atom->x;
    double *_mass = atom->mass;
    int *mask = atom->mask;
    int *type = atom->type;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      spring_energy += 0.5 * _mass[type[i]] * fbond * lam[universe->iworld] *
          (x[i][0] * x[i][0] + x[i][1] * x[i][1] + x[i][2] * x[i][2]);
    }
  } else if (method != PIMD) {
    error->universe_all(FLERR, fmt::format("Unknown method parameter for fix {}", style));
  }
  reduce_bead_and_total(spring_energy, se_bead, total_spring_energy);
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::compute_pote()
{
  pe_bead = 0.0;
  pote = 0.0;
  c_pe->compute_scalar();
  pe_bead = c_pe->scalar;
  double pot_energy_partition = pe_bead / universe->procs_per_world[universe->iworld];
  MPI_Allreduce(&pot_energy_partition, &pote, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::compute_tote()
{
  tote = totke + pote + total_spring_energy;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::compute_t_prim()
{
  t_prim = 1.5 * group->count(igroup) * np * force->boltz * temp -
      total_spring_energy * inverse_np;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::compute_t_vir()
{
  t_vir = -0.5 * inverse_np * vir_;
  t_cv = 1.5 * group->count(igroup) * force->boltz * temp - 0.5 * inverse_np * centroid_vir;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::compute_p_prim()
{
  double inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
  p_prim = static_cast<double>(group->count(igroup)) * np * force->boltz * temp * inv_volume -
      (2.0 / 3.0) * inv_volume * total_spring_energy;
  p_prim *= force->nktv2p;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::compute_p_cv()
{
  double inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
  p_md = (1.0 / 3.0) * inv_volume * (totke + vir);
  if ((method == NMPIMD || method == CMD) && universe->iworld == 0) {
    p_cv = (1.0 / 3.0) * inv_volume * ((2.0 * ke_bead - centroid_vir) * force->nktv2p + vir) / np;
  } else if (method == PIMD) {
    p_cv = (1.0 / 3.0) * inv_volume * ((2.0 * totke / np - centroid_vir) * force->nktv2p + vir) / np;
  }
  MPI_Bcast(&p_cv, 1, MPI_DOUBLE, 0, universe->uworld);
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::write_restart(FILE *fp)
{
  int nsize = base_restart_size() + subclass_restart_size();
  double *list;
  memory->create(list, nsize, "FixPIMDNVE:list");
  int n = pack_base_restart(list);
  pack_subclass_restart(list, n);
  if (comm->me == 0) {
    int size = nsize * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    if (nsize) fwrite(list, sizeof(double), nsize, fp);
  }
  memory->destroy(list);
}

/* ---------------------------------------------------------------------- */

int FixPIMDNVE::base_restart_size() const
{
  return 0;
}

/* ---------------------------------------------------------------------- */

int FixPIMDNVE::pack_base_restart(double *) const
{
  return 0;
}

/* ---------------------------------------------------------------------- */

int FixPIMDNVE::unpack_base_restart(const double *)
{
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixPIMDNVE::restart(char *buf)
{
  auto *list = (double *) buf;
  unpack_subclass_restart(list, unpack_base_restart(list));
}

/* ---------------------------------------------------------------------- */

double FixPIMDNVE::compute_vector(int n)
{
  const int prefix = nuclear_vector_size();
  if (n < prefix) return compute_nuclear_vector(n);
  return compute_subclass_vector(n - prefix);
}

/* ---------------------------------------------------------------------- */

int FixPIMDNVE::nuclear_vector_size() const
{
  return 10;
}

/* ---------------------------------------------------------------------- */

double FixPIMDNVE::compute_nuclear_vector(int n) const
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
  return 0.0;
}
