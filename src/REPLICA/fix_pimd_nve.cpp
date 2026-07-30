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

enum { PHYSICAL, NORMAL };
enum { BAOAB, OBABO };
enum { SINGLE_PROC, MULTI_PROC };

void FixPIMDNVE::init_defaults()
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
}

bool FixPIMDNVE::parse_common_keyword(int narg, char **arg, int &i)
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
    if (fmass < 0.0 || fmass > np) error->all(FLERR, "Invalid fmass value for fix {}", style);
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

void FixPIMDNVE::parse_arguments(int narg, char **arg, const KeywordParser &subclass_parser)
{
  for (int i = 3; i < narg;) {
    if (parse_common_keyword(narg, arg, i)) continue;
    if (subclass_parser && subclass_parser(narg, arg, i)) continue;
    error->all(FLERR, "Unknown keyword {} for fix {}", arg[i], style);
  }
}

FixPIMDNVE::FixPIMDNVE(LAMMPS *lmp, int narg, char **arg, bool) :
    Fix(lmp, narg, arg), mass(nullptr), rootworld(MPI_COMM_NULL), plansend(nullptr),
    planrecv(nullptr), tagsend(nullptr), tagrecv(nullptr), bufsend(nullptr), bufrecv(nullptr),
    bufbeads(nullptr), bufsorted(nullptr), bufsortedall(nullptr), tagsendall(nullptr),
    tagrecvall(nullptr), bufsendall(nullptr), bufrecvall(nullptr), counts(nullptr),
    displacements(nullptr), lam(nullptr), M_x2xp(nullptr), M_xp2x(nullptr), modeindex(nullptr),
    _omega_k(nullptr), Lan_s(nullptr), Lan_c(nullptr), xc(nullptr), xcall(nullptr),
    x_unwrap(nullptr), id_pe(nullptr), id_press(nullptr), c_pe(nullptr), c_press(nullptr)
{
  init_defaults();

  if (domain->dimension != 3)
    error->universe_all(FLERR, fmt::format("Fix {} requires a 3d system", style));
  if (narg < 3) utils::missing_cmd_args(FLERR, std::string("fix ") + style, error);
}

FixPIMDNVE::FixPIMDNVE(LAMMPS *lmp, int narg, char **arg) : FixPIMDNVE(lmp, narg, arg, true)
{
  parse_arguments(narg, arg, {});
  if (method == CMD) error->all(FLERR, "Fix pimd/nve does not support method cmd");
  finish_constructor_setup();
}

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
  if (rootworld != MPI_COMM_NULL) MPI_Comm_free(&rootworld);
}

int FixPIMDNVE::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

void FixPIMDNVE::init()
{
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

  omega_np = np / (hbar * beta) * sqrt(force->mvv2e);
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

  setup_subclass_state();
  t_prim = t_vir = t_cv = p_cv = p_md = 0.0;
}

void FixPIMDNVE::setup(int vflag)
{
  if (method == NMPIMD || method == CMD) {
    begin_normal_mode_coordinate_propagation();
  } else if (method == PIMD) {
    unmap_coordinates(atom->x, atom->image);
    inter_replica_comm(atom->x);
    spring_force();
  } else {
    error->universe_all(FLERR, fmt::format("Unknown method parameter for fix {}", style));
  }
  after_force_transform_hook();
  collect_xc();
  compute_spring_energy();
  compute_t_prim();
  compute_p_prim();
  if (method == NMPIMD || method == CMD)
    finalize_setup_normal_mode_coordinates();
  else
    remap_coordinates(atom->x, atom->image);

  post_force(vflag);
  compute_totke();
  end_of_step();
  schedule_common_computes();
}

void FixPIMDNVE::initial_integrate(int /*vflag*/)
{
  b_step();
  if (method == NMPIMD || method == CMD) {
    begin_normal_mode_coordinate_propagation();
    propagate_normal_mode_coordinate_halfstep();
    propagate_normal_mode_coordinate_halfstep();
    finalize_normal_mode_coordinate_propagation();
  } else if (method == PIMD) {
    unmap_coordinates(atom->x, atom->image);
    q_step();
    q_step();
    collect_xc();
    remap_coordinates(atom->x, atom->image);
  } else {
    error->universe_all(FLERR, fmt::format("Unknown method parameter for fix {}", style));
  }
}

void FixPIMDNVE::final_integrate()
{
  b_step();
}

void FixPIMDNVE::post_force(int /*flag*/)
{
  prepare_common_virial_state();
  compute_vir();
  compute_xf_vir();
  compute_cvir();
  compute_t_vir();

  if (method == PIMD) {
    unmap_coordinates(atom->x, atom->image);
    inter_replica_comm(atom->x);
    spring_force();
    compute_spring_energy();
    compute_t_prim();
    remap_coordinates(atom->x, atom->image);
  }

  compute_pote();
  if (method == NMPIMD || method == CMD) prepare_normal_mode_forces();
  after_force_transform_hook();

  schedule_common_computes();
}

void FixPIMDNVE::end_of_step()
{
  compute_totke();
  compute_p_cv();
  compute_tote();
}

void FixPIMDNVE::setup_subclass_state() {}

void FixPIMDNVE::after_force_transform_hook() {}

void FixPIMDNVE::unmap_coordinates(double **coords, imageint *image)
{
  if (!mapflag) return;

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) domain->unmap(coords[i], image[i]);
}

void FixPIMDNVE::remap_coordinates(double **coords, imageint *image)
{
  if (!mapflag) return;

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) domain->unmap_inv(coords[i], image[i]);
}

double **FixPIMDNVE::normal_mode_transform_buffer()
{
  if (cmode == SINGLE_PROC) return bufsortedall;
  return bufbeads;
}

void FixPIMDNVE::forward_normal_mode_transform(double **ptr)
{
  inter_replica_comm(ptr);
  nmpimd_transform(normal_mode_transform_buffer(), ptr, M_x2xp[universe->iworld]);
}

void FixPIMDNVE::backward_normal_mode_transform(double **ptr)
{
  inter_replica_comm(ptr);
  nmpimd_transform(normal_mode_transform_buffer(), ptr, M_xp2x[universe->iworld]);
}

void FixPIMDNVE::finalize_setup_normal_mode_coordinates()
{
  backward_normal_mode_transform(atom->x);
  remap_coordinates(atom->x, atom->image);
}

void FixPIMDNVE::begin_normal_mode_coordinate_propagation()
{
  unmap_coordinates(atom->x, atom->image);
  forward_normal_mode_transform(atom->x);
}

void FixPIMDNVE::propagate_normal_mode_coordinate_halfstep()
{
  qc_step();
  a_step();
}

void FixPIMDNVE::finalize_normal_mode_coordinate_propagation()
{
  collect_xc();
  compute_spring_energy();
  compute_t_prim();
  compute_p_prim();
  backward_normal_mode_transform(atom->x);
  remap_coordinates(atom->x, atom->image);
}

void FixPIMDNVE::prepare_common_virial_state()
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
}

void FixPIMDNVE::prepare_normal_mode_forces()
{
  forward_normal_mode_transform(atom->f);
}

void FixPIMDNVE::schedule_common_computes()
{
  c_pe->addstep(update->ntimestep + 1);
  c_press->addstep(update->ntimestep + 1);
}

int FixPIMDNVE::subclass_restart_size() const
{
  return 0;
}

int FixPIMDNVE::pack_subclass_restart(double *, int n) const
{
  return n;
}

int FixPIMDNVE::unpack_subclass_restart(const double *, int n)
{
  return n;
}

int FixPIMDNVE::subclass_vector_size() const
{
  return 0;
}

double FixPIMDNVE::compute_subclass_vector(int) const
{
  return 0.0;
}

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

void FixPIMDNVE::b_step()
{
  apply_force_velocity_kick();
}

void FixPIMDNVE::apply_force_velocity_kick()
{
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

  double omega_np_dt_half = omega_np * update->dt * 0.5;
  if (fmmode == PHYSICAL) {
    for (int i = 0; i < np; i++) {
      _omega_k[i] = omega_np * sqrt(lam[i]) / sqrt(fmass);
      Lan_c[i] = cos(sqrt(lam[i]) * omega_np_dt_half);
      Lan_s[i] = sin(sqrt(lam[i]) * omega_np_dt_half);
    }
  } else if (fmmode == NORMAL) {
    for (int i = 0; i < np; i++) {
      _omega_k[i] = omega_np / sqrt(fmass);
      Lan_c[i] = cos(omega_np_dt_half);
      Lan_s[i] = sin(omega_np_dt_half);
    }
  } else {
    error->universe_all(FLERR, "Unknown fmmode setting; only physical and normal are supported!");
  }

  int iworld = universe->iworld;
  for (int i = 1; i <= atom->ntypes; i++) {
    mass[i] = atom->mass[i] * fmass;
    if (iworld != 0 && fmmode == NORMAL) mass[i] *= lam[iworld];
  }
}

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

void FixPIMDNVE::comm_init()
{
  if (np != universe->nworlds)
    error->all(FLERR, "Fix {}: np must equal universe->nworlds", style);

  int nlocal = atom->nlocal;
  if (cmode == SINGLE_PROC) {
    memory->destroy(counts);
    memory->destroy(displacements);
    memory->create(counts, nreplica, "FixPIMDNVE:counts");
    memory->create(displacements, nreplica, "FixPIMDNVE:displacements");
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
    int isend = ireplica + i + 1;
    if (isend >= nreplica) isend -= nreplica;

    int irecv = ireplica - (i + 1);
    if (irecv < 0) irecv += nreplica;

    plansend[i] = universe->root_proc[isend];
    planrecv[i] = universe->root_proc[irecv];
    modeindex[i] = irecv;
  }

  x_next = (universe->iworld + 1 + universe->nworlds) % universe->nworlds;
  x_last = (universe->iworld - 1 + universe->nworlds) % universe->nworlds;
}

void FixPIMDNVE::reallocate_xc()
{
  maxxc = atom->nmax;
  memory->destroy(xc);
  memory->create(xc, maxxc, 3, "FixPIMDNVE:xc");
}

void FixPIMDNVE::reallocate_x_unwrap()
{
  maxunwrap = atom->nmax;
  memory->destroy(x_unwrap);
  memory->create(x_unwrap, maxunwrap, 3, "FixPIMDNVE:x_unwrap");
}

void FixPIMDNVE::reallocate()
{
  maxlocal = atom->nmax;
  ntotal = atom->natoms;
  if (cmode == SINGLE_PROC) {
    memory->destroy(bufsorted);
    memory->destroy(bufsortedall);
    memory->create(bufsorted, ntotal, 3, "FixPIMDNVE:bufsorted");
    memory->create(bufsortedall, nreplica * ntotal, 3, "FixPIMDNVE:bufsortedall");
  } else {
    memory->destroy(bufsend);
    memory->destroy(bufrecv);
    memory->destroy(tagsend);
    memory->destroy(tagrecv);
    memory->destroy(bufbeads);
    memory->create(bufsend, maxlocal, 3, "FixPIMDNVE:bufsend");
    memory->create(bufrecv, maxlocal, 3, "FixPIMDNVE:bufrecv");
    memory->create(tagsend, maxlocal, "FixPIMDNVE:tagsend");
    memory->create(tagrecv, maxlocal, "FixPIMDNVE:tagrecv");
    memory->create(bufbeads, nreplica, maxlocal * 3, "FixPIMDNVE:bufbeads");

    memory->destroy(tagsendall);
    memory->destroy(tagrecvall);
    memory->destroy(bufsendall);
    memory->destroy(bufrecvall);
    memory->destroy(counts);
    memory->destroy(displacements);
    memory->create(tagsendall, ntotal, "FixPIMDNVE:tagsendall");
    memory->create(tagrecvall, ntotal, "FixPIMDNVE:tagrecvall");
    memory->create(bufsendall, ntotal, 3, "FixPIMDNVE:bufsendall");
    memory->create(bufrecvall, ntotal, 3, "FixPIMDNVE:bufrecvall");
    memory->create(counts, nprocs, "FixPIMDNVE:counts_multi");
    memory->create(displacements, nprocs, "FixPIMDNVE:displacements_multi");
  }
}

void FixPIMDNVE::inter_replica_comm(double **ptr)
{
  MPI_Request requests[2];
  MPI_Status statuses[2];
  if (atom->nmax > maxlocal) reallocate();
  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;

  if (cmode == SINGLE_PROC) {
    int m = 0;
    for (int i = 0; i < nlocal; i++) {
      tagint tagtmp = tag[i];
      bufsorted[tagtmp - 1][0] = ptr[i][0];
      bufsorted[tagtmp - 1][1] = ptr[i][1];
      bufsorted[tagtmp - 1][2] = ptr[i][2];
      m++;
    }
    MPI_Allgatherv(bufsorted[0], 3 * m, MPI_DOUBLE, bufsortedall[0], counts, displacements,
                   MPI_DOUBLE, universe->uworld);
  } else {
    for (int i = 0; i < nlocal; i++) {
      bufbeads[ireplica][3 * i + 0] = ptr[i][0];
      bufbeads[ireplica][3 * i + 1] = ptr[i][1];
      bufbeads[ireplica][3 * i + 2] = ptr[i][2];
    }
    int m = 0;
    for (int i = 0; i < nlocal; i++) {
      tagsend[m] = tag[i];
      bufsend[m][0] = ptr[i][0];
      bufsend[m][1] = ptr[i][1];
      bufsend[m][2] = ptr[i][2];
      m++;
    }
    MPI_Gather(&m, 1, MPI_INT, counts, 1, MPI_INT, 0, world);
    displacements[0] = 0;
    for (int i = 0; i < nprocs - 1; i++) displacements[i + 1] = displacements[i] + counts[i];
    MPI_Gatherv(tagsend, m, MPI_LMP_TAGINT, tagsendall, counts, displacements, MPI_LMP_TAGINT, 0,
                world);
    for (int i = 0; i < nprocs; i++) counts[i] *= 3;
    for (int i = 0; i < nprocs - 1; i++) displacements[i + 1] = displacements[i] + counts[i];
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
      for (int i = 0; i < ntotal; i++) {
        int mapped = atom->map(tagrecvall[i]);
        if (mapped < 0 || mapped >= nlocal) continue;
        bufbeads[modeindex[iplan]][3 * mapped + 0] = bufrecvall[i][0];
        bufbeads[modeindex[iplan]][3 * mapped + 1] = bufrecvall[i][1];
        bufbeads[modeindex[iplan]][3 * mapped + 2] = bufrecvall[i][2];
      }
    }
  }
}

void FixPIMDNVE::remove_com_motion()
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

double FixPIMDNVE::estimator_atom_count(bool restrict_group) const
{
  return restrict_group ? static_cast<double>(group->count(igroup)) : static_cast<double>(atom->natoms);
}

double FixPIMDNVE::local_kinetic_energy_sum(bool restrict_group) const
{
  double kine = 0.0;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  for (int i = 0; i < nlocal; i++) {
    if (restrict_group && !(mask[i] & groupbit)) continue;
    for (int j = 0; j < 3; j++) kine += 0.5 * mass[type[i]] * atom->v[i][j] * atom->v[i][j];
  }
  return kine * force->mvv2e;
}

double FixPIMDNVE::local_normal_mode_spring_energy_sum(bool restrict_group) const
{
  double energy = 0.0;
  double **x = atom->x;
  double *_mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (restrict_group && !(mask[i] & groupbit)) continue;
    energy += 0.5 * _mass[type[i]] * fbond * lam[universe->iworld] *
        (x[i][0] * x[i][0] + x[i][1] * x[i][1] + x[i][2] * x[i][2]);
  }
  return energy;
}

double FixPIMDNVE::local_xf_virial_sum(bool restrict_group) const
{
  double xf = 0.0;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  for (int i = 0; i < nlocal; i++) {
    if (restrict_group && !(mask[i] & groupbit)) continue;
    for (int j = 0; j < 3; j++) xf += x_unwrap[i][j] * atom->f[i][j];
  }
  return xf;
}

double FixPIMDNVE::local_centroid_virial_sum(bool restrict_group) const
{
  double xcf = 0.0;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  for (int i = 0; i < nlocal; i++) {
    if (restrict_group && !(mask[i] & groupbit)) continue;
    for (int j = 0; j < 3; j++) xcf += (x_unwrap[i][j] - xc[i][j]) * atom->f[i][j];
  }
  return xcf;
}

void FixPIMDNVE::reduce_bead_and_total(double local_value, double &bead_value, double &total_value) const
{
  MPI_Allreduce(&local_value, &bead_value, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&bead_value, &total_value, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  total_value /= universe->procs_per_world[universe->iworld];
}

double FixPIMDNVE::reduce_partition_scalar(double partition_scalar) const
{
  double total_scalar = 0.0;
  MPI_Allreduce(&partition_scalar, &total_scalar, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  return total_scalar;
}

void FixPIMDNVE::compute_xf_vir()
{
  vir_ = 0.0;
  double xf = local_xf_virial_sum(false);
  MPI_Allreduce(&xf, &vir_, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
}

void FixPIMDNVE::compute_cvir()
{
  centroid_vir = 0.0;
  double xcf = local_centroid_virial_sum(false);
  MPI_Allreduce(&xcf, &centroid_vir, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
}

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

void FixPIMDNVE::compute_totke()
{
  totke = ke_bead = 0.0;
  double kine = local_kinetic_energy_sum(false);
  reduce_bead_and_total(kine, ke_bead, totke);
}

void FixPIMDNVE::compute_spring_energy()
{
  total_spring_energy = se_bead = 0.0;
  if (method == NMPIMD || method == CMD) {
    spring_energy = local_normal_mode_spring_energy_sum(false);
  } else if (method != PIMD) {
    error->universe_all(FLERR, fmt::format("Unknown method parameter for fix {}", style));
  }
  reduce_bead_and_total(spring_energy, se_bead, total_spring_energy);
}

void FixPIMDNVE::compute_pote()
{
  pe_bead = 0.0;
  pote = 0.0;
  c_pe->compute_scalar();
  pe_bead = c_pe->scalar;
  double pot_energy_partition = pe_bead / universe->procs_per_world[universe->iworld];
  pote = reduce_partition_scalar(pot_energy_partition);
}

void FixPIMDNVE::compute_tote()
{
  tote = totke + pote + total_spring_energy;
}

void FixPIMDNVE::compute_t_prim()
{
  t_prim = 1.5 * estimator_atom_count(false) * np * force->boltz * temp -
      total_spring_energy * inverse_np;
}

void FixPIMDNVE::compute_t_vir()
{
  t_vir = -0.5 * inverse_np * vir_;
  t_cv = 1.5 * estimator_atom_count(false) * force->boltz * temp - 0.5 * inverse_np * centroid_vir;
}

void FixPIMDNVE::compute_p_prim()
{
  double inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
  p_prim = estimator_atom_count(false) * np * force->boltz * temp * inv_volume -
      (2.0 / 3.0) * inv_volume * total_spring_energy;
  p_prim *= force->nktv2p;
}

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

void FixPIMDNVE::write_restart(FILE *fp)
{
  int nsize = size_restart_global();
  double *list;
  memory->create(list, nsize, "FixPIMDNVE:list");
  pack_restart_data(list);
  if (comm->me == 0) {
    int size = nsize * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    if (nsize) fwrite(list, sizeof(double), nsize, fp);
  }
  memory->destroy(list);
}

int FixPIMDNVE::size_restart_global()
{
  return base_restart_size() + subclass_restart_size();
}

int FixPIMDNVE::pack_restart_data(double *list)
{
  int n = pack_base_restart(list);
  return pack_subclass_restart(list, n);
}

int FixPIMDNVE::base_restart_size() const
{
  return 0;
}

int FixPIMDNVE::pack_base_restart(double *) const
{
  return 0;
}

int FixPIMDNVE::unpack_base_restart(const double *)
{
  return 0;
}

void FixPIMDNVE::restart(char *buf)
{
  auto *list = (double *) buf;
  unpack_subclass_restart(list, unpack_base_restart(list));
}

double FixPIMDNVE::compute_vector(int n)
{
  const int prefix = nuclear_vector_size();
  if (n < prefix) return compute_nuclear_vector(n);
  return compute_subclass_vector(n - prefix);
}

int FixPIMDNVE::nuclear_vector_size() const
{
  return 10;
}

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
