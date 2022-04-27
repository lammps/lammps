/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Package      FixPIMD
   Purpose      Quantum Path Integral Algorithm for Quantum Chemistry
   Copyright    Voth Group @ University of Chicago
   Authors      Chris Knight & Yuxing Peng (yuxing at uchicago.edu)

   Updated      Oct-01-2011
   Version      1.0
------------------------------------------------------------------------- */

#include "fix_pimd.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum { PIMD, NMPIMD, CMD };

/* ---------------------------------------------------------------------- */

FixPIMD::FixPIMD(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  max_nsend = 0;
  tag_send = nullptr;
  buf_send = nullptr;

  max_nlocal = 0;
  buf_recv = nullptr;
  buf_beads = nullptr;

  size_plan = 0;
  plan_send = plan_recv = nullptr;

  M_x2xp = M_xp2x = M_f2fp = M_fp2f = nullptr;
  lam = nullptr;
  mode_index = nullptr;

  mass = nullptr;

  array_atom = nullptr;
  nhc_eta = nullptr;
  nhc_eta_dot = nullptr;
  nhc_eta_dotdot = nullptr;
  nhc_eta_mass = nullptr;

  method = PIMD;
  fmass = 1.0;
  nhc_temp = 298.15;
  nhc_nchain = 2;
  sp = 1.0;

  for (int i = 3; i < narg - 1; i += 2) {
    if (strcmp(arg[i], "method") == 0) {
      if (strcmp(arg[i + 1], "pimd") == 0)
        method = PIMD;
      else if (strcmp(arg[i + 1], "nmpimd") == 0)
        method = NMPIMD;
      else if (strcmp(arg[i + 1], "cmd") == 0)
        method = CMD;
      else
        error->universe_all(FLERR, "Unknown method parameter for fix pimd");
    } else if (strcmp(arg[i], "fmass") == 0) {
      fmass = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (fmass < 0.0 || fmass > 1.0)
        error->universe_all(FLERR, "Invalid fmass value for fix pimd");
    } else if (strcmp(arg[i], "sp") == 0) {
      sp = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (fmass < 0.0) error->universe_all(FLERR, "Invalid sp value for fix pimd");
    } else if (strcmp(arg[i], "temp") == 0) {
      nhc_temp = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (nhc_temp < 0.0) error->universe_all(FLERR, "Invalid temp value for fix pimd");
    } else if (strcmp(arg[i], "nhc") == 0) {
      nhc_nchain = utils::inumeric(FLERR, arg[i + 1], false, lmp);
      if (nhc_nchain < 2) error->universe_all(FLERR, "Invalid nhc value for fix pimd");
    } else
      error->universe_all(FLERR, fmt::format("Unknown keyword {} for fix pimd", arg[i]));
  }

  /* Initiation */

  size_peratom_cols = 12 * nhc_nchain + 3;

  nhc_offset_one_1 = 3 * nhc_nchain;
  nhc_offset_one_2 = 3 * nhc_nchain + 3;
  nhc_size_one_1 = sizeof(double) * nhc_offset_one_1;
  nhc_size_one_2 = sizeof(double) * nhc_offset_one_2;

  restart_peratom = 1;
  peratom_flag = 1;
  peratom_freq = 1;

  global_freq = 1;
  vector_flag = 1;
  size_vector = 2;
  extvector = 1;
  comm_forward = 3;

  atom->add_callback(Atom::GROW);       // Call LAMMPS to allocate memory for per-atom array
  atom->add_callback(Atom::RESTART);    // Call LAMMPS to re-assign restart-data for per-atom array

  grow_arrays(atom->nmax);

  // some initilizations

  nhc_ready = false;
}

/* ---------------------------------------------------------------------- */
FixPIMD::~FixPIMD()
{
  delete[] mass;
  atom->delete_callback(id, Atom::GROW);
  atom->delete_callback(id, Atom::RESTART);

  memory->destroy(M_x2xp);
  memory->destroy(M_xp2x);
  memory->destroy(M_f2fp);
  memory->destroy(M_fp2f);
  memory->sfree(lam);

  if (buf_beads)
    for (int i = 0; i < np; i++) memory->sfree(buf_beads[i]);
  delete[] buf_beads;
  delete[] plan_send;
  delete[] plan_recv;
  delete[] mode_index;

  memory->sfree(tag_send);
  memory->sfree(buf_send);
  memory->sfree(buf_recv);

  memory->destroy(array_atom);
  memory->destroy(nhc_eta);
  memory->destroy(nhc_eta_dot);
  memory->destroy(nhc_eta_dotdot);
  memory->destroy(nhc_eta_mass);
}

/* ---------------------------------------------------------------------- */
int FixPIMD::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::init()
{
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, "Fix pimd requires an atom map, see atom_modify");

  if (universe->me == 0 && universe->uscreen)
    fprintf(universe->uscreen, "Fix pimd initializing Path-Integral ...\n");

  // prepare the constants

  np = universe->nworlds;
  inverse_np = 1.0 / np;

  /* The first solution for the force constant, using SI units

  const double Boltzmann = 1.3806488E-23;    // SI unit: J/K
  const double Plank     = 6.6260755E-34;    // SI unit: m^2 kg / s

  double hbar = Plank / ( 2.0 * MY_PI ) * sp;
  double beta = 1.0 / ( Boltzmann * input.nh_temp);

  // - P / ( beta^2 * hbar^2)   SI unit: s^-2
  double _fbond = -1.0 / (beta*beta*hbar*hbar) * input.nbeads;

  // convert the units: s^-2 -> (kcal/mol) / (g/mol) / (A^2)
  fbond = _fbond * 4.184E+26;

  */

  /* The current solution, using LAMMPS internal real units */

  const double Boltzmann = force->boltz;
  const double Plank = force->hplanck;

  double hbar = Plank / (2.0 * MY_PI);
  double beta = 1.0 / (Boltzmann * nhc_temp);
  double _fbond = 1.0 * np / (beta * beta * hbar * hbar);

  omega_np = sqrt(np) / (hbar * beta) * sqrt(force->mvv2e);
  fbond = -_fbond * force->mvv2e;

  if (universe->me == 0)
    printf("Fix pimd -P/(beta^2 * hbar^2) = %20.7lE (kcal/mol/A^2)\n\n", fbond);

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  comm_init();

  mass = new double[atom->ntypes + 1];

  if (method == CMD || method == NMPIMD)
    nmpimd_init();
  else
    for (int i = 1; i <= atom->ntypes; i++) mass[i] = atom->mass[i] / np * fmass;

  if (!nhc_ready) nhc_init();
}

/* ---------------------------------------------------------------------- */

void FixPIMD::setup(int vflag)
{
  if (universe->me == 0 && universe->uscreen)
    fprintf(universe->uscreen, "Setting up Path-Integral ...\n");

  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPIMD::initial_integrate(int /*vflag*/)
{
  nhc_update_v();
  nhc_update_x();
}

/* ---------------------------------------------------------------------- */

void FixPIMD::final_integrate()
{
  nhc_update_v();
}

/* ---------------------------------------------------------------------- */

void FixPIMD::post_force(int /*flag*/)
{
  for (int i = 0; i < atom->nlocal; i++)
    for (int j = 0; j < 3; j++) atom->f[i][j] /= np;

  comm_exec(atom->x);
  spring_force();

  if (method == CMD || method == NMPIMD) {
    /* forward comm for the force on ghost atoms */

    nmpimd_fill(atom->f);

    /* inter-partition comm */

    comm_exec(atom->f);

    /* normal-mode transform */

    nmpimd_transform(buf_beads, atom->f, M_f2fp[universe->iworld]);
  }
}

/* ----------------------------------------------------------------------
   Nose-Hoover Chains
------------------------------------------------------------------------- */

void FixPIMD::nhc_init()
{
  double tau = 1.0 / omega_np;
  double KT = force->boltz * nhc_temp;

  double mass0 = KT * tau * tau;
  int max = 3 * atom->nlocal;

  for (int i = 0; i < max; i++) {
    for (int ichain = 0; ichain < nhc_nchain; ichain++) {
      nhc_eta[i][ichain] = 0.0;
      nhc_eta_dot[i][ichain] = 0.0;
      nhc_eta_dot[i][ichain] = 0.0;
      nhc_eta_dotdot[i][ichain] = 0.0;
      nhc_eta_mass[i][ichain] = mass0;
      if ((method == CMD || method == NMPIMD) && universe->iworld == 0)
        ;
      else
        nhc_eta_mass[i][ichain] *= fmass;
    }

    nhc_eta_dot[i][nhc_nchain] = 0.0;

    for (int ichain = 1; ichain < nhc_nchain; ichain++)
      nhc_eta_dotdot[i][ichain] = (nhc_eta_mass[i][ichain - 1] * nhc_eta_dot[i][ichain - 1] *
                                       nhc_eta_dot[i][ichain - 1] * force->mvv2e -
                                   KT) /
          nhc_eta_mass[i][ichain];
  }

  // Zero NH acceleration for CMD

  if (method == CMD && universe->iworld == 0)
    for (int i = 0; i < max; i++)
      for (int ichain = 0; ichain < nhc_nchain; ichain++) nhc_eta_dotdot[i][ichain] = 0.0;

  nhc_ready = true;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::nhc_update_x()
{
  int n = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;

  if (method == CMD || method == NMPIMD) {
    nmpimd_fill(atom->v);
    comm_exec(atom->v);

    /* borrow the space of atom->f to store v in cartisian */

    v = atom->f;
    nmpimd_transform(buf_beads, v, M_xp2x[universe->iworld]);
  }

  for (int i = 0; i < n; i++) {
    x[i][0] += dtv * v[i][0];
    x[i][1] += dtv * v[i][1];
    x[i][2] += dtv * v[i][2];
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::nhc_update_v()
{
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

  t_sys = 0.0;
  if (method == CMD && universe->iworld == 0) return;

  double expfac;
  int nmax = 3 * atom->nlocal;
  double KT = force->boltz * nhc_temp;
  double kecurrent, t_current;

  double dthalf = 0.5 * update->dt;
  double dt4 = 0.25 * update->dt;
  double dt8 = 0.125 * update->dt;

  for (int i = 0; i < nmax; i++) {
    int iatm = i / 3;
    int idim = i % 3;

    double *vv = v[iatm];

    kecurrent = mass[type[iatm]] * vv[idim] * vv[idim] * force->mvv2e;
    t_current = kecurrent / force->boltz;

    double *eta = nhc_eta[i];
    double *eta_dot = nhc_eta_dot[i];
    double *eta_dotdot = nhc_eta_dotdot[i];

    eta_dotdot[0] = (kecurrent - KT) / nhc_eta_mass[i][0];

    for (int ichain = nhc_nchain - 1; ichain > 0; ichain--) {
      expfac = exp(-dt8 * eta_dot[ichain + 1]);
      eta_dot[ichain] *= expfac;
      eta_dot[ichain] += eta_dotdot[ichain] * dt4;
      eta_dot[ichain] *= expfac;
    }

    expfac = exp(-dt8 * eta_dot[1]);
    eta_dot[0] *= expfac;
    eta_dot[0] += eta_dotdot[0] * dt4;
    eta_dot[0] *= expfac;

    // Update particle velocities half-step

    double factor_eta = exp(-dthalf * eta_dot[0]);
    vv[idim] *= factor_eta;

    t_current *= (factor_eta * factor_eta);
    kecurrent = force->boltz * t_current;
    eta_dotdot[0] = (kecurrent - KT) / nhc_eta_mass[i][0];

    for (int ichain = 0; ichain < nhc_nchain; ichain++) eta[ichain] += dthalf * eta_dot[ichain];

    eta_dot[0] *= expfac;
    eta_dot[0] += eta_dotdot[0] * dt4;
    eta_dot[0] *= expfac;

    for (int ichain = 1; ichain < nhc_nchain; ichain++) {
      expfac = exp(-dt8 * eta_dot[ichain + 1]);
      eta_dot[ichain] *= expfac;
      eta_dotdot[ichain] =
          (nhc_eta_mass[i][ichain - 1] * eta_dot[ichain - 1] * eta_dot[ichain - 1] - KT) /
          nhc_eta_mass[i][ichain];
      eta_dot[ichain] += eta_dotdot[ichain] * dt4;
      eta_dot[ichain] *= expfac;
    }

    t_sys += t_current;
  }

  t_sys /= nmax;
}

/* ----------------------------------------------------------------------
   Normal Mode PIMD
------------------------------------------------------------------------- */

void FixPIMD::nmpimd_init()
{
  memory->create(M_x2xp, np, np, "fix_feynman:M_x2xp");
  memory->create(M_xp2x, np, np, "fix_feynman:M_xp2x");
  memory->create(M_f2fp, np, np, "fix_feynman:M_f2fp");
  memory->create(M_fp2f, np, np, "fix_feynman:M_fp2f");

  lam = (double *) memory->smalloc(sizeof(double) * np, "FixPIMD::lam");

  // Set up  eigenvalues

  lam[0] = 0.0;
  if (np % 2 == 0) lam[np - 1] = 4.0 * np;

  for (int i = 2; i <= np / 2; i++) {
    lam[2 * i - 3] = lam[2 * i - 2] = 2.0 * np * (1.0 - 1.0 * cos(2.0 * MY_PI * (i - 1) / np));
  }

  // Set up eigenvectors for non-degenerated modes

  for (int i = 0; i < np; i++) {
    M_x2xp[0][i] = 1.0 / np;
    if (np % 2 == 0) M_x2xp[np - 1][i] = 1.0 / np * pow(-1.0, i);
  }

  // Set up eigenvectors for degenerated modes

  for (int i = 0; i < (np - 1) / 2; i++)
    for (int j = 0; j < np; j++) {
      M_x2xp[2 * i + 1][j] = sqrt(2.0) * cos(2.0 * MY_PI * (i + 1) * j / np) / np;
      M_x2xp[2 * i + 2][j] = -sqrt(2.0) * sin(2.0 * MY_PI * (i + 1) * j / np) / np;
    }

  // Set up Ut

  for (int i = 0; i < np; i++)
    for (int j = 0; j < np; j++) {
      M_xp2x[i][j] = M_x2xp[j][i] * np;
      M_f2fp[i][j] = M_x2xp[i][j] * np;
      M_fp2f[i][j] = M_xp2x[i][j];
    }

  // Set up masses

  int iworld = universe->iworld;

  for (int i = 1; i <= atom->ntypes; i++) {
    mass[i] = atom->mass[i];

    if (iworld) {
      mass[i] *= lam[iworld];
      mass[i] *= fmass;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::nmpimd_fill(double **ptr)
{
  comm_ptr = ptr;
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

void FixPIMD::nmpimd_transform(double **src, double **des, double *vector)
{
  int n = atom->nlocal;
  int m = 0;

  for (int i = 0; i < n; i++)
    for (int d = 0; d < 3; d++) {
      des[i][d] = 0.0;
      for (int j = 0; j < np; j++) { des[i][d] += (src[j][m] * vector[j]); }
      m++;
    }
}

/* ---------------------------------------------------------------------- */

void FixPIMD::spring_force()
{
  spring_energy = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double *_mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double *xlast = buf_beads[x_last];
  double *xnext = buf_beads[x_next];

  for (int i = 0; i < nlocal; i++) {
    double delx1 = xlast[0] - x[i][0];
    double dely1 = xlast[1] - x[i][1];
    double delz1 = xlast[2] - x[i][2];
    xlast += 3;
    domain->minimum_image(delx1, dely1, delz1);

    double delx2 = xnext[0] - x[i][0];
    double dely2 = xnext[1] - x[i][1];
    double delz2 = xnext[2] - x[i][2];
    xnext += 3;
    domain->minimum_image(delx2, dely2, delz2);

    double ff = fbond * _mass[type[i]];

    double dx = delx1 + delx2;
    double dy = dely1 + dely2;
    double dz = delz1 + delz2;

    f[i][0] -= (dx) *ff;
    f[i][1] -= (dy) *ff;
    f[i][2] -= (dz) *ff;

    spring_energy += (dx * dx + dy * dy + dz * dz);
  }
}

/* ----------------------------------------------------------------------
   Comm operations
------------------------------------------------------------------------- */

void FixPIMD::comm_init()
{
  if (size_plan) {
    delete[] plan_send;
    delete[] plan_recv;
  }

  if (method == PIMD) {
    size_plan = 2;
    plan_send = new int[2];
    plan_recv = new int[2];
    mode_index = new int[2];

    int rank_last = universe->me - comm->nprocs;
    int rank_next = universe->me + comm->nprocs;
    if (rank_last < 0) rank_last += universe->nprocs;
    if (rank_next >= universe->nprocs) rank_next -= universe->nprocs;

    plan_send[0] = rank_next;
    plan_send[1] = rank_last;
    plan_recv[0] = rank_last;
    plan_recv[1] = rank_next;

    mode_index[0] = 0;
    mode_index[1] = 1;
    x_last = 1;
    x_next = 0;
  } else {
    size_plan = np - 1;
    plan_send = new int[size_plan];
    plan_recv = new int[size_plan];
    mode_index = new int[size_plan];

    for (int i = 0; i < size_plan; i++) {
      plan_send[i] = universe->me + comm->nprocs * (i + 1);
      if (plan_send[i] >= universe->nprocs) plan_send[i] -= universe->nprocs;

      plan_recv[i] = universe->me - comm->nprocs * (i + 1);
      if (plan_recv[i] < 0) plan_recv[i] += universe->nprocs;

      mode_index[i] = (universe->iworld + i + 1) % (universe->nworlds);
    }

    x_next = (universe->iworld + 1 + universe->nworlds) % (universe->nworlds);
    x_last = (universe->iworld - 1 + universe->nworlds) % (universe->nworlds);
  }

  if (buf_beads) {
    for (int i = 0; i < np; i++) delete[] buf_beads[i];
    delete[] buf_beads;
  }

  buf_beads = new double *[np];
  for (int i = 0; i < np; i++) buf_beads[i] = nullptr;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::comm_exec(double **ptr)
{
  int nlocal = atom->nlocal;

  if (nlocal > max_nlocal) {
    max_nlocal = nlocal + 200;
    int size = sizeof(double) * max_nlocal * 3;
    buf_recv = (double *) memory->srealloc(buf_recv, size, "FixPIMD:x_recv");

    for (int i = 0; i < np; i++)
      buf_beads[i] = (double *) memory->srealloc(buf_beads[i], size, "FixPIMD:x_beads[i]");
  }

  // copy local positions

  memcpy(buf_beads[universe->iworld], &(ptr[0][0]), sizeof(double) * nlocal * 3);

  // go over comm plans

  for (int iplan = 0; iplan < size_plan; iplan++) {
    // sendrecv nlocal

    int nsend;

    MPI_Sendrecv(&(nlocal), 1, MPI_INT, plan_send[iplan], 0, &(nsend), 1, MPI_INT, plan_recv[iplan],
                 0, universe->uworld, MPI_STATUS_IGNORE);

    // allocate arrays

    if (nsend > max_nsend) {
      max_nsend = nsend + 200;
      tag_send =
          (tagint *) memory->srealloc(tag_send, sizeof(tagint) * max_nsend, "FixPIMD:tag_send");
      buf_send =
          (double *) memory->srealloc(buf_send, sizeof(double) * max_nsend * 3, "FixPIMD:x_send");
    }

    // send tags

    MPI_Sendrecv(atom->tag, nlocal, MPI_LMP_TAGINT, plan_send[iplan], 0, tag_send, nsend,
                 MPI_LMP_TAGINT, plan_recv[iplan], 0, universe->uworld, MPI_STATUS_IGNORE);

    // wrap positions

    double *wrap_ptr = buf_send;
    int ncpy = sizeof(double) * 3;

    for (int i = 0; i < nsend; i++) {
      int index = atom->map(tag_send[i]);

      if (index < 0) {
        auto mesg = fmt::format("Atom {} is missing at world [{}] rank [{}] "
                                "required by rank [{}] ({}, {}, {}).\n",
                                tag_send[i], universe->iworld, comm->me, plan_recv[iplan],
                                atom->tag[0], atom->tag[1], atom->tag[2]);
        error->universe_one(FLERR, mesg);
      }

      memcpy(wrap_ptr, ptr[index], ncpy);
      wrap_ptr += 3;
    }

    // sendrecv x

    MPI_Sendrecv(buf_send, nsend * 3, MPI_DOUBLE, plan_recv[iplan], 0, buf_recv, nlocal * 3,
                 MPI_DOUBLE, plan_send[iplan], 0, universe->uworld, MPI_STATUS_IGNORE);

    // copy x

    memcpy(buf_beads[mode_index[iplan]], buf_recv, sizeof(double) * nlocal * 3);
  }
}

/* ---------------------------------------------------------------------- */

int FixPIMD::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;

  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = comm_ptr[j][0];
    buf[m++] = comm_ptr[j][1];
    buf[m++] = comm_ptr[j][2];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    comm_ptr[i][0] = buf[m++];
    comm_ptr[i][1] = buf[m++];
    comm_ptr[i][2] = buf[m++];
  }
}

/* ----------------------------------------------------------------------
   Memory operations
------------------------------------------------------------------------- */

double FixPIMD::memory_usage()
{
  return (double) atom->nmax * size_peratom_cols * sizeof(double);
}

/* ---------------------------------------------------------------------- */

void FixPIMD::grow_arrays(int nmax)
{
  if (nmax == 0) return;
  int count = nmax * 3;

  memory->grow(array_atom, nmax, size_peratom_cols, "FixPIMD::array_atom");
  memory->grow(nhc_eta, count, nhc_nchain, "FixPIMD::nh_eta");
  memory->grow(nhc_eta_dot, count, nhc_nchain + 1, "FixPIMD::nh_eta_dot");
  memory->grow(nhc_eta_dotdot, count, nhc_nchain, "FixPIMD::nh_eta_dotdot");
  memory->grow(nhc_eta_mass, count, nhc_nchain, "FixPIMD::nh_eta_mass");
}

/* ---------------------------------------------------------------------- */

void FixPIMD::copy_arrays(int i, int j, int /*delflag*/)
{
  int i_pos = i * 3;
  int j_pos = j * 3;

  memcpy(nhc_eta[j_pos], nhc_eta[i_pos], nhc_size_one_1);
  memcpy(nhc_eta_dot[j_pos], nhc_eta_dot[i_pos], nhc_size_one_2);
  memcpy(nhc_eta_dotdot[j_pos], nhc_eta_dotdot[i_pos], nhc_size_one_1);
  memcpy(nhc_eta_mass[j_pos], nhc_eta_mass[i_pos], nhc_size_one_1);
}

/* ---------------------------------------------------------------------- */

int FixPIMD::pack_exchange(int i, double *buf)
{
  int offset = 0;
  int pos = i * 3;

  memcpy(buf + offset, nhc_eta[pos], nhc_size_one_1);
  offset += nhc_offset_one_1;
  memcpy(buf + offset, nhc_eta_dot[pos], nhc_size_one_2);
  offset += nhc_offset_one_2;
  memcpy(buf + offset, nhc_eta_dotdot[pos], nhc_size_one_1);
  offset += nhc_offset_one_1;
  memcpy(buf + offset, nhc_eta_mass[pos], nhc_size_one_1);
  offset += nhc_offset_one_1;

  return size_peratom_cols;
}

/* ---------------------------------------------------------------------- */

int FixPIMD::unpack_exchange(int nlocal, double *buf)
{
  int offset = 0;
  int pos = nlocal * 3;

  memcpy(nhc_eta[pos], buf + offset, nhc_size_one_1);
  offset += nhc_offset_one_1;
  memcpy(nhc_eta_dot[pos], buf + offset, nhc_size_one_2);
  offset += nhc_offset_one_2;
  memcpy(nhc_eta_dotdot[pos], buf + offset, nhc_size_one_1);
  offset += nhc_offset_one_1;
  memcpy(nhc_eta_mass[pos], buf + offset, nhc_size_one_1);
  offset += nhc_offset_one_1;

  return size_peratom_cols;
}

/* ---------------------------------------------------------------------- */

int FixPIMD::pack_restart(int i, double *buf)
{
  int offset = 0;
  int pos = i * 3;
  // pack buf[0] this way because other fixes unpack it
  buf[offset++] = size_peratom_cols + 1;

  memcpy(buf + offset, nhc_eta[pos], nhc_size_one_1);
  offset += nhc_offset_one_1;
  memcpy(buf + offset, nhc_eta_dot[pos], nhc_size_one_2);
  offset += nhc_offset_one_2;
  memcpy(buf + offset, nhc_eta_dotdot[pos], nhc_size_one_1);
  offset += nhc_offset_one_1;
  memcpy(buf + offset, nhc_eta_mass[pos], nhc_size_one_1);
  offset += nhc_offset_one_1;

  return size_peratom_cols + 1;
}

/* ---------------------------------------------------------------------- */

void FixPIMD::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int>(extra[nlocal][m]);
  m++;

  int pos = nlocal * 3;

  memcpy(nhc_eta[pos], extra[nlocal] + m, nhc_size_one_1);
  m += nhc_offset_one_1;
  memcpy(nhc_eta_dot[pos], extra[nlocal] + m, nhc_size_one_2);
  m += nhc_offset_one_2;
  memcpy(nhc_eta_dotdot[pos], extra[nlocal] + m, nhc_size_one_1);
  m += nhc_offset_one_1;
  memcpy(nhc_eta_mass[pos], extra[nlocal] + m, nhc_size_one_1);
  m += nhc_offset_one_1;

  nhc_ready = true;
}

/* ---------------------------------------------------------------------- */

int FixPIMD::maxsize_restart()
{
  return size_peratom_cols + 1;
}

/* ---------------------------------------------------------------------- */

int FixPIMD::size_restart(int /*nlocal*/)
{
  return size_peratom_cols + 1;
}

/* ---------------------------------------------------------------------- */

double FixPIMD::compute_vector(int n)
{
  if (n == 0) { return spring_energy; }
  if (n == 1) { return t_sys; }
  return 0.0;
}
