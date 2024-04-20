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
   Contributing author: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#include "pair_coul_long_cs_gpu.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "kspace.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

using namespace LAMMPS_NS;

static constexpr double EWALD_F =  1.12837917;
static constexpr double EWALD_P =  9.95473818e-1;
static constexpr double B0      = -0.1335096380159268;
static constexpr double B1      = -2.57839507e-1;
static constexpr double B2      = -1.37203639e-1;
static constexpr double B3      = -8.88822059e-3;
static constexpr double B4      = -5.80844129e-3;
static constexpr double B5      =  1.14652755e-1;

static constexpr double EPSILON = 1.0e-20;
static constexpr double EPS_EWALD = 1.0e-6;
static constexpr double EPS_EWALD_SQR = 1.0e-12;

// External functions from cuda library for atom decomposition

int clcs_gpu_init(const int ntypes, double **scale, const int nlocal, const int nall,
                  const int max_nbors, const int maxspecial, const double cell_size, int &gpu_mode,
                  FILE *screen, double host_cut_coulsq, double *host_special_coul,
                  const double qqrd2e, const double g_ewald);
void clcs_gpu_reinit(const int ntypes, double **scale);
void clcs_gpu_clear();
int **clcs_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int &host_start, int **ilist, int **jnum,
                         const double cpu_time, bool &success, double *host_q, double *boxlo,
                         double *prd);
void clcs_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, int &host_start,
                      const double cpu_time, bool &success, double *host_q, const int nlocal,
                      double *boxlo, double *prd);
double clcs_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairCoulLongCSGPU::PairCoulLongCSGPU(LAMMPS *lmp) : PairCoulLongCS(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  cpu_time = 0.0;
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairCoulLongCSGPU::~PairCoulLongCSGPU()
{
  clcs_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairCoulLongCSGPU::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    double sublo[3], subhi[3];
    if (domain->triclinic == 0) {
      sublo[0] = domain->sublo[0];
      sublo[1] = domain->sublo[1];
      sublo[2] = domain->sublo[2];
      subhi[0] = domain->subhi[0];
      subhi[1] = domain->subhi[1];
      subhi[2] = domain->subhi[2];
    } else {
      domain->bbox(domain->sublo_lamda, domain->subhi_lamda, sublo, subhi);
    }
    inum = atom->nlocal;
    firstneigh = clcs_gpu_compute_n(neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi,
                                    atom->tag, atom->nspecial, atom->special, eflag, vflag,
                                    eflag_atom, vflag_atom, host_start, &ilist, &numneigh, cpu_time,
                                    success, atom->q, domain->boxlo, domain->prd);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    clcs_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh, firstneigh,
                     eflag, vflag, eflag_atom, vflag_atom, host_start, cpu_time, success, atom->q,
                     atom->nlocal, domain->boxlo, domain->prd);
  }
  if (!success) error->one(FLERR, "Insufficient memory on accelerator");

  if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0)
    neighbor->build_topology();
  if (host_start < inum) {
    cpu_time = platform::walltime();
    cpu_compute(host_start, inum, eflag, vflag, ilist, numneigh, firstneigh);
    cpu_time = platform::walltime() - cpu_time;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCoulLongCSGPU::init_style()
{
  cut_respa = nullptr;

  if (!atom->q_flag) error->all(FLERR, "Pair style coul/long/cs/gpu requires atom attribute q");

  // Call init_one calculation make sure scale is correct
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) { init_one(i, j); }
    }
  }
  double cell_size = cut_coul + neighbor->skin;

  cut_coulsq = cut_coul * cut_coul;

  // ensure use of KSpace long-range solver, set g_ewald

  if (force->kspace == nullptr) error->all(FLERR, "Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables

  if (ncoultablebits) init_tables(cut_coul, cut_respa);

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;
  int mnf = 5e-2 * neighbor->oneatom;
  int success = clcs_gpu_init(atom->ntypes + 1, scale, atom->nlocal, atom->nlocal + atom->nghost,
                              mnf, maxspecial, cell_size, gpu_mode, screen, cut_coulsq,
                              force->special_coul, force->qqrd2e, g_ewald);

  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

void PairCoulLongCSGPU::reinit()
{
  Pair::reinit();

  clcs_gpu_reinit(atom->ntypes + 1, scale);
}

/* ---------------------------------------------------------------------- */

double PairCoulLongCSGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + clcs_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairCoulLongCSGPU::cpu_compute(int start, int inum, int eflag, int /* vflag */, int *ilist,
                                    int *numneigh, int **firstneigh)
{
  int i, j, ii, jj, jnum, itable, itype, jtype;
  double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, ecoul, fpair;
  double fraction, table;
  double r, r2inv, forcecoul, factor_coul;
  double grij, expm2, prefactor, t, erfc, u;
  int *jlist;
  double rsq;

  ecoul = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  double *special_coul = force->special_coul;
  double qqrd2e = force->qqrd2e;

  // loop over neighbors of my atoms

  for (ii = start; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cut_coulsq) {
        rsq +=
            EPSILON;    // Add Epsilon for case: r = 0; Interaction must be removed by special bond;
        r2inv = 1.0 / rsq;
        if (!ncoultablebits || rsq <= tabinnersq) {
          r = sqrt(rsq);
          prefactor = qqrd2e * scale[itype][jtype] * qtmp * q[j];
          if (factor_coul < 1.0) {
            // When bonded parts are being calculated a minimal distance (EPS_EWALD)
            // has to be added to the prefactor and erfc in order to make the
            // used approximation functions for the Ewald correction valid
            grij = g_ewald * (r + EPS_EWALD);
            expm2 = exp(-grij * grij);
            t = 1.0 / (1.0 + EWALD_P * grij);
            u = 1.0 - t;
            erfc = t * (1. + u * (B0 + u * (B1 + u * (B2 + u * (B3 + u * (B4 + u * B5)))))) * expm2;
            prefactor /= (r + EPS_EWALD);
            forcecoul = prefactor * (erfc + EWALD_F * grij * expm2 - (1.0 - factor_coul));
            // Additionally r2inv needs to be accordingly modified since the later
            // scaling of the overall force shall be consistent
            r2inv = 1.0 / (rsq + EPS_EWALD_SQR);
          } else {
            grij = g_ewald * r;
            expm2 = exp(-grij * grij);
            t = 1.0 / (1.0 + EWALD_P * grij);
            u = 1.0 - t;
            erfc = t * (1. + u * (B0 + u * (B1 + u * (B2 + u * (B3 + u * (B4 + u * B5)))))) * expm2;
            prefactor /= r;
            forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);
          }
        } else {
          union_int_float_t rsq_lookup;
          rsq_lookup.f = rsq;
          itable = rsq_lookup.i & ncoulmask;
          itable >>= ncoulshiftbits;
          fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
          table = ftable[itable] + fraction * dftable[itable];
          forcecoul = scale[itype][jtype] * qtmp * q[j] * table;
          if (factor_coul < 1.0) {
            table = ctable[itable] + fraction * dctable[itable];
            prefactor = scale[itype][jtype] * qtmp * q[j] * table;
            forcecoul -= (1.0 - factor_coul) * prefactor;
          }
        }

        fpair = forcecoul * r2inv;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;

        if (eflag) {
          if (rsq < cut_coulsq) {
            if (!ncoultablebits || rsq <= tabinnersq)
              ecoul = prefactor * erfc;
            else {
              table = etable[itable] + fraction * detable[itable];
              ecoul = scale[itype][jtype] * qtmp * q[j] * table;
            }
            if (factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
          } else
            ecoul = 0.0;
        }

        if (evflag) ev_tally_full(i, 0.0, ecoul, fpair, delx, dely, delz);
      }
    }
  }
}
