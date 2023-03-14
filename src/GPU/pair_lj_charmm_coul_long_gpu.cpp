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
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "pair_lj_charmm_coul_long_gpu.h"

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

#define EWALD_F 1.12837917
#define EWALD_P 0.3275911
#define A1 0.254829592
#define A2 -0.284496736
#define A3 1.421413741
#define A4 -1.453152027
#define A5 1.061405429

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int crml_gpu_init(const int ntypes, double cut_bothsq, double **host_lj1, double **host_lj2,
                  double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                  const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen, double host_cut_ljsq,
                  double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                  const double g_ewald, const double cut_lj_innersq, const double denom_lj,
                  double **epsilon, double **sigma, const bool mix_arithmetic);
void crml_gpu_clear();
int **crml_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int &host_start, int **ilist, int **jnum,
                         const double cpu_time, bool &success, double *host_q, double *boxlo,
                         double *prd);
void crml_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, int &host_start,
                      const double cpu_time, bool &success, double *host_q, const int nlocal,
                      double *boxlo, double *prd);
double crml_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairLJCharmmCoulLongGPU::PairLJCharmmCoulLongGPU(LAMMPS *lmp) :
    PairLJCharmmCoulLong(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  reinitflag = 0;
  cpu_time = 0.0;
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairLJCharmmCoulLongGPU::~PairLJCharmmCoulLongGPU()
{
  crml_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulLongGPU::compute(int eflag, int vflag)
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
    firstneigh = crml_gpu_compute_n(neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi,
                                    atom->tag, atom->nspecial, atom->special, eflag, vflag,
                                    eflag_atom, vflag_atom, host_start, &ilist, &numneigh, cpu_time,
                                    success, atom->q, domain->boxlo, domain->prd);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    crml_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh, firstneigh,
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

void PairLJCharmmCoulLongGPU::init_style()
{
  cut_respa = nullptr;

  if (!atom->q_flag)
    error->all(FLERR, "Pair style lj/charmm/coul/long/gpu requires atom attribute q");

  // Repeat cutsq calculation because done after call to init_style

  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) init_one(i, j);
    }
  }

  cut_lj_innersq = cut_lj_inner * cut_lj_inner;
  cut_ljsq = cut_lj * cut_lj;
  cut_coulsq = cut_coul * cut_coul;
  cut_bothsq = MAX(cut_ljsq, cut_coulsq);

  denom_lj =
      (cut_ljsq - cut_lj_innersq) * (cut_ljsq - cut_lj_innersq) * (cut_ljsq - cut_lj_innersq);

  double cell_size = sqrt(cut_bothsq) + neighbor->skin;

  // ensure use of KSpace long-range solver, set g_ewald

  if (force->kspace == nullptr) error->all(FLERR, "Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables

  if (ncoultablebits) init_tables(cut_coul, cut_respa);

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;

  bool arithmetic = true;
  for (int i = 1; i < atom->ntypes + 1; i++)
    for (int j = i + 1; j < atom->ntypes + 1; j++) {
      if (epsilon[i][j] != sqrt(epsilon[i][i] * epsilon[j][j])) arithmetic = false;
      if (sigma[i][j] != 0.5 * (sigma[i][i] + sigma[j][j])) arithmetic = false;
    }

  int mnf = 5e-2 * neighbor->oneatom;
  int success =
      crml_gpu_init(atom->ntypes + 1, cut_bothsq, lj1, lj2, lj3, lj4, offset, force->special_lj,
                    atom->nlocal, atom->nlocal + atom->nghost, mnf, maxspecial, cell_size, gpu_mode,
                    screen, cut_ljsq, cut_coulsq, force->special_coul, force->qqrd2e, g_ewald,
                    cut_lj_innersq, denom_lj, epsilon, sigma, arithmetic);
  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

double PairLJCharmmCoulLongGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + crml_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulLongGPU::cpu_compute(int start, int inum, int eflag, int /* vflag */,
                                          int *ilist, int *numneigh, int **firstneigh)
{
  int i, j, ii, jj, jnum, itype, jtype, itable;
  double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, evdwl, ecoul, fpair;
  double fraction, table;
  double r, r2inv, r6inv, forcecoul, forcelj, factor_coul, factor_lj;
  double grij, expm2, prefactor, t, erfc;
  double philj, switch1, switch2;
  int *jlist;
  double rsq;

  evdwl = ecoul = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
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
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cut_bothsq) {
        r2inv = 1.0 / rsq;

        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            r = sqrt(rsq);
            grij = g_ewald * r;
            expm2 = exp(-grij * grij);
            t = 1.0 / (1.0 + EWALD_P * grij);
            erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
            prefactor = qqrd2e * qtmp * q[j] / r;
            forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);
            if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
            table = ftable[itable] + fraction * dftable[itable];
            forcecoul = qtmp * q[j] * table;
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction * dctable[itable];
              prefactor = qtmp * q[j] * table;
              forcecoul -= (1.0 - factor_coul) * prefactor;
            }
          }
        } else
          forcecoul = 0.0;

        if (rsq < cut_ljsq) {
          r6inv = r2inv * r2inv * r2inv;
          jtype = type[j];
          forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
          if (rsq > cut_lj_innersq) {
            switch1 = (cut_ljsq - rsq) * (cut_ljsq - rsq) *
                (cut_ljsq + 2.0 * rsq - 3.0 * cut_lj_innersq) / denom_lj;
            switch2 = 12.0 * rsq * (cut_ljsq - rsq) * (rsq - cut_lj_innersq) / denom_lj;
            philj = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]);
            forcelj = forcelj * switch1 + philj * switch2;
          }
        } else
          forcelj = 0.0;

        fpair = (forcecoul + factor_lj * forcelj) * r2inv;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;

        if (eflag) {
          if (rsq < cut_coulsq) {
            if (!ncoultablebits || rsq <= tabinnersq)
              ecoul = prefactor * erfc;
            else {
              table = etable[itable] + fraction * detable[itable];
              ecoul = qtmp * q[j] * table;
            }
            if (factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
          } else
            ecoul = 0.0;

          if (rsq < cut_ljsq) {
            evdwl = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]);
            if (rsq > cut_lj_innersq) {
              switch1 = (cut_ljsq - rsq) * (cut_ljsq - rsq) *
                  (cut_ljsq + 2.0 * rsq - 3.0 * cut_lj_innersq) / denom_lj;
              evdwl *= switch1;
            }
            evdwl *= factor_lj;
          } else
            evdwl = 0.0;
        }

        if (evflag) ev_tally_full(i, evdwl, ecoul, fpair, delx, dely, delz);
      }
    }
  }
}
