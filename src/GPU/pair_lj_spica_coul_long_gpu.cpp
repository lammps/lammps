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
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "pair_lj_spica_coul_long_gpu.h"

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

int spical_gpu_init(const int ntypes, double **cutsq, int **lj_type, double **host_lj1,
                    double **host_lj2, double **host_lj3, double **host_lj4, double **offset,
                    double *special_lj, const int nlocal, const int nall, const int max_nbors,
                    const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                    double **host_cut_ljsq, double host_cut_coulsq, double *host_special_coul,
                    const double qqrd2e, const double g_ewald);
void spical_gpu_clear();
int **spical_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                           int *host_type, double *sublo, double *subhi, tagint *tag,
                           int **nspecial, tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start, int **ilist,
                           int **jnum, const double cpu_time, bool &success, double *host_q,
                           double *boxlo, double *prd);
void spical_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                        int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                        const bool vflag, const bool eatom, const bool vatom, int &host_start,
                        const double cpu_time, bool &success, double *host_q, const int nlocal,
                        double *boxlo, double *prd);
double spical_gpu_bytes();

#include "lj_spica_common.h"

using namespace LJSPICAParms;

/* ---------------------------------------------------------------------- */

PairLJSPICACoulLongGPU::PairLJSPICACoulLongGPU(LAMMPS *lmp) :
    PairLJSPICACoulLong(lmp), gpu_mode(GPU_FORCE)
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

PairLJSPICACoulLongGPU::~PairLJSPICACoulLongGPU()
{
  spical_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairLJSPICACoulLongGPU::compute(int eflag, int vflag)
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
    firstneigh = spical_gpu_compute_n(neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi,
                                      atom->tag, atom->nspecial, atom->special, eflag, vflag,
                                      eflag_atom, vflag_atom, host_start, &ilist, &numneigh,
                                      cpu_time, success, atom->q, domain->boxlo, domain->prd);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    spical_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh, firstneigh,
                       eflag, vflag, eflag_atom, vflag_atom, host_start, cpu_time, success, atom->q,
                       atom->nlocal, domain->boxlo, domain->prd);
  }
  if (!success) error->one(FLERR, "Insufficient memory on accelerator");

  if (host_start < inum) {
    cpu_time = platform::walltime();
    if (evflag) {
      if (eflag)
        cpu_compute<1, 1>(host_start, inum, ilist, numneigh, firstneigh);
      else
        cpu_compute<1, 0>(host_start, inum, ilist, numneigh, firstneigh);
    } else
      cpu_compute<0, 0>(host_start, inum, ilist, numneigh, firstneigh);
    cpu_time = platform::walltime() - cpu_time;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJSPICACoulLongGPU::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR, "Pair style lj/spica/coul/long/gpu requires atom attribute q");

  // Repeat cutsq calculation because done after call to init_style
  double maxcut = -1.0;
  double cut;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        cut = init_one(i, j);
        cut *= cut;
        if (cut > maxcut) maxcut = cut;
        cutsq[i][j] = cutsq[j][i] = cut;
      } else
        cutsq[i][j] = cutsq[j][i] = 0.0;
    }
  }
  double cell_size = sqrt(maxcut) + neighbor->skin;

  cut_coulsq = cut_coul * cut_coul;

  // insure use of KSpace long-range solver, set g_ewald

  if (force->kspace == nullptr) error->all(FLERR, "Pair style is incompatible with KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables

  if (ncoultablebits) init_tables(cut_coul, nullptr);

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;
  int mnf = 5e-2 * neighbor->oneatom;
  int success = spical_gpu_init(atom->ntypes + 1, cutsq, lj_type, lj1, lj2, lj3, lj4, offset,
                                force->special_lj, atom->nlocal, atom->nlocal + atom->nghost, mnf,
                                maxspecial, cell_size, gpu_mode, screen, cut_ljsq, cut_coulsq,
                                force->special_coul, force->qqrd2e, g_ewald);
  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

double PairLJSPICACoulLongGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + spical_gpu_bytes();
}

/* ---------------------------------------------------------------------- */
template <int EVFLAG, int EFLAG>
void PairLJSPICACoulLongGPU::cpu_compute(int start, int inum, int *ilist, int *numneigh,
                                         int **firstneigh)
{
  int i, j, ii, jj;
  double qtmp, xtmp, ytmp, ztmp;
  double r2inv, forcecoul, forcelj, factor_coul, factor_lj;

  const double *const *const x = atom->x;
  double *const *const f = atom->f;
  const double *const q = atom->q;
  const int *const type = atom->type;
  const double *const special_coul = force->special_coul;
  const double *const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  double fxtmp, fytmp, fztmp;

  // loop over neighbors of my atoms

  for (ii = start; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp = fytmp = fztmp = 0.0;

    const int itype = type[i];
    const int *const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx * delx + dely * dely + delz * delz;
      const int jtype = type[j];

      double evdwl = 0.0;
      double ecoul = 0.0;
      double fpair = 0.0;

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0 / rsq;
        const int ljt = lj_type[itype][jtype];

        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            const double r = sqrt(rsq);
            const double grij = g_ewald * r;
            const double expm2 = exp(-grij * grij);
            const double t = 1.0 / (1.0 + EWALD_P * grij);
            const double erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
            const double prefactor = qqrd2e * qtmp * q[j] / r;
            forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);
            if (EFLAG) ecoul = prefactor * erfc;
            if (factor_coul < 1.0) {
              forcecoul -= (1.0 - factor_coul) * prefactor;
              if (EFLAG) ecoul -= (1.0 - factor_coul) * prefactor;
            }
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            int itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            const double fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
            const double table = ftable[itable] + fraction * dftable[itable];
            forcecoul = qtmp * q[j] * table;
            if (EFLAG) {
              const double table2 = etable[itable] + fraction * detable[itable];
              ecoul = qtmp * q[j] * table2;
            }
            if (factor_coul < 1.0) {
              const double table2 = ctable[itable] + fraction * dctable[itable];
              const double prefactor = qtmp * q[j] * table2;
              forcecoul -= (1.0 - factor_coul) * prefactor;
              if (EFLAG) ecoul -= (1.0 - factor_coul) * prefactor;
            }
          }
        } else {
          forcecoul = 0.0;
          ecoul = 0.0;
        }

        if (rsq < cut_ljsq[itype][jtype]) {

          if (ljt == LJ12_4) {
            const double r4inv = r2inv * r2inv;
            forcelj = r4inv * (lj1[itype][jtype] * r4inv * r4inv - lj2[itype][jtype]);

            if (EFLAG)
              evdwl = r4inv * (lj3[itype][jtype] * r4inv * r4inv - lj4[itype][jtype]) -
                  offset[itype][jtype];

          } else if (ljt == LJ9_6) {
            const double r3inv = r2inv * sqrt(r2inv);
            const double r6inv = r3inv * r3inv;
            forcelj = r6inv * (lj1[itype][jtype] * r3inv - lj2[itype][jtype]);
            if (EFLAG)
              evdwl =
                  r6inv * (lj3[itype][jtype] * r3inv - lj4[itype][jtype]) - offset[itype][jtype];

          } else if (ljt == LJ12_6) {
            const double r6inv = r2inv * r2inv * r2inv;
            forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
            if (EFLAG)
              evdwl =
                  r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];

          } else if (ljt == LJ12_5) {
            const double r5inv = r2inv * r2inv * sqrt(r2inv);
            const double r7inv = r5inv * r2inv;
            forcelj = r5inv * (lj1[itype][jtype] * r7inv - lj2[itype][jtype]);
            if (EFLAG)
              evdwl =
                  r5inv * (lj3[itype][jtype] * r7inv - lj4[itype][jtype]) - offset[itype][jtype];
          }

          if (EFLAG) evdwl *= factor_lj;

        } else {
          forcelj = 0.0;
          evdwl = 0.0;
        }

        fpair = (forcecoul + factor_lj * forcelj) * r2inv;

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;

        if (EVFLAG) ev_tally_full(i, evdwl, ecoul, fpair, delx, dely, delz);
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}
