/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Dac Nguyen (Northwestern)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_born_coul_long_cs_gpu.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "integrate.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "neigh_request.h"
#include "universe.h"
#include "update.h"
#include "domain.h"
#include "kspace.h"
#include "gpu_extra.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define EWALD_F   1.12837917
#define EWALD_P   9.95473818e-1
#define B0       -0.1335096380159268
#define B1       -2.57839507e-1
#define B2       -1.37203639e-1
#define B3       -8.88822059e-3
#define B4       -5.80844129e-3
#define B5        1.14652755e-1

#define EPSILON 1.0e-20
#define EPS_EWALD 1.0e-6
#define EPS_EWALD_SQR 1.0e-12

// External functions from cuda library for atom decomposition

int bornclcs_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv,
                    double **host_born1, double **host_born2,
                    double **host_born3, double **host_a,
                    double **host_c, double **host_d,
                    double **sigma, double **offset, double *special_lj,
                    const int inum, const int nall, const int max_nbors,
                    const int maxspecial, const double cell_size,
                    int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                    double host_cut_coulsq, double *host_special_coul,
                    const double qqrd2e, const double g_ewald);
void bornclcs_gpu_clear();
int** bornclcs_gpu_compute_n(const int ago, const int inum_full, const int nall,
                           double **host_x, int *host_type, double *sublo,
                           double *subhi, tagint *tag, int **nspecial,
                           tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start,
                           int **ilist, int **jnum,  const double cpu_time,
                           bool &success, double *host_q, double *boxlo,
                           double *prd);
void bornclcs_gpu_compute(const int ago, const int inum_full, const int nall,
                        double **host_x, int *host_type, int *ilist, int *numj,
                        int **firstneigh, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        const double cpu_time, bool &success, double *host_q,
                        const int nlocal, double *boxlo, double *prd);
double bornclcs_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairBornCoulLongCSGPU::PairBornCoulLongCSGPU(LAMMPS *lmp) :
  PairBornCoulLongCS(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  reinitflag = 0;
  cpu_time = 0.0;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairBornCoulLongCSGPU::~PairBornCoulLongCSGPU()
{
  bornclcs_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairBornCoulLongCSGPU::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = bornclcs_gpu_compute_n(neighbor->ago, inum, nall, atom->x,
                                      atom->type, domain->sublo, domain->subhi,
                                      atom->tag, atom->nspecial, atom->special,
                                      eflag, vflag, eflag_atom, vflag_atom,
                                      host_start, &ilist, &numneigh, cpu_time,
                                      success, atom->q, domain->boxlo,
                                      domain->prd);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    bornclcs_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type,
                       ilist, numneigh, firstneigh, eflag, vflag, eflag_atom,
                       vflag_atom, host_start, cpu_time, success, atom->q,
                       atom->nlocal, domain->boxlo, domain->prd);
  }
  if (!success)
    error->one(FLERR,"Insufficient memory on accelerator");

  if (host_start<inum) {
    cpu_time = MPI_Wtime();
    cpu_compute(host_start, inum, eflag, vflag, ilist, numneigh, firstneigh);
    cpu_time = MPI_Wtime() - cpu_time;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBornCoulLongCSGPU::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,
      "Pair style born/coul/long/cs/gpu requires atom attribute q");
  if (force->newton_pair)
    error->all(FLERR,
       "Cannot use newton pair with born/coul/long/cs/gpu pair style");

  // Repeat cutsq calculation because done after call to init_style
  double maxcut = -1.0;
  double cut;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        cut = init_one(i,j);
        cut *= cut;
        if (cut > maxcut)
          maxcut = cut;
        cutsq[i][j] = cutsq[j][i] = cut;
      } else
        cutsq[i][j] = cutsq[j][i] = 0.0;
    }
  }
  double cell_size = sqrt(maxcut) + neighbor->skin;

  cut_coulsq = cut_coul * cut_coul;

  // insure use of KSpace long-range solver, set g_ewald

  if (force->kspace == NULL)
    error->all(FLERR,"Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  int maxspecial=0;
  if (atom->molecular)
    maxspecial=atom->maxspecial;
  int success = bornclcs_gpu_init(atom->ntypes+1, cutsq,  rhoinv,
                                born1, born2, born3, a, c, d, sigma,
                                offset, force->special_lj, atom->nlocal,
                                  atom->nlocal+atom->nghost, 300, maxspecial,
                                   cell_size, gpu_mode, screen, cut_ljsq,
                                cut_coulsq, force->special_coul,
                                force->qqrd2e, g_ewald);

  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE) {
    int irequest = neighbor->request(this,instance_me);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  }
}

/* ---------------------------------------------------------------------- */

double PairBornCoulLongCSGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + bornclcs_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairBornCoulLongCSGPU::cpu_compute(int start, int inum, int eflag,
                                      int /* vflag */, int *ilist,
                                      int *numneigh, int **firstneigh)
{
  int i,j,ii,jj,jnum,itable,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double fraction,table;
  double r,rsq,rexp,r2inv,r6inv,forcecoul,forceborn,factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc,u;
  int *jlist;

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
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {

        if (rsq < cut_coulsq) {
          rsq += EPSILON; // Add Epsilon for case: r = 0; Interaction must be removed by special bond;
          r2inv = 1.0/rsq;
          if (!ncoultablebits || rsq <= tabinnersq) {
            r = sqrt(rsq);
            prefactor = qqrd2e * qtmp*q[j];
            if (factor_coul < 1.0) {
              // When bonded parts are being calculated a minimal distance (EPS_EWALD)
              // has to be added to the prefactor and erfc in order to make the
              // used approximation functions for the Ewald correction valid
              grij = g_ewald * (r+EPS_EWALD);
              expm2 = exp(-grij*grij);
              t = 1.0 / (1.0 + EWALD_P*grij);
              u = 1.0 - t;
              erfc = t * (1.+u*(B0+u*(B1+u*(B2+u*(B3+u*(B4+u*B5)))))) * expm2;
              prefactor /= (r+EPS_EWALD);
              forcecoul = prefactor * (erfc + EWALD_F*grij*expm2 - (1.0-factor_coul));
              // Additionally r2inv needs to be accordingly modified since the later
              // scaling of the overall force shall be consistent
              r2inv = 1.0/(rsq + EPS_EWALD_SQR);
            } else {
              grij = g_ewald * r;
              expm2 = exp(-grij*grij);
              t = 1.0 / (1.0 + EWALD_P*grij);
              u = 1.0 - t;
              erfc = t * (1.+u*(B0+u*(B1+u*(B2+u*(B3+u*(B4+u*B5)))))) * expm2;
              prefactor /= r;
              forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
            }
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
            table = ftable[itable] + fraction*dftable[itable];
            forcecoul = qtmp*q[j] * table;
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction*dctable[itable];
              prefactor = qtmp*q[j] * table;
              forcecoul -= (1.0-factor_coul)*prefactor;
            }
          }

          forcecoul *= r2inv;

        } else forcecoul = 0;

        r2inv = 1.0/rsq;
        r = sqrt(rsq);
        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          rexp = exp((sigma[itype][jtype]-r)*rhoinv[itype][jtype]);
          forceborn = born1[itype][jtype]*r*rexp - born2[itype][jtype]*r6inv
            + born3[itype][jtype]*r2inv*r6inv;
        } else forceborn = 0.0;

        fpair = forcecoul + factor_lj*forceborn * r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;

        if (eflag) {
          if (rsq < cut_coulsq) {
            ecoul = prefactor*erfc;
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;
          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = a[itype][jtype]*rexp - c[itype][jtype]*r6inv
              + d[itype][jtype]*r6inv*r2inv - offset[itype][jtype];
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (evflag) ev_tally_full(i,evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
  }
}
