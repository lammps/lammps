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
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "pair_lj_sdk_coul_long_gpu.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "integrate.h"
#include "memory.h"
#include "error.h"
#include "neigh_request.h"
#include "universe.h"
#include "update.h"
#include "domain.h"
#include "string.h"
#include "kspace.h"
#include "gpu_extra.h"

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

// External functions from cuda library for atom decomposition

int cmml_gpu_init(const int ntypes, double **cutsq, int **lj_type,
                  double **host_lj1, double **host_lj2, double **host_lj3,
                  double **host_lj4, double **offset, double *special_lj,
                  const int nlocal, const int nall, const int max_nbors,
                  const int maxspecial, const double cell_size, int &gpu_mode,
                  FILE *screen, double **host_cut_ljsq, double host_cut_coulsq,
                  double *host_special_coul, const double qqrd2e,
                  const double g_ewald);
void cmml_gpu_clear();
int ** cmml_gpu_compute_n(const int ago, const int inum, const int nall,
                          double **host_x, int *host_type, double *sublo,
                          double *subhi, int *tag, int **nspecial,
                          int **special, const bool eflag, const bool vflag,
                          const bool eatom, const bool vatom, int &host_start,
                          int **ilist, int **jnum, const double cpu_time,
                          bool &success, double *host_q, double *boxlo,
                          double *prd);
void cmml_gpu_compute(const int ago, const int inum, const int nall,
                      double **host_x, int *host_type, int *ilist, int *numj,
                      int **firstneigh, const bool eflag, const bool vflag,
                      const bool eatom, const bool vatom, int &host_start,
                      const double cpu_time, bool &success, double *host_q,
                      const int nlocal, double *boxlo, double *prd);
double cmml_gpu_bytes();

#include "lj_sdk_common.h"

using namespace LAMMPS_NS;
using namespace LJSDKParms;

/* ---------------------------------------------------------------------- */

PairLJSDKCoulLongGPU::PairLJSDKCoulLongGPU(LAMMPS *lmp) :
  PairLJSDKCoulLong(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  cpu_time = 0.0;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairLJSDKCoulLongGPU::~PairLJSDKCoulLongGPU()
{
  cmml_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairLJSDKCoulLongGPU::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = cmml_gpu_compute_n(neighbor->ago, inum, nall, atom->x,
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
    cmml_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type,
                     ilist, numneigh, firstneigh, eflag, vflag, eflag_atom,
                     vflag_atom, host_start, cpu_time, success, atom->q,
                     atom->nlocal, domain->boxlo, domain->prd);
  }
  if (!success)
    error->one(FLERR,"Insufficient memory on accelerator");

  if (host_start<inum) {
    cpu_time = MPI_Wtime();
    if (evflag) {
      if (eflag) cpu_compute<1,1>(host_start, inum, ilist, numneigh, firstneigh);
      else cpu_compute<1,0>(host_start, inum, ilist, numneigh, firstneigh);
    } else cpu_compute<0,0>(host_start, inum, ilist, numneigh, firstneigh);
    cpu_time = MPI_Wtime() - cpu_time;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJSDKCoulLongGPU::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style lj/sdk/coul/long/gpu requires atom attribute q");
  if (force->newton_pair)
    error->all(FLERR,"Cannot use newton pair with lj/sdk/coul/long/gpu pair style");

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
    error->all(FLERR,"Pair style is incompatible with KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables

  if (ncoultablebits) init_tables();

  int maxspecial=0;
  if (atom->molecular)
    maxspecial=atom->maxspecial;
  int success = cmml_gpu_init(atom->ntypes+1, cutsq, lj_type, lj1, lj2, lj3,
                              lj4, offset, force->special_lj, atom->nlocal,
                              atom->nlocal+atom->nghost, 300, maxspecial,
                              cell_size, gpu_mode, screen, cut_ljsq,
                              cut_coulsq, force->special_coul,
                              force->qqrd2e, g_ewald);
  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE) {
    int irequest = neighbor->request(this);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJSDKCoulLongGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + cmml_gpu_bytes();
}

/* ---------------------------------------------------------------------- */
template <int EVFLAG, int EFLAG>
void PairLJSDKCoulLongGPU::cpu_compute(int start, int inum, int *ilist,
                                       int *numneigh, int **firstneigh)
{
  int i,j,ii,jj;
  double qtmp,xtmp,ytmp,ztmp;
  double r,rsq,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;

  const double * const * const x = atom->x;
  double * const * const f = atom->f;
  const double * const q = atom->q;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_coul = force->special_coul;
  const double * const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  double fxtmp,fytmp,fztmp;

  // loop over neighbors of my atoms

  for (ii = start; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp=fytmp=fztmp=0.0;

    const int itype = type[i];
    const int * const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx*delx + dely*dely + delz*delz;
      const int jtype = type[j];

      double evdwl = 0.0;
      double ecoul = 0.0;
      double fpair = 0.0;

      if (rsq < cutsq[itype][jtype]) {
          r2inv = 1.0/rsq;
        const int ljt = lj_type[itype][jtype];

        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            const double r = sqrt(rsq);
            const double grij = g_ewald * r;
            const double expm2 = exp(-grij*grij);
            const double t = 1.0 / (1.0 + EWALD_P*grij);
            const double erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
            const double prefactor = qqrd2e * qtmp*q[j]/r;
            forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
            if (EFLAG) ecoul = prefactor*erfc;
            if (factor_coul < 1.0) {
              forcecoul -= (1.0-factor_coul)*prefactor;
              if (EFLAG) ecoul -= (1.0-factor_coul)*prefactor;
            }
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            int itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            const double fraction = (rsq_lookup.f - rtable[itable]) *
                                     drtable[itable];
            const double table = ftable[itable] + fraction*dftable[itable];
            forcecoul = qtmp*q[j] * table;
            if (EFLAG) {
              const double table2 = etable[itable] + fraction*detable[itable];
              ecoul = qtmp*q[j] * table2;
            }
            if (factor_coul < 1.0) {
              const double table2 = ctable[itable] + fraction*dctable[itable];
              const double prefactor = qtmp*q[j] * table2;
              forcecoul -= (1.0-factor_coul)*prefactor;
              if (EFLAG) ecoul -= (1.0-factor_coul)*prefactor;
            }
          }
        } else {
          forcecoul = 0.0;
          ecoul = 0.0;
        }


        if (rsq < cut_ljsq[itype][jtype]) {

          if (ljt == LJ12_4) {
            const double r4inv=r2inv*r2inv;
            forcelj = r4inv*(lj1[itype][jtype]*r4inv*r4inv
                             - lj2[itype][jtype]);

            if (EFLAG)
              evdwl = r4inv*(lj3[itype][jtype]*r4inv*r4inv
                             - lj4[itype][jtype]) - offset[itype][jtype];

          } else if (ljt == LJ9_6) {
            const double r3inv = r2inv*sqrt(r2inv);
            const double r6inv = r3inv*r3inv;
            forcelj = r6inv*(lj1[itype][jtype]*r3inv
                             - lj2[itype][jtype]);
            if (EFLAG)
              evdwl = r6inv*(lj3[itype][jtype]*r3inv
                             - lj4[itype][jtype]) - offset[itype][jtype];

          } else if (ljt == LJ12_6) {
            const double r6inv = r2inv*r2inv*r2inv;
            forcelj = r6inv*(lj1[itype][jtype]*r6inv
                             - lj2[itype][jtype]);
            if (EFLAG)
              evdwl = r6inv*(lj3[itype][jtype]*r6inv
                             - lj4[itype][jtype]) - offset[itype][jtype];
          }

          if (EFLAG) evdwl *= factor_lj;

        } else {
          forcelj=0.0;
          evdwl = 0.0;
        }

        fpair = (forcecoul + factor_lj*forcelj) * r2inv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;

        if (EVFLAG) ev_tally_full(i,evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}
