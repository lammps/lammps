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
   Contributing author: Trung Dac Nguyen (ORNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_coul_dsf_gpu.h"
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
#include "gpu_extra.h"

#define MY_PIS 1.77245385090551602729
#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int cdsf_gpu_init(const int ntypes, const int nlocal, const int nall,
                  const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen,
                  const double host_cut_coulsq,
                  double *host_special_coul, const double qqrd2e,
                  const double e_shift, const double f_shift,
                  const double alpha);
void cdsf_gpu_clear();
int ** cdsf_gpu_compute_n(const int ago, const int inum,
                          const int nall, double **host_x, int *host_type,
                          double *sublo, double *subhi, tagint *tag, int **nspecial,
                          tagint **special, const bool eflag, const bool vflag,
                          const bool eatom, const bool vatom, int &host_start,
                          int **ilist, int **jnum, const double cpu_time,
                          bool &success, double *host_q, double *boxlo,
                          double *prd);
void cdsf_gpu_compute(const int ago, const int inum,
                      const int nall, double **host_x, int *host_type,
                      int *ilist, int *numj, int **firstneigh,
                      const bool eflag, const bool vflag, const bool eatom,
                      const bool vatom, int &host_start, const double cpu_time,
                      bool &success, double *host_q, const int nlocal,
                      double *boxlo, double *prd);
double cdsf_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairCoulDSFGPU::PairCoulDSFGPU(LAMMPS *lmp) : PairCoulDSF(lmp),
  gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  reinitflag = 0;
  cpu_time = 0.0;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairCoulDSFGPU::~PairCoulDSFGPU()
{
  cdsf_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairCoulDSFGPU::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = cdsf_gpu_compute_n(neighbor->ago, inum, nall, atom->x,
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
    cdsf_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type,
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

void PairCoulDSFGPU::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style coul/dsf/gpu requires atom attribute q");

  if (force->newton_pair)
    error->all(FLERR,"Cannot use newton pair with coul/dsf/gpu pair style");

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
  double erfcc = erfc(alpha*cut_coul);
  double erfcd = exp(-alpha*alpha*cut_coul*cut_coul);
  f_shift = -(erfcc/cut_coulsq + 2.0/MY_PIS*alpha*erfcd/cut_coul);
  e_shift = erfcc/cut_coul - f_shift*cut_coul;

  int maxspecial=0;
  if (atom->molecular)
    maxspecial=atom->maxspecial;
  int success = cdsf_gpu_init(atom->ntypes+1, atom->nlocal,
                              atom->nlocal+atom->nghost, 300, maxspecial,
                              cell_size, gpu_mode, screen, cut_coulsq,
                              force->special_coul, force->qqrd2e, e_shift,
                              f_shift, alpha);
  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE) {
    int irequest = neighbor->request(this,instance_me);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  }
}

/* ---------------------------------------------------------------------- */

double PairCoulDSFGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + cdsf_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairCoulDSFGPU::cpu_compute(int start, int inum, int eflag,
                                 int /* vflag */, int *ilist,
                                 int *numneigh, int **firstneigh)
{
  int i,j,ii,jj,jnum;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double r,rsq,r2inv,forcecoul,factor_coul;
  double prefactor,erfcc,erfcd,t;
  int *jlist;

  ecoul = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double qqrd2e = force->qqrd2e;

  // loop over neighbors of my atoms

  for (ii = start; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (evflag) {
      double e_self = -(e_shift/2.0 + alpha/MY_PIS) * qtmp*qtmp*qqrd2e;
      ev_tally(i,i,nlocal,0,0.0,e_self,0.0,0.0,0.0,0.0);
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cut_coulsq) {
        r2inv = 1.0/rsq;
        r = sqrt(rsq);
        prefactor = qqrd2e*qtmp*q[j]/r;
        erfcd = exp(-alpha*alpha*r*r);
        t = 1.0 / (1.0 + EWALD_P*alpha*r);
        erfcc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * erfcd;
        forcecoul = prefactor * (erfcc/r + 2.0*alpha/MY_PIS * erfcd +
          r*f_shift) * r;
        if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;

        fpair = forcecoul * r2inv;
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;

        if (eflag) {
          if (rsq < cut_coulsq) {
            ecoul = prefactor * (erfcc - r*e_shift - rsq*f_shift);
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;
        }

        if (evflag) ev_tally_full(i,0.0,ecoul,fpair,delx,dely,delz);
      }
    }
  }
}
