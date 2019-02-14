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
   Contributing author: Axel Kohlmeyer (Temple)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_coul_long_gpu.h"
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
#include "kspace.h"
#include "gpu_extra.h"

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int cl_gpu_init(const int ntypes, double **scale,
                const int nlocal, const int nall, const int max_nbors,
                const int maxspecial, const double cell_size, int &gpu_mode,
                FILE *screen, double host_cut_coulsq, double *host_special_coul,
                const double qqrd2e, const double g_ewald);
void cl_gpu_reinit(const int ntypes, double **scale);
void cl_gpu_clear();
int ** cl_gpu_compute_n(const int ago, const int inum,
                        const int nall, double **host_x, int *host_type,
                        double *sublo, double *subhi, tagint *tag,
                        int **nspecial, tagint **special, const bool eflag,
                        const bool vflag, const bool eatom, const bool vatom,
                        int &host_start, int **ilist, int **jnum,
                        const double cpu_time, bool &success, double *host_q,
                        double *boxlo, double *prd);
void cl_gpu_compute(const int ago, const int inum, const int nall,
                    double **host_x, int *host_type, int *ilist, int *numj,
                    int **firstneigh, const bool eflag, const bool vflag,
                    const bool eatom, const bool vatom, int &host_start,
                    const double cpu_time, bool &success, double *host_q,
                    const int nlocal, double *boxlo, double *prd);
double cl_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairCoulLongGPU::PairCoulLongGPU(LAMMPS *lmp) :
  PairCoulLong(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  cpu_time = 0.0;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairCoulLongGPU::~PairCoulLongGPU()
{
  cl_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairCoulLongGPU::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = cl_gpu_compute_n(neighbor->ago, inum, nall, atom->x,
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
    cl_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type,
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

void PairCoulLongGPU::init_style()
{
  cut_respa = NULL;

  if (!atom->q_flag)
    error->all(FLERR,"Pair style coul/long/gpu requires atom attribute q");
  if (force->newton_pair)
    error->all(FLERR,"Cannot use newton pair with coul/long/gpu pair style");

  // Call init_one calculation make sure scale is correct
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        init_one(i,j);
      }
    }
  }
  double cell_size = cut_coul + neighbor->skin;

  cut_coulsq = cut_coul * cut_coul;

  // insure use of KSpace long-range solver, set g_ewald

  if (force->kspace == NULL)
    error->all(FLERR,"Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables

  if (ncoultablebits) init_tables(cut_coul,cut_respa);

  int maxspecial=0;
  if (atom->molecular)
    maxspecial=atom->maxspecial;
  int success = cl_gpu_init(atom->ntypes+1, scale,
                            atom->nlocal, atom->nlocal+atom->nghost, 300,
                            maxspecial, cell_size, gpu_mode, screen, cut_coulsq,
                            force->special_coul, force->qqrd2e, g_ewald);

  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE) {
    int irequest = neighbor->request(this,instance_me);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulLongGPU::reinit()
{
  Pair::reinit();

  cl_gpu_reinit(atom->ntypes+1, scale);
}

/* ---------------------------------------------------------------------- */

double PairCoulLongGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + cl_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairCoulLongGPU::cpu_compute(int start, int inum, int eflag,
                                  int /* vflag */, int *ilist,
                                  int *numneigh, int **firstneigh)
{
  int i,j,ii,jj,jnum,itable;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double fraction,table;
  double r,r2inv,forcecoul,factor_coul;
  double grij,expm2,prefactor,t,erfc;
  int *jlist;
  double rsq;

  ecoul = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
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

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      r2inv = 1.0/rsq;

      if (rsq < cut_coulsq) {
        if (!ncoultablebits || rsq <= tabinnersq) {
          r = sqrt(rsq);
          grij = g_ewald * r;
          expm2 = exp(-grij*grij);
          t = 1.0 / (1.0 + EWALD_P*grij);
          erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
          prefactor = qqrd2e * qtmp*q[j]/r;
          forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
          if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
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

        fpair = forcecoul * r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;

        if (eflag) {
          if (rsq < cut_coulsq) {
            if (!ncoultablebits || rsq <= tabinnersq)
              ecoul = prefactor*erfc;
            else {
              table = etable[itable] + fraction*detable[itable];
              ecoul = qtmp*q[j] * table;
            }
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;
        }

        if (evflag) ev_tally_full(i,0.0,ecoul,fpair,delx,dely,delz);
      }
    }
  }
}
