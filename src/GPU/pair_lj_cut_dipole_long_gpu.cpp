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

#include "pair_lj_cut_dipole_long_gpu.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
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
#include "gpu_extra.h"

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

using namespace LAMMPS_NS;
using namespace MathConst;

// External functions from cuda library for atom decomposition

int dplj_gpu_init(const int ntypes, double **cutsq, double **host_lj1,
                 double **host_lj2, double **host_lj3, double **host_lj4,
                 double **offset, double *special_lj, const int nlocal,
                 const int nall, const int max_nbors, const int maxspecial,
                 const double cell_size, int &gpu_mode, FILE *screen,
                 double **host_cut_ljsq, const double host_cut_coulsq,
                 double *host_special_coul, const double qqrd2e, const double g_ewald);
void dplj_gpu_clear();
int ** dplj_gpu_compute_n(const int ago, const int inum,
                         const int nall, double **host_x, int *host_type,
                         double *sublo, double *subhi, tagint *tag,
                         int **nspecial, tagint **special, const bool eflag,
                         const bool vflag, const bool eatom, const bool vatom,
                         int &host_start, int **ilist, int **jnum,
                         const double cpu_time, bool &success,
                         double *host_q, double **host_mu,
                         double *boxlo, double *prd);
void dplj_gpu_compute(const int ago, const int inum,
                     const int nall, double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh,
                     const bool eflag, const bool vflag, const bool eatom,
                     const bool vatom, int &host_start, const double cpu_time,
                     bool &success, double *host_q, double **host_mu,
                     const int nlocal, double *boxlo, double *prd);
double dplj_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairLJCutDipoleLongGPU::PairLJCutDipoleLongGPU(LAMMPS *lmp) : PairLJCutDipoleLong(lmp),
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

PairLJCutDipoleLongGPU::~PairLJCutDipoleLongGPU()
{
  dplj_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairLJCutDipoleLongGPU::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = dplj_gpu_compute_n(neighbor->ago, inum, nall, atom->x,
                                   atom->type, domain->sublo, domain->subhi,
                                   atom->tag, atom->nspecial, atom->special,
                                   eflag, vflag, eflag_atom, vflag_atom,
                                   host_start, &ilist, &numneigh, cpu_time,
                                   success, atom->q, atom->mu, domain->boxlo,
                                   domain->prd);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    dplj_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type,
                    ilist, numneigh, firstneigh, eflag, vflag, eflag_atom,
                    vflag_atom, host_start, cpu_time, success, atom->q,
                    atom->mu, atom->nlocal, domain->boxlo, domain->prd);
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

void PairLJCutDipoleLongGPU::init_style()
{
  if (!atom->q_flag || !atom->mu_flag || !atom->torque_flag)
    error->all(FLERR,"Pair dipole/cut/gpu requires atom attributes q, mu, torque");

  if (force->newton_pair)
    error->all(FLERR,"Cannot use newton pair with dipole/cut/gpu pair style");

  if (strcmp(update->unit_style,"electron") == 0)
    error->all(FLERR,"Cannot (yet) use 'electron' units with dipoles");

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

  // setup force tables

  if (ncoultablebits) init_tables(cut_coul,NULL);

  int maxspecial=0;
  if (atom->molecular)
    maxspecial=atom->maxspecial;
  int success = dplj_gpu_init(atom->ntypes+1, cutsq, lj1, lj2, lj3, lj4,
                             offset, force->special_lj, atom->nlocal,
                             atom->nlocal+atom->nghost, 300, maxspecial,
                             cell_size, gpu_mode, screen, cut_ljsq, cut_coulsq,
                             force->special_coul, force->qqrd2e, g_ewald);
  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE) {
    int irequest = neighbor->request(this,instance_me);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJCutDipoleLongGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + dplj_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairLJCutDipoleLongGPU::cpu_compute(int start, int inum, int eflag, int vflag,
                                   int *ilist, int *numneigh,
                                   int **firstneigh)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,r,rinv,r2inv,r6inv;
  double forcecoulx,forcecouly,forcecoulz,fforce;
  double tixcoul,tiycoul,tizcoul;
  double fx,fy,fz,fdx,fdy,fdz,fax,fay,faz;
  double pdotp,pidotr,pjdotr,pre1,pre2,pre3;
  double grij,expm2,t,erfc;
  double g0,g1,g2,b0,b1,b2,b3,d0,d1,d2,d3;
  double zdix,zdiy,zdiz,zdjx,zdjy,zdjz,zaix,zaiy,zaiz,zajx,zajy,zajz;
  double g0b1_g1b2_g2b3,g0d1_g1d2_g2d3;
  double forcelj,factor_coul,factor_lj,facm1;
  double evdwl,ecoul;
  int *jlist;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  double **mu = atom->mu;
  double **torque = atom->torque;
  int *type = atom->type;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  double qqrd2e = force->qqrd2e;

  pre1 = 2.0 * g_ewald / MY_PIS;
  pre2 = 4.0 * pow(g_ewald,3.0) / MY_PIS;
  pre3 = 8.0 * pow(g_ewald,5.0) / MY_PIS;

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
        r2inv = 1.0/rsq;
        rinv = sqrt(r2inv);

        if (rsq < cut_coulsq) {
          r = sqrt(rsq);
          grij = g_ewald * r;
          expm2 = exp(-grij*grij);
          t = 1.0 / (1.0 + EWALD_P*grij);
          erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

          pdotp = mu[i][0]*mu[j][0] + mu[i][1]*mu[j][1] + mu[i][2]*mu[j][2];
          pidotr = mu[i][0]*delx + mu[i][1]*dely + mu[i][2]*delz;
          pjdotr = mu[j][0]*delx + mu[j][1]*dely + mu[j][2]*delz;

          g0 = qtmp*q[j];
          g1 = qtmp*pjdotr - q[j]*pidotr + pdotp;
          g2 = -pidotr*pjdotr;

          if (factor_coul > 0.0) {
            b0 = erfc * rinv;
            b1 = (b0 + pre1*expm2) * r2inv;
            b2 = (3.0*b1 + pre2*expm2) * r2inv;
            b3 = (5.0*b2 + pre3*expm2) * r2inv;

            g0b1_g1b2_g2b3 = g0*b1 + g1*b2 + g2*b3;
            fdx = delx * g0b1_g1b2_g2b3 -
              b1 * (qtmp*mu[j][0] - q[j]*mu[i][0]) +
              b2 * (pjdotr*mu[i][0] + pidotr*mu[j][0]);
            fdy = dely * g0b1_g1b2_g2b3 -
              b1 * (qtmp*mu[j][1] - q[j]*mu[i][1]) +
              b2 * (pjdotr*mu[i][1] + pidotr*mu[j][1]);
            fdz = delz * g0b1_g1b2_g2b3 -
              b1 * (qtmp*mu[j][2] - q[j]*mu[i][2]) +
              b2 * (pjdotr*mu[i][2] + pidotr*mu[j][2]);

            zdix = delx * (q[j]*b1 + b2*pjdotr) - b1*mu[j][0];
            zdiy = dely * (q[j]*b1 + b2*pjdotr) - b1*mu[j][1];
            zdiz = delz * (q[j]*b1 + b2*pjdotr) - b1*mu[j][2];
            zdjx = delx * (-qtmp*b1 + b2*pidotr) - b1*mu[i][0];
            zdjy = dely * (-qtmp*b1 + b2*pidotr) - b1*mu[i][1];
            zdjz = delz * (-qtmp*b1 + b2*pidotr) - b1*mu[i][2];

            if (factor_coul < 1.0) {
              fdx *= factor_coul;
              fdy *= factor_coul;
              fdz *= factor_coul;
              zdix *= factor_coul;
              zdiy *= factor_coul;
              zdiz *= factor_coul;
              zdjx *= factor_coul;
              zdjy *= factor_coul;
              zdjz *= factor_coul;
            }
          } else {
            fdx = fdy = fdz = 0.0;
            zdix = zdiy = zdiz = 0.0;
            zdjx = zdjy = zdjz = 0.0;
          }

          if (factor_coul < 1.0) {
            d0 = (erfc - 1.0) * rinv;
            d1 = (d0 + pre1*expm2) * r2inv;
            d2 = (3.0*d1 + pre2*expm2) * r2inv;
            d3 = (5.0*d2 + pre3*expm2) * r2inv;

            g0d1_g1d2_g2d3 = g0*d1 + g1*d2 + g2*d3;
            fax = delx * g0d1_g1d2_g2d3 -
              d1 * (qtmp*mu[j][0] - q[j]*mu[i][0]) +
              d2 * (pjdotr*mu[i][0] + pidotr*mu[j][0]);
            fay = dely * g0d1_g1d2_g2d3 -
              d1 * (qtmp*mu[j][1] - q[j]*mu[i][1]) +
              d2 * (pjdotr*mu[i][1] + pidotr*mu[j][1]);
            faz = delz * g0d1_g1d2_g2d3 -
              d1 * (qtmp*mu[j][2] - q[j]*mu[i][2]) +
              d2 * (pjdotr*mu[i][2] + pidotr*mu[j][2]);

            zaix = delx * (q[j]*d1 + d2*pjdotr) - d1*mu[j][0];
            zaiy = dely * (q[j]*d1 + d2*pjdotr) - d1*mu[j][1];
            zaiz = delz * (q[j]*d1 + d2*pjdotr) - d1*mu[j][2];
            zajx = delx * (-qtmp*d1 + d2*pidotr) - d1*mu[i][0];
            zajy = dely * (-qtmp*d1 + d2*pidotr) - d1*mu[i][1];
            zajz = delz * (-qtmp*d1 + d2*pidotr) - d1*mu[i][2];

            if (factor_coul > 0.0) {
              facm1 = 1.0 - factor_coul;
              fax *= facm1;
              fay *= facm1;
              faz *= facm1;
              zaix *= facm1;
              zaiy *= facm1;
              zaiz *= facm1;
              zajx *= facm1;
              zajy *= facm1;
              zajz *= facm1;
            }
          } else {
            fax = fay = faz = 0.0;
            zaix = zaiy = zaiz = 0.0;
            zajx = zajy = zajz = 0.0;
          }

          forcecoulx = fdx + fax;
          forcecouly = fdy + fay;
          forcecoulz = fdz + faz;

          tixcoul = mu[i][1]*(zdiz + zaiz) - mu[i][2]*(zdiy + zaiy);
          tiycoul = mu[i][2]*(zdix + zaix) - mu[i][0]*(zdiz + zaiz);
          tizcoul = mu[i][0]*(zdiy + zaiy) - mu[i][1]*(zdix + zaix);
        } else {
          forcecoulx = forcecouly = forcecoulz = 0.0;
          tixcoul = tiycoul = tizcoul = 0.0;
        }

        // LJ interaction

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          fforce = factor_lj * forcelj*r2inv;
        } else fforce = 0.0;

        // total force

        fx = qqrd2e*forcecoulx + delx*fforce;
        fy = qqrd2e*forcecouly + dely*fforce;
        fz = qqrd2e*forcecoulz + delz*fforce;

        // force & torque accumulation

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        torque[i][0] += qqrd2e*tixcoul;
        torque[i][1] += qqrd2e*tiycoul;
        torque[i][2] += qqrd2e*tizcoul;

        if (eflag) {
          if (rsq < cut_coulsq && factor_coul > 0.0) {
            ecoul = qqrd2e*(b0*g0 + b1*g1 + b2*g2);
            if (factor_coul < 1.0) {
              ecoul *= factor_coul;
              ecoul += (1-factor_coul) * qqrd2e * (d0*g0 + d1*g1 + d2*g2);
            }
          } else ecoul = 0.0;

          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
              offset[itype][jtype];
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (evflag) ev_tally_xyz_full(i,evdwl,ecoul,fx,fy,fz,delx,dely,delz);
      }
    }
  }
}
