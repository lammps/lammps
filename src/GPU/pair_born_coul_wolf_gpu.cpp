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
   Contributing authors: Trung Dac Nguyen (ORNL)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "pair_born_coul_wolf_gpu.h"
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
#include "string.h"
#include "gpu_extra.h"

using namespace LAMMPS_NS;
using namespace MathConst;

// External functions from cuda library for atom decomposition

int borncw_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv,
                    double **host_born1, double **host_born2, 
                    double **host_born3, double **host_a, double **host_c, 
                    double **host_d, double **sigma, double **offset, 
                    double *special_lj, const int inum,
                    const int nall, const int max_nbors, const int maxspecial,
                    const double cell_size, int &gpu_mode, FILE *screen,
                    double **host_cut_ljsq, double host_cut_coulsq,
                    double *host_special_coul, const double qqrd2e,
                    const double alf, const double e_shift, const double f_shift);
void borncw_gpu_clear();
int ** borncw_gpu_compute_n(const int ago, const int inum_full, const int nall,
                            double **host_x, int *host_type, double *sublo,
                            double *subhi, tagint *tag, int **nspecial, 
                            tagint **special, const bool eflag, const bool vflag,
                            const bool eatom, const bool vatom, int &host_start,
                            int **ilist, int **jnum, const double cpu_time,
                            bool &success, double *host_q, double *boxlo,
                            double *prd);
void borncw_gpu_compute(const int ago, const int inum_full, const int nall,
                        double **host_x, int *host_type, int *ilist, int *numj,
                        int **firstneigh, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, int &host_start,
                        const double cpu_time, bool &success, double *host_q,
                        const int nlocal, double *boxlo, double *prd);
double borncw_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairBornCoulWolfGPU::PairBornCoulWolfGPU(LAMMPS *lmp) : PairBornCoulWolf(lmp), 
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

PairBornCoulWolfGPU::~PairBornCoulWolfGPU()
{
  borncw_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairBornCoulWolfGPU::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;
  
  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;
  
  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = borncw_gpu_compute_n(neighbor->ago, inum, nall,
                                      atom->x, atom->type, domain->sublo,
                                      domain->subhi, atom->tag, atom->nspecial,
                                      atom->special, eflag, vflag, eflag_atom,
                                      vflag_atom, host_start, 
                                      &ilist, &numneigh, cpu_time, success, 
                                      atom->q, domain->boxlo, domain->prd);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    borncw_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type,
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

void PairBornCoulWolfGPU::init_style()
{
  if (force->newton_pair) 
    error->all(FLERR,
      "Cannot use newton pair with born/coul/wolf/gpu pair style");

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

  double e_shift = erfc(alf*cut_coul)/cut_coul;
  double f_shift = -(e_shift+ 2.0*alf/MY_PIS * exp(-alf*alf*cut_coul*cut_coul)) / 
    cut_coul; 

  int maxspecial=0;
  if (atom->molecular)
    maxspecial=atom->maxspecial;
  int success = borncw_gpu_init(atom->ntypes+1, cutsq, rhoinv, 
                                born1, born2, born3, a, c, d, sigma, offset, 
                                force->special_lj, atom->nlocal,
                                atom->nlocal+atom->nghost, 300, maxspecial,
                                cell_size, gpu_mode, screen, cut_ljsq,
                                cut_coulsq, force->special_coul, force->qqrd2e, 
                                alf, e_shift, f_shift);
  GPU_EXTRA::check_flag(success,error,world);

  if (gpu_mode == GPU_FORCE) {
    int irequest = neighbor->request(this);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
  } 
}

/* ---------------------------------------------------------------------- */

double PairBornCoulWolfGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + borncw_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairBornCoulWolfGPU::cpu_compute(int start, int inum, int eflag, int vflag,
                                      int *ilist, int *numneigh,
                                      int **firstneigh) {
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,qtmp,delx,dely,delz,evdwl,ecoul,fpair;
  double rsq,r2inv,r6inv,forcecoul,forceborn,factor_coul,factor_lj;
  double erfcc,erfcd,v_sh,dvdrr,e_self,qisq;
  double prefactor;
  double r,rexp;
  int *jlist;

  evdwl = ecoul = 0.0;
  
  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  double qqrd2e = force->qqrd2e;
  
  double e_shift = erfc(alf*cut_coul)/cut_coul;
  double f_shift = -(e_shift+ 2.0*alf/MY_PIS * exp(-alf*alf*cut_coul*cut_coul)) / 
    cut_coul; 

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

    qisq = qtmp*qtmp;
    e_self = -(e_shift/2.0 + alf/MY_PIS) * qisq*qqrd2e;
    if (evflag) ev_tally(i,i,nlocal,0,0.0,e_self,0.0,0.0,0.0,0.0);

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
  
        if (rsq < cut_coulsq) {
          r = sqrt(rsq);
          prefactor = qqrd2e*qtmp*q[j]/r;
          erfcc = erfc(alf*r); 
          erfcd = exp(-alf*alf*r*r);
          v_sh = (erfcc - e_shift*r) * prefactor; 
          dvdrr = (erfcc/rsq + 2.0*alf/MY_PIS * erfcd/r) + f_shift;
          forcecoul = dvdrr*rsq*prefactor;
          if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
        } else forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          r = sqrt(rsq);
          rexp = exp((sigma[itype][jtype]-r)*rhoinv[itype][jtype]);
          forceborn = born1[itype][jtype]*r*rexp - born2[itype][jtype]*r6inv + 
            born3[itype][jtype]*r2inv*r6inv;
        } else forceborn = 0.0;
	
        fpair = (factor_coul*forcecoul + factor_lj*forceborn) * r2inv;
	
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;

        if (eflag) {
          if (rsq < cut_coulsq) {
            ecoul = v_sh;
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;
          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = a[itype][jtype]*rexp - c[itype][jtype]*r6inv +
              d[itype][jtype]*r6inv*r2inv - offset[itype][jtype];
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (evflag) ev_tally_full(i,evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
  }
}
