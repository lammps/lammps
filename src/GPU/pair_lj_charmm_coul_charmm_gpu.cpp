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

#include "pair_lj_charmm_coul_charmm_gpu.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int crm_gpu_init(const int ntypes, double cut_bothsq, double **host_lj1, double **host_lj2,
                 double **host_lj3, double **host_lj4, double *special_lj, const int nlocal,
                 const int nall, const int max_nbors, const int maxspecial, const double cell_size,
                 int &gpu_mode, FILE *screen, double host_cut_ljsq, double host_cut_coulsq,
                 double *host_special_coul, const double qqrd2e, const double cut_lj_innersq,
                 const double cut_coul_innersq, const double denom_lj, const double denom_coul,
                 double **epsilon, double **sigma, const bool mix_arithmetic);
void crm_gpu_clear();
int **crm_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                        int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag, const bool eatom,
                        const bool vatom, int &host_start, int **ilist, int **jnum,
                        const double cpu_time, bool &success, double *host_q, double *boxlo,
                        double *prd);
void crm_gpu_compute(const int ago, const int inum, const int nall, double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start, const double cpu_time,
                     bool &success, double *host_q, const int nlocal, double *boxlo, double *prd);
double crm_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairLJCharmmCoulCharmmGPU::PairLJCharmmCoulCharmmGPU(LAMMPS *lmp) :
    PairLJCharmmCoulCharmm(lmp), gpu_mode(GPU_FORCE)
{
  reinitflag = 0;
  cpu_time = 0.0;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairLJCharmmCoulCharmmGPU::~PairLJCharmmCoulCharmmGPU()
{
  crm_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulCharmmGPU::compute(int eflag, int vflag)
{
  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    inum = atom->nlocal;
    firstneigh = crm_gpu_compute_n(neighbor->ago, inum, nall, atom->x, atom->type, domain->sublo,
                                   domain->subhi, atom->tag, atom->nspecial, atom->special, eflag,
                                   vflag, eflag_atom, vflag_atom, host_start, &ilist, &numneigh,
                                   cpu_time, success, atom->q, domain->boxlo, domain->prd);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    crm_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh, firstneigh,
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

void PairLJCharmmCoulCharmmGPU::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR, "Pair style lj/charmm/coul/long/gpu requires atom attribute q");

  // Repeated cutsq calculation in init_one() is required for GPU package

  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) init_one(i, j);
    }
  }

  cut_lj_innersq = cut_lj_inner * cut_lj_inner;
  cut_coul_innersq = cut_coul_inner * cut_coul_inner;
  cut_ljsq = cut_lj * cut_lj;
  cut_coulsq = cut_coul * cut_coul;
  cut_bothsq = MAX(cut_ljsq, cut_coulsq);

  denom_lj =
      (cut_ljsq - cut_lj_innersq) * (cut_ljsq - cut_lj_innersq) * (cut_ljsq - cut_lj_innersq);
  denom_lj = 1.0 / denom_lj;

  denom_coul = (cut_coulsq - cut_coul_innersq) * (cut_coulsq - cut_coul_innersq) *
      (cut_coulsq - cut_coul_innersq);
  denom_coul = 1.0 / denom_coul;

  double cell_size = sqrt(cut_bothsq) + neighbor->skin;

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
      crm_gpu_init(atom->ntypes + 1, cut_bothsq, lj1, lj2, lj3, lj4, force->special_lj,
                   atom->nlocal, atom->nlocal + atom->nghost, mnf, maxspecial, cell_size, gpu_mode,
                   screen, cut_ljsq, cut_coulsq, force->special_coul, force->qqrd2e, cut_lj_innersq,
                   cut_coul_innersq, denom_lj, denom_coul, epsilon, sigma, arithmetic);
  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

double PairLJCharmmCoulCharmmGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + crm_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulCharmmGPU::cpu_compute(int start, int inum, int eflag, int /* vflag */,
                                            int *ilist, int *numneigh, int **firstneigh)
{
  int i, j, ii, jj, jnum, itype, jtype;
  double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, evdwl, ecoul, fpair;
  double rsq, r2inv, r6inv, forcecoul, forcelj, factor_coul, factor_lj;
  double philj, switch1, switch2;
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
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cut_bothsq) {
        r2inv = 1.0 / rsq;

        if (rsq < cut_coulsq) {
          forcecoul = qqrd2e * qtmp * q[j] * sqrt(r2inv);
          if (rsq > cut_coul_innersq) {
            switch1 = (cut_coulsq - rsq) * (cut_coulsq - rsq) *
                (cut_coulsq + 2.0 * rsq - 3.0 * cut_coul_innersq) * denom_coul;
            forcecoul *= switch1;
          }
        } else
          forcecoul = 0.0;

        if (rsq < cut_ljsq) {
          r6inv = r2inv * r2inv * r2inv;
          jtype = type[j];
          forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
          if (rsq > cut_lj_innersq) {
            switch1 = (cut_ljsq - rsq) * (cut_ljsq - rsq) *
                (cut_ljsq + 2.0 * rsq - 3.0 * cut_lj_innersq) * denom_lj;
            switch2 = 12.0 * rsq * (cut_ljsq - rsq) * (rsq - cut_lj_innersq) * denom_lj;
            philj = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]);
            forcelj = forcelj * switch1 + philj * switch2;
          }
        } else
          forcelj = 0.0;

        fpair = (factor_coul * forcecoul + factor_lj * forcelj) * r2inv;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;

        if (eflag) {
          if (rsq < cut_coulsq) {
            ecoul = qqrd2e * qtmp * q[j] * sqrt(r2inv);
            if (rsq > cut_coul_innersq) {
              switch1 = (cut_coulsq - rsq) * (cut_coulsq - rsq) *
                  (cut_coulsq + 2.0 * rsq - 3.0 * cut_coul_innersq) * denom_coul;
              ecoul *= switch1;
            }
            ecoul *= factor_coul;
          } else
            ecoul = 0.0;

          if (rsq < cut_ljsq) {
            evdwl = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]);
            if (rsq > cut_lj_innersq) {
              switch1 = (cut_ljsq - rsq) * (cut_ljsq - rsq) *
                  (cut_ljsq + 2.0 * rsq - 3.0 * cut_lj_innersq) * denom_lj;
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
