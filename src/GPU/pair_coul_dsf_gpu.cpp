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
   Contributing author: Trung Dac Nguyen (ORNL)
------------------------------------------------------------------------- */

#include "pair_coul_dsf_gpu.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "gpu_extra.h"
#include "math_const.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace EwaldConst;
using MathConst::MY_PIS;

// External functions from cuda library for atom decomposition

int cdsf_gpu_init(const int ntypes, const int nlocal, const int nall, const int max_nbors,
                  const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                  const double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                  const double e_shift, const double f_shift, const double alpha);
void cdsf_gpu_clear();
int **cdsf_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int &host_start, int **ilist, int **jnum,
                         const double cpu_time, bool &success, double *host_q, double *boxlo,
                         double *prd);
void cdsf_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, int &host_start,
                      const double cpu_time, bool &success, double *host_q, const int nlocal,
                      double *boxlo, double *prd);
double cdsf_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairCoulDSFGPU::PairCoulDSFGPU(LAMMPS *lmp) : PairCoulDSF(lmp), gpu_mode(GPU_FORCE)
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

PairCoulDSFGPU::~PairCoulDSFGPU()
{
  cdsf_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairCoulDSFGPU::compute(int eflag, int vflag)
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
    firstneigh = cdsf_gpu_compute_n(neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi,
                                    atom->tag, atom->nspecial, atom->special, eflag, vflag,
                                    eflag_atom, vflag_atom, host_start, &ilist, &numneigh, cpu_time,
                                    success, atom->q, domain->boxlo, domain->prd);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    cdsf_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh, firstneigh,
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

void PairCoulDSFGPU::init_style()
{
  if (!atom->q_flag) error->all(FLERR, "Pair style coul/dsf/gpu requires atom attribute q");

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
  double erfcc = erfc(alpha * cut_coul);
  double erfcd = exp(-alpha * alpha * cut_coul * cut_coul);
  f_shift = -(erfcc / cut_coulsq + 2.0 / MY_PIS * alpha * erfcd / cut_coul);
  e_shift = erfcc / cut_coul - f_shift * cut_coul;

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;
  int mnf = 5e-2 * neighbor->oneatom;
  int success = cdsf_gpu_init(atom->ntypes + 1, atom->nlocal, atom->nlocal + atom->nghost, mnf,
                              maxspecial, cell_size, gpu_mode, screen, cut_coulsq,
                              force->special_coul, force->qqrd2e, e_shift, f_shift, alpha);
  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

double PairCoulDSFGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + cdsf_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairCoulDSFGPU::cpu_compute(int start, int inum, int eflag, int /* vflag */, int *ilist,
                                 int *numneigh, int **firstneigh)
{
  int i, j, ii, jj, jnum;
  double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, ecoul, fpair;
  double r, rsq, r2inv, forcecoul, factor_coul;
  double prefactor, erfcc, erfcd, t;
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
      double e_self = -(e_shift / 2.0 + alpha / MY_PIS) * qtmp * qtmp * qqrd2e;
      ev_tally(i, i, nlocal, 0, 0.0, e_self, 0.0, 0.0, 0.0, 0.0);
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cut_coulsq) {
        r2inv = 1.0 / rsq;
        r = sqrt(rsq);
        prefactor = qqrd2e * qtmp * q[j] / r;
        erfcd = exp(-alpha * alpha * r * r);
        t = 1.0 / (1.0 + EWALD_P * alpha * r);
        erfcc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * erfcd;
        forcecoul = prefactor * (erfcc / r + 2.0 * alpha / MY_PIS * erfcd + r * f_shift) * r;
        if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;

        fpair = forcecoul * r2inv;
        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;

        if (eflag) {
          if (rsq < cut_coulsq) {
            ecoul = prefactor * (erfcc - r * e_shift - rsq * f_shift);
            if (factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
          } else
            ecoul = 0.0;
        }

        if (evflag) ev_tally_full(i, 0.0, ecoul, fpair, delx, dely, delz);
      }
    }
  }
}
