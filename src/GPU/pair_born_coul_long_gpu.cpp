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

#include "pair_born_coul_long_gpu.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "gpu_extra.h"
#include "kspace.h"
#include "math_const.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace EwaldConst;

// External functions from cuda library for atom decomposition

int borncl_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv, double **host_born1,
                    double **host_born2, double **host_born3, double **host_a, double **host_c,
                    double **host_d, double **sigma, double **offset, double *special_lj,
                    const int inum, const int nall, const int max_nbors, const int maxspecial,
                    const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                    double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                    const double g_ewald);
void borncl_gpu_clear();
int **borncl_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                           int *host_type, double *sublo, double *subhi, tagint *tag,
                           int **nspecial, tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int &host_start, int **ilist,
                           int **jnum, const double cpu_time, bool &success, double *host_q,
                           double *boxlo, double *prd);
void borncl_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                        int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                        const bool vflag, const bool eatom, const bool vatom, int &host_start,
                        const double cpu_time, bool &success, double *host_q, const int nlocal,
                        double *boxlo, double *prd);
double borncl_gpu_bytes();

/* ---------------------------------------------------------------------- */

PairBornCoulLongGPU::PairBornCoulLongGPU(LAMMPS *lmp) : PairBornCoulLong(lmp), gpu_mode(GPU_FORCE)
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

PairBornCoulLongGPU::~PairBornCoulLongGPU()
{
  borncl_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairBornCoulLongGPU::compute(int eflag, int vflag)
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
    firstneigh = borncl_gpu_compute_n(neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi,
                                      atom->tag, atom->nspecial, atom->special, eflag, vflag,
                                      eflag_atom, vflag_atom, host_start, &ilist, &numneigh,
                                      cpu_time, success, atom->q, domain->boxlo, domain->prd);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    borncl_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh, firstneigh,
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

void PairBornCoulLongGPU::init_style()
{
  if (!atom->q_flag) error->all(FLERR, "Pair style born/coul/long/gpu requires atom attribute q");

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

  // ensure use of KSpace long-range solver, set g_ewald

  if (force->kspace == nullptr) error->all(FLERR, "Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables

  if (ncoultablebits) init_tables(cut_coul, cut_respa);

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;
  int mnf = 5e-2 * neighbor->oneatom;
  int success = borncl_gpu_init(
      atom->ntypes + 1, cutsq, rhoinv, born1, born2, born3, a, c, d, sigma, offset,
      force->special_lj, atom->nlocal, atom->nlocal + atom->nghost, mnf, maxspecial, cell_size,
      gpu_mode, screen, cut_ljsq, cut_coulsq, force->special_coul, force->qqrd2e, g_ewald);

  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

double PairBornCoulLongGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + borncl_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairBornCoulLongGPU::cpu_compute(int start, int inum, int eflag, int /* vflag */, int *ilist,
                                      int *numneigh, int **firstneigh)
{
  int i, j, ii, jj, jnum, itype, jtype;
  double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, evdwl, ecoul, fpair;
  double r, rexp, r2inv, r6inv, forcecoul, forceborn, factor_coul, factor_lj;
  double grij, expm2, prefactor, t, erfc;
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
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0 / rsq;
        r = sqrt(rsq);

        if (rsq < cut_coulsq) {
          grij = g_ewald * r;
          expm2 = exp(-grij * grij);
          t = 1.0 / (1.0 + EWALD_P * grij);
          erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
          prefactor = qqrd2e * qtmp * q[j] / r;
          forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);
          if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;
        } else
          forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv * r2inv * r2inv;
          rexp = exp((sigma[itype][jtype] - r) * rhoinv[itype][jtype]);
          forceborn = born1[itype][jtype] * r * rexp - born2[itype][jtype] * r6inv +
              born3[itype][jtype] * r2inv * r6inv;
        } else
          forceborn = 0.0;

        fpair = (forcecoul + factor_lj * forceborn) * r2inv;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;

        if (eflag) {
          if (rsq < cut_coulsq) {
            ecoul = prefactor * erfc;
            if (factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
          } else
            ecoul = 0.0;
          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = a[itype][jtype] * rexp - c[itype][jtype] * r6inv +
                d[itype][jtype] * r6inv * r2inv - offset[itype][jtype];
            evdwl *= factor_lj;
          } else
            evdwl = 0.0;
        }

        if (evflag) ev_tally_full(i, evdwl, ecoul, fpair, delx, dely, delz);
      }
    }
  }
}
