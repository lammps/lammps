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

#include "pair_soft_gpu.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "math_const.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int soft_gpu_init(const int ntypes, double **cutsq, double **prefactor, double **cut,
                  double *special_lj, const int nlocal, const int nall, const int max_nbors,
                  const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen);
void soft_gpu_reinit(const int ntypes, double **cutsq, double **host_prefactor, double **host_cut);
void soft_gpu_clear();
int **soft_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int &host_start, int **ilist, int **jnum,
                         const double cpu_time, bool &success);
void soft_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, int &host_start,
                      const double cpu_time, bool &success);
double soft_gpu_bytes();

using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairSoftGPU::PairSoftGPU(LAMMPS *lmp) : PairSoft(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  cpu_time = 0.0;
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairSoftGPU::~PairSoftGPU()
{
  soft_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairSoftGPU::compute(int eflag, int vflag)
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
    firstneigh =
        soft_gpu_compute_n(neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi, atom->tag,
                           atom->nspecial, atom->special, eflag, vflag, eflag_atom, vflag_atom,
                           host_start, &ilist, &numneigh, cpu_time, success);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    soft_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh, firstneigh,
                     eflag, vflag, eflag_atom, vflag_atom, host_start, cpu_time, success);
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

void PairSoftGPU::init_style()
{

  // Repeat cutsq calculation because done after call to init_style
  double maxcut = -1.0;
  double mcut;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        mcut = init_one(i, j);
        mcut *= mcut;
        if (mcut > maxcut) maxcut = mcut;
        cutsq[i][j] = cutsq[j][i] = mcut;
      } else
        cutsq[i][j] = cutsq[j][i] = 0.0;
    }
  }
  double cell_size = sqrt(maxcut) + neighbor->skin;

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;
  int mnf = 5e-2 * neighbor->oneatom;
  int success =
      soft_gpu_init(atom->ntypes + 1, cutsq, prefactor, cut, force->special_lj, atom->nlocal,
                    atom->nlocal + atom->nghost, mnf, maxspecial, cell_size, gpu_mode, screen);
  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

void PairSoftGPU::reinit()
{
  Pair::reinit();

  soft_gpu_reinit(atom->ntypes + 1, cutsq, prefactor, cut);
}

/* ---------------------------------------------------------------------- */

double PairSoftGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + soft_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairSoftGPU::cpu_compute(int start, int inum, int eflag, int /* vflag */, int *ilist,
                              int *numneigh, int **firstneigh)
{
  int i, j, ii, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double r, rsq, arg, factor_lj;
  int *jlist;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  double *special_lj = force->special_lj;

  // loop over neighbors of my atoms

  for (ii = start; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        arg = MY_PI * r / cut[itype][jtype];
        if (r > 0.0)
          fpair = factor_lj * prefactor[itype][jtype] * sin(arg) * MY_PI / cut[itype][jtype] / r;
        else
          fpair = 0.0;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;

        if (eflag) evdwl = factor_lj * prefactor[itype][jtype] * (1.0 + cos(arg));

        if (evflag) ev_tally_full(i, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }
}
