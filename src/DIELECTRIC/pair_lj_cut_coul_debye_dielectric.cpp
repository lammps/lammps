/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#include "pair_lj_cut_coul_debye_dielectric.h"

#include "atom.h"
#include "atom_vec_dielectric.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;
using MathConst::MY_PIS;

static constexpr double EPSILON = 1.0e-6;

/* ---------------------------------------------------------------------- */

PairLJCutCoulDebyeDielectric::PairLJCutCoulDebyeDielectric(LAMMPS *_lmp) : PairLJCutCoulDebye(_lmp)
{
  efield = nullptr;
  epot = nullptr;
  nmax = 0;
  no_virial_fdotr_compute = 1;
}

/* ---------------------------------------------------------------------- */

PairLJCutCoulDebyeDielectric::~PairLJCutCoulDebyeDielectric()
{
  memory->destroy(efield);
  memory->destroy(epot);
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulDebyeDielectric::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double qtmp, etmp, xtmp, ytmp, ztmp, delx, dely, delz, evdwl, ecoul;
  double fpair_i;
  double rsq, r2inv, r6inv, forcecoul, forcelj, factor_coul, factor_lj, efield_i, epot_i;
  double r, rinv, screening;
  int *ilist, *jlist, *numneigh, **firstneigh;

  if (atom->nmax > nmax) {
    memory->destroy(efield);
    memory->destroy(epot);
    nmax = atom->nmax;
    memory->create(efield, nmax, 3, "pair:efield");
    memory->create(epot, nmax, "pair:epot");
  }

  evdwl = ecoul = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q_scaled;
  double *eps = atom->epsilon;
  double **norm = atom->mu;
  double *curvature = atom->curvature;
  double *area = atom->area;
  int *type = atom->type;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    etmp = eps[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // self term Eq. (55) for I_{ii} and Eq. (52) and in Barros et al
    double curvature_threshold = sqrt(area[i]);
    if (curvature[i] < curvature_threshold) {
      double sf = curvature[i] / (4.0 * MY_PIS * curvature_threshold) * area[i] * q[i];
      efield[i][0] = sf * norm[i][0];
      efield[i][1] = sf * norm[i][1];
      efield[i][2] = sf * norm[i][2];
    } else {
      efield[i][0] = efield[i][1] = efield[i][2] = 0;
    }

    epot[i] = 0;

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

        if (rsq < cut_coulsq[itype][jtype] && rsq > EPSILON) {
          r = sqrt(rsq);
          rinv = 1.0 / r;
          screening = exp(-kappa * r);
          efield_i = qqrd2e * q[j] * screening * (kappa + rinv);
          forcecoul = qtmp * efield_i;
          epot_i = efield_i;
        } else
          epot_i = efield_i = forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv * r2inv * r2inv;
          forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
        } else
          forcelj = 0.0;

        fpair_i = (factor_coul * etmp * forcecoul + factor_lj * forcelj) * r2inv;
        f[i][0] += delx * fpair_i;
        f[i][1] += dely * fpair_i;
        f[i][2] += delz * fpair_i;

        efield_i *= (factor_coul * etmp * r2inv);
        efield[i][0] += delx * efield_i;
        efield[i][1] += dely * efield_i;
        efield[i][2] += delz * efield_i;
        epot[i] += epot_i;

        if (eflag) {
          if (rsq < cut_coulsq[itype][jtype]) {
            ecoul = factor_coul * qqrd2e * qtmp * q[j] * 0.5 * (etmp + eps[j]) * rinv * screening;
          } else
            ecoul = 0.0;
          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];
            evdwl *= factor_lj;
          } else
            evdwl = 0.0;
        }

        if (evflag) ev_tally_full(i, evdwl, ecoul, fpair_i, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutCoulDebyeDielectric::init_style()
{
  avec = dynamic_cast<AtomVecDielectric *>(atom->style_match("dielectric"));
  if (!avec) error->all(FLERR, "Pair lj/cut/coul/debye/dielectric requires atom style dielectric");

  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

double PairLJCutCoulDebyeDielectric::single(int i, int j, int itype, int jtype, double rsq,
                                            double factor_coul, double factor_lj, double &fforce)
{
  double r2inv, r6inv, forcecoul, forcelj, phicoul, ei, ej, philj;
  double r, rinv, screening;
  double *q = atom->q_scaled;
  double *eps = atom->epsilon;

  r2inv = 1.0 / rsq;
  if (rsq < cut_coulsq[itype][jtype]) {
    r = sqrt(rsq);
    rinv = 1.0 / r;
    screening = exp(-kappa * r);
    forcecoul = force->qqrd2e * q[i] * q[j] * screening * (kappa + rinv) * eps[i];
  } else
    forcecoul = 0.0;
  if (rsq < cut_ljsq[itype][jtype]) {
    r6inv = r2inv * r2inv * r2inv;
    forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
  } else
    forcelj = 0.0;
  fforce = (factor_coul * forcecoul + factor_lj * forcelj) * r2inv;

  double eng = 0.0;
  if (eps[i] == 1)
    ei = 0;
  else
    ei = eps[i];
  if (eps[j] == 1)
    ej = 0;
  else
    ej = eps[j];
  if (rsq < cut_coulsq[itype][jtype]) {
    phicoul = force->qqrd2e * q[i] * q[j] * rinv * screening;
    phicoul *= 0.5 * (ei + ej);
    eng += factor_coul * phicoul;
  }
  if (rsq < cut_ljsq[itype][jtype]) {
    philj = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];
    eng += factor_lj * philj;
  }

  return eng;
}
