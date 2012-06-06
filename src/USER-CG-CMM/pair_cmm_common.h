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
   Common functionality for the CMM coarse grained MD potentials.
   Contributing author: Axel Kohlmeyer <akohlmey@gmail.com>
------------------------------------------------------------------------- */

#ifndef LMP_PAIR_CMM_COMMON_H
#define LMP_PAIR_CMM_COMMON_H

#include "pair.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "respa.h"
#include "update.h"

#include "cg_cmm_parms.h"

#include "math.h"

namespace LAMMPS_NS {

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define EWALD_A1  0.254829592
#define EWALD_A2 -0.284496736
#define EWALD_A3  1.421413741
#define EWALD_A4 -1.453152027
#define EWALD_A5  1.061405429

  class PairCMMCommon : public Pair , public CGCMMParms {
    public:

    PairCMMCommon(class LAMMPS *);
    virtual ~PairCMMCommon();

    virtual void settings(int, char **);
    virtual void coeff(int, char **);
    virtual void init_style();
    virtual void init_list(int, class NeighList *);
    virtual double init_one(int, int);

    virtual void write_restart(FILE *);
    virtual void read_restart(FILE *);
    virtual void write_restart_settings(FILE *);
    virtual void read_restart_settings(FILE *);

    virtual double memory_usage();

    protected:

    // coarse grain flags
    int **cg_type;

    // lennard jones parameters
    double cut_lj_global, **cut, **cut_lj, **cut_ljsq;
    double **epsilon, **sigma;
    double **lj1, **lj2, **lj3, **lj4, **offset;

    // coulomb parameters
    int allocated_coul; // 0/1 = whether coulomb arrays are allocated
    double cut_coul_global, cut_coulsq_global, kappa, g_ewald;
    double **cut_coul, **cut_coulsq;

    // tables
    double tabinnersq;
    double *rtable,*drtable,*ftable,*dftable,*ctable,*dctable;
    double *etable,*detable,*ptable,*dptable,*vtable,*dvtable;
    int ncoulshiftbits,ncoulmask;

    // r-RESPA parameters
    double *cut_respa;

    // methods
    virtual void allocate();

    private:

    // disable default constructor
    PairCMMCommon();

    protected:
    // general optimizeable real space loops
    template < const int EVFLAG, const int EFLAG,
      const int NEWTON_PAIR, const int COUL_TYPE >
      void eval_verlet();
    template < const int NEWTON_PAIR, const int COUL_TYPE >
      void eval_inner();
    template < const int NEWTON_PAIR, const int COUL_TYPE >
      void eval_middle();
    template < const int EVFLAG, const int EFLAG, const int VFLAG,
      const int NEWTON_PAIR, const int COUL_TYPE >
      void eval_outer();

    // this one is not performance critical... no template needed.
    double eval_single(int, int, int, int, int,
                     double, double, double, double &);
  };

/* ---------------------------------------------------------------------- */
/* this is the inner heart of the CG potentials. */
#define CG_LJ_INNER(eflag,fvar)                                         \
  fvar=factor_lj;                                                       \
  if (eflag) evdwl=factor_lj;                                           \
                                                                        \
  if (cgt == CG_LJ12_4) {                                               \
    const double r4inv=r2inv*r2inv;                                     \
                                                                        \
    fvar *= r4inv*(lj1[itype][jtype]*r4inv*r4inv                        \
                    - lj2[itype][jtype]);                               \
                                                                        \
    if (eflag) {                                                        \
      evdwl *= r4inv*(lj3[itype][jtype]*r4inv*r4inv                     \
                      - lj4[itype][jtype]) - offset[itype][jtype];      \
    }                                                                   \
  } else if (cgt == CG_LJ9_6) {                                         \
    const double r3inv = r2inv*sqrt(r2inv);                             \
    const double r6inv = r3inv*r3inv;                                   \
    fvar *= r6inv*(lj1[itype][jtype]*r3inv                              \
                    - lj2[itype][jtype]);                               \
    if (eflag) {                                                        \
      evdwl *= r6inv*(lj3[itype][jtype]*r3inv                           \
                      - lj4[itype][jtype]) - offset[itype][jtype];      \
    }                                                                   \
  } else if (cgt == CG_LJ12_6) {                                        \
    const double r6inv = r2inv*r2inv*r2inv;                             \
    fvar *= r6inv*(lj1[itype][jtype]*r6inv                              \
                    - lj2[itype][jtype]);                               \
    if (eflag) {                                                        \
      evdwl *= r6inv*(lj3[itype][jtype]*r6inv                           \
                      - lj4[itype][jtype]) - offset[itype][jtype];      \
    }                                                                   \
  } else {                                                              \
    /* do nothing. this is a "cannot happen(TM)" case */                \
    ;                                                                   \
  }

#define CG_LJ_ENERGY(eflag)                                             \
  if (eflag) {                                                          \
    evdwl=factor_lj;                                                    \
                                                                        \
    if (cgt == CG_LJ12_4) {                                             \
      const double r4inv=r2inv*r2inv;                                   \
      evdwl *= r4inv*(lj3[itype][jtype]*r4inv*r4inv                     \
                      - lj4[itype][jtype]) - offset[itype][jtype];      \
    } else if (cgt == CG_LJ9_6) {                                       \
      const double r3inv = r2inv*sqrt(r2inv);                           \
      const double r6inv = r3inv*r3inv;                                 \
      evdwl *= r6inv*(lj3[itype][jtype]*r3inv                           \
                      - lj4[itype][jtype]) - offset[itype][jtype];      \
    } else if (cgt == CG_LJ12_6) {                                      \
      const double r6inv = r2inv*r2inv*r2inv;                           \
      evdwl *= r6inv*(lj3[itype][jtype]*r6inv                           \
                      - lj4[itype][jtype]) - offset[itype][jtype];      \
    } else {                                                            \
      /* do nothing. this is a "cannot happen(TM)" case */              \
      ;                                                                 \
    }                                                                   \
  }                                                                     \



  template < const int EVFLAG, const int EFLAG,
    const int NEWTON_PAIR, const int COUL_TYPE >
    void PairCMMCommon::eval_verlet()
  {
    double ** const x = atom->x;
    double ** const f = atom->f;
    const double  * const q = atom->q;
    const int * const type = atom->type;
    const int nlocal = atom->nlocal;
    const double * const special_lj = force->special_lj;
    const double * const special_coul = force->special_coul;
    const double qqrd2e = force->qqrd2e;
    double factor_lj,factor_coul;

    const int inum = list->inum;
    const int * const ilist = list->ilist;
    const int * const numneigh = list->numneigh;
    int * const * const firstneigh = list->firstneigh;

    // loop over neighbors of my atoms

    int ii,jj;
    for (ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      double qtmp = (COUL_TYPE != CG_COUL_NONE) ? q[i] : 0.0;

      const int itype = type[i];
      const int * const jlist = firstneigh[i];
      const int jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        int j2 = jlist[jj];
        factor_lj = special_lj[sbmask(j2)];
        factor_coul = special_coul[sbmask(j2)];
        const int j = j2 & NEIGHMASK;

        const double delx = xtmp - x[j][0];
        const double dely = ytmp - x[j][1];
        const double delz = ztmp - x[j][2];
        const double rsq = delx*delx + dely*dely + delz*delz;
        const int jtype = type[j];

        double evdwl = 0.0;
        double ecoul = 0.0;
        double fpair = 0.0;

        const double r2inv = 1.0/rsq;
        const int cgt=cg_type[itype][jtype];

        if (rsq < cutsq[itype][jtype]) {
          if (COUL_TYPE == CG_COUL_NONE) {
            CG_LJ_INNER(EFLAG,fpair);
            fpair *= r2inv;
          } else {
            double forcelj  = 0.0;
            double forcecoul = 0.0;

            if (rsq < cut_ljsq[itype][jtype]) {
              CG_LJ_INNER(EFLAG,forcelj);
            }

            // coulomb with cutoff and screening
            if ((COUL_TYPE == CG_COUL_CUT) || (COUL_TYPE == CG_COUL_DEBYE)) {
              if (rsq < cut_coulsq[itype][jtype]) {
                double r=sqrt(rsq);
                double qscreen=exp(-kappa*r);
                forcecoul = factor_coul * qqrd2e
                  * qtmp * q[j] * qscreen * (kappa + 1.0/r);
                if (EFLAG) ecoul=factor_coul*qqrd2e
                  * qtmp*q[j] * qscreen / r;
              }
            }

            if (COUL_TYPE == CG_COUL_LONG) {
              if (rsq < cut_coulsq_global) {
                if (!ncoultablebits || rsq <= tabinnersq) {
                  const double r = sqrt(rsq);

                  const double grij = g_ewald * r;
                  const double expm2 = exp(-grij*grij);
                  const double t = 1.0 / (1.0 + EWALD_P*grij);
                  const double erfc = t * (EWALD_A1+t*(EWALD_A2+t*(EWALD_A3+t*(EWALD_A4+t*EWALD_A5)))) * expm2;
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
                  const double fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
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
              }
            }
            fpair = (forcecoul + forcelj) * r2inv;
          }
          f[i][0] += delx*fpair;
          f[i][1] += dely*fpair;
          f[i][2] += delz*fpair;
          if (NEWTON_PAIR || j < nlocal) {
            f[j][0] -= delx*fpair;
            f[j][1] -= dely*fpair;
            f[j][2] -= delz*fpair;
          }
          if (EVFLAG) ev_tally(i,j,nlocal,NEWTON_PAIR,
                               evdwl,ecoul,fpair,delx,dely,delz);
        }
      }
    }
    if (vflag_fdotr) virial_fdotr_compute();
  }

/* ---------------------------------------------------------------------- */

  template < const int NEWTON_PAIR, const int COUL_TYPE >
    void PairCMMCommon::eval_inner()
  {
    double ** const x = atom->x;
    double ** const f = atom->f;
    const double  * const q = atom->q;
    const int * const type = atom->type;
    const int nlocal = atom->nlocal;
    const double * const special_lj = force->special_lj;
    const double * const special_coul = force->special_coul;
    const double qqrd2e = force->qqrd2e;
    double factor_lj,factor_coul;

    const int inum = listinner->inum;
    const int * const ilist = listinner->ilist;
    const int * const numneigh = listinner->numneigh;
    int * const * const firstneigh = listinner->firstneigh;

    const double cut_out_on = cut_respa[0];
    const double cut_out_off = cut_respa[1];

    const double cut_out_diff = cut_out_off - cut_out_on;
    const double cut_out_on_sq = cut_out_on*cut_out_on;
    const double cut_out_off_sq = cut_out_off*cut_out_off;

    // loop over neighbors of my atoms

    int ii,jj;
    for (ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      double qtmp = (COUL_TYPE != CG_COUL_NONE) ? q[i] : 0.0;
      const int itype = type[i];
      const int * const jlist = firstneigh[i];
      const int jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        int j2 = jlist[jj];
        factor_lj = special_lj[sbmask(j2)];
        factor_coul = special_coul[sbmask(j2)];
        const int j = j2 & NEIGHMASK;

        const double delx = xtmp - x[j][0];
        const double dely = ytmp - x[j][1];
        const double delz = ztmp - x[j][2];
        const double rsq = delx*delx + dely*dely + delz*delz;
        const int jtype = type[j];

        double evdwl = 0.0;
        double ecoul = 0.0;
        double fpair = 0.0;

        const double r2inv = 1.0/rsq;
        const int cgt=cg_type[itype][jtype];

        if (rsq < cut_out_off_sq) {
          if (COUL_TYPE == CG_COUL_NONE) {
            CG_LJ_INNER(0,fpair);
            fpair *= r2inv;
            if (rsq > cut_out_on_sq) {
              const double rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
              fpair *= 1.0 - rsw*rsw*(3.0 - 2.0*rsw);
            }
          } else {
            double forcelj  = 0.0;
            double forcecoul = 0.0;

            if (rsq < cut_ljsq[itype][jtype]) {
              CG_LJ_INNER(0,forcelj);
            }

            forcecoul = qqrd2e * qtmp*q[j]*sqrt(r2inv);
            if (factor_coul < 1.0) forcecoul -= (1.0 -factor_coul)*forcecoul;

            fpair = (forcecoul + forcelj) * r2inv;
            if (rsq > cut_out_on_sq) {
              const double rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
              fpair *= 1.0 - rsw*rsw*(3.0 - 2.0*rsw);
            }
          }

          f[i][0] += delx*fpair;
          f[i][1] += dely*fpair;
          f[i][2] += delz*fpair;
          if (NEWTON_PAIR || j < nlocal) {
            f[j][0] -= delx*fpair;
            f[j][1] -= dely*fpair;
            f[j][2] -= delz*fpair;
          }
        }
      }
    }
  }

/* ---------------------------------------------------------------------- */

  template < const int NEWTON_PAIR, const int COUL_TYPE >
    void PairCMMCommon::eval_middle()
  {
    double ** const x = atom->x;
    double ** const f = atom->f;
    const double  * const q = atom->q;
    const int * const type = atom->type;
    const int nlocal = atom->nlocal;
    const double * const special_lj = force->special_lj;
    const double * const special_coul = force->special_coul;
    const double qqrd2e = force->qqrd2e;
    double factor_lj,factor_coul;

    const int inum = listmiddle->inum;
    const int * const ilist = listmiddle->ilist;
    const int * const numneigh = listmiddle->numneigh;
    int * const * const firstneigh = listmiddle->firstneigh;

    const double cut_in_off = cut_respa[0];
    const double cut_in_on = cut_respa[1];
    const double cut_out_on = cut_respa[2];
    const double cut_out_off = cut_respa[3];

    const double cut_in_diff = cut_in_on - cut_in_off;
    const double cut_out_diff = cut_out_off - cut_out_on;
    const double cut_in_off_sq = cut_in_off*cut_in_off;
    const double cut_in_on_sq = cut_in_on*cut_in_on;
    const double cut_out_on_sq = cut_out_on*cut_out_on;
    const double cut_out_off_sq = cut_out_off*cut_out_off;

    // loop over neighbors of my atoms

    int ii,jj;
    for (ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      double qtmp = (COUL_TYPE != CG_COUL_NONE) ? q[i] : 0.0;
      const int itype = type[i];
      const int * const jlist = firstneigh[i];
      const int jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        int j2 = jlist[jj];
        factor_lj = special_lj[sbmask(j2)];
        factor_coul = special_coul[sbmask(j2)];
        const int j = j2 & NEIGHMASK;

        const double delx = xtmp - x[j][0];
        const double dely = ytmp - x[j][1];
        const double delz = ztmp - x[j][2];
        const double rsq = delx*delx + dely*dely + delz*delz;
        const int jtype = type[j];

        double evdwl = 0.0;
        double ecoul = 0.0;
        double fpair = 0.0;

        const double r2inv = 1.0/rsq;
        const int cgt=cg_type[itype][jtype];

        if (rsq < cut_out_off_sq && rsq > cut_in_off_sq) {
          if (COUL_TYPE == CG_COUL_NONE) {
            CG_LJ_INNER(0,fpair);
            fpair *= r2inv;
            if (rsq < cut_in_on_sq) {
              const double rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff;
              fpair *= rsw*rsw*(3.0 - 2.0*rsw);
            }
            if (rsq > cut_out_on_sq) {
              const double rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
              fpair *= 1.0 + rsw*rsw*(2.0*rsw - 3.0);
            }
          } else {
            double forcelj  = 0.0;
            double forcecoul = 0.0;

            if (rsq < cut_ljsq[itype][jtype]) {
              CG_LJ_INNER(0,forcelj);
            }

            forcecoul = qqrd2e * qtmp*q[j]*sqrt(r2inv);
            if (factor_coul < 1.0) forcecoul -= (1.0 -factor_coul)*forcecoul;

            fpair = (forcecoul + forcelj) * r2inv;
            if (rsq < cut_in_on_sq) {
              const double rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff;
              fpair *= rsw*rsw*(3.0 - 2.0*rsw);
            }
            if (rsq > cut_out_on_sq) {
              const double rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
              fpair *= 1.0 + rsw*rsw*(2.0*rsw - 3.0);
            }
          }

          f[i][0] += delx*fpair;
          f[i][1] += dely*fpair;
          f[i][2] += delz*fpair;
          if (NEWTON_PAIR || j < nlocal) {
            f[j][0] -= delx*fpair;
            f[j][1] -= dely*fpair;
            f[j][2] -= delz*fpair;
          }
        }
      }
    }
  }

/* ---------------------------------------------------------------------- */

  template < const int EVFLAG, const int EFLAG, const int VFLAG,
    const int NEWTON_PAIR, const int COUL_TYPE >
    void PairCMMCommon::eval_outer()
  {
    double ** const x = atom->x;
    double ** const f = atom->f;
    const double  * const q = atom->q;
    const int * const type = atom->type;
    const int nlocal = atom->nlocal;
    const double * const special_lj = force->special_lj;
    const double * const special_coul = force->special_coul;
    const double qqrd2e = force->qqrd2e;
    double factor_lj,factor_coul;

    const int inum = listouter->inum;
    const int * const ilist = listouter->ilist;
    const int * const numneigh = listouter->numneigh;
    int * const * const firstneigh = listouter->firstneigh;

    const double cut_in_off = cut_respa[2];
    const double cut_in_on = cut_respa[3];

    const double cut_in_diff = cut_in_on - cut_in_off;
    const double cut_in_off_sq = cut_in_off*cut_in_off;
    const double cut_in_on_sq = cut_in_on*cut_in_on;

    // loop over neighbors of my atoms

    int ii,jj;
    for (ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      double qtmp = (COUL_TYPE != CG_COUL_NONE) ? q[i] : 0.0;
      const int itype = type[i];
      const int * const jlist = firstneigh[i];
      const int jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        int j2 = jlist[jj];
        factor_lj = special_lj[sbmask(j2)];
        factor_coul = special_coul[sbmask(j2)];
        const int j = j2 & NEIGHMASK;

        const double delx = xtmp - x[j][0];
        const double dely = ytmp - x[j][1];
        const double delz = ztmp - x[j][2];
        const double rsq = delx*delx + dely*dely + delz*delz;
        const int jtype = type[j];

        double evdwl = 0.0;
        double ecoul = 0.0;
        double fpair = 0.0;

        const double r2inv = 1.0/rsq;
        const int cgt=cg_type[itype][jtype];

        if (rsq < cutsq[itype][jtype]) {
          if (COUL_TYPE == CG_COUL_NONE) {
            double forcelj=0.0;

            if (rsq > cut_in_off_sq) {
              CG_LJ_INNER(0,forcelj);
              fpair = forcelj*r2inv;
              if (rsq < cut_in_on_sq) {
                const double rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff;
                fpair *= rsw*rsw*(3.0 - 2.0*rsw);
              }

              f[i][0] += delx*fpair;
              f[i][1] += dely*fpair;
              f[i][2] += delz*fpair;
              if (NEWTON_PAIR || j < nlocal) {
                f[j][0] -= delx*fpair;
                f[j][1] -= dely*fpair;
                f[j][2] -= delz*fpair;
              }
            }

            CG_LJ_ENERGY(EFLAG);

            if (VFLAG) {
              if (rsq <= cut_in_off_sq) {
                CG_LJ_INNER(0,fpair);
                fpair *= r2inv;
              } else if (rsq < cut_in_on_sq) {
                fpair = forcelj*r2inv;
              }
            }

            if (EVFLAG) ev_tally(i,j,nlocal,NEWTON_PAIR,
                                 evdwl,ecoul,fpair,delx,dely,delz);
          } else {
            double forcelj  = 0.0;
            double forcecoul = 0.0;

            if (rsq < cut_ljsq[itype][jtype]) {
              CG_LJ_INNER(EFLAG,forcelj);
            }

            // coulomb with cutoff and screening
            if ((COUL_TYPE == CG_COUL_CUT) || (COUL_TYPE == CG_COUL_DEBYE)) {
              if (rsq < cut_coulsq[itype][jtype]) {
                double r=sqrt(rsq);
                double qscreen=exp(-kappa*r);
                forcecoul = factor_coul * qqrd2e
                  * qtmp * q[j] * qscreen * (kappa + 1.0/r);
                if (EFLAG) ecoul=factor_coul*qqrd2e
                  * qtmp*q[j] * qscreen / r;
              }
            }

            if (COUL_TYPE == CG_COUL_LONG) {
              if (rsq < cut_coulsq_global) {
                if (!ncoultablebits || rsq <= tabinnersq) {
                  const double r = sqrt(rsq);

                  const double grij = g_ewald * r;
                  const double expm2 = exp(-grij*grij);
                  const double t = 1.0 / (1.0 + EWALD_P*grij);
                  const double erfc = t * (EWALD_A1+t*(EWALD_A2+t*(EWALD_A3+t*(EWALD_A4+t*EWALD_A5)))) * expm2;
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
                  const double fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
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
              }
            }
            fpair = (forcecoul + forcelj) * r2inv;
            f[i][0] += delx*fpair;
            f[i][1] += dely*fpair;
            f[i][2] += delz*fpair;
            if (NEWTON_PAIR || j < nlocal) {
              f[j][0] -= delx*fpair;
              f[j][1] -= dely*fpair;
              f[j][2] -= delz*fpair;
            }
            if (EVFLAG) ev_tally(i,j,nlocal,NEWTON_PAIR,
                                 evdwl,ecoul,fpair,delx,dely,delz);
          }
        }
      }
    }
  }
/* ------------------------------------------------------------------------ */

}

#undef EWALD_F
#undef EWALD_P
#undef EWALD_A1
#undef EWALD_A2
#undef EWALD_A3
#undef EWALD_A4
#undef EWALD_A5
#endif
