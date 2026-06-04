/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_PAIR_FUNCTOR_TIP4P_H
#define LMP_PAIR_FUNCTOR_TIP4P_H

// Dedicated FUNCTOR driver for TIP4P-style pair styles: a van der Waals
// evaluator EVAL combined with long-range Coulomb where the Coulomb charge of
// the O atom sits on a virtual "M" site.  Reimplements
// src/KSPACE/pair_lj_cut_tip4p_long.cpp.  Unlike the policy-based PairFunctor,
// the M-site geometry changes the neighbor-loop structure (charge sites differ
// from atom positions, forces on the M site are redistributed to O and H, and
// the virial cannot be computed as F dot r), so this is a separate base that
// composes EVAL (for the LJ term) and a CoulLong policy member (for the
// real-space Coulomb math + tables) rather than reusing PairFunctor's kernel.
//
// Only core headers are needed, so this builds without the KSPACE package; at
// run time it requires a kspace_style plus bond and angle styles (for the M-site
// geometry).

#include "pair.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "functor_coul_long.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <string>

namespace LAMMPS_NS {

template <class EVAL> class PairFunctorTIP4P : public Pair {
 protected:
  using Coeff = typename EVAL::Coeff;
  using Param = typename EVAL::Param;

  double cut_global;        // global LJ cutoff
  Coeff **coeffs;           // raw per-type-pair vdW coefficients
  Param *params;            // derived per-type-pair vdW parameters (flat, aligned)
  int nparams;              // matrix stride
  functor::CoulLong coul;   // long-range Coulomb math + tables (global cutoff)
  typename EVAL::Global gvars;    // evaluator style-global parameters (empty for LJ)

  int typeO, typeH, typeB, typeA;
  std::string typeO_str, typeH_str, typeB_str, typeA_str;
  double qdist;             // O-M distance
  double alpha;             // M-site geometric factor

  int nmax;
  int **hneigh;             // [i][0,1] = the two H of O atom i; [2] = newsite valid flag
  double **newsite;         // M-site position for each O atom

  void compute_newsite(const double *xO, const double *xH1, const double *xH2, double *xM) const
  {
    const double delx1 = xH1[0] - xO[0], dely1 = xH1[1] - xO[1], delz1 = xH1[2] - xO[2];
    const double delx2 = xH2[0] - xO[0], dely2 = xH2[1] - xO[1], delz2 = xH2[2] - xO[2];
    xM[0] = xO[0] + alpha * 0.5 * (delx1 + delx2);
    xM[1] = xO[1] + alpha * 0.5 * (dely1 + dely2);
    xM[2] = xO[2] + alpha * 0.5 * (delz1 + delz2);
  }

 public:
  PairFunctorTIP4P(LAMMPS *lmp) :
      Pair(lmp), cut_global(0.0), coeffs(nullptr), params(nullptr), nparams(0), typeO(0), typeH(0),
      typeB(0), typeA(0), qdist(0.0), alpha(0.0), nmax(0), hneigh(nullptr), newsite(nullptr)
  {
    tip4pflag = 1;
    ewaldflag = pppmflag = 1;
    single_enable = 0;
    respa_enable = 0;
    writedata = 0;    // FUNCTOR styles do not write pair coeffs to data files
    no_virial_fdotr_compute = 1;    // virial computed explicitly (bonded H far from O)
  }

  ~PairFunctorTIP4P() override
  {
    if (copymode) return;
    if (allocated) {
      memory->destroy(setflag);
      memory->destroy(cutsq);
      memory->destroy(coeffs);
      delete[] params;
    }
    memory->destroy(hneigh);
    memory->destroy(newsite);
  }

  void allocate()
  {
    allocated = 1;
    const int n = atom->ntypes + 1;
    memory->create(setflag, n, n, "pair:setflag");
    for (int i = 1; i < n; i++)
      for (int j = i; j < n; j++) setflag[i][j] = 0;
    memory->create(cutsq, n, n, "pair:cutsq");
    memory->create(coeffs, n, n, "pair:coeffs");
    nparams = n;
    params = new Param[(std::size_t) n * n];
    coul.allocate(n);
  }

  // -----------------------------------------------------------------
  // compute (M-site neighbor loop; reimplements PairLJCutTIP4PLong::compute)
  // -----------------------------------------------------------------

  void compute(int eflag, int vflag) override
  {
    double evdwl, ecoul;
    int vlist[6];
    double v[6], fd[3], fO[3], fH[3];

    evdwl = ecoul = 0.0;
    ev_init(eflag, vflag);

    const int nlocal = atom->nlocal;
    const int nall = nlocal + atom->nghost;

    if (atom->nmax > nmax) {
      nmax = atom->nmax;
      memory->destroy(hneigh);
      memory->create(hneigh, nmax, 3, "pair:hneigh");
      memory->destroy(newsite);
      memory->create(newsite, nmax, 3, "pair:newsite");
    }
    if (neighbor->ago == 0)
      for (int i = 0; i < nall; i++) hneigh[i][0] = -1;
    for (int i = 0; i < nall; i++) hneigh[i][2] = 0;

    double **f = atom->f;
    double **x = atom->x;
    double *q = atom->q;
    tagint *tag = atom->tag;
    int *type = atom->type;
    double *special_coul = force->special_coul;
    double *special_lj = force->special_lj;
    const int newton_pair = force->newton_pair;

    const double qqrd2e = coul.qqrd2e;
    const double g_ewald = coul.g_ewald;
    const int ncoultablebits = coul.ncoultablebits;
    const double tabinnersq = coul.tabinnersq;
    const double cut_coul = coul.cut_coul_global;
    const double cut_coulsq = cut_coul * cut_coul;
    const double cut_coulsqplus = (cut_coul + 2.0 * qdist) * (cut_coul + 2.0 * qdist);

    const int inum = list->inum;
    int *ilist = list->ilist;
    int *numneigh = list->numneigh;
    int **firstneigh = list->firstneigh;

    for (int ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      const double qtmp = q[i];
      const double xtmp = x[i][0], ytmp = x[i][1], ztmp = x[i][2];
      const int itype = type[i];
      const Param *pi = params + (std::size_t) itype * nparams;

      int iH1 = -1, iH2 = -1;
      double *x1;
      if (itype == typeO) {
        if (hneigh[i][0] < 0) {
          iH1 = atom->map(tag[i] + 1);
          iH2 = atom->map(tag[i] + 2);
          if (iH1 == -1 || iH2 == -1) error->one(FLERR, "TIP4P hydrogen is missing");
          if (type[iH1] != typeH || type[iH2] != typeH)
            error->one(FLERR, "TIP4P hydrogen has incorrect atom type");
          iH1 = domain->closest_image(i, iH1);
          iH2 = domain->closest_image(i, iH2);
          compute_newsite(x[i], x[iH1], x[iH2], newsite[i]);
          hneigh[i][0] = iH1;
          hneigh[i][1] = iH2;
          hneigh[i][2] = 1;
        } else {
          iH1 = hneigh[i][0];
          iH2 = hneigh[i][1];
          if (hneigh[i][2] == 0) {
            hneigh[i][2] = 1;
            compute_newsite(x[i], x[iH1], x[iH2], newsite[i]);
          }
        }
        x1 = newsite[i];
      } else
        x1 = x[i];

      int *jlist = firstneigh[i];
      const int jnum = numneigh[i];

      for (int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        const double factor_lj = special_lj[sbmask(j)];
        const double factor_coul = special_coul[sbmask(j)];
        j &= NEIGHMASK;

        double delx = xtmp - x[j][0];
        double dely = ytmp - x[j][1];
        double delz = ztmp - x[j][2];
        double rsq = delx * delx + dely * dely + delz * delz;
        const int jtype = type[j];
        const Param &p = pi[jtype];

        // LJ interaction based on the true atom separation

        if (rsq < p.cutsq) {
          const auto lj = EVAL::template pair<true>(rsq, p, factor_lj);
          const double forcelj = lj.fpair;

          f[i][0] += delx * forcelj;
          f[i][1] += dely * forcelj;
          f[i][2] += delz * forcelj;
          f[j][0] -= delx * forcelj;
          f[j][1] -= dely * forcelj;
          f[j][2] -= delz * forcelj;

          evdwl = eflag ? lj.energy : 0.0;
          if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, forcelj, delx, dely, delz);
        }

        // Coulomb interaction based on the (possibly M-site-shifted) separation

        if (rsq < cut_coulsqplus) {
          int jH1 = -1, jH2 = -1;
          if (itype == typeO || jtype == typeO) {
            double *x2;
            if (jtype == typeO) {
              if (hneigh[j][0] < 0) {
                jH1 = atom->map(tag[j] + 1);
                jH2 = atom->map(tag[j] + 2);
                if (jH1 == -1 || jH2 == -1) error->one(FLERR, "TIP4P hydrogen is missing");
                if (type[jH1] != typeH || type[jH2] != typeH)
                  error->one(FLERR, "TIP4P hydrogen has incorrect atom type");
                jH1 = domain->closest_image(j, jH1);
                jH2 = domain->closest_image(j, jH2);
                compute_newsite(x[j], x[jH1], x[jH2], newsite[j]);
                hneigh[j][0] = jH1;
                hneigh[j][1] = jH2;
                hneigh[j][2] = 1;
              } else {
                jH1 = hneigh[j][0];
                jH2 = hneigh[j][1];
                if (hneigh[j][2] == 0) {
                  hneigh[j][2] = 1;
                  compute_newsite(x[j], x[jH1], x[jH2], newsite[j]);
                }
              }
              x2 = newsite[j];
            } else
              x2 = x[j];

            delx = x1[0] - x2[0];
            dely = x1[1] - x2[1];
            delz = x1[2] - x2[2];
            rsq = delx * delx + dely * dely + delz * delz;
          }

          if (rsq < cut_coulsq) {
            const double r2inv = 1.0 / rsq;
            double forcecoul, prefactor = 0.0, erfc = 0.0;
            double fraction = 0.0;
            int itable = 0;
            using namespace EwaldConst;
            if (!ncoultablebits || rsq <= tabinnersq) {
              const double r = sqrt(rsq);
              const double grij = g_ewald * r;
              const double expm2 = exp(-grij * grij);
              const double t = 1.0 / (1.0 + EWALD_P * grij);
              erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
              prefactor = qqrd2e * qtmp * q[j] / r;
              forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);
              if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;
            } else {
              Pair::union_int_float_t rsq_lookup;
              rsq_lookup.f = rsq;
              itable = rsq_lookup.i & coul.ncoulmask;
              itable >>= coul.ncoulshiftbits;
              fraction = ((double) rsq_lookup.f - coul.rtable[itable]) * coul.drtable[itable];
              double table = coul.ftable[itable] + fraction * coul.dftable[itable];
              forcecoul = qtmp * q[j] * table;
              if (factor_coul < 1.0) {
                table = coul.ctable[itable] + fraction * coul.dctable[itable];
                prefactor = qtmp * q[j] * table;
                forcecoul -= (1.0 - factor_coul) * prefactor;
              }
            }

            const double cforce = forcecoul * r2inv;

            // apply Coulomb force; if i or j is an O atom the force acts on the
            // M site and is partitioned over O and the two H atoms
            // (Feenstra, J Comp Chem 20, 786 (1999)), preserving force + torque

            int n = 0, key = 0;

            if (itype != typeO) {
              f[i][0] += delx * cforce;
              f[i][1] += dely * cforce;
              f[i][2] += delz * cforce;
              if (vflag) {
                v[0] = x[i][0] * delx * cforce;
                v[1] = x[i][1] * dely * cforce;
                v[2] = x[i][2] * delz * cforce;
                v[3] = x[i][0] * dely * cforce;
                v[4] = x[i][0] * delz * cforce;
                v[5] = x[i][1] * delz * cforce;
              }
              vlist[n++] = i;
            } else {
              key++;
              fd[0] = delx * cforce;
              fd[1] = dely * cforce;
              fd[2] = delz * cforce;
              fO[0] = fd[0] * (1 - alpha);
              fO[1] = fd[1] * (1 - alpha);
              fO[2] = fd[2] * (1 - alpha);
              fH[0] = 0.5 * alpha * fd[0];
              fH[1] = 0.5 * alpha * fd[1];
              fH[2] = 0.5 * alpha * fd[2];
              f[i][0] += fO[0];
              f[i][1] += fO[1];
              f[i][2] += fO[2];
              f[iH1][0] += fH[0];
              f[iH1][1] += fH[1];
              f[iH1][2] += fH[2];
              f[iH2][0] += fH[0];
              f[iH2][1] += fH[1];
              f[iH2][2] += fH[2];
              if (vflag) {
                const double *xH1 = x[iH1];
                const double *xH2 = x[iH2];
                v[0] = x[i][0] * fO[0] + xH1[0] * fH[0] + xH2[0] * fH[0];
                v[1] = x[i][1] * fO[1] + xH1[1] * fH[1] + xH2[1] * fH[1];
                v[2] = x[i][2] * fO[2] + xH1[2] * fH[2] + xH2[2] * fH[2];
                v[3] = x[i][0] * fO[1] + xH1[0] * fH[1] + xH2[0] * fH[1];
                v[4] = x[i][0] * fO[2] + xH1[0] * fH[2] + xH2[0] * fH[2];
                v[5] = x[i][1] * fO[2] + xH1[1] * fH[2] + xH2[1] * fH[2];
              }
              vlist[n++] = i;
              vlist[n++] = iH1;
              vlist[n++] = iH2;
            }

            if (jtype != typeO) {
              f[j][0] -= delx * cforce;
              f[j][1] -= dely * cforce;
              f[j][2] -= delz * cforce;
              if (vflag) {
                v[0] -= x[j][0] * delx * cforce;
                v[1] -= x[j][1] * dely * cforce;
                v[2] -= x[j][2] * delz * cforce;
                v[3] -= x[j][0] * dely * cforce;
                v[4] -= x[j][0] * delz * cforce;
                v[5] -= x[j][1] * delz * cforce;
              }
              vlist[n++] = j;
            } else {
              key += 2;
              fd[0] = -delx * cforce;
              fd[1] = -dely * cforce;
              fd[2] = -delz * cforce;
              fO[0] = fd[0] * (1 - alpha);
              fO[1] = fd[1] * (1 - alpha);
              fO[2] = fd[2] * (1 - alpha);
              fH[0] = 0.5 * alpha * fd[0];
              fH[1] = 0.5 * alpha * fd[1];
              fH[2] = 0.5 * alpha * fd[2];
              f[j][0] += fO[0];
              f[j][1] += fO[1];
              f[j][2] += fO[2];
              f[jH1][0] += fH[0];
              f[jH1][1] += fH[1];
              f[jH1][2] += fH[2];
              f[jH2][0] += fH[0];
              f[jH2][1] += fH[1];
              f[jH2][2] += fH[2];
              if (vflag) {
                const double *xH1 = x[jH1];
                const double *xH2 = x[jH2];
                v[0] += x[j][0] * fO[0] + xH1[0] * fH[0] + xH2[0] * fH[0];
                v[1] += x[j][1] * fO[1] + xH1[1] * fH[1] + xH2[1] * fH[1];
                v[2] += x[j][2] * fO[2] + xH1[2] * fH[2] + xH2[2] * fH[2];
                v[3] += x[j][0] * fO[1] + xH1[0] * fH[1] + xH2[0] * fH[1];
                v[4] += x[j][0] * fO[2] + xH1[0] * fH[2] + xH2[0] * fH[2];
                v[5] += x[j][1] * fO[2] + xH1[1] * fH[2] + xH2[1] * fH[2];
              }
              vlist[n++] = j;
              vlist[n++] = jH1;
              vlist[n++] = jH2;
            }

            if (eflag) {
              if (!ncoultablebits || rsq <= tabinnersq) {
                ecoul = prefactor * erfc;
              } else {
                const double table = coul.etable[itable] + fraction * coul.detable[itable];
                ecoul = qtmp * q[j] * table;
              }
              if (factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
            } else
              ecoul = 0.0;

            if (evflag) ev_tally_tip4p(key, vlist, v, ecoul, alpha);
          }
        }
      }
    }
  }

  // -----------------------------------------------------------------
  // setup / coefficients
  // -----------------------------------------------------------------

  void settings(int narg, char **arg) override
  {
    if (narg < 6 || narg > 7) error->all(FLERR, "Illegal pair_style command");
    typeO_str = arg[0];
    typeH_str = arg[1];
    typeB_str = arg[2];
    typeA_str = arg[3];
    qdist = utils::numeric(FLERR, arg[4], false, lmp);
    cut_global = utils::numeric(FLERR, arg[5], false, lmp);
    coul.cut_coul_global = (narg == 7) ? utils::numeric(FLERR, arg[6], false, lmp) : cut_global;

    if (allocated) {
      for (int i = 1; i <= atom->ntypes; i++)
        for (int j = i; j <= atom->ntypes; j++)
          if (setflag[i][j]) coeffs[i][j].cut = cut_global;
    }
  }

  void coeff(int narg, char **arg) override
  {
    // resolve the TIP4P atom/bond/angle types (skipped on restart, where the
    // type strings are empty and the integer types are already set)
    if (typeO_str.size() > 0) {
      typeO = utils::expand_type_int(FLERR, typeO_str, Atom::ATOM, lmp, true);
      typeH = utils::expand_type_int(FLERR, typeH_str, Atom::ATOM, lmp, true);
      typeB = utils::expand_type_int(FLERR, typeB_str, Atom::BOND, lmp, true);
      typeA = utils::expand_type_int(FLERR, typeA_str, Atom::ANGLE, lmp, true);
    }

    if (narg < 2) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
    if (!allocated) allocate();

    int ilo, ihi, jlo, jhi;
    utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
    utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

    int nconsumed = 0;
    Coeff c = EVAL::parse(narg, arg, lmp, cut_global, nconsumed);
    if (narg > 2 + nconsumed)
      error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));

    int count = 0;
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo, i); j <= jhi; j++) {
        coeffs[i][j] = c;
        setflag[i][j] = 1;
        count++;
      }
    }
    if (count == 0)
      error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
  }

  void init_style() override
  {
    if (atom->tag_enable == 0)
      error->all(FLERR, "Pair style {} requires atom IDs", force->pair_style);
    if (!force->newton_pair)
      error->all(FLERR, "Pair style {} requires newton pair on", force->pair_style);
    if (!atom->q_flag)
      error->all(FLERR, "Pair style {} requires atom attribute q", force->pair_style);
    if (force->bond == nullptr) error->all(FLERR, "Must use a bond style with TIP4P potential");
    if (force->angle == nullptr) error->all(FLERR, "Must use an angle style with TIP4P potential");

    neighbor->add_request(this);
    coul.init_style(this, lmp);    // qqrd2e, kspace check + g_ewald, build tables

    const double theta = force->angle->equilibrium_angle(typeA);
    const double blen = force->bond->equilibrium_distance(typeB);
    alpha = qdist / (cos(0.5 * theta) * blen);

    const double mincut = coul.cut_coul_global + qdist + blen + neighbor->skin;
    if (comm->get_comm_cutoff() < mincut) {
      if (comm->me == 0)
        error->warning(FLERR, "Increasing communication cutoff to {:.8} for TIP4P pair style",
                       mincut);
      comm->cutghostuser = mincut;
    }
  }

  double init_one(int i, int j) override
  {
    if (setflag[i][j] == 0) {
      if constexpr (EVAL::has_mixing)
        coeffs[i][j] = EVAL::mix(coeffs[i][i], coeffs[j][j], this);
      else
        error->all(FLERR, "All pair coeffs are not set");
    }

    const std::size_t ij = (std::size_t) i * nparams + j;
    const std::size_t ji = (std::size_t) j * nparams + i;
    params[ij] = EVAL::derive(coeffs[i][j], offset_flag, gvars);
    params[ji] = params[ij];

    // no LJ term for any interaction involving a water H atom
    if (i == typeH || j == typeH) {
      params[ij].cutsq = 0.0;
      params[ji].cutsq = 0.0;
    }

    // neighbor/communication cutoff includes the O-M offset for Coulomb
    return MAX(coeffs[i][j].cut, coul.cut_coul_global + 2.0 * qdist);
  }

  // expose the Coulomb cutoff and TIP4P geometry to the KSpace (pppm/tip4p) style
  void *extract(const char *str, int &dim) override
  {
    dim = 0;
    if (strcmp(str, "qdist") == 0) return (void *) &qdist;
    if (strcmp(str, "typeO") == 0) return (void *) &typeO;
    if (strcmp(str, "typeH") == 0) return (void *) &typeH;
    if (strcmp(str, "typeA") == 0) return (void *) &typeA;
    if (strcmp(str, "typeB") == 0) return (void *) &typeB;
    if (strcmp(str, "cut_coul") == 0) return (void *) &coul.cut_coul_global;
    return nullptr;
  }

  // -----------------------------------------------------------------
  // restart I/O
  // -----------------------------------------------------------------

  void write_restart(FILE *fp) override
  {
    write_restart_settings(fp);
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++) {
        fwrite(&setflag[i][j], sizeof(int), 1, fp);
        if (setflag[i][j]) fwrite(&coeffs[i][j], sizeof(Coeff), 1, fp);
      }
  }

  void read_restart(FILE *fp) override
  {
    read_restart_settings(fp);
    allocate();
    const int me = comm->me;
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++) {
        if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
        MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
        if (setflag[i][j]) {
          if (me == 0) utils::sfread(FLERR, &coeffs[i][j], sizeof(Coeff), 1, fp, nullptr, error);
          MPI_Bcast(&coeffs[i][j], sizeof(Coeff), MPI_BYTE, 0, world);
        }
      }
  }

  void write_restart_settings(FILE *fp) override
  {
    fwrite(&typeO, sizeof(int), 1, fp);
    fwrite(&typeH, sizeof(int), 1, fp);
    fwrite(&typeB, sizeof(int), 1, fp);
    fwrite(&typeA, sizeof(int), 1, fp);
    fwrite(&qdist, sizeof(double), 1, fp);
    fwrite(&cut_global, sizeof(double), 1, fp);
    fwrite(&coul.cut_coul_global, sizeof(double), 1, fp);
    fwrite(&offset_flag, sizeof(int), 1, fp);
    fwrite(&mix_flag, sizeof(int), 1, fp);
    fwrite(&tail_flag, sizeof(int), 1, fp);
    fwrite(&ncoultablebits, sizeof(int), 1, fp);
    fwrite(&tabinner, sizeof(double), 1, fp);
  }

  void read_restart_settings(FILE *fp) override
  {
    const int me = comm->me;
    if (me == 0) {
      utils::sfread(FLERR, &typeO, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, &typeH, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, &typeB, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, &typeA, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, &qdist, sizeof(double), 1, fp, nullptr, error);
      utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
      utils::sfread(FLERR, &coul.cut_coul_global, sizeof(double), 1, fp, nullptr, error);
      utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, &tail_flag, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, &ncoultablebits, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, &tabinner, sizeof(double), 1, fp, nullptr, error);
    }
    MPI_Bcast(&typeO, 1, MPI_INT, 0, world);
    MPI_Bcast(&typeH, 1, MPI_INT, 0, world);
    MPI_Bcast(&typeB, 1, MPI_INT, 0, world);
    MPI_Bcast(&typeA, 1, MPI_INT, 0, world);
    MPI_Bcast(&qdist, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&coul.cut_coul_global, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
    MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
    MPI_Bcast(&tail_flag, 1, MPI_INT, 0, world);
    MPI_Bcast(&ncoultablebits, 1, MPI_INT, 0, world);
    MPI_Bcast(&tabinner, 1, MPI_DOUBLE, 0, world);

    // the type strings only exist on the original run; clear them so coeff()
    // keeps the integer types read above
    typeO_str.clear();
    typeH_str.clear();
    typeB_str.clear();
    typeA_str.clear();
  }
};

}    // namespace LAMMPS_NS

#endif
