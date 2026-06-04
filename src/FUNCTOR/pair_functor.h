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

#ifndef LMP_PAIR_FUNCTOR_H
#define LMP_PAIR_FUNCTOR_H

// Templated driver for simple, pairwise-additive pair styles.
//
//   template <class EVAL, class COUL> class PairFunctor : public Pair
//
// EVAL is a per-potential "evaluator" functor (see functor_evaluators.h) and
// COUL is a compile-time Coulomb policy (see functor_coul_policies.h).  The
// driver owns all of the neighbor-loop boilerplate and bakes in the three
// optimizations used by the OPT/OPENMP packages:
//
//   1. atom->x / atom->f are accessed as a contiguous array of vec3_t structs;
//   2. the force on the outer-loop atom is accumulated in local scalars and
//      written back to the global array exactly once per atom;
//   3. the inner loop is compiled separately for each combination of the
//      energy/virial, energy, Newton, and Coulomb-table flags via a template,
//      so those branches are eliminated from the hot path.
//
// In addition the per-type-pair parameters are stored as a contiguous matrix of
// EVAL::Param structs (params[itype][jtype]) and read through a const reference,
// rather than as a fan of separate double** arrays.  Unlike the OPT package,
// this matrix is the native store (filled once in init_one), so there is no
// per-step rebuild.
//
// A concrete pair style is a thin subclass binding EVAL and COUL together, e.g.
//   class PairLJCutFunctor : public PairFunctor<EvaluatorLJCut, CoulNone> {...};

#include "pair.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "functor_vec.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstddef>
#include <cstdio>

namespace LAMMPS_NS {

template <class EVAL, class COUL> class PairFunctor : public Pair {
 protected:
  using Coeff = typename EVAL::Coeff;
  using Param = typename EVAL::Param;

  double cut_global;
  Coeff **coeffs;    // raw per-type-pair coefficients (for mixing + restart)
  Param *params;     // derived parameters as a flat, 64-byte-aligned matrix of structs;
  int nparams;       // entry (itype,jtype) is params[itype*nparams + jtype]
  COUL coul;         // Coulomb policy state (empty for CoulNone)

 public:
  PairFunctor(LAMMPS *lmp) :
      Pair(lmp), cut_global(0.0), coeffs(nullptr), params(nullptr), nparams(0)
  {
    single_enable = 1;
    respa_enable = 0;
    writedata = 0;
    born_matrix_enable = 0;
  }

  ~PairFunctor() override
  {
    if (copymode) return;
    if (allocated) {
      memory->destroy(setflag);
      memory->destroy(cutsq);
      memory->destroy(coeffs);
      delete[] params;
    }
  }

  // -----------------------------------------------------------------
  // compute: dispatch to the templated kernel (OPT-style branch table)
  // -----------------------------------------------------------------

  void compute(int eflag, int vflag) override
  {
    ev_init(eflag, vflag);
    const int newton = force->newton_pair;

    if constexpr (COUL::has_table) {
      const int ctable = (ncoultablebits != 0);
      if (evflag) {
        if (eflag) {
          if (newton) {
            if (ctable) eval<1, 1, 1, 1>(); else eval<1, 1, 1, 0>();
          } else {
            if (ctable) eval<1, 1, 0, 1>(); else eval<1, 1, 0, 0>();
          }
        } else {
          if (newton) {
            if (ctable) eval<1, 0, 1, 1>(); else eval<1, 0, 1, 0>();
          } else {
            if (ctable) eval<1, 0, 0, 1>(); else eval<1, 0, 0, 0>();
          }
        }
      } else {
        if (newton) {
          if (ctable) eval<0, 0, 1, 1>(); else eval<0, 0, 1, 0>();
        } else {
          if (ctable) eval<0, 0, 0, 1>(); else eval<0, 0, 0, 0>();
        }
      }
    } else {
      if (evflag) {
        if (eflag) {
          if (newton) eval<1, 1, 1, 0>(); else eval<1, 1, 0, 0>();
        } else {
          if (newton) eval<1, 0, 1, 0>(); else eval<1, 0, 0, 0>();
        }
      } else {
        if (newton) eval<0, 0, 1, 0>(); else eval<0, 0, 0, 0>();
      }
    }

    if (vflag_fdotr) virial_fdotr_compute();
  }

  template <int EVFLAG, int EFLAG, int NEWTON_PAIR, int CTABLE> void eval()
  {
    double **x = atom->x;
    double **f = atom->f;
    const int *_noalias type = atom->type;
    const int nlocal = atom->nlocal;
    const double *_noalias special_lj = force->special_lj;

    auto *_noalias xx = (functor::vec3_t *) x[0];
    auto *_noalias ff = (functor::vec3_t *) f[0];

    [[maybe_unused]] const double *_noalias q = nullptr;
    [[maybe_unused]] const double *_noalias special_coul = nullptr;
    [[maybe_unused]] double cut_coulsq = 0.0;
    if constexpr (COUL::needs_charge) {
      q = atom->q;
      special_coul = force->special_coul;
    }
    if constexpr (COUL::has_coul) cut_coulsq = coul.cut_coulsq;

    const int inum = list->inum;
    const int *_noalias ilist = list->ilist;
    const int *_noalias numneigh = list->numneigh;
    const int *const *_noalias firstneigh = list->firstneigh;

    for (int ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      const double xtmp = xx[i].x;
      const double ytmp = xx[i].y;
      const double ztmp = xx[i].z;
      const int itype = type[i];
      const Param *_noalias pi = params + (std::size_t) itype * nparams;
      [[maybe_unused]] double qi = 0.0;
      if constexpr (COUL::needs_charge) qi = q[i];

      const int *_noalias jlist = firstneigh[i];
      const int jnum = numneigh[i];

      double tmpfx = 0.0, tmpfy = 0.0, tmpfz = 0.0;

      for (int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        const int sb = sbmask(j);
        j &= NEIGHMASK;

        const double delx = xtmp - xx[j].x;
        const double dely = ytmp - xx[j].y;
        const double delz = ztmp - xx[j].z;
        const double rsq = delx * delx + dely * dely + delz * delz;
        const int jtype = type[j];
        const Param &p = pi[jtype];

        // p.cutsq is the van der Waals (evaluator) cutoff.  For CoulNone it is
        // the only cutoff; for a Coulomb policy the outer guard is the larger of
        // the vdW and Coulomb cutoffs, and each term is gated separately.
        // NB: a separate sb==0 "special bond" fast path (as in OPT) was tried
        // here and measured *slower* for this driver; do not re-add it without
        // re-benchmarking.

        double fpair = 0.0, evdwl = 0.0, ecoul = 0.0;
        bool within;

        if constexpr (!COUL::has_coul) {
          within = (rsq < p.cutsq);
          if (within) {
            const auto v = EVAL::template pair<EFLAG>(rsq, p, special_lj[sb]);
            fpair = v.fpair;
            evdwl = v.energy;
          }
        } else {
          const double outer = (p.cutsq > cut_coulsq) ? p.cutsq : cut_coulsq;
          within = (rsq < outer);
          if (within) {
            fpair = 0.0;
            evdwl = 0.0;
            if (rsq < p.cutsq) {
              const auto v = EVAL::template pair<EFLAG>(rsq, p, special_lj[sb]);
              fpair = v.fpair;
              evdwl = v.energy;
            }
            if (rsq < cut_coulsq) {
              const auto c =
                  coul.template eval_coul<EFLAG, CTABLE>(rsq, qi, q[j], special_coul[sb]);
              fpair += c.fpair;
              ecoul = c.energy;
            }
          }
        }

        if (within) {
          tmpfx += delx * fpair;
          tmpfy += dely * fpair;
          tmpfz += delz * fpair;
          if (NEWTON_PAIR || j < nlocal) {
            ff[j].x -= delx * fpair;
            ff[j].y -= dely * fpair;
            ff[j].z -= delz * fpair;
          }

          if (EVFLAG)
            ev_tally(i, j, nlocal, NEWTON_PAIR, evdwl, ecoul, fpair, delx, dely, delz);
        }
      }

      ff[i].x += tmpfx;
      ff[i].y += tmpfy;
      ff[i].z += tmpfz;
    }
  }

  // -----------------------------------------------------------------
  // setup / coefficients
  // -----------------------------------------------------------------

  virtual void allocate()
  {
    allocated = 1;
    const int n = atom->ntypes + 1;

    memory->create(setflag, n, n, "pair:setflag");
    for (int i = 1; i < n; i++)
      for (int j = i; j < n; j++) setflag[i][j] = 0;

    memory->create(cutsq, n, n, "pair:cutsq");
    memory->create(coeffs, n, n, "pair:coeffs");

    // the hot parameter matrix is a flat, contiguous block of cache-line-sized
    // Param structs (alignas(64)) so that params[itype*nparams+jtype] reads
    // exactly one cache line.  Param's over-alignment makes plain new[] use the
    // aligned operator new[], which is portable (unlike posix_memalign) and
    // guarantees the 64-byte base alignment in all builds.
    nparams = n;
    params = new Param[(std::size_t) n * n];
  }

  void settings(int narg, char **arg) override
  {
    if constexpr (!COUL::has_coul) {
      if (narg != 1) error->all(FLERR, "Illegal pair_style command");
      cut_global = utils::numeric(FLERR, arg[0], false, lmp);
    } else {
      // combined vdW + Coulomb: "cut_lj [cut_coul]"; cut_coul defaults to cut_lj
      if (narg < 1 || narg > 2) error->all(FLERR, "Illegal pair_style command");
      cut_global = utils::numeric(FLERR, arg[0], false, lmp);
      const double cut_coul = (narg == 2) ? utils::numeric(FLERR, arg[1], false, lmp) : cut_global;
      coul.cut_coulsq = cut_coul * cut_coul;
    }

    if (allocated) {
      for (int i = 1; i <= atom->ntypes; i++)
        for (int j = i; j <= atom->ntypes; j++)
          if (setflag[i][j]) coeffs[i][j].cut = cut_global;
    }
  }

  void coeff(int narg, char **arg) override
  {
    if (narg < 2) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
    if (!allocated) allocate();

    int ilo, ihi, jlo, jhi;
    utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
    utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

    Coeff c = EVAL::parse(narg, arg, lmp, cut_global);

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
    if constexpr (COUL::needs_charge) {
      if (!atom->q_flag)
        error->all(FLERR, "Pair style {} requires atom attribute q", EVAL::name());
    }
    if constexpr (COUL::has_coul) coul.init_style(lmp);
    neighbor->add_request(this);
  }

  double init_one(int i, int j) override
  {
    if (setflag[i][j] == 0) {
      if constexpr (EVAL::has_mixing)
        coeffs[i][j] = EVAL::mix(coeffs[i][i], coeffs[j][j], this);
      else
        error->all(FLERR, "All pair coeffs are not set");
    }

    const Param p = EVAL::derive(coeffs[i][j], offset_flag);
    params[(std::size_t) i * nparams + j] = p;
    params[(std::size_t) j * nparams + i] = p;

    // the neighbor/communication cutoff for this pair is the larger of the vdW
    // (evaluator) cutoff and the Coulomb cutoff
    double cut = coeffs[i][j].cut;
    if constexpr (COUL::has_coul) {
      const double cut_coul = sqrt(coul.cut_coulsq);
      if (cut_coul > cut) cut = cut_coul;
    }
    return cut;
  }

  // -----------------------------------------------------------------
  // single() for compute group/group, pressure, pair_write
  // -----------------------------------------------------------------

  double single(int i, int j, int itype, int jtype, double rsq, double factor_coul,
                double factor_lj, double &fforce) override
  {
    const Param &p = params[(std::size_t) itype * nparams + jtype];
    auto [fpair, evdwl] = EVAL::template pair<true>(rsq, p, factor_lj);
    double ecoul = 0.0;

    if constexpr (COUL::has_coul) {
      auto [fcoul, ec] = coul.single_coul(lmp, i, j, rsq, factor_coul);
      fpair += fcoul;
      ecoul = ec;
    } else {
      (void) i;
      (void) j;
      (void) factor_coul;
    }

    fforce = fpair;
    return evdwl + ecoul;
  }

  // -----------------------------------------------------------------
  // restart I/O (serialize the raw Coeff structs)
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
    fwrite(&cut_global, sizeof(double), 1, fp);
    fwrite(&offset_flag, sizeof(int), 1, fp);
    fwrite(&mix_flag, sizeof(int), 1, fp);
    fwrite(&tail_flag, sizeof(int), 1, fp);
    if constexpr (COUL::has_coul) coul.write_restart(fp);
  }

  void read_restart_settings(FILE *fp) override
  {
    const int me = comm->me;
    if (me == 0) {
      utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
      utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, &tail_flag, sizeof(int), 1, fp, nullptr, error);
    }
    MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
    MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
    MPI_Bcast(&tail_flag, 1, MPI_INT, 0, world);
    if constexpr (COUL::has_coul) coul.read_restart(fp, lmp);
  }
};

}    // namespace LAMMPS_NS

#endif
