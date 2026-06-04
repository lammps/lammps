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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/charmm/coul/long/functor,PairLJCharmmCoulLongFunctor);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CHARMM_COUL_LONG_FUNCTOR_H
#define LMP_PAIR_LJ_CHARMM_COUL_LONG_FUNCTOR_H

#include "evaluator_lj_charmm.h"
#include "functor_coul_long.h"
#include "pair_functor.h"

#include "memory.h"

#include <cmath>
#include <cstring>

namespace LAMMPS_NS {

// The pairwise kernel comes from PairFunctor<EvaluatorLJCharmm, CoulLong>.  In
// addition the CHARMM force field computes the 1-4 (dihedral) special LJ/Coulomb
// interactions in dihedral_style charmm, which reads the 1-4 LJ coefficients and
// the "implicit" flag from the pair style via extract().  Those per-type-pair
// 1-4 arrays (lj14_*) are derived here from the eps14/sigma14 coefficients and
// exposed through extract(); they are not used by the pairwise kernel itself.

class PairLJCharmmCoulLongFunctor :
    public PairFunctor<functor::EvaluatorLJCharmm, functor::CoulLong> {
  using Base = PairFunctor<functor::EvaluatorLJCharmm, functor::CoulLong>;

 public:
  PairLJCharmmCoulLongFunctor(class LAMMPS *);

  ~PairLJCharmmCoulLongFunctor() override
  {
    if (copymode) return;
    if (allocated) {
      memory->destroy(lj14_1);
      memory->destroy(lj14_2);
      memory->destroy(lj14_3);
      memory->destroy(lj14_4);
    }
  }

  void allocate() override
  {
    Base::allocate();
    const int n = atom->ntypes + 1;
    memory->create(lj14_1, n, n, "pair:lj14_1");
    memory->create(lj14_2, n, n, "pair:lj14_2");
    memory->create(lj14_3, n, n, "pair:lj14_3");
    memory->create(lj14_4, n, n, "pair:lj14_4");
  }

  double init_one(int i, int j) override
  {
    const double cut = Base::init_one(i, j);

    // 1-4 LJ coefficients for dihedral_style charmm (eps14/sigma14 were mixed
    // into coeffs[i][j] by the base init_one when the pair was not set)
    const double e14 = coeffs[i][j].eps14;
    const double s14 = coeffs[i][j].sigma14;
    lj14_1[i][j] = lj14_1[j][i] = 48.0 * e14 * pow(s14, 12.0);
    lj14_2[i][j] = lj14_2[j][i] = 24.0 * e14 * pow(s14, 6.0);
    lj14_3[i][j] = lj14_3[j][i] = 4.0 * e14 * pow(s14, 12.0);
    lj14_4[i][j] = lj14_4[j][i] = 4.0 * e14 * pow(s14, 6.0);
    return cut;
  }

  void *extract(const char *str, int &dim) override
  {
    dim = 2;
    if (strcmp(str, "lj14_1") == 0) return (void *) lj14_1;
    if (strcmp(str, "lj14_2") == 0) return (void *) lj14_2;
    if (strcmp(str, "lj14_3") == 0) return (void *) lj14_3;
    if (strcmp(str, "lj14_4") == 0) return (void *) lj14_4;
    dim = 0;
    if (strcmp(str, "implicit") == 0) return (void *) &implicit;
    return Base::extract(str, dim);
  }

 protected:
  double **lj14_1, **lj14_2, **lj14_3, **lj14_4;
  int implicit;
};

}    // namespace LAMMPS_NS

#endif
#endif
