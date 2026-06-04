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

#ifndef LMP_PAIR_LJ_CHARMM_COUL_LONG_FUNCTOR_IMPL_H
#define LMP_PAIR_LJ_CHARMM_COUL_LONG_FUNCTOR_IMPL_H

#include "atom.h"
#include "memory.h"

#include <cmath>
#include <cstring>

namespace LAMMPS_NS {

// Shared CHARMM 1-4 glue for the lj/charmm/coul/long/functor style, written once
// and parameterized on the functor driver base (BASE is the serial
// PairFunctor<EvaluatorLJCharmm, CoulLong> or the threaded
// PairFunctorOMP<EvaluatorLJCharmm, CoulLong>).  The pairwise kernel comes
// entirely from BASE; this mixin only adds the per-type-pair 1-4 (dihedral) LJ
// arrays (lj14_*) and the "implicit" flag that dihedral_style charmm reads back
// from the pair style through extract().  They are derived here from the
// eps14/sigma14 coefficients and are not used by the pairwise kernel itself.
// CHARMM also defaults to arithmetic (Lorentz-Berthelot) mixing.

template <class BASE> class PairLJCharmmCoulLongFunctorImpl : public BASE {
 public:
  PairLJCharmmCoulLongFunctorImpl(class LAMMPS *lmp) : BASE(lmp)
  {
    implicit = 0;
    lj14_1 = lj14_2 = lj14_3 = lj14_4 = nullptr;
    this->mix_flag = BASE::ARITHMETIC;    // CHARMM uses arithmetic mixing by default
  }

  ~PairLJCharmmCoulLongFunctorImpl() override
  {
    if (this->copymode) return;
    if (this->allocated) {
      this->memory->destroy(lj14_1);
      this->memory->destroy(lj14_2);
      this->memory->destroy(lj14_3);
      this->memory->destroy(lj14_4);
    }
  }

  void allocate() override
  {
    BASE::allocate();
    const int n = this->atom->ntypes + 1;
    this->memory->create(lj14_1, n, n, "pair:lj14_1");
    this->memory->create(lj14_2, n, n, "pair:lj14_2");
    this->memory->create(lj14_3, n, n, "pair:lj14_3");
    this->memory->create(lj14_4, n, n, "pair:lj14_4");
  }

  double init_one(int i, int j) override
  {
    const double cut = BASE::init_one(i, j);

    // 1-4 LJ coefficients for dihedral_style charmm (eps14/sigma14 were mixed
    // into coeffs[i][j] by the base init_one when the pair was not set)
    const double e14 = this->coeffs[i][j].eps14;
    const double s14 = this->coeffs[i][j].sigma14;
    lj14_1[i][j] = lj14_1[j][i] = 48.0 * e14 * std::pow(s14, 12.0);
    lj14_2[i][j] = lj14_2[j][i] = 24.0 * e14 * std::pow(s14, 6.0);
    lj14_3[i][j] = lj14_3[j][i] = 4.0 * e14 * std::pow(s14, 12.0);
    lj14_4[i][j] = lj14_4[j][i] = 4.0 * e14 * std::pow(s14, 6.0);
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
    return BASE::extract(str, dim);
  }

 protected:
  double **lj14_1, **lj14_2, **lj14_3, **lj14_4;
  int implicit;
};

}    // namespace LAMMPS_NS

#endif
