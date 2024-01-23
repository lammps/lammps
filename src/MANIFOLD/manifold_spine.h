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

#ifndef LMP_MANIFOLD_SPINE_H
#define LMP_MANIFOLD_SPINE_H

#include "manifold.h"

namespace LAMMPS_NS {

namespace user_manifold {

  // A dendritic spine approximation:
  class manifold_spine : public manifold {
   public:
    enum { NPARAMS = 5 };    // Number of parameters.
    manifold_spine(LAMMPS *lmp, int, char **);
    double g(const double *x) override;
    void n(const double *x, double *nn) override;
    double g_and_n(const double *x, double *nn) override;

    static const char *type() { return "spine"; }
    const char *id() override { return type(); }

    static int expected_argc() { return NPARAMS; }
    int nparams() override { return NPARAMS; }

   protected:
    int power;
  };

  class manifold_spine_two : public manifold_spine {
   public:
    manifold_spine_two(LAMMPS *lmp, int, char **);

    static const char *type() { return "spine/two"; }
    const char *id() override { return type(); }
  };
}    // namespace user_manifold

}    // namespace LAMMPS_NS

#endif    // LMP_MANIFOLD_SPINE_H
