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

#ifndef LMP_MANIFOLD_DUMBBELL_H
#define LMP_MANIFOLD_DUMBBELL_H

#include "manifold.h"

namespace LAMMPS_NS {

namespace user_manifold {

  // A dendritic dumbbell approximation:
  class manifold_dumbbell : public manifold {
   public:
    enum { NPARAMS = 4 };    // Number of parameters.
    manifold_dumbbell(LAMMPS *lmp, int, char **);
    double g(const double *x) override;
    void n(const double *x, double *nn) override;

    static const char *type() { return "dumbbell"; }
    const char *id() override { return type(); }

    static int expected_argc() { return NPARAMS; }
    int nparams() override { return NPARAMS; }
  };
}    // namespace user_manifold

}    // namespace LAMMPS_NS

#endif    // LMP_MANIFOLD_DUMBBELL_H
