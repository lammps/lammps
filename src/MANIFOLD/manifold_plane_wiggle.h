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

#ifndef LMP_MANIFOLD_PLANE_WIGGLE_H
#define LMP_MANIFOLD_PLANE_WIGGLE_H

#include "manifold.h"

namespace LAMMPS_NS {

namespace user_manifold {

  // A 2D wiggly/wave-y plane (Like z = A cos(kx))
  class manifold_plane_wiggle : public manifold {
   public:
    enum { NPARAMS = 2 };    // Number of parameters.
    manifold_plane_wiggle(LAMMPS *lmp, int, char **);
    double g(const double *x) override;
    void n(const double *x, double *n) override;
    static const char *type() { return "plane/wiggle"; }
    const char *id() override { return type(); }
    static int expected_argc() { return NPARAMS; }
    int nparams() override { return NPARAMS; }
  };
}    // namespace user_manifold

}    // namespace LAMMPS_NS

#endif    // LMP_MANIFOLD_PLANE_WIGGLE_H
