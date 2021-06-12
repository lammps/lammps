/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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
    virtual ~manifold_plane_wiggle() {}
    virtual double g(const double *x);
    virtual void n(const double *x, double *n);
    static const char *type() { return "plane/wiggle"; }
    virtual const char *id() { return type(); }
    static int expected_argc() { return NPARAMS; }
    virtual int nparams() { return NPARAMS; }
  };
}    // namespace user_manifold

}    // namespace LAMMPS_NS

#endif    // LMP_MANIFOLD_PLANE_WIGGLE_H
