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

#ifndef LMP_MANIFOLD_SPHERE_H
#define LMP_MANIFOLD_SPHERE_H

#include "manifold.h"

namespace LAMMPS_NS {

namespace user_manifold {

  // A sphere:
  class manifold_sphere : public manifold {
   public:
    enum { NPARAMS = 1 };
    manifold_sphere(LAMMPS *lmp, int, char **) : manifold(lmp) {}

    virtual ~manifold_sphere() {}
    virtual double g(const double *x)
    {
      double R = params[0];
      double r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
      return r2 - R * R;
    }

    virtual double g_and_n(const double *x, double *nn)
    {
      double R = params[0];
      double r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
      nn[0] = 2 * x[0];
      nn[1] = 2 * x[1];
      nn[2] = 2 * x[2];

      return r2 - R * R;
    }

    virtual void n(const double *x, double *nn)
    {
      nn[0] = 2 * x[0];
      nn[1] = 2 * x[1];
      nn[2] = 2 * x[2];
    }

    virtual void H(double * /*x*/, double h[3][3])
    {
      h[0][1] = h[0][2] = h[1][0] = h[1][2] = h[2][0] = h[2][1] = 0.0;
      h[0][0] = h[1][1] = h[2][2] = 2.0;
    }

    static const char *type() { return "sphere"; }
    virtual const char *id() { return type(); }
    static int expected_argc() { return NPARAMS; }
    virtual int nparams() { return NPARAMS; }
  };
}    // namespace user_manifold

}    // namespace LAMMPS_NS

#endif    // LMP_MANIFOLD_SPHERE_H
