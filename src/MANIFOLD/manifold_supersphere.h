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

#ifndef LMP_MANIFOLD_SUPERSPHERE_H
#define LMP_MANIFOLD_SUPERSPHERE_H

#include "manifold.h"

namespace LAMMPS_NS {

namespace user_manifold {

  // A sphere:
  class manifold_supersphere : public manifold {
   public:
    enum { NPARAMS = 2 };
    manifold_supersphere(LAMMPS *lmp, int, char **) : manifold(lmp) {}

    virtual ~manifold_supersphere() {}

    double my_sign(double a) { return (a > 0) - (a < 0); }

    virtual double g(const double *x)
    {
      double R = params[0];
      double q = params[1];
      double xx = fabs(x[0]);
      double yy = fabs(x[1]);
      double zz = fabs(x[2]);

      double rr = pow(xx, q) + pow(yy, q) + pow(zz, q);

      return rr - pow(R, q);
    }

    virtual void n(const double *x, double *nn)
    {
      double q = params[1];
      double xx = fabs(x[0]);
      double yy = fabs(x[1]);
      double zz = fabs(x[2]);

      nn[0] = q * my_sign(x[0]) * pow(xx, q - 1);
      nn[1] = q * my_sign(x[1]) * pow(yy, q - 1);
      nn[2] = q * my_sign(x[2]) * pow(zz, q - 1);
    }

    static const char *type() { return "supersphere"; }
    virtual const char *id() { return type(); }
    static int expected_argc() { return NPARAMS; }
    virtual int nparams() { return NPARAMS; }
  };
}    // namespace user_manifold

}    // namespace LAMMPS_NS

#endif
