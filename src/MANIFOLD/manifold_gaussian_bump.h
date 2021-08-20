/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   -----------------------------------------------------------------------

   This file is a part of the MANIFOLD package.

   This package allows LAMMPS to perform MD simulations of particles
   constrained on a manifold (i.e., a 2D subspace of the 3D simulation
   box). It achieves this using the RATTLE constraint algorithm applied
   to single-particle constraint functions g(xi,yi,zi) = 0 and their
   derivative (i.e. the normal of the manifold) n = grad(g).

   It is very easy to add your own manifolds to the current zoo
   (we now have sphere, a dendritic spine approximation, a 2D plane (for
   testing purposes) and a wave-y plane.
   See the README file for more info.

   Stefan Paquay, spaquay@brandeis.edu
   Brandeis University, Waltham, MA, USA.

   This package was mainly developed at
   Applied Physics/Theory of Polymers and Soft Matter,
   Eindhoven University of Technology (TU/e), The Netherlands

   Thanks to Remy Kusters at TU/e for testing.

   This software is distributed under the GNU General Public License.

------------------------------------------------------------------------- */

#ifndef LMP_MANIFOLD_GAUSSIAN_BUMP_H
#define LMP_MANIFOLD_GAUSSIAN_BUMP_H

#include "manifold.h"

namespace LAMMPS_NS {

namespace user_manifold {

  // A Gaussian bump with a smoothed decay to flat 2D.
  class manifold_gaussian_bump : public manifold {
   public:
    enum { NPARAMS = 4 };
    manifold_gaussian_bump(class LAMMPS *, int, char **);
    virtual ~manifold_gaussian_bump();

    virtual double g(const double *);
    virtual void n(const double *, double *);

    // Variant of g that computes n at the same time.
    virtual double g_and_n(const double *x, double *nn);

    static const char *type() { return "gaussian_bump"; }
    virtual const char *id() { return type(); }

    virtual int nparams() { return NPARAMS; }
    virtual void post_param_init();

   private:
    // Some private constants:
    double AA, ll, ll2, f_at_rc, fp_at_rc;
    double rc1, rc2, rc12, rc22, dr, inv_dr;

    // Stuff for the look-up table:
    double lut_x0, lut_x1;
    int lut_Nbins;
    double lut_dx;
    double *lut_z;
    double *lut_zp;

    double gaussian_bump(double) const;
    double gaussian_bump_x2(double) const;
    double gaussian_bump_der(double) const;

    void make_lut();
    double lut_get_z(double rr) const;
    double lut_get_zp(double rr) const;
    void lut_get_z_and_zp(double rr, double &zz, double &zzp) const;

    void test_lut();

    double taper(double);
    double taper_der(double);
  };
}    // namespace user_manifold

}    // namespace LAMMPS_NS

#endif    // LMP_MANIFOLD_GAUSSIAN_BUMP_H
