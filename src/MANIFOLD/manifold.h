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

   Stefan Paquay, stefanpaquay@gmail.com
   Applied Physics/Theory of Polymers and Soft Matter,
   Eindhoven University of Technology (TU/e), The Netherlands

   Thanks to Remy Kusters at TU/e for testing.

   This software is distributed under the GNU General Public License.

------------------------------------------------------------------------- */

#ifndef LMP_MANIFOLD_H
#define LMP_MANIFOLD_H

#include "pointers.h"
#include <cmath>

namespace LAMMPS_NS {
namespace user_manifold {

  // Abstract base class.
  class manifold : protected Pointers {
   public:
    manifold(class LAMMPS *lmp) : Pointers(lmp), params(nullptr) {}
    ~manifold() override { delete[] params; }
    virtual double g(const double *) = 0;
    virtual void n(const double *, double *) = 0;

    // Variant of g that computes n at the same time.
    virtual double g_and_n(const double *x, double *nn)
    {
      this->n(x, nn);
      return g(x);
    }

    virtual const char *id() = 0;

    virtual void set_atom_id(tagint /*a_id*/) {}
    virtual int nparams() = 0;
    // double *get_params() { return params; };

    // Overload if any initialization depends on params:
    virtual void post_param_init() {}
    virtual void checkup() {}    // Some diagnostics...

    double *params;
  };

  // Some utility functions that are templated, so I implement them
  // here in the header.
  template <unsigned int size> inline double infnorm(double *vect)
  {
    double largest = fabs(vect[0]);
    for (unsigned int i = 1; i < size; ++i) {
      double c = fabs(vect[i]);
      largest = (c > largest) ? c : largest;
    }
    return largest;
  }

  inline double dot(double *a, double *b)
  {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }

}    // namespace user_manifold

}    // namespace LAMMPS_NS

#endif    // LMP_MANIFOLD_H
