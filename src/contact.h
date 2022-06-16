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

#ifndef LMP_CONTACT
_H
#define LMP_CONTACT
_H

namespace LAMMPS_NS {
namespace Contact {

  class ContactModel {
   public:
    ContactModel();
    void touch_JKR(int);
    void normal_JKR(double&);
    void normal_DMT(double&);
    void normal_Hooke(double&);
    double normal_damping(double, double);
    double critical_normal(double, double);

    int normal_model, damping_model, tangential_model;
    int roll_model, twist_model;
    double E, G, poisson, damp, coh;

   private:
     double a, knfac;
     ContactGeometry geom;
  };


  class ContactGeometry {
   public:
    ContactGeometry();
    void add_data();
    double r, rinv, rsq, Reff, radsum, delta, dR;
  };
}    // namespace Contact
}    // namespace LAMMPS_NS
#endif
