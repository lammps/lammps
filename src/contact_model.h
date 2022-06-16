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

#ifndef LMP_CONTACT_MODEL_H
#define LMP_CONTACT_MODEL_H

namespace LAMMPS_NS {
namespace Contact_Model {

  class ContactModel {
   public:
    ContactModel();
    void set_strain(const double, const double);
    void step_deform(const double, const double);
    bool reduce();
    void get_box(double[3][3], double);
    void get_rot(double[3][3]);
    void get_inverse_cob(int[3][3]);

   private:
     double a, knfac;
  };
}    // namespace Contact_Model
}    // namespace LAMMPS_NS
#endif
