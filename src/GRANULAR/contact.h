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
    void reset_contact();
    bool check_contact();
    void prep_contact();
    void calculate_forces(double *, double *, double *, double *);
    double pulloff_distance(double, double);

    int normal_model, damping_model, tangential_model;
    int roll_model, twist_model;
    int limit_damping;
    double cutoff_type;
    double Emod, poisson;                    // variables used in defining mixed interactions
    double k_norm, gamma_norm, cohesion;     // normal_coeffs
    double k_tang, gamma_tang, mu_tang;      // tangential_coeffs - wutang?
    double k_roll, gamma_roll, mu_roll;      // roll_coeffs
    double k_twist, gamma_twist, mu_twist;   // twist_coeffs

    double radi, radj, meff, dt;
    double xi[3], xj[3], vi[3], vj[3], omegai[3], omegaj[3];
    int history_update, roll_history_index, twist_history_index;

   private:
    double a, knfac, Fncrit, Fscrit, Frcrit, damp_normal_prefactor;
    double fs[3], fr[3], ft[3];
    double dx[3], nx[3], r, rsq, rinv, Reff, radsum, delta;
    double vr[3], vn[3], vnnr, vt[3], wr[3], vtr[3], vrel;
    double magtwist, magtortwist;
    bool touch;

    int prep_flag, check_flag;
    int mindlin_rescale, mindlin_force;

    void touch_JKR(int);
    double normal_JKR();
    double normal_DMT();
    double normal_Hertz();
    double normal_Hooke();
    double normal_damping();
    void tangential_no_history();
    void tangential_history(double *);
    void tangential_mindlin(double *);
    void rolling(double *);
    void twisting_marshall(double *);
    void twisting_SDS(double *);

  };

}    // namespace Contact
}    // namespace LAMMPS_NS
#endif
