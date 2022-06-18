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

#ifndef LMP_CONTACT_H
#define LMP_CONTACT_H

namespace Contact {

  enum {HOOKE, HERTZ, HERTZ_MATERIAL, DMT, JKR};
  enum {VELOCITY, MASS_VELOCITY, VISCOELASTIC, TSUJI};
  enum {TANGENTIAL_NOHISTORY, TANGENTIAL_HISTORY,
        TANGENTIAL_MINDLIN, TANGENTIAL_MINDLIN_RESCALE,
        TANGENTIAL_MINDLIN_FORCE, TANGENTIAL_MINDLIN_RESCALE_FORCE};
  enum {TWIST_NONE, TWIST_SDS, TWIST_MARSHALL};
  enum {ROLL_NONE, ROLL_SDS};

  #define PI27SQ 266.47931882941264802866    // 27*PI**2
  #define THREEROOT3 5.19615242270663202362  // 3*sqrt(3)
  #define SIXROOT6 14.69693845669906728801   // 6*sqrt(6)
  #define INVROOT6 0.40824829046386307274    // 1/sqrt(6)
  #define FOURTHIRDS (4.0/3.0)               // 4/3
  #define ONETHIRD (1.0/3.0)                 // 1/3
  #define THREEQUARTERS 0.75                 // 3/4

  #define EPSILON 1e-10

  class ContactModel {
   public:
    ContactModel();
    void reset_contact();
    bool check_contact();
    void prep_contact();
    void calculate_forces(double *, double *, double *, double *);
    double calculate_heat();
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
    double conductivity;

    double radi, radj, meff, dt, Ti, Tj;
    double *xi, *xj, *vi, *vj, *omegai, *omegaj;
    int history_update, roll_history_index, twist_history_index;

    double fs[3], fr[3], ft[3], magtortwist;

   private:
    double a, knfac, Fntot, Fncrit, Fscrit, Frcrit, damp_normal_prefactor;
    double dx[3], nx[3], r, rsq, rinv, Reff, radsum, delta, dR;
    double vr[3], vn[3], vnnr, vt[3], wr[3], vtr[3], vrl[3], relrot[3], vrel;
    double magtwist;
    bool touch;

    int prep_flag, check_flag;
    int mindlin_rescale, mindlin_force;

    bool touch_JKR(int);
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

#endif
