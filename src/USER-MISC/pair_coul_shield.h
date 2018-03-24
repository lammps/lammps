/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(coul/shield,PairCoulShield)

#else

#ifndef LMP_PAIR_COUL_SHIELD_H
#define LMP_PAIR_COUL_SHIELD_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCoulShield : public Pair {
 public:
  PairCoulShield(class LAMMPS *);
  virtual ~PairCoulShield();

  virtual void compute(int, int);

  virtual void settings(int, char **);
  virtual void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
  virtual void write_restart(FILE *);
  virtual void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);

  virtual double single(int, int, int, int, double, double, double, double &);

 protected:
  double cut_global;
  double **cut;
  double **sigmae, **offset;
  //double a_eps, b_eps, eps_s;
  int tap_flag;

  void allocate();

/* ----Calculate the long-range cutoff term */
  inline double calc_Tap(double r_ij, double Rcut) {
    double Tap,r;
    //int Tap_coeff[8] = {1,0,0,0,-35,84,-70,20};
    double Tap_coeff[8] = {1.0,0.0,0.0,0.0,-35.0,84.0,-70.0,20.0};

    r = r_ij/Rcut;
    Tap = 0.0;

    Tap = Tap_coeff[7] * r + Tap_coeff[6];
    Tap = Tap * r  + Tap_coeff[5];
    Tap = Tap * r  + Tap_coeff[4];
    Tap = Tap * r  + Tap_coeff[3];
    Tap = Tap * r  + Tap_coeff[2];
    Tap = Tap * r  + Tap_coeff[1];
    Tap = Tap * r  + Tap_coeff[0];

    return(Tap);
  }

 /* ----Calculate the derivatives of long-range cutoff term */
  inline double calc_dTap(double r_ij, double Rcut) {
    double dTap,r;
    double Tap_coeff[8] = {1.0,0.0,0.0,0.0,-35.0,84.0,-70.0,20.0};

    r = r_ij/Rcut;
    dTap = 0.0;

    dTap = 7.0*Tap_coeff[7] * r + 6.0*Tap_coeff[6];
    dTap = dTap * r  + 5.0*Tap_coeff[5];
    dTap = dTap * r  + 4.0*Tap_coeff[4];
    dTap = dTap * r  + 3.0*Tap_coeff[3];
    dTap = dTap * r  + 2.0*Tap_coeff[2];
    dTap = dTap * r  + Tap_coeff[1];
    dTap = dTap/Rcut;

    return(dTap);
  }
};

}

#endif
#endif
