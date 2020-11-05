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

#ifdef FIX_CLASS
FixStyle(esfriction, FixESFriction)
#else

#ifndef LMP_FIX_ESFRICTION_H
#define LMP_FIX_ESFRICTION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixESFriction : public Fix {
 public:
  FixESFriction(class LAMMPS *, int, char **);
  virtual ~FixESFriction();
  int setmask();
  void init();
  void post_run();
  void post_force(int);
  double get_ehat();
  double get_kenergy(double, double, double, double);

  double const ev_to_joules = 1.6022e-19;
  double const joules_to_ev = 1.0 / ev_to_joules;
  double const Nav = 6.0221e23; // Avogadro's constant
  double const bohr_radius = 5.2917721067e-1 ; // Angstroms
  double const h_cross = 1.0546e-34; // planck's constant/(2pi) in J-s
  double const m_e = 9.1090e-31; // electron mass in Kg
  double const e_charge = 1.2; // sqrt(eV-nm)
  double const pi = 3.14159265359;

 protected:
  int my_rank; // process rank
  // input quantities
  double Z1, Z2, A1, A2; // atomic number and mass number of PKA & target
  double Ne_val; // no. of electrons in outermost shell; valance electrons
  double ES_cutoff; // Energy cutoff for applying electronic stopping
  double Vol_molar; // volume per mole
  double pka_vx, pka_vy, pka_vz; // pka velocities
  // computed quantities
  double beta; // coefficient of friction due to ES (computed)
  double E_pka; // energy of PKA
  double ES_loss_NRT, ES_loss_LS; // ES loss using NRT & LS models
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
