/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing authors: Steven E. Strong and Joel D. Eaves
   Joel.Eaves@Colorado.edu
   ------------------------------------------------------------------------- */
#ifdef FIX_CLASS
// clang-format off
FixStyle(flow/gauss,FixFlowGauss);
// clang-format on
#else

#ifndef LMP_FIX_FLOWGAUSS_H
#define LMP_FIX_FLOWGAUSS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixFlowGauss : public Fix {
 public:
  FixFlowGauss(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  double compute_scalar();
  double compute_vector(int n);

 protected:
  int dimension;
  bool flow[3];       //flag if each direction is conserved
  double a_app[3];    //applied acceleration
  double mTot;        //total mass of constrained group
  double f_tot[3];    //total applied force
  double pe_tot;      //total added energy
  double dt;          //timestep
  bool workflag;      //if calculate work done by fix
  int ilevel_respa;
};

}    // namespace LAMMPS_NS

#endif
#endif
