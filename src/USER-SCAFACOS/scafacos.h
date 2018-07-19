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

#ifdef KSPACE_CLASS

KSpaceStyle(scafacos,Scafacos)

#else

#ifndef LMP_SCAFACOS_H
#define LMP_SCAFACOS_H

#include "kspace.h"
#include "fcs.h"

namespace LAMMPS_NS {

class Scafacos : public KSpace {
 public:
  Scafacos(class LAMMPS *, int, char **);
  ~Scafacos();
  void init();
  void setup();
  void compute(int, int);
  int modify_param(int, char **);
  double memory_usage();

 private: 
  int me;

  char *method;
  double tolerance;
  double *epot,**efield;
  int tolerance_type;
  int initialized,maxatom;

  FCS fcs;                // ScaFaCoS handle
  FCSResult result;       // result for each ScaFaCoS call

  // simulation state: box, natoms
  // so ScaFaCoS can detect if changes, e.g. for NPT

  fcs_float old_box_x[3],old_box_y[3],old_box_z[3];
  fcs_float old_origin[3];
  fcs_int old_periodicity[3];
  fcs_int old_natoms;

  void check_result(FCSResult);
  void setup_handle();
  bool box_has_changed();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
