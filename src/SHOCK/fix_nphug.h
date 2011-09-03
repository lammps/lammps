/* ----------------------------------------------------------------------
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

FixStyle(nphug,FixNPHug)

#else

#ifndef LMP_FIX_NPHUG_H
#define LMP_FIX_NPHUG_H

#include "fix_nh.h"

namespace LAMMPS_NS {

class FixNPHug : public FixNH {
 public:
  FixNPHug(class LAMMPS *, int, char **);
  ~FixNPHug();
  void init();
  void setup(int);
  int modify_param(int, char **);
  int pack_restart_data(double *); // pack restart data
  int size_restart();              // return size
  void restart(char *);
 
 private:
  class Compute *pe;               // PE compute pointer

  void compute_temp_target();
  double compute_vector(int);
  double compute_etotal();
  double compute_vol();
  double compute_hugoniot();
  double compute_us();
  double compute_up();

  char *id_pe;
  int peflag;
  int v0_set,p0_set,e0_set;
  double v0,p0,e0,rho0;
  int idir;
  int uniaxial;
};

}

#endif
#endif
