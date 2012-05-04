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

#ifndef LMP_KSPACE_H
#define LMP_KSPACE_H

#include "pointers.h"

namespace LAMMPS_NS {

class KSpace : protected Pointers {
  friend class ThrOMP;
  friend class FixOMP;
 public:
  double energy;                  // accumulated energy
  double virial[6];               // accumlated virial
  double *eatom,**vatom;          // accumulated per-atom energy/virial

  double g_ewald;
  int nx_pppm,ny_pppm,nz_pppm;
 
  int compute_flag;               // 0 if skip compute()

  KSpace(class LAMMPS *, int, char **);
  virtual ~KSpace();
  void modify_params(int, char **);
  void *extract(const char *);
  void compute_dummy(int, int);

  // general child-class methods

  virtual void init() = 0;
  virtual void setup() = 0;
  virtual void compute(int, int) = 0;
  virtual void timing(int, double &, double &) {}
  virtual double memory_usage() {return 0.0;}

 protected:
  int gridflag,gewaldflag,differentiation_flag;
  int order;
  int slabflag;
  int suffix_flag;                  // suffix compatibility flag
  double scale;
  double slab_volfactor;
 
  double accuracy;                  // accuracy of KSpace solver (force units)
  double accuracy_absolute;         // user-specifed accuracy in force units
  double accuracy_relative;         // user-specified dimensionless accuracy
                                    // accurary = acc_rel * two_charge_force
  double two_charge_force;          // force in user units of two point
                                    // charges separated by 1 Angstrom

  int evflag,evflag_atom;
  int eflag_either,eflag_global,eflag_atom;
  int vflag_either,vflag_global,vflag_atom;
  int maxeatom,maxvatom;

  void ev_setup(int, int);
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Bad kspace_modify slab parameter

Kspace_modify value for the slab/volume keyword must be >= 2.0.

W: Kspace_modify slab param < 2.0 may cause unphysical behavior

The kspace_modify slab parameter should be larger to insure periodic
grids padded with empty space do not overlap.

*/
