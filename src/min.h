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

#ifndef MIN_H
#define MIN_H

#include "pointers.h"

namespace LAMMPS_NS {

class Min : protected Pointers {
 public:
  double einitial,efinal,eprevious;
  double gnorm2_init,gnorminf_init,gnorm2_final,gnorminf_final;
  int niter,neval;
  double dmin,dmax;
  int linestyle,lineiter;

  Min(class LAMMPS *);
  virtual ~Min();
  virtual void init() = 0;
  virtual void run() = 0;
  virtual double memory_usage() {return 0.0;}

  void modify_params(int, char **);

 protected:
  int eflag,vflag;            // flags for energy/virial computation
  int virial_style;           // compute virial explicitly or implicitly

  int nelist_atom;                  // # of PE,virial computes to check
  int nvlist_global,nvlist_atom;
  class Compute **elist_atom;       // list of PE,virial Computes
  class Compute **vlist_global;
  class Compute **vlist_atom;

  void ev_setup();
  void ev_set(int);
};

}

#endif
