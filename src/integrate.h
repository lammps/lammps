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

#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Integrate : protected Pointers {
 public:
  Integrate(class LAMMPS *, int, char **);
  virtual ~Integrate();
  virtual void init() = 0;
  virtual void setup() = 0;
  virtual void setup_minimal(int) = 0;
  virtual void run(int) = 0;
  virtual void cleanup() {}
  virtual void reset_dt() {}
  virtual double memory_usage() {return 0.0;}

 protected:
  int eflag,vflag;                  // flags for energy/virial computation
  int virial_style;                 // compute virial explicitly or implicitly

  int nelist_global,nelist_atom;    // # of PE,virial computes to check
  int nvlist_global,nvlist_atom;
  class Compute **elist_global;     // lists of PE,virial Computes
  class Compute **elist_atom;
  class Compute **vlist_global;
  class Compute **vlist_atom;

  void ev_setup();
  void ev_set(int);
};

}

#endif
