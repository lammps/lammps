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

#ifndef LMP_FIX_HYPER_H
#define LMP_FIX_HYPER_H

#include "fix.h"

namespace LAMMPS_NS {

class FixHyper : public Fix {
 public:
  FixHyper(class LAMMPS *, int, char **);
  virtual ~FixHyper() {}
  void *extract(const char *, int &);

  // must be provided by child class

  virtual void init_hyper() = 0;
  virtual void build_bond_list(int) = 0;
  virtual double query(int) = 0;

 protected:
  int hyperflag;
};

}

#endif

/* ERROR/WARNING messages:

*/
