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
// clang-format off
FixStyle(update/special/bonds,FixUpdateSpecialBonds)
// clang-format on
#else

#ifndef LMP_FIX_UPDATE_SPECIAL_BONDS_H
#define LMP_FIX_UPDATE_SPECIAL_BONDS_H

#include "fix.h"

#include <set>
#include <utility>

namespace LAMMPS_NS {

class FixUpdateSpecialBonds : public Fix {
 public:
  FixUpdateSpecialBonds(class LAMMPS *, int, char **);
  ~FixUpdateSpecialBonds();
  int setmask();
  void setup(int);
  void pre_exchange();
  void pre_force(int);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);  
  void add_broken_bond(int,int);

 protected:
  std::set <std::pair<tagint, tagint>> broken_pairs;
  inline int sbmask(int j) const {
    return j >> SBBITS & 3;
  }
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix censor/bonded/pairs command

Self-explanatory.

E: Fix censor/bonded/pairs requires Newton bond off

Self-explanatory.

E: Fix censor/bonded/pairs requires atom bonds

Self-explanatory.

E: Fix censor/bonded/pairs must be used without special bonds

Self-explanatory. Look at the atom modify special command.

E: Fix censor/bonded/pairs requires special_bonds = 0,0,0

Self-explanatory.

*/
