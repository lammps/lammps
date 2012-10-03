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

PairStyle(lj/cut/tip4p/long/opt,PairLJCutTIP4PLongOpt)

#else

#ifndef LMP_PAIR_LJ_CUT_TIP4P_LONG_OPT_H
#define LMP_PAIR_LJ_CUT_TIP4P_LONG_OPT_H

#include "pair_lj_cut_tip4p_long.h"

namespace LAMMPS_NS {

  class PairLJCutTIP4PLongOpt : public PairLJCutTIP4PLong {
 public:
    PairLJCutTIP4PLongOpt(class LAMMPS *);
    virtual ~PairLJCutTIP4PLongOpt() {};

  virtual void compute(int, int);
  virtual double memory_usage();

 protected:
  template < const int, const int, const int, const int >
  void eval();
  void compute_newsite_opt(const double *, const double *,
                           const double *, double *) const;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: TIP4P hydrogen is missing

The TIP4P pairwise computation failed to find the correct H atom
within a water molecule.

E: TIP4P hydrogen has incorrect atom type

The TIP4P pairwise computation found an H atom whose type does not
agree with the specified H type.

*/
