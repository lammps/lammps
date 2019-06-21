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

#ifndef LMP_ASA_DATA_H
#define LMP_ASA_DATA_H

#include "asa_user.h"
#include "pointers.h"
#include "pair_cac.h"
#include "atom_vec_cac.h"
using namespace std;

namespace LAMMPS_NS {

class Asa_Data : protected Pointers {
 public:
  Asa_Data(class LAMMPS *, class PairCAC *);
  Asa_Data(class LAMMPS *, class AtomVecCAC *);
 ~Asa_Data();
  
  asacg_parm *cgParm;
  asa_parm *asaParm;
  asa_objective *Objective;

  AtomVecCAC *avec_pointer;
  PairCAC *pair_pointer;

  int class_flag; //0 for Pair CAC calling and 1 for Atom Vec calling, increment as new ones added
  int call_asa_cg(double *x,double *lo,double *hi, ASA_INT n,
    double grad_tol, double (*valgrad) (asa_objective *), double *Work, ASA_INT *iWork);
  double myvalue(asa_objective *asa);
  void mygrad(asa_objective *asa);
  double myvalue_surfmin(asa_objective *asa);
  void mygrad_surfmin(asa_objective *asa);
  double myvalue_neigh_check(asa_objective *asa);
  void mygrad_neigh_check(asa_objective *asa);  
  void allocate();
  
};

}

#endif

/* ERROR/WARNING messages:

E: Non-numeric positions - simulation unstable

UNDOCUMENTED

*/
