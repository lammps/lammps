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

PairStyle(e3b,PairE3B)

#else

#ifndef LMP_PAIR_E3B_H
#define LMP_PAIR_E3B_H

#include "pair.h"

namespace LAMMPS_NS {

class PairE3B : public Pair {
 public:
  PairE3B(class LAMMPS *);
  virtual ~PairE3B();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  virtual double init_one(int, int);
  virtual void init_style();

protected:
  //potential parameters
  int typeO;
  double ea,eb,ec;      //three body energies
  double k3;            //three body exponential decay (units inverse length)
  double rs,rc3,rc2;    //rs: switching cutuff, rc3: cutoff for 3-body
  double e2,k2;         //2-body energy and exp decay
  double cutmax;   //max cutoff of all interactions
  double rc2sq,rc3sq,rc3deltaSq;
  double sc_denom,sc_num;

  //list of indexes of Os and Hs in each pair
  int pairmax,pairPerAtom;      // size of pair list
  int **pairO,***pairH;         // pair lists
  double ***exps,****del3,***fpair3,*sumExp;
  int maxID;  //size of global sumExp array
  size_t nbytes;     //size of sumExp array in bytes
  int natoms;        //to make sure number of atoms is constant

  virtual void allocate();
  void allocateE3B();
  bool allocatedE3B;
  //for reading settings from pair_style input
  bool checkKeyword(const char *,const char *,const int, const int);
  void checkInputs(const double &bondL);
  void presetParam(const int flag,bool &repeatFlag,double &bondL);
  tagint find_maxID();
};
}

#endif
#endif
