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

PairStyle(list,PairList)

#else

#ifndef LMP_PAIR_LIST_H
#define LMP_PAIR_LIST_H

#include "pair.h"

namespace LAMMPS_NS {

class PairList : public Pair {
 public:
  PairList(class LAMMPS *);
  virtual ~PairList();

  virtual void compute(int, int);
  virtual void settings(int, char **);
  virtual void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
  virtual double memory_usage();

 protected:
  void allocate();

  enum { NONE=0, HARM, MORSE, LJ126 };

  // potential specific parameters
  struct harm_p  { double k, r0;          };
  struct morse_p { double d0, alpha, r0;  };
  struct lj126_p { double epsilon, sigma; };

  union parm_u { 
    struct harm_p harm;
    struct morse_p morse;
    struct lj126_p lj126;
  };

  typedef struct {
    int id1,id2;        // global atom ids
    double cutsq;       // cutoff**2 for this pair
    double offset;      // energy offset
    union parm_u parm;  // parameters for style
  } list_parm_t;    

 protected:
  double cut_global;    // global cutoff distance
  int *style;           // list of styles for pair interactions
  list_parm_t *params;  // lisf of pair interaction parameters
  int npairs;           // # of atom pairs in global list
  int check_flag;       // 1 if checking for missing pairs
};

}

#endif
#endif

/* ERROR/WARNING messages:


*/
