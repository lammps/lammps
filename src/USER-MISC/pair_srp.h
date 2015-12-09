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

PairStyle(srp,PairSRP)

#else

#ifndef LMP_PAIR_SRP_H
#define LMP_PAIR_SRP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSRP : public Pair {
 public:
  PairSRP(class LAMMPS *);
  virtual ~PairSRP();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  virtual void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  virtual void write_data(FILE *);
  virtual void write_data_all(FILE *);
  virtual void write_restart(FILE *);
  virtual void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);

 private:
  inline void onetwoexclude(int* &, int &, int* &, int* &, int** &);
  inline void remapBonds(int &);
  void allocate();
  void getMinDist(double** &, double &, double &, double &, double &,
                  double &, int &, int &, int &, int &);
  bool min, midpoint;
  double **cut;
  double **a0;
  double **srp;
  double cut_global;
  int bptype;
  int btype;
  class Fix *f_srp;
  char *fix_id;
  int exclude,maxcount;
  int **segment;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
