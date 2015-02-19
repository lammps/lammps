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

#ifdef PAIR_CLASS

PairStyle(coul/streitz,PairCoulStreitz)

#else

#ifndef LMP_PAIR_COUL_Streitz_H
#define LMP_PAIR_COUL_Streitz_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCoulStreitz : public Pair {
 public:
  PairCoulStreitz(class LAMMPS *);
  virtual ~PairCoulStreitz();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double memory_usage();
  virtual void *extract(const char *, int &);

 protected:
  struct Param {
    double chi, eta, gamma, zeta, zcore;
    int ielement;
  };

  int nmax;
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  int *elem2param;            // mapping from element triplets to parameters
  int *map;                     // mapping from atom types to elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  double precision;
  Param *params;                // parameter set for an I-J-K interaction

  // Kspace parameters
  int kspacetype;
  double cut_coul, cut_coulsq;
  double *cut_respa;
  double **scale;
  
  // Wolf
  double g_wolf, woself, dwoself;
  
  // Ewald
  double g_ewald;

  // QEq
  double *qeq_x, *qeq_j, *qeq_g, *qeq_z, *qeq_c;

  void allocate();
  virtual void read_file(char *);
  void setup();
  double self(Param *, double);
  void coulomb_integral_wolf(double, double, double, double &, double &,
        double &, double &);
  void wolf_sum(double, double, double, double, double, double, double, 
		  double, double &, double &);
  void coulomb_integral_ewald(double, double, double, double &, double &,
        double &, double &);
  void ewald_sum(double, double, double, double, double, double, double, 
		  double, double &, double &, double);

};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
