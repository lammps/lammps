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

#ifndef PAIR_TERSOFF_H
#define PAIR_TERSOFF_H

#include "pair.h"

namespace LAMMPS_NS {

class PairTersoff : public Pair {
 public:
  PairTersoff(class LAMMPS *);
  virtual ~PairTersoff();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

 protected:
  struct Param {
    double lam1,lam2,lam3;
    double c,d,h;
    double gamma,powerm;
    double powern,beta;
    double biga,bigb,bigd,bigr;
    double cut,cutsq;
    double c1,c2,c3,c4;
    int ielement,jelement,kelement;
    int powermint;
    double Z_i,Z_j;
    double ZBLcut,ZBLexpscale;
  };
  
  double PI,PI2,PI4;
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  int ***elem2param;            // mapping from element triplets to parameters
  int *map;                     // mapping from atom types to elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  Param *params;                // parameter set for an I-J-K interaction

  void allocate();
  virtual void read_file(char *);
  void setup();
  virtual void repulsive(Param *, double, double &, int, double &);
  double zeta(Param *, double, double, double *, double *);
  void force_zeta(Param *, double, double, double &, double &, int, double &);
  void attractive(Param *, double, double, double, double *, double *,
		  double *, double *, double *);

  double ters_fc(double, Param *);
  double ters_fc_d(double, Param *);
  virtual double ters_fa(double, Param *);
  virtual double ters_fa_d(double, Param *);
  double ters_bij(double, Param *);
  double ters_bij_d(double, Param *);
  double ters_gijk(double, Param *);
  double ters_gijk_d(double, Param *);
  void ters_zetaterm_d(double, double *, double, double *, double,
			       double *, double *, double *, Param *);
  void costheta_d(double *, double, double *, double,
		  double *, double *, double *);

  // vector functions, inline for efficiency

  inline double vec3_dot(double *x, double *y) {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  }
  inline void vec3_add(double *x, double *y, double *z) {
    z[0] = x[0]+y[0];  z[1] = x[1]+y[1];  z[2] = x[2]+y[2];
  }
  inline void vec3_scale(double k, double *x, double *y) {
    y[0] = k*x[0];  y[1] = k*x[1];  y[2] = k*x[2];
  }
  inline void vec3_scaleadd(double k, double *x, double *y, double *z) {
    z[0] = k*x[0]+y[0];  z[1] = k*x[1]+y[1];  z[2] = k*x[2]+y[2];
  }
};

}

#endif
