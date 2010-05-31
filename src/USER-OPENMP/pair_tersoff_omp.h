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

PairStyle(tersoff/omp,PairTersoffOMP)

#else

#ifndef LMP_PAIR_TERSOFF_OMP_H
#define LMP_PAIR_TERSOFF_OMP_H

#include "pair_omp.h"

namespace LAMMPS_NS {

class PairTersoffOMP : public PairOMP {
 public:
  PairTersoffOMP(class LAMMPS *);
  virtual ~PairTersoffOMP();
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

  double ters_fc_d(double, Param *);
  virtual double ters_fa(double, Param *);
  virtual double ters_fa_d(double, Param *);
  double ters_bij(double, Param *);
  double ters_bij_d(double, Param *);

  void ters_zetaterm_d(double, double *, double, double *, double,
			       double *, double *, double *, Param *);
  void costheta_d(double *, double, double *, double,
		  double *, double *, double *);

  // small methods that are called a *lot* and have no side effects.
  inline double ters_fc(const double r, const Param *param) const
    {
      double ters_R = param->bigr;
      double ters_D = param->bigd;

      if (r < ters_R-ters_D) return 1.0;
      if (r > ters_R+ters_D) return 0.0;
      return 0.5*(1.0 - sin(PI2*(r - ters_R)/ters_D));
    };

  inline double ters_gijk(const double costheta, const Param *param) const {
    const double ters_c = param->c;
    const double ters_d = param->d;
    const double ters_r = ters_c/ters_d;
    const double ters_h = param->h - costheta;

    return param->gamma*(1.0 + (ters_r*ters_r) - (ters_c*ters_c)
			 / ((ters_d*ters_d) + (ters_h*ters_h)));
  };

  inline double ters_gijk_d(const double costheta, const Param *param) const {

    const double ters_c = param->c;
    const double ters_d = param->d;
    const double ters_h = param->h - costheta;

    const double numerator = -2.0 * (ters_c*ters_c) * ters_h;
    const double denom = (ters_d*ters_d) + (ters_h*ters_h);

    return param->gamma*numerator/(denom*denom);
  };

  // vector functions, inline for efficiency

  inline double vec3_dot(double *x, double *y) const {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  }
  inline void vec3_add(double *x, double *y, double *z) const {
    z[0] = x[0]+y[0];  z[1] = x[1]+y[1];  z[2] = x[2]+y[2];
  }
  inline void vec3_scale(double k, double *x, double *y) const {
    y[0] = k*x[0];  y[1] = k*x[1];  y[2] = k*x[2];
  }
  inline void vec3_scaleadd(double k, double *x, double *y, double *z) const {
    z[0] = k*x[0]+y[0];  z[1] = k*x[1]+y[1];  z[2] = k*x[2]+y[2];
  }
};

}

#endif
#endif
