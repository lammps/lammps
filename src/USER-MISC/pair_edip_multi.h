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

PairStyle(edip/multi,PairEDIPMulti)

#else

#ifndef LMP_PAIR_EDIP_MULTI_H
#define LMP_PAIR_EDIP_MULTI_H

#include "pair.h"

namespace LAMMPS_NS {

class PairEDIPMulti : public Pair {
 public:
  PairEDIPMulti(class LAMMPS *);
  virtual ~PairEDIPMulti();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();

 protected:
  struct Param {
    double A, B;//coefficients for pair interaction I-J
    double cutoffA;//cut-off distance for pair interaction I-J
    double cutoffC;//lower cut-off distance for calculating Z_I
    double alpha;//coefficient for calculating Z_I
    double beta;//attractive term for pair I-J
    double sigma;//cut-off coefficient for pair I-J
    double rho;//pair I-J
    double gamma;//coefficient for three-body interaction I-J-K
    double eta, lambda;//coefficients for function h(l,Z)
    double mu, Q0;//coefficients for function Q(Z)
    double u1, u2, u3, u4;//coefficients for function tau(Z)
    double cutsq;
    int ielement,jelement,kelement;
  };

  double *preForceCoord;

  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  int ***elem2param;            // mapping from element triplets to parameters
  int *map;                     // mapping from atom types to elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  Param *params;                // parameter set for an I-J-K interaction

  // max number of interaction per atom for f(Z) environment potential

  static const int leadDimInteractionList = 64;

  void allocate();
  void allocatePreLoops(void);
  void deallocatePreLoops(void);

  void read_file(char *);
  void setup();

  void edip_pair(double, double, Param *, double &, double &, double &);
  void edip_fc(double, Param *, double &, double &);
  void edip_fcut2(double, Param *, double &, double &);
  void edip_tau(double, Param *, double &, double &);
  void edip_h(double, double, Param *, double &, double &, double &);
  void edip_fcut3(double, Param *, double &, double &);

  double vec3_dot(double x[3], double y[3])
  {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  }

  void vec3_add(double k1, double x[3], double k2, double y[3], double *z)
  {
    z[0] = k1 * x[0] + k2 * y[0];
    z[1] = k1 * x[1] + k2 * y[1];
    z[2] = k1 * x[2] + k2 * y[2];
  }

  //dr_ij=r_j - r_i
  //dr_ik=r_k - r_i
  void costheta_d(double *dr_ij, double r_ij, double *dr_ik, double r_ik,
                  double *dri, double *drj, double *drk)
  {
    double costheta;

    costheta = vec3_dot(dr_ij, dr_ik) / r_ij / r_ik;
    vec3_add(1 / r_ij / r_ik, dr_ik, -costheta / r_ij / r_ij, dr_ij, drj);
    vec3_add(1 / r_ij / r_ik, dr_ij, -costheta / r_ik / r_ik, dr_ik, drk);
    vec3_add(-1, drj, -1, drk, dri);
  }

};

}

#endif
#endif
