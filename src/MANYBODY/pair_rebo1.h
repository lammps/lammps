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

PairStyle(rebo1,PairREBO1)

#else

#ifndef LMP_PAIR_REBO1_H
#define LMP_PAIR_REBO1_H

#include "pair_tersoff.h"

// See "Humbird and Graves, JCP 120, 2405 (2004) for equations

namespace LAMMPS_NS {

class PairREBO1 : public PairTersoff {
 public:
  PairREBO1(class LAMMPS *);
  virtual ~PairREBO1() {}
  void compute(int, int);

 protected:

  // Library file containing the spline coefficient
  void read_lib();

  // Moliere repulsive, Eq. A7-8
  double moliere_aB; 
  double moliere_c[3];
  double moliere_d[3];
  double firsov_const;

  void read_file(char *);

  // Repulsive term, Eq. A3
  void repulsive(Param *, double, double &, int, double &);

  // Cutoff function, Eq. A6
  virtual double ters_fc(double, Param *);
  virtual double ters_fc_d(double, Param *);

  // Attraction term, Eq. A4
  double ters_fa(double, Param *);
  double ters_fa_d(double, Param *);

  // Zeta term, Eq. A13
  virtual double zeta(Param *, Param *, double, double, double *, double *);
  virtual void force_zeta(Param *, double, double, double &,
                          double &, int, double &, int, int);
  virtual void ters_zetaterm_d(Param *, Param *, double, double *, double, double *, double,
                               double *, double *, double *);

  // Bij term, Eq. A12
  virtual double ters_bij(double, Param *, int, int);
  virtual double ters_bij_d(double, Param *, int, int);

  // Angular term; Eq. A13; inlined functions for efficiency

  inline double ters_gijk(const double cos_theta,
                          const Param * const param) const {
    const double ters_c = param->c;
    const double ters_d = param->d;
    const double hcth = param->h - cos_theta;

    return ters_c + ters_d * hcth * hcth;
  }

  // replace the inline double with void below, since derivative of 
  // cos_theta is not -sin_theta

  /*
  inline double ters_gijk_d(const double cos_theta,
                            const Param * const param) const {
    const double ters_c = param->c;
    const double ters_d = param->d;
    const double hcth = param->h - cos_theta;
    const double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
    return 2.0 * ters_d * hcth * sin_theta;
  }
  */

  void ters_gijk_d(Param*, double *, double, double *, double,
                  double *, double *, double *);

  // Attractive term, Eq. A4
  void attractive(Param *, Param *, double, double, double, double *, double *,
                  double *, double *, double *);

  // Coordination terms, Hij of Eq. A12
  void count_neigh();
  int nmax;
  double coordenergy[5][402], coordforce[5][402], coordnumber[5][402];
  double *NCl, *NSi;

  // communication functions
  int pack_flag;
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int pack_forward_comm(int , int *, double *, int, int *);
  void unpack_forward_comm(int , int , double *);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Pair rebo1 requires metal or real units

This is a current restriction of this pair potential.

E: Cannot open Tersoff potential file %s

The specified potential file cannot be opened.  Check that the path
and name are correct.

E: Incorrect format in Tersoff potential file

Incorrect number of words per line in the potential file.

E: Illegal Tersoff parameter

One or more of the coefficients defined in the potential file is
invalid.

*/
