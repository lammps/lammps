/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(tersoff,PairTersoff);
// clang-format on
#else

#ifndef LMP_PAIR_TERSOFF_H
#define LMP_PAIR_TERSOFF_H

#include "pair.h"

namespace LAMMPS_NS {

class PairTersoff : public Pair {
 public:
  PairTersoff(class LAMMPS *);
  ~PairTersoff() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  template <int SHIFT_FLAG, int EVFLAG, int EFLAG, int VFLAG_ATOM> void eval();

  static constexpr int NPARAMS_PER_LINE = 17;

  struct Param {
    double lam1, lam2, lam3;
    double c, d, h;
    double gamma, powerm;
    double powern, beta;
    double biga, bigb, bigd, bigr;
    double cut, cutsq;
    double c1, c2, c3, c4;
    int ielement, jelement, kelement;
    int powermint;
    double Z_i, Z_j;    // added for TersoffZBL
    double ZBLcut, ZBLexpscale;
    double c5, ca1, ca4;    // added for TersoffMOD
    double powern_del;
    double c0;    // added for TersoffMODC
  };

 protected:
  Param *params;      // parameter set for an I-J-K interaction
  double cutmax;      // max cutoff for all elements
  int maxshort;       // size of short neighbor list array
  int *neighshort;    // short neighbor list array

  int shift_flag;    // flag to turn on/off shift
  double shift;      // negative change in equilibrium bond length

  virtual void allocate();
  virtual void read_file(char *);
  virtual void setup_params();
  virtual void repulsive(Param *, double, double &, int, double &);
  virtual double zeta(Param *, double, double, double *, double *);
  virtual void force_zeta(Param *, double, double, double &, double &, int, double &);
  void attractive(Param *, double, double, double, double *, double *, double *, double *,
                  double *);

  virtual double ters_fc(double, Param *);
  virtual double ters_fc_d(double, Param *);
  virtual double ters_fa(double, Param *);
  virtual double ters_fa_d(double, Param *);
  virtual double ters_bij(double, Param *);
  virtual double ters_bij_d(double, Param *);

  virtual void ters_zetaterm_d(double, double *, double, double, double *, double, double, double *,
                               double *, double *, Param *);
  void costheta_d(double *, double, double *, double, double *, double *, double *);

  // inlined functions for efficiency

  inline double ters_gijk(const double costheta, const Param *const param) const
  {
    const double ters_c = param->c * param->c;
    const double ters_d = param->d * param->d;
    const double hcth = param->h - costheta;

    return param->gamma * (1.0 + ters_c / ters_d - ters_c / (ters_d + hcth * hcth));
  }

  inline double ters_gijk_d(const double costheta, const Param *const param) const
  {
    const double ters_c = param->c * param->c;
    const double ters_d = param->d * param->d;
    const double hcth = param->h - costheta;
    const double numerator = -2.0 * ters_c * hcth;
    const double denominator = 1.0 / (ters_d + hcth * hcth);
    return param->gamma * numerator * denominator * denominator;
  }
};

}    // namespace LAMMPS_NS

#endif
#endif
