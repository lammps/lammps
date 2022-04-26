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
PairStyle(comb,PairComb);
// clang-format on
#else

#ifndef LMP_PAIR_COMB_H
#define LMP_PAIR_COMB_H

#include "pair.h"

namespace LAMMPS_NS {

class PairComb : public Pair {
 public:
  PairComb(class LAMMPS *);
  ~PairComb() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;

  virtual double yasu_char(double *, int &);
  double enegtot;

  static constexpr int NPARAMS_PER_LINE = 49;

 protected:
  struct Param {
    double lam11, lam12, lam21, lam22;
    double c, d, h;
    double gamma, powerm;
    double powern, beta;
    double biga1, biga2, bigb1, bigb2;
    double bigd, bigr;
    double cut, cutsq;
    double c1, c2, c3, c4;
    double plp1, plp3, plp6, a123, aconf;
    double rlm1, rlm2;
    double romiga, romigb, romigc, romigd, addrep;
    double QU1, QL1, DU1, DL1, Qo1, dQ1, aB1, bB1, nD1, bD1;
    double QU2, QL2, DU2, DL2, Qo2, dQ2, aB2, bB2, nD2, bD2;
    double chi, dj, dk, dl, dm, esm1, esm2, cmn1, cmn2, cml1, cml2;
    double coulcut, lcut, lcutsq, hfocor;
    int ielement, jelement, kelement;
    int powermint;
  };

  double cutmax;    // max cutoff for all elements
  double precision;
  Param *params;    // parameter set for an I-J-K interaction

  int nmax;
  double *qf;

  double *esm, **fafb, **dfafb, **ddfafb, **phin, **dphin, **erpaw;
  double *charge;
  int **intype, *typeno;
  int *NCo, cor_flag, cuo_flag, cuo_flag1, cuo_flag2;
  double **bbij;

  int pgsize;                   // size of neighbor page
  int oneatom;                  // max # of neighbors for one atom
  int *sht_num, **sht_first;    // short-range neighbor list
  MyPage<int> *ipage;           // neighbor list pages
  double cutmin;

  void allocate();
  virtual void read_file(char *);
  void setup_params();
  virtual void repulsive(Param *, double, double &, int, double &, double, double);
  double zeta(Param *, double, double, double *, double *);
  void force_zeta(Param *, int, int, int, double, double, double, double, double &, double &,
                  double &);
  void attractive(Param *, double, double, double, double *, double *, double *, double *,
                  double *);
  double elp(Param *, double, double, double *, double *);
  void flp(Param *, double, double, double *, double *, double *, double *, double *);
  double comb_fc(double, Param *);
  double comb_fc_d(double, Param *);
  double comb_fc2(double);
  double comb_fc2_d(double);
  double comb_fc3(double);
  double comb_fc3_d(double);
  virtual double comb_fa(double, Param *, double, double);
  virtual double comb_fa_d(double, Param *, double, double);
  double comb_bij(double, Param *);
  double comb_bij_d(double, Param *);

  inline double comb_gijk(const double costheta, const Param *const param) const
  {
    const double comb_c = param->c * param->c;
    const double comb_d = param->d * param->d;
    const double hcth = param->h - costheta;

    return param->gamma * (1.0 + comb_c / comb_d - comb_c / (comb_d + hcth * hcth));
  }

  inline double comb_gijk_d(const double costheta, const Param *const param) const
  {
    const double comb_c = param->c * param->c;
    const double comb_d = param->d * param->d;
    const double hcth = param->h - costheta;
    const double numerator = -2.0 * comb_c * hcth;
    const double denominator = 1.0 / (comb_d + hcth * hcth);
    return param->gamma * numerator * denominator * denominator;
  }

  void comb_zetaterm_d(double, double *, double, double *, double, double *, double *, double *,
                       Param *);
  void costheta_d(double *, double, double *, double, double *, double *, double *);
  double self(Param *, double, double);
  void sm_table();
  void potal_calc(double &, double &, double &);
  void tri_point(double, int &, int &, int &, double &, double &, double &, int &);
  void direct(int, int, int, int, double, double, double, double, double, double, double, double,
              double, double &, double &);
  void field(Param *, double, double, double, double &, double &);
  double qfo_self(Param *, double, double);
  void qfo_short(Param *, int, int, double, double, double, double &, double &);
  void qfo_direct(int, int, int, int, double, double, double, double, double, double &);
  void qfo_field(Param *, double, double, double, double &, double &);
  void qsolve(double *);
  void Over_cor(Param *, double, int, double &, double &);
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  void Short_neigh();
};

}    // namespace LAMMPS_NS

#endif
#endif
