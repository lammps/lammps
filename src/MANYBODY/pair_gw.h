/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(gw,PairGW);
// clang-format on
#else

#ifndef LMP_PAIR_GW_H
#define LMP_PAIR_GW_H

#include "pair.h"

namespace LAMMPS_NS {

class PairGW : public Pair {
 public:
  PairGW(class LAMMPS *);
  ~PairGW() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  static constexpr int NPARAMS_PER_LINE = 17;

 protected:
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
    double Z_i, Z_j;
    double ZBLcut, ZBLexpscale;
  };

  Param *params;    // parameter set for an I-J-K interaction
  double cutmax;    // max cutoff for all elements

  int **pages;     // neighbor list pages
  int maxlocal;    // size of numneigh, firstneigh arrays
  int maxpage;     // # of pages currently allocated
  int pgsize;      // size of neighbor page
  int oneatom;     // max # of neighbors for one atom

  int *GW_numneigh;       // # of pair neighbors for each atom
  int **GW_firstneigh;    // ptr to 1st neighbor of each atom

  void GW_neigh();
  void add_pages(int howmany = 1);

  void allocate();
  virtual void read_file(char *);
  void setup_params();
  virtual void repulsive(Param *, double, double &, int, double &);
  double zeta(Param *, double, double, double *, double *);
  virtual void force_zeta(Param *, double, double, double &, double &, int, double &);
  void attractive(Param *, double, double, double, double *, double *, double *, double *,
                  double *);

  double gw_fc(double, Param *);
  double gw_fc_d(double, Param *);
  virtual double gw_fa(double, Param *);
  virtual double gw_fa_d(double, Param *);
  double gw_bij(double, Param *);
  double gw_bij_d(double, Param *);

  void gw_zetaterm_d(double, double *, double, double *, double, double *, double *, double *,
                     Param *);
  void costheta_d(double *, double, double *, double, double *, double *, double *);

  // inlined functions for efficiency

  inline double gw_gijk(const double costheta, const Param *const param) const
  {
    const double gw_c = param->c * param->c;
    const double gw_d = param->d * param->d;
    const double hcth = param->h - costheta;

    //printf("gw_gijk: gw_c=%f gw_d=%f hcth=%f=%f-%f\n", gw_c, gw_d, hcth, param->h, costheta);

    return param->gamma * (1.0 + gw_c / gw_d - gw_c / (gw_d + hcth * hcth));
  }

  inline double gw_gijk_d(const double costheta, const Param *const param) const
  {
    const double gw_c = param->c * param->c;
    const double gw_d = param->d * param->d;
    const double hcth = param->h - costheta;
    const double numerator = -2.0 * gw_c * hcth;
    const double denominator = 1.0 / (gw_d + hcth * hcth);
    return param->gamma * numerator * denominator * denominator;
  }
};

}    // namespace LAMMPS_NS

#endif
#endif
