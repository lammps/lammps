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
PairStyle(polymorphic,PairPolymorphic);
// clang-format on
#else

#ifndef LMP_PAIR_POLYMORPHIC_H
#define LMP_PAIR_POLYMORPHIC_H

#include "pair.h"

namespace LAMMPS_NS {
// forward declaration
class TabularFunction;

class PairPolymorphic : public Pair {
 public:
  PairPolymorphic(class LAMMPS *);
  ~PairPolymorphic() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

 protected:
  struct PairParameters {
    double cut;
    double cutsq;
    double xi;
    TabularFunction *U;
    TabularFunction *V;
    TabularFunction *W;
    TabularFunction *F;
    PairParameters();
    ~PairParameters();
  };

  struct TripletParameters {
    TabularFunction *P;
    TabularFunction *G;
    TripletParameters();
    ~TripletParameters();
  };

  double epsilon;
  int eta;
  int nx, nr, ng;    // table sizes
  double maxX;

  // parameter sets
  PairParameters *pairParameters;          // for I-J interaction
  TripletParameters *tripletParameters;    // for I-J-K interaction

  int neighsize, numneighV, numneighW, numneighW1;
  int *firstneighV, *firstneighW, *firstneighW1;
  double *delxV, *delyV, *delzV, *drV;
  double *delxW, *delyW, *delzW, *drW;

  double cutmax;    // max cutoff for all elements
  double cutmaxsq;
  int npair, ntriple;
  int *match;

  void allocate();

  virtual void read_file(char *);
  void setup_params();
#if defined(LMP_POLYMORPHIC_WRITE_TABLES)
  void write_tables(int);
#endif
  void attractive(PairParameters *, PairParameters *, TripletParameters *, double, double, double,
                  double *, double *, double *, double *, double *);

  void ters_zetaterm_d(double, double *, double, double *, double, double *, double *, double *,
                       PairParameters *, PairParameters *, TripletParameters *);
  void costheta_d(double *, double, double *, double, double *, double *, double *);
};
}    // namespace LAMMPS_NS

#endif
#endif
