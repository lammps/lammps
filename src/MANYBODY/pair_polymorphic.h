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
  virtual ~PairPolymorphic();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

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

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style polymorphic requires atom IDs

This is a requirement to use the polymorphic potential.

E: Pair style polymorphic requires newton pair on

See the newton command.  This is a restriction to use the polymorphic
potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open polymorphic potential file %s

The specified polymorphic potential file cannot be opened.  Check that
the path and name are correct.

E: Incorrect number of elements in potential file

Self-explanatory.

E: Element not defined in potential file

The specified element is not in the potential file.

E: Potential file incompatible with this pair style version

UNDOCUMENTED

E: Error reading potential file header

UNDOCUMENTED

*/
