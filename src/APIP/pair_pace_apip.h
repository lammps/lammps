/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*
This file is a modified version of src/ML-PACE/pair_pace.h.

Original copyright:
Copyright 2021 Yury Lysogorskiy^1, Cas van der Oord^2, Anton Bochkarev^1,
 Sarath Menon^1, Matteo Rinaldi^1, Thomas Hammerschmidt^1, Matous Mrovec^1,
 Aidan Thompson^3, Gabor Csanyi^2, Christoph Ortner^4, Ralf Drautz^1

^1: Ruhr-University Bochum, Bochum, Germany
^2: University of Cambridge, Cambridge, United Kingdom
^3: Sandia National Laboratories, Albuquerque, New Mexico, USA
^4: University of British Columbia, Vancouver, BC, Canada
*/

//
// Originally created by Lysogorskiy Yury on 27.02.20.
// Adaptive precision added by David Immel in 2025
// (d.immel@fz-juelich.de, FZJ, Germany).
//

#ifdef PAIR_CLASS
// clang-format off
PairStyle(pace/apip,PairPACEAPIP);
// clang-format on
#else

#ifndef LMP_PAIR_PACE_APIP_H
#define LMP_PAIR_PACE_APIP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairPACEAPIP : public Pair {
 public:
  PairPACEAPIP(class LAMMPS *);
  ~PairPACEAPIP() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void setup() override;

  void *extract(const char *, int &) override;
  void *extract_peratom(const char *, int &) override;

 protected:
  struct ACEImpl *aceimpl;
  int nmax_corerep;

  virtual void allocate();
  double *corerep_factor;    //per-atom core-rep factor (= 1 - fcut)
  int flag_corerep_factor;

  double **scale;
  bool recursive;    // "recursive" option for ACERecursiveEvaluator

  int chunksize;

  // start of adaptive-precision modifications by DI
  virtual double *get_e_ref_ptr();
  virtual double compute_factor_lambda(double);
  virtual int check_abort_condition(double *, double *, int *, int);

  bool lambda_thermostat;    // true/false there is one/no fix lambda_thermostat

  void calculate_time_per_atom();

  // stats required for load balancing
  int n_computations_accumulated;    // number of accumulated computations
  double time_wall_accumulated;      // accumulated compute time
  double time_per_atom;              // average time of one computation
  // end of adaptive-precision modifications by DI
};
}    // namespace LAMMPS_NS

#endif
#endif
