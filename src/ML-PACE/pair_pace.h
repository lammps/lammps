/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*
Copyright 2021 Yury Lysogorskiy^1, Cas van der Oord^2, Anton Bochkarev^1,
 Sarath Menon^1, Matteo Rinaldi^1, Thomas Hammerschmidt^1, Matous Mrovec^1,
 Aidan Thompson^3, Gabor Csanyi^2, Christoph Ortner^4, Ralf Drautz^1

^1: Ruhr-University Bochum, Bochum, Germany
^2: University of Cambridge, Cambridge, United Kingdom
^3: Sandia National Laboratories, Albuquerque, New Mexico, USA
^4: University of British Columbia, Vancouver, BC, Canada
*/

//
// Created by Lysogorskiy Yury on 27.02.20.
//

#ifdef PAIR_CLASS
// clang-format off
PairStyle(pace,PairPACE);
// clang-format on
#else

#ifndef LMP_PAIR_PACE_H
#define LMP_PAIR_PACE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairPACE : public Pair {
 public:
  PairPACE(class LAMMPS *);
  ~PairPACE() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  void *extract(const char *, int &) override;
  void *extract_peratom(const char *, int &) override;

 protected:
  struct ACEImpl *aceimpl;
  int nmax_corerep = 0;

  virtual void allocate();
  double *corerep_factor;                    //per-atom core-rep factor (= 1 - fcut)
  int flag_corerep_factor;

  double **scale;
  bool recursive;    // "recursive" option for ACERecursiveEvaluator

  int chunksize;
};
}    // namespace LAMMPS_NS

#endif
#endif
