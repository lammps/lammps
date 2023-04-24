/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*
Copyright 2022 Yury Lysogorskiy^1, Anton Bochkarev^1, Matous Mrovec^1, Ralf Drautz^1

^1: Ruhr-University Bochum, Bochum, Germany
 */

//
// Created by Lysogorskiy Yury on 1.01.22.
//

#ifdef PAIR_CLASS
// clang-format off
PairStyle(pace/extrapolation,PairPACEExtrapolation)
// clang-format on
#else

#ifndef LMP_PAIR_PACE_AL_H
#define LMP_PAIR_PACE_AL_H

#include "pair.h"
#include <vector>

namespace LAMMPS_NS {

class PairPACEExtrapolation : public Pair {
 public:
  PairPACEExtrapolation(class LAMMPS *);
  ~PairPACEExtrapolation() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void *extract(const char *, int &) override;
  void *extract_peratom(const char *, int &) override;

 protected:
  struct ACEALImpl *aceimpl;
  int nmax;

  virtual void allocate();
  std::vector<std::string> element_names;    // list of elements (used by dump pace/extrapolation)
  double *extrapolation_grade_gamma;         //per-atom gamma value

  int flag_compute_extrapolation_grade;

  double **scale;

  int chunksize;
};

}    // namespace LAMMPS_NS

#endif
#endif
