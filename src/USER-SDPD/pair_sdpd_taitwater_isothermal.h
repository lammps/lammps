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
PairStyle(sdpd/taitwater/isothermal,PairSDPDTaitwaterIsothermal);
// clang-format on
#else

#ifndef LMP_PAIR_SDPD_TAITWATER_MORRIS_ISOTHERMAL_H
#define LMP_PAIR_SDPD_TAITWATER_MORRIS_ISOTHERMAL_H

#include "pair.h"
#ifdef USE_ZEST
#include "zest.hpp"
#include <random>
#endif

namespace LAMMPS_NS {

class PairSDPDTaitwaterIsothermal : public Pair {
 public:
  PairSDPDTaitwaterIsothermal(class LAMMPS *);
  virtual ~PairSDPDTaitwaterIsothermal();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);
  virtual void init_style();

 protected:
  double viscosity, temperature;
  double *rho0, *soundspeed, *B;
  double **cut;

  void allocate();

  unsigned int seed;
#ifdef USE_ZEST
  std::mt19937_64 generator;
  Ziggurat<zest::StandardNormal, std::mt19937_64> gaussian;
#else
  class RanMars *random;
#endif
};

}    // namespace LAMMPS_NS

#endif
#endif
