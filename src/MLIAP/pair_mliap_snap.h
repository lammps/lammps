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

PairStyle(mliap/snap,PairMLIAPSNAP)

#else

#ifndef LMP_PAIR_MLIAP_SNAP_H
#define LMP_PAIR_MLIAP_SNAP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairMLIAPSNAP : public Pair {
public:
  PairMLIAPSNAP(class LAMMPS *);
  ~PairMLIAPSNAP();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
  virtual double memory_usage();

  double rcutfac, quadraticflag; // declared public to workaround gcc 4.9
  int ncoeff;                    //  compiler bug, manifest in KOKKOS package

protected:
  virtual void allocate();
  void read_files(char *, char *);
  inline int equal(double* x,double* y);
  inline double dist2(double* x,double* y);

  void compute_bispectrum();

  double* atomenergy;           // energies for all atoms in list
  double** beta;                // betas for all atoms in list
  double** bispectrum;          // bispectrum components for all atoms in list
  int beta_max;                 // length of beta
  class MLIAPModelSNAP* model;
  class MLIAPDescriptorSNAP* descriptor;
};

}

#endif
#endif

