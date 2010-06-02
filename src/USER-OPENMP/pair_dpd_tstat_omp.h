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

PairStyle(dpd/tstat/omp,PairDPDTstatOMP)

#else

#ifndef LMP_PAIR_DPD_TSTAT_OMP_H
#define LMP_PAIR_DPD_TSTAT_OMP_H

#include "pair_dpd_omp.h"

namespace LAMMPS_NS {

class PairDPDTstatOMP : public PairDPDOMP {
 public:
  PairDPDTstatOMP(class LAMMPS *);
  ~PairDPDTstatOMP(){};
  virtual void compute(int, int);
  virtual void settings(int, char **);
  virtual void coeff(int, char **);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);

 protected:
  double t_start,t_stop;
  template <int, int, int> void eval_tstat();
};

}

#endif
#endif
