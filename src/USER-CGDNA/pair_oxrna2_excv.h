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

PairStyle(oxrna2/excv,PairOxrna2Excv)

#else

#ifndef LMP_PAIR_OXRNA2_EXCV_H
#define LMP_PAIR_OXRNA2_EXCV_H

#include "pair_oxdna_excv.h"

namespace LAMMPS_NS {

class PairOxrna2Excv : public PairOxdnaExcv {
 public:
  PairOxrna2Excv(class LAMMPS *);
  virtual ~PairOxrna2Excv();
  virtual void compute_interaction_sites(double *, double *,
    double *, double *, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

*/
