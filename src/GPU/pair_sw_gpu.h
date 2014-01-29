/* ----------------------------------------------------------------------
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

PairStyle(sw/gpu,PairSWGPU)

#else

#ifndef LMP_PAIR_SW_GPU_H
#define LMP_PAIR_SW_GPU_H

#include "pair_sw.h"

namespace LAMMPS_NS {

class PairSWGPU : public PairSW {
 public:
  PairSWGPU(class LAMMPS *);
  ~PairSWGPU();
  void compute(int, int);
  double init_one(int, int);
  void init_style();

 enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 protected:
  void allocate();

  int gpu_mode;
  double cpu_time;
  int *gpulist;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Insufficient memory on accelerator

UNDOCUMENTED

E: Pair style sw/gpu requires atom IDs

UNDOCUMENTED

E: Pair style sw/gpu requires newton pair off

UNDOCUMENTED

E: All pair coeffs are not set

UNDOCUMENTED

U: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

U: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

U: No matching element in ADP potential file

The ADP potential file does not contain elements that match the
requested elements.

U: Cannot open ADP potential file %s

The specified ADP potential file cannot be opened.  Check that the
path and name are correct.

U: Incorrect element names in ADP potential file

The element names in the ADP file do not match those requested.

*/
