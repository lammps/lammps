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
PairStyle(eam/alloy/gpu,PairEAMAlloyGPU);
// clang-format on
#else

#ifndef LMP_PAIR_EAM_ALLOY_GPU_H
#define LMP_PAIR_EAM_ALLOY_GPU_H

#include "pair_eam.h"

namespace LAMMPS_NS {

class PairEAMAlloyGPU : public PairEAM {
 public:
  PairEAMAlloyGPU(class LAMMPS *);
  virtual ~PairEAMAlloyGPU();
  void coeff(int, char **);
  void compute(int, int);
  void init_style();
  double single(int, int, int, int, double, double, double, double &);
  double memory_usage();
  void *extract(const char *, int &) { return nullptr; }

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);

  enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 protected:
  void read_file(char *);
  void file2array();

  int gpu_mode;
  double cpu_time;
  void *fp_pinned;
  bool fp_single;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Insufficient memory on accelerator

There is insufficient memory on one of the devices specified for the gpu
package

E: Cannot use newton pair with eam/alloy/gpu pair style

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: No matching element in EAM potential file

The EAM potential file does not contain elements that match the
requested elements.

E: Cannot open EAM potential file %s

The specified EAM potential file cannot be opened.  Check that the
path and name are correct.

E: Incorrect element names in EAM potential file

The element names in the EAM file do not match those requested.

*/
