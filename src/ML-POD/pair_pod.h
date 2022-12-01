/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(pod,PairPOD);
// clang-format on
#else

#ifndef LMP_PAIR_POD_H
#define LMP_PAIR_POD_H

#include "pair.h"

namespace LAMMPS_NS {

class PairPOD : public Pair {
 public:
  PairPOD(class LAMMPS *);
  ~PairPOD() override;
  void compute(int, int) override;

  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;

  int dim;    // typically 3

  double *gd;             // global linear descriptors
  double *gdall;          // global linear descriptors summed over all MPI ranks
  double *podcoeff;       // POD coefficients
  double *newpodcoeff;    // normalized POD coefficients
  double *energycoeff;    // energy coefficients
  double *forcecoeff;     // force coefficients

  void estimate_tempmemory();
  void free_tempmemory();
  void allocate_tempmemory();

  void lammpsNeighPairs(double **x, int **firstneigh, int *atomtype, int *map, int *numneigh,
                        int i);

 protected:
  int nablockmax;    // maximum number of atoms per computation block
  int nij;           //  number of atom pairs
  int nijmax;        // maximum number of atom pairs
  int szd;           // size of tmpmem

  class MLPOD *podptr;

  // temporary arrays for computation blocks

  double *tmpmem;      // temporary memory
  int *typeai;         // types of atoms I only
  int *numneighsum;    // cumulative sum for an array of numbers of neighbors
  double *rij;         // (xj - xi) for all pairs (I, J)
  int *idxi;           // storing linear indices for all pairs (I, J)
  int *ai;             // IDs of atoms I for all pairs (I, J)
  int *aj;             // IDs of atoms J for all pairs (I, J)
  int *ti;             // types of atoms I for all pairs (I, J)
  int *tj;             // types of atoms J  for all pairs (I, J)

  bool peratom_warn;    // print warning about missing per-atom energies or stresses
};

}    // namespace LAMMPS_NS

#endif
#endif
