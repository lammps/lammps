/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Paul Lafourcade (CEA-DAM-DIF, Arpajon, France)
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(slcsa/atom,ComputeSLCSAAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_SLCSA_ATOM_H
#define LMP_COMPUTE_SLCSA_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSLCSAAtom : public Compute {
 public:
  ComputeSLCSAAtom(class LAMMPS *, int, char **);
  ~ComputeSLCSAAtom() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  //  double memory_usage() override;
  int compute_ncomps(int);
  int argmax(double *, int);

 private:
  struct value_t {
    int which;         // type of data: COMPUTE, FIX, VARIABLE
    int argindex;      // 1-based index if data is vector, else 0
    std::string id;    // compute/fix/variable ID
    union {
      class Compute *c;
      class Fix *f;
      int v;
    } val;
  };
  value_t descriptorval;
  int nmax;
  int ncols;
  int nevery;
  int ncomps;
  int nclasses;
  const char *database_mean_descriptor_file;
  const char *lda_scalings_file;
  const char *lr_decision_file;
  const char *lr_bias_file;
  const char *covmat_file;
  const char *maha_file;
  class NeighList *list;

  // LDA dimension reduction
  double **lda_scalings;
  double *database_mean_descriptor;

  // LR classification
  double *lr_bias;
  double **lr_decision;

  // Mahalanobis distance calculation
  double ***icov_list;
  double **mean_projected_descriptors;
  double *maha_thresholds;

  // Per-atom local arrays
  double *full_descriptor;
  double *projected_descriptor;
  double *scores;
  double *probas;
  double *prodright;
  double *dmaha;

  // Output array
  double **classification;
};

}    // namespace LAMMPS_NS

#endif
#endif
