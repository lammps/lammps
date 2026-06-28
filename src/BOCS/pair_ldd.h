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
/* ------------------------------------------------------
    This file is part of the BOCS package for LAMMPS.
    Contributed by Michael R. DeLyser, mrd5285@psu.edu
    and Maria C. Lesniewski, mjl6766@psu.edu
    The Pennsylvania State University
   ------------------------------------------------------ */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(ldd,PairLdd);
// clang-format on
#else

#ifndef LMP_PAIR_LDD_H
#define LMP_PAIR_LDD_H

#include "pair.h"

namespace LAMMPS_NS {

//Forward Declarations
class LddIndicator;
class LddPotential;

// the registration tables below use LMP_REGISTRY_CONST (see lmptype.h): const,
// except in a GPU-enabled Kokkos build where the host-only factory function
// pointers must not be shadowed into device memory.

typedef LddIndicator *(*IndicatorCreator)(LAMMPS *);
struct LddIndicatorInfo {
  const char *name;
  IndicatorCreator creator;
};
extern LMP_REGISTRY_CONST LddIndicatorInfo ldd_indicator_table[];
extern LMP_REGISTRY_CONST int num_ldd_indicator;

typedef LddPotential *(*PotentialCreator)(LAMMPS *);
struct LddPotentialInfo {
  const char *name;
  PotentialCreator creator;
};
extern LMP_REGISTRY_CONST LddPotentialInfo ldd_potential_table[];
extern LMP_REGISTRY_CONST int num_ldd_potential;

class PairLdd : public Pair {
 public:
  PairLdd(class LAMMPS *);
  ~PairLdd() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void *extract_peratom(const char *, int &) override;
  void coeff_ldd(int si, int sj, int narg, char **arg);

  class LddIndicator *new_indicator(const std::string &);
  class LddPotential *new_potential(const std::string &);

  // Functions to calculate the local densities, the gradients,
  // and the associated energies
  void LDD_calculate_LDs();
  void LDD_calculate_energies();

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 protected:
  // all per-species-pair settings below are indexed by species (0..nelements-1);
  // atom types are mapped to species via the inherited map[] array.
  bool **self_interaction;    // self-interaction flag for each species pair
  bool **ignore_pair;         // ignored species pairs
  bool *ignore_me;            // species with no active interaction as a central atom
  bool **bGradient;           // species pairs that also carry a gradient interaction

  LddIndicator ***Inds;         // indicator function w(r) for each ordered species pair
  LddPotential ***Potls;        // U_rho potential for each ordered species pair
  LddPotential ***GradPotls;    // U_grad potential for each ordered species pair

  // per-atom local-density data owned by the pair style (recomputed every step,
  // not stored in data/restart files).  Communicated via the pack/unpack methods
  // below and exposed to the rest of LAMMPS through extract_peratom() / fix pair.
  int nmax;                   // current allocated length of the per-atom arrays
  double **local_density;     // local density of each type around an atom
  double **grad_density;      // gradient of local density (3 components per type)
  double **ld_energy;         // u_{b|t_I}(rho_I) per type
  double **ld_grad_energy;    // u_{\nabla b|t_I}(rho_I) per type
  double *total_energy;       // sum of all local-density energies on an atom

  void allocate();
  void allocate_species();
  void grow_peratom();
  void ErrorDoubleKeyword(const char *);
  void ErrorNumKeywordArgs(const char *, const char *);
  void read_file(char *filename);    // reads ldd potential file, executes coeff_ldd per entry

};

}    // namespace LAMMPS_NS

#endif
#endif
