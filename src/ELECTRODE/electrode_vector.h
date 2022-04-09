/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meißner (TUHH)
------------------------------------------------------------------------- */

#include "pointers.h"

#include <map>

namespace LAMMPS_NS {

class ElectrodeVector : protected Pointers {
 public:
  ElectrodeVector(class LAMMPS *, int, double);
  ~ElectrodeVector();
  void setup(const std::map<tagint, int> &, class Pair *, class NeighList *);
  void compute_vector();
  double *vector;
  int igroup;

 private:
  int groupbit;
  bigint ngroup;
  double **cutsq;
  double g_ewald, eta;
  std::map<tagint, int> tag_to_iele;
  std::vector<bigint> mpos;
  class Pair *pair;
  class NeighList *list;
  class ElectrodeKSpace *electrode_kspace;

  void update_mpos();

  void pair_contribution();
  double calc_erfc(double);

  double setup_time_total;
  double reduce_time_total;
  double kspace_time_total;
  double pair_time_total;
  double boundary_time_total;
  double b_time_total;
  double alloc_time_total;
  double mpos_time_total;
};

}    // namespace LAMMPS_NS
