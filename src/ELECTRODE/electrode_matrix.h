/* ----------------------------------------------------------------------
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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meißner (TUHH)
------------------------------------------------------------------------- */

#ifndef LMP_ELECTRODE_MATRIX_H
#define LMP_ELECTRODE_MATRIX_H

#include "pointers.h"

#include <map>
#include <unordered_map>

namespace LAMMPS_NS {

class ElectrodeMatrix : protected Pointers {
 public:
  ElectrodeMatrix(class LAMMPS *, int, double);
  void setup(const std::unordered_map<tagint, int> &, class Pair *, class NeighList *);
  void setup_tf(const std::map<int, double> &);
  void setup_eta(int);
  void compute_array(double **, bool);
  int igroup;

 private:
  int groupbit;
  bigint ngroup;
  double **cutsq;
  double g_ewald, eta;
  bool tfflag;
  bool etaflag;
  int eta_index;
  std::map<int, double> tf_types;
  std::unordered_map<tagint, int> tag_to_iele;
  bool assigned;
  std::vector<bigint> mpos;
  class Pair *pair;
  class NeighList *list;
  class ElectrodeKSpace *electrode_kspace;

  void update_mpos();
  void pair_contribution(double **);
  void self_contribution(double **);
  void tf_contribution(double **);
};

}    // namespace LAMMPS_NS

#endif
