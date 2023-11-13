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

#ifndef LMP_ELECTRODE_VECTOR_H
#define LMP_ELECTRODE_VECTOR_H

#include "pointers.h"
#include <map>

namespace LAMMPS_NS {

class ElectrodeVector : protected Pointers {
 public:
  ElectrodeVector(class LAMMPS *, int, int, double, bool);
  ~ElectrodeVector() override;
  void setup(class Pair *, class NeighList *, bool);
  void setup_tf(const std::map<int, double> &);
  void setup_eta(int);
  void compute_vector(double *);
  int igroup, source_group;

 private:
  bool invert_source;
  int groupbit, source_grpbit;
  bigint ngroup;
  double **cutsq;
  double g_ewald, eta;
  bool tfflag;
  bool etaflag;
  int eta_index;
  std::map<int, double> tf_types;
  class Pair *pair;
  class NeighList *list;
  class ElectrodeKSpace *electrode_kspace;

  void pair_contribution(double *);
  void self_contribution(double *);
  void tf_contribution(double *);

  double kspace_time_total;
  double pair_time_total;
  double boundary_time_total;
  double b_time_total;

  bool timer_flag;
};

}    // namespace LAMMPS_NS

#endif
