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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meissner (TUHH)
------------------------------------------------------------------------- */

#ifndef LMP_ELECTRODE_MAT_CG_H
#define LMP_ELECTRODE_MAT_CG_H

#include "electrode_cg.h"
#include <unordered_map>

namespace LAMMPS_NS {

class ElectrodeMatCG : public ElectrodeCG {
 public:
  // ChargeSolver methods
  ElectrodeMatCG(class LAMMPS *);
  ~ElectrodeMatCG() noexcept;    // TODO why do we need noexcept here
  void update_solver(std::vector<tagint>, std::vector<int>) override;
  double memory_use() override;

  //setup
  void setup_solver(double, std::unordered_map<tagint, int>, int);
  void set_elastance(int, double **);

 private:
  int n_mat;
  bool matrix_set;
  double **elastance;
  std::vector<double> qele_world;
  std::unordered_map<tagint, int> tag_to_iele;    // inverse of global taglist:
  std::vector<int> iele_local;                    // electrode IDs owned by me

  std::vector<double> ele_ele_interaction(const std::vector<double> &) override;
};

}    // namespace LAMMPS_NS

#endif

