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

#ifndef LMP_GRAN_SUB_MOD_H
#define LMP_GRAN_SUB_MOD_H

#include "granular_model.h"
#include "pointers.h"    // IWYU pragma: export

namespace LAMMPS_NS {
namespace Granular_NS {

class GranSubMod : protected Pointers {
 public:
  GranSubMod(class GranularModel *, class LAMMPS *);
  virtual ~GranSubMod();

  int num_coeffs;
  double *coeffs;
  void read_restart();
  virtual void mix_coeffs(double*, double*);
  virtual void coeffs_to_local() {};
  virtual void init() {}; // called after all sub models + coeffs defined

  void allocate_coeffs();
  std::string name;

  int size_history;
  int nondefault_history_transfer;
  double *transfer_history_factor;

  int history_index;
  int beyond_contact;  // If the sub model contact extends beyond overlap
  int allow_cohesion;  // If the sub model works with a cohesive normal force
  int area_flag;       // If the sub model requires area

  GranularModel *gm;

 protected:
  int allocated;

  double mix_stiffnessE(double, double, double, double);
  double mix_stiffnessG(double, double, double, double);
  double mix_stiffnessE_wall(double, double);
  double mix_stiffnessG_wall(double, double);
  double mix_geom(double, double);
};

}    // namespace GranularModel
}    // namespace LAMMPS_NS

#endif /* GRAN_SUB_MOD_H */
