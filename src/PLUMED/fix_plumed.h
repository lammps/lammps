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

#ifdef FIX_CLASS
// clang-format off
FixStyle(plumed,FixPlumed);
// clang-format on
#else

#ifndef LMP_FIX_PLUMED_H
#define LMP_FIX_PLUMED_H

#include "fix.h"

#include <vector>

// forward declaration
namespace PLMD {
class Plumed;
}

namespace LAMMPS_NS {

class FixPlumed : public Fix {
 public:
  FixPlumed(class LAMMPS *, int, char **);
  ~FixPlumed() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double compute_scalar() override;
  void *extract(const char *, int &) override;
  void reset_dt() override;
  int modify_param(int narg, char **arg) override;
  double memory_usage() override;

 private:
  enum {
    PATH_INTEGRAL_OFF,
    PATH_INTEGRAL_CENTROID,
    PATH_INTEGRAL_BEAD_MEAN,
    PATH_INTEGRAL_BEAD_DENSITY
  };

  int path_contraction_flag;
  double path_contraction;
  std::vector<double> contracted_positions, contracted_forces, contraction_mean;

  PLMD::Plumed *p;                    // pointer to plumed object
  class Fix *pimd_fix;                // fix providing the Cartesian PIMD centroid
  int nlocal;                         // number of atoms local to this process
  int natoms;                         // total number of atoms
  int path_integral_mode;             // path-integral coupling mode
  int plumed_active;                  // this partition runs a PLUMED instance
  double centroid_force_scale;        // chain-rule force scale supplied by the PIMD fix
  double bead_density_force_scale;    // 1/P scaling for shared bead-density forces
  int *gatindex;                      // array of atom indexes local to this process
  double *masses;                     // array of masses for local atoms
  double *charges;                    // array of charges for local atoms
  double *centroid_coordinates;       // tag-ordered Cartesian centroid coordinates
  double *centroid_positions;         // local centroid coordinates passed to PLUMED
  double *centroid_forces;            // local centroid forces returned by PLUMED
  double *centroid_forces_all;        // tag-ordered centroid forces shared by all beads
  int *centroid_virial_pending;       // asks the NMPIMD fix to collect this step's bias virial
  int *bead_bias_virial_pending;      // asks NMPIMD to transform bead forces and refresh virial
  double *forces_before_plumed;       // local forces before bead-density bias
  int nlevels_respa;                  // this is something to enable respa
  double bias;                        // output bias potential
  class Compute *c_pe;                // Compute for the energy
  class Compute *c_press;             // Compute for the pressure
  int plumedNeedsEnergy;              // Flag to trigger calculation of the
                                      // energy and virial
  char *id_pe, *id_press;             // ID for potential energy and pressure compute
  char *nonfinite_trace_prefix;       // opt-in diagnostic output prefix
  char *id_pimd;                      // ID for the coupled fix pimd/langevin

  void check_path_integral_compatibility();
  void check_normal_mode_post_force_order();
  void update_atom_data();
  void post_force_centroid();
  void prepare_contracted_coordinates();
  void pullback_contracted_forces();
  void trace_nonfinite_state(const char *, const double *, bool);
  void trace_centroid_nonfinite(const char *);
};

};    // namespace LAMMPS_NS

#endif
#endif
