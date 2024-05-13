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

#ifndef LMP_GRANULAR_MODEL_H
#define LMP_GRANULAR_MODEL_H

#include "pointers.h"    // IWYU pragma: export

namespace LAMMPS_NS {
namespace Granular_NS {

enum SubModelType {
  NORMAL = 0,
  DAMPING,
  TANGENTIAL,
  ROLLING,
  TWISTING,
  HEAT,
  NSUBMODELS
};

enum ContactType {
  PAIR = 0,
  WALL = 1,
  WALLREGION = 2
};

// forward declarations
class GranSubMod;
class GranSubModNormal;
class GranSubModDamping;
class GranSubModTangential;
class GranSubModRolling;
class GranSubModTwisting;
class GranSubModHeat;

class GranularModel : protected Pointers {
 public:
  GranularModel(class LAMMPS *);
  ~GranularModel() override;
  void init();
  bool check_contact();
  void calculate_forces();
  double pulloff_distance(double, double);

  int add_sub_model(char **, int, int, SubModelType);
  int define_classic_model(char **, int, int);
  void construct_sub_model(std::string, SubModelType);
  int mix_coeffs(GranularModel*, GranularModel*);

  void write_restart(FILE *);
  void read_restart(FILE *);

  // Sub models
  GranSubModNormal *normal_model;
  GranSubModDamping *damping_model;
  GranSubModTangential *tangential_model;
  GranSubModRolling *rolling_model;
  GranSubModTwisting *twisting_model;
  GranSubModHeat *heat_model;
  GranSubMod *sub_models[NSUBMODELS];

  // Extra options
  int beyond_contact, limit_damping, history_update;
  ContactType contact_type;

  // History variables
  int size_history, nondefault_history_transfer;
  double *transfer_history_factor;
  double *history;

  // Contact properties/output
  double Fnormal, forces[3], torquesi[3], torquesj[3], dq;

  double radi, radj, meff, dt, Ti, Tj, contact_radius;
  double Fntot, magtortwist;

  double *xi, *xj, *vi, *vj, *omegai, *omegaj;
  double fs[3], fr[3], ft[3];

  double dx[3], nx[3], r, rsq, rinv, Reff, radsum, delta, dR;
  double vr[3], vn[3], vnnr, vt[3], wr[3], vtr[3], vrl[3], relrot[3], vrel;
  double magtwist;
  bool touch;

 protected:
  int rolling_defined, twisting_defined, heat_defined; // Flag optional sub models
  int classic_model;                                   // Flag original pair/gran calculations
  int contact_radius_flag;                             // Flag whether contact radius is needed

  int nclass;

  typedef class GranSubMod *(*GranSubModCreator)(class GranularModel *, class LAMMPS *);
  GranSubModCreator *gran_sub_mod_class;
  char **gran_sub_mod_names;
  int *gran_sub_mod_types;
};

}    // namespace Granular_NS
}    // namespace LAMMPS_NS

#endif
