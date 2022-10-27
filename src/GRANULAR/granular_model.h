/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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

#define EPSILON 1e-10
#define NSUBMODELS 6

enum SubmodelType {
  NORMAL = 0,
  DAMPING = 1,
  TANGENTIAL = 2,
  ROLLING = 3,
  TWISTING = 4,
  HEAT = 5
}; // Relative order matters since some derive coeffs from others

enum ContactType {
  PAIR = 0,
  WALL = 1,
  WALLREGION = 2
};

// forward declaration
class GSM;
class GSMNormal;
class GSMDamping;
class GSMTangential;
class GSMRolling;
class GSMTwisting;
class GSMHeat;

class GranularModel : protected Pointers {
 public:
  GranularModel(class LAMMPS *);
  ~GranularModel();
  void init();
  bool check_contact();
  void calculate_forces();
  double pulloff_distance(double, double);

  int init_submodel(char **, int, int, SubmodelType);
  int init_classic_model(char **, int, int);
  void construct_submodel(std::string, SubmodelType);
  int mix_coeffs(GranularModel*, GranularModel*);

  void write_restart(FILE *);
  void read_restart(FILE *);

  // Sub models
  GSMNormal *normal_model;
  GSMDamping *damping_model;
  GSMTangential *tangential_model;
  GSMRolling *rolling_model;
  GSMTwisting *twisting_model;
  GSMHeat *heat_model;
  GSM *sub_models[NSUBMODELS];  // Need to resize if we add more model flavors

  // Extra options
  int beyond_contact, limit_damping, history_update;
  ContactType contact_type;

  // History variables
  int size_history, nondefault_history_transfer;
  double *transfer_history_factor;
  double *history;

  // Contact properties/output
  double forces[3], torquesi[3], torquesj[3], dq;

  double radi, radj, meff, dt, Ti, Tj, area;
  double Fntot, magtortwist;

  int i, j;
  double *xi, *xj, *vi, *vj, *omegai, *omegaj;
  double fs[3], fr[3], ft[3];

  double dx[3], nx[3], r, rsq, rinv, Reff, radsum, delta, dR;
  double vr[3], vn[3], vnnr, vt[3], wr[3], vtr[3], vrl[3], relrot[3], vrel;
  double magtwist;
  bool touch;

 protected:
  int rolling_defined, twisting_defined, heat_defined; // Used to quickly skip undefined submodels
  int classic_model;
};

}    // namespace Granular_NS
}    // namespace LAMMPS_NS

#endif
