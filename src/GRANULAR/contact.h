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

#ifndef LMP_CONTACT_H
#define LMP_CONTACT_H

#include "pointers.h"    // IWYU pragma: export

namespace LAMMPS_NS {
namespace Contact {

#define EPSILON 1e-10
#define NSUBMODELS 6

enum ModelType {
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
class NormalModel;
class DampingModel;
class TangentialModel;
class RollingModel;
class TwistingModel;
class HeatModel;
class SubModel;

class ContactModel : protected Pointers {
 public:
  ContactModel(class LAMMPS *);
  ~ContactModel();
  void init();
  bool check_contact();
  void prep_contact();
  void calculate_forces();
  double pulloff_distance(double, double);

  void init_model(std::string, ModelType);
  int init_classic_model(char **, int, int);
  int mix_coeffs(ContactModel*, ContactModel*);

  void write_restart(FILE *);
  void read_restart(FILE *);

  // Sub models
  NormalModel *normal_model;
  DampingModel *damping_model;
  TangentialModel *tangential_model;
  RollingModel *rolling_model;
  TwistingModel *twisting_model;
  HeatModel *heat_model;
  SubModel *sub_models[NSUBMODELS];  // Need to resize if we add more model flavors

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

}    // namespace Contact
}    // namespace LAMMPS_NS

#endif
