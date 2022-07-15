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

#include "pointers.h"
#include "normal_contact_models.h"
#include "tangential_contact_models.h"
#include "damping_contact_models.h"
#include "rolling_contact_models.h"
#include "twisting_contact_models.h"
#include "heat_models.h"

using namespace LAMMPS_NS;

namespace Contact {

#define EPSILON 1e-10

class ContactModel : protected Pointers {
public:
  ContactModel();
  ~ContactModel();
  int init();
  bool check_contact();
  void reset_contact();
  void calculate_forces();
  double calculate_heat();
  double pulloff_distance(double, double);

  void init_model(char*, int);

  void mix_coeffs(ContactModel*, ContactModel*);

  void write_restart(FILE *);
  void read_restart(FILE *);

  // Sub models
  NormalModel normal_model;
  DampingModel damping_model;
  TangentialModel tangential_model;
  RollingModel rolling_model;
  TwistingModel twisting_model;
  HeatModel heat_model;
  SubModel *sub_models[6];  // Need to resize if we add more model flavors

  // Extra options
  int beyond_contact, limit_damping;
  double cutoff_type;

  // Contact properties/output
  double *forces, *torquesi, *torquesj;
  double *history;

  double radi, radj, meff, dt, Ti, Tj, area;
  double Fntot, magtortwist;

  double *xi, *xj, *vi, *vj, *omegai, *omegaj;
  double fs[3], fr[3], ft[3];

private:
  double dx[3], nx[3], r, rsq, rinv, Reff, radsum, delta, dR;
  double vr[3], vn[3], vnnr, vt[3], wr[3], vtr[3], vrl[3], relrot[3], vrel;
  double magtwist;
  bool touch;

  int prep_flag, check_flag;
};

}    // namespace Contact

#endif
