// clang-format off
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
   This class contains a series of tools for DEM contacts
   Multiple models can be defined and used to calculate forces
   and torques based on contact geometry

   Contributing authors:
   Dan Bolintineanu (SNL), Joel Clemmer (SNL)
----------------------------------------------------------------------- */

#include "contact.h"
#include "contact_sub_models.h"
#include "contact_normal_models.h"
#include "contact_tangential_models.h"
#include "contact_damping_models.h"
#include "contact_rolling_models.h"
#include "contact_twisting_models.h"
#include "contact_heat_models.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace Contact;
using namespace MathExtra;

ContactModel::ContactModel(LAMMPS *lmp) : Pointers(lmp)
{
  limit_damping = 0;
  beyond_contact = 0;
  nondefault_history_transfer = 0;

  contact_type = PAIR;

  reset_contact();

  normal_model = nullptr;
  damping_model = nullptr;
  tangential_model = nullptr;
  rolling_model = nullptr;
  twisting_model = nullptr;
  heat_model = nullptr;

  for (int i = 0; i < NSUBMODELS; i++) sub_models[i] = nullptr;
  transfer_history_factor = nullptr;
}

/* ---------------------------------------------------------------------- */

ContactModel::~ContactModel()
{
  delete [] transfer_history_factor;

  delete normal_model;
  delete damping_model;
  delete tangential_model;
  delete rolling_model;
  delete twisting_model;
  delete heat_model;
}

/* ---------------------------------------------------------------------- */

void ContactModel::init_model(std::string model_name, ModelType model_type)
{
  if (model_type == NORMAL) {
    delete normal_model;
    if (model_name == "none") normal_model = new NormalNone(lmp);
    else if (model_name == "hooke") normal_model = new NormalHooke(lmp);
    else if (model_name == "hertz") normal_model = new NormalHertz(lmp);
    else if (model_name == "hertz/material") normal_model = new NormalHertzMaterial(lmp);
    else if (model_name == "dmt") normal_model = new NormalDMT(lmp);
    else if (model_name == "jkr") normal_model = new NormalJKR(lmp);
    else error->all(FLERR, "Normal model name {} not recognized", model_name);
    sub_models[model_type] = normal_model;

  } else if (model_type == TANGENTIAL) {
    delete tangential_model;
    if (model_name == "none") tangential_model = new TangentialNone(lmp);
    else if (model_name == "linear_nohistory") tangential_model = new TangentialLinearNoHistory(lmp);
    else if (model_name == "linear_history") tangential_model = new TangentialLinearHistory(lmp);
    else if (model_name == "mindlin") tangential_model = new TangentialMindlin(lmp);
    else if (model_name == "mindlin/force") tangential_model = new TangentialMindlinForce(lmp);
    else if (model_name == "mindlin_rescale") tangential_model = new TangentialMindlinRescale(lmp);
    else if (model_name == "mindlin_rescale/force") tangential_model = new TangentialMindlinRescaleForce(lmp);
    else error->all(FLERR, "Tangential model name {} not recognized", model_name);
    sub_models[model_type] = tangential_model;

  } else if (model_type == DAMPING) {
    delete damping_model;
    if (model_name == "none") damping_model = new DampingNone(lmp);
    else if (model_name == "velocity") damping_model = new DampingVelocity(lmp);
    else if (model_name == "mass_velocity") damping_model = new DampingMassVelocity(lmp);
    else if (model_name == "viscoelastic") damping_model = new DampingViscoelastic(lmp);
    else if (model_name == "tsuji") damping_model = new DampingTsuji(lmp);
    else error->all(FLERR, "Damping model name {} not recognized", model_name);
    sub_models[model_type] = damping_model;

  } else if (model_type == ROLLING) {
    delete rolling_model;
    rolling_defined = 1;
    if (model_name == "none") {
      rolling_model = new RollingNone(lmp);
      rolling_defined = 0;
    } else if (model_name == "sds") rolling_model = new RollingSDS(lmp);
    else error->all(FLERR, "Rolling model name {} not recognized", model_name);
    sub_models[model_type] = rolling_model;

  } else if (model_type == TWISTING) {
    delete twisting_model;
    twisting_defined = 1;
    if (model_name == "none") {
      twisting_model = new TwistingNone(lmp);
      twisting_defined = 0;
    } else if (model_name == "sds") twisting_model = new TwistingSDS(lmp);
    else if (model_name == "marshall") twisting_model = new TwistingMarshall(lmp);
    else error->all(FLERR, "Twisting model name {} not recognized", model_name);
    sub_models[model_type] = twisting_model;

  } else if (model_type == HEAT) {
    delete heat_model;
    heat_defined = 1;
    if (model_name == "none") {
      heat_model = new HeatNone(lmp);
      heat_defined = 0;
    } else if (model_name == "area") heat_model = new HeatArea(lmp);
    else error->all(FLERR, "Heat model name not {} recognized", model_name);
    sub_models[model_type] = heat_model;
  } else {
    error->all(FLERR, "Illegal model type {}", model_type);
  }

  sub_models[model_type]->name.assign(model_name);
  sub_models[model_type]->contact = this;
  sub_models[model_type]->allocate_coeffs();
}

/* ---------------------------------------------------------------------- */

int ContactModel::init_classic_model(char **arg, int iarg, int narg)
{
  double kn, kt, gamman, gammat, xmu;

  if (iarg + 6 >= narg)
    error->all(FLERR,"Insufficient arguments provided for classic gran model command");

  kn = utils::numeric(FLERR,arg[iarg + 1],false,lmp);
  if (strcmp(arg[iarg + 2],"NULL") == 0) kt = kn * 2.0 / 7.0;
  else kt = utils::numeric(FLERR,arg[iarg + 2],false,lmp);

  gamman = utils::numeric(FLERR,arg[iarg + 3],false,lmp);
  if (strcmp(arg[iarg + 4],"NULL") == 0) gammat = 0.5 * gamman;
  else gammat = utils::numeric(FLERR,arg[iarg + 4],false,lmp);

  xmu = utils::numeric(FLERR,arg[iarg + 5],false,lmp);
  int dampflag = utils::inumeric(FLERR,arg[iarg + 6],false,lmp);
  if (dampflag == 0) gammat = 0.0;

  if (kn < 0.0 || kt < 0.0 || gamman < 0.0 || gammat < 0.0 ||
      xmu < 0.0 || xmu > 10000.0 || dampflag < 0 || dampflag > 1)
    error->all(FLERR,"Illegal classic gran model command");

  if (strcmp(arg[iarg],"hooke") == 0) {
    init_model("hooke", NORMAL);
    init_model("linear_nohistory", TANGENTIAL);
    init_model("mass_velocity", DAMPING);
  } else if (strcmp(arg[iarg],"hooke/history") == 0) {
    init_model("hooke", NORMAL);
    init_model("linear_history", TANGENTIAL);
    init_model("mass_velocity", DAMPING);
  } else if (strcmp(arg[iarg],"hertz/history") == 0) {
    // convert Kn and Kt from pressure units to force/distance^2 if Hertzian
    kn /= force->nktv2p;
    kt /= force->nktv2p;
    init_model("hertz", NORMAL);
    init_model("mindlin", TANGENTIAL);
    init_model("viscoelastic", DAMPING);
  } else error->all(FLERR,"Invalid classic gran model");

  // ensure additional models are undefined
  init_model("none", ROLLING);
  init_model("none", TWISTING);
  init_model("none", HEAT);

  // manually parse coeffs
  normal_model->coeffs[0] = kn;
  normal_model->coeffs[1] = gamman;
  tangential_model->coeffs[0] = kt;
  tangential_model->coeffs[1] = gammat;
  tangential_model->coeffs[2] = xmu;

  normal_model->coeffs_to_local();
  tangential_model->coeffs_to_local();
  damping_model->coeffs_to_local();

  iarg += 7;
  return iarg;
}

/* ---------------------------------------------------------------------- */

void ContactModel::init()
{
  int i, j;
  for (i = 0; i < NSUBMODELS; i++)
    if (!sub_models[i]) init_model("none", (ModelType) i);

  // Must have valid normal, damping, and tangential models
  if (normal_model->name == "none") error->all(FLERR, "Must specify normal contact model");
  if (damping_model->name == "none") error->all(FLERR, "Must specify damping contact model");
  if (tangential_model->name == "none") error->all(FLERR, "Must specify tangential contact model");

  int size_cumulative;
  size_history = 0;
  for (i = 0; i < NSUBMODELS; i++) {
    if (sub_models[i]->nondefault_history_transfer)
      nondefault_history_transfer = 1;
    if (sub_models[i]->beyond_contact)
      beyond_contact = 1;
    size_history += sub_models[i]->size_history;
  }

  if (nondefault_history_transfer) {
    transfer_history_factor = new double[size_history];

    for (i = 0; i < size_history; i++) {
      // Find which model owns this history value
      size_cumulative = 0;
      for (j = 0; j < NSUBMODELS; j++) {
        if (size_cumulative + sub_models[j]->size_history > i) break;
        size_cumulative += sub_models[j]->size_history;
      }

      // Check if model has nondefault transfers, if so copy its array
      transfer_history_factor[i] = -1;
      if (j != NSUBMODELS) {
        if (sub_models[j]->nondefault_history_transfer) {
          transfer_history_factor[i] = sub_models[j]->transfer_history_factor[i - size_cumulative];
        }
      }
    }
  }

  for (i = 0; i < NSUBMODELS; i++) sub_models[i]->init();
}

/* ---------------------------------------------------------------------- */

int ContactModel::mix_coeffs(ContactModel *c1, ContactModel *c2)
{
  int i;
  for (i = 0; i < NSUBMODELS; i++) {
    if (c1->sub_models[i]->name != c2->sub_models[i]->name) return i;

    init_model(c1->sub_models[i]->name, (ModelType) i);
    sub_models[i]->mix_coeffs(c1->sub_models[i]->coeffs, c2->sub_models[i]->coeffs);
  }

  limit_damping = MAX(c1->limit_damping, c2->limit_damping);

  return -1;
}

/* ---------------------------------------------------------------------- */

void ContactModel::write_restart(FILE *fp)
{
  int num_char, num_coeffs;

  for (int i = 0; i < NSUBMODELS; i++) {
    num_char = sub_models[i]->name.length();
    num_coeffs = sub_models[i]->num_coeffs;
    fwrite(&num_char, sizeof(int), 1, fp);
    fwrite(sub_models[i]->name.data(), sizeof(char), num_char, fp);
    fwrite(&num_coeffs, sizeof(int), 1, fp);
    fwrite(sub_models[i]->coeffs, sizeof(double), num_coeffs, fp);
  }
}

/* ---------------------------------------------------------------------- */

void ContactModel::read_restart(FILE *fp)
{
  int num_char, num_coeff;

  for (int i = 0; i < NSUBMODELS; i++) {
    if (comm->me == 0)
      utils::sfread(FLERR, &num_char, sizeof(int), 1, fp, nullptr, error);
    MPI_Bcast(&num_char, 1, MPI_INT, 0, world);

    std::string model_name (num_char, ' ');
    if (comm->me == 0)
      utils::sfread(FLERR, const_cast<char*>(model_name.data()), sizeof(char),num_char, fp, nullptr, error);
    MPI_Bcast(const_cast<char*>(model_name.data()), num_char, MPI_CHAR, 0, world);

    init_model(model_name, (ModelType) i);

    if (comm->me == 0) {
      utils::sfread(FLERR, &num_coeff, sizeof(int), 1, fp, nullptr, error);
      if (num_coeff != sub_models[i]->num_coeffs)
        error->one(FLERR, "Invalid contact model written to restart file");
    }
    MPI_Bcast(&num_coeff, 1, MPI_INT, 0, world);

    if (comm->me == 0) {
      utils::sfread(FLERR, sub_models[i]->coeffs, sizeof(double), num_coeff, fp, nullptr, error);
    }
    MPI_Bcast(sub_models[i]->coeffs, num_coeff, MPI_DOUBLE, 0, world);
    sub_models[i]->coeffs_to_local();
  }
}

/* ---------------------------------------------------------------------- */

bool ContactModel::check_contact()
{
  if (contact_type == WALL) {
    // Used by fix_wall_gran.cpp
    //   radj = radius of wall
    //   dx already provided
    rsq = lensq3(dx);
    radsum = radi;
    if (radj == 0) Reff = radi;
    else Reff = radi * radj / (radi + radj);
  } else if (contact_type == WALLREGION) {
    // Used by fix_wall_gran_region.cpp
    //   radj = radius of wall
    //   dx and r already provided
    rsq = r * r;
    radsum = radi;
    if (radj == 0) Reff = radi;
    else Reff = radi * radj / (radi + radj);
  } else {
    sub3(xi, xj, dx);
    rsq = lensq3(dx);
    radsum = radi + radj;
    Reff = radi * radj / radsum;
  }

  touch = normal_model->touch();
  return touch;
}

/* ---------------------------------------------------------------------- */

void ContactModel::prep_contact()
{
  double temp[3];

  // Standard geometric quantities
  r = sqrt(rsq);
  rinv = 1.0 / r;
  delta = radsum - r;
  dR = delta * Reff;
  scale3(rinv, dx, nx);

  // relative translational velocity
  sub3(vi, vj, vr);

  // normal component
  vnnr = dot3(vr, nx); //v_R . n
  scale3(vnnr, nx, vn);

  // tangential component
  sub3(vr, vn, vt);

  // relative rotational velocity
  scaleadd3(radi, omegai, radj, omegaj, wr);

  // relative tangential velocities
  cross3(wr, nx, temp);
  sub3(vt, temp, vtr);
  vrel = len3(vtr);

  if (rolling_defined || twisting_defined)
    sub3(omegai, omegaj, relrot);

  if (rolling_defined) {
    // rolling velocity, see eq. 31 of Wang et al, Particuology v 23, p 49 (2015)
    // this is different from the Marshall papers, which use the Bagi/Kuhn formulation
    // for rolling velocity (see Wang et al for why the latter is wrong)
    vrl[0] = Reff * (relrot[1] * nx[2] - relrot[2] * nx[1]);
    vrl[1] = Reff * (relrot[2] * nx[0] - relrot[0] * nx[2]);
    vrl[2] = Reff * (relrot[0] * nx[1] - relrot[1] * nx[0]);
  }

  if (twisting_defined) {
    // omega_T (eq 29 of Marshall)
    magtwist = dot3(relrot, nx);
  }
}

/* ---------------------------------------------------------------------- */

void ContactModel::calculate_forces()
{
  // calculate forces/torques

  forces[0] = 0.0;
  double Fne, Fdamp;
  area = normal_model->calculate_area();
  normal_model->set_knfac();
  Fne = normal_model->calculate_forces();

  Fdamp = damping_model->calculate_forces();
  Fntot = Fne + Fdamp;

  normal_model->set_fncrit(); // Needed for tangential, rolling, twisting

  if (limit_damping && Fntot < 0.0) Fntot = 0.0;

  tangential_model->calculate_forces();
  if (rolling_defined) rolling_model->calculate_forces();
  if (twisting_defined) twisting_model->calculate_forces();

  // sum contributions

  scale3(Fntot, nx, forces);
  add3(forces, fs, forces);

  //May need to rethink eventually tris..
  cross3(nx, fs, torquesi);

  double dist_to_contact = radi - 0.5 * delta;
  scale3(dist_to_contact, torquesi);

  if (contact_type == PAIR) {
    copy3(torquesi, torquesj);
    dist_to_contact = radj - 0.5 * delta;
    scale3(dist_to_contact, torquesj);
  }

  double torroll[3];
  if (rolling_defined) {
    cross3(nx, fr, torroll);
    scale3(Reff, torroll);
    add3(torquesi, torroll, torquesi);
    if (contact_type == PAIR) sub3(torquesj, torroll, torquesj);
  }

  double tortwist[3];
  if (twisting_defined) {
    scale3(magtortwist, nx, tortwist);
    add3(torquesi, tortwist, torquesi);
    if (contact_type == PAIR) sub3(torquesj, tortwist, torquesj);
  }

  if (heat_defined) {
    dq = heat_model->calculate_heat();
  }
}

/* ----------------------------------------------------------------------
   compute pull-off distance (beyond contact) for a given radius and atom type
   use temporary variables since this does not use a specific contact geometry
------------------------------------------------------------------------- */

double ContactModel::pulloff_distance(double radi, double radj)
{
  return normal_model->pulloff_distance(radi, radj);
}
