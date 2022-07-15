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
-------------------------------------------------------------------------

   This class contains a series of tools for DEM contacts
   Multiple models can be defined and used to calculate forces
   and torques based on contact geometry
*/

#include "contact.h"
#include "math_const.h"
#include "math_extra.h"
#include "pointers.h"
#include "comm.h"
#include "error.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathExtra;
using namespace MathConst;

namespace Contact {

enum {NORMAL, TANGENTIAL, DAMPING, ROLLING, TWISTING, HEAT}

ContactModel::ContactModel()
{
  limit_damping = 0;
  beyond_contact = 0;
  cutoff_type = 0.0;
  normal_model = nullptr;
  tangential_model = nullptr;
  damping_model = nullptr;
  rolling_model = nullptr;
  twisting_model = nullptr;
  reset_contact();

  sub_models = {nullptr};
}

/* ---------------------------------------------------------------------- */

void ContactModel::init_model(char *model_name, int model_type)
{
  if (model_type == NORMAL) {
    if (strcmp(model_name, "hooke") == 0) normal_model = new NormalHooke();
    else if (strcmp(model_name, "hertz") == 0) normal_model = new NormalHertz();
    else if (strcmp(model_name, "hertz/material") == 0) normal_model = new NormalHertzMaterial();
    else if (strcmp(model_name, "dmt") == 0) normal_model = new NormalDMT();
    else if (strcmp(model_name, "jkr") == 0) normal_model = new NormalJKR();
    else error->all(FLERR, "Normal model name not recognized");
    sub_models[model_type] = &normal_model;

  } else if (model_type == TANGENTIAL) {
    if (strcmp(model_name, "linear_nohistory") == 0) tangential_model = new TangentialLinearNoHistory();
    else if (strcmp(model_name, "linear_history") == 0) tangential_model = new TangentialLinearHistory();
    else if (strcmp(model_name, "mindlin") == 0) tangential_model = new TangentialMindlin();
    else if (strcmp(model_name, "mindlin/force") == 0) tangential_model = new TangentialMindlinForce();
    else if (strcmp(model_name, "mindlin_rescale") == 0) tangential_model = new TangentialMindlinRescale();
    else if (strcmp(model_name, "mindlin_rescale/force") == 0) tangential_model = new TangentialMindlinRescaleForce();
    else error->all(FLERR, "Tangential model name not recognized");
    sub_models[model_type] = &tangential_model;

  } else if (model_type == DAMPING) {
    if (strcmp(model_name, "velocity") == 0) damping_model = new DampingVelocity();
    else if (strcmp(model_name, "mass_velocity") == 0) damping_model = new DampingMassVelocity();
    else if (strcmp(model_name, "viscoelastic") == 0) damping_model = new DampingViscoelastic();
    else if (strcmp(model_name, "tsuji") == 0) damping_model = new DampingTsuji();
    else error->all(FLERR, "Damping model name not recognized");
    sub_models[model_type] = &damping_model;

  } else if (model_type == ROLLING) {
    if (strcmp(model_name, "none") == 0) delete rolling_model;
    else if (strcmp(model_name, "sds") == 0) rolling_model = new RollingSDS();
    else error->all(FLERR, "Rolling model name not recognized");
    sub_models[model_type] = &rolling_model;

  } else if (model_type == TWISTING) {
    if (strcmp(model_name, "none") == 0) delete twisting_model;
    else if (strcmp(model_name, "sds") == 0) twisting_model = new TwistingSDS();
    else if (strcmp(model_name, "marshall") == 0) twisting_model = new TwistingMarshall();
    else error->all(FLERR, "Twisting model name not recognized");
    sub_models[model_type] = &twisting_model;

  } else if (model_type == HEAT) {
    if (strcmp(model_name, "none") == 0) delete heat_model;
    else if (strcmp(model_name, "area") == 0) heat_model = new HeatArea();
    else error->all(FLERR, "Heat model name not recognized");
    sub_models[model_type] = &heat_model;
  }

  sub_models[model_type]->name.assign(model_name);
  sub_models[model_type]->contact = *this;
  sub_models[model_type]->allocate_coeffs();
}

/* ---------------------------------------------------------------------- */

void ContactModel::mix_coeffs(ContactModel *c1, ContactModel *c2)
{
  for (int i = 0; i < 5; i++)
    sum_model[i]->mix_coeffs(c1->submodel[i], c2->submodel[i]);

  limit_damping = MAX(c1->limit_damping, c2->limit_damping);
  cutoff_type = MAX(c1->cutoff_type, c2->cutoff_type);
}

/* ---------------------------------------------------------------------- */

void ContactModel::write_restart(FILE *fp)
{
  int num_char = -1;

  for (int i = 0; i < 5; i++) {
    if (sub_models[i]) {
      sub_models[i]->write_restart(fp);
    } else {
      fwrite(&num_char, sizeof(int), 1, fp);
    }
  }
}

/* ---------------------------------------------------------------------- */

void ContactModel::read_restart(FILE *fp)
{
  int num_char;

  for (int i = 0; i < 5; i++) {
    if (comm->me == 0)
      utils::sfread(FLERR,&num_char,sizeof(int),1,fp,nullptr,error);
    MPI_BCast(&num_char, 1, MPI_INT, 0, world);

    if (num_char != -1) {
      std::string model_name(num_char);
      if (comm->me == 0)
        utils::sfread(FLERR,const_cast<char*>(model_name.data()),sizeof(char),num_char,fp,nullptr,  error);
      MPI_BCast(const_cast<char*>(model_name.data()), num_char, MPI_CHAR, world);

      init_model(const_cast<char*>(model_name.data(), i));
      sub_models[i]->read_restart();
    }
  }
}

/* ---------------------------------------------------------------------- */

void ContactModel::reset_contact()
{
  radi = radj = Fntot = Fncrit = magtortwist = 0.0;
  xi = xj = vi = vj = omegai = omegaj = nullptr;
  forces = torquesi = torquesj = history = nullptr;

  prep_flag = check_flag = 0;
  touch = false;
}

/* ---------------------------------------------------------------------- */

bool ContactModel::check_contact()
{
  check_flag = 1;

  sub3(xi, xj, dx);
  rsq = lensq3(dx);
  radsum = radi + radj;
  Reff = radi * radj / radsum;

  touch = false;
  if (normal_model.beyond_contact) normal_model.touch(touch);
  else touch = (rsq < radsum*radsum);

  return touch;
}

/* ---------------------------------------------------------------------- */

void ContactModel::prep_contact()
{
  prep_flag = 1;
  // If it hasn't already been done, test if the contact exists
  if (check_flag != 1) touch = check_contact();
  if (!touch) return;

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

  if (roll_model || twist_model)
    sub3(omegai, omegaj, relrot);

  if (rolling_model) {
    // rolling velocity, see eq. 31 of Wang et al, Particuology v 23, p 49 (2015)
    // this is different from the Marshall papers, which use the Bagi/Kuhn formulation
    // for rolling velocity (see Wang et al for why the latter is wrong)
    vrl[0] = Reff * (relrot[1] * nx[2] - relrot[2] * nx[1]);
    vrl[1] = Reff * (relrot[2] * nx[0] - relrot[0] * nx[2]);
    vrl[2] = Reff * (relrot[0] * nx[1] - relrot[1] * nx[0]);
  }

  if (twisting_model) {
    // omega_T (eq 29 of Marshall)
    magtwist = dot3(relrot, nx);
  }
}

/* ---------------------------------------------------------------------- */

void ContactModel::calculate_forces()
{
  // If it hasn't already been done, run prep calculations
  if (prep_flag != 1) prep_contact();
  if (!touch) {
    forces[0] = forces[1] = forces[2] = 0.0;
    return;
  }

  //**********************************************
  // normal forces
  //**********************************************

  double Fne, Fdamp;
  Fne = normal_model.calculate_forces();
  Fdamp = damping_model.calculate_forces();

  Fntot = Fne + Fdamp;
  if (limit_damping && Fntot < 0.0) Fntot = 0.0;
  normal_model.set_fncrit();

  tangential_model.calculate_forces();
  if (rolling_model) rolling_model.calculate_forces();
  if (twisting_model) twisting_model.calculate_forces();

  //**********************************************
  // sum contributions
  //**********************************************

  scale3(Fntot, nx, forces);
  add3(forces, fs, forces);

  //May need to rethink this for use with walls (and eventually tris)..
  cross3(nx, fs, torquesi);
  copy3(torquesi, torquesj);

  double dist_to_contact = radi-0.5*delta;
  scale3(dist_to_contact, torquesi);
  dist_to_contact = radj-0.5*delta;
  scale3(dist_to_contact, torquesj);

  double torroll[3];
  if (rolling_model) {
    cross3(nx, fr, torroll);
    scale3(Reff, torroll);
    add3(torquesi, torroll, torquesi);
    sub3(torquesj, torroll, torquesj);
  }

  double tortwist[3];
  if (twisting_model) {
    scale3(magtortwist, nx, tortwist);
    add3(torquesi, tortwist, torquesi);
    sub3(torquesj, tortwist, torquesj);
  }
}

/* ---------------------------------------------------------------------- */

double ContactModel::calculate_heat()
{
  return heat_model.calculate_heat();
}

/* ----------------------------------------------------------------------
   compute pull-off distance (beyond contact) for a given radius and atom type
   use temporary variables since this does not use a specific contact geometry
------------------------------------------------------------------------- */

double ContactModel::pulloff_distance(double radi, double radj)
{
  return normal_model.pulloff_distance(radi, radj);
}

}
