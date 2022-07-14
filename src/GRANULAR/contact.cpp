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

ContactModel::ContactModel()
{
  limit_damping = 0;
  cutoff_type = 0.0;
  normal_model = nullptr;
  tangential_model = nullptr;
  damping_model = nullptr;
  rolling_model = nullptr;
  twisting_model = nullptr;
  reset_contact();
}

/* ---------------------------------------------------------------------- */
void ContactModel::init_normal(char *model_name){
  if (strcmp(model_name, "hooke") == 0) normal_model = new Hooke();
  else if (strcmp(model_name, "hertz") == 0) normal_model = new Hertz(0);
  else if (strcmp(model_name, "hertz/material") == 0) normal_model = new Hertz(1);
  //...
  else error->all(FLERR, "Normal model name not recognized");

  normal_model->model_name.assign(model_name);
  normal_model->contact = *this;
  normal_model->allocate_coeffs();
}

/* ---------------------------------------------------------------------- */
void ContactModel::init_damping(char *model_name){
  if (strcmp(model_name, "velocity") == 0) tangential_model = new Velocity(*this);
  //...
  else error->all(FLERR, "Damping model name not recognized");
  tangential_model->model_name.assign(model_name);
  tangential_model->contact = *this;
  tangential_model->allocate_coeffs();
}

/* ---------------------------------------------------------------------- */
void ContactModel::init_tangential(char *model_name){
  if (strcmp(model_name, "linear") == 0) tangential_model = new LinearNohistory(*this);
  //...
  else error->all(FLERR, "Tangential model name not recognized");
  damping_model->model_name.assign(model_name);
  damping_model->contact = *this;
  damping_model->allocate_coeffs();
}


/* ---------------------------------------------------------------------- */
// .... same for rolling, twisting


void ContactModel::write_restart(FILE *fp){
  normal_model->write_restart(fp);
  tangential_model->write_restart(fp);
  damping_model->write_restart(fp);
  if (rolling_model){
    rolling_model->write_restart(fp);
  }
  else{
    int num_char = -1;
    fwrite(&num_char, sizeof(int), 1, fp);
  }
  if (twisting_model){
    twisting_model->write_restart(fp);
  }
  else{
    int num_char = -1;
    fwrite(&num_char, sizeof(int), 1, fp);
  }
}

void ContactModel::read_restart(FILE *fp){
  int num_char;
  //Normal model
  if (comm->me == 0){
    utils::sfread(FLERR,&num_char,sizeof(int),1,fp,nullptr,error);
  }
  MPI_BCast(&num_char, 1, MPI_INT, 0, world);
  std::string model_name(num_char);
  if (comm->me == 0){
    utils::sfread(FLERR,const_cast<char*>(model_name.data()),sizeof(char),num_char,fp,nullptr,error);
  }
  MPI_BCast(const_cast<char*>(model_name.data()), num_char, MPI_CHAR, world);
  init_normal(const_cast<char*>(model_name.data()));
  normal_model->read_restart();

  //Tangential model
  if (comm->me == 0){
    utils::sfread(FLERR,&num_char,sizeof(int),1,fp,nullptr,error);
  }
  MPI_BCast(&num_char, 1, MPI_INT, 0, world);
  std::string model_name(num_char);
  if (comm->me == 0){
    utils::sfread(FLERR,const_cast<char*>(model_name.data()),sizeof(char),num_char,fp,nullptr,error);
  }
  init_tangential(const_cast<char*>(model_name.data()));
  tangential_model->read_restart();

  //Damping
  if (comm->me == 0){
    utils::sfread(FLERR,&num_char,sizeof(int),1,fp,nullptr,error);
  }
  MPI_BCast(&num_char, 1, MPI_INT, 0, world);
  std::string model_name(num_char);
  if (comm->me == 0){
    utils::sfread(FLERR,const_cast<char*>(model_name.data()),sizeof(char),num_char,fp,nullptr,error);
  }
  init_tangential(const_cast<char*>(model_name.data()));
  damping_model->read_restart();

  //Optional (rolling, twisting) - only if num_char is > 0.

}

/* ---------------------------------------------------------------------- */

void ContactModel::reset_contact()
{
  radi = radj = 0.0;
  xi = xj = vi = vj = omegai = omegaj = nullptr;

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
  Reff = radi*radj/radsum;

  touch = false;
  if (normal_model->beyond_contact) normal_model->touch(touch);
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

  if (roll_model) {
    // rolling velocity, see eq. 31 of Wang et al, Particuology v 23, p 49 (2015)
    // this is different from the Marshall papers, which use the Bagi/Kuhn formulation
    // for rolling velocity (see Wang et al for why the latter is wrong)
    vrl[0] = Reff * (relrot[1] * nx[2] - relrot[2] * nx[1]);
    vrl[1] = Reff * (relrot[2] * nx[0] - relrot[0] * nx[2]);
    vrl[2] = Reff * (relrot[0] * nx[1] - relrot[1] * nx[0]);
  }

  if (twist_model) {
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
  Fne = normal_model->calculate_forces();
  Fdamp = damping_model->calculate_forces();

  Fntot = Fne + Fdamp;
  if (limit_damping && Fntot < 0.0) Fntot = 0.0;
  normal_model->set_fncrit();

  tangential_model->calculate_forces();
  if (roll_model) roll_model->calculate_forces();
  if (twist_model) twist_model->calculate_forces();

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
  if (roll_model) {
    cross3(nx, fr, torroll);
    scale3(Reff, torroll);
    add3(torquesi, torroll, torquesi);
    sub3(torquesj, torroll, torquesj);
  }

  double tortwist[3];
  if (twist_model) {
    scale3(magtortwist, nx, tortwist);
    add3(torquesi, tortwist, torquesi);
    sub3(torquesj, tortwist, torquesj);
  }
}

/* ---------------------------------------------------------------------- */

double ContactModel::calculate_heat()
{
  double dT = Ti - Tj;
  return conductivity * a * dT;
}

/* ---------------------------------------------------------------------- */



/* ---------------------------------------------------------------------- */

void ContactModel::tangential_no_history()
{
  double gamma_scaled = gamma_tang * damp_normal_prefactor;
  double fsmag, Ft;

  // classic pair gran/hooke (no history)
  fsmag = gamma_scaled * vrel;
  if (vrel != 0.0) Ft = MIN(Fscrit,fsmag) / vrel;
  else Ft = 0.0;

  Ft = -Ft;
  scale3(Ft, vtr, fs);
}

/* ---------------------------------------------------------------------- */

void ContactModel::tangential_history(double *history)
{
  double gamma_scaled = gamma_tang * damp_normal_prefactor;
  double k = k_tang;
  int frame_update = 0;
  double magfs, rsht, shrmag, prjmag, temp_dbl, temp_array[3];

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (history_update) {
    rsht = dot3(history, nx);
    frame_update = fabs(rsht) * k > EPSILON * Fscrit;

    if (frame_update) {
      shrmag = len3(history);
      // projection
      scale3(rsht, nx, history);
      // also rescale to preserve magnitude
      prjmag = len3(history);
      if (prjmag > 0) temp_dbl = shrmag / prjmag;
      else temp_dbl = 0;
      scale3(temp_dbl, history);
    }

    // update history
    // tangential force
    // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
    temp_dbl = k * dt;
    scale3(temp_dbl, vtr, temp_array);
    sub3(history, temp_array, history);
  }

  // tangential forces = history + tangential velocity damping
  temp_dbl = -gamma_norm;
  scale3(temp_dbl, vtr, fs);

  // rescale frictional displacements and forces if needed
  magfs = len3(fs);
  if (magfs > Fscrit) {
    shrmag = len3(history);
    if (shrmag != 0.0) {
      temp_dbl = Fscrit / magfs;
      scale3(temp_dbl, fs, history);
      scale3(gamma_norm, vtr, temp_array);
      add3(history, temp_array, history);
      temp_dbl = Fscrit / magfs;
      scale3(temp_dbl, fs);
    } else {
      zero3(fs);
    }
  }
}

/* ---------------------------------------------------------------------- */

void ContactModel::tangential_mindlin(double *history)
{
  double k_scaled, gamma_scaled, magfs, rsht, shrmag, prjmag, temp_dbl;
  double temp_array[3];
  int frame_update = 0;

  gamma_scaled = gamma_tang * damp_normal_prefactor;
  k_scaled = k_tang * a;
  if (mindlin_rescale) {
    // on unloading, rescale the shear displacements/force
    if (a < history[3]) {
      temp_dbl = a / history[3];
      scale3(temp_dbl, history);
    }
  }

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (history_update) {
    rsht = dot3(history, nx);
    if (mindlin_force)
      frame_update = fabs(rsht) > EPSILON * Fscrit;
    else
      frame_update = fabs(rsht) * k_scaled > EPSILON * Fscrit;

    if (frame_update) {
      shrmag = len3(history);
      // projection
      scale3(rsht, nx, history);
      // also rescale to preserve magnitude
      prjmag = len3(history);
      if (prjmag > 0) temp_dbl = shrmag / prjmag;
      else temp_dbl = 0;
      scale3(temp_dbl, history);
    }

    // update history
    if (mindlin_force) {
      // tangential force
      // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
      temp_dbl = -k_scaled * dt;
      scale3(temp_dbl, vtr, temp_array);
    } else {
      scale3(dt, vtr, temp_array);
    }
    add3(history, temp_array, history);

    if (mindlin_rescale) history[3] = a;
  }

  // tangential forces = history + tangential velocity damping
  temp_dbl = -gamma_scaled;
  scale3(temp_dbl, vtr, fs);

  if (! mindlin_force) {
    scale3(k_scaled, history, temp_array);
    add3(fs, temp_array, fs);
  }

  // rescale frictional displacements and forces if needed
  magfs = len3(fs);
  if (magfs > Fscrit) {
    shrmag = len3(history);
    if (shrmag != 0.0) {
      temp_dbl = Fscrit / magfs;
      scale3(temp_dbl, fs, history);
      scale3(gamma_tang, vtr, temp_array);
      add3(history, temp_array, history);
      if (! mindlin_force) {
        temp_dbl = -1.0 / k_tang;
        scale3(temp_dbl, history);
      }
      temp_dbl = Fscrit / magfs;
      scale3(temp_dbl, fs);
    } else {
      zero3(fs);
    }
  }
}

/* ---------------------------------------------------------------------- */

void ContactModel::rolling(double *history)
{
  int rhist0, rhist1, rhist2, frameupdate;
  double rolldotn, rollmag, prjmag, magfr, hist_temp[3], temp_dbl, temp_array[3];

  rhist0 = roll_history_index;
  rhist1 = rhist0 + 1;
  rhist2 = rhist1 + 1;

  Frcrit = mu_roll * Fncrit;

  if (history_update) {
    hist_temp[0] = history[rhist0];
    hist_temp[1] = history[rhist1];
    hist_temp[2] = history[rhist2];
    rolldotn = dot3(hist_temp, nx);

    frameupdate = fabs(rolldotn)*k_roll > EPSILON*Frcrit;
    if (frameupdate) { // rotate into tangential plane
      rollmag = len3(hist_temp);
      // projection
      temp_dbl = -rolldotn;
      scale3(temp_dbl, nx, temp_array);
      sub3(hist_temp, temp_array, hist_temp);

      // also rescale to preserve magnitude
      prjmag = len3(hist_temp);
      if (prjmag > 0) temp_dbl = rollmag / prjmag;
      else temp_dbl = 0;
      scale3(temp_dbl, hist_temp);
    }
    scale3(dt, vrl, temp_array);
    add3(hist_temp, temp_array, hist_temp);
  }

  scaleadd3(k_roll, hist_temp, gamma_roll, vrl, fr);
  negate3(fr);

  // rescale frictional displacements and forces if needed

  magfr = len3(fr);
  if (magfr > Frcrit) {
    rollmag = len3(hist_temp);
    if (rollmag != 0.0) {
      temp_dbl = -Frcrit / (magfr * k_roll);
      scale3(temp_dbl, fr, temp_array);
      add3(hist_temp, temp_array, hist_temp);

      temp_dbl = -gamma_roll/k_roll;
      scale3(temp_dbl, vrl, temp_array);
      add3(hist_temp, temp_array, hist_temp);

      temp_dbl = Frcrit / magfr;
      scale3(temp_dbl, fr);
    } else {
      zero3(fr);
    }
  }

  history[rhist0] = hist_temp[0];
  history[rhist1] = hist_temp[1];
  history[rhist2] = hist_temp[2];
}

/* ---------------------------------------------------------------------- */

void ContactModel::twisting_marshall(double *history)
{
  // Overwrite twist coefficients with derived values
  k_twist = 0.5 * k_tang * a * a; // eq 32 of Marshall paper
  gamma_twist = 0.5 * gamma_tang * a * a;
  mu_twist = TWOTHIRDS * a * mu_tang;

  twisting_SDS(history);
}

/* ---------------------------------------------------------------------- */

void ContactModel::twisting_SDS(double *history)
{
  double signtwist, Mtcrit;

  if (history_update) {
    history[twist_history_index] += magtwist * dt;
  }

  magtortwist = -k_twist * history[twist_history_index] - gamma_twist*magtwist; // M_t torque (eq 30)
  signtwist = (magtwist > 0) - (magtwist < 0);
  Mtcrit = mu_twist * Fncrit; // critical torque (eq 44)
  if (fabs(magtortwist) > Mtcrit) {
    history[twist_history_index] = (Mtcrit * signtwist - gamma_twist * magtwist) / k_twist;
    magtortwist = -Mtcrit * signtwist; // eq 34
  }
}

/* ----------------------------------------------------------------------
   compute pull-off distance (beyond contact) for a given radius and atom type
   use temporary variables since this does not use a specific contact geometry
------------------------------------------------------------------------- */

double ContactModel::pulloff_distance(double radi, double radj)
{
  double Ecaled, a_tmp, Reff_tmp;

  if (normal_model != JKR) return radi+radj;

  Reff_tmp = radi * radj / (radi + radj);
  if (Reff_tmp <= 0) return 0;
  Ecaled = k_norm * THREEQUARTERS;
  a_tmp = cbrt(9 * MY_PI * cohesion * Reff_tmp * Reff_tmp / (4 * Ecaled));
  return a_tmp * a_tmp / Reff_tmp - 2 * sqrt(MY_PI * cohesion * a_tmp / Ecaled);
}

}
