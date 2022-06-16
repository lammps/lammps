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

#include <cmath>
#include "contact.h"

namespace LAMMPS_NS {
  namespace Contact{

enum {HOOKE, HERTZ, HERTZ_MATERIAL, DMT, JKR};
enum {VELOCITY, MASS_VELOCITY, VISCOELASTIC, TSUJI};
enum {TANGENTIAL_NOHISTORY, TANGENTIAL_HISTORY,
      TANGENTIAL_MINDLIN, TANGENTIAL_MINDLIN_RESCALE,
      TANGENTIAL_MINDLIN_FORCE, TANGENTIAL_MINDLIN_RESCALE_FORCE};
enum {TWIST_NONE, TWIST_SDS, TWIST_MARSHALL};
enum {ROLL_NONE, ROLL_SDS};

#define PI27SQ 266.47931882941264802866    // 27*PI**2
#define THREEROOT3 5.19615242270663202362  // 3*sqrt(3)
#define SIXROOT6 14.69693845669906728801   // 6*sqrt(6)
#define INVROOT6 0.40824829046386307274    // 1/sqrt(6)
#define FOURTHIRDS (4.0/3.0)                 // 4/3
#define THREEQUARTERS 0.75                 // 3/4

#define EPSILON 1e-10

ContactModel::ContactModel()
{
  k_norm = cohesion = gamma_norm = 0.0;
  k_tang = gamma_tang = mu_tang = 0.0;
  k_roll = gamma_roll = mu_roll = 0.0;
  k_twist = gamma_twist = mu_twist = 0.0;

  limit_damping = 0;
  cutoff_type = 0.0;
  reset_contact();
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

  MathExtra::sub3(xi, xj, dx);
  rsq = MathExtra::lensq3(dx);
  radsum = radi + radj;
  Reff = radi*radj/radsum;

  touch = false;
  if (normal_model == JKR) touch = touch_JKR(touch);
  else touch = (rsq < radsum*radsum);

  return touch
}

/* ---------------------------------------------------------------------- */

void ContactModel::prep_contact()
{
  prep_flag = 1;

  // If it hasn't already been done, test if the contact exists
  if (check_flag != 1) touch = check_contact();
  if (!touch) return;

  // Set flags
  mindlin_rescale = mindlin_force = 0;
  if (tangential_model == TANGENTIAL_MINDLIN_RESCALE ||
      tangential_model == TANGENTIAL_MINDLIN_RESCALE_FORCE)
    mindlin_rescale = 1
  if (tangential_model == TANGENTIAL_MINDLIN_FORCE ||
      tangential_model == TANGENTIAL_MINDLIN_RESCALE_FORCE)
    mindlin_force = 1

  double temp[3];

  // Standard geometric quantities
  r = sqrt(rsq);
  rinv = 1.0/r;
  delta = radsum - r;
  dR = delta*Reff
  MathExtra::scale3(rinv, dx, nx);

  // relative translational velocity
  MathExtra::sub3(v[i], v[j], vr);

  // normal component
  vnnr = MathExtra::dot3(vr, nx); //v_R . n
  MathExtra::scale3(vnnr, nx, vn);

  // tangential component
  MathExtra::sub3(vr, vn, vt);

  // relative rotational velocity
  MathExtra::scaleadd3(radi, omegai, radj, omegaj, wr);

  // relative tangential velocities
  MathExtra::cross3(wr, nx, temp);
  MathExtra::sub3(vt, temp, vtr);
  vrel = MathExtra::len(vtr);

  if (roll_model != NONE || twist_model != NONE)
    MathExtra::sub3(omega[i], omega[j], relrot);

  if (roll_model != NONE) {
    // rolling velocity, see eq. 31 of Wang et al, Particuology v 23, p 49 (2015)
    // this is different from the Marshall papers, which use the Bagi/Kuhn formulation
    // for rolling velocity (see Wang et al for why the latter is wrong)
    vrl[0] = Reff * (relrot[1] * n[2] - relrot[2] * nx[1]);
    vrl[1] = Reff * (relrot[2] * n[0] - relrot[0] * nx[2]);
    vrl[2] = Reff * (relrot[0] * n[1] - relrot[1] * nx[0]);
  }

  if (twist_model != NONE) {
    // omega_T (eq 29 of Marshall)
    magtwist = MathExtra::dot3(relrot, nx);
  }
}


/* ---------------------------------------------------------------------- */

void ContactModel::calculate_forces(double *forces, double *torquesi, double *torquesj, double *history)
{
  // If it hasn't already been done, run prep calculations
  if (prep_flag != 1) prep_contact();
  if (!touch) {
    forces[0] = forces[1] = forces[2] = 0.0;
    return
  }

  //**********************************************
  // normal forces
  //**********************************************

  // Also calculates: a, knfac, Fncrit
  double Fne;
  if (normal_model == JKR) {
    Fne = normal_JKR();
  } else if (normal_model == DMT) {
    Fne = normal_DMT();
  } else if (normal_model == HERTZ || normal_model == HERTZ_MATERIAL) {
    Fne = normal_hertz();
  } else {
    Fne = normal_hooke();
  }

  // NOTE: consider restricting Hooke to only have
  // 'velocity' as an option for damping?
  // Also calculates: damp_normal_prefactor
  double Fdamp = normal_damping();

  double Fntot = Fne + Fdamp;
  if (limit_damping && (Fntot < 0.0)) Fntot = 0.0;

  //**********************************************
  // tangential force, including history effects
  //**********************************************
  // For linear, mindlin, mindlin_rescale:
  // history = cumulative tangential displacement
  //
  // For mindlin/force, mindlin_rescale/force:
  // history = cumulative tangential elastic force

  Fncrit = critical_normal(Fne, Fntot, geom, model);
  Fscrit = mu_tang * Fncrit;

  if (tangential_model == TANGENTIAL_NOHISTORY) {
    tangential_no_history();
  } else if (tangential_model == TANGENTIAL_HISTORY) {
    tangential_history(history);
  }
  tangential_forces(history);

  //**********************************************
  // rolling force
  //**********************************************

  if (roll_model != ROLL_NONE)
    rolling(history);

  //****************************************
  // twisting torque, including history effects
  //****************************************

  if (twist_model == TWIST_MARSHALL) {
    twist_marshall(history);
  } else if (twist_model == TWIST_SDS) {
    twist_SDS(history);
  }

  //**********************************************
  // sum contributions
  //**********************************************

  MathExtra::scale3(Fntot, nx, forces);
  MathExtra::add3(forces, fs, forces);

  MathExtra::cross3(nx, fs, torquesi);
  MathExtra::copy3(torquesi, torquesj);

  double dist_to_contact = radi-0.5*delta;
  MathExtra::scale3(dist_to_contact, torquesi);
  dist_to_contact = radj-0.5*delta;
  MathExtra::scale3(dist_to_contact, torquesj);

  double torroll[3];
  if (roll_model != ROLL_NONE) {
    MathExtra::cross3(nx, fr, torroll);
    MathExtra::scale3(Reff, torroll);
    MathExtra::add3(torquesi, torroll, torquesi);
    MathExtra::sub3(torquesj, torroll, torquesj);
  }

  double tortwist[3];
  if (twist_model != NONE_TWIST) {
    MathExtra::scale3(magtortwist, nx, tortwist);
    MathExtra::add3(torquesi, tortwist, torquesi)
    MathExtra::sub3(torquesj, tortwist, torquesj)
  }
}

/* ---------------------------------------------------------------------- */

double ContactModel::touch_JKR(int touch)
{
  double Escaled, R2, delta_pulloff, dist_pulloff;
  bool touchflag;

  Escaled = E*THREEQUARTERS;
  if (touch) {
    R2 = Reff * Reff;
    a = cbrt(9.0 * MY_PI * cohesion * R2 / (4 * Escaled));
    delta_pulloff = a * a / Reff - 2 * sqrt(MY_PI * cohesion * a / Escaled);
    dist_pulloff = radsum - delta_pulloff;
    touchflag = (rsq < dist_pulloff * dist_pulloff);
  } else {
    touchflag = (rsq < radsum * radsum);
  }
  return touchflag;
}

/* ---------------------------------------------------------------------- */

double ContactModel::normal_JKR()
{
  double R2, dR2, t0, t1, t2, t3, t4, t5, t6;
  double sqrt1, sqrt2, sqrt3, a, a2, F_pulloff;

  R2 = Reff * Reff;
  dR2 = dR * dR;
  t0 = cohesion * cohesion * R2 * R2 * E;
  t1 = PI27SQ*t0;
  t2 = 8 * dR * dR2 * E * E *  E;
  t3 = 4 * dR2 * E;

  // in case sqrt(0) < 0 due to precision issues
  sqrt1 = MAX(0, t0 * (t1 + 2 * t2));
  t4 = cbrt(t1 + t2 + THREEROOT3 * MY_PI * sqrt(sqrt1));
  t5 = t3 / t4 + t4 / E;
  sqrt2 = MAX(0, 2 * dR + t5);
  t6 = sqrt(sqrt2);
  sqrt3 = MAX(0, 4 * dR - t5 + SIXROOT6 * cohesion * MY_PI * R2 / (E * t6));
  a = INVROOT6 * (t6 + sqrt(sqrt3));
  a2 = a * a;
  double Fne = E * a * a2 / Reff - MY_2PI * a2 * sqrt(4 * cohesion * E / (MY_PI * a));
  double F_pulloff = 3 * MY_PI * cohesion * Reff;

  knfac = E * a;
  Fncrit = fabs(Fne + 2*F_pulloff);
  return Fne;
}

/* ---------------------------------------------------------------------- */

double ContactModel::normal_DMT(double &Fne)
{
  a = sqrt(dR);
  double Fne = a * E * delta;
  Fne -= 4 * MY_PI * cohesion * Reff;
  F_pulloff = 4 * MY_PI * cohesion * Reff;

  knfac = E * a;
  Fncrit = fabs(Fne + 2*F_pulloff);
  return Fne;
}

/* ---------------------------------------------------------------------- */

double  ContactModel::normal_Hertz(double &Fne)
{
  a = sqrt(dR);
  double Fne = E * delta * a;

  knfac = E * a;
  return Fne;
}

/* ---------------------------------------------------------------------- */

double  ContactModel::normal_Hooke(double &Fne)
{
  a = sqrt(dR);
  double Fne = E * delta;

  knfac = E;
  return Fne;
}

/* ---------------------------------------------------------------------- */

double ContactModel::normal_damping()
{
  double damp_normal;
  if (damping_model == VELOCITY) {
    damp_normal = 1;
  } else if (damping_model == MASS_VELOCITY) {
    damp_normal = meff;
  } else if (damping_model == VISCOELASTIC) {
    damp_normal = a * meff;
  } else if (damping_model == TSUJI) {
    damp_normal = sqrt(meff * knfac);
  } else damp_normal = 0.0;

  damp_normal_prefactor = gamma_norm * damp_normal;
  return -damp_normal_prefactor * vnnr;
}

/* ---------------------------------------------------------------------- */

double ContactModel::tangential_no_history()
{
  double gamma_scaled = gamma_tang * damp_normal_prefactor;
  double fsmag, Ft;

  // classic pair gran/hooke (no history)
  fsmag = gamma_scaled * vrel;
  if (vrel != 0.0) Ft = MIN(Fscrit,fsmag) / vrel;
  else Ft = 0.0;

  Ft = -Ft;
  MathExtra::scale3(Ft, vtr, fs);
}

/* ---------------------------------------------------------------------- */

double ContactModel::tangential_forces(double *history)
{
  double gamma_scaled = gamma_tang * damp_normal_prefactor;
  double k = k_tang;
  int frame_update = 0;
  double fsmag, rsht, shrmag, prjmag, temp_dbl, temp_array[3];

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (history_update) {
    rsht = MathExtra::dot(history, nx);
    frameupdate = fabs(rsht) * k > EPSILON * Fscrit;

    if (frameupdate) {
      shrmag = MathExtra::len3(history);
      // projection
      MathExtra::scale3(rsht, nx, history);
      // also rescale to preserve magnitude
      prjmag = MathExtra::len3(history);
      if (prjmag > 0) temp_double = shrmag / prjmag;
      else temp_double = 0;
      MathExtra::scale3(temp_double, history);
    }

    // update history
    // tangential force
    // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
    temp_dbl = k * dt;
    MathExtra::scale3(temp_dbl, vtr, temp_array);
    MathExtra::sub3(history, temp_array, history);
  }

  // tangential forces = history + tangential velocity damping
  temp_double = -gamma;
  MathExtra::scale3(temp_double, vtr, fs);

  // rescale frictional displacements and forces if needed
  fsmag = MathExtra::len3(fs);
  if (fsmag > Fscrit) {
    shrmag = MathExtra::len3(history);
    if (shrmag != 0.0) {
      temp_double = Fscrit / magfs;
      MathExtra::scale3(temp_double, fs, history);
      MathExtra::scale3(gamma, vtr, temp_array);
      MathExtra::add3(history, temp_array, history);
      temp_double = Fscrit / magfs;
      MathExtra::scale3(temp_double, fs);
    } else {
      MathExtra::zero3(fs);
    }
  }
}

/* ---------------------------------------------------------------------- */

double ContactModel::tangential_mindlin(double *history)
{
  double k_scaled, gamma_scaled, fsmag, rsht, shrmag, prjmag, temp_dbl;
  double temp_array[3];
  int frame_update = 0;

  gamma_scaled = gamma_tang * damp_normal_prefactor;
  k_scaled = k_tange * a;
  if (mindlin_rescale) {
    // on unloading, rescale the shear displacements/force
    if (a < history[3]) {
      temp_double = a / history[3];
      MathExtra::scale3(temp_double, history);
    }
  }

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (history_update) {
    rsht = MathExtra::dot(history, nx);
    if (mindlin_force)
      frameupdate = fabs(rsht) > EPSILON * Fscrit;
    else
      frameupdate = fabs(rsht) * k_scaled > EPSILON * Fscrit;

    if (frameupdate) {
      shrmag = MathExtra::len3(history);
      // projection
      MathExtra::scale3(rsht, nx, history);
      // also rescale to preserve magnitude
      prjmag = MathExtra::len3(history);
      if (prjmag > 0) temp_double = shrmag / prjmag;
      else temp_double = 0;
      MathExtra::scale3(temp_double, history);
    }

    // update history
    if (mindlin_force) {
      // tangential force
      // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
      temp_dbl = -k_scaled * dt;
      MathExtra::scale3(temp_dbl, vtr, temp_array);
    } else {
      MathExtra::scale3(dt, vtr, temp_array);
    }
    MathExtra::add3(history, temp_array, history);

    if (mindlin_rescale) history[3] = a;
  }

  // tangential forces = history + tangential velocity damping
  temp_double = -gamma_scaled;
  MathExtra::scale3(temp_double, vtr, fs);

  if (! mindlin_force) {
    MathExtra::scale3(k_scaled, history, temp_array);
    MathExtra::add3(fs, temp_array, fs);
  }

  // rescale frictional displacements and forces if needed
  fsmag = MathExtra::len3(fs);
  if (fsmag > Fscrit) {
    shrmag = MathExtra::len3(history);
    if (shrmag != 0.0) {
      temp_double = Fscrit / magfs;
      MathExtra::scale3(temp_double, fs, history);
      MathExtra::scale3(gamma, vtr, temp_array);
      MathExtra::add3(history, temp_array, history);
      if (! mindlin_force) {
        temp_double = -1.0 / k;
        MathExtra::scale3(temp_double, history);
      }
      temp_double = Fscrit / magfs;
      MathExtra::scale3(temp_double, fs);
    } else {
      MathExtra::zero3(fs);
    }
  }
}

/* ---------------------------------------------------------------------- */

double ContactModel::rolling(double *history)
{
  int rhist0, rhist1, rhist2, frameupdate;
  double rolldotn, rollmag, prjmag, magfr, hist_temp[3], temp_array[3];

  rhist0 = roll_history_index;
  rhist1 = rhist0 + 1;
  rhist2 = rhist1 + 1;

  Frcrit = mu_roll * Fncrit;

  if (history_update) {
    hist_temp[0] = history[rhist0];
    hist_temp[1] = history[rhist1];
    hist_temp[2] = history[rhist2];
    rolldotn = MathExtra::dot3(hist_temp, nx);

    frameupdate = fabs(rolldotn)*k_roll > EPSILON*Frcrit;
    if (frameupdate) { // rotate into tangential plane
      rollmag = MathExtra::len3(temp);
      // projection
      temp_double = -rolldotn;
      MathExtra::scale3(temp_double, nx, temp);

      // also rescale to preserve magnitude
      prjmag = MathExtra::len3(hist_temp);
      if (prjmag > 0) temp_double = rollmag / prjmag;
      else temp_double = 0;
      MathExtra::scale3(temp_double, hist_temp);
    }
    MathExtra::scale3(dt, vrl, temp_array);
    MathExtra::add3(hist_temp, temp_array, hist_temp)
  }

  MathExtra::scaleadd3(k_roll, hist_tepm, gamma_roll, vrl, fr);
  MathExtra::negate3(fr);

  // rescale frictional displacements and forces if needed

  magfr = MathExtra::len3(fr);
  if (magfr > Frcrit) {
    rollmag = MathExtra::len3(temp);
    if (rollmag != 0.0) {
      temp_double = -Frcrit / (magfr * k_roll);
      MathExtra::scale3(temp_double, fr, temp_array);
      MathExtra::add3(hist_temp, temp_array, hist_temp);

      temp_double = -gamma_roll/k_roll;
      MathExtra::scale3(temp_double, vrl, temp_array);
      MathExtra::add3(hist_temp, temp_array, hist_temp)

      temp_double = Frcrit / magfr;
      MathExtra::scale3(temp_double, fr);
    } else {
      MathExtra::zero3(fr);
    }
  }

  history[rhist0] = hist_temp[0];
  history[rhist1] = hist_temp[1];
  history[rhist2] = hist_temp[2];
}

/* ---------------------------------------------------------------------- */

double ContactModel::twisting_marshall(double *history)
{
  // Overwrite twist coefficients with derived values
  k_twist = 0.5 * k_tangential * a * a; // eq 32 of Marshall paper
  gamma_twist = 0.5 * gamma_tangential * a * a;
  mu_twist = TWOTHIRDS * a * mu_tangential;

  twisting_SDS(history);
}

/* ---------------------------------------------------------------------- */

double ContactModel::twisting_SDS(double *history)
{
  double signtwist, Mtcrit;

  if (historyupdate) {
    history[twist_history_index] += magtwist * dt;
  }

  magtortwist = -k_twist * history[twist_history_index] - damp_twist*magtwist; // M_t torque (eq 30)
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
  Ecaled = E * THREEQUARTERS;
  a_tmp = cbrt(9 * MY_PI * cohesion * Reff_tmp * Reff_tmp / (4 * Ecaled));
  return a_tmp * a_tmp / Reff_tmp - 2 * sqrt(MY_PI * cohesion * a_tmp / Ecaled);
}

}}
