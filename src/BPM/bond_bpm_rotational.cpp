/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Joel Clemmer (SNL)
                        Gabriel Alkuino (SyracuseU)
                        Teng Zhang (SyracuseU)
------------------------------------------------------------------------- */

#include "bond_bpm_rotational.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_bond_history.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

static constexpr double EPSILON = 1e-10;

enum { DEM, DERIVATIVE };
enum { AVERAGE, PARTICLE };

using namespace LAMMPS_NS;
using MathConst::MY_SQRT2;
using MathConst::MY_PI;

/* ---------------------------------------------------------------------- */

static double acos_limit(double c)
{
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;
  return acos(c);
}

static double wrap_2pi(double t)
{
  if (t > MY_PI) t -= 2.0 * MY_PI;
  if (t < -MY_PI) t += 2.0 * MY_PI;
  return t;
}

/* ----------------------------------------------------------------------
   Average of two quaternions using shortest path SLERP.
------------------------------------------------------------------------- */
static void quat_average(double *qA, double *qB, double *q_out)
{
  double qB_use[4] = {qB[0], qB[1], qB[2], qB[3]};
  double dot = qA[0] * qB[0] + qA[1] * qB[1] + qA[2] * qB[2] + qA[3] * qB[3];
  if (dot < 0.0) {
    qB_use[0] = -qB[0]; qB_use[1] = -qB[1];
    qB_use[2] = -qB[2]; qB_use[3] = -qB[3];
    dot = -dot;
  }
  double s = 1.0 / sqrt(2.0 * (1.0 + dot));
  q_out[0] = (qA[0] + qB_use[0]) * s;
  q_out[1] = (qA[1] + qB_use[1]) * s;
  q_out[2] = (qA[2] + qB_use[2]) * s;
  q_out[3] = (qA[3] + qB_use[3]) * s;
  MathExtra::qnormalize(q_out);
}

/* ----------------------------------------------------------------------
   Decompose quaternion into swing and twist using Wang's formula
------------------------------------------------------------------------- */
static void decompose_swing_twist(double *quat, double *swing, double *twist)
{
  double tmp = sqrt(quat[0] * quat[0] + quat[3] * quat[3]);
  double psi;
  if (tmp < EPSILON) {
    psi = 0.0;
  } else {
    psi = 2.0 * acos_limit(quat[0] / tmp);
  }
  if (quat[3] < 0.0) psi = -psi;

  double c = cos(psi / 2.0);
  double s = sin(psi / 2.0);
  twist[0] = c; twist[1] = 0.0; twist[2] = 0.0; twist[3] = s;

  double twist_inv[4];
  MathExtra::qconjugate(twist, twist_inv);
  MathExtra::quatquat(twist_inv, quat, swing);
  MathExtra::qnormalize(swing);
}

/* ---------------------------------------------------------------------- */

BondBPMRotational::BondBPMRotational(LAMMPS *_lmp) :
    BondBPM(_lmp), Kr(nullptr), Ks(nullptr), Kt(nullptr), Kb(nullptr), gr(nullptr),
    gs(nullptr), gt(nullptr), gb(nullptr), Fcr(nullptr), Fcs(nullptr), Tct(nullptr),
    Tcb(nullptr)
{
  partial_flag = 1;
  smooth_flag = 1;
  normalize_flag = 0;
  writedata = 0;

  frame_style = AVERAGE;
  damping_style = DERIVATIVE;

  update_flag = 1;
  nhistory = 7;
  id_fix_bond_history = utils::strdup("HISTORY_BPM_ROTATIONAL");

  single_extra = 7;
  svector = new double[7];
}

/* ---------------------------------------------------------------------- */

BondBPMRotational::~BondBPMRotational()
{
  delete[] svector;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(Kr);
    memory->destroy(Ks);
    memory->destroy(Kt);
    memory->destroy(Kb);
    memory->destroy(Fcr);
    memory->destroy(Fcs);
    memory->destroy(Tct);
    memory->destroy(Tcb);
    memory->destroy(gr);
    memory->destroy(gs);
    memory->destroy(gt);
    memory->destroy(gb);
  }
}

/* ----------------------------------------------------------------------
  Store data for a single bond - if bond added after LAMMPS init
------------------------------------------------------------------------- */

double BondBPMRotational::store_bond(int n, int i, int j)
{
  double delx, dely, delz, r, rinv;
  double **x = atom->x;
  tagint *tag = atom->tag;
  double **bondstore = fix_bond_history->bondstore;

  // ri points from atom with smaller tag to larger tag
  if (tag[i] < tag[j]) {
    delx = x[j][0] - x[i][0];
    dely = x[j][1] - x[i][1];
    delz = x[j][2] - x[i][2];
  } else {
    delx = x[i][0] - x[j][0];
    dely = x[i][1] - x[j][1];
    delz = x[i][2] - x[j][2];
  }

  r = sqrt(delx * delx + dely * dely + delz * delz);
  rinv = 1.0 / r;

  bondstore[n][0] = r;
  bondstore[n][1] = delx * rinv;
  bondstore[n][2] = dely * rinv;
  bondstore[n][3] = delz * rinv;
  if (damping_style == DERIVATIVE) {
    bondstore[n][4] = 0.0;  // gamma
    bondstore[n][5] = 0.0;  // theta
    bondstore[n][6] = 0.0;  // psi
  }

  if (i < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[i]; m++) {
      if (atom->bond_atom[i][m] == tag[j]) {
        fix_bond_history->update_atom_value(i, m, 0, r);
        fix_bond_history->update_atom_value(i, m, 1, delx * rinv);
        fix_bond_history->update_atom_value(i, m, 2, dely * rinv);
        fix_bond_history->update_atom_value(i, m, 3, delz * rinv);
        if (damping_style == DERIVATIVE)
          for (int a = 4; a < 7; a++)
            fix_bond_history->update_atom_value(i, m, a, 0.0);
      }
    }
  }

  if (j < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[j]; m++) {
      if (atom->bond_atom[j][m] == tag[i]) {
        fix_bond_history->update_atom_value(j, m, 0, r);
        fix_bond_history->update_atom_value(j, m, 1, delx * rinv);
        fix_bond_history->update_atom_value(j, m, 2, dely * rinv);
        fix_bond_history->update_atom_value(j, m, 3, delz * rinv);
        if (damping_style == DERIVATIVE)
          for (int a = 4; a < 7; a++)
            fix_bond_history->update_atom_value(j, m, a, 0.0);
      }
    }
  }
  return r;
}

/* ----------------------------------------------------------------------
  Store data for all bonds called once
------------------------------------------------------------------------- */

void BondBPMRotational::store_data()
{
  int i, j, m, type;
  double delx, dely, delz, r, rinv;
  double **x = atom->x;
  int **bond_type = atom->bond_type;
  tagint *tag = atom->tag;

  for (i = 0; i < atom->nlocal; i++) {
    for (m = 0; m < atom->num_bond[i]; m++) {
      type = bond_type[i][m];
      if (type <= 0) continue;

      j = atom->map(atom->bond_atom[i][m]);
      if (j == -1) error->one(FLERR, "Atom missing in BPM bond");

      // ri points from smaller tag to larger tag
      if (tag[i] < tag[j]) {
        delx = x[j][0] - x[i][0];
        dely = x[j][1] - x[i][1];
        delz = x[j][2] - x[i][2];
      } else {
        delx = x[i][0] - x[j][0];
        dely = x[i][1] - x[j][1];
        delz = x[i][2] - x[j][2];
      }

      domain->minimum_image(FLERR, delx, dely, delz);
      r = sqrt(delx * delx + dely * dely + delz * delz);
      rinv = 1.0 / r;

      fix_bond_history->update_atom_value(i, m, 0, r);
      fix_bond_history->update_atom_value(i, m, 1, delx * rinv);
      fix_bond_history->update_atom_value(i, m, 2, dely * rinv);
      fix_bond_history->update_atom_value(i, m, 3, delz * rinv);

      if (damping_style == DERIVATIVE)
        for (int a = 4; a < 7; a++)
          fix_bond_history->update_atom_value(i, m, a, 0.0);
    }
  }
}

/* ----------------------------------------------------------------------
  Construction of elastic forces + derivative-based damping based on
  Alkuino et al. 2026
------------------------------------------------------------------------- */

double BondBPMRotational::average_frame_forces(int i1, int i2, int type,
                                                double *ri, double *rf,
                                                double *force1on2, double *torque1on2,
                                                double *torque2on1, double &ebond,
                                                double *bondstore)
{
  double **quat = atom->quat;
  double **v = atom->v;

  double Kr_type = Kr[type];
  double Ks_type = Ks[type];

  double ri_norm = MathExtra::len3(ri);
  if (normalize_flag) {
    double ri_inv = 1.0 / ri_norm;
    Kr_type *= ri_inv;
    Ks_type *= ri_inv;
  }

  // Store prior values for damping
  double gamma_prior, theta_prior, psi_prior;
  if (damping_style == DERIVATIVE) {
    gamma_prior = bondstore[4];
    theta_prior = bondstore[5];
    psi_prior = bondstore[6];
  }

  // Current quaternions
  double q1f[4], q2f[4];
  q1f[0] = quat[i1][0]; q1f[1] = quat[i1][1]; q1f[2] = quat[i1][2]; q1f[3] = quat[i1][3];
  q2f[0] = quat[i2][0]; q2f[1] = quat[i2][1]; q2f[2] = quat[i2][2]; q2f[3] = quat[i2][3];

  // qcf: global to C frame
  double qcf[4], qcf_inv[4];
  quat_average(q1f, q2f, qcf);
  MathExtra::qconjugate(qcf, qcf_inv);

  // Radial force in C frame
  double rc[3];
  MathExtra::quatrotvec(qcf_inv, rf, rc);
  double rc_norm = MathExtra::len3(rc);
  double Fr_c[3];
  double Fr_elastic_mag = Kr_type * (rc_norm - ri_norm);
  double Fr_mag = Fr_elastic_mag;

  // Shear force in C frame
  double ri_dot_rc = MathExtra::dot3(ri, rc);
  double c = ri_dot_rc / (ri_norm * rc_norm);
  double gamma = acos_limit(c);

  double rixrc[3], s_hat[3], t_hat[3];
  MathExtra::cross3(ri, rc, rixrc);
  double rixrc_mag = MathExtra::len3(rixrc);

  double Fst_c[3] = {0.0, 0.0, 0.0};
  double Tst_c[3] = {0.0, 0.0, 0.0};

  double rc_hat[3];
  MathExtra::normalize3(rc, rc_hat);
  double Fs_elastic_mag = 0.0;
  double Fs_mag = 0.0;
  if (rixrc_mag > EPSILON) {
    MathExtra::normalize3(rixrc, t_hat);
    MathExtra::cross3(t_hat, rc_hat, s_hat);
    MathExtra::norm3(s_hat);
    Fs_elastic_mag = Ks_type * rc_norm * gamma;
    Fs_mag = Fs_elastic_mag;
  }

  // Orientation of atom2 w.r.t. C frame
  double qu0[4];
  MathExtra::quatquat(qcf_inv, q2f, qu0);
  MathExtra::qnormalize(qu0);

  // Quaternion that rotates z to rc direction
  double qm[4];
  double temp = rc_norm + rc[2];
  if (temp < 0.0) temp = 0.0;
  qm[0] = MY_SQRT2 * 0.5 * sqrt(temp / rc_norm);

  temp = sqrt(rc[0]*rc[0] + rc[1]*rc[1]);
  if (temp > EPSILON) {
    double temp2 = rc_norm - rc[2];
    if (temp2 < 0.0) temp2 = 0.0;
    double factor = -MY_SQRT2 * 0.5 * sqrt(temp2 / rc_norm) / temp;
    qm[1] = factor * rc[1];
    qm[2] = -factor * rc[0];
  } else {
    qm[1] = 0.0;
    qm[2] = 0.0;
  }
  qm[3] = 0.0;
  if (qm[0] == 0.0 && qm[1] == 0.0 && qm[2] == 0.0) {
    qm[2] = 1.0;
  }
  MathExtra::qnormalize(qm);

  // qu = qm^-1 * qu0 * qm
  double qm_inv[4], qu[4], qtmp[4];
  MathExtra::qconjugate(qm, qm_inv);
  MathExtra::quatquat(qm_inv, qu0, qtmp);
  MathExtra::quatquat(qtmp, qm, qu);
  MathExtra::qnormalize(qu);

  // Swing-twist decomposition
  double swing2[4], twist2[4];
  decompose_swing_twist(qu, swing2, twist2);

  double qu_inv[4];
  MathExtra::qconjugate(qu, qu_inv);
  double swing1[4], twist1[4];
  decompose_swing_twist(qu_inv, swing1, twist1);

  // Combined twist = twist2 * twist1^-1
  double twist1_inv[4], twist[4];
  MathExtra::qconjugate(twist1, twist1_inv);
  MathExtra::quatquat(twist2, twist1_inv, twist);
  double psi = 2.0 * acos_limit(twist[0]);

  double taxis[3];
  if (fabs(psi) > EPSILON) {
    double twist_vec[3] = {twist[1], twist[2], twist[3]};
    MathExtra::normalize3(twist_vec, taxis);
  } else {
    taxis[0] = 0.0; taxis[1] = 0.0; taxis[2] = 1.0;
  }

  // Combined swing = swing2 * swing1^-1
  double swing1_inv[4], swing[4];
  MathExtra::qconjugate(swing1, swing1_inv);
  MathExtra::quatquat(swing2, swing1_inv, swing);
  double theta = 2.0 * acos_limit(swing[0]);

  double baxis[3];
  if (fabs(theta) > EPSILON) {
    double swing_vec[3] = {swing[1], swing[2], swing[3]};
    MathExtra::normalize3(swing_vec, baxis);
  } else {
    baxis[0] = 0.0; baxis[1] = 0.0; baxis[2] = 1.0;
  }

  // Torques in C' frame with damping
  double Tt_p[3], Tb_p[3];
  double Tt_elastic_mag = Kt[type] * psi;
  double Tt_mag = Tt_elastic_mag;
  double Tb_elastic_mag = Kb[type] * theta;
  double Tb_mag = Tb_elastic_mag;

  // Damping forces, if relevant

  if (damping_style == DERIVATIVE) {
    double dt_inv = 1.0 / update->dt;

    // radial
    double dv[3], dv_c[3];
    MathExtra::sub3(v[i2], v[i1], dv);
    MathExtra::quatrotvec(qcf_inv, dv, dv_c);
    double dvdotr = MathExtra::dot3(dv_c, rc_hat);
    Fr_mag += gr[type] * dvdotr;

    // shear
    if (rixrc_mag > EPSILON) {
      double dgamma = wrap_2pi(gamma - gamma_prior) * dt_inv;
      Fs_mag += gs[type] * rc_norm * dgamma;
    }

    // twisting
    double dpsi = wrap_2pi(psi - psi_prior) * dt_inv;
    Tt_mag += gt[type] * dpsi;

    // bending
    double dtheta = wrap_2pi(theta - theta_prior) * dt_inv;
    Tb_mag += gb[type] * dtheta;

    bondstore[4] = gamma;
    bondstore[5] = theta;
    bondstore[6] = psi;
  }

  // Calculate net radial/shear forces (after damping)

  MathExtra::scale3(Fr_mag, rc_hat, Fr_c);

  if (rixrc_mag > EPSILON) {
    MathExtra::scale3(Fs_mag, s_hat, Fst_c);
    double Tst_mag = 0.5 * rc_norm * Fs_mag;
    MathExtra::scale3(Tst_mag, t_hat, Tst_c);
  }

  MathExtra::scale3(Tt_mag, taxis, Tt_p);
  MathExtra::scale3(Tb_mag, baxis, Tb_p);

  // Rotate torques to C frame
  double Tb_c[3], Tt_c[3];
  MathExtra::quatrotvec(qm, Tb_p, Tb_c);
  MathExtra::quatrotvec(qm, Tt_p, Tt_c);

  // Total forces and torques in C frame
  double F1_c[3], T1_c[3], T2_c[3], Ttmp[3];
  MathExtra::add3(Fr_c, Fst_c, F1_c);
  MathExtra::add3(Tb_c, Tt_c, Ttmp);
  MathExtra::add3(Tst_c, Ttmp, T1_c);
  MathExtra::sub3(Tst_c, Ttmp, T2_c);

  // Rotate to global frame
  double force2on1[3];
  MathExtra::quatrotvec(qcf, F1_c, force2on1);
  MathExtra::scale3(-1.0, force2on1, force1on2);

  MathExtra::quatrotvec(qcf, T1_c, torque2on1);
  MathExtra::quatrotvec(qcf, T2_c, torque1on2);

  // Breaking criterion (use undamped magnitudes)
  double breaking = fabs(Fr_elastic_mag) / Fcr[type] + Fs_elastic_mag / Fcs[type] +
                    Tb_elastic_mag / Tcb[type] + Tt_elastic_mag / Tct[type];
  if (breaking < 0.0) breaking = 0.0;

  // Approximate bond energy:
  ebond = 0.5 * (Kr_type * (rc_norm - ri_norm) * (rc_norm - ri_norm) +
	      Ks_type * rc_norm * rc_norm * gamma * gamma +
        Kt[type] * psi * psi +
        Kb[type] * theta * theta);

  return breaking;
}

/* ----------------------------------------------------------------------
  Construction of elastic forces + derivative-based damping based on:
    1) Y. Wang Acta Geotechnica 2009
    2) P. Mora & Y. Wang Advances in Geomcomputing 2009
------------------------------------------------------------------------- */

double BondBPMRotational::particle_frame_forces(int i1, int i2, int type,
                                            double *r0_neg, double *r,
                                            double *force1on2, double *torque1on2,
                                            double *torque2on1, double &ebond,
                                            double *bondstore)
{
  double breaking, temp, r0_dot_rb, c;
  double cos_phi, sin_phi;
  double mag_in_plane, mag_out_plane;
  double Fs_mag, Tt_mag, Tb_mag;

  double q1[4], q2[4];
  double q2inv[4], mq[4], mqinv[4], qp21[4], q21[4], qtmp[4];
  double r0[3], rb[3], rb_x_r0[3], s[3], t[3];
  double Fr, Fn[3], Fs[3], Fsr[3], Fsrt[3], Fsrp[3], F_rot[3];
  double Tsr[3], Tst[3], Tb[3], Tt[3], Tbp[3], Ttp[3], Tsrp[3], T_rot[3];

  // Store prior values for damping
  double gamma_prior, theta_prior, psi_prior;
  if (damping_style == DERIVATIVE) {
    gamma_prior = bondstore[4];
    theta_prior = bondstore[5];
    psi_prior = bondstore[6];
  }

  double **quat = atom->quat;
  double r_mag = MathExtra::len3(r);
  double r_mag_inv = 1.0 / r_mag;
  double r0_mag = MathExtra::len3(r0_neg);
  double r0_mag_inv = 1.0 / r0_mag;
  double Kr_type = Kr[type];
  double Ks_type = Ks[type];
  if (normalize_flag) {
    Kr_type *= r0_mag_inv;
    Ks_type *= r0_mag_inv;
  }

  q1[0] = quat[i1][0];
  q1[1] = quat[i1][1];
  q1[2] = quat[i1][2];
  q1[3] = quat[i1][3];

  q2[0] = quat[i2][0];
  q2[1] = quat[i2][1];
  q2[2] = quat[i2][2];
  q2[3] = quat[i2][3];

  // Calculate normal forces, rb = bond vector in particle 1's frame
  MathExtra::qconjugate(q2, q2inv);
  MathExtra::quatrotvec(q2inv, r, rb);
  MathExtra::negate3(rb); // Note, reverse of Mora & Wang and Alkuino et al.
  MathExtra::scale3(-1.0, r0_neg, r0);
  Fr = Kr_type * (r_mag - r0_mag);

  // Calculate forces due to tangential displacements (no rotation)
  r0_dot_rb = MathExtra::dot3(r0, rb);
  c = r0_dot_rb * r_mag_inv * r0_mag_inv;
  double gamma = acos_limit(c);

  MathExtra::cross3(rb, r0, rb_x_r0);
  MathExtra::cross3(rb, rb_x_r0, s);
  MathExtra::norm3(s);

  Fs_mag = Ks_type * r_mag * gamma;

  // Calculate torque due to tangential displacements
  MathExtra::cross3(r0, rb, t);
  MathExtra::norm3(t);

  // Relative rotation force/torque
  // Use representation of X'Y'Z' rotations from Wang, Mora 2009
  temp = r_mag + rb[2];
  if (temp < 0.0) temp = 0.0;
  mq[0] = MY_SQRT2 * 0.5 * sqrt(temp * r_mag_inv);

  temp = sqrt(rb[0] * rb[0] + rb[1] * rb[1]);
  if (temp != 0.0) {
    mq[1] = -MY_SQRT2 * 0.5 / temp;
    temp = r_mag - rb[2];
    if (temp < 0.0) temp = 0.0;
    mq[1] *= sqrt(temp * r_mag_inv);
    mq[2] = -mq[1];
    mq[1] *= rb[1];
    mq[2] *= rb[0];
  } else {
    // If aligned along z axis, x,y terms zero (r_mag-rb[2] = 0)
    mq[1] = 0.0;
    mq[2] = 0.0;
  }
  mq[3] = 0.0;

  // qp21 = opposite of r^\circ_21 in Wang
  // q21 = opposite of r_21 in Wang
  MathExtra::quatquat(q2inv, q1, qp21);
  MathExtra::qconjugate(mq, mqinv);
  MathExtra::quatquat(mqinv, qp21, qtmp);
  MathExtra::quatquat(qtmp, mq, q21);

  double psi;
  temp = sqrt(q21[0] * q21[0] + q21[3] * q21[3]);
  if (temp != 0.0) {
    psi = 2.0 * acos_limit(q21[0] / temp);
  } else {
    psi = 0.0;
  }

  // Map negative rotations
  if (q21[3] < 0.0)    // sin = q21[3]/temp
    psi = -psi;

  if (q21[3] == 0.0) psi = 0.0;

  c = q21[0] * q21[0] - q21[1] * q21[1] - q21[2] * q21[2] + q21[3] * q21[3];
  double theta = acos_limit(c);

  // Separately calculate magnitude of quaternion in x-y and out of x-y planes
  // to avoid dividing by zero
  mag_out_plane = (q21[0] * q21[0] + q21[3] * q21[3]);
  mag_in_plane = (q21[1] * q21[1] + q21[2] * q21[2]);

  if (mag_in_plane == 0.0) {
    // No rotation => no bending/shear torque or extra shear force
    // achieve by setting cos/sin = 0
    cos_phi = 0.0;
    sin_phi = 0.0;
  } else if (mag_out_plane == 0.0) {
    // Calculate angle in plane
    cos_phi = q21[2] / sqrt(mag_in_plane);
    sin_phi = -q21[1] / sqrt(mag_in_plane);
  } else {
    // Default equations in Mora, Wang 2009
    cos_phi = q21[1] * q21[3] + q21[0] * q21[2];
    sin_phi = q21[2] * q21[3] - q21[0] * q21[1];

    cos_phi /= sqrt(mag_out_plane * mag_in_plane);
    sin_phi /= sqrt(mag_out_plane * mag_in_plane);
  }

  Tbp[0] = -Kb[type] * theta * sin_phi;
  Tbp[1] = Kb[type] * theta * cos_phi;
  Tbp[2] = 0.0;

  Ttp[0] = 0.0;
  Ttp[1] = 0.0;
  Ttp[2] = Kt[type] * psi;

  Fsrp[0] = -0.5 * Ks_type * r_mag * theta * cos_phi;
  Fsrp[1] = -0.5 * Ks_type * r_mag * theta * sin_phi;
  Fsrp[2] = 0.0;

  Tsrp[0] = 0.25 * Ks_type * r_mag * r_mag * theta * sin_phi;
  Tsrp[1] = -0.25 * Ks_type * r_mag * r_mag * theta * cos_phi;
  Tsrp[2] = 0.0;

  // Damping forces, if relevant

  if (damping_style == DERIVATIVE) {
    double dt_inv = 1.0 / update->dt;

    // radial
    double **v = atom->v;
    double dv[3], rhat[3];
    MathExtra::sub3(v[i2], v[i1], dv);
    MathExtra::scale3(r_mag_inv, r, rhat);
    double dvdotr = MathExtra::dot3(dv, rhat);
    Fr += gr[type] * dvdotr;

    // shear
    double dgamma = wrap_2pi(gamma - gamma_prior) * dt_inv;
    Fs_mag += gs[type] * r_mag * dgamma;

    // twisting
    double dpsi = wrap_2pi(psi - psi_prior) * dt_inv;
    Ttp[2] += gt[type] * dpsi;

    // bending
    double dtheta = wrap_2pi(theta - theta_prior) * dt_inv;
    Tbp[0] += -gb[type] * dtheta * sin_phi;
    Tbp[1] += gb[type] * dtheta * cos_phi;
    Fsrp[0] += -0.5 * gs[type] * r_mag * dtheta * cos_phi;
    Fsrp[1] += -0.5 * gs[type] * r_mag * dtheta * sin_phi;
    Tsrp[0] += 0.25 * gs[type] * r_mag * r_mag * dtheta * sin_phi;
    Tsrp[1] += -0.25 * gs[type] * r_mag * r_mag * dtheta * cos_phi;

    bondstore[4] = gamma;
    bondstore[5] = theta;
    bondstore[6] = psi;
  }

  // Calculate net radial/shear forces (after damping)

  MathExtra::scale3(Fr * r_mag_inv, rb, Fn);
  MathExtra::scale3(Fr * r_mag_inv, rb, F_rot);
  MathExtra::scale3(Fs_mag, s, Fsrt);
  double Tst_mag = 0.5 * r_mag * Fs_mag;
  MathExtra::scale3(Tst_mag, t, Tst);

  // Rotate forces/torques back to 1st particle's frame

  MathExtra::quatrotvec(mq, Fsrp, Fsr);
  MathExtra::quatrotvec(mq, Tsrp, Tsr);
  MathExtra::quatrotvec(mq, Tbp, Tb);
  MathExtra::quatrotvec(mq, Ttp, Tt);

  // Sum forces and calculate magnitudes
  F_rot[0] = Fn[0] + Fsrt[0] + Fsr[0];
  F_rot[1] = Fn[1] + Fsrt[1] + Fsr[1];
  F_rot[2] = Fn[2] + Fsrt[2] + Fsr[2];
  MathExtra::quatrotvec(q2, F_rot, force1on2);

  T_rot[0] = Tst[0] + Tsr[0] + Tt[0] + Tb[0];
  T_rot[1] = Tst[1] + Tsr[1] + Tt[1] + Tb[1];
  T_rot[2] = Tst[2] + Tsr[2] + Tt[2] + Tb[2];
  MathExtra::quatrotvec(q2, T_rot, torque1on2);

  T_rot[0] = Tst[0] + Tsr[0] - Tt[0] - Tb[0];
  T_rot[1] = Tst[1] + Tsr[1] - Tt[1] - Tb[1];
  T_rot[2] = Tst[2] + Tsr[2] - Tt[2] - Tb[2];
  MathExtra::quatrotvec(q2, T_rot, torque2on1);

  MathExtra::add3(Fsrt, Fsr, Fs);
  Fs_mag = MathExtra::len3(Fs);
  Tt_mag = MathExtra::len3(Tt);
  Tb_mag = MathExtra::len3(Tb);

  // Includes damping contribution (unlike in average_frame_forces) for backwards compatibility
  breaking = Fr / Fcr[type] + Fs_mag / Fcs[type] + Tb_mag / Tcb[type] + Tt_mag / Tct[type];
  if (breaking < 0.0) breaking = 0.0;

  // Approximate bond energy:
  ebond = 0.5 * (Kr_type * (r_mag - r0_mag) * (r_mag - r0_mag) +
                 Ks_type * r_mag * r_mag * gamma * gamma +
                 Kt[type] * psi * psi +
                 Kb[type] * theta * theta);

  return breaking;

}

/* ----------------------------------------------------------------------
  Calculate damping using formulation in
        Y. Wang, F. Alonso-Marroquin, W. Guo 2015
  Note: n points towards 1 vs pointing towards 2
---------------------------------------------------------------------- */

void BondBPMRotational::dem_damping_forces(int i1, int i2, int type, double *r,
                                             double *force1on2, double *torque1on2, double *torque2on1)
{
  double v1dotr, v2dotr, w1dotr, w2dotr;
  double rhat[3], s1[3], s2[3], tdamp[3], tmp[3];
  double vn1[3], vn2[3], vt1[3], vt2[3], vroll[3];
  double wxn1[3], wxn2[3], wn1[3], wn2[3];

  // Note: 1 <-> 2 swapped from older 2025 version (before average frame and derivative damping)
  MathExtra::normalize3(r, rhat);

  double **v = atom->v;
  double **omega = atom->omega;

  // Damp normal velocity difference
  v1dotr = MathExtra::dot3(v[i1], rhat);
  v2dotr = MathExtra::dot3(v[i2], rhat);

  MathExtra::scale3(v1dotr, rhat, vn1);
  MathExtra::scale3(v2dotr, rhat, vn2);

  MathExtra::sub3(vn1, vn2, tmp);
  MathExtra::scale3(gr[type], tmp);
  MathExtra::add3(force1on2, tmp, force1on2);

  // Damp tangential objective velocities
  MathExtra::sub3(v[i1], vn1, vt1);
  MathExtra::sub3(v[i2], vn2, vt2);

  MathExtra::sub3(vt2, vt1, tmp);
  MathExtra::scale3(0.5, tmp);

  MathExtra::cross3(omega[i1], r, s1);
  MathExtra::scale3(0.5, s1);
  MathExtra::sub3(s1, tmp, s1);    // Eq 12

  MathExtra::cross3(omega[i2], r, s2);
  MathExtra::scale3(-0.5, s2);
  MathExtra::add3(s2, tmp, s2);    // Eq 13

  MathExtra::sub3(s1, s2, tmp);
  MathExtra::scale3(gs[type], tmp);
  MathExtra::add3(force1on2, tmp, force1on2);

  // Apply corresponding torque
  MathExtra::cross3(r, tmp, tdamp);
  MathExtra::scale3(-0.5, tdamp);
  MathExtra::add3(torque1on2, tdamp, torque1on2);
  MathExtra::add3(torque2on1, tdamp, torque2on1);

  // Damp rolling
  MathExtra::cross3(omega[i1], rhat, wxn1);
  MathExtra::cross3(omega[i2], rhat, wxn2);
  MathExtra::sub3(wxn1, wxn2, vroll);    // Eq. 31
  MathExtra::cross3(r, vroll, tdamp);

  MathExtra::scale3(0.5 * gb[type], tdamp);
  MathExtra::add3(torque1on2, tdamp, torque1on2);
  MathExtra::scale3(-1.0, tdamp);
  MathExtra::add3(torque2on1, tdamp, torque2on1);

  // Damp twist
  w1dotr = MathExtra::dot3(omega[i1], rhat);
  w2dotr = MathExtra::dot3(omega[i2], rhat);

  MathExtra::scale3(w1dotr, rhat, wn1);
  MathExtra::scale3(w2dotr, rhat, wn2);

  MathExtra::sub3(wn1, wn2, tdamp);    // Eq. 38
  MathExtra::scale3(0.5 * gt[type], tdamp);
  MathExtra::add3(torque1on2, tdamp, torque1on2);
  MathExtra::scale3(-1.0, tdamp);
  MathExtra::add3(torque2on1, tdamp, torque2on1);
}

/* ---------------------------------------------------------------------- */

void BondBPMRotational::compute(int eflag, int vflag)
{
  pre_compute();

  int i1, i2, itmp, n, type;
  double ri_norm, ebond, breaking, smooth;
  double rf[3], ri[3], force1on2[3], torque1on2[3], torque2on1[3];

  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  tagint *tag = atom->tag;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  double **bondstore = fix_bond_history->bondstore;

  for (n = 0; n < nbondlist; n++) {
    if (bondlist[n][2] <= 0) continue;

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    ri_norm = bondstore[n][0];

    // Order so i1 has lower tag
    if (tag[i2] < tag[i1]) {
      itmp = i1;
      i1 = i2;
      i2 = itmp;
    }

    if (ri_norm < EPSILON || std::isnan(ri_norm)) ri_norm = store_bond(n, i1, i2);

    // ri points from i1 (lower tag) to i2 (higher tag)
    ri[0] = bondstore[n][1] * ri_norm;
    ri[1] = bondstore[n][2] * ri_norm;
    ri[2] = bondstore[n][3] * ri_norm;

    // rf = x[i2] - x[i1]
    MathExtra::sub3(x[i2], x[i1], rf);

    if (frame_style == AVERAGE) {
      breaking = average_frame_forces(i1, i2, type, ri, rf,
                                  force1on2, torque1on2, torque2on1, ebond, bondstore[n]);
    } else {
      breaking = particle_frame_forces(i1, i2, type, ri, rf,
                                force1on2, torque1on2, torque2on1, ebond, bondstore[n]);
    }

    if (damping_style == DEM)
      dem_damping_forces(i1, i2, type, rf, force1on2, torque1on2, torque2on1);

    if ((breaking >= 1.0) && break_flag) {
      bondlist[n][2] = 0;
      process_broken(i1, i2);
      continue;
    }

    if (smooth_flag) {
      smooth = breaking * breaking;
      smooth = 1.0 - smooth * smooth;
    } else {
      smooth = 1.0;
    }

    MathExtra::scale3(smooth, force1on2);

    if (newton_bond || i1 < nlocal) {
      f[i1][0] -= force1on2[0];
      f[i1][1] -= force1on2[1];
      f[i1][2] -= force1on2[2];

      MathExtra::scale3(smooth, torque2on1);
      torque[i1][0] += torque2on1[0];
      torque[i1][1] += torque2on1[1];
      torque[i1][2] += torque2on1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += force1on2[0];
      f[i2][1] += force1on2[1];
      f[i2][2] += force1on2[2];

      MathExtra::scale3(smooth, torque1on2);
      torque[i2][0] += torque1on2[0];
      torque[i2][1] += torque1on2[1];
      torque[i2][2] += torque1on2[2];
    }

    if (evflag)
      ev_tally_xyz(i1, i2, nlocal, newton_bond, ebond, -force1on2[0], -force1on2[1],
                   -force1on2[2], rf[0], rf[1], rf[2]);
  }

  post_compute();
}

/* ---------------------------------------------------------------------- */

void BondBPMRotational::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(Kr, np1, "bond:Kr");
  memory->create(Ks, np1, "bond:Ks");
  memory->create(Kt, np1, "bond:Kt");
  memory->create(Kb, np1, "bond:Kb");
  memory->create(Fcr, np1, "bond:Fcr");
  memory->create(Fcs, np1, "bond:Fcs");
  memory->create(Tct, np1, "bond:Tct");
  memory->create(Tcb, np1, "bond:Tcb");
  memory->create(gr, np1, "bond:gr");
  memory->create(gs, np1, "bond:gs");
  memory->create(gt, np1, "bond:gt");
  memory->create(gb, np1, "bond:gb");

  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
   args: Kr Ks Kt Kb Fcr Fcs Tct Tcb gr gs gt gb
------------------------------------------------------------------------- */

void BondBPMRotational::coeff(int narg, char **arg)
{
  if (narg != 13) error->all(FLERR, "Incorrect args for bond coefficients: "
                            "expected Kr Ks Kt Kb Fcr Fcs Tct Tcb gr gs gt gb");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double Kr_one = utils::numeric(FLERR, arg[1], false, lmp);
  double Ks_one = utils::numeric(FLERR, arg[2], false, lmp);
  double Kt_one = utils::numeric(FLERR, arg[3], false, lmp);
  double Kb_one = utils::numeric(FLERR, arg[4], false, lmp);
  double Fcr_one = utils::numeric(FLERR, arg[5], false, lmp);
  double Fcs_one = utils::numeric(FLERR, arg[6], false, lmp);
  double Tct_one = utils::numeric(FLERR, arg[7], false, lmp);
  double Tcb_one = utils::numeric(FLERR, arg[8], false, lmp);
  double gr_one = utils::numeric(FLERR, arg[9], false, lmp);
  double gs_one = utils::numeric(FLERR, arg[10], false, lmp);
  double gt_one = utils::numeric(FLERR, arg[11], false, lmp);
  double gb_one = utils::numeric(FLERR, arg[12], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    Kr[i] = Kr_one;
    Ks[i] = Ks_one;
    Kt[i] = Kt_one;
    Kb[i] = Kb_one;
    Fcr[i] = Fcr_one;
    Fcs[i] = Fcs_one;
    Tct[i] = Tct_one;
    Tcb[i] = Tcb_one;
    gr[i] = gr_one;
    gs[i] = gs_one;
    gt[i] = gt_one;
    gb[i] = gb_one;
    setflag[i] = 1;
    count++;

    if (Fcr[i] / Kr[i] > max_stretch) max_stretch = Fcr[i] / Kr[i];
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   check for correct settings and create fix
------------------------------------------------------------------------- */

void BondBPMRotational::init_style()
{
  // History: [0]=ri_mag, [1-3]=ri_hat
  //  if derivative damping: [4]=gamma, [5]=theta, [6]=psi
  if (damping_style == DEM)
    nhistory = 4;
  else
    nhistory = 7;

  BondBPM::init_style();

  if (!atom->quat_flag || !atom->radius_flag || !atom->omega_flag)
    error->all(FLERR, "Bond bpm/rotational requires atom style bpm/sphere");
  if (comm->ghost_velocity == 0)
    error->all(FLERR, "Bond bpm/rotational requires ghost atoms store velocity");

  if (domain->dimension == 2)
    error->warning(FLERR, "Bond style bpm/rotational not intended for 2d use");
}

/* ---------------------------------------------------------------------- */

void BondBPMRotational::settings(int narg, char **arg)
{
  BondBPM::settings(narg, arg);

  int iarg;
  for (std::size_t i = 0; i < leftover_iarg.size(); i++) {
    iarg = leftover_iarg[i];
    if (strcmp(arg[iarg], "smooth") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "bond_style bpm/rotational smooth", error);
      smooth_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      i += 1;
    } else if (strcmp(arg[iarg], "normalize") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "bond_style bpm/rotational normalize", error);
      normalize_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      i += 1;
    } else if (strcmp(arg[iarg], "frame") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "bond_style bpm/rotational frame", error);
      if (strcmp(arg[iarg + 1], "particle") == 0)
        frame_style = PARTICLE;
      else if (strcmp(arg[iarg + 1], "average") == 0)
        frame_style = AVERAGE;
      else
        error->all(FLERR, "Unknown frame style {}", arg[iarg + 1]);
      i += 1;
    } else if (strcmp(arg[iarg], "damping") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "bond_style bpm/rotational damping", error);
      if (strcmp(arg[iarg + 1], "dem") == 0)
        damping_style = DEM;
      else if (strcmp(arg[iarg + 1], "derivative") == 0)
        damping_style = DERIVATIVE;
      else
        error->all(FLERR, "Unknown damping style {}", arg[iarg + 1]);
      i += 1;
    } else {
      error->all(FLERR, "Illegal bond bpm command, invalid argument {}", arg[iarg]);
    }
  }

  if (smooth_flag && !break_flag)
    error->all(FLERR, "Illegal bond bpm command, must turn off smoothing with break no option");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondBPMRotational::write_restart(FILE *fp)
{
  BondBPM::write_restart(fp);
  write_restart_settings(fp);

  fwrite(&Kr[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&Ks[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&Kt[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&Kb[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&Fcr[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&Fcs[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&Tct[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&Tcb[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gr[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gs[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gt[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gb[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondBPMRotational::read_restart(FILE *fp)
{
  BondBPM::read_restart(fp);
  read_restart_settings(fp);
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &Kr[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &Ks[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &Kt[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &Kb[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &Fcr[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &Fcs[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &Tct[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &Tcb[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gr[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gs[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gt[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gb[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&Kr[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Ks[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Kt[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Kb[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Fcr[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Fcs[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Tct[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Tcb[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gr[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gs[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gt[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gb[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void BondBPMRotational::write_restart_settings(FILE *fp)
{
  fwrite(&smooth_flag, sizeof(int), 1, fp);
  fwrite(&normalize_flag, sizeof(int), 1, fp);
  fwrite(&damping_style, sizeof(int), 1, fp);
  fwrite(&frame_style, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void BondBPMRotational::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &smooth_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &normalize_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &damping_style, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &frame_style, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&smooth_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&normalize_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&damping_style, 1, MPI_INT, 0, world);
  MPI_Bcast(&frame_style, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

double BondBPMRotational::single(int type, double rsq, int i, int j, double &fforce)
{
  if (type <= 0) return 0.0;

  int flipped = 0;
  if (atom->tag[j] < atom->tag[i]) {
    int itmp = i;
    i = j;
    j = itmp;
    flipped = 1;
  }

  double ri_norm = 0.0;
  double ri[3], rf[3], bondstore[7];
  for (int n = 0; n < atom->num_bond[i]; n++) {
    if (atom->bond_atom[i][n] == atom->tag[j]) {
      ri_norm = fix_bond_history->get_atom_value(i, n, 0);
      ri[0] = fix_bond_history->get_atom_value(i, n, 1) * ri_norm;
      ri[1] = fix_bond_history->get_atom_value(i, n, 2) * ri_norm;
      ri[2] = fix_bond_history->get_atom_value(i, n, 3) * ri_norm;
      if (damping_style == DERIVATIVE) {
        bondstore[4] = fix_bond_history->get_atom_value(i, n, 4);
        bondstore[5] = fix_bond_history->get_atom_value(i, n, 5);
        bondstore[6] = fix_bond_history->get_atom_value(i, n, 6);
      }
    }
  }

  double **x = atom->x;
  MathExtra::sub3(x[j], x[i], rf);

  double force1on2[3], torque1on2[3], torque2on1[3], ebond;

  // Note, derivative damping will update bondstore, but this doesn't affect real values
  double breaking;
  if (frame_style == AVERAGE) {
    breaking = average_frame_forces(i, j, type, ri, rf,
                                force1on2, torque1on2, torque2on1, ebond, bondstore);
  } else {
    breaking = particle_frame_forces(i, j, type, ri, rf,
                              force1on2, torque1on2, torque2on1, ebond, bondstore);
  }

  if (damping_style == DEM)
    dem_damping_forces(i, j, type, rf, force1on2, torque1on2, torque2on1);

  double rf_hat[3];
  double rf_mag = sqrt(rsq);
  MathExtra::scale3(1.0 / rf_mag, rf, rf_hat);
  fforce = MathExtra::dot3(force1on2, rf_hat);

  double smooth = 1.0;
  if (smooth_flag) {
    smooth = breaking * breaking;
    smooth = 1.0 - smooth * smooth;
    fforce *= smooth;
  }

  MathExtra::scale3(smooth, force1on2);
  svector[0] = ri_norm;
  if (flipped) {
    svector[1] = -ri[0];
    svector[2] = -ri[1];
    svector[3] = -ri[2];
    svector[4] = -force1on2[0];
    svector[5] = -force1on2[1];
    svector[6] = -force1on2[2];
  } else {
    svector[1] = ri[0];
    svector[2] = ri[1];
    svector[3] = ri[2];
    svector[4] = force1on2[0];
    svector[5] = force1on2[1];
    svector[6] = force1on2[2];
  }

  return ebond;
}
