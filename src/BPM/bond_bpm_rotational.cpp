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
#include "modify.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

static constexpr double EPSILON = 1e-10;

using namespace LAMMPS_NS;
using MathConst::MY_SQRT2;
/* ---------------------------------------------------------------------- */

static double acos_limit(double c)
{
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;
  return acos(c);
}

/* ---------------------------------------------------------------------- */

BondBPMRotational::BondBPMRotational(LAMMPS *_lmp) :
    BondBPM(_lmp), Kr(nullptr), Ks(nullptr), Kt(nullptr), Kb(nullptr), gnorm(nullptr),
    gslide(nullptr), groll(nullptr), gtwist(nullptr), Fcr(nullptr), Fcs(nullptr), Tct(nullptr),
    Tcb(nullptr)
{
  partial_flag = 1;
  smooth_flag = 1;
  normalize_flag = 0;

  nhistory = 4;

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
    memory->destroy(gnorm);
    memory->destroy(gslide);
    memory->destroy(groll);
    memory->destroy(gtwist);
  }
}

/* ----------------------------------------------------------------------
  Store data for a single bond - if bond added after LAMMPS init (e.g. pour)
------------------------------------------------------------------------- */

double BondBPMRotational::store_bond(int n, int i, int j)
{
  double delx, dely, delz, r, rinv;
  double **x = atom->x;
  tagint *tag = atom->tag;
  double **bondstore = fix_bond_history->bondstore;

  if (tag[i] < tag[j]) {
    delx = x[i][0] - x[j][0];
    dely = x[i][1] - x[j][1];
    delz = x[i][2] - x[j][2];
  } else {
    delx = x[j][0] - x[i][0];
    dely = x[j][1] - x[i][1];
    delz = x[j][2] - x[i][2];
  }

  r = sqrt(delx * delx + dely * dely + delz * delz);
  rinv = 1.0 / r;

  bondstore[n][0] = r;
  bondstore[n][1] = delx * rinv;
  bondstore[n][2] = dely * rinv;
  bondstore[n][3] = delz * rinv;

  if (i < atom->nlocal) {
    for (int m = 0; m < atom->num_bond[i]; m++) {
      if (atom->bond_atom[i][m] == tag[j]) {
        fix_bond_history->update_atom_value(i, m, 0, r);
        fix_bond_history->update_atom_value(i, m, 1, delx * rinv);
        fix_bond_history->update_atom_value(i, m, 2, dely * rinv);
        fix_bond_history->update_atom_value(i, m, 3, delz * rinv);
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

      //Skip if bond was turned off
      if (type < 0) continue;

      // map to find index n for tag
      j = atom->map(atom->bond_atom[i][m]);
      if (j == -1) error->one(FLERR, "Atom missing in BPM bond");

      // Save orientation as pointing towards small tag
      if (tag[i] < tag[j]) {
        delx = x[i][0] - x[j][0];
        dely = x[i][1] - x[j][1];
        delz = x[i][2] - x[j][2];
      } else {
        delx = x[j][0] - x[i][0];
        dely = x[j][1] - x[i][1];
        delz = x[j][2] - x[i][2];
      }

      // Get closest image in case bonded with ghost
      domain->minimum_image(delx, dely, delz);
      r = sqrt(delx * delx + dely * dely + delz * delz);
      rinv = 1.0 / r;

      fix_bond_history->update_atom_value(i, m, 0, r);
      fix_bond_history->update_atom_value(i, m, 1, delx * rinv);
      fix_bond_history->update_atom_value(i, m, 2, dely * rinv);
      fix_bond_history->update_atom_value(i, m, 3, delz * rinv);
    }
  }

  fix_bond_history->post_neighbor();
}

/* ----------------------------------------------------------------------
  Calculate forces using formulation in:
    1) Y. Wang Acta Geotechnica 2009
    2) P. Mora & Y. Wang Advances in Geomcomputing 2009
---------------------------------------------------------------------- */

double BondBPMRotational::elastic_forces(int i1, int i2, int type, double r_mag, double r0_mag,
                                         double r_mag_inv, double * /*rhat*/, double *r, double *r0,
                                         double *force1on2, double *torque1on2, double *torque2on1)
{
  double breaking, temp, r0_dot_rb, c, gamma;
  double psi, theta, cos_phi, sin_phi;
  double mag_in_plane, mag_out_plane;
  double Fs_mag, Tt_mag, Tb_mag;

  double q1[4], q2[4];
  double q2inv[4], mq[4], mqinv[4], qp21[4], q21[4], qtmp[4];
  double rb[3], rb_x_r0[3], s[3], t[3];
  double Fr, Fs[3], Fsp[3], F_rot[3], Ftmp[3];
  double Ts[3], Tb[3], Tt[3], Tbp[3], Ttp[3], Tsp[3], T_rot[3], Ttmp[3];

  double **quat = atom->quat;
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
  Fr = Kr_type * (r_mag - r0_mag);

  MathExtra::scale3(Fr * r_mag_inv, rb, F_rot);

  // Calculate forces due to tangential displacements (no rotation)
  r0_dot_rb = MathExtra::dot3(r0, rb);
  c = r0_dot_rb * r_mag_inv * r0_mag_inv;
  gamma = acos_limit(c);

  MathExtra::cross3(rb, r0, rb_x_r0);
  MathExtra::cross3(rb, rb_x_r0, s);
  MathExtra::norm3(s);

  MathExtra::scale3(Ks_type * r_mag * gamma, s, Fs);

  // Calculate torque due to tangential displacements
  MathExtra::cross3(r0, rb, t);
  MathExtra::norm3(t);

  MathExtra::scale3(0.5 * r_mag * Ks_type * r_mag * gamma, t, Ts);

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
  theta = acos_limit(c);

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

  Fsp[0] = -0.5 * Ks_type * r_mag * theta * cos_phi;
  Fsp[1] = -0.5 * Ks_type * r_mag * theta * sin_phi;
  Fsp[2] = 0.0;

  Tsp[0] = 0.25 * Ks_type * r_mag * r_mag * theta * sin_phi;
  Tsp[1] = -0.25 * Ks_type * r_mag * r_mag * theta * cos_phi;
  Tsp[2] = 0.0;

  // Rotate forces/torques back to 1st particle's frame

  MathExtra::quatrotvec(mq, Fsp, Ftmp);
  MathExtra::quatrotvec(mq, Tsp, Ttmp);
  for (int m = 0; m < 3; m++) {
    Fs[m] += Ftmp[m];
    Ts[m] += Ttmp[m];
  }

  MathExtra::quatrotvec(mq, Tbp, Tb);
  MathExtra::quatrotvec(mq, Ttp, Tt);

  // Sum forces and calculate magnitudes
  F_rot[0] += Fs[0];
  F_rot[1] += Fs[1];
  F_rot[2] += Fs[2];
  MathExtra::quatrotvec(q2, F_rot, force1on2);

  T_rot[0] = Ts[0] + Tt[0] + Tb[0];
  T_rot[1] = Ts[1] + Tt[1] + Tb[1];
  T_rot[2] = Ts[2] + Tt[2] + Tb[2];
  MathExtra::quatrotvec(q2, T_rot, torque1on2);

  T_rot[0] = Ts[0] - Tt[0] - Tb[0];
  T_rot[1] = Ts[1] - Tt[1] - Tb[1];
  T_rot[2] = Ts[2] - Tt[2] - Tb[2];
  MathExtra::quatrotvec(q2, T_rot, torque2on1);

  Fs_mag = MathExtra::len3(Fs);
  Tt_mag = MathExtra::len3(Tt);
  Tb_mag = MathExtra::len3(Tb);

  breaking = Fr / Fcr[type] + Fs_mag / Fcs[type] + Tb_mag / Tcb[type] + Tt_mag / Tct[type];
  if (breaking < 0.0) breaking = 0.0;

  return breaking;
}

/* ----------------------------------------------------------------------
  Calculate damping using formulation in
        Y. Wang, F. Alonso-Marroquin, W. Guo 2015
  Note: n points towards 1 vs pointing towards 2
---------------------------------------------------------------------- */

void BondBPMRotational::damping_forces(int i1, int i2, int type, double *rhat, double *r,
                                       double *force1on2, double *torque1on2, double *torque2on1)
{
  double v1dotr, v2dotr, w1dotr, w2dotr;
  double s1[3], s2[3], tdamp[3], tmp[3];
  double vn1[3], vn2[3], vt1[3], vt2[3], vroll[3];
  double wxn1[3], wxn2[3], wn1[3], wn2[3];

  double **v = atom->v;
  double **omega = atom->omega;

  // Damp normal velocity difference
  v1dotr = MathExtra::dot3(v[i1], rhat);
  v2dotr = MathExtra::dot3(v[i2], rhat);

  MathExtra::scale3(v1dotr, rhat, vn1);
  MathExtra::scale3(v2dotr, rhat, vn2);

  MathExtra::sub3(vn1, vn2, tmp);
  MathExtra::scale3(gnorm[type], tmp);
  MathExtra::add3(force1on2, tmp, force1on2);

  // Damp tangential objective velocities
  MathExtra::sub3(v[i1], vn1, vt1);
  MathExtra::sub3(v[i2], vn2, vt2);

  MathExtra::sub3(vt2, vt1, tmp);
  MathExtra::scale3(0.5, tmp);

  MathExtra::cross3(omega[i1], r, s1);
  MathExtra::scale3(-0.5, s1);
  MathExtra::sub3(s1, tmp, s1);    // Eq 12

  MathExtra::cross3(omega[i2], r, s2);
  MathExtra::scale3(0.5, s2);
  MathExtra::add3(s2, tmp, s2);    // Eq 13

  MathExtra::sub3(s1, s2, tmp);
  MathExtra::scale3(gslide[type], tmp);
  MathExtra::add3(force1on2, tmp, force1on2);

  // Apply corresponding torque
  MathExtra::cross3(r, tmp, tdamp);
  MathExtra::scale3(0.5, tdamp);
  MathExtra::add3(torque1on2, tdamp, torque1on2);
  MathExtra::add3(torque2on1, tdamp, torque2on1);

  // Damp rolling
  MathExtra::cross3(omega[i1], rhat, wxn1);
  MathExtra::cross3(omega[i2], rhat, wxn2);
  MathExtra::sub3(wxn1, wxn2, vroll);    // Eq. 31
  MathExtra::cross3(r, vroll, tdamp);

  MathExtra::scale3(0.5 * groll[type], tdamp);
  MathExtra::add3(torque1on2, tdamp, torque1on2);
  MathExtra::scale3(-1.0, tdamp);
  MathExtra::add3(torque2on1, tdamp, torque2on1);

  // Damp twist
  w1dotr = MathExtra::dot3(omega[i1], rhat);
  w2dotr = MathExtra::dot3(omega[i2], rhat);

  MathExtra::scale3(w1dotr, rhat, wn1);
  MathExtra::scale3(w2dotr, rhat, wn2);

  MathExtra::sub3(wn1, wn2, tdamp);    // Eq. 38
  MathExtra::scale3(0.5 * gtwist[type], tdamp);
  MathExtra::add3(torque1on2, tdamp, torque1on2);
  MathExtra::scale3(-1.0, tdamp);
  MathExtra::add3(torque2on1, tdamp, torque2on1);
}

/* ---------------------------------------------------------------------- */

void BondBPMRotational::compute(int eflag, int vflag)
{

  if (!fix_bond_history->stored_flag) {
    fix_bond_history->stored_flag = true;
    store_data();
  }

  if (hybrid_flag)
    fix_bond_history->compress_history();

  int i1, i2, itmp, n, type;
  double r[3], r0[3], rhat[3];
  double rsq, r0_mag, r_mag, r_mag_inv;
  double breaking, smooth;
  double force1on2[3], torque1on2[3], torque2on1[3];

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

    // skip bond if already broken
    if (bondlist[n][2] <= 0) continue;

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    r0_mag = bondstore[n][0];

    // Ensure pair is always ordered such that r0 points in
    // a consistent direction and to ensure numerical operations
    // are identical to minimize the possibility that a bond straddling
    // an mpi grid (newton off) doesn't break on one proc but not the other
    if (tag[i2] < tag[i1]) {
      itmp = i1;
      i1 = i2;
      i2 = itmp;
    }

    // If bond hasn't been set - should be initialized to zero
    if (r0_mag < EPSILON || std::isnan(r0_mag)) r0_mag = store_bond(n, i1, i2);

    r0[0] = bondstore[n][1];
    r0[1] = bondstore[n][2];
    r0[2] = bondstore[n][3];
    MathExtra::scale3(r0_mag, r0);

    // Note this is the reverse of Mora & Wang
    MathExtra::sub3(x[i1], x[i2], r);

    rsq = MathExtra::lensq3(r);
    r_mag = sqrt(rsq);
    r_mag_inv = 1.0 / r_mag;
    MathExtra::scale3(r_mag_inv, r, rhat);

    // ------------------------------------------------------//
    //  Calculate forces, check if bond breaks
    // ------------------------------------------------------//

    breaking = elastic_forces(i1, i2, type, r_mag, r0_mag, r_mag_inv, rhat, r, r0, force1on2,
                              torque1on2, torque2on1);

    if (breaking >= 1.0) {
      bondlist[n][2] = 0;
      process_broken(i1, i2);
      continue;
    }

    damping_forces(i1, i2, type, rhat, r, force1on2, torque1on2, torque2on1);

    if (smooth_flag) {
      smooth = breaking * breaking;
      smooth = 1.0 - smooth * smooth;
    } else {
      smooth = 1.0;
    }

    // ------------------------------------------------------//
    //  Apply forces and torques to particles
    // ------------------------------------------------------//

    if (newton_bond || i1 < nlocal) {
      f[i1][0] -= force1on2[0] * smooth;
      f[i1][1] -= force1on2[1] * smooth;
      f[i1][2] -= force1on2[2] * smooth;

      torque[i1][0] += torque2on1[0] * smooth;
      torque[i1][1] += torque2on1[1] * smooth;
      torque[i1][2] += torque2on1[2] * smooth;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += force1on2[0] * smooth;
      f[i2][1] += force1on2[1] * smooth;
      f[i2][2] += force1on2[2] * smooth;

      torque[i2][0] += torque1on2[0] * smooth;
      torque[i2][1] += torque1on2[1] * smooth;
      torque[i2][2] += torque1on2[2] * smooth;
    }

    if (evflag)
      ev_tally_xyz(i1, i2, nlocal, newton_bond, 0.0, -force1on2[0] * smooth, -force1on2[1] * smooth,
                   -force1on2[2] * smooth, r[0], r[1], r[2]);
  }

  if (hybrid_flag)
    fix_bond_history->uncompress_history();
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
  memory->create(gnorm, np1, "bond:gnorm");
  memory->create(gslide, np1, "bond:gslide");
  memory->create(groll, np1, "bond:groll");
  memory->create(gtwist, np1, "bond:gtwist");

  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondBPMRotational::coeff(int narg, char **arg)
{
  if (narg != 13) error->all(FLERR, "Incorrect args for bond coefficients");
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
  double gnorm_one = utils::numeric(FLERR, arg[9], false, lmp);
  double gslide_one = utils::numeric(FLERR, arg[10], false, lmp);
  double groll_one = utils::numeric(FLERR, arg[11], false, lmp);
  double gtwist_one = utils::numeric(FLERR, arg[12], false, lmp);

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
    gnorm[i] = gnorm_one;
    gslide[i] = gslide_one;
    groll[i] = groll_one;
    gtwist[i] = gtwist_one;
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
      if (iarg + 1 > narg) error->all(FLERR, "Illegal bond bpm command, missing option for smooth");
      smooth_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      i += 1;
    } else if (strcmp(arg[iarg], "normalize") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal bond bpm command, missing option for normalize");
      normalize_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      i += 1;
    } else {
      error->all(FLERR, "Illegal bond bpm command, invalid argument {}", arg[iarg]);
    }
  }
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
  fwrite(&gnorm[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gslide[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&groll[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&gtwist[1], sizeof(double), atom->nbondtypes, fp);
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
    utils::sfread(FLERR, &gnorm[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gslide[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &groll[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &gtwist[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&Kr[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Ks[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Kt[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Kb[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Fcr[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Fcs[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Tct[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&Tcb[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gnorm[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gslide[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&groll[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gtwist[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void BondBPMRotational::write_restart_settings(FILE *fp)
{
  fwrite(&smooth_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void BondBPMRotational::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) utils::sfread(FLERR, &smooth_flag, sizeof(int), 1, fp, nullptr, error);
  MPI_Bcast(&smooth_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

double BondBPMRotational::single(int type, double rsq, int i, int j, double &fforce)
{
  // Not yet enabled
  if (type <= 0) return 0.0;

  int flipped = 0;
  if (atom->tag[j] < atom->tag[i]) {
    int itmp = i;
    i = j;
    j = itmp;
    flipped = 1;
  }

  double r0_mag, r_mag, r_mag_inv;
  double r0[3], r[3], rhat[3];
  for (int n = 0; n < atom->num_bond[i]; n++) {
    if (atom->bond_atom[i][n] == atom->tag[j]) {
      r0_mag = fix_bond_history->get_atom_value(i, n, 0);
      r0[0] = fix_bond_history->get_atom_value(i, n, 1);
      r0[1] = fix_bond_history->get_atom_value(i, n, 2);
      r0[2] = fix_bond_history->get_atom_value(i, n, 3);
    }
  }

  double **x = atom->x;
  MathExtra::scale3(r0_mag, r0);
  MathExtra::sub3(x[i], x[j], r);

  r_mag = sqrt(rsq);
  r_mag_inv = 1.0 / r_mag;
  MathExtra::scale3(r_mag_inv, r, rhat);

  double force1on2[3], torque1on2[3], torque2on1[3];
  double breaking = elastic_forces(i, j, type, r_mag, r0_mag, r_mag_inv, rhat, r, r0, force1on2,
                                   torque1on2, torque2on1);
  damping_forces(i, j, type, rhat, r, force1on2, torque1on2, torque2on1);
  fforce = MathExtra::dot3(force1on2, rhat);
  fforce *= -1;

  double smooth = 1.0;
  if (smooth_flag) {
    smooth = breaking * breaking;
    smooth = 1.0 - smooth * smooth;
    fforce *= smooth;
  }

  // set single_extra quantities

  svector[0] = r0_mag;
  if (flipped) {
    svector[1] = -r0[0];
    svector[2] = -r0[1];
    svector[3] = -r0[2];
    svector[4] = force1on2[0] * smooth;
    svector[5] = force1on2[1] * smooth;
    svector[6] = force1on2[2] * smooth;
  } else {
    svector[1] = r0[0];
    svector[2] = r0[1];
    svector[3] = r0[2];
    svector[4] = -force1on2[0] * smooth;
    svector[5] = -force1on2[1] * smooth;
    svector[6] = -force1on2[2] * smooth;
  }

  return 0.0;
}
