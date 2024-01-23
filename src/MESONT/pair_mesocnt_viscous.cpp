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
   Contributing author: Philipp Kloza (University of Cambridge)
                        pak37@cam.ac.uk
------------------------------------------------------------------------- */

#include "pair_mesocnt_viscous.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathExtra;
using MathConst::MY_PI;

#define SELF_CUTOFF 3
#define RHOMIN 10.0

#define QUAD_FINF 129
#define QUAD_FSEMI 10

/* ---------------------------------------------------------------------- */

void PairMesoCNTViscous::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  int i, j, k, i1, i2, j1, j2;
  int endflag, endindex;
  int clen, numchain;
  int *end, *nchain;
  int **chain;
  double fend, lp, scale, sumw, sumw_inv;
  double evdwl, evdwl_chain;
  double vtot, fvisc_tot;
  double *r1, *r2, *q1, *q2, *qe;
  double *vr1, *vr2, *vq1, *vq2;
  double vr[3], vp1[3], vp2[3], vp[3], vrel[3], fvisc[3];
  double ftotal[3], ftorque[3], torque[3], delr1[3], delr2[3], delqe[3];
  double t1[3], t2[3], t3[3];
  double dr1_sumw[3], dr2_sumw[3];
  double dr1_w[3], dr2_w[3], dq1_w[3], dq2_w[3];
  double fgrad_r1_p1[3], fgrad_r1_p2[3], fgrad_r2_p1[3], fgrad_r2_p2[3];
  double fgrad_q_p1[3], fgrad_q_p2[3];
  double q1_dr1_w[3][3], q1_dr2_w[3][3], q2_dr1_w[3][3], q2_dr2_w[3][3];
  double dr1_p1[3][3], dr1_p2[3][3], dr2_p1[3][3], dr2_p2[3][3];
  double dq_p1[3][3], dq_p2[3][3];
  double temp[3][3];

  evdwl = 0.0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  int flag, cols;
  int buckled_index = atom->find_custom("buckled", flag, cols);

  // update bond neighbor list when necessary

  if (update->ntimestep == neighbor->lastcall) {
    if (neigh_flag)
      bond_neigh_topo();
    else
      bond_neigh_id();
  }

  // iterate over all bonds

  for (i = 0; i < nbondlist; i++) {
    i1 = bondlist[i][0];
    i2 = bondlist[i][1];

    r1 = x[i1];
    r2 = x[i2];

    vr1 = v[i1];
    vr2 = v[i2];
    add3(vr1, vr2, vr);
    scale3(0.5, vr);

    numchain = numchainlist[i];
    end = endlist[i];
    nchain = nchainlist[i];
    chain = chainlist[i];

    // iterate over all neighbouring chains

    for (j = 0; j < numchain; j++) {
      clen = nchain[j];
      if (clen < 2) continue;

      // check if segment-segment interactions are necessary

      endflag = end[j];
      int buckled = 0;
      if (buckled_index != -1) {
        for (k = 0; k < clen; k++) {
          if (atom->ivector[buckled_index][chain[j][k]]) {
            buckled = 1;
            break;
          }
        }
      }

      if (segment_flag || buckled || j == selfid[i] || endflag == 3) {

        zero3(vp1);
        zero3(vp2);
        sumw = 0.0;

        for (k = 0; k < clen - 1; k++) {

          // segment-segment interaction

          // exclude SELF_CUTOFF neighbors in self-chain

          if (j == selfid[i]) {
            int min11 = abs(k - selfpos[i][0]);
            int min12 = abs(k - selfpos[i][1]);
            int min21 = abs(k + 1 - selfpos[i][0]);
            int min22 = abs(k + 1 - selfpos[i][1]);
            int min = min11;
            if (min12 < min) min = min12;
            if (min21 < min) min = min21;
            if (min22 < min) min = min22;

            if (min < SELF_CUTOFF) continue;
          }

          j1 = chain[j][k];
          j2 = chain[j][k + 1];
          j1 &= NEIGHMASK;
          j2 &= NEIGHMASK;
          q1 = x[j1];
          q2 = x[j2];
          vq1 = v[j1];
          vq2 = v[j2];

          // weighted velocity for friction

          w[k] = weight(r1, r2, q1, q2);
          if (w[k] == 0.0) continue;
          sumw += w[k];
          scaleadd3(w[k], vq1, vp1, vp1);
          scaleadd3(w[k], vq2, vp2, vp2);

          // check if orientation of segment needs to be flipped to prevent overlap

          geometry(r1, r2, q1, q2, q1, p, m, param, basis);
          if (param[0] > cutoff) continue;

          double calpha = cos(param[1]);
          double salpha = sin(param[1]);
          double hsq = param[0] * param[0];

          double ceta = calpha * param[4];
          double seta = salpha * param[4];
          double dsq1 = hsq + seta * seta;
          if (ceta < param[2]) {
            double dceta = ceta - param[2];
            dsq1 += dceta * dceta;
          } else if (ceta > param[3]) {
            double dceta = ceta - param[3];
            dsq1 += dceta * dceta;
          }

          ceta = calpha * param[5];
          seta = salpha * param[5];

          double dsq2 = hsq + seta * seta;
          if (ceta < param[2]) {
            double dceta = ceta - param[2];
            dsq2 += dceta * dceta;
          } else if (ceta > param[3]) {
            double dceta = ceta - param[3];
            dsq2 += dceta * dceta;
          }

          if (dsq1 > cutoffsq && dsq2 > cutoffsq) continue;

          int jj1, jj2;

          // do flip if necessary

          if (dsq1 < dsq2) {
            jj1 = j1;
            jj2 = j2;

          } else {
            if (param[1] > MY_PI)
              param[1] -= MY_PI;
            else
              param[1] += MY_PI;

            double eta2 = -param[5];
            param[5] = -param[4];
            param[4] = eta2;
            param[6] = eta2;

            negate3(m);

            jj1 = j2;
            jj2 = j1;
          }

          // first force contribution

          fsemi(param, evdwl, fend, flocal);

          if (evdwl == 0.0) continue;

          // transform to global coordinate system

          matvec(basis[0], basis[1], basis[2], flocal[0], fglobal[0]);
          matvec(basis[0], basis[1], basis[2], flocal[1], fglobal[1]);

          // forces acting on segment 2

          add3(fglobal[0], fglobal[1], ftotal);
          scaleadd3(fend, m, ftotal, ftotal);
          scale3(-0.5, ftotal);

          sub3(r1, p, delr1);
          sub3(r2, p, delr2);
          cross3(delr1, fglobal[0], t1);
          cross3(delr2, fglobal[1], t2);
          add3(t1, t2, torque);

          cross3(torque, m, ftorque);
          lp = param[5] - param[4];
          scale3(1.0 / lp, ftorque);

          add3(ftotal, ftorque, fglobal[2]);
          sub3(ftotal, ftorque, fglobal[3]);

          // add forces to nodes

          scaleadd3(0.5, fglobal[0], f[i1], f[i1]);
          scaleadd3(0.5, fglobal[1], f[i2], f[i2]);
          scaleadd3(0.5, fglobal[2], f[jj1], f[jj1]);
          scaleadd3(0.5, fglobal[3], f[jj2], f[jj2]);
          scaleadd3(0.5 * fend, m, f[jj1], f[jj1]);

          // add energy

          if (eflag_either) {
            evdwl = 0.5 * evdwl;
            double evdwl_atom = 0.25 * evdwl;
            if (eflag_global) eng_vdwl += evdwl;
            if (eflag_atom) {
              eatom[i1] += evdwl_atom;
              eatom[i2] += evdwl_atom;
              eatom[jj1] += evdwl_atom;
              eatom[jj2] += evdwl_atom;
            }
          }

          // second force contribution

          param[6] += lp;

          fsemi(param, evdwl, fend, flocal);

          if (evdwl == 0.0) continue;

          // transform to global coordinate system

          matvec(basis[0], basis[1], basis[2], flocal[0], fglobal[0]);
          matvec(basis[0], basis[1], basis[2], flocal[1], fglobal[1]);

          // forces acting on segment 2

          add3(fglobal[0], fglobal[1], ftotal);
          scaleadd3(fend, m, ftotal, ftotal);
          scale3(-0.5, ftotal);

          scaleadd3(lp, m, p, p);

          sub3(r1, p, delr1);
          sub3(r2, p, delr2);
          cross3(delr1, fglobal[0], t1);
          cross3(delr2, fglobal[1], t2);
          add3(t1, t2, torque);

          cross3(torque, m, ftorque);
          scale3(1.0 / lp, ftorque);

          add3(ftotal, ftorque, fglobal[2]);
          sub3(ftotal, ftorque, fglobal[3]);

          // add forces to nodes

          scaleadd3(-0.5, fglobal[0], f[i1], f[i1]);
          scaleadd3(-0.5, fglobal[1], f[i2], f[i2]);
          scaleadd3(-0.5, fglobal[2], f[jj2], f[jj2]);
          scaleadd3(0.5, fglobal[3], f[jj1], f[jj1]);
          sub3(f[jj2], fglobal[3], f[jj2]);
          scaleadd3(-0.5 * fend, m, f[jj2], f[jj2]);

          // add energy

          if (eflag_either) {
            evdwl = -0.5 * evdwl;
            double evdwl_atom = 0.25 * evdwl;
            if (eflag_global) eng_vdwl += evdwl;
            if (eflag_atom) {
              eatom[i1] += evdwl_atom;
              eatom[i2] += evdwl_atom;
              eatom[jj1] += evdwl_atom;
              eatom[jj2] += evdwl_atom;
            }
          }
        }

        if (sumw == 0.0) continue;
        sumw_inv = 1.0 / sumw;

        // mean chain velocity and relative velocity

        add3(vp1, vp2, vp);
        scale3(0.5 * sumw_inv, vp);
        sub3(vp, vr, vrel);
        vtot = dot3(vrel, basis[2]);

        // friction forces

        if (vtot == 0.0)
          fvisc_tot = 0.0;
        else {
          fvisc_tot = fvisc_max / (1.0 + exp(-kvisc * (fabs(vtot) - vvisc))) - fvisc_shift;
          fvisc_tot *= 0.25 * (1 - s(sumw)) * (param[3] - param[2]);
          if (vtot < 0) fvisc_tot = -fvisc_tot;
        }
        scale3(fvisc_tot, basis[2], fvisc);

        // add friction forces to current segment

        add3(fvisc, f[i1], f[i1]);
        add3(fvisc, f[i2], f[i2]);

        // add friction forces to neighbor chain

        for (k = 0; k < clen - 1; k++) {
          if (w[k] == 0.0) continue;
          j1 = chain[j][k];
          j2 = chain[j][k + 1];
          j1 &= NEIGHMASK;
          j2 &= NEIGHMASK;
          scale = w[k] * sumw_inv;

          scaleadd3(-scale, fvisc, f[j1], f[j1]);
          scaleadd3(-scale, fvisc, f[j2], f[j2]);
        }
      } else {

        // segment-chain interaction

        // assign end position

        endflag = end[j];
        if (endflag == 1) {
          endindex = chain[j][0];
          qe = x[endindex];
        } else if (endflag == 2) {
          endindex = chain[j][clen - 1];
          qe = x[endindex];
        }

        // compute substitute straight (semi-)infinite CNT

        zero3(p1);
        zero3(p2);
        zero3(dr1_sumw);
        zero3(dr2_sumw);
        zero3(vp1);
        zero3(vp2);
        zeromat3(q1_dr1_w);
        zeromat3(q2_dr1_w);
        zeromat3(q1_dr2_w);
        zeromat3(q2_dr2_w);
        for (k = 0; k < clen; k++) {
          wnode[k] = 0.0;
          zero3(dq_w[k]);
          zeromat3(q1_dq_w[k]);
          zeromat3(q2_dq_w[k]);
        }
        sumw = 0.0;

        for (k = 0; k < clen - 1; k++) {
          j1 = chain[j][k];
          j2 = chain[j][k + 1];
          j1 &= NEIGHMASK;
          j2 &= NEIGHMASK;
          q1 = x[j1];
          q2 = x[j2];
          vq1 = v[j1];
          vq2 = v[j2];

          weight(r1, r2, q1, q2, w[k], dr1_w, dr2_w, dq1_w, dq2_w);

          if (w[k] == 0.0) {
            if (endflag == 1 && k == 0)
              endflag = 0;
            else if (endflag == 2 && k == clen - 2)
              endflag = 0;
            continue;
          }

          sumw += w[k];
          wnode[k] += w[k];
          wnode[k + 1] += w[k];

          scaleadd3(w[k], q1, p1, p1);
          scaleadd3(w[k], q2, p2, p2);

          // weighted velocity for friction

          scaleadd3(w[k], vq1, vp1, vp1);
          scaleadd3(w[k], vq2, vp2, vp2);

          // weight gradient terms

          add3(dr1_w, dr1_sumw, dr1_sumw);
          add3(dr2_w, dr2_sumw, dr2_sumw);

          outer3(q1, dr1_w, temp);
          plus3(temp, q1_dr1_w, q1_dr1_w);
          outer3(q2, dr1_w, temp);
          plus3(temp, q2_dr1_w, q2_dr1_w);
          outer3(q1, dr2_w, temp);
          plus3(temp, q1_dr2_w, q1_dr2_w);
          outer3(q2, dr2_w, temp);
          plus3(temp, q2_dr2_w, q2_dr2_w);

          add3(dq1_w, dq_w[k], dq_w[k]);
          add3(dq2_w, dq_w[k + 1], dq_w[k + 1]);

          outer3(q1, dq1_w, temp);
          plus3(temp, q1_dq_w[k], q1_dq_w[k]);
          outer3(q1, dq2_w, temp);
          plus3(temp, q1_dq_w[k + 1], q1_dq_w[k + 1]);
          outer3(q2, dq1_w, temp);
          plus3(temp, q2_dq_w[k], q2_dq_w[k]);
          outer3(q2, dq2_w, temp);
          plus3(temp, q2_dq_w[k + 1], q2_dq_w[k + 1]);
        }

        if (sumw == 0.0) continue;

        sumw_inv = 1.0 / sumw;
        scale3(sumw_inv, p1);
        scale3(sumw_inv, p2);

        // compute geometry and forces

        if (endflag == 0) {

          // infinite CNT case

          geometry(r1, r2, p1, p2, nullptr, p, m, param, basis);

          if (param[0] > cutoff) continue;
          if (param[2] >= 0 || param[3] <= 0) {
            double salpha = sin(param[1]);
            double sxi1 = salpha * param[2];
            double sxi2 = salpha * param[3];
            double hsq = param[0] * param[0];

            if (sxi1 * sxi1 + hsq > cutoffsq && sxi2 * sxi2 + hsq > cutoffsq) continue;
          }

          finf(param, evdwl, flocal);

        } else {

          // semi-infinite CNT case
          // endflag == 1: CNT end at start of chain
          // endflag == 2: CNT end at end of chain

          if (endflag == 1)
            geometry(r1, r2, p1, p2, qe, p, m, param, basis);
          else
            geometry(r1, r2, p2, p1, qe, p, m, param, basis);

          if (param[0] > cutoff) continue;
          if (param[2] >= 0 || param[3] <= 0) {
            double hsq = param[0] * param[0];
            double calpha = cos(param[1]);
            double etamin = calpha * param[2];
            double dsq1;
            if (etamin < param[6]) {
              dsq1 = hsq + pow(sin(param[1]) * param[6], 2);
              dsq1 += pow(param[2] - calpha * param[6], 2);
            } else
              dsq1 = hsq + pow(sin(param[1]) * param[2], 2);

            etamin = calpha * param[3];
            double dsq2;
            if (etamin < param[6]) {
              dsq2 = hsq + pow(sin(param[1]) * param[6], 2);
              dsq2 += pow(param[3] - calpha * param[6], 2);
            } else
              dsq2 = hsq + pow(sin(param[1]) * param[3], 2);

            if (dsq1 > cutoffsq && dsq2 > cutoffsq) continue;
          }

          fsemi(param, evdwl, fend, flocal);
        }

        if (evdwl == 0.0) continue;
        evdwl *= 0.5;

        // transform to global coordinate system

        matvec(basis[0], basis[1], basis[2], flocal[0], fglobal[0]);
        matvec(basis[0], basis[1], basis[2], flocal[1], fglobal[1]);

        // mean chain velocity and relative velocity

        add3(vp1, vp2, vp);
        scale3(0.5 * sumw_inv, vp);
        sub3(vp, vr, vrel);
        vtot = dot3(vrel, basis[2]);

        // friction forces

        if (vtot == 0.0)
          fvisc_tot = 0.0;
        else {
          fvisc_tot = fvisc_max / (1.0 + exp(-kvisc * (fabs(vtot) - vvisc))) - fvisc_shift;
          fvisc_tot *= 0.25 * (1 - s(sumw)) * (param[3] - param[2]);
          if (vtot < 0) fvisc_tot = -fvisc_tot;
        }
        scale3(fvisc_tot, basis[2], fvisc);

        // add friction forces to current segment

        add3(fvisc, f[i1], f[i1]);
        add3(fvisc, f[i2], f[i2]);

        // forces acting on approximate chain

        add3(fglobal[0], fglobal[1], ftotal);
        if (endflag) scaleadd3(fend, m, ftotal, ftotal);
        scale3(-0.5, ftotal);

        sub3(r1, p, delr1);
        sub3(r2, p, delr2);
        cross3(delr1, fglobal[0], t1);
        cross3(delr2, fglobal[1], t2);
        add3(t1, t2, torque);

        // additional torque contribution from chain end

        if (endflag) {
          sub3(qe, p, delqe);
          cross3(delqe, m, t3);
          scale3(fend, t3);
          add3(t3, torque, torque);
        }

        cross3(torque, m, ftorque);
        lp = param[5] - param[4];
        scale3(1.0 / lp, ftorque);

        if (endflag == 2) {
          add3(ftotal, ftorque, fglobal[3]);
          sub3(ftotal, ftorque, fglobal[2]);
        } else {
          add3(ftotal, ftorque, fglobal[2]);
          sub3(ftotal, ftorque, fglobal[3]);
        }

        scale3(0.5, fglobal[0]);
        scale3(0.5, fglobal[1]);
        scale3(0.5, fglobal[2]);
        scale3(0.5, fglobal[3]);

        // weight gradient terms acting on current segment

        outer3(p1, dr1_sumw, temp);
        minus3(q1_dr1_w, temp, dr1_p1);
        outer3(p2, dr1_sumw, temp);
        minus3(q2_dr1_w, temp, dr1_p2);
        outer3(p1, dr2_sumw, temp);
        minus3(q1_dr2_w, temp, dr2_p1);
        outer3(p2, dr2_sumw, temp);
        minus3(q2_dr2_w, temp, dr2_p2);

        transpose_matvec(dr1_p1, fglobal[2], fgrad_r1_p1);
        transpose_matvec(dr1_p2, fglobal[3], fgrad_r1_p2);
        transpose_matvec(dr2_p1, fglobal[2], fgrad_r2_p1);
        transpose_matvec(dr2_p2, fglobal[3], fgrad_r2_p2);

        // add forces to nodes in current segment

        add3(fglobal[0], f[i1], f[i1]);
        add3(fglobal[1], f[i2], f[i2]);

        scaleadd3(sumw_inv, fgrad_r1_p1, f[i1], f[i1]);
        scaleadd3(sumw_inv, fgrad_r1_p2, f[i1], f[i1]);
        scaleadd3(sumw_inv, fgrad_r2_p1, f[i2], f[i2]);
        scaleadd3(sumw_inv, fgrad_r2_p2, f[i2], f[i2]);

        // add forces in approximate chain

        for (k = 0; k < clen - 1; k++) {
          if (w[k] == 0.0) continue;
          j1 = chain[j][k];
          j2 = chain[j][k + 1];
          j1 &= NEIGHMASK;
          j2 &= NEIGHMASK;
          scale = w[k] * sumw_inv;
          scaleadd3(scale, fglobal[2], f[j1], f[j1]);
          scaleadd3(scale, fglobal[3], f[j2], f[j2]);

          // friction forces

          scaleadd3(-scale, fvisc, f[j1], f[j1]);
          scaleadd3(-scale, fvisc, f[j2], f[j2]);
        }

        // weight gradient terms acting on approximate chain
        // iterate over nodes instead of segments

        for (k = 0; k < clen; k++) {
          if (wnode[k] == 0.0) continue;
          j1 = chain[j][k];
          j1 &= NEIGHMASK;

          outer3(p1, dq_w[k], temp);
          minus3(q1_dq_w[k], temp, dq_p1);
          outer3(p2, dq_w[k], temp);
          minus3(q2_dq_w[k], temp, dq_p2);

          transpose_matvec(dq_p1, fglobal[2], fgrad_q_p1);
          transpose_matvec(dq_p2, fglobal[3], fgrad_q_p2);

          scaleadd3(sumw_inv, fgrad_q_p1, f[j1], f[j1]);
          scaleadd3(sumw_inv, fgrad_q_p2, f[j1], f[j1]);
        }

        // force on node at CNT end

        if (endflag) scaleadd3(0.5 * fend, m, f[endindex], f[endindex]);

        // compute energy

        if (eflag_either) {
          if (eflag_global) eng_vdwl += evdwl;
          if (eflag_atom) {
            eatom[i1] += 0.25 * evdwl;
            eatom[i2] += 0.25 * evdwl;
            for (k = 0; k < clen - 1; k++) {
              if (w[k] == 0.0) continue;
              j1 = chain[j][k];
              j2 = chain[j][k + 1];
              j1 &= NEIGHMASK;
              j2 &= NEIGHMASK;
              evdwl_chain = 0.5 * w[k] * sumw_inv * evdwl;
              eatom[j1] += evdwl_chain;
              eatom[j2] += evdwl_chain;
            }
          }
        }
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMesoCNTViscous::coeff(int narg, char **arg)
{
  if (narg < 6) utils::missing_cmd_args(FLERR, "pair_coeff", error);
  read_file(arg[2]);

  fvisc_max = utils::numeric(FLERR, arg[3], false, lmp);
  kvisc = utils::numeric(FLERR, arg[4], false, lmp);
  vvisc = utils::numeric(FLERR, arg[5], false, lmp);
  fvisc_shift = fvisc_max / (1.0 + exp(kvisc * vvisc));

  nend_types = narg - 6;

  if (!allocated) allocate();

  // end atom types
  for (int i = 6; i < narg; i++) end_types[i - 6] = utils::inumeric(FLERR, arg[i], false, lmp);

  // units, eV to energy unit conversion
  ang = force->angstrom;
  ang_inv = 1.0 / ang;
  if (strcmp(update->unit_style, "real") == 0)
    eunit = 23.06054966;
  else if (strcmp(update->unit_style, "metal") == 0)
    eunit = 1.0;
  else if (strcmp(update->unit_style, "si") == 0)
    eunit = 1.6021765e-19;
  else if (strcmp(update->unit_style, "cgs") == 0)
    eunit = 1.6021765e-12;
  else if (strcmp(update->unit_style, "electron") == 0)
    eunit = 3.674932248e-2;
  else if (strcmp(update->unit_style, "micro") == 0)
    eunit = 1.6021765e-4;
  else if (strcmp(update->unit_style, "nano") == 0)
    eunit = 1.6021765e2;
  else
    error->all(FLERR, "Pair style mesocnt/viscous does not support {} units", update->unit_style);
  funit = eunit * ang_inv;

  // potential variables
  sig = sig_ang * ang;
  r = r_ang * ang;
  rsq = r * r;
  d = 2.0 * r;
  d_ang = 2.0 * r_ang;
  rc = 3.0 * sig;
  cutoff = rc + d;
  cutoffsq = cutoff * cutoff;
  cutoff_ang = cutoff * ang_inv;
  cutoffsq_ang = cutoff_ang * cutoff_ang;
  comega = 0.275 * (1.0 - 1.0 / (1.0 + 0.59 * r_ang));
  ctheta = 0.35 + 0.0226 * (r_ang - 6.785);

  // compute spline coefficients
  spline_coeff(uinf_data, uinf_coeff, delh_uinf, uinf_points);
  spline_coeff(gamma_data, gamma_coeff, delh_gamma, gamma_points);
  spline_coeff(phi_data, phi_coeff, delh_phi, delpsi_phi, phi_points);
  spline_coeff(usemi_data, usemi_coeff, delh_usemi, delxi_usemi, usemi_points);

  memory->destroy(uinf_data);
  memory->destroy(gamma_data);
  memory->destroy(phi_data);
  memory->destroy(usemi_data);

  // compute Gauss-Legendre quadrature nodes and weights
  gl_init_nodes(QUAD_FINF, gl_nodes_finf);
  gl_init_nodes(QUAD_FSEMI, gl_nodes_fsemi);
  gl_init_weights(QUAD_FINF, gl_nodes_finf, gl_weights_finf);
  gl_init_weights(QUAD_FSEMI, gl_nodes_fsemi, gl_weights_fsemi);

  int ntypes = atom->ntypes;
  for (int i = 1; i <= ntypes; i++)
    for (int j = i; j <= ntypes; j++) setflag[i][j] = 1;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairMesoCNTViscous::init_style()
{
  if (atom->tag_enable == 0) error->all(FLERR, "Pair style mesocnt/viscous requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR, "Pair style mesocnt/viscous requires newton pair on");
  if (force->special_lj[1] == 0.0 || force->special_lj[2] == 0.0 || force->special_lj[3] == 0.0)
    error->all(FLERR, "Pair mesocnt/viscous requires all special_bond lj values to be non-zero");
  if (comm->ghost_velocity == 0)
    error->all(FLERR, "Pair mesocnt/viscous requires ghost atoms store velocity");

  // need a full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   weight for averaged friction from CNT chain
------------------------------------------------------------------------- */

inline double PairMesoCNTViscous::weight(const double *r1, const double *r2, const double *p1,
                                         const double *p2)
{
  double dr, dp, rhoc, rhomin, rho;
  double r[3], p[3];

  add3(r1, r2, r);
  add3(p1, p2, p);

  dr = sqrt(0.25 * distsq3(r1, r2) + rsq);
  dp = sqrt(0.25 * distsq3(p1, p2) + rsq);
  rhoc = dr + dp + rc;
  rhomin = RHOMIN * ang;
  rho = 0.5 * sqrt(distsq3(r, p));

  return s((rho - rhomin) / (rhoc - rhomin));
}

/* ----------------------------------------------------------------------
   weight for substitute CNT chain
   computes gradients with respect to positions
------------------------------------------------------------------------- */

inline void PairMesoCNTViscous::weight(const double *r1, const double *r2, const double *p1,
                                       const double *p2, double &w, double *dr1_w, double *dr2_w,
                                       double *dp1_w, double *dp2_w)
{
  double dr, dp, rhoc, rhomin, rho, frac, arg, factor;
  double r[3], p[3];
  double dr_rho[3], dr_rhoc[3], dp_rhoc[3];

  add3(r1, r2, r);
  add3(p1, p2, p);
  scale3(0.5, r);
  scale3(0.5, p);

  dr = sqrt(0.25 * distsq3(r1, r2) + rsq);
  dp = sqrt(0.25 * distsq3(p1, p2) + rsq);
  rhoc = dr + dp + rc;
  rhomin = RHOMIN * ang;
  rho = sqrt(distsq3(r, p));

  frac = 1.0 / (rhoc - rhomin);
  arg = frac * (rho - rhomin);
  w = s(arg);

  if (w == 0.0 || w == 1.0) {
    zero3(dr1_w);
    zero3(dr2_w);
    zero3(dp1_w);
    zero3(dp2_w);
  } else {
    factor = ds(arg) * frac;

    sub3(r, p, dr_rho);
    sub3(r1, r2, dr_rhoc);
    sub3(p1, p2, dp_rhoc);
    scale3(0.5 / rho, dr_rho);
    scale3(0.25 / dr, dr_rhoc);
    scale3(0.25 / dp, dp_rhoc);

    scaleadd3(-arg, dr_rhoc, dr_rho, dr1_w);
    scaleadd3(arg, dr_rhoc, dr_rho, dr2_w);
    negate3(dr_rho);
    scaleadd3(-arg, dp_rhoc, dr_rho, dp1_w);
    scaleadd3(arg, dp_rhoc, dr_rho, dp2_w);
    scale3(factor, dr1_w);
    scale3(factor, dr2_w);
    scale3(factor, dp1_w);
    scale3(factor, dp2_w);
  }
}
