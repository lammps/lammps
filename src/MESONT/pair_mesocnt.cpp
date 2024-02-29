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

#include "pair_mesocnt.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"
#include "update.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

using namespace LAMMPS_NS;
using namespace MathExtra;
using MathConst::MY_2PI;
using MathConst::MY_PI;

static constexpr int SELF_CUTOFF = 3;
static constexpr double SMALL = 1.0e-6;
static constexpr double SWITCH = 1.0e-4;
static constexpr double RHOMIN = 10.0;

static constexpr int QUAD_FINF = 129;
static constexpr int QUAD_FSEMI = 10;

static constexpr int BISECTION_STEPS = 1000000;
static constexpr double BISECTION_EPS = 1.0e-15;

/* ---------------------------------------------------------------------- */

PairMesoCNT::PairMesoCNT(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  respa_enable = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  no_virial_fdotr_compute = 0;
  writedata = 0;
  ghostneigh = 0;

  comm_forward = 3;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairMesoCNT::~PairMesoCNT()
{
  if (allocated) {
    memory->destroy(cutsq);
    memory->destroy(setflag);

    memory->destroy(end_types);

    memory->destroy(uinf_coeff);
    memory->destroy(gamma_coeff);
    memory->destroy(phi_coeff);
    memory->destroy(usemi_coeff);

    memory->destroy(numchainlist);
    memory->destroy(nchainlist);
    memory->destroy(endlist);
    memory->destroy(chainlist);

    memory->destroy(selfid);
    memory->destroy(selfpos);

    memory->destroy(w);
    memory->destroy(wnode);
    memory->destroy(dq_w);
    memory->destroy(q1_dq_w);
    memory->destroy(q2_dq_w);

    memory->destroy(param);

    memory->destroy(flocal);
    memory->destroy(fglobal);
    memory->destroy(basis);

    memory->destroy(gl_nodes_finf);
    memory->destroy(gl_nodes_fsemi);
    memory->destroy(gl_weights_finf);
    memory->destroy(gl_weights_fsemi);
  }
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  int i, j, k, i1, i2, j1, j2;
  int endflag, endindex;
  int clen, numchain;
  int *end, *nchain;
  int **chain;
  double fend, lp, scale, sumw, sumw_inv;
  double evdwl, evdwl_chain;
  double *r1, *r2, *q1, *q2, *qe;
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

          if (dsq1 < dsq2) {
            jj1 = j1;
            jj2 = j2;
          } else {
            if (param[1] > MY_PI)
              param[1] -= MY_PI;
            else
              param[1] += MY_PI;

            double temp = -param[5];
            param[5] = -param[4];
            param[4] = temp;
            param[6] = temp;

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
      } else {

        // segment-chain interaction

        // assign end position

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

/* ---------------------------------------------------------------------- */

void PairMesoCNT::allocate()
{
  allocated = 1;
  int ntypes = atom->ntypes;
  int np1 = ntypes + 1;
  int init_size = 1;

  memory->create(cutsq, np1, np1, "pair:cutsq");
  memory->create(setflag, np1, np1, "pair:setflag");
  for (int i = 1; i <= ntypes; i++)
    for (int j = i; j <= ntypes; j++) setflag[i][j] = 0;

  memory->create(end_types, nend_types, "pair:end_types");

  memory->create(uinf_coeff, uinf_points, 4, "pair:uinf_coeff");
  memory->create(gamma_coeff, gamma_points, 4, "pair:gamma_coeff");
  memory->create(phi_coeff, phi_points, phi_points, 4, 4, "pair:phi_coeff");
  memory->create(usemi_coeff, usemi_points, usemi_points, 4, 4, "pair:usemi_coeff");

  memory->create(numchainlist, init_size, "pair:numchainlist");
  memory->create(nchainlist, init_size, init_size, "pair:nchainlist");
  memory->create(endlist, init_size, init_size, "pair:endlist");
  memory->create(chainlist, init_size, init_size, init_size, "pair:chainlist");

  memory->create(selfid, init_size, "pair:selfid");
  memory->create(selfpos, init_size, 2, "pair:selfpos");

  memory->create(w, init_size, "pair:w");
  memory->create(wnode, init_size, "pair:wnode");
  memory->create(dq_w, init_size, 3, "pair:dq_w");
  memory->create(q1_dq_w, init_size, 3, 3, "pair:q1_dq_w");
  memory->create(q2_dq_w, init_size, 3, 3, "pair:q2_dq_w");

  memory->create(param, 7, "pair:param");

  memory->create(flocal, 2, 3, "pair:flocal");
  memory->create(fglobal, 4, 3, "pair:fglobal");
  memory->create(basis, 3, 3, "pair:basis");

  memory->create(gl_nodes_finf, QUAD_FINF, "pair:gl_nodes_finf");
  memory->create(gl_nodes_fsemi, QUAD_FSEMI, "pair:gl_nodes_fsemi");
  memory->create(gl_weights_finf, QUAD_FINF, "pair:gl_weights_finf");
  memory->create(gl_weights_fsemi, QUAD_FSEMI, "pair:gl_weights_fsemi");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMesoCNT::settings(int narg, char **arg)
{
  if (narg < 1)
    utils::missing_cmd_args(FLERR, "pair_coeff", error);
  else if (narg > 3)
    error->all(FLERR, "Too many arguments in pair_style mesocnt command");

  neigh_cutoff = utils::numeric(FLERR, arg[0], false, lmp);

  // segment flag (default false)

  if (narg > 1) {
    if (strcmp(arg[1], "segment") == 0)
      segment_flag = true;
    else if (strcmp(arg[1], "chain") == 0)
      segment_flag = false;
    else
      error->all(
          FLERR,
          "Unknown second argument {} in pair_style mesocnt command, must be 'chain' or 'segment'",
          arg[1]);
  } else
    segment_flag = false;

  // neigh flag (default false)

  if (narg > 2) {
    if (strcmp(arg[2], "topology") == 0)
      neigh_flag = true;
    else if (strcmp(arg[2], "id") == 0)
      neigh_flag = false;
    else
      error->all(
          FLERR,
          "Unknown third argument {} in pair_style mesocnt command, must be 'id' or 'topology'",
          arg[2]);
  } else
    neigh_flag = false;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMesoCNT::coeff(int narg, char **arg)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "pair_coeff", error);
  read_file(arg[2]);

  nend_types = narg - 3;

  if (!allocated) allocate();

  // end atom types
  for (int i = 3; i < narg; i++) end_types[i - 3] = utils::inumeric(FLERR, arg[i], false, lmp);

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
    error->all(FLERR, "Pair style mesocnt does not support {} units", update->unit_style);
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

void PairMesoCNT::init_style()
{
  if (atom->tag_enable == 0) error->all(FLERR, "Pair style mesocnt requires atom IDs");
  if (force->newton_pair == 0) error->all(FLERR, "Pair style mesocnt requires newton pair on");
  if (force->special_lj[1] == 0.0 || force->special_lj[2] == 0.0 || force->special_lj[3] == 0.0)
    error->all(FLERR, "Pair mesocnt requires special_bond lj x y z to have non-zero x, y and z");

  // need a full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMesoCNT::init_one(int /* i */, int /* j */)
{
  return neigh_cutoff;
}

/* ----------------------------------------------------------------------
   update bond neighbor lists based on atom and mol IDs
------------------------------------------------------------------------- */

void PairMesoCNT::bond_neigh_id()
{
  int nlocal = atom->nlocal;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  int *numneigh = list->numneigh;
  int *numneigh_max;
  memory->create(numneigh_max, nbondlist, "pair:numneigh_max");

  for (int i = 0; i < nbondlist; i++) {
    int i1 = bondlist[i][0];
    int i2 = bondlist[i][1];
    int numneigh1, numneigh2;

    // prevent ghost atom with undefined neighbors

    if (i1 > nlocal - 1)
      numneigh1 = 0;
    else
      numneigh1 = numneigh[i1];
    if (i2 > nlocal - 1)
      numneigh2 = 0;
    else
      numneigh2 = numneigh[i2];
    numneigh_max[i] = numneigh1 + numneigh2 + 2;
  }

  // create temporary arrays for chain creation

  memory->create(reduced_nlist, nbondlist, "pair:reduced_nlist");
  memory->create_ragged(reduced_neighlist, nbondlist, numneigh_max, "pair:reduced_neighlist");

  // reduce neighbors to common list and find longest common list size

  for (int i = 0; i < nbondlist; i++) {
    int i1 = bondlist[i][0];
    int i2 = bondlist[i][1];

    neigh_common(i1, i2, reduced_nlist[i], reduced_neighlist[i]);

    // sort list according to atom-id

    sort(reduced_neighlist[i], reduced_nlist[i]);
  }

  // resize chain arrays

  memory->destroy(numchainlist);
  memory->destroy(nchainlist);
  memory->destroy(endlist);
  memory->destroy(chainlist);

  memory->grow(selfid, nbondlist, "pair:selfid");
  memory->grow(selfpos, nbondlist, 2, "pair:selfpos");

  // count neighbor chains per bond

  memory->create(numchainlist, nbondlist, "pair:numchainlist");

  int numchain_max = 0;
  for (int i = 0; i < nbondlist; i++) {
    numchainlist[i] = count_chains(reduced_neighlist[i], reduced_nlist[i]);
    if (numchain_max < numchainlist[i]) numchain_max = numchainlist[i];
  }

  // count neighbor chain lengths per bond

  memory->create_ragged(nchainlist, nbondlist, numchainlist, "pair:nchainlist");

  for (int i = 0; i < nbondlist; i++)
    chain_lengths(reduced_neighlist[i], reduced_nlist[i], nchainlist[i]);

  // create connected neighbor chains and check for ends
  // MEMORY: prevent zero-size array creation for chainlist

  memory->create_ragged(endlist, nbondlist, numchainlist, "pair:endlist");
  if (numchain_max > 0)
    memory->create_ragged(chainlist, nbondlist, numchainlist, nchainlist, "pair:chainlist");
  else
    memory->create(chainlist, 1, 1, 1, "pair:chainlist");

  int nchain_max = 0;
  for (int i = 0; i < nbondlist; i++) {
    int *reduced_neigh = reduced_neighlist[i];
    int *end = endlist[i];
    int *nchain = nchainlist[i];
    int **chain = chainlist[i];

    // set up connected chains and check for ends

    chain_split(reduced_neigh, reduced_nlist[i], nchain, chain, end);

    // find longest chain

    for (int j = 0; j < numchainlist[i]; j++)
      if (nchain_max < nchain[j]) nchain_max = nchain[j];

    // find selfid and selfpos

    tagint tag1 = atom->tag[bondlist[i][0]];
    tagint tag2 = atom->tag[bondlist[i][1]];

    bool found1 = false;
    bool found2 = false;

    for (int j = 0; j < numchainlist[i]; j++) {
      for (int k = 0; k < nchain[j]; k++) {
        tagint temp_tag = atom->tag[chain[j][k]];
        if (tag1 == temp_tag) {
          selfid[i] = j;
          selfpos[i][0] = k;
          found1 = true;
        } else if (tag2 == temp_tag) {
          selfid[i] = j;
          selfpos[i][1] = k;
          found2 = true;
        }

        if (found1 && found2) break;
      }
      if (found1 && found2) break;
    }
  }

  // resize potential arrays

  memory->grow(w, nchain_max, "pair:w");
  memory->grow(wnode, nchain_max, "pair:wnode");
  memory->grow(dq_w, nchain_max, 3, "pair:dq_w");
  memory->grow(q1_dq_w, nchain_max, 3, 3, "pair:q1_dq_w");
  memory->grow(q2_dq_w, nchain_max, 3, 3, "pair:q2_dq_w");

  // destroy temporary arrays for chain creation

  memory->destroy(numneigh_max);
  memory->destroy(reduced_neighlist);
  memory->destroy(reduced_nlist);
}

/* ----------------------------------------------------------------------
   update bond neighbor lists based on bond topology
------------------------------------------------------------------------- */

void PairMesoCNT::bond_neigh_topo()
{
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  int *type = atom->type;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  comm->forward_comm(this);

  // create version of atom->special with local ids and correct images

  int atom1, atom2;
  int **special_local;
  memory->create(special_local, nlocal + nghost, 2, "pair:special_local");

  for (int i = 0; i < nlocal + nghost; i++) {
    atom1 = atom->map(special[i][0]);
    special_local[i][0] = domain->closest_image(i, atom1);
    if (nspecial[i][0] == 1)
      special_local[i][1] = -1;
    else {
      atom2 = atom->map(special[i][1]);
      special_local[i][1] = domain->closest_image(i, atom2);
    }
  }

  int *numneigh = list->numneigh;
  int *numneigh_max;
  memory->create(numneigh_max, nbondlist, "pair:numneigh_max");

  for (int i = 0; i < nbondlist; i++) {
    int i1 = bondlist[i][0];
    int i2 = bondlist[i][1];
    int numneigh1, numneigh2;

    // prevent ghost atom with undefined neighbors

    if (i1 > nlocal - 1)
      numneigh1 = 0;
    else
      numneigh1 = numneigh[i1];
    if (i2 > nlocal - 1)
      numneigh2 = 0;
    else
      numneigh2 = numneigh[i2];
    numneigh_max[i] = numneigh1 + numneigh2 + 2;
  }

  // create temporary arrays for chain creation

  memory->create(reduced_nlist, nbondlist, "pair:reduced_nlist");
  memory->create_ragged(reduced_neighlist, nbondlist, numneigh_max, "pair:reduced_neighlist");

  // reduce neighbors to common list and find longest common list size

  for (int i = 0; i < nbondlist; i++) {
    int i1 = bondlist[i][0];
    int i2 = bondlist[i][1];
    neigh_common(i1, i2, reduced_nlist[i], reduced_neighlist[i]);
  }

  // resize chain arrays

  memory->destroy(numchainlist);
  memory->destroy(nchainlist);
  memory->destroy(endlist);
  memory->destroy(chainlist);

  memory->grow(selfid, nbondlist, "pair:selfid");
  memory->grow(selfpos, nbondlist, 2, "pair:selfpos");

  // split neighbor list into neighbor chains based on bond topology

  int **chainid, **chainpos;
  memory->create_ragged(chainid, nbondlist, reduced_nlist, "pair:chainid");
  memory->create_ragged(chainpos, nbondlist, reduced_nlist, "pair:chainpos");
  memory->create(numchainlist, nbondlist, "pair:numchainlist");

  bool empty_neigh = true;
  for (int i = 0; i < nbondlist; i++) {

    numchainlist[i] = 0;
    if (reduced_nlist[i] == 0) continue;

    // map local ids to reduced neighbor list ids

    std::unordered_map<int, int> reduced_map;
    for (int j = 0; j < reduced_nlist[i]; j++) {
      reduced_map[reduced_neighlist[i][j]] = j;
      chainid[i][j] = -1;
    }

    // assign chain ids and positions

    for (int j = 0; j < reduced_nlist[i]; j++) {

      // skip if atom is already assigned to a chain

      if (chainid[i][j] != -1) continue;

      // iterate along bonded atoms in both directions

      chainid[i][j] = numchainlist[i];
      chainpos[i][j] = 0;

      if (reduced_neighlist[i][j] == bondlist[i][0]) {
        selfid[i] = numchainlist[i];
        selfpos[i][0] = 0;
      } else if (reduced_neighlist[i][j] == bondlist[i][1]) {
        selfid[i] = numchainlist[i];
        selfpos[i][1] = 0;
      }

      int curr_local, next_local;
      int curr_reduced, next_reduced;

      // down the chain: k = 0; up the chain: k = 1

      for (int k = 0; k < 2; k++) {
        curr_local = reduced_neighlist[i][j];
        next_local = special_local[curr_local][k];

        while (next_local != -1) {
          try {
            curr_reduced = reduced_map.at(curr_local);
            next_reduced = reduced_map.at(next_local);
          } catch (const std::out_of_range &) {
            break;
          }

          chainid[i][next_reduced] = numchainlist[i];
          if (k == 0)
            chainpos[i][next_reduced] = chainpos[i][curr_reduced] - 1;
          else
            chainpos[i][next_reduced] = chainpos[i][curr_reduced] + 1;
          if (special_local[next_local][0] != curr_local) {
            curr_local = next_local;
            next_local = special_local[next_local][0];
          } else {
            curr_local = next_local;
            next_local = special_local[next_local][1];
          }

          if (curr_local == bondlist[i][0]) {
            selfid[i] = numchainlist[i];
            selfpos[i][0] = chainpos[i][next_reduced];
          } else if (curr_local == bondlist[i][1]) {
            selfid[i] = numchainlist[i];
            selfpos[i][1] = chainpos[i][next_reduced];
          }
        }
      }
      numchainlist[i]++;
    }
    if (numchainlist[i]) empty_neigh = false;
  }

  memory->destroy(special_local);

  // count neighbor chain lengths per bond

  int **chainpos_min, **chainpos_max;
  memory->create_ragged(chainpos_min, nbondlist, numchainlist, "pair:chainpos_min");
  memory->create_ragged(chainpos_max, nbondlist, numchainlist, "pair:chainpos_max");
  memory->create_ragged(nchainlist, nbondlist, numchainlist, "pair:nchainlist");

  int nchain_max = 0;
  for (int i = 0; i < nbondlist; i++) {
    for (int j = 0; j < numchainlist[i]; j++) {
      chainpos_min[i][j] = 0;
      chainpos_max[i][j] = 0;
    }

    for (int j = 0; j < reduced_nlist[i]; j++) {
      int cid = chainid[i][j];
      int cpos = chainpos[i][j];
      if (cpos < chainpos_min[i][cid]) chainpos_min[i][cid] = cpos;
      if (cpos > chainpos_max[i][cid]) chainpos_max[i][cid] = cpos;
    }

    for (int j = 0; j < numchainlist[i]; j++) {
      int clen = chainpos_max[i][j] - chainpos_min[i][j] + 1;
      nchainlist[i][j] = clen;
      if (clen > nchain_max) nchain_max = clen;
    }
  }

  // create connected neighbor chains and check for ends
  // MEMORY: prevent zero-size array creation for chainlist

  memory->create_ragged(endlist, nbondlist, numchainlist, "pair:endlist");
  if (empty_neigh)
    memory->create(chainlist, 1, 1, 1, "pair:chainlist");
  else
    memory->create_ragged(chainlist, nbondlist, numchainlist, nchainlist, "pair:chainlist");

  for (int i = 0; i < nbondlist; i++) {

    // shift selfpos

    selfpos[i][0] -= chainpos_min[i][selfid[i]];
    selfpos[i][1] -= chainpos_min[i][selfid[i]];

    // sort atoms into chain lists

    for (int j = 0; j < reduced_nlist[i]; j++) {
      int cid = chainid[i][j];
      int cpos = chainpos[i][j] - chainpos_min[i][cid];
      chainlist[i][cid][cpos] = reduced_neighlist[i][j];
    }

    // check for ends

    for (int j = 0; j < numchainlist[i]; j++) {
      int clen = nchainlist[i][j];
      int cstart = chainlist[i][j][0];
      int cend = chainlist[i][j][clen - 1];

      bool estart = match_end(type[cstart]);
      bool eend = match_end(type[cend]);

      if (estart && eend)
        endlist[i][j] = 3;
      else if (estart)
        endlist[i][j] = 1;
      else if (eend)
        endlist[i][j] = 2;
      else
        endlist[i][j] = 0;
    }
  }

  // destroy remaining temporary arrays for chain creation

  memory->destroy(reduced_neighlist);
  memory->destroy(reduced_nlist);

  memory->destroy(chainpos_min);
  memory->destroy(chainpos_max);
  memory->destroy(chainid);
  memory->destroy(chainpos);

  memory->destroy(numneigh_max);

  // resize potential arrays

  memory->grow(w, nchain_max, "pair:w");
  memory->grow(wnode, nchain_max, "pair:wnode");
  memory->grow(dq_w, nchain_max, 3, "pair:dq_w");
  memory->grow(q1_dq_w, nchain_max, 3, 3, "pair:q1_dq_w");
  memory->grow(q2_dq_w, nchain_max, 3, 3, "pair:q2_dq_w");
}

/* ----------------------------------------------------------------------
   extract common neighbor list for bond
------------------------------------------------------------------------- */

void PairMesoCNT::neigh_common(int i1, int i2, int &numred, int *redlist)
{
  int nlocal = atom->nlocal;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int numneigh1, numneigh2;
  int *neighlist1, *neighlist2;

  numred = 0;

  // prevent ghost atom with undefined neighbors

  if (i1 > nlocal - 1 && i2 > nlocal - 1) {
    return;
  } else if (i1 > nlocal - 1) {
    numneigh2 = numneigh[i2];
    neighlist2 = firstneigh[i2];
    redlist[numred++] = i2;
    for (int j = 0; j < numneigh2; j++) redlist[numred++] = neighlist2[j] & NEIGHMASK;
    return;
  } else if (i2 > nlocal - 1) {
    numneigh1 = numneigh[i1];
    neighlist1 = firstneigh[i1];
    redlist[numred++] = i1;
    for (int j = 0; j < numneigh1; j++) redlist[numred++] = neighlist1[j] & NEIGHMASK;
    return;
  } else {
    numneigh1 = numneigh[i1];
    numneigh2 = numneigh[i2];
    neighlist1 = firstneigh[i1];
    neighlist2 = firstneigh[i2];

    for (int j = 0; j < numneigh1; j++) neighlist1[j] &= NEIGHMASK;
    for (int j = 0; j < numneigh2; j++) neighlist2[j] &= NEIGHMASK;

    // sort and unify of neighbor lists

    std::sort(neighlist1, neighlist1 + numneigh1);
    std::sort(neighlist2, neighlist2 + numneigh2);
    std::vector<int> temp;
    std::set_union(neighlist1, neighlist1 + numneigh1, neighlist2, neighlist2 + numneigh2,
                   std::back_inserter(temp));

    for (const auto &j : temp) redlist[numred++] = j;

    return;
  }
}

/* ----------------------------------------------------------------------
   count neighbor chains of given bond
------------------------------------------------------------------------- */

int PairMesoCNT::count_chains(int *redlist, int numred)
{
  if (numred == 0) return 0;

  tagint *tag = atom->tag;
  tagint *mol = atom->molecule;
  int count = 1;

  // split neighbor list and count chains

  for (int j = 0; j < numred - 1; j++) {
    int j1 = redlist[j];
    int j2 = redlist[j + 1];
    if (tag[j2] - tag[j1] != 1 || mol[j1] != mol[j2]) count++;
  }

  return count;
}

/* ----------------------------------------------------------------------
   count lengths of neighbor chains of given bond
------------------------------------------------------------------------- */

void PairMesoCNT::chain_lengths(int *redlist, int numred, int *nchain)
{
  if (numred == 0) return;

  tagint *tag = atom->tag;
  tagint *mol = atom->molecule;
  int clen = 0;
  int cid = 0;

  // split neighbor list into connected chains

  for (int j = 0; j < numred - 1; j++) {
    int j1 = redlist[j];
    int j2 = redlist[j + 1];
    clen++;
    if (tag[j2] - tag[j1] != 1 || mol[j1] != mol[j2]) {
      nchain[cid++] = clen;
      clen = 0;
    }
  }
  clen++;
  nchain[cid++] = clen;
}

/* ----------------------------------------------------------------------
   split neighbors into chains and identify ends
------------------------------------------------------------------------- */

void PairMesoCNT::chain_split(int *redlist, int numred, int *nchain, int **chain, int *end)
{
  if (numred == 0) return;

  tagint *tag = atom->tag;
  tagint *mol = atom->molecule;
  int *type = atom->type;
  int clen = 0;
  int cid = 0;

  // split neighbor list into connected chains

  for (int j = 0; j < numred - 1; j++) {
    int j1 = redlist[j];
    int j2 = redlist[j + 1];
    chain[cid][clen++] = j1;
    if (tag[j2] - tag[j1] != 1 || mol[j1] != mol[j2]) {
      cid++;
      clen = 0;
    }
  }
  chain[cid][clen++] = redlist[numred - 1];
  cid++;

  // check for chain ends

  for (int j = 0; j < cid; j++) {
    int clen = nchain[j];
    int cstart = chain[j][0];
    int cend = chain[j][clen - 1];

    bool estart = match_end(type[cstart]);
    bool eend = match_end(type[cend]);

    if (estart && eend)
      end[j] = 3;
    else if (estart)
      end[j] = 1;
    else if (eend)
      end[j] = 2;
    else
      end[j] = 0;
  }
}

/* ----------------------------------------------------------------------
   insertion sort list according to corresponding atom ID
------------------------------------------------------------------------- */

void PairMesoCNT::sort(int *list, int size)
{
  int i, j, temp1, temp2;
  tagint *tag = atom->tag;
  for (i = 1; i < size; i++) {
    j = i;
    temp1 = list[j - 1];
    temp2 = list[j];
    while (tag[temp1] > tag[temp2]) {
      list[j] = temp1;
      list[j - 1] = temp2;
      j--;
      if (j == 0) break;
      temp1 = list[j - 1];
      temp2 = list[j];
    }
  }
}

/* ---------------------------------------------------------------------- */

int PairMesoCNT::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = ubuf(atom->nspecial[j][0]).d;
    buf[m++] = ubuf(atom->special[j][0]).d;
    if (atom->nspecial[j][0] == 1)
      buf[m++] = ubuf(-1).d;
    else
      buf[m++] = ubuf(atom->special[j][1]).d;
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairMesoCNT::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    atom->nspecial[i][0] = (int) ubuf(buf[m++]).i;
    atom->special[i][0] = (tagint) ubuf(buf[m++]).i;
    if (atom->nspecial[i][0] > 1) atom->special[i][1] = (tagint) ubuf(buf[m]).i;
    m++;
  }
}

/* ----------------------------------------------------------------------
   check if type is in end_types
------------------------------------------------------------------------- */

bool PairMesoCNT::match_end(int type)
{
  for (int i = 0; i < nend_types; i++)
    if (type == end_types[i]) return true;

  return false;
}

/* ----------------------------------------------------------------------
   read mesocnt potential file
------------------------------------------------------------------------- */

void PairMesoCNT::read_file(const char *file)
{
  if (comm->me == 0) {

    // open file and skip first line

    PotentialFileReader reader(lmp, file, "mesocnt");
    reader.skip_line();

    // parse global parameters: 4x integer then 4x double

    try {
      auto values = reader.next_values(4);
      uinf_points = values.next_int();
      gamma_points = values.next_int();
      phi_points = values.next_int();
      usemi_points = values.next_int();

      values = reader.next_values(4);
      r_ang = values.next_double();
      sig_ang = values.next_double();
      delta1 = values.next_double();
      delta2 = values.next_double();
    } catch (std::exception &e) {
      error->one(FLERR, "Error parsing mesocnt potential file header: {}", e.what());
    }

    // allocate and read potential tables

    memory->create(uinf_data, uinf_points, "pair:uinf_data");
    memory->create(gamma_data, gamma_points, "pair:gamma_data");
    memory->create(phi_data, phi_points, phi_points, "pair:phi_data");
    memory->create(usemi_data, usemi_points, usemi_points, "pair:usemi_data");

    read_data(reader, uinf_data, hstart_uinf, delh_uinf, uinf_points);
    read_data(reader, gamma_data, hstart_gamma, delh_gamma, gamma_points);
    read_data(reader, phi_data, hstart_phi, psistart_phi, delh_phi, delpsi_phi, phi_points);
    read_data(reader, usemi_data, hstart_usemi, xistart_usemi, delh_usemi, delxi_usemi,
              usemi_points);
  }

  MPI_Bcast(&uinf_points, 1, MPI_INT, 0, world);
  MPI_Bcast(&gamma_points, 1, MPI_INT, 0, world);
  MPI_Bcast(&phi_points, 1, MPI_INT, 0, world);
  MPI_Bcast(&usemi_points, 1, MPI_INT, 0, world);
  MPI_Bcast(&r_ang, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&sig_ang, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&delta1, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&delta2, 1, MPI_DOUBLE, 0, world);

  // allocate table arrays on other MPI ranks

  if (comm->me != 0) {
    memory->create(uinf_data, uinf_points, "pair:uinf_data");
    memory->create(gamma_data, gamma_points, "pair:gamma_data");
    memory->create(phi_data, phi_points, phi_points, "pair:phi_data");
    memory->create(usemi_data, usemi_points, usemi_points, "pair:usemi_data");
  }

  // broadcast tables

  MPI_Bcast(&hstart_uinf, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&hstart_gamma, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&hstart_phi, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&psistart_phi, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&hstart_usemi, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&xistart_usemi, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&delh_uinf, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&delh_gamma, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&delh_phi, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&delpsi_phi, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&delh_usemi, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&delxi_usemi, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(uinf_data, uinf_points, MPI_DOUBLE, 0, world);
  MPI_Bcast(gamma_data, gamma_points, MPI_DOUBLE, 0, world);
  MPI_Bcast(&phi_data[0][0], phi_points * phi_points, MPI_DOUBLE, 0, world);
  MPI_Bcast(&usemi_data[0][0], usemi_points * usemi_points, MPI_DOUBLE, 0, world);
}

/* ----------------------------------------------------------------------
   read 1D data file. Only called from MPI rank 0
------------------------------------------------------------------------- */

void PairMesoCNT::read_data(PotentialFileReader &reader, double *data, double &xstart, double &dx,
                            int ninput)
{
  // read values from file

  int serror = 0;
  double x, xtemp, dxtemp;

  for (int i = 0; i < ninput; i++) {
    try {
      auto values = reader.next_values(2);

      if (i > 0) xtemp = x;
      x = values.next_double();
      data[i] = values.next_double();

      if (i == 0) {
        xstart = x;
      } else {
        dxtemp = x - xtemp;
        if (i == 1) dx = dxtemp;
        if (fabs(dxtemp - dx) / dx > SMALL) serror++;
      }
    } catch (std::exception &e) {
      error->one(FLERR, "Error parsing data for mesocnt potential: {}", e.what());
    }

    // warn if spacing between data points is not constant

    if (serror)
      error->warning(FLERR, "{} spacings in first column were different from first", serror);
  }
}

/* ----------------------------------------------------------------------
   read 2D data file only called from MPI rank 0
------------------------------------------------------------------------- */

void PairMesoCNT::read_data(PotentialFileReader &reader, double **data, double &xstart,
                            double &ystart, double &dx, double &dy, int ninput)
{
  // read values from file

  int sxerror = 0;
  int syerror = 0;
  double x, y, xtemp, ytemp, dxtemp, dytemp;

  for (int i = 0; i < ninput; i++) {
    try {
      if (i > 0) xtemp = x;
      for (int j = 0; j < ninput; j++) {
        if (j > 0) ytemp = y;

        auto values = reader.next_values(3);
        x = values.next_double();
        y = values.next_double();
        data[i][j] = values.next_double();

        if (i == 0 && j == 0) ystart = y;
        if (j > 0) {
          dytemp = y - ytemp;
          if (j == 1) dy = dytemp;
          if (fabs(dytemp - dy) / dy > SMALL) syerror++;
        }
      }
      if (i == 0) {
        xstart = x;
      } else {
        dxtemp = x - xtemp;
        if (i == 1) dx = dxtemp;
        if (fabs(dxtemp - dx) / dx > SMALL) sxerror++;
      }
    } catch (std::exception &e) {
      error->one(FLERR, "Error parsing data for mesocnt potential: {}", e.what());
    }
  }

  // warn if spacing between data points is not constant

  if (sxerror)
    error->warning(FLERR, "{} spacings in first column were different from first", sxerror);
  if (syerror)
    error->warning(FLERR, "{} spacings in second column were different from first", syerror);
}

/* ----------------------------------------------------------------------
   compute cubic spline coefficients
------------------------------------------------------------------------- */

void PairMesoCNT::spline_coeff(double *data, double **coeff, double dx, int size)
{
  double *u = data;
  double **g = coeff;
  int n = size;

  double d, *p, *bprime, *dprime, **b;
  memory->create(p, n, "pair:p");
  memory->create(b, n, n, "pair:b");
  memory->create(bprime, n, "pair:bprime");
  memory->create(dprime, n, "pair:dprime");

  double dx_inv = 1.0 / dx;
  double dxsq_inv = dx_inv * dx_inv;
  double dxcb_inv = dx_inv * dxsq_inv;

  double ax[4][4] = {{1, 0, 0, 0},
                     {0, 1, 0, 0},
                     {-3 * dxsq_inv, -2 * dx_inv, 3 * dxsq_inv, -dx_inv},
                     {2 * dxcb_inv, dxsq_inv, -2 * dxcb_inv, dxsq_inv}};

  // compute finite difference derivatives at boundaries

  p[0] = (u[1] - u[0]) * dx_inv;
  p[n - 1] = (u[n - 1] - u[n - 2]) * dx_inv;

  // compute derivatives inside domain

  for (int i = 1; i < n - 1; i++) {
    if (i > 1) b[i][i - 1] = dx;
    b[i][i] = 4 * dx;
    if (i < n - 2) b[i][i + 1] = dx;
  }
  bprime[1] = b[1][1];
  for (int i = 2; i < n - 1; i++) bprime[i] = b[i][i] - b[i][i - 1] * b[i - 1][i] / bprime[i - 1];

  for (int i = 1; i < n - 1; i++) {
    d = 3 * (u[i + 1] - u[i - 1]);
    if (i == 1) d -= dx * p[i - 1];
    if (i == n - 2) d -= dx * p[i + 1];
    dprime[i] = d;
    if (i != 1) dprime[i] -= b[i][i - 1] * dprime[i - 1] / bprime[i - 1];
  }

  p[n - 2] = dprime[n - 2] / bprime[n - 2];
  for (int i = n - 3; i > 0; i--) p[i] = (dprime[i] - b[i][i + 1] * p[i + 1]) / bprime[i];

  // compute spline coefficients

  for (int i = 1; i < n; i++) {
    for (int j = 0; j < 4; j++) g[i][j] = 0;

    double k[4] = {u[i - 1], p[i - 1], u[i], p[i]};

    for (int j = 0; j < 4; j++)
      for (int l = 0; l < 4; l++) g[i][j] += ax[j][l] * k[l];
  }

  memory->destroy(p);
  memory->destroy(b);
  memory->destroy(bprime);
  memory->destroy(dprime);
}

/* ----------------------------------------------------------------------
   compute bicubic spline coefficients
------------------------------------------------------------------------- */

void PairMesoCNT::spline_coeff(double **data, double ****coeff, double dx, double dy, int size)
{
  double **u = data;
  double ****g = coeff;
  int n = size;

  double d, *bprime, *dprime, **p, **q, **s, **b;
  memory->create(p, n, n, "pair:p");
  memory->create(q, n, n, "pair:q");
  memory->create(s, n, n, "pair:s");
  memory->create(b, n, n, "pair:b");
  memory->create(bprime, n, "pair:bprime");
  memory->create(dprime, n, "pair:dprime");

  double dx_inv = 1.0 / dx;
  double dy_inv = 1.0 / dy;
  double dxsq_inv = dx_inv * dx_inv;
  double dysq_inv = dy_inv * dy_inv;
  double dxcb_inv = dx_inv * dxsq_inv;
  double dycb_inv = dy_inv * dysq_inv;

  double ax[4][4] = {{1, 0, 0, 0},
                     {0, 1, 0, 0},
                     {-3 * dxsq_inv, -2 * dx_inv, 3 * dxsq_inv, -dx_inv},
                     {2 * dxcb_inv, dxsq_inv, -2 * dxcb_inv, dxsq_inv}};
  double ay[4][4] = {{1, 0, 0, 0},
                     {0, 1, 0, 0},
                     {-3 * dysq_inv, -2 * dy_inv, 3 * dysq_inv, -dy_inv},
                     {2 * dycb_inv, dysq_inv, -2 * dycb_inv, dysq_inv}};

  // compute finite difference derivatives at boundaries

  for (int i = 0; i < n; i++) {
    p[0][i] = (u[1][i] - u[0][i]) * dx_inv;
    p[n - 1][i] = (u[n - 1][i] - u[n - 2][i]) * dx_inv;
  }

  for (int i = 0; i < n; i++) {
    q[i][0] = (u[i][1] - u[i][0]) * dy_inv;
    q[i][n - 1] = (u[i][n - 1] - u[i][n - 2]) * dy_inv;
  }

  s[0][0] = (p[0][1] - p[0][0]) * dy_inv;
  s[0][n - 1] = (p[0][n - 1] - p[0][n - 2]) * dy_inv;
  s[n - 1][0] = (p[n - 1][1] - p[n - 1][0]) * dy_inv;
  s[n - 1][n - 1] = (p[n - 1][n - 1] - p[n - 1][n - 2]) * dy_inv;

  // compute derivatives inside domain

  // sweep in x

  for (int i = 1; i < n - 1; i++) {
    if (i > 1) b[i][i - 1] = dx;
    b[i][i] = 4 * dx;
    if (i < n - 2) b[i][i + 1] = dx;
  }
  bprime[1] = b[1][1];
  for (int i = 2; i < n - 1; i++) bprime[i] = b[i][i] - b[i][i - 1] * b[i - 1][i] / bprime[i - 1];

  // compute p

  for (int j = 0; j < n; j++) {
    for (int i = 1; i < n - 1; i++) {
      d = 3 * (u[i + 1][j] - u[i - 1][j]);
      if (i == 1) d -= dx * p[i - 1][j];
      if (i == n - 2) d -= dx * p[i + 1][j];
      dprime[i] = d;
      if (i != 1) dprime[i] -= b[i][i - 1] * dprime[i - 1] / bprime[i - 1];
    }

    p[n - 2][j] = dprime[n - 2] / bprime[n - 2];
    for (int i = n - 3; i > 0; i--) p[i][j] = (dprime[i] - b[i][i + 1] * p[i + 1][j]) / bprime[i];
  }

  // compute s

  for (int j = 0; j < n; j += n - 1) {
    for (int i = 1; i < n - 1; i++) {
      d = 3 * (q[i + 1][j] - q[i - 1][j]);
      if (i == 1) d -= dx * s[i - 1][j];
      if (i == n - 2) d -= dx * s[i + 1][j];
      dprime[i] = d;
      if (i != 1) dprime[i] -= b[i][i - 1] * dprime[i - 1] / bprime[i - 1];
    }

    s[n - 2][j] = dprime[n - 2] / bprime[n - 2];
    for (int i = n - 3; i > 0; i--) s[i][j] = (dprime[i] - b[i][i + 1] * s[i + 1][j]) / bprime[i];
  }

  // sweep in y

  for (int i = 1; i < n - 1; i++) {
    if (i > 1) b[i][i - 1] = dy;
    b[i][i] = 4 * dy;
    if (i < n - 2) b[i][i + 1] = dy;
  }
  bprime[1] = b[1][1];
  for (int i = 2; i < n - 1; i++) bprime[i] = b[i][i] - b[i][i - 1] * b[i - 1][i] / bprime[i - 1];

  // compute q

  for (int i = 0; i < n; i++) {
    for (int j = 1; j < n - 1; j++) {
      d = 3 * (u[i][j + 1] - u[i][j - 1]);
      if (j == 1) d -= dy * q[i][j - 1];
      if (j == n - 2) d -= dy * q[i][j + 1];
      dprime[j] = d;
      if (j != 1) dprime[j] -= b[j][j - 1] * dprime[j - 1] / bprime[j - 1];
    }

    q[i][n - 2] = dprime[n - 2] / bprime[n - 2];
    for (int j = n - 3; j > 0; j--) q[i][j] = (dprime[j] - b[j][j + 1] * q[i][j + 1]) / bprime[j];
  }

  // compute s

  for (int i = 0; i < n; i++) {
    for (int j = 1; j < n - 1; j++) {
      d = 3 * (p[i][j + 1] - p[i][j - 1]);
      if (j == 1) d -= dy * s[i][j - 1];
      if (j == n - 2) d -= dy * s[i][j + 1];
      dprime[j] = d;
      if (j != 1) dprime[j] -= b[j][j - 1] * dprime[j - 1] / bprime[j - 1];
    }

    s[i][n - 2] = dprime[n - 2] / bprime[n - 2];
    for (int j = n - 3; j > 0; j--) s[i][j] = (dprime[j] - b[j][j + 1] * s[i][j + 1]) / bprime[j];
  }

  for (int i = 1; i < n; i++)
    for (int j = 1; j < n; j++) {
      for (int l = 0; l < 4; l++)
        for (int m = 0; m < 4; m++) g[i][j][l][m] = 0;

      double k[4][4] = {{u[i - 1][j - 1], q[i - 1][j - 1], u[i - 1][j], q[i - 1][j]},
                        {p[i - 1][j - 1], s[i - 1][j - 1], p[i - 1][j], s[i - 1][j]},
                        {u[i][j - 1], q[i][j - 1], u[i][j], q[i][j]},
                        {p[i][j - 1], s[i][j - 1], p[i][j], s[i][j]}};

      for (int l = 0; l < 4; l++)
        for (int m = 0; m < 4; m++)
          for (int n = 0; n < 4; n++)
            for (int o = 0; o < 4; o++) g[i][j][l][m] += ax[l][n] * k[n][o] * ay[m][o];
    }

  memory->destroy(p);
  memory->destroy(q);
  memory->destroy(s);
  memory->destroy(b);
  memory->destroy(bprime);
  memory->destroy(dprime);
}

/* ----------------------------------------------------------------------
   cubic spline evaluation
------------------------------------------------------------------------- */

inline double PairMesoCNT::spline(double x, double xstart, double dx, double **coeff,
                                  int coeff_size)
{
  int i = ceil((x - xstart) / dx);

  // linear extrapolation

  if (i < 1) {
    return coeff[1][0] + coeff[1][1] * (x - xstart);

    // constant extrapolation

  } else if (i > coeff_size - 1) {
    i = coeff_size - 1;
    x = xstart + (coeff_size - 1) * dx;
  }

  // cubic interpolation

  double xlo = xstart + (i - 1) * dx;
  double xbar = x - xlo;

  return coeff[i][0] + xbar * (coeff[i][1] + xbar * (coeff[i][2] + xbar * coeff[i][3]));
}

/* ----------------------------------------------------------------------
   cubic spline derivative
------------------------------------------------------------------------- */

inline double PairMesoCNT::dspline(double x, double xstart, double dx, double **coeff,
                                   int coeff_size)
{
  int i = ceil((x - xstart) / dx);

  // constant extrapolation

  if (i < 1) {
    return coeff[1][1];

    // constant extrapolation

  } else if (i > coeff_size - 1) {
    i = coeff_size - 1;
    x = xstart + (coeff_size - 1) * dx;
  }

  // cubic interpolation

  double xlo = xstart + (i - 1) * dx;
  double xbar = x - xlo;

  return coeff[i][1] + xbar * (2 * coeff[i][2] + 3 * xbar * coeff[i][3]);
}

/* ----------------------------------------------------------------------
   bicubic spline evaluation
------------------------------------------------------------------------- */

inline double PairMesoCNT::spline(double x, double y, double xstart, double ystart, double dx,
                                  double dy, double ****coeff, int coeff_size)
{
  int i = ceil((x - xstart) / dx);
  int j = ceil((y - ystart) / dy);

  // constant extrapolation

  if (i < 1) {
    i = 1;
    x = xstart;
  } else if (i > coeff_size - 1) {
    i = coeff_size - 1;
    x = xstart + (coeff_size - 1) * dx;
  }

  if (j < 1) {
    j = 1;
    y = ystart;
  } else if (j > coeff_size - 1) {
    j = coeff_size - 1;
    y = ystart + (coeff_size - 1) * dy;
  }

  // cubic interpolation

  double xlo = xstart + (i - 1) * dx;
  double ylo = ystart + (j - 1) * dy;
  double xbar = x - xlo;
  double ybar = y - ylo;

  double y0 = coeff[i][j][0][0] +
      ybar * (coeff[i][j][0][1] + ybar * (coeff[i][j][0][2] + ybar * (coeff[i][j][0][3])));
  double y1 = coeff[i][j][1][0] +
      ybar * (coeff[i][j][1][1] + ybar * (coeff[i][j][1][2] + ybar * (coeff[i][j][1][3])));
  double y2 = coeff[i][j][2][0] +
      ybar * (coeff[i][j][2][1] + ybar * (coeff[i][j][2][2] + ybar * (coeff[i][j][2][3])));
  double y3 = coeff[i][j][3][0] +
      ybar * (coeff[i][j][3][1] + ybar * (coeff[i][j][3][2] + ybar * (coeff[i][j][3][3])));

  return y0 + xbar * (y1 + xbar * (y2 + xbar * y3));
}

/* ----------------------------------------------------------------------
   bicubic spline partial x derivative
------------------------------------------------------------------------- */

inline double PairMesoCNT::dxspline(double x, double y, double xstart, double ystart, double dx,
                                    double dy, double ****coeff, int coeff_size)
{
  int i = ceil((x - xstart) / dx);
  int j = ceil((y - ystart) / dy);

  // constant extrapolation

  if (i < 1) {
    i = 1;
    x = xstart;
  } else if (i > coeff_size - 1) {
    i = coeff_size - 1;
    x = xstart + (coeff_size - 1) * dx;
  }

  if (j < 1) {
    j = 1;
    y = ystart;
  } else if (j > coeff_size - 1) {
    j = coeff_size - 1;
    y = ystart + (coeff_size - 1) * dy;
  }

  // cubic interpolation

  double xlo = xstart + (i - 1) * dx;
  double ylo = ystart + (j - 1) * dy;
  double xbar = x - xlo;
  double ybar = y - ylo;

  double y1 = coeff[i][j][1][0] +
      ybar * (coeff[i][j][1][1] + ybar * (coeff[i][j][1][2] + ybar * (coeff[i][j][1][3])));
  double y2 = coeff[i][j][2][0] +
      ybar * (coeff[i][j][2][1] + ybar * (coeff[i][j][2][2] + ybar * (coeff[i][j][2][3])));
  double y3 = coeff[i][j][3][0] +
      ybar * (coeff[i][j][3][1] + ybar * (coeff[i][j][3][2] + ybar * (coeff[i][j][3][3])));

  return y1 + xbar * (2 * y2 + 3 * xbar * y3);
}

/* ----------------------------------------------------------------------
   bicubic spline partial y derivative
------------------------------------------------------------------------- */

inline double PairMesoCNT::dyspline(double x, double y, double xstart, double ystart, double dx,
                                    double dy, double ****coeff, int coeff_size)
{
  int i = ceil((x - xstart) / dx);
  int j = ceil((y - ystart) / dy);

  // constant extrapolation

  if (i < 1) {
    i = 1;
    x = xstart;
  } else if (i > coeff_size - 1) {
    i = coeff_size - 1;
    x = xstart + (coeff_size - 1) * dx;
  }

  if (j < 1) {
    j = 1;
    y = ystart;
  } else if (j > coeff_size - 1) {
    j = coeff_size - 1;
    y = ystart + (coeff_size - 1) * dy;
  }

  // cubic interpolation

  double xlo = xstart + (i - 1) * dx;
  double ylo = ystart + (j - 1) * dy;
  double xbar = x - xlo;
  double ybar = y - ylo;

  double y0 = coeff[i][j][0][1] + ybar * (2 * coeff[i][j][0][2] + 3 * ybar * coeff[i][j][0][3]);
  double y1 = coeff[i][j][1][1] + ybar * (2 * coeff[i][j][1][2] + 3 * ybar * coeff[i][j][1][3]);
  double y2 = coeff[i][j][2][1] + ybar * (2 * coeff[i][j][2][2] + 3 * ybar * coeff[i][j][2][3]);
  double y3 = coeff[i][j][3][1] + ybar * (2 * coeff[i][j][3][2] + 3 * ybar * coeff[i][j][3][3]);

  return y0 + xbar * (y1 + xbar * (y2 + xbar * y3));
}

/* ----------------------------------------------------------------------
   compute local geometric parameters
------------------------------------------------------------------------- */

void PairMesoCNT::geometry(const double *r1, const double *r2, const double *p1, const double *p2,
                           const double *qe, double *p, double *m, double *param, double **basis)
{
  double r[3], delr[3], l[3], rbar[3], pbar[3], delrbar[3];
  double psil[3], psim[3], dell_psim[3], delpsil_m[3];
  double delr1[3], delr2[3], delp1[3], delp2[3], delpqe[3];

  double *ex = basis[0];
  double *ey = basis[1];
  double *ez = basis[2];

  add3(r1, r2, r);
  scale3(0.5, r);
  add3(p1, p2, p);
  scale3(0.5, p);

  sub3(p, r, delr);

  sub3(r2, r1, l);
  norm3(l);
  sub3(p2, p1, m);
  norm3(m);

  double psi = dot3(l, m);
  if (psi > 1.0)
    psi = 1.0;
  else if (psi < -1.0)
    psi = -1.0;
  double denom = 1.0 - psi * psi;

  copy3(l, psil);
  scale3(psi, psil);
  copy3(m, psim);
  scale3(psi, psim);

  double rhoe, etae, taur, taup;
  if (qe) {
    sub3(p, qe, delpqe);
    rhoe = dot3(delpqe, m);
  } else
    rhoe = 0;

  // parallel case

  if (denom < SWITCH) {
    taur = dot3(delr, l) - rhoe * psi;
    taup = -rhoe;
    etae = 0;

    // non-parallel case

  } else {
    double frac = 1.0 / denom;
    sub3(l, psim, dell_psim);
    sub3(psil, m, delpsil_m);
    taur = dot3(delr, dell_psim) * frac;
    taup = dot3(delr, delpsil_m) * frac;
    etae = -rhoe - taup;
  }

  scaleadd3(taur, l, r, rbar);
  scaleadd3(taup, m, p, pbar);
  sub3(pbar, rbar, delrbar);

  double h = len3(delrbar);
  if (h > cutoff) {
    param[0] = h;
    return;
  }
  if (h * ang_inv < SMALL) h = SMALL * ang;

  copy3(delrbar, ex);
  copy3(l, ez);
  scale3(1.0 / h, ex);
  cross3(ez, ex, ey);

  double alpha;
  alpha = (dot3(m, ey) < 0) ? acos(psi) : MY_2PI - acos(psi);

  sub3(r1, rbar, delr1);
  sub3(r2, rbar, delr2);
  sub3(p1, pbar, delp1);
  sub3(p2, pbar, delp2);
  double xi1 = dot3(delr1, l);
  double xi2 = dot3(delr2, l);
  double eta1 = dot3(delp1, m);
  double eta2 = dot3(delp2, m);

  param[0] = h;
  param[1] = alpha;
  param[2] = xi1;
  param[3] = xi2;
  param[4] = eta1;
  param[5] = eta2;
  param[6] = etae;
}

/* ----------------------------------------------------------------------
   weight for substitute CNT chain
   computes gradients with respect to positions
------------------------------------------------------------------------- */

inline void PairMesoCNT::weight(const double *r1, const double *r2, const double *p1,
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

/* ----------------------------------------------------------------------
   forces for infinite CNT case
------------------------------------------------------------------------- */

void PairMesoCNT::finf(const double *param, double &evdwl, double **f)
{
  double h = param[0] * ang_inv;
  double alpha = param[1];
  double xi1 = param[2] * ang_inv;
  double xi2 = param[3] * ang_inv;

  double sin_alpha = sin(alpha);
  double sin_alphasq = sin_alpha * sin_alpha;

  // parallel case

  if (sin_alphasq < SWITCH) {
    double ubar = spline(h, hstart_uinf, delh_uinf, uinf_coeff, uinf_points);
    double delxi = xi2 - xi1;
    f[0][0] = 0.5 * delxi * dspline(h, hstart_uinf, delh_uinf, uinf_coeff, uinf_points) * funit;
    f[1][0] = f[0][0];
    f[0][1] = 0;
    f[1][1] = 0;
    f[0][2] = ubar * funit;
    f[1][2] = -f[0][2];
    evdwl = ubar * delxi * eunit;

    // non-parallel case

  } else {
    double sin_alpha_inv = 1.0 / sin_alpha;
    double sin_alphasq_inv = sin_alpha_inv * sin_alpha_inv;
    double cos_alpha = cos(alpha);
    double cot_alpha = cos_alpha * sin_alpha_inv;

    double omega = 1.0 / (1.0 - comega * sin_alphasq);
    double c1 = omega * sin_alpha;
    double c1_inv = 1.0 / c1;
    double domega = 2.0 * comega * cos_alpha * c1 * omega;

    double gamma_orth = spline(h, hstart_gamma, delh_gamma, gamma_coeff, gamma_points);
    double dgamma_orth = dspline(h, hstart_gamma, delh_gamma, gamma_coeff, gamma_points);
    double gamma = 1.0 + (gamma_orth - 1.0) * sin_alphasq;
    double gamma_inv = 1.0 / gamma;
    double dalpha_gamma = 2.0 * (gamma_orth - 1.0) * sin_alpha * cos_alpha;
    double dh_gamma = dgamma_orth * sin_alphasq;

    double zeta1 = xi1 * c1;
    double zeta2 = xi2 * c1;

    double smooth = s5((h - d_ang - delta1) / (delta2 - delta1));
    double dsmooth = ds5((h - d_ang - delta1) / (delta2 - delta1));
    double g = d_ang + delta2;
    double hsq = h * h;

    double zetaminbar;
    if (h >= g)
      zetaminbar = 0;
    else
      zetaminbar = sqrt(g * g - hsq);
    double zetamin = smooth * zetaminbar;
    double zetamax = sqrt(cutoffsq_ang - hsq);

    double dzetaminbar;
    if (h >= g)
      dzetaminbar = 0;
    else
      dzetaminbar = -h / zetaminbar;
    double dzetamin = dzetaminbar * smooth + zetaminbar * dsmooth / (delta2 - delta1);
    double dzetamax = -h / zetamax;

    double delzeta1 = fabs(zeta1) - zetamin;
    double delzeta2 = fabs(zeta2) - zetamin;

    double zeta_range_inv = 1.0 / (zetamax - zetamin);
    double psi1 = delzeta1 * zeta_range_inv;
    double psi2 = delzeta2 * zeta_range_inv;

    double delta_phi, delta_dh_phi;
    double dzeta_phi1, dzeta_phi2;

    // if psi1 or psi2 are out of interpolation range,
    // use numerical integration to calculate phi and its derivatives directly rather than using splines

    if (psi1 < 0 || psi2 < 0) {
      error->warning(FLERR,
                     "Segment - infinite chain interaction outside of interpolation range. "
                     "Attempting numerical integration, but performance may be poor and simulation "
                     "likely unstable.");

      double scale = 0.5 * (zeta2 - zeta1);
      double shift = 0.5 * (zeta1 + zeta2);

      delta_phi = 0.0;
      delta_dh_phi = 0.0;

      for (int i = 0; i < QUAD_FINF; i++) {
        double zeta = scale * gl_nodes_finf[i] + shift;
        double spline_arg = sqrt(hsq + zeta * zeta);
        delta_phi += gl_weights_finf[i] *
            spline(spline_arg, hstart_uinf, delh_uinf, uinf_coeff, uinf_points);
        delta_dh_phi += gl_weights_finf[i] * h *
            dspline(spline_arg, hstart_uinf, delh_uinf, uinf_coeff, uinf_points) / spline_arg;
      }

      delta_phi *= scale;
      delta_dh_phi *= scale;

      dzeta_phi1 =
          spline(sqrt(hsq + zeta1 * zeta1), hstart_uinf, delh_uinf, uinf_coeff, uinf_points);
      dzeta_phi2 =
          spline(sqrt(hsq + zeta2 * zeta2), hstart_uinf, delh_uinf, uinf_coeff, uinf_points);
    } else {
      double phi1 =
          spline(h, psi1, hstart_phi, psistart_phi, delh_phi, delpsi_phi, phi_coeff, phi_points);
      double phi2 =
          spline(h, psi2, hstart_phi, psistart_phi, delh_phi, delpsi_phi, phi_coeff, phi_points);
      double dh_phibar1 =
          dxspline(h, psi1, hstart_phi, psistart_phi, delh_phi, delpsi_phi, phi_coeff, phi_points);
      double dh_phibar2 =
          dxspline(h, psi2, hstart_phi, psistart_phi, delh_phi, delpsi_phi, phi_coeff, phi_points);
      double dpsi_phibar1 =
          dyspline(h, psi1, hstart_phi, psistart_phi, delh_phi, delpsi_phi, phi_coeff, phi_points);
      double dpsi_phibar2 =
          dyspline(h, psi2, hstart_phi, psistart_phi, delh_phi, delpsi_phi, phi_coeff, phi_points);

      double dzeta_range = dzetamax - dzetamin;
      double dh_psi1 = -zeta_range_inv * (dzetamin + dzeta_range * psi1);
      double dh_psi2 = -zeta_range_inv * (dzetamin + dzeta_range * psi2);
      double dh_phi1 = dh_phibar1 + dpsi_phibar1 * dh_psi1;
      double dh_phi2 = dh_phibar2 + dpsi_phibar2 * dh_psi2;

      dzeta_phi1 = dpsi_phibar1 * zeta_range_inv;
      dzeta_phi2 = dpsi_phibar2 * zeta_range_inv;

      if (zeta1 < 0) {
        phi1 = -phi1;
        dh_phi1 = -dh_phi1;
      }
      if (zeta2 < 0) {
        phi2 = -phi2;
        dh_phi2 = -dh_phi2;
      }

      delta_phi = phi2 - phi1;
      delta_dh_phi = dh_phi2 - dh_phi1;
    }

    double delta_dzeta_phi = dzeta_phi2 - dzeta_phi1;

    double c2 = gamma * c1_inv;
    double u = c2 * delta_phi;
    double c3 = u * gamma_inv;

    double dh_u = dh_gamma * c3 + c2 * delta_dh_phi;
    double dalpha_u = dalpha_gamma * c3 +
        c1_inv * (domega * sin_alpha + omega * cos_alpha) *
            (gamma * (xi2 * dzeta_phi2 - xi1 * dzeta_phi1) - u);

    double lr_inv = 1.0 / (xi2 - xi1);
    double cx = h * gamma * sin_alphasq_inv * delta_dzeta_phi;
    double cy = gamma * cot_alpha * delta_dzeta_phi;

    f[0][0] = lr_inv * (xi2 * dh_u - cx) * funit;
    f[1][0] = lr_inv * (-xi1 * dh_u + cx) * funit;
    f[0][1] = lr_inv * (dalpha_u - xi2 * cy) * funit;
    f[1][1] = lr_inv * (-dalpha_u + xi1 * cy) * funit;
    f[0][2] = gamma * dzeta_phi1 * funit;
    f[1][2] = -gamma * dzeta_phi2 * funit;
    evdwl = u * eunit;
  }
}

/* ----------------------------------------------------------------------
   forces for semi-infinite CNT case
------------------------------------------------------------------------- */

void PairMesoCNT::fsemi(const double *param, double &evdwl, double &fend, double **f)
{
  double h = param[0] * ang_inv;
  double alpha = param[1];
  double xi1 = param[2] * ang_inv;
  double xi2 = param[3] * ang_inv;
  double etae = param[6] * ang_inv;

  double sin_alpha = sin(alpha);
  double sin_alphasq = sin_alpha * sin_alpha;
  double cos_alpha = cos(alpha);

  double omega = 1.0 / (1.0 - comega * sin_alphasq);
  double omegasq = omega * omega;
  double domega = 2 * comega * sin_alpha * cos_alpha * omegasq;

  double theta = 1.0 - ctheta * sin_alphasq;
  double dtheta = -2 * ctheta * sin_alpha * cos_alpha;

  double c1 = omega * sin_alpha;
  double c1sq = c1 * c1;
  double c2 = theta * etae;

  double gamma_orth = spline(h, hstart_gamma, delh_gamma, gamma_coeff, gamma_points);
  double dgamma_orth = dspline(h, hstart_gamma, delh_gamma, gamma_coeff, gamma_points);
  double gamma = 1.0 + (gamma_orth - 1) * sin_alphasq;
  double gamma_inv = 1.0 / gamma;
  double dalpha_gamma = 2 * (gamma_orth - 1) * sin_alpha * cos_alpha;
  double dh_gamma = dgamma_orth * sin_alphasq;

  double scale = 0.5 * (xi2 - xi1);
  double shift = 0.5 * (xi1 + xi2);
  double c3 = gamma * scale;

  double jh = 0;
  double jh1 = 0;
  double jh2 = 0;
  double jxi = 0;
  double jxi1 = 0;
  double ubar = 0;

  for (int i = 0; i < QUAD_FSEMI; i++) {
    double xibar = scale * gl_nodes_fsemi[i] + shift;
    double g = xibar * c1;
    double hbar = sqrt(h * h + g * g);
    double thetabar = xibar * cos_alpha - c2;

    double c = gl_weights_fsemi[i];

    double u = c *
        spline(hbar, thetabar, hstart_usemi, xistart_usemi, delh_usemi, delxi_usemi, usemi_coeff,
               usemi_points);
    double uh;
    if (hbar == 0)
      uh = 0;
    else
      uh = c / hbar *
          dxspline(hbar, thetabar, hstart_usemi, xistart_usemi, delh_usemi, delxi_usemi,
                   usemi_coeff, usemi_points);
    double uxi = c *
        dyspline(hbar, thetabar, hstart_usemi, xistart_usemi, delh_usemi, delxi_usemi, usemi_coeff,
                 usemi_points);

    double uh1 = xibar * uh;
    jh += uh;
    jh1 += uh1;
    jh2 += xibar * uh1;
    jxi += uxi;
    jxi1 += xibar * uxi;
    ubar += u;
  }

  jh *= c3;
  jh1 *= c3;
  jh2 *= c3;
  jxi *= c3;
  jxi1 *= c3;
  ubar *= c3;

  double c4 = gamma_inv * ubar;
  double dh_ubar = dh_gamma * c4 + h * jh;
  double dalpha_ubar = dalpha_gamma * c4 + c1 * (domega * sin_alpha + omega * cos_alpha) * jh2 -
      sin_alpha * jxi1 - dtheta * etae * jxi;

  double cx = h * (omegasq * jh1 + cos_alpha * ctheta * jxi);
  double cy = sin_alpha * (cos_alpha * omegasq * jh1 + (ctheta - 1) * jxi);
  double cz1 = c1sq * jh1 + cos_alpha * jxi;
  double cz2 = c1sq * jh2 + cos_alpha * jxi1;

  double l_inv = 1.0 / (xi2 - xi1);
  f[0][0] = l_inv * (xi2 * dh_ubar - cx) * funit;
  f[1][0] = l_inv * (cx - xi1 * dh_ubar) * funit;
  f[0][1] = l_inv * (dalpha_ubar - xi2 * cy) * funit;
  f[1][1] = l_inv * (xi1 * cy - dalpha_ubar) * funit;
  f[0][2] = l_inv * (cz2 + ubar - xi2 * cz1) * funit;
  f[1][2] = l_inv * (xi1 * cz1 - cz2 - ubar) * funit;
  evdwl = ubar * eunit;

  fend = theta * jxi * funit;
}

/* ----------------------------------------------------------------------
   n-th Legendre polynomial
------------------------------------------------------------------------- */

double PairMesoCNT::legendre(int n, double x)
{
  if (n == 0) return 1.0;
  if (n == 1) return x;

  std::vector<double> lcache(n + 1, 0.0);

  lcache[0] = 1.0;
  lcache[1] = x;

  for (int i = 2; i <= n; i++)
    lcache[i] = ((2 * i - 1) * x * lcache[i - 1] - (i - 1) * lcache[i - 2]) / i;

  return lcache[n];
}

/* ----------------------------------------------------------------------
   find all roots of Legendre polynomial
------------------------------------------------------------------------- */

void PairMesoCNT::gl_init_nodes(int quad, double *gl_nodes)
{
  int k_start, k_end, k_offset;
  if (quad % 2) {
    k_start = 1;
    k_end = (quad - 1) / 2 + 1;
    k_offset = 2;
    gl_nodes[k_end - 1] = 0.0;
  } else {
    k_start = 0;
    k_end = quad / 2;
    k_offset = 1;
  }

  int root = 0;
  for (int k = k_start; k < k_end; k++) {

    double theta = (ceil(0.5 * quad) - 0.25 - k) * MY_PI / (quad + 0.5);
    double a = cos((ceil(0.5 * quad) - k) * MY_PI / (quad + 1.0));
    double b = cos(theta);
    double c;

    // perform bisection

    int iter = 0;

    do {
      c = 0.5 * (a + b);
      if (legendre(quad, c) == 0.0) break;
      if (legendre(quad, a) * legendre(quad, c) < 0)
        b = c;
      else
        a = c;
      iter++;
    } while (fabs(a - b) >= BISECTION_EPS && iter <= BISECTION_STEPS);

    gl_nodes[k_end + root] = c;
    gl_nodes[k_end - root - k_offset] = -c;
    root++;
  }
}

/* ----------------------------------------------------------------------
   find all Gauss-Legendre quadrature weights
------------------------------------------------------------------------- */

void PairMesoCNT::gl_init_weights(int quad, double *gl_nodes, double *gl_weights)
{
  for (int i = 0; i < quad; i++) {
    double x = gl_nodes[i];
    double dlegendre = quad * (x * legendre(quad, x) - legendre(quad - 1, x)) / (x * x - 1.0);

    gl_weights[i] = 2.0 / ((1.0 - x * x) * dlegendre * dlegendre);
  }
}
