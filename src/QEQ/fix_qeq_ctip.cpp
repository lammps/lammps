// clang-format off
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
   Contributing author: Gabriel Plummer (NASA)
------------------------------------------------------------------------- */

#include "fix_qeq_ctip.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace EwaldConst;
using namespace MathConst;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixQEqCTIP::FixQEqCTIP(LAMMPS *lmp, int narg, char **arg) :
    FixQEq(lmp, narg, arg), reff(nullptr), reffsq(nullptr), reff4(nullptr), reff7(nullptr),
    s2d_self(nullptr), shield(nullptr), shieldcu(nullptr), reffc(nullptr), reffcsq(nullptr),
    reffc4(nullptr), reffc7(nullptr), s2d_shift(nullptr), f_shift(nullptr), e_shift(nullptr)
{
  cdamp = 0.30;
  maxrepeat = 10;

  // optional args
  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "cdamp") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix qeq/ctip cdamp", error);
      cdamp = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "maxrepeat") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix qeq/ctip maxrepeat", error);
      maxrepeat = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "warn") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix qeq/ctip warn", error);
      maxwarn = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Unknown fix qeq/ctip keyword: {}", arg[iarg]);
  }

  extract_ctip();
}

FixQEqCTIP::~FixQEqCTIP()
{
  delete[] reff;
  delete[] reffsq;
  delete[] reff4;
  delete[] reff7;
  delete[] s2d_self;
  memory->destroy(shield);
  memory->destroy(shieldcu);
  memory->destroy(reffc);
  memory->destroy(reffcsq);
  memory->destroy(reffc4);
  memory->destroy(reffc7);
  memory->destroy(s2d_shift);
  memory->destroy(f_shift);
  memory->destroy(e_shift);
}

/* ---------------------------------------------------------------------- */

void FixQEqCTIP::init()
{
  FixQEq::init();

  neighbor->add_request(this, NeighConst::REQ_FULL);

  int ntypes = atom->ntypes;
  memory->create(shld, ntypes + 1, ntypes + 1, "qeq:shielding");

  delete[] reff;
  delete[] reffsq;
  delete[] reff4;
  delete[] reff7;
  delete[] s2d_self;

  reff = new double[ntypes];
  reffsq = new double[ntypes];
  reff4 = new double[ntypes];
  reff7 = new double[ntypes];
  s2d_self = new double[ntypes];

  double r = cutoff;
  double rsq = r * r;
  double r6 = rsq * rsq * rsq;

  double erfcd_cut = exp(-cdamp * cdamp * rsq);
  double t_cut = 1.0 / (1.0 + EWALD_P * cdamp * r);
  double erfcc_cut =
      (t_cut * (A1 + t_cut * (A2 + t_cut * (A3 + t_cut * (A4 + t_cut * A5)))) * erfcd_cut) / r;

  for (int elt1 = 0; elt1 < ntypes; elt1++) {
    reff[elt1] = std::cbrt(rsq * r + 1.0 / (gamma[elt1 + 1] * gamma[elt1 + 1] * gamma[elt1 + 1]));
    reffsq[elt1] = reff[elt1] * reff[elt1];
    reff4[elt1] = reffsq[elt1] * reffsq[elt1];
    reff7[elt1] = reff4[elt1] * reffsq[elt1] * reff[elt1];
    s2d_self[elt1] = 2.0 * force->qqr2e *
        (1.5 * erfcc_cut + 2.0 * cdamp / MY_PIS * erfcd_cut +
         cdamp * cdamp * cdamp / MY_PIS * rsq * erfcd_cut + 0.5 / reff[elt1] - 1.5 / r +
         r6 / reff7[elt1] + cdamp / MY_PIS);
  }

  memory->destroy(shield);
  memory->destroy(shieldcu);
  memory->destroy(reffc);
  memory->destroy(reffcsq);
  memory->destroy(reffc4);
  memory->destroy(reffc7);
  memory->destroy(s2d_shift);
  memory->destroy(f_shift);
  memory->destroy(e_shift);

  memory->create(shield, ntypes, ntypes, "qeq:shield");
  memory->create(shieldcu, ntypes, ntypes, "qeq:shieldcu");
  memory->create(reffc, ntypes, ntypes, "qeq:reffc");
  memory->create(reffcsq, ntypes, ntypes, "qeq:reffcsq");
  memory->create(reffc4, ntypes, ntypes, "qeq:reffc4");
  memory->create(reffc7, ntypes, ntypes, "qeq:reffc7");
  memory->create(s2d_shift, ntypes, ntypes, "qeq:s2d_shift");
  memory->create(f_shift, ntypes, ntypes, "qeq:f_shift");
  memory->create(e_shift, ntypes, ntypes, "qeq:e_shift");

  double cutoffsq = cutoff * cutoff;
  double cutoffcu = cutoffsq * cutoff;
  double cutoff4 = cutoffsq * cutoffsq;
  double cdampcu = cdamp * cdamp * cdamp;

  for (int elt1 = 0; elt1 < atom->ntypes; elt1++) {
    for (int elt2 = 0; elt2 < atom->ntypes; elt2++) {
      shield[elt1][elt2] = sqrt(gamma[elt1 + 1] * gamma[elt2 + 1]);
      shieldcu[elt1][elt2] = shield[elt1][elt2] * shield[elt1][elt2] * shield[elt1][elt2];
      reffc[elt1][elt2] = cbrt(cutoffcu + 1 / shieldcu[elt1][elt2]);
      reffcsq[elt1][elt2] = reffc[elt1][elt2] * reffc[elt1][elt2];
      reffc4[elt1][elt2] = reffcsq[elt1][elt2] * reffcsq[elt1][elt2];
      reffc7[elt1][elt2] = reffc4[elt1][elt2] * reffcsq[elt1][elt2] * reffc[elt1][elt2];
      s2d_shift[elt1][elt2] = 2.0 * erfcc_cut / cutoffcu +
          4.0 * cdamp / MY_PIS * erfcd_cut / cutoffsq + 4.0 * cdampcu / MY_PIS * erfcd_cut -
          2 / cutoffcu + 4 * cutoff4 / reffc7[elt1][elt2] - 2 * cutoff / reffc4[elt1][elt2];
      f_shift[elt1][elt2] = erfcc_cut / cutoffsq + 2.0 * cdamp / MY_PIS * erfcd_cut / cutoff -
          1 / cutoffsq + cutoffsq / reffc4[elt1][elt2];
      e_shift[elt1][elt2] = erfcc_cut / cutoff + 1 / reffc[elt1][elt2] - 1 / cutoff;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixQEqCTIP::extract_ctip()
{
  Pair *pair = force->pair_match("^coul/ctip",0);
  if (pair == nullptr) error->all(FLERR,"No pair style coul/ctip for fix qeq/ctip");
  int tmp;
  chi = (double *) pair->extract("chi",tmp);
  eta = (double *) pair->extract("eta",tmp);
  gamma = (double *) pair->extract("gamma",tmp);
  zeta = (double *) pair->extract("zeta",tmp);
  zcore = (double *) pair->extract("zcore",tmp);
  qmin = (double *) pair->extract("qmin",tmp);
  qmax = (double *) pair->extract("qmax",tmp);
  omega = (double *) pair->extract("omega",tmp);
  if (chi == nullptr || eta == nullptr || gamma == nullptr || zeta == nullptr ||
      zcore == nullptr || qmin == nullptr || qmax == nullptr || omega == nullptr)
    error->all(FLERR,  "Fix qeq/ctip could not extract all params from pair style coul/ctip");

}

/* ---------------------------------------------------------------------- */

void FixQEqCTIP::pre_force(int /*vflag*/)
{

  int i,n;

  if (update->ntimestep % nevery) return;

  nlocal = atom->nlocal;

  if (atom->nmax > nmax) reallocate_storage();

  if (nlocal > n_cap*DANGER_ZONE || m_fill > m_cap*DANGER_ZONE)
    reallocate_matrix();

  for (i=1; i <= maxrepeat; i++) {
       init_matvec();
       matvecs = CG(b_s, s);         // CG on s - parallel
       matvecs += CG(b_t, t);        // CG on t - parallel
       matvecs /= 2;
       n=calculate_check_Q();
       MPI_Allreduce(&n, &nout, 1, MPI_INT, MPI_SUM, world);
       if (nout == 0) break;
  }

  if (i > maxrepeat && comm->me == 0)
    error->all(FLERR,"Fix qeq some charges not bound within the domain");

  if (force->kspace) force->kspace->qsum_qsq();
}

/* ---------------------------------------------------------------------- */

void FixQEqCTIP::init_matvec()
{
  compute_H();

  int inum, ii, i;
  int *ilist;
  double *q = atom->q, qi;
  int *type = atom->type;

  double r = cutoff;
  double rsq = r*r;

  double erfcd_cut = exp(-cdamp * cdamp * rsq);
  double t_cut = 1.0 / (1.0 + EWALD_P * cdamp * r);

  inum = list->inum;
  ilist = list->ilist;

  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {

      qi=q[i];
      if (qi < qmin[type[i]]) {
        Hdia_inv[i] = 1. / (eta[type[i]]+2*omega[type[i]]-s2d_self[type[i]-1]);
        b_s[i]      = -((chi[type[i]]-2*qmin[type[i]]*omega[type[i]]) + chizj[i]);
      } else if (qi < qmax[type[i]]) {
        Hdia_inv[i] = 1. / (eta[type[i]]-s2d_self[type[i]-1]);
        b_s[i]      = -(chi[type[i]] + chizj[i]);
      } else {
        Hdia_inv[i] = 1. / (eta[type[i]]+2*omega[type[i]]-s2d_self[type[i]-1]);
        b_s[i]      = -((chi[type[i]]-2*qmax[type[i]]*omega[type[i]]) + chizj[i]);
      }

      b_t[i]      = -1.0;
      t[i] = t_hist[i][2] + 3 * (t_hist[i][0] - t_hist[i][1]);
      s[i] = 4*(s_hist[i][0]+s_hist[i][2])-(6*s_hist[i][1]+s_hist[i][3]);
    }
  }

  pack_flag = 2;
  comm->forward_comm(this); //Dist_vector(s);
  pack_flag = 3;
  comm->forward_comm(this); //Dist_vector(t);
}

/* ---------------------------------------------------------------------- */

void FixQEqCTIP::compute_H()
{
  int inum, jnum, *ilist, *jlist, *numneigh, **firstneigh;
  int i, j, ii, jj;
  double dx, dy, dz, r_sqr, r, reff;
  double cutoffsq, erfcd_cut, t_cut;
  double erfcc, erfcd, t;

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  cutoffsq = cutoff * cutoff;
  erfcd_cut = exp(-cdamp * cdamp * cutoffsq);
  t_cut = 1.0 / (1.0 + EWALD_P * cdamp * cutoff);

  // fill in the H matrix
  m_fill = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      jlist = firstneigh[i];
      jnum = numneigh[i];
      H.firstnbr[i] = m_fill;

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        dx = x[j][0] - x[i][0];
        dy = x[j][1] - x[i][1];
        dz = x[j][2] - x[i][2];
        r_sqr = dx*dx + dy*dy + dz*dz;

        if (r_sqr <= cutoff_sq) {
          H.jlist[m_fill] = j;
          r = sqrt(r_sqr);
          reff = cbrt(r_sqr * r + 1 / shieldcu[type[i]-1][type[j]-1]);
          erfcd = exp(-cdamp * cdamp * r_sqr);
          t = 1.0 / (1.0 + EWALD_P * cdamp * r);
          erfcc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * erfcd;
          H.val[m_fill] = 0.5*force->qqr2e*(erfcc/r+1/reff-1/r-e_shift[type[i]-1][type[j]-1]+f_shift[type[i]-1][type[j]-1]*(r-cutoff)-s2d_shift[type[i]-1][type[j]-1]*0.5*(r-cutoff)*(r-cutoff));
          m_fill++;
        }
      }
      H.numnbrs[i] = m_fill - H.firstnbr[i];
    }
  }

  if (m_fill >= H.m)
    error->all(FLERR,"Fix qeq/ctip has insufficient H matrix size: m_fill={} H.m={}\n",m_fill, H.m);
}

/* ---------------------------------------------------------------------- */

void FixQEqCTIP::sparse_matvec(sparse_matrix *A, double *x, double *b)
{
  int i, j, itr_j;
  double *q=atom->q, qi;
  int *type = atom->type;

  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;

  for (i = 0; i < nlocal; ++i) {
    if (atom->mask[i] & groupbit) {
      qi=q[i];
      if (qi < qmin[type[i]]) {
        b[i] = (eta[type[i]]+2*omega[type[i]]-s2d_self[type[i]-1])*x[i];
      } else if (qi < qmax[type[i]]) {
        b[i] = (eta[type[i]]-s2d_self[type[i]-1]) * x[i];
      } else {
        b[i] = (eta[type[i]]+2*omega[type[i]]-s2d_self[type[i]-1])*x[i];
      }
    }
  }

  for (i = nlocal; i < nall; ++i) {
    if (atom->mask[i] & groupbit)
      b[i] = 0;
  }

  for (i = 0; i < nlocal; ++i) {
    if (atom->mask[i] & groupbit) {
      for (itr_j=A->firstnbr[i]; itr_j<A->firstnbr[i]+A->numnbrs[i]; itr_j++) {
        j = A->jlist[itr_j];
        b[i] += A->val[itr_j] * x[j];
        b[j] += A->val[itr_j] * x[i];
      }
    }
  }

}

/* ---------------------------------------------------------------------- */

int FixQEqCTIP::calculate_check_Q()
{
  int i, k, inum, ii;
  int *ilist;
  double u, s_sum, t_sum;
  double *q = atom->q;
  int *type = atom->type;
  double qi_old,qi_new;
  double qi_check1,qi_check2;
  double qi_check3;
  int n;

  inum = list->inum;
  ilist = list->ilist;

  s_sum = parallel_vector_acc( s, inum );
  t_sum = parallel_vector_acc( t, inum);
  u = s_sum / t_sum;

  n = 0;
  for( ii = 0; ii < inum; ++ii ) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      qi_old = q[i];
      q[i] = s[i] - u * t[i];

      for( k = 4; k > 0; --k ) {
        s_hist[i][k] = s_hist[i][k-1];
        t_hist[i][k] = t_hist[i][k-1];
      }
      s_hist[i][0] = s[i];
      t_hist[i][0] = t[i];

      qi_new = q[i];
      qi_check1=(qi_new-qmin[type[i]])*(qi_old-qmin[type[i]]);
      qi_check2=(qi_new-qmax[type[i]])*(qi_old-qmax[type[i]]);
      if ( qi_check1 < 0.0 || qi_check2 < 0.0 ) {
        qi_check3=abs(qi_new-qi_old);
        if (qi_check3 > tolerance) n++;
      }
    }
  }

  pack_flag = 4;
  comm->forward_comm( this ); //Dist_vector( atom->q );

  return n;
}
