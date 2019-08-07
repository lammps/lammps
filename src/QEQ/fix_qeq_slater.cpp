/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Ray Shan (Sandia)
------------------------------------------------------------------------- */

#include "fix_qeq_slater.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "force.h"
#include "group.h"
#include "pair.h"
#include "kspace.h"
#include "respa.h"
#include "math_const.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixQEqSlater::FixQEqSlater(LAMMPS *lmp, int narg, char **arg) :
  FixQEq(lmp, narg, arg)
{
  alpha = 0.20;

  // optional arg
  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"alpha") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix qeq/slater command");
      alpha = atof(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix qeq/slater command");
  }

  if (streitz_flag) extract_streitz();
}

/* ---------------------------------------------------------------------- */

void FixQEqSlater::init()
{
  if (!atom->q_flag)
    error->all(FLERR,"Fix qeq/slater requires atom attribute q");

  ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR,"Fix qeq/slater group has no atoms");

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  int ntypes = atom->ntypes;
  for (int i = 1; i <= ntypes; i++) {
    if (zeta[i] == 0.0)
      error->all(FLERR,"Invalid param file for fix qeq/slater");
  }

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixQEqSlater::extract_streitz()
{
  Pair *pair = force->pair_match("coul/streitz",1);
  if (pair == NULL) error->all(FLERR,"No pair coul/streitz for fix qeq/slater");
  int tmp;
  chi = (double *) pair->extract("chi",tmp);
  eta = (double *) pair->extract("eta",tmp);
  gamma = (double *) pair->extract("gamma",tmp);
  zeta = (double *) pair->extract("zeta",tmp);
  zcore = (double *) pair->extract("zcore",tmp);
  if (chi == NULL || eta == NULL || gamma == NULL
                  || zeta == NULL || zcore == NULL)
    error->all(FLERR,
        "Fix qeq/slater could not extract params from pair coul/streitz");

}

/* ---------------------------------------------------------------------- */

void FixQEqSlater::pre_force(int /*vflag*/)
{
  if (update->ntimestep % nevery) return;

  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;

  if (atom->nmax > nmax) reallocate_storage();

  if (nlocal > n_cap*DANGER_ZONE || m_fill > m_cap*DANGER_ZONE)
    reallocate_matrix();

  init_matvec();
  matvecs = CG(b_s, s);         // CG on s - parallel
  matvecs += CG(b_t, t);        // CG on t - parallel
  calculate_Q();

  if (force->kspace) force->kspace->qsum_qsq();
}

/* ---------------------------------------------------------------------- */

void FixQEqSlater::init_matvec()
{
  compute_H();

  int inum, ii, i;
  int *ilist;

  inum = list->inum;
  ilist = list->ilist;

  for( ii = 0; ii < inum; ++ii ) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      Hdia_inv[i] = 1. / eta[ atom->type[i] ];
      b_s[i]      = -( chi[atom->type[i]] + chizj[i] );
      b_t[i]      = -1.0;
      t[i] = t_hist[i][2] + 3 * ( t_hist[i][0] - t_hist[i][1] );
      s[i] = 4*(s_hist[i][0]+s_hist[i][2])-(6*s_hist[i][1]+s_hist[i][3]);
    }
  }

  pack_flag = 2;
  comm->forward_comm_fix(this); //Dist_vector( s );
  pack_flag = 3;
  comm->forward_comm_fix(this); //Dist_vector( t );
}

/* ---------------------------------------------------------------------- */

void FixQEqSlater::compute_H()
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double r, rsq, delr[3];
  double zei, zej, zj, zjtmp;

  int *type = atom->type;
  double **x = atom->x;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  m_fill = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    zei = zeta[itype];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    H.firstnbr[i] = m_fill;
    zjtmp = 0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      jtype = type[j];
      zej = zeta[jtype];
      zj = zcore[jtype];

      delr[0] = x[i][0] - x[j][0];
      delr[1] = x[i][1] - x[j][1];
      delr[2] = x[i][2] - x[j][2];
      rsq = delr[0]*delr[0] + delr[1]*delr[1] + delr[2]*delr[2];

      if (rsq > cutoff_sq) continue;

      r = sqrt(rsq);
      H.jlist[m_fill] = j;
      H.val[m_fill] = calculate_H(zei, zej, zj, r, zjtmp);
      m_fill++;
    }
    H.numnbrs[i] = m_fill - H.firstnbr[i];
    chizj[i] = zjtmp;
  }

  if (m_fill >= H.m) {
    char str[128];
    sprintf(str,"H matrix size has been exceeded: m_fill=%d H.m=%d\n",
             m_fill, H.m );
    error->warning(FLERR,str);
    error->all(FLERR,"Fix qeq/slater has insufficient QEq matrix size");
  }

}

/* ---------------------------------------------------------------------- */

double FixQEqSlater::calculate_H(double zei, double zej, double zj,
                double r, double &zjtmp)
{
  double rinv = 1.0/r;

  double exp2zir = exp(-2.0*zei*r);
  double zei2 = zei*zei;
  double zei4 = zei2*zei2;
  double zei6 = zei2*zei4;

  double exp2zjr = exp(-2.0*zej*r);
  double zej2 = zej*zej;
  double zej4 = zej2*zej2;
  double zej6 = zej2*zej4;

  double sm1 = 11.0/8.0;
  double sm2 = 3.00/4.0;
  double sm3 = 1.00/6.0;

  double erfcr = erfc(alpha*r);
  double qqrd2e = force->qqrd2e;

  double etmp1, etmp2;
  double e1, e2, e3, e4;
  double ci_jfi, ci_fifj;

  e1 = e2 = e3 = e4 = 0.0;
  etmp1 = etmp2 = 0.0;

  ci_jfi = -zei*exp2zir - rinv*exp2zir;

  if (zei == zej) {
    ci_fifj = -exp2zir*(rinv + zei*(sm1 + sm2*zei*r + sm3*zei2*r*r));
  } else {
    e1 = zei*zej4/((zei+zej)*(zei+zej)*(zei-zej)*(zei-zej));
    e2 = zej*zei4/((zei+zej)*(zei+zej)*(zej-zei)*(zej-zei));
    e3 = (3.0*zei2*zej4-zej6) /
         ((zei+zej)*(zei+zej)*(zei+zej)*(zei-zej)*(zei-zej)*(zei-zej));
    e4 = (3.0*zej2*zei4-zei6) /
         ((zei+zej)*(zei+zej)*(zei+zej)*(zej-zei)*(zej-zei)*(zej-zei));
    ci_fifj = -exp2zir*(e1+e3/r) - exp2zjr*(e2+e4/r);
  }

  etmp1 = 1.00 * (ci_jfi - ci_fifj);
  etmp2 = 0.50 * (ci_fifj + erfcr*rinv);

  zjtmp += qqrd2e * zj * etmp1;
  return qqrd2e * etmp2;

}

/* ---------------------------------------------------------------------- */

double FixQEqSlater::calculate_H_wolf(double zei, double zej, double zj,
                double r, double &zjtmp)
{
  double rinv = 1.0/r;

  double exp2zir = exp(-2.0*zei*r);
  double zei2 = zei*zei;
  double zei4 = zei2*zei2;
  double zei6 = zei2*zei4;

  double exp2zjr = exp(-2.0*zej*r);
  double zej2 = zej*zej;
  double zej4 = zej2*zej2;
  double zej6 = zej2*zej4;

  double sm1 = 11.0/8.0;
  double sm2 = 3.00/4.0;
  double sm3 = 1.00/6.0;
  double e1, e2, e3, e4;

  double rc = cutoff;
  double rcinv = 1.0/rc;
  double rcinv2 = rcinv*rcinv;
  double exp2zirsh = exp(-2.0*zei*rc);
  double exp2zjrsh = exp(-2.0*zej*rc);

  double eshift, fshift, ci_jfi, ci_fifj;
  double etmp1, etmp2, etmp3;

  double a = alpha;
  double erfcr = erfc(a*r);
  double erfcrc = erfc(a*rc);

  double qqrd2e = force->qqrd2e;

  etmp1 = etmp2 = etmp3 = 0.0;
  e1 = e2 = e3 = e4 = 0.0;

  eshift = -zei*exp2zirsh - rcinv*exp2zirsh;
  fshift = 2.0*zei2*exp2zirsh + rcinv2*exp2zirsh + 2.0*zei*rcinv*exp2zirsh;

  ci_jfi = -zei*exp2zir - rinv*exp2zir - eshift - (r-rc)*fshift;

  if (zei == zej) {
    eshift = -exp2zirsh*(rcinv + zei*(sm1 + sm2*zei*rc + sm3*zei2*rc*rc));
    ci_fifj = -exp2zir*(rinv + zei*(sm1 + sm2*zei*r + sm3*zei2*r*r))
              - eshift - (r-rc)*fshift;
  } else {
    e1 = zei*zej4/((zei+zej)*(zei+zej)*(zei-zej)*(zei-zej));
    e2 = zej*zei4/((zei+zej)*(zei+zej)*(zej-zei)*(zej-zei));
    e3 = (3.0*zei2*zej4-zej6) /
         ((zei+zej)*(zei+zej)*(zei+zej)*(zei-zej)*(zei-zej)*(zei-zej));
    e4 = (3.0*zej2*zei4-zei6) /
         ((zei+zej)*(zei+zej)*(zei+zej)*(zej-zei)*(zej-zei)*(zej-zei));

    eshift = -exp2zirsh*(e1+e3/rc) - exp2zjrsh*(e2+e4/rc);
    ci_fifj = -exp2zir*(e1+e3/r) - exp2zjr*(e2+e4/r)
              - eshift - (r-rc)*fshift;
  }

  etmp1 = erfcr/r - erfcrc/rc;
  etmp2 = 1.00 * (ci_jfi - ci_fifj);
  etmp3 = 0.50 * (etmp1 + ci_fifj);

  zjtmp += qqrd2e * zj * etmp2;
  return qqrd2e * etmp3;

}

/* ---------------------------------------------------------------------- */

void FixQEqSlater::sparse_matvec( sparse_matrix *A, double *x, double *b )
{
  int i, j, itr_j;

  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;

  double r = cutoff;
  double woself = 0.50*erfc(alpha*r)/r + alpha/MY_PIS;

  for( i = 0; i < nlocal; ++i ) {
    if (atom->mask[i] & groupbit)
      b[i] = (eta[atom->type[i]] - 2.0*force->qqr2e*woself) * x[i];
  }

  for( i = nlocal; i < nall; ++i ) {
    if (atom->mask[i] & groupbit)
      b[i] = 0;
  }

  for( i = 0; i < nlocal; ++i ) {
    if (atom->mask[i] & groupbit) {
      for( itr_j=A->firstnbr[i]; itr_j<A->firstnbr[i]+A->numnbrs[i]; itr_j++) {
        j = A->jlist[itr_j];
        b[i] += A->val[itr_j] * x[j];
        b[j] += A->val[itr_j] * x[i];
      }
    }
  }

}
