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

#include "fix_qeq_shielded.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "force.h"
#include "group.h"
#include "pair.h"
#include "kspace.h"
#include "respa.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixQEqShielded::FixQEqShielded(LAMMPS *lmp, int narg, char **arg) :
  FixQEq(lmp, narg, arg) {
  if (reax_flag) extract_reax();
}

/* ---------------------------------------------------------------------- */

void FixQEqShielded::init()
{
  if (!atom->q_flag)
    error->all(FLERR,"Fix qeq/shielded requires atom attribute q");

  ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR,"Fix qeq/shielded group has no atoms");

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  int ntypes = atom->ntypes;
  memory->create(shld,ntypes+1,ntypes+1,"qeq:shielding");

  init_shielding();

  int i;
  for (i = 1; i <= ntypes; i++) {
    if (gamma[i] == 0.0)
      error->all(FLERR,"Invalid param file for fix qeq/shielded");
  }

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

}

/* ---------------------------------------------------------------------- */

void FixQEqShielded::extract_reax()
{
  Pair *pair = force->pair_match("reax/c",1);
  if (pair == NULL) error->all(FLERR,"No pair reax/c for fix qeq/shielded");
  int tmp;
  chi = (double *) pair->extract("chi",tmp);
  eta = (double *) pair->extract("eta",tmp);
  gamma = (double *) pair->extract("gamma",tmp);
  if (chi == NULL || eta == NULL || gamma == NULL)
    error->all(FLERR,
        "Fix qeq/slater could not extract params from pair reax/c");
}


/* ---------------------------------------------------------------------- */

void FixQEqShielded::init_shielding()
{
  int i,j;
  double d7, swa2, swa3, swb2, swb3;

  int ntypes = atom->ntypes;
  for( i = 1; i <= ntypes; ++i )
    for( j = 1; j <= ntypes; ++j )
      shld[i][j] = pow( gamma[i] * gamma[j], -1.5 );

  if (fabs(swa) > 0.01 && comm->me == 0)
    error->warning(FLERR,"Fix qeq has non-zero lower Taper radius cutoff");
  if (swb < 0)
    error->all(FLERR, "Fix qeq has negative upper Taper radius cutoff");
  else if (swb < 5 && comm->me == 0)
    error->warning(FLERR,"Fix qeq has very low Taper radius cutoff");

  d7 = pow( swb - swa, 7 );
  swa2 = swa*swa;
  swa3 = swa2*swa;
  swb2 = swb*swb;
  swb3 = swb2*swb;

  Tap[7] =  20.0 / d7;
  Tap[6] = -70.0 * (swa + swb) / d7;
  Tap[5] =  84.0 * (swa2 + 3.0*swa*swb + swb2) / d7;
  Tap[4] = -35.0 * (swa3 + 9.0*swa2*swb + 9.0*swa*swb2 + swb3 ) / d7;
  Tap[3] = 140.0 * (swa3*swb + 3.0*swa2*swb2 + swa*swb3 ) / d7;
  Tap[2] =-210.0 * (swa3*swb2 + swa2*swb3) / d7;
  Tap[1] = 140.0 * swa3 * swb3 / d7;
  Tap[0] = (-35.0*swa3*swb2*swb2 + 21.0*swa2*swb3*swb2 -
            7.0*swa*swb3*swb3 + swb3*swb3*swb ) / d7;
}

/* ---------------------------------------------------------------------- */

void FixQEqShielded::pre_force(int /*vflag*/)
{
  if (update->ntimestep % nevery) return;

  nlocal = atom->nlocal;

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

void FixQEqShielded::init_matvec()
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

void FixQEqShielded::compute_H()
{
  int inum, jnum, *ilist, *jlist, *numneigh, **firstneigh;
  int i, j, ii, jj;
  double **x;
  double dx, dy, dz, r_sqr, r;

  int *type = atom->type;
  x = atom->x;
  int *mask = atom->mask;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // fill in the H matrix
  m_fill = 0;
  r_sqr = 0;
  for( ii = 0; ii < inum; ii++ ) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      jlist = firstneigh[i];
      jnum = numneigh[i];
      H.firstnbr[i] = m_fill;

      for( jj = 0; jj < jnum; jj++ ) {
        j = jlist[jj];
        j &= NEIGHMASK;

        dx = x[j][0] - x[i][0];
        dy = x[j][1] - x[i][1];
        dz = x[j][2] - x[i][2];
        r_sqr = dx*dx + dy*dy + dz*dz;

        if (r_sqr <= cutoff_sq) {
          H.jlist[m_fill] = j;
          r = sqrt(r_sqr);
          H.val[m_fill] = 0.5 * calculate_H( r, shld[type[i]][type[j]] );
          m_fill++;
        }
      }
      H.numnbrs[i] = m_fill - H.firstnbr[i];
    }
  }

  if (m_fill >= H.m) {
    char str[128];
    sprintf(str,"H matrix size has been exceeded: m_fill=%d H.m=%d\n",
             m_fill, H.m );
    error->warning(FLERR,str);
    error->all(FLERR,"Fix qeq/shielded has insufficient QEq matrix size");
  }
}

/* ---------------------------------------------------------------------- */

double FixQEqShielded::calculate_H( double r, double gamma )
{
  double Taper, denom;

  Taper = Tap[7] * r + Tap[6];
  Taper = Taper * r + Tap[5];
  Taper = Taper * r + Tap[4];
  Taper = Taper * r + Tap[3];
  Taper = Taper * r + Tap[2];
  Taper = Taper * r + Tap[1];
  Taper = Taper * r + Tap[0];

  denom = r * r * r + gamma;
  denom = pow(denom,0.3333333333333);

  return Taper * EV_TO_KCAL_PER_MOL / denom;
}
