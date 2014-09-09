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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fix_qeq_point.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "force.h"
#include "group.h"
#include "respa.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixQEqPoint::FixQEqPoint(LAMMPS *lmp, int narg, char **arg) :
  FixQEq(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

void FixQEqPoint::init()
{
  if (!atom->q_flag) error->all(FLERR,"Fix qeq/point requires atom attribute q");

  ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR,"Fix qeq/point group has no atoms");

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
  neighbor->requests[irequest]->half = 1;
  neighbor->requests[irequest]->full = 0;

  int ntypes = atom->ntypes;
  memory->create(shld,ntypes+1,ntypes+1,"qeq:shileding");

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

}

/* ---------------------------------------------------------------------- */

void FixQEqPoint::pre_force(int vflag)
{
  if (update->ntimestep % nevery) return;

  n = atom->nlocal;
  N = atom->nlocal + atom->nghost;

  if( atom->nmax > nmax ) reallocate_storage();

  if( n > n_cap*DANGER_ZONE || m_fill > m_cap*DANGER_ZONE )
    reallocate_matrix();

  init_matvec();
  matvecs = CG(b_s, s);    	// CG on s - parallel
  matvecs += CG(b_t, t); 	// CG on t - parallel
  calculate_Q();

}

/* ---------------------------------------------------------------------- */

void FixQEqPoint::init_matvec()
{
  compute_H();

  int nn, ii, i;
  int *ilist;

  nn = list->inum;
  ilist = list->ilist;

  for( ii = 0; ii < nn; ++ii ) {
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

void FixQEqPoint::compute_H()
{
  int inum, jnum, *ilist, *jlist, *numneigh, **firstneigh;
  int i, j, ii, jj, flag;
  double **x, SMALL = 0.0001;
  double dx, dy, dz, r_sqr, r;

  tagint *tag = atom->tag;
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

        dx = x[j][0] - x[i][0];
        dy = x[j][1] - x[i][1];
        dz = x[j][2] - x[i][2];
        r_sqr = dx*dx + dy*dy + dz*dz;

        flag = 0;
        if (r_sqr <= cutoff_sq) {
          if (j < n) flag = 1;
          else if (tag[i] < tag[j]) flag = 1;
          else if (tag[i] == tag[j]) {
            if (dz > SMALL) flag = 1;
            else if (fabs(dz) < SMALL) {
              if (dy > SMALL) flag = 1;
              else if (fabs(dy) < SMALL && dx > SMALL)
                flag = 1;
	    }
	  }
	}

        if( flag ) {
          H.jlist[m_fill] = j;
	  r = sqrt(r_sqr);
          H.val[m_fill] = 1.0/r;
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
    error->all(FLERR,"Fix qeq/point has insufficient QEq matrix size");
  }
}

/* ---------------------------------------------------------------------- */
