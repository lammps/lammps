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
   Contributing author: Ray Shan (Sandia)
------------------------------------------------------------------------- */

#include "fix_qeq_point.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixQEqPoint::FixQEqPoint(LAMMPS *lmp, int narg, char **arg) : FixQEq(lmp, narg, arg)
{
  if (narg == 10) {
    if (strcmp(arg[8], "warn") == 0) {
      maxwarn = utils::logical(FLERR, arg[9], false, lmp);
    } else
      error->all(FLERR, "Illegal fix qeq/point command");
  } else if (narg > 8)
    error->all(FLERR, "Illegal fix qeq/point command");
}

/* ---------------------------------------------------------------------- */

void FixQEqPoint::init()
{
  FixQEq::init();

  neighbor->add_request(this, NeighConst::REQ_FULL);

  int ntypes = atom->ntypes;
  memory->create(shld, ntypes + 1, ntypes + 1, "qeq:shielding");
}

// clang-format off

/* ---------------------------------------------------------------------- */

void FixQEqPoint::pre_force(int /*vflag*/)
{
  if (update->ntimestep % nevery) return;

  nlocal = atom->nlocal;

  if (atom->nmax > nmax) reallocate_storage();

  if (nlocal > n_cap*DANGER_ZONE || m_fill > m_cap*DANGER_ZONE)
    reallocate_matrix();

  init_matvec();
  matvecs = CG(b_s, s);         // CG on s - parallel
  matvecs += CG(b_t, t);        // CG on t - parallel
  matvecs /= 2;
  calculate_Q();

  if (force->kspace) force->kspace->qsum_qsq();
}

/* ---------------------------------------------------------------------- */

void FixQEqPoint::init_matvec()
{
  compute_H();

  int inum, ii, i;
  int *ilist;

  inum = list->inum;
  ilist = list->ilist;

  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      Hdia_inv[i] = 1. / eta[atom->type[i]];
      b_s[i]      = -(chi[atom->type[i]] + chizj[i]);
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

void FixQEqPoint::compute_H()
{
  int inum, jnum, *ilist, *jlist, *numneigh, **firstneigh;
  int i, j, ii, jj;
  double **x;
  double dx, dy, dz, r_sqr, r;

  x = atom->x;
  int *mask = atom->mask;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

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
          H.val[m_fill] = 0.5/r;
          m_fill++;
        }
      }
      H.numnbrs[i] = m_fill - H.firstnbr[i];
    }
  }

  if (m_fill >= H.m)
    error->all(FLERR,"Fix qeq/point has insufficient H matrix "
                                 "size: m_fill={} H.m={}\n",m_fill, H.m);
}

/* ---------------------------------------------------------------------- */
