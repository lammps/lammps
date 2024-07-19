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

#include "fix_qeq_shielded.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixQEqShielded::FixQEqShielded(LAMMPS *lmp, int narg, char **arg) : FixQEq(lmp, narg, arg)
{
  if (narg == 10) {
    if (strcmp(arg[8], "warn") == 0) {
      maxwarn = utils::logical(FLERR, arg[9], false, lmp);
    } else
      error->all(FLERR, "Illegal fix qeq/shielded command");
  } else if (narg > 8)
    error->all(FLERR, "Illegal fix qeq/shielded command");
  if (reax_flag) extract_reax();
}

/* ---------------------------------------------------------------------- */

void FixQEqShielded::init()
{
  FixQEq::init();

  neighbor->add_request(this, NeighConst::REQ_FULL);

  const int ntypes = atom->ntypes;
  memory->create(shld, ntypes + 1, ntypes + 1, "qeq:shielding");

  init_shielding();

  // check if valid parameters for all atom types in the fix group are provided
  const int *type = atom->type;
  const int *mask = atom->mask;
  int tmp = 0, tmp_all = 0;
  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      if (gamma[type[i]] == 0.0)
        tmp = type[i];
    }
  }
  MPI_Allreduce(&tmp, &tmp_all, 1, MPI_INT, MPI_MAX, world);
  if (tmp_all)
    error->all(FLERR, "Invalid QEq parameters for atom type {} provided", tmp_all);
}

/* ---------------------------------------------------------------------- */

void FixQEqShielded::extract_reax()
{
  const int nlocal = atom->nlocal;
  const int *mask = atom->mask;
  const int *type = atom->type;

  Pair *pair = force->pair_match("^reaxff", 0);
  if (pair == nullptr) error->all(FLERR, "No reaxff pair style for fix qeq/shielded");
  int tmp, tmp_all;
  chi = (double *) pair->extract("chi", tmp);
  eta = (double *) pair->extract("eta", tmp);
  gamma = (double *) pair->extract("gamma", tmp);
  if ((chi == nullptr) || (eta == nullptr) || (gamma == nullptr))
    error->all(FLERR, "Fix qeq/shielded could not extract all QEq parameters from pair reaxff");

  tmp = tmp_all = 0;
  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      if ((chi[type[i]] == 0.0) && (eta[type[i]] == 0.0) && (gamma[type[i]] == 0.0))
        tmp = type[i];
    }
  }
  MPI_Allreduce(&tmp, &tmp_all, 1, MPI_INT, MPI_MAX, world);
  if (tmp_all)
    error->all(FLERR, "No QEq parameters for atom type {} provided by pair reaxff", tmp_all);
}

// clang-format off
/* ---------------------------------------------------------------------- */

void FixQEqShielded::init_shielding()
{
  int i,j;
  double d7, swa2, swa3, swb2, swb3;

  int ntypes = atom->ntypes;
  for (i = 1; i <= ntypes; ++i)
    for (j = 1; j <= ntypes; ++j)
      shld[i][j] = pow(gamma[i] * gamma[j], -1.5);

  if (fabs(swa) > 0.01 && comm->me == 0)
    error->warning(FLERR,"Fix qeq has non-zero lower Taper radius cutoff");
  if (swb < 0)
    error->all(FLERR, "Fix qeq has negative upper Taper radius cutoff");
  else if (swb < 5 && comm->me == 0)
    error->warning(FLERR,"Fix qeq has very low Taper radius cutoff");

  d7 = pow(swb - swa, 7);
  swa2 = swa*swa;
  swa3 = swa2*swa;
  swb2 = swb*swb;
  swb3 = swb2*swb;

  Tap[7] =  20.0 / d7;
  Tap[6] = -70.0 * (swa + swb) / d7;
  Tap[5] =  84.0 * (swa2 + 3.0*swa*swb + swb2) / d7;
  Tap[4] = -35.0 * (swa3 + 9.0*swa2*swb + 9.0*swa*swb2 + swb3) / d7;
  Tap[3] = 140.0 * (swa3*swb + 3.0*swa2*swb2 + swa*swb3) / d7;
  Tap[2] =-210.0 * (swa3*swb2 + swa2*swb3) / d7;
  Tap[1] = 140.0 * swa3 * swb3 / d7;
  Tap[0] = (-35.0*swa3*swb2*swb2 + 21.0*swa2*swb3*swb2 -
            7.0*swa*swb3*swb3 + swb3*swb3*swb) / d7;
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
  matvecs /= 2;
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
          H.val[m_fill] = 0.5 * calculate_H(r, shld[type[i]][type[j]]);
          m_fill++;
        }
      }
      H.numnbrs[i] = m_fill - H.firstnbr[i];
    }
  }

  if (m_fill >= H.m)
    error->all(FLERR,"Fix qeq/shielded has insufficient H matrix "
                                 "size: m_fill={} H.m={}\n",m_fill,H.m);
}

/* ---------------------------------------------------------------------- */

double FixQEqShielded::calculate_H(double r, double gamma)
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
  denom = pow(denom,1.0/3.0);

  return Taper * EV_TO_KCAL_PER_MOL / denom;
}
