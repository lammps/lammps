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
   Contributing author:  G. Ziegenhain, gerolf@ziegenhain.com
                         Copyright (C) 2007
   Updated algorithm by: Brian Barnes, brian.c.barnes11.civ@mail.mil
------------------------------------------------------------------------- */

#include "compute_ackland_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <utility>

using namespace LAMMPS_NS;

enum{UNKNOWN,BCC,FCC,HCP,ICO};

/* ---------------------------------------------------------------------- */

ComputeAcklandAtom::ComputeAcklandAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if ((narg < 3) || (narg > 5))
    error->all(FLERR,"Illegal compute ackland/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  structure = nullptr;
  maxneigh = 0;
  legacy = 0;
  distsq = nullptr;
  nearest = nullptr;
  nearest_n0 = nullptr;
  nearest_n1 = nullptr;

  int iarg = 3;
  while (narg > iarg) {
    if (strcmp("legacy",arg[iarg]) == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid compute ackland/atom command");
      legacy = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    }
  }
}

/* ---------------------------------------------------------------------- */

ComputeAcklandAtom::~ComputeAcklandAtom()
{
  memory->destroy(structure);
  memory->destroy(distsq);
  memory->destroy(nearest);
  memory->destroy(nearest_n0);
  memory->destroy(nearest_n1);
}

/* ---------------------------------------------------------------------- */

void ComputeAcklandAtom::init()
{
  // need an occasional full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"ackland/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute ackland/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeAcklandAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeAcklandAtom::compute_peratom()
{
  int i,j,ii,jj,k,n,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int chi[8];

  invoked_peratom = update->ntimestep;

  // grow structure array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(structure);
    nmax = atom->nmax;
    memory->create(structure,nmax,"compute/ackland/atom:ackland");
    vector_atom = structure;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute structure parameter for each atom in group
  // use full neighbor list

  double **x = atom->x;
  int *mask = atom->mask;
  double cutsq = force->pair->cutforce * force->pair->cutforce;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // ensure distsq and nearest arrays are long enough

      if (jnum > maxneigh) {
        memory->destroy(distsq);
        memory->destroy(nearest);
        memory->destroy(nearest_n0);
        memory->destroy(nearest_n1);
        maxneigh = jnum;
        memory->create(distsq,maxneigh,"compute/ackland/atom:distsq");
        memory->create(nearest,maxneigh,"compute/ackland/atom:nearest");
        memory->create(nearest_n0,maxneigh,"compute/ackland/atom:nearest_n0");
        memory->create(nearest_n1,maxneigh,"compute/ackland/atom:nearest_n1");
      }

      // loop over list of all neighbors within force cutoff
      // distsq[] = distance sq to each
      // nearest[] = atom indices of neighbors

      n = 0;
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq) {
          distsq[n] = rsq;
          nearest[n++] = j;
        }
      }

      // Select 6 nearest neighbors

      select2(6,n,distsq,nearest);

      // Mean squared separation

      double r0_sq = 0.;
      for (j = 0; j < 6; j++)
        r0_sq += distsq[j];
      r0_sq /= 6.;

      // n0 near neighbors with: distsq<1.45*r0_sq
      // n1 near neighbors with: distsq<1.55*r0_sq

      double n0_dist_sq = 1.45*r0_sq,
        n1_dist_sq = 1.55*r0_sq;
      int n0 = 0, n1 = 0;
      for (j = 0; j < n; j++) {
        if (distsq[j] < n1_dist_sq) {
          nearest_n1[n1++] = nearest[j];
          if (distsq[j] < n0_dist_sq) {
            nearest_n0[n0++] = nearest[j];
          }
        }
      }

      // Evaluate all angles <(r_ij,rik) forall n0 particles with:
      // distsq < 1.45*r0_sq

      double bond_angle;
      double norm_j, norm_k;
      chi[0] = chi[1] = chi[2] = chi[3] = chi[4] = chi[5] = chi[6] = chi[7] = 0;
      double x_ij, y_ij, z_ij, x_ik, y_ik, z_ik;
      for (j = 0; j < n0; j++) {
        x_ij = x[i][0]-x[nearest_n0[j]][0];
        y_ij = x[i][1]-x[nearest_n0[j]][1];
        z_ij = x[i][2]-x[nearest_n0[j]][2];
        norm_j = sqrt (x_ij*x_ij + y_ij*y_ij + z_ij*z_ij);
        if (norm_j <= 0.) continue;
        for (k = j+1; k < n0; k++) {
          x_ik = x[i][0]-x[nearest_n0[k]][0];
          y_ik = x[i][1]-x[nearest_n0[k]][1];
          z_ik = x[i][2]-x[nearest_n0[k]][2];
          norm_k = sqrt (x_ik*x_ik + y_ik*y_ik + z_ik*z_ik);
          if (norm_k <= 0.)
            continue;

          bond_angle = (x_ij*x_ik + y_ij*y_ik + z_ij*z_ik) / (norm_j*norm_k);

          // Histogram for identifying the relevant peaks

          if (bond_angle < -0.945) chi[0]++;
          else if (bond_angle < -0.915) chi[1]++;
          else if (bond_angle < -0.755) chi[2]++;
          else if (bond_angle < -0.195) chi[3]++;
          else if (bond_angle < 0.195) chi[4]++;
          else if (bond_angle < 0.245) chi[5]++;
          else if (bond_angle < 0.795) chi[6]++;
          else chi[7]++;
        }
      }
      if (legacy) {

        // This is the original implementation by Gerolf Ziegenhain
        // Deviations from the different lattice structures

        double delta_bcc = 0.35*chi[4]/(double)(chi[5]+chi[6]-chi[4]);
        double delta_cp = fabs(1.-(double)chi[6]/24.);
        double delta_fcc = 0.61*(fabs((double)(chi[0]+chi[1]-6.))+
                                 (double)chi[2])/6.0;
        double delta_hcp = (fabs((double)chi[0]-3.)+
                            fabs((double)chi[0]+(double)chi[1]+
                                 (double)chi[2]+(double)chi[3]-9.0))/12.0;

        // Identification of the local structure according to the reference

        if (chi[0] == 7)       { delta_bcc = 0.; }
        else if (chi[0] == 6)  { delta_fcc = 0.; }
        else if (chi[0] <= 3)  { delta_hcp = 0.; }

        if (chi[7] > 0.)
          structure[i] = UNKNOWN;
        else
          if (chi[4] < 3.)
            {
              if (n1 > 13 || n1 < 11)
                structure[i] = UNKNOWN;
              else
                structure[i] = ICO;
            } else
            if (delta_bcc <= delta_cp)
              {
                if (n1 < 11)
                  structure[i] = UNKNOWN;
                else
                  structure[i] = BCC;
              } else
              if (n1 > 12 || n1 < 11)
                structure[i] = UNKNOWN;
              else
                if (delta_fcc < delta_hcp)
                  structure[i] = FCC;
                else
                  structure[i] = HCP;

      } else {

        // This is the updated implementation by Brian Barnes

        if (chi[7] > 0 || n0 < 11) structure[i] = UNKNOWN;
        else if (chi[0] == 7) structure[i] = BCC;
        else if (chi[0] == 6) structure[i] = FCC;
        else if (chi[0] == 3) structure[i] = HCP;
        else {
          // Deviations from the different lattice structures

          double delta_cp = fabs(1.-(double)chi[6]/24.);

          // ensure we do not get divide by zero
          // and if we will, make delta_bcc irrelevant
          double delta_bcc = delta_cp + 1.0;
          int chi56m4 = chi[5]+chi[6]-chi[4];

          // note that chi[7] presumed zero
          if (chi56m4 != 0) delta_bcc = 0.35*chi[4]/(double)chi56m4;

          double delta_fcc = 0.61*(fabs((double)(chi[0]+chi[1]-6))
                                   +(double)chi[2])/6.0;

          double delta_hcp = (fabs((double)chi[0]-3.)
                              +fabs((double)chi[0]
                                    +(double)chi[1]
                                    +(double)chi[2]
                                    +(double)chi[3]
                                    -9.0))/12.0;

          // Identification of the local structure according to the reference

          if (delta_bcc >= 0.1 && delta_cp >= 0.1 && delta_fcc >= 0.1
              && delta_hcp >= 0.1) structure[i] = UNKNOWN;

          // not part of Ackland-Jones 2006; included for backward compatibility
          if (chi[4] < 3. && n1 == 12) structure[i] = ICO;

          else {
            if (delta_bcc <= delta_cp && n1 > 10 && n1 < 13) structure[i] = BCC;
            else {
              if (n0 > 12) structure[i] = UNKNOWN;
              else {
                if (delta_fcc < delta_hcp) structure[i] = FCC;
                else
                  structure[i] = HCP;
              }
            }
          }
        }
      }
    } else structure[i] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   2 select routines from Numerical Recipes (slightly modified)
   find k smallest values in array of length n
   2nd routine sorts auxiliary array at same time
------------------------------------------------------------------------- */

void ComputeAcklandAtom::select(int k, int n, double *arr)
  {
  int i,ir,j,l,mid;
  double a;

  arr--;
  l = 1;
  ir = n;
  while (true) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) std::swap(arr[l],arr[ir]);
      return;
    } else {
      mid=(l+ir) >> 1;
      std::swap(arr[mid],arr[l+1]);
      if (arr[l] > arr[ir]) std::swap(arr[l],arr[ir]);
      if (arr[l+1] > arr[ir]) std::swap(arr[l+1],arr[ir]);
      if (arr[l] > arr[l+1]) std::swap(arr[l],arr[l+1]);
      i = l+1;
      j = ir;
      a = arr[l+1];
      while (true) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        std::swap(arr[i],arr[j]);
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      if (j >= k) ir = j-1;
      if (j <= k) l = i;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeAcklandAtom::select2(int k, int n, double *arr, int *iarr)
{
  int i,ir,j,l,mid,ia;
  double a;

  arr--;
  iarr--;
  l = 1;
  ir = n;
  while (true) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
        std::swap(arr[l],arr[ir]);
        std::swap(iarr[l],iarr[ir]);
      }
      return;
    } else {
      mid=(l+ir) >> 1;
      std::swap(arr[mid],arr[l+1]);
      std::swap(iarr[mid],iarr[l+1]);
      if (arr[l] > arr[ir]) {
        std::swap(arr[l],arr[ir]);
        std::swap(iarr[l],iarr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
        std::swap(arr[l+1],arr[ir]);
        std::swap(iarr[l+1],iarr[ir]);
      }
      if (arr[l] > arr[l+1]) {
        std::swap(arr[l],arr[l+1]);
        std::swap(iarr[l],iarr[l+1]);
      }
      i = l+1;
      j = ir;
      a = arr[l+1];
      ia = iarr[l+1];
      while (true) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        std::swap(arr[i],arr[j]);
        std::swap(iarr[i],iarr[j]);
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      iarr[l+1] = iarr[j];
      iarr[j] = ia;
      if (j >= k) ir = j-1;
      if (j <= k) l = i;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeAcklandAtom::memory_usage()
{
  double bytes = (double)nmax * sizeof(double);
  return bytes;
}
