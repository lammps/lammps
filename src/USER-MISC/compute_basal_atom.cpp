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
   Contributing author: C.D. Barrett, cdb333@cavs.msstate.edu
                        Copyright (C) 2013
------------------------------------------------------------------------- */


#include "string.h"
#include "compute_basal_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include <math.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeBasalAtom::ComputeBasalAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute basal/atom command");

  peratom_flag = 1;
  size_peratom_cols = 3;

  nmax = 0;
  BPV = NULL;
  maxneigh = 0;
  distsq = NULL;
  nearest = NULL;
  nearest_n0 = NULL;
  nearest_n1 = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeBasalAtom::~ComputeBasalAtom()
{
  memory->destroy(BPV);
  memory->destroy(distsq);
  memory->destroy(nearest);
  memory->destroy(nearest_n0);
  memory->destroy(nearest_n1);
}

/* ---------------------------------------------------------------------- */

void ComputeBasalAtom::init()
{
  // need an occasional full neighbor list

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count1 = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"basal/atom") == 0) count1++;
  if (count1 > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute basal/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeBasalAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeBasalAtom::compute_peratom()
{
  int i,j,ii,jj,k,n,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,var5,var6,var7;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int chi[8];
  int value;
  int count;
  int k2[3];
  int j1[3];
  double x4[3],y4[3],z4[3],x5[3],y5[3],z5[3],x6[3],y6[3],z6[3];
  double x7[3],y7[3],z7[3];

  invoked_peratom = update->ntimestep;

  // grow structure array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(BPV);
    nmax = atom->nmax;
    memory->create(BPV,nmax,3,"basal/atom:basal");
    array_atom = BPV;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list->index);

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
      	memory->create(distsq,maxneigh,"compute/basal/atom:distsq");
      	memory->create(nearest,maxneigh,"compute/basal/atom:nearest");
	memory->create(nearest_n0,maxneigh,"compute/basal/atom:nearest_n0");
	memory->create(nearest_n1,maxneigh,"compute/basal/atom:nearest_n1");
      }
      // neighbor selection is identical to ackland/atom algorithm

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

      double r0_sq = 0.0;
      for (j = 0; j < 6; j++) r0_sq += distsq[j];
      r0_sq /= 6.0;

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

      // Evaluate all angles <(r_ij,rik) forall n0 particles with: distsq<1.45*r0_sq
      double bond_angle;
      double norm_j, norm_k;
      chi[0] = chi[1] = chi[2] = chi[3] = chi[4] = chi[5] = chi[6] = chi[7] = 0;
      double x_ij, y_ij, z_ij, x_ik, y_ik, z_ik,x3[n0],y3[n0],z3[n0],
        xmean5, ymean5, zmean5, xmean6, ymean6, zmean6, xmean7, ymean7, zmean7;
      for (j = 0; j < n0; j++) {
	x_ij = x[i][0]-x[nearest_n0[j]][0];
	y_ij = x[i][1]-x[nearest_n0[j]][1];
	z_ij = x[i][2]-x[nearest_n0[j]][2];
	norm_j = sqrt (x_ij*x_ij + y_ij*y_ij + z_ij*z_ij);
	if (norm_j <= 0.) {continue;}
	for (k = j+1; k < n0; k++) {
	  x_ik = x[i][0]-x[nearest_n0[k]][0];
	  y_ik = x[i][1]-x[nearest_n0[k]][1];
	  z_ik = x[i][2]-x[nearest_n0[k]][2];
	  norm_k = sqrt (x_ik*x_ik + y_ik*y_ik + z_ik*z_ik);
	  if (norm_k <= 0.) {continue;}
	  bond_angle = (x_ij*x_ik + y_ij*y_ik + z_ij*z_ik) / (norm_j*norm_k);
	  //find all bond angles that are about 180 degrees
	  if (-1. <= bond_angle && bond_angle < -0.945) { 
		x3[chi[0]] = x_ik - x_ij;
		y3[chi[0]] = y_ik - y_ij;
		z3[chi[0]] = z_ik - z_ij;
                chi[0]++;
 	  }
	}
      }
      // for atoms that have 2 or 3 ~180 bond angles:
      if (2 == chi[0] || 3 == chi[0]) {
          count = value = 0;
      	  if (chi[0] == 2) {
            k2[0] = 0;
            j1[0] = 1;
          }
          else {
            k2[0] = 0;
            k2[1] = 0;
            k2[2] = 1;
            j1[0]=1;
            j1[1]=2;
            j1[2]=2;
          }
          xmean5 = ymean5 = zmean5 = xmean6 = ymean6 = zmean6 = xmean7 = ymean7 = zmean7 = 0;
	  for (j = 0; j < chi[0]; j++) {
            for (k = j+1; k < chi[0]; k++) {
	       //get cross products
               x4[count] = y3[j1[count]]*z3[k2[count]]-y3[k2[count]]*z3[j1[count]];
               y4[count] = z3[j1[count]]*x3[k2[count]]-z3[k2[count]]*x3[j1[count]];
               z4[count] = x3[j1[count]]*y3[k2[count]]-x3[k2[count]]*y3[j1[count]];
	       //get all sign combinations of cross products
               x5[count] = x4[count]*copysign(1.0,x4[count]);
               y5[count] = y4[count]*copysign(1.0,x4[count]);
               z5[count] = z4[count]*copysign(1.0,x4[count]);
               x6[count] = x4[count]*copysign(1.0,y4[count]);
               y6[count] = y4[count]*copysign(1.0,y4[count]);
               z6[count] = z4[count]*copysign(1.0,y4[count]);
               x7[count] = x4[count]*copysign(1.0,z4[count]);
               y7[count] = y4[count]*copysign(1.0,z4[count]);
               z7[count] = z4[count]*copysign(1.0,z4[count]);
	       //get average cross products
               xmean5 = xmean5 + x5[count];
               ymean5 = ymean5 + y5[count];
               zmean5 = zmean5 + z5[count];
               xmean6 = xmean6 + x6[count];
               ymean6 = ymean6 + y6[count];
               zmean6 = zmean6 + z6[count];
               xmean7 = xmean7 + x7[count];
               ymean7 = ymean7 + y7[count];
               zmean6 = zmean6 + z7[count];
               count++;
            }
          }
          xmean5 = xmean5/count;
          xmean6 = xmean6/count;
          xmean7 = xmean7/count;
          ymean5 = ymean5/count;
          ymean6 = ymean6/count;
          ymean7 = ymean7/count;
          zmean5 = zmean5/count;
          zmean6 = zmean6/count;
          zmean7 = zmean7/count;
          var5 = var6 = var7 = 0.0;
	  //find standard deviations
          for (j=0;j<count;j++){
            var5 = var5 + x5[j]*x5[j]-2*x5[j]*xmean5+xmean5*xmean5+y5[j]*y5[j]-2*y5[j]*ymean5+ymean5*ymean5+z5[j]*z5[j]-2*z5[j]*zmean5+zmean5*zmean5;
            var6 = var6 + x6[j]*x6[j]-2*x6[j]*xmean6+xmean6*xmean6+y6[j]*y6[j]-2*y6[j]*ymean6+ymean6*ymean6+z6[j]*z6[j]-2*z6[j]*zmean6+zmean6*zmean6;
            var7 = var7 + x7[j]*x7[j]-2*x7[j]*xmean7+xmean7*xmean7+y7[j]*y7[j]-2*y7[j]*ymean7+ymean7*ymean7+z7[j]*z7[j]-2*z7[j]*zmean7+zmean7*zmean7;
          }
          //select sign combination with minimum standard deviation
          if (var5 < var6) {
              if (var5 < var7) { value = 0;}
              else {value = 2;}
          }
          else if (var6 < var7) {value = 1;}
          else {value = 2;}
	  //BPV is average of cross products of all neighbor vectors which are part of 180 degree angles
          BPV[i][0] = 0;
          BPV[i][1] = 0;
          BPV[i][2] = 0;
          for (k=0;k<count;k++) {
           if (value == 0){
               BPV[i][0] = BPV[i][0]+x5[k];
               BPV[i][1] = BPV[i][1]+y5[k];
               BPV[i][2] = BPV[i][2]+z5[k];
           }
           else if (value == 1) {
               BPV[i][0] = BPV[i][0]+x6[k];
               BPV[i][1] = BPV[i][1]+y6[k];
               BPV[i][2] = BPV[i][2]+z6[k];
           }
           else {
               BPV[i][0] = BPV[i][0]+x7[k];
               BPV[i][1] = BPV[i][1]+y7[k];
               BPV[i][2] = BPV[i][2]+z7[k];
           }
          }
      }
      //for atoms with more than three 180 degree bond angles:
      else if (chi[0] > 3) {
          double x44[3], y44[3], z44[3], S0;
          int l, m;
          count = value = 0;
          S0 = 100000;
          k2[0] = 0;
          k2[1] = 0;
          k2[2] = 1;
          j1[0]=1;
          j1[1]=2;
          j1[2]=2;
	  //algorithm is as above, but now all combinations of three 180 degree angles are compared, and the combination with minimum standard deviation is chosen
          for (j=0; j<chi[0]; j++) {
              for (k=j+1; k<chi[0]; k++) {
                  for (l=k+1; l<chi[0]; l++) {
                      if (k >= chi[0] || l >= chi[0]) continue;
		      //get unique combination of three neighbor vectors
                      x4[0] = x3[j];
                      x4[1] = x3[k];
                      x4[2] = x3[l];
                      y4[0] = y3[j];
                      y4[1] = y3[k];
                      y4[2] = y3[l];
                      z4[0] = z3[j];
                      z4[1] = z3[k];
                      z4[2] = z3[l];
                      xmean5 = ymean5 = zmean5 = xmean6 = ymean6 = zmean6 = xmean7 = ymean7 = zmean7 = 0;
                      for (m=0;m<3;m++) {
			//get cross products
                        x44[m] = y4[j1[m]]*z4[k2[m]]-y4[k2[m]]*z4[j1[m]];
                        y44[m] = z4[j1[m]]*x4[k2[m]]-z4[k2[m]]*x4[j1[m]];
                        z44[m] = x4[j1[m]]*y4[k2[m]]-x4[k2[m]]*y4[j1[m]];
                        x5[m] = x44[m]*copysign(1.0,x44[m]);
                        y5[m] = y44[m]*copysign(1.0,x44[m]);
                        z5[m] = z44[m]*copysign(1.0,x44[m]);
                        x6[m] = x44[m]*copysign(1.0,y44[m]);
                        y6[m] = y44[m]*copysign(1.0,y44[m]);
                        z6[m] = z44[m]*copysign(1.0,y44[m]);
                        x7[m] = x44[m]*copysign(1.0,z44[m]);
                        y7[m] = y44[m]*copysign(1.0,z44[m]);
                        z7[m] = z44[m]*copysign(1.0,z44[m]);
			//get average cross products
                        xmean5 = xmean5 + x5[m];
                        ymean5 = ymean5 + y5[m];
                        zmean5 = zmean5 + z5[m];
                        xmean6 = xmean6 + x6[m];
                        ymean6 = ymean6 + y6[m];
                        zmean6 = zmean6 + z6[m];
                        xmean7 = xmean7 + x7[m];
                        ymean7 = ymean7 + y7[m];
                        zmean6 = zmean6 + z7[m];
                      }
                      xmean5 = xmean5/3;
                      xmean6 = xmean6/3;
                      xmean7 = xmean7/3;
                      ymean5 = ymean5/3;
                      ymean6 = ymean6/3;
                      ymean7 = ymean7/3;
                      zmean5 = zmean5/3;
                      zmean6 = zmean6/3;
                      zmean7 = zmean7/3;
                      var5 = var6 = var7 = 0;
		      //get standard deviations
                      for (m=0;m<3;m++){
                            var5 = var5 + x5[m]*x5[m]-2*x5[m]*xmean5+xmean5*xmean5+y5[m]*y5[m]-2*y5[m]*ymean5+ymean5*ymean5+z5[m]*z5[m]-2*z5[m]*zmean5+zmean5*zmean5;
                            var6 = var6 + x6[m]*x6[m]-2*x6[m]*xmean6+xmean6*xmean6+y6[m]*y6[m]-2*y6[m]*ymean6+ymean6*ymean6+z6[m]*z6[m]-2*z6[m]*zmean6+zmean6*zmean6;
                            var7 = var7 + x7[m]*x7[m]-2*x7[m]*xmean7+xmean7*xmean7+y7[m]*y7[m]-2*y7[m]*ymean7+ymean7*ymean7+z7[m]*z7[m]-2*z7[m]*zmean7+zmean7*zmean7;
                      }
		      //choose minimum standard deviation
                      if (var5 < S0) {
                          S0 = var5;
                          BPV[i][0] = (x5[0]+x5[1]+x5[2])/3;
                          BPV[i][1] = (y5[0]+y5[1]+x5[2])/3;
                          BPV[i][2] = (z5[0]+z5[1]+z5[2])/3;
                      }
                      if (var6 < S0) {
                          S0 = var6;
                          BPV[i][0] = (x6[0]+x6[1]+x6[2])/3;
                          BPV[i][1] = (y6[0]+y6[1]+x6[2])/3;
                          BPV[i][2] = (z6[0]+z6[1]+z6[2])/3;
                      }
                      if (var7 < S0) {
                          S0 = var7;
                          BPV[i][0] = (x7[0]+x7[1]+x7[2])/3;
                          BPV[i][1] = (y7[0]+y7[1]+x7[2])/3;
                          BPV[i][2] = (z7[0]+z7[1]+z7[2])/3;
                      }
                  }
              }
          }
      }
      //if there are less than two ~180 degree bond angles, the algorithm returns null
      else BPV[i][0] = BPV[i][1] = BPV[i][2] = 0.0;

      //normalize BPV:
      double Mag = sqrt(BPV[i][0]*BPV[i][0] + 
                        BPV[i][1]*BPV[i][1] + BPV[i][2]*BPV[i][2]);
      if (Mag > 0){
        BPV[i][0] = BPV[i][0]/Mag;
        BPV[i][1] = BPV[i][1]/Mag;
        BPV[i][2] = BPV[i][2]/Mag;
      }
    } else BPV[i][0] = BPV[i][1] = BPV[i][2] = 0.0;
  }
}
/* ----------------------------------------------------------------------
   2 select routines from Numerical Recipes (slightly modified)
   find k smallest values in array of length n
   2nd routine sorts auxiliary array at same time
------------------------------------------------------------------------- */

#define SWAP(a,b)   tmp = a; a = b; b = tmp;
#define ISWAP(a,b) itmp = a; a = b; b = itmp;

void ComputeBasalAtom::select(int k, int n, double *arr)
  {
  int i,ir,j,l,mid;
  double a,tmp;

  arr--;
  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir])
      }
      return;
    } else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1])
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir])
      }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir])
      }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1])
      }
      i = l+1;
      j = ir;
      a = arr[l+1];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j])
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      if (j >= k) ir = j-1;
      if (j <= k) l = i;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeBasalAtom::select2(int k, int n, double *arr, int *iarr)
{
  int i,ir,j,l,mid,ia,itmp;
  double a,tmp;

  arr--;
  iarr--;
  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir])
	ISWAP(iarr[l],iarr[ir])
      }
      return;
    } else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1])
      ISWAP(iarr[mid],iarr[l+1])
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir])
	ISWAP(iarr[l],iarr[ir])
      }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir])
	ISWAP(iarr[l+1],iarr[ir])
      }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1])
	ISWAP(iarr[l],iarr[l+1])
      }
      i = l+1;
      j = ir;
      a = arr[l+1];
      ia = iarr[l+1];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j])
	ISWAP(iarr[i],iarr[j])
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

double ComputeBasalAtom::memory_usage()
{
  double bytes = 3*nmax * sizeof(double);
  return bytes;
}
