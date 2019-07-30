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
   Contributing author: Wengen Ouyang (Tel Aviv University)
   e-mail: w.g.ouyang at gmail dot com
   based on previous versions by Jaap Kroes

   This is a complete version of the potential described in
   [Kolmogorov & Crespi, Phys. Rev. B 71, 235415 (2005)]
------------------------------------------------------------------------- */

#include "pair_kolmogorov_crespi_full.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <mpi.h>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "my_page.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define DELTA 4
#define PGDELTA 1

/* ---------------------------------------------------------------------- */

PairKolmogorovCrespiFull::PairKolmogorovCrespiFull(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  one_coeff = 1;

  // initialize element to parameter maps
  nelements = 0;
  elements = NULL;
  nparams = maxparam = 0;
  params = NULL;
  elem2param = NULL;
  cutKCsq = NULL;
  map = NULL;

  nmax = 0;
  maxlocal = 0;
  KC_numneigh = NULL;
  KC_firstneigh = NULL;
  ipage = NULL;
  pgsize = oneatom = 0;

  normal = NULL;
  dnormal = NULL;
  dnormdri = NULL;

  // always compute energy offset
  offset_flag = 1;

  // set comm size needed by this Pair
  comm_forward = 39;
  tap_flag = 0;
}

/* ---------------------------------------------------------------------- */

PairKolmogorovCrespiFull::~PairKolmogorovCrespiFull()
{
  memory->destroy(KC_numneigh);
  memory->sfree(KC_firstneigh);
  delete [] ipage;
  memory->destroy(normal);
  memory->destroy(dnormal);
  memory->destroy(dnormdri);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(offset);
  }

  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  memory->destroy(params);
  memory->destroy(elem2param);
  memory->destroy(cutKCsq);
  if (allocated) delete [] map;
}

/* ---------------------------------------------------------------------- */

void PairKolmogorovCrespiFull::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,k,l,kk,ll;
  tagint itag,jtag;
  double prodnorm1,prodnorm2,fkcx,fkcy,fkcz;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair,fpair1,fpair2;
  double rsq,r,rhosq1,rhosq2,exp0,exp1,exp2,r2inv,r6inv,r8inv,Tap,dTap,Vkc;
  double frho1,frho2,sumC1,sumC2,sumC11,sumC22,sumCff,fsum,rdsq1,rdsq2;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *KC_neighs_i,*KC_neighs_j;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double dprodnorm1[3] = {0.0, 0.0, 0.0};
  double dprodnorm2[3] = {0.0, 0.0, 0.0};
  double fp1[3] = {0.0, 0.0, 0.0};
  double fp2[3] = {0.0, 0.0, 0.0};
  double fprod1[3] = {0.0, 0.0, 0.0};
  double fprod2[3] = {0.0, 0.0, 0.0};
  double fk[3] = {0.0, 0.0, 0.0};
  double fl[3] = {0.0, 0.0, 0.0};
  double delkj[3] = {0.0, 0.0, 0.0};
  double delli[3] = {0.0, 0.0, 0.0};

  inum = list->inum;	
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  // Build full neighbor list
  KC_neigh();
  // Calculate the normals
  calc_normal();

  // communicate the normal vector and its derivatives
  comm->forward_comm_pair(this);

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      jtag = tag[j];

      // two-body interactions from full neighbor list, skip half of them
      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      // only include the interation between different layers
      if (rsq < cutsq[itype][jtype] && atom->molecule[i] != atom->molecule[j]) {

        int iparam_ij = elem2param[map[itype]][map[jtype]];
        Param& p = params[iparam_ij];

        r = sqrt(rsq);
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        r8inv = r2inv*r6inv;
	// turn on/off taper function
	if (tap_flag) {
	  Tap = calc_Tap(r,sqrt(cutsq[itype][jtype]));
	  dTap = calc_dTap(r,sqrt(cutsq[itype][jtype]));
	} else {Tap = 1.0; dTap = 0.0;}

        // Calculate the transverse distance
        // note that rho_ij does not equal to rho_ji except when normals are all along z
        prodnorm1 = normal[i][0]*delx + normal[i][1]*dely + normal[i][2]*delz;
        prodnorm2 = normal[j][0]*delx + normal[j][1]*dely + normal[j][2]*delz;
        rhosq1 = rsq - prodnorm1*prodnorm1;  // rho_ij
        rhosq2 = rsq - prodnorm2*prodnorm2;  // rho_ji
        rdsq1 = rhosq1*p.delta2inv; // (rho_ij/delta)^2
        rdsq2 = rhosq2*p.delta2inv; // (rho_ji/delta)^2

        // store exponents
        exp0 = exp(-p.lambda*(r-p.z0));
        exp1 = exp(-rdsq1);
        exp2 = exp(-rdsq2);

        sumC1 = p.C0 + p.C2*rdsq1 + p.C4*rdsq1*rdsq1;
        sumC2 = p.C0 + p.C2*rdsq2 + p.C4*rdsq2*rdsq2;
        sumC11 = (p.C2 + 2.0*p.C4*rdsq1)*p.delta2inv;
        sumC22 = (p.C2 + 2.0*p.C4*rdsq2)*p.delta2inv;
        frho1 = exp1*sumC1;
        frho2 = exp2*sumC2;
        sumCff = p.C + frho1 + frho2;
	Vkc = -p.A*p.z06*r6inv + exp0*sumCff;

        // derivatives
        fpair = -6.0*p.A*p.z06*r8inv + p.lambda*exp0/r*sumCff;
        fpair1 = 2.0*exp0*exp1*(p.delta2inv*sumC1 - sumC11);
        fpair2 = 2.0*exp0*exp2*(p.delta2inv*sumC2 - sumC22);
        fsum = fpair + fpair1 + fpair2;
        // derivatives of the product of rij and ni, the result is a vector
        dprodnorm1[0] = dnormdri[0][0][i]*delx + dnormdri[1][0][i]*dely + dnormdri[2][0][i]*delz;
        dprodnorm1[1] = dnormdri[0][1][i]*delx + dnormdri[1][1][i]*dely + dnormdri[2][1][i]*delz;
        dprodnorm1[2] = dnormdri[0][2][i]*delx + dnormdri[1][2][i]*dely + dnormdri[2][2][i]*delz;
        // derivatives of the product of rji and nj, the result is a vector
        dprodnorm2[0] = dnormdri[0][0][j]*delx + dnormdri[1][0][j]*dely + dnormdri[2][0][j]*delz;
        dprodnorm2[1] = dnormdri[0][1][j]*delx + dnormdri[1][1][j]*dely + dnormdri[2][1][j]*delz;
        dprodnorm2[2] = dnormdri[0][2][j]*delx + dnormdri[1][2][j]*dely + dnormdri[2][2][j]*delz;
        fp1[0] = prodnorm1*normal[i][0]*fpair1;
        fp1[1] = prodnorm1*normal[i][1]*fpair1;
        fp1[2] = prodnorm1*normal[i][2]*fpair1;
        fp2[0] = prodnorm2*normal[j][0]*fpair2;
        fp2[1] = prodnorm2*normal[j][1]*fpair2;
        fp2[2] = prodnorm2*normal[j][2]*fpair2;
        fprod1[0] = prodnorm1*dprodnorm1[0]*fpair1;
        fprod1[1] = prodnorm1*dprodnorm1[1]*fpair1;
        fprod1[2] = prodnorm1*dprodnorm1[2]*fpair1;
        fprod2[0] = prodnorm2*dprodnorm2[0]*fpair2;
        fprod2[1] = prodnorm2*dprodnorm2[1]*fpair2;
        fprod2[2] = prodnorm2*dprodnorm2[2]*fpair2;
        fkcx = (delx*fsum - fp1[0] - fp2[0])*Tap - Vkc*dTap*delx/r;
        fkcy = (dely*fsum - fp1[1] - fp2[1])*Tap - Vkc*dTap*dely/r;
        fkcz = (delz*fsum - fp1[2] - fp2[2])*Tap - Vkc*dTap*delz/r;

        f[i][0] += fkcx - fprod1[0]*Tap;
        f[i][1] += fkcy - fprod1[1]*Tap;
        f[i][2] += fkcz - fprod1[2]*Tap;
        f[j][0] -= fkcx + fprod2[0]*Tap;
        f[j][1] -= fkcy + fprod2[1]*Tap;
        f[j][2] -= fkcz + fprod2[2]*Tap;

	// calculate the forces acted on the neighbors of atom i from atom j
	KC_neighs_i = KC_firstneigh[i];
	for (kk = 0; kk < KC_numneigh[i]; kk++) {
	  k = KC_neighs_i[kk];
          if (k == i) continue;
          // derivatives of the product of rij and ni respect to rk, k=0,1,2, where atom k is the neighbors of atom i
          dprodnorm1[0] = dnormal[0][0][kk][i]*delx + dnormal[1][0][kk][i]*dely + dnormal[2][0][kk][i]*delz;
          dprodnorm1[1] = dnormal[0][1][kk][i]*delx + dnormal[1][1][kk][i]*dely + dnormal[2][1][kk][i]*delz;
          dprodnorm1[2] = dnormal[0][2][kk][i]*delx + dnormal[1][2][kk][i]*dely + dnormal[2][2][kk][i]*delz;
          fk[0] = (-prodnorm1*dprodnorm1[0]*fpair1)*Tap;
          fk[1] = (-prodnorm1*dprodnorm1[1]*fpair1)*Tap;
          fk[2] = (-prodnorm1*dprodnorm1[2]*fpair1)*Tap;
          f[k][0] += fk[0];
          f[k][1] += fk[1];
          f[k][2] += fk[2];
          delkj[0] = x[k][0] - x[j][0];
          delkj[1] = x[k][1] - x[j][1];
          delkj[2] = x[k][2] - x[j][2];
          if (evflag) ev_tally_xyz(k,j,nlocal,newton_pair,0.0,0.0,fk[0],fk[1],fk[2],delkj[0],delkj[1],delkj[2]);
	}

	// calculate the forces acted on the neighbors of atom j from atom i
	KC_neighs_j = KC_firstneigh[j];
	for (ll = 0; ll < KC_numneigh[j]; ll++) {
	  l = KC_neighs_j[ll];
          if (l == j) continue;
          // derivatives of the product of rji and nj respect to rl, l=0,1,2, where atom l is the neighbors of atom j
          dprodnorm2[0] = dnormal[0][0][ll][j]*delx + dnormal[1][0][ll][j]*dely + dnormal[2][0][ll][j]*delz;
          dprodnorm2[1] = dnormal[0][1][ll][j]*delx + dnormal[1][1][ll][j]*dely + dnormal[2][1][ll][j]*delz;
          dprodnorm2[2] = dnormal[0][2][ll][j]*delx + dnormal[1][2][ll][j]*dely + dnormal[2][2][ll][j]*delz;
          fl[0] = (-prodnorm2*dprodnorm2[0]*fpair2)*Tap;
          fl[1] = (-prodnorm2*dprodnorm2[1]*fpair2)*Tap;
          fl[2] = (-prodnorm2*dprodnorm2[2]*fpair2)*Tap;
          f[l][0] += fl[0];
          f[l][1] += fl[1];
          f[l][2] += fl[2];
          delli[0] = x[l][0] - x[i][0];
          delli[1] = x[l][1] - x[i][1];
          delli[2] = x[l][2] - x[i][2];
          if (evflag) ev_tally_xyz(l,i,nlocal,newton_pair,0.0,0.0,fl[0],fl[1],fl[2],delli[0],delli[1],delli[2]);
	}

        if (eflag) {
          if (tap_flag) evdwl = Tap*Vkc;
          else  evdwl = Vkc - offset[itype][jtype];
        }

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,evdwl,0,fkcx,fkcy,fkcz,delx,dely,delz);
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   Calculate the normals for each atom
------------------------------------------------------------------------- */
void PairKolmogorovCrespiFull::calc_normal()
{
  int i,j,ii,jj,inum,jnum;
  int cont,id,ip,m;
  double nn,xtp,ytp,ztp,delx,dely,delz,nn2;
  int *ilist,*jlist;
  double pv12[3],pv31[3],pv23[3],n1[3],dni[3],dnn[3][3],vet[3][3],dpvdri[3][3];
  double dn1[3][3][3],dpv12[3][3][3],dpv23[3][3][3],dpv31[3][3][3];

  double **x = atom->x;

  // grow normal array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(normal);
    memory->destroy(dnormal);
    memory->destroy(dnormdri);
    nmax = atom->nmax;
    memory->create(normal,nmax,3,"KolmogorovCrespiFull:normal");
    memory->create(dnormdri,3,3,nmax,"KolmogorovCrespiFull:dnormdri");
    memory->create(dnormal,3,3,3,nmax,"KolmogorovCrespiFull:dnormal");
  }

  inum = list->inum;	
  ilist = list->ilist;
  //Calculate normals
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtp = x[i][0];
    ytp = x[i][1];
    ztp = x[i][2];

    //   Initialize the arrays
    for (id = 0; id < 3; id++){
      pv12[id] = 0.0;
      pv31[id] = 0.0;
      pv23[id] = 0.0;
      n1[id] = 0.0;
      dni[id] = 0.0;
      normal[i][id] = 0.0;
      for (ip = 0; ip < 3; ip++){
        vet[ip][id] = 0.0;
        dnn[ip][id] = 0.0;
        dpvdri[ip][id] = 0.0;
        dnormdri[ip][id][i] = 0.0;
        for (m = 0; m < 3; m++){
          dpv12[ip][id][m] = 0.0;
          dpv31[ip][id][m] = 0.0;
          dpv23[ip][id][m] = 0.0;
          dn1[ip][id][m] = 0.0;
          dnormal[ip][id][m][i] = 0.0;
        }
      }
    }

    cont = 0;
    jlist = KC_firstneigh[i];
    jnum = KC_numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = x[j][0] - xtp;
      dely = x[j][1] - ytp;
      delz = x[j][2] - ztp;
      vet[cont][0] = delx;
      vet[cont][1] = dely;
      vet[cont][2] = delz;
      cont++;
    }

    if (cont <= 1) {
      normal[i][0] = 0.0;
      normal[i][1] = 0.0;
      normal[i][2] = 1.0;
      // derivatives of normal vector is zero
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          dnormdri[id][ip][i] = 0.0;
          for (m = 0; m < 3; m++){
            dnormal[id][ip][m][i] = 0.0;
          }
        }
      }
    }
    else if (cont == 2) {
      // for the atoms at the edge who has only two neighbor atoms
      pv12[0] = vet[0][1]*vet[1][2] - vet[1][1]*vet[0][2];
      pv12[1] = vet[0][2]*vet[1][0] - vet[1][2]*vet[0][0];
      pv12[2] = vet[0][0]*vet[1][1] - vet[1][0]*vet[0][1];
      dpvdri[0][0] = 0.0;
      dpvdri[0][1] = vet[0][2]-vet[1][2];
      dpvdri[0][2] = vet[1][1]-vet[0][1];
      dpvdri[1][0] = vet[1][2]-vet[0][2];
      dpvdri[1][1] = 0.0;
      dpvdri[1][2] = vet[0][0]-vet[1][0];
      dpvdri[2][0] = vet[0][1]-vet[1][1];
      dpvdri[2][1] = vet[1][0]-vet[0][0];
      dpvdri[2][2] = 0.0;

      // derivatives respect to the first neighbor, atom k
      dpv12[0][0][0] =  0.0;
      dpv12[0][1][0] =  vet[1][2];
      dpv12[0][2][0] = -vet[1][1];
      dpv12[1][0][0] = -vet[1][2];
      dpv12[1][1][0] =  0.0;
      dpv12[1][2][0] =  vet[1][0];
      dpv12[2][0][0] =  vet[1][1];
      dpv12[2][1][0] = -vet[1][0];
      dpv12[2][2][0] =  0.0;

      // derivatives respect to the second neighbor, atom l
      dpv12[0][0][1] =  0.0;
      dpv12[0][1][1] = -vet[0][2];
      dpv12[0][2][1] =  vet[0][1];
      dpv12[1][0][1] =  vet[0][2];
      dpv12[1][1][1] =  0.0;
      dpv12[1][2][1] = -vet[0][0];
      dpv12[2][0][1] = -vet[0][1];
      dpv12[2][1][1] =  vet[0][0];
      dpv12[2][2][1] =  0.0;

      // derivatives respect to the third neighbor, atom n
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          dpv12[id][ip][2] = 0.0;
        }
      }

      n1[0] = pv12[0];
      n1[1] = pv12[1];
      n1[2] = pv12[2];
      // the magnitude of the normal vector
      nn2 = n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2];
      nn = sqrt(nn2);
      if (nn == 0) error->one(FLERR,"The magnitude of the normal vector is zero");
      // the unit normal vector
      normal[i][0] = n1[0]/nn;
      normal[i][1] = n1[1]/nn;
      normal[i][2] = n1[2]/nn;
      // derivatives of nn, dnn:3x1 vector
      dni[0] = (n1[0]*dpvdri[0][0] + n1[1]*dpvdri[1][0] + n1[2]*dpvdri[2][0])/nn;
      dni[1] = (n1[0]*dpvdri[0][1] + n1[1]*dpvdri[1][1] + n1[2]*dpvdri[2][1])/nn;
      dni[2] = (n1[0]*dpvdri[0][2] + n1[1]*dpvdri[1][2] + n1[2]*dpvdri[2][2])/nn;
      // derivatives of unit vector ni respect to ri, the result is 3x3 matrix
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          dnormdri[id][ip][i] = dpvdri[id][ip]/nn - n1[id]*dni[ip]/nn2;
        }
      }

      // derivatives of non-normalized normal vector, dn1:3x3x3 array
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          for (m = 0; m < 3; m++){
            dn1[id][ip][m] = dpv12[id][ip][m];
          }
        }
      }
      // derivatives of nn, dnn:3x3 vector
      // dnn[id][m]: the derivative of nn respect to r[id][m], id,m=0,1,2
      // r[id][m]: the id's component of atom m
      for (m = 0; m < 3; m++){
        for (id = 0; id < 3; id++){
          dnn[id][m] = (n1[0]*dn1[0][id][m] + n1[1]*dn1[1][id][m] + n1[2]*dn1[2][id][m])/nn;
        }
      }
      // dnormal[id][ip][m][i]: the derivative of normal[id] respect to r[ip][m], id,ip=0,1,2
      // for atom m, which is a neighbor atom of atom i, m=0,jnum-1
      for (m = 0; m < 3; m++){
        for (id = 0; id < 3; id++){
          for (ip = 0; ip < 3; ip++){
            dnormal[id][ip][m][i] = dn1[id][ip][m]/nn - n1[id]*dnn[ip][m]/nn2;
          }
        }
      }
    }
//##############################################################################################

    else if(cont == 3) {
      // for the atoms at the edge who has only two neighbor atoms
      pv12[0] = vet[0][1]*vet[1][2] - vet[1][1]*vet[0][2];
      pv12[1] = vet[0][2]*vet[1][0] - vet[1][2]*vet[0][0];
      pv12[2] = vet[0][0]*vet[1][1] - vet[1][0]*vet[0][1];
      // derivatives respect to the first neighbor, atom k
      dpv12[0][0][0] =  0.0;
      dpv12[0][1][0] =  vet[1][2];
      dpv12[0][2][0] = -vet[1][1];
      dpv12[1][0][0] = -vet[1][2];
      dpv12[1][1][0] =  0.0;
      dpv12[1][2][0] =  vet[1][0];
      dpv12[2][0][0] =  vet[1][1];
      dpv12[2][1][0] = -vet[1][0];
      dpv12[2][2][0] =  0.0;
      // derivatives respect to the second neighbor, atom l
      dpv12[0][0][1] =  0.0;
      dpv12[0][1][1] = -vet[0][2];
      dpv12[0][2][1] =  vet[0][1];
      dpv12[1][0][1] =  vet[0][2];
      dpv12[1][1][1] =  0.0;
      dpv12[1][2][1] = -vet[0][0];
      dpv12[2][0][1] = -vet[0][1];
      dpv12[2][1][1] =  vet[0][0];
      dpv12[2][2][1] =  0.0;

      // derivatives respect to the third neighbor, atom n
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          dpv12[id][ip][2] = 0.0;
        }
      }

      pv31[0] = vet[2][1]*vet[0][2] - vet[0][1]*vet[2][2];
      pv31[1] = vet[2][2]*vet[0][0] - vet[0][2]*vet[2][0];
      pv31[2] = vet[2][0]*vet[0][1] - vet[0][0]*vet[2][1];
      // derivatives respect to the first neighbor, atom k
      dpv31[0][0][0] =  0.0;
      dpv31[0][1][0] = -vet[2][2];
      dpv31[0][2][0] =  vet[2][1];
      dpv31[1][0][0] =  vet[2][2];
      dpv31[1][1][0] =  0.0;
      dpv31[1][2][0] = -vet[2][0];
      dpv31[2][0][0] = -vet[2][1];
      dpv31[2][1][0] =  vet[2][0];
      dpv31[2][2][0] =  0.0;
      // derivatives respect to the third neighbor, atom n
      dpv31[0][0][2] =  0.0;
      dpv31[0][1][2] =  vet[0][2];
      dpv31[0][2][2] = -vet[0][1];
      // derivatives of pv13[1] to rn
      dpv31[1][0][2] = -vet[0][2];
      dpv31[1][1][2] =  0.0;
      dpv31[1][2][2] =  vet[0][0];
      // derivatives of pv13[2] to rn
      dpv31[2][0][2] =  vet[0][1];
      dpv31[2][1][2] = -vet[0][0];
      dpv31[2][2][2] =  0.0;

      // derivatives respect to the second neighbor, atom l
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          dpv31[id][ip][1] = 0.0;
        }
      }

      pv23[0] = vet[1][1]*vet[2][2] - vet[2][1]*vet[1][2];
      pv23[1] = vet[1][2]*vet[2][0] - vet[2][2]*vet[1][0];
      pv23[2] = vet[1][0]*vet[2][1] - vet[2][0]*vet[1][1];
      // derivatives respect to the second neighbor, atom k
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          dpv23[id][ip][0] = 0.0;
        }
      }
      // derivatives respect to the second neighbor, atom l
      dpv23[0][0][1] =  0.0;
      dpv23[0][1][1] =  vet[2][2];
      dpv23[0][2][1] = -vet[2][1];
      dpv23[1][0][1] = -vet[2][2];
      dpv23[1][1][1] =  0.0;
      dpv23[1][2][1] =  vet[2][0];
      dpv23[2][0][1] =  vet[2][1];
      dpv23[2][1][1] = -vet[2][0];
      dpv23[2][2][1] =  0.0;
      // derivatives respect to the third neighbor, atom n
      dpv23[0][0][2] =  0.0;
      dpv23[0][1][2] = -vet[1][2];
      dpv23[0][2][2] =  vet[1][1];
      dpv23[1][0][2] =  vet[1][2];
      dpv23[1][1][2] =  0.0;
      dpv23[1][2][2] = -vet[1][0];
      dpv23[2][0][2] = -vet[1][1];
      dpv23[2][1][2] =  vet[1][0];
      dpv23[2][2][2] =  0.0;

//############################################################################################
      // average the normal vectors by using the 3 neighboring planes
      n1[0] = (pv12[0] + pv31[0] + pv23[0])/cont;
      n1[1] = (pv12[1] + pv31[1] + pv23[1])/cont;
      n1[2] = (pv12[2] + pv31[2] + pv23[2])/cont;
      // the magnitude of the normal vector
      nn2 = n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2];
      nn = sqrt(nn2);
      if (nn == 0) error->one(FLERR,"The magnitude of the normal vector is zero");
      // the unit normal vector
      normal[i][0] = n1[0]/nn;
      normal[i][1] = n1[1]/nn;
      normal[i][2] = n1[2]/nn;

      // for the central atoms, dnormdri is always zero
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          dnormdri[id][ip][i] = 0.0;
        }
      } // end of derivatives of normals respect to atom i

      // derivatives of non-normalized normal vector, dn1:3x3x3 array
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
          for (m = 0; m < 3; m++){
            dn1[id][ip][m] = (dpv12[id][ip][m] + dpv23[id][ip][m] + dpv31[id][ip][m])/cont;
          }
        }
      }
      // derivatives of nn, dnn:3x3 vector
      // dnn[id][m]: the derivative of nn respect to r[id][m], id,m=0,1,2
      // r[id][m]: the id's component of atom m
      for (m = 0; m < 3; m++){
        for (id = 0; id < 3; id++){
          dnn[id][m] = (n1[0]*dn1[0][id][m] + n1[1]*dn1[1][id][m] + n1[2]*dn1[2][id][m])/nn;
        }
      }
      // dnormal[id][ip][m][i]: the derivative of normal[id] respect to r[ip][m], id,ip=0,1,2
      // for atom m, which is a neighbor atom of atom i, m=0,jnum-1
      for (m = 0; m < 3; m++){
        for (id = 0; id < 3; id++){
          for (ip = 0; ip < 3; ip++){
            dnormal[id][ip][m][i] = dn1[id][ip][m]/nn - n1[id]*dnn[ip][m]/nn2;
          }
        }
      }
    }
    else {
      error->one(FLERR,"There are too many neighbors for calculating normals");
    }

//##############################################################################################
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairKolmogorovCrespiFull::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style kolmolgorov/crespi/full requires newton pair on");
  if (!atom->molecule_flag)
    error->all(FLERR,"Pair style kolmolgorov/crespi/full requires atom attribute molecule");

  // need a full neighbor list, including neighbors of ghosts

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->ghost = 1;

  // local KC neighbor list
  // create pages if first time or if neighbor pgsize/oneatom has changed

  int create = 0;
  if (ipage == NULL) create = 1;
  if (pgsize != neighbor->pgsize) create = 1;
  if (oneatom != neighbor->oneatom) create = 1;

  if (create) {
    delete [] ipage;
    pgsize = neighbor->pgsize;
    oneatom = neighbor->oneatom;

    int nmypage= comm->nthreads;
    ipage = new MyPage<int>[nmypage];
    for (int i = 0; i < nmypage; i++)
      ipage[i].init(oneatom,pgsize,PGDELTA);
  }
}


/* ----------------------------------------------------------------------
 create neighbor list from main neighbor list for calculating the normals
------------------------------------------------------------------------- */

void PairKolmogorovCrespiFull::KC_neigh()
{
  int i,j,ii,jj,n,allnum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *neighptr;

  double **x = atom->x;
  int *type = atom->type;

  if (atom->nmax > maxlocal) {
    maxlocal = atom->nmax;
    memory->destroy(KC_numneigh);
    memory->sfree(KC_firstneigh);
    memory->create(KC_numneigh,maxlocal,"KolmogorovCrespiFull:numneigh");
    KC_firstneigh = (int **) memory->smalloc(maxlocal*sizeof(int *),
                                             "KolmogorovCrespiFull:firstneigh");
  }

  allnum = list->inum + list->gnum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // store all KC neighs of owned and ghost atoms
  // scan full neighbor list of I

  ipage->reset();

  for (ii = 0; ii < allnum; ii++) {
    i = ilist[ii];

    n = 0;
    neighptr = ipage->vget();

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = map[type[i]];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq != 0 && rsq < cutKCsq[itype][jtype] && atom->molecule[i] == atom->molecule[j]) {
        neighptr[n++] = j;
      }
    }

    KC_firstneigh[i] = neighptr;
    KC_numneigh[i] = n;
    if (n > 3) error->one(FLERR,"There are too many neighbors for some atoms, please check your configuration");
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairKolmogorovCrespiFull::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(offset,n+1,n+1,"pair:offset");
  map = new int[atom->ntypes+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairKolmogorovCrespiFull::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2) error->all(FLERR,"Illegal pair_style command");
  if (strcmp(force->pair_style,"hybrid/overlay")!=0)
    error->all(FLERR,"ERROR: requires hybrid/overlay pair_style");

  cut_global = force->numeric(FLERR,arg[0]);
  if (narg == 2) tap_flag = force->numeric(FLERR,arg[1]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairKolmogorovCrespiFull::coeff(int narg, char **arg)
{
  int i,j,n;

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) elements[i] = NULL;

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }


  read_file(arg[2]);

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        cut[i][j] = cut_global;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairKolmogorovCrespiFull::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  if (!offset_flag)
    error->all(FLERR,"Must use 'pair_modify shift yes' with this pair style");

  if (offset_flag && (cut[i][j] > 0.0)) {
    int iparam_ij = elem2param[map[i]][map[j]];
    Param& p = params[iparam_ij];
    offset[i][j] = -p.A*pow(p.z0/cut[i][j],6);
  } else offset[i][j] = 0.0;
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   read Kolmogorov-Crespi potential file
------------------------------------------------------------------------- */

void PairKolmogorovCrespiFull::read_file(char *filename)
{
  int params_per_line = 12;
  char **words = new char*[params_per_line+1];
  memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential(filename);
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open KC potential file %s",filename);
      error->one(FLERR,str);
    }
  }

  // read each line out of file, skipping blank lines or leading '#'
  // store line of params if all 3 element tags are in element list

  int i,j,n,m,nwords,ielement,jelement;
  char line[MAXLINE],*ptr;
  int eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR,"Insufficient format in KC potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    // ielement,jelement = 1st args
    // if these 2 args are in element list, then parse this line
    // else skip to next line (continue)

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                          "pair:params");
    }

    params[nparams].ielement = ielement;
    params[nparams].jelement = jelement;
    params[nparams].z0       = atof(words[2]);
    params[nparams].C0       = atof(words[3]);
    params[nparams].C2       = atof(words[4]);
    params[nparams].C4       = atof(words[5]);
    params[nparams].C        = atof(words[6]);
    params[nparams].delta    = atof(words[7]);
    params[nparams].lambda   = atof(words[8]);
    params[nparams].A        = atof(words[9]);
    // S provides a convenient scaling of all energies
    params[nparams].S        = atof(words[10]);
    params[nparams].rcut     = atof(words[11]);

    // energies in meV further scaled by S
    double meV = 1.0e-3*params[nparams].S;
    params[nparams].C *= meV;
    params[nparams].A *= meV;
    params[nparams].C0 *= meV;
    params[nparams].C2 *= meV;
    params[nparams].C4 *= meV;

    // precompute some quantities
    params[nparams].delta2inv = pow(params[nparams].delta,-2);
    params[nparams].z06 = pow(params[nparams].z0,6);

    nparams++;
    //if(nparams >= pow(atom->ntypes,3)) break;
  }
  memory->destroy(elem2param);
  memory->destroy(cutKCsq);
  memory->create(elem2param,nelements,nelements,"pair:elem2param");
  memory->create(cutKCsq,nelements,nelements,"pair:cutKCsq");
  for (i = 0; i < nelements; i++) {
    for (j = 0; j < nelements; j++) {
      n = -1;
      for (m = 0; m < nparams; m++) {
        if (i == params[m].ielement && j == params[m].jelement) {
          if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
          n = m;
        }
      }
      if (n < 0) error->all(FLERR,"Potential file is missing an entry");
      elem2param[i][j] = n;
      cutKCsq[i][j] = params[n].rcut*params[n].rcut;
    }
  }
  delete [] words;
}

/* ---------------------------------------------------------------------- */

double PairKolmogorovCrespiFull::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                         double /*factor_coul*/, double factor_lj,
                         double &fforce)
{
  double r,r2inv,r6inv,r8inv,forcelj,philj;
  double Tap,dTap,Vkc,fpair;

  int iparam_ij = elem2param[map[itype]][map[jtype]];
  Param& p = params[iparam_ij];

  r = sqrt(rsq);
  // turn on/off taper function
  if (tap_flag) {
    Tap = calc_Tap(r,sqrt(cutsq[itype][jtype]));
    dTap = calc_dTap(r,sqrt(cutsq[itype][jtype]));
  } else {Tap = 1.0; dTap = 0.0;}

  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;
  r8inv = r2inv*r6inv;

  Vkc = -p.A*p.z06*r6inv;
  // derivatives
  fpair = -6.0*p.A*p.z06*r8inv;
  forcelj = fpair;
  fforce = factor_lj*(forcelj*Tap - Vkc*dTap/r);

  if (tap_flag) philj = Vkc*Tap;
  else philj = Vkc - offset[itype][jtype];
  return factor_lj*philj;
}

/* ---------------------------------------------------------------------- */

int PairKolmogorovCrespiFull::pack_forward_comm(int n, int *list, double *buf,
                               int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m,l,ip,id;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = normal[j][0];
    buf[m++] = normal[j][1];
    buf[m++] = normal[j][2];
    buf[m++] = dnormdri[0][0][j];
    buf[m++] = dnormdri[0][1][j];
    buf[m++] = dnormdri[0][2][j];
    buf[m++] = dnormdri[1][0][j];
    buf[m++] = dnormdri[1][1][j];
    buf[m++] = dnormdri[1][2][j];
    buf[m++] = dnormdri[2][0][j];
    buf[m++] = dnormdri[2][1][j];
    buf[m++] = dnormdri[2][2][j];
    for (l = 0; l < 3; l++){
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
	  buf[m++] = dnormal[id][ip][l][j];
        }
      }
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairKolmogorovCrespiFull::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last,l,ip,id;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    normal[i][0] = buf[m++];
    normal[i][1] = buf[m++];
    normal[i][2] = buf[m++];
    dnormdri[0][0][i] = buf[m++];
    dnormdri[0][1][i] = buf[m++];
    dnormdri[0][2][i] = buf[m++];
    dnormdri[1][0][i] = buf[m++];
    dnormdri[1][1][i] = buf[m++];
    dnormdri[1][2][i] = buf[m++];
    dnormdri[2][0][i] = buf[m++];
    dnormdri[2][1][i] = buf[m++];
    dnormdri[2][2][i] = buf[m++];
    for (l = 0; l < 3; l++){
      for (id = 0; id < 3; id++){
        for (ip = 0; ip < 3; ip++){
	  dnormal[id][ip][l][i] = buf[m++];
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
