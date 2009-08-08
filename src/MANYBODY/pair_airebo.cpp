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
   Contributing author: Ase Henry (MIT)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include "pair_airebo.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define MAXLINE 1024
#define TOL 1.0e-9
#define PI 3.14159265
#define PGDELTA 1

/* ---------------------------------------------------------------------- */

PairAIREBO::PairAIREBO(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  one_coeff = 1;

  maxlocal = 0;
  REBO_numneigh = NULL;
  REBO_firstneigh = NULL;
  closestdistsq = NULL;
  maxpage = 0;
  pages = NULL;
  nC = nH = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairAIREBO::~PairAIREBO()
{
  memory->sfree(REBO_numneigh);
  memory->sfree(REBO_firstneigh);
  for (int i = 0; i < maxpage; i++) memory->sfree(pages[i]);
  memory->sfree(pages);
  memory->sfree(closestdistsq);
  memory->sfree(nC);
  memory->sfree(nH);

  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);

    memory->destroy_2d_double_array(cutljsq);
    memory->destroy_2d_double_array(lj1);
    memory->destroy_2d_double_array(lj2);
    memory->destroy_2d_double_array(lj3);
    memory->destroy_2d_double_array(lj4);
    delete [] map;
  }
}

/* ---------------------------------------------------------------------- */

void PairAIREBO::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  REBO_neigh();
  FREBO(eflag,vflag);
  if (ljflag) FLJ(eflag,vflag);
  if (torflag) TORSION(eflag,vflag);
  
  if (vflag_fdotr) virial_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairAIREBO::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  // only sized by C,H = 2 types

  cutljsq = memory->create_2d_double_array(2,2,"pair:cutljsq");
  lj1 = memory->create_2d_double_array(2,2,"pair:lj1");
  lj2 = memory->create_2d_double_array(2,2,"pair:lj2");
  lj3 = memory->create_2d_double_array(2,2,"pair:lj3");
  lj4 = memory->create_2d_double_array(2,2,"pair:lj4");

  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairAIREBO::settings(int narg, char **arg)
{
  if (narg != 1 && narg != 3) error->all("Illegal pair_style command");

  cutlj = force->numeric(arg[0]);

  ljflag = torflag = 1;
  if (narg == 3) {
    ljflag = force->inumeric(arg[1]);
    torflag = force->inumeric(arg[2]);
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairAIREBO::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all("Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all("Incorrect args for pair coefficients");

  // read args that map atom types to C and H
  // map[i] = which element (0,1) the Ith atom type is, -1 if NULL

  for (int i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    } else if (strcmp(arg[i],"C") == 0) {
      map[i-2] = 0;
    } else if (strcmp(arg[i],"H") == 0) {
      map[i-2] = 1;
    } else error->all("Incorrect args for pair coefficients");
  }

  // read potential file and initialize fitting splines

  read_file(arg[2]);
  spline_init();

  // clear setflag since coeff() called once with I,J = * *

  int n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
	setflag[i][j] = 1;
	count++;
      }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairAIREBO::init_style()
{
  if (atom->tag_enable == 0)
    error->all("Pair style AIREBO requires atom IDs");
  if (force->newton_pair == 0)
    error->all("Pair style AIREBO requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  // local REBO neighbor list memory

  pgsize = neighbor->pgsize;
  oneatom = neighbor->oneatom;
  if (maxlocal == 0) add_pages(0);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairAIREBO::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all("All pair coeffs are not set");

  // convert to C,H types

  int ii = map[i];
  int jj = map[j];

  // use C-C values for these cutoffs since C atoms are biggest

  // cut3rebo = 3 REBO distances

  cut3rebo = 3.0 * rcmax[0][0];

  // cutljrebosq = furthest distance from an owned atom a ghost atom can be
  //               to need its REBO neighs computed
  // interaction = M-K-I-J-L-N with I = owned and J = ghost
  //   this insures N is in the REBO neigh list of L
  //   since I-J < rcLJmax and J-L < rmax

  double cutljrebo = rcLJmax[0][0] + rcmax[0][0];
  cutljrebosq = cutljrebo * cutljrebo;

  // cutmax = furthest distance from an owned atom
  //          at which another atom will feel force, i.e. the ghost cutoff
  // for REBO term in potential:
  //   interaction = M-K-I-J-L-N with I = owned and J = ghost
  //   I to N is max distance = 3 REBO distances
  // for LJ term in potential:
  //   short interaction = M-K-I-J-L-N with I = owned, J = ghost, I-J < rcLJmax
  //   rcLJmax + 2*rcmax, since I-J < rcLJmax and J-L,L-N = REBO distances
  //   long interaction = I-J with I = owned and J = ghost
  //   cutlj*sigma, since I-J < LJ cutoff

  double cutmax = cut3rebo;
  if (ljflag) {
    cutmax = MAX(cutmax,rcLJmax[0][0] + 2.0*rcmax[0][0]);
    cutmax = MAX(cutmax,cutlj*sigma[0][0]);
  }

  cutsq[i][j] = cutmax * cutmax;
  cutljsq[ii][jj] = cutlj*sigma[ii][jj] * cutlj*sigma[ii][jj];
  lj1[ii][jj] = 48.0 * epsilon[ii][jj] * pow(sigma[ii][jj],12.0);
  lj2[ii][jj] = 24.0 * epsilon[ii][jj] * pow(sigma[ii][jj],6.0);
  lj3[ii][jj] = 4.0 * epsilon[ii][jj] * pow(sigma[ii][jj],12.0);
  lj4[ii][jj] = 4.0 * epsilon[ii][jj] * pow(sigma[ii][jj],6.0);

  cutsq[j][i] = cutsq[i][j];
  cutljsq[jj][ii] = cutljsq[ii][jj];
  lj1[jj][ii] = lj1[ii][jj];
  lj2[jj][ii] = lj2[ii][jj];
  lj3[jj][ii] = lj3[ii][jj];
  lj4[jj][ii] = lj4[ii][jj];

  return cutmax;
}

/* ----------------------------------------------------------------------
   create REBO neighbor list from main neighbor list
   REBO neighbor list stores neighbors of ghost atoms
------------------------------------------------------------------------- */

void PairAIREBO::REBO_neigh()
{
  int i,j,ii,jj,m,n,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,dS;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *neighptr;

  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (nall > maxlocal) {
    maxlocal = atom->nmax;
    memory->sfree(REBO_numneigh);
    memory->sfree(REBO_firstneigh);
    memory->sfree(closestdistsq);
    memory->sfree(nC);
    memory->sfree(nH);
    REBO_numneigh = (int *)
      memory->smalloc(maxlocal*sizeof(int),"AIREBO:numneigh");
    REBO_firstneigh = (int **)
      memory->smalloc(maxlocal*sizeof(int *),"AIREBO:firstneigh");
    closestdistsq = (double *) 
      memory->smalloc(maxlocal*sizeof(double),"AIREBO:closestdistsq");
    nC = (double *) memory->smalloc(maxlocal*sizeof(double),"AIREBO:nC");
    nH = (double *) memory->smalloc(maxlocal*sizeof(double),"AIREBO:nH");
  }
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // initialize ghost atom references to -1 and closest distance = cutljrebosq

  for (i = nlocal; i < nall; i++) {
    REBO_numneigh[i] = -1;
    closestdistsq[i] = cutljrebosq;
  }

  // store all REBO neighs of owned atoms
  // scan full neighbor list of I
  // if J is ghost and within LJ cutoff:
  //   flag it via REBO_numneigh so its REBO neighbors will be stored below
  //   REBO requires neighbors of neighbors of i,j in each i,j LJ interaction

  npage = 0;
  int npnt = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == maxpage) add_pages(npage);
    }
    neighptr = &pages[npage][npnt];
    n = 0;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = map[type[i]];
    nC[i] = nH[i] = 0.0;
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jtype = map[type[j]];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < rcmaxsq[itype][jtype]) {
	neighptr[n++] = j;
	if (jtype == 0)
	  nC[i] += Sp(sqrt(rsq),rcmin[itype][jtype],rcmax[itype][jtype],dS);
	else
	  nH[i] += Sp(sqrt(rsq),rcmin[itype][jtype],rcmax[itype][jtype],dS);
      }
      if (j >= nlocal && rsq < closestdistsq[j]) {
	REBO_numneigh[j] = i;
	closestdistsq[j] = rsq;
      }
    }

    REBO_firstneigh[i] = neighptr;
    REBO_numneigh[i] = n;
    npnt += n;
    if (npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }

  // store REBO neighs of ghost atoms I that have been flagged in REBO_numneigh
  // find by scanning full neighbor list of owned atom M = closest neigh of I

  for (i = nlocal; i < nall; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == maxpage) add_pages(npage);
    }
    neighptr = &pages[npage][npnt];
    n = 0;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = map[type[i]];
    nC[i] = nH[i] = 0.0;
    m = REBO_numneigh[i];
    if (m < 0) {
      REBO_firstneigh[i] = neighptr;
      REBO_numneigh[i] = 0;
      continue;
    }

    jtype = map[type[m]];
    delx = xtmp - x[m][0];
    dely = ytmp - x[m][1];
    delz = ztmp - x[m][2];
    rsq = delx*delx + dely*dely + delz*delz;

    // add M as neigh of I if close enough

    if (rsq < rcmaxsq[itype][jtype]) {
      neighptr[n++] = m;
      if (jtype == 0)
	nC[i] += Sp(sqrt(rsq),rcmin[itype][jtype],rcmax[itype][jtype],dS);
      else
	nH[i] += Sp(sqrt(rsq),rcmin[itype][jtype],rcmax[itype][jtype],dS);
    }

    jlist = firstneigh[m];
    jnum = numneigh[m];

    // add M's neighbors as neighs of I if close enough
    // skip I itself

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      if (j == i) continue;
      jtype = map[type[j]];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < rcmaxsq[itype][jtype]) {
	neighptr[n++] = j;
	if (jtype == 0)
	  nC[i] += Sp(sqrt(rsq),rcmin[itype][jtype],rcmax[itype][jtype],dS);
	else
	  nH[i] += Sp(sqrt(rsq),rcmin[itype][jtype],rcmax[itype][jtype],dS);
      }
    }

    REBO_firstneigh[i] = neighptr;
    REBO_numneigh[i] = n;
    npnt += n;
    if (npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }
}

/* ----------------------------------------------------------------------
   REBO forces and energy
------------------------------------------------------------------------- */

void PairAIREBO::FREBO(int eflag, int vflag)
{
  int i,j,k,m,ii,inum,itype,jtype;
  double delx,dely,delz,evdwl,fpair;
  double rsq,rij,wij;
  double Qij,Aij,alphaij,VR,pre,dVRdi,VA,term,bij,dVAdi,dVA;
  double dwij,del[3];
  int *ilist,*REBO_neighs;

  evdwl = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;

  // two-body interactions from REBO neighbor list, skip half of them

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = map[type[i]];
    REBO_neighs = REBO_firstneigh[i];

    for (k = 0; k < REBO_numneigh[i]; k++) {
      j = REBO_neighs[k];
      if (tag[i] > tag[j]) continue;
      jtype = map[type[j]];

      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      rij = sqrt(rsq);
      wij = Sp(rij,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
      if (wij <= TOL) continue;
      
      Qij = Q[itype][jtype];
      Aij = A[itype][jtype];
      alphaij = alpha[itype][jtype];

      VR = wij*(1.0+(Qij/rij)) * Aij*exp(-alphaij*rij);
      pre = wij*Aij * exp(-alphaij*rij);
      dVRdi = pre * ((-alphaij)-(Qij/rsq)-(Qij*alphaij/rij));
      dVRdi += VR/wij * dwij;

      VA = dVA = 0.0;
      for (m = 0; m < 3; m++) {
	term = -wij * BIJc[itype][jtype][m] * exp(-Beta[itype][jtype][m]*rij);
	VA += term;
	dVA += -Beta[itype][jtype][m] * term;
      }
      dVA += VA/wij * dwij;
      del[0] = delx;
      del[1] = dely;
      del[2] = delz;
      bij = bondorder(i,j,del,rij,VA,f,vflag_atom);
      dVAdi = bij*dVA;

      fpair = -(dVRdi+dVAdi) / rij;
      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (eflag) evdwl = VR + bij*VA;
      if (evflag) ev_tally(i,j,nlocal,newton_pair,
			   evdwl,0.0,fpair,delx,dely,delz);
    }
  }
}

/* ----------------------------------------------------------------------
   compute LJ forces and energy
   find 3- and 4-step paths between atoms I,J via REBO neighbor lists
------------------------------------------------------------------------- */

void PairAIREBO::FLJ(int eflag, int vflag)
{
  int i,j,k,m,ii,jj,kk,mm,inum,jnum,itype,jtype,ktype,mtype;
  int atomi,atomj,atomk,atomm;
  int testpath,npath,done;
  double evdwl,fpair;
  double rsq,best,wik,wkm,cij,rij,dwij,dwik,dwkj,dwkm,dwmj;
  double delij[3],rijsq,delik[3],rik,deljk[3];
  double rkj,wkj,dC,VLJ,dVLJ,VA,Str,dStr,Stb;
  double delkm[3],rkm,deljm[3],rmj,wmj,r2inv,r6inv,scale,delscale[3];
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *REBO_neighs_i,*REBO_neighs_k;
  double delikS[3],deljkS[3],delkmS[3],deljmS[3],delimS[3];
  double rikS,rkjS,rkmS,rmjS,wikS,dwikS;
  double wkjS,dwkjS,wkmS,dwkmS,wmjS,dwmjS;
  double fpair1,fpair2,fpair3;
  double fi[3],fj[3],fk[3],fm[3];

  // I-J interaction from full neighbor list
  // skip 1/2 of interactions since only consider each pair once

  evdwl = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  int *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = map[type[i]];
    atomi = i;
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      atomj = j;
      if (tag[i] > tag[j]) continue;
      jtype = map[type[j]];

      delij[0] = x[i][0] - x[j][0];
      delij[1] = x[i][1] - x[j][1];
      delij[2] = x[i][2] - x[j][2];
      rijsq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];

      // if outside of LJ cutoff, skip
      // if outside of 4-path cutoff, best = 0.0, no need to test paths
      // if outside of 2-path cutoff but inside 4-path cutoff,
      //   best = 0.0, test 3-,4-paths
      // if inside 2-path cutoff, best = wij, only test 3-,4-paths if best < 1
	
      if (rijsq >= cutljsq[itype][jtype]) continue;
      rij = sqrt(rijsq);
      if (rij >= cut3rebo) {
	best = 0.0;
	testpath = 0;
      } else if (rij >= rcmax[itype][jtype]) {
	best = 0.0;
	testpath = 1;
      } else {
	best = Sp(rij,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
	npath = 2;
	if (best < 1.0) testpath = 1;
	else testpath = 0;
      }

      done = 0;
      if (testpath) {

	// test all 3-body paths = I-K-J
	// I-K interactions come from atom I's REBO neighbors
	// if wik > current best, compute wkj
	// if best = 1.0, done

	REBO_neighs_i = REBO_firstneigh[i];
	for (kk = 0; kk < REBO_numneigh[i] && done==0; kk++) {
	  k = REBO_neighs_i[kk];
	  if (k == j) continue;
	  ktype = map[type[k]];

	  delik[0] = x[i][0] - x[k][0];
	  delik[1] = x[i][1] - x[k][1];
	  delik[2] = x[i][2] - x[k][2];
	  rsq = delik[0]*delik[0] + delik[1]*delik[1] + delik[2]*delik[2];
	  if (rsq < rcmaxsq[itype][ktype]) {
	    rik = sqrt(rsq);
	    wik = Sp(rik,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
	  } else wik = 0.0;

	  if (wik > best) {
	    deljk[0] = x[j][0] - x[k][0];
	    deljk[1] = x[j][1] - x[k][1];
	    deljk[2] = x[j][2] - x[k][2];
	    rsq = deljk[0]*deljk[0] + deljk[1]*deljk[1] + deljk[2]*deljk[2];
	    if (rsq < rcmaxsq[ktype][jtype]) {
	      rkj = sqrt(rsq);
	      wkj = Sp(rkj,rcmin[ktype][jtype],rcmax[ktype][jtype],dwkj);
	      if (wik*wkj > best) {
		best = wik*wkj;
		npath = 3;
 		atomk = k;
            	delikS[0] = delik[0];
            	delikS[1] = delik[1];
            	delikS[2] = delik[2];
	    	rikS = rik;
	    	wikS = wik;
	    	dwikS = dwik;
            	deljkS[0] = deljk[0];
            	deljkS[1] = deljk[1];
            	deljkS[2] = deljk[2];
	    	rkjS = rkj;
	    	wkjS = wkj;
	    	dwkjS = dwkj;
		if (best == 1.0) {
		  done = 1; 
		  break;
		}
	      }
	    }

	    // test all 4-body paths = I-K-M-J
	    // K-M interactions come from atom K's REBO neighbors
	    // if wik*wkm > current best, compute wmj
	    // if best = 1.0, done

	    REBO_neighs_k = REBO_firstneigh[k];
	    for (mm = 0; mm < REBO_numneigh[k] && done==0; mm++) {
	      m = REBO_neighs_k[mm];
	      if (m == i || m == j) continue;
	      mtype = map[type[m]];
	      delkm[0] = x[k][0] - x[m][0];
	      delkm[1] = x[k][1] - x[m][1];
	      delkm[2] = x[k][2] - x[m][2];
	      rsq = delkm[0]*delkm[0] + delkm[1]*delkm[1] + delkm[2]*delkm[2];
	      if (rsq < rcmaxsq[ktype][mtype]) {
		rkm = sqrt(rsq);
		wkm = Sp(rkm,rcmin[ktype][mtype],rcmax[ktype][mtype],dwkm);
	      } else wkm = 0.0;

	      if (wik*wkm > best) {
		deljm[0] = x[j][0] - x[m][0];
		deljm[1] = x[j][1] - x[m][1];
		deljm[2] = x[j][2] - x[m][2];
		rsq = deljm[0]*deljm[0] + deljm[1]*deljm[1] +
		  deljm[2]*deljm[2];
		if (rsq < rcmaxsq[mtype][jtype]) {
		  rmj = sqrt(rsq);
		  wmj = Sp(rmj,rcmin[mtype][jtype],rcmax[mtype][jtype],dwmj);
		  if (wik*wkm*wmj > best) {
		    best = wik*wkm*wmj;
		    npath = 4;
		    atomk = k;
		    delikS[0] = delik[0];
	            delikS[1] = delik[1];
        	    delikS[2] = delik[2];
		    rikS = rik;
		    wikS = wik;
		    dwikS = dwik;
		    atomm = m;
	            delkmS[0] = delkm[0];
        	    delkmS[1] = delkm[1];
            	    delkmS[2] = delkm[2];
		    rkmS = rkm;
		    wkmS = wkm;
		    dwkmS = dwkm;
	            deljmS[0] = deljm[0];
        	    deljmS[1] = deljm[1];
           	    deljmS[2] = deljm[2];
		    rmjS = rmj;
		    wmjS = wmj;
		    dwmjS = dwmj;
		    if (best == 1.0) {
		      done = 1;
		      break;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
		
      cij = 1.0 - best;
      if (cij == 0.0) continue;

      // compute LJ forces and energy
  
      r2inv = 1.0/rijsq;
      r6inv = r2inv*r2inv*r2inv;
      VLJ = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);
      dVLJ = -r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]) / rij;
      
      Str = Sp2(rij,rcLJmin[itype][jtype],rcLJmax[itype][jtype],dStr);
      VA = Str*cij*VLJ;
      if (Str > 0.0) {
	scale = rcmin[itype][jtype] / rij;
	delscale[0] = scale * delij[0];
	delscale[1] = scale * delij[1];
	delscale[2] = scale * delij[2];
	Stb = bondorderLJ(i,j,delscale,rcmin[itype][jtype],VA,
			  delij,rij,f,vflag_atom);
      } else Stb = 0.0;

      fpair = -(dStr * (Stb*cij*VLJ - cij*VLJ) +
		dVLJ * (Str*Stb*cij + cij - Str*cij)) / rij;
      
      f[i][0] += delij[0]*fpair;
      f[i][1] += delij[1]*fpair;
      f[i][2] += delij[2]*fpair;
      f[j][0] -= delij[0]*fpair;
      f[j][1] -= delij[1]*fpair;
      f[j][2] -= delij[2]*fpair;
      
      if (eflag) evdwl = VA*Stb + (1.0-Str)*cij*VLJ;
      if (evflag) ev_tally(i,j,nlocal,newton_pair,
			   evdwl,0.0,fpair,delij[0],delij[1],delij[2]);

      if (cij < 1.0) {
	dC = Str*Stb*VLJ + (1.0-Str)*VLJ;
	if (npath == 2) {
	  fpair = dC*dwij / rij;
	  f[atomi][0] += delij[0]*fpair;
	  f[atomi][1] += delij[1]*fpair;
	  f[atomi][2] += delij[2]*fpair;
	  f[atomj][0] -= delij[0]*fpair;
	  f[atomj][1] -= delij[1]*fpair;
	  f[atomj][2] -= delij[2]*fpair;

	  if (vflag_atom) v_tally2(atomi,atomj,fpair,delij);

	} else if (npath == 3) {
	  fpair1 = dC*dwikS*wkjS / rikS;
	  fi[0] = delikS[0]*fpair1;
	  fi[1] = delikS[1]*fpair1;
	  fi[2] = delikS[2]*fpair1;
	  fpair2 = dC*wikS*dwkjS / rkjS;
	  fj[0] = deljkS[0]*fpair2;
	  fj[1] = deljkS[1]*fpair2;
	  fj[2] = deljkS[2]*fpair2;

	  f[atomi][0] += fi[0];
	  f[atomi][1] += fi[1];
	  f[atomi][2] += fi[2];
	  f[atomj][0] += fj[0];
	  f[atomj][1] += fj[1];
	  f[atomj][2] += fj[2];
	  f[atomk][0] -= fi[0] + fj[0];
	  f[atomk][1] -= fi[1] + fj[1];
	  f[atomk][2] -= fi[2] + fj[2];

	  if (vflag_atom) 
	    v_tally3(atomi,atomj,atomk,fi,fj,delikS,deljkS);

	} else {
	  fpair1 = dC*dwikS*wkmS*wmjS / rikS;
	  fi[0] = delikS[0]*fpair1;
	  fi[1] = delikS[1]*fpair1;
	  fi[2] = delikS[2]*fpair1;

	  fpair2 = dC*wikS*dwkmS*wmjS / rkmS;
	  fk[0] = delkmS[0]*fpair2 - fi[0];
	  fk[1] = delkmS[1]*fpair2 - fi[1];
	  fk[2] = delkmS[2]*fpair2 - fi[2];

	  fpair3 = dC*wikS*wkmS*dwmjS / rmjS;
	  fj[0] = deljmS[0]*fpair3;
	  fj[1] = deljmS[1]*fpair3;
	  fj[2] = deljmS[2]*fpair3;

	  fm[0] = -delkmS[0]*fpair2 - fj[0];
	  fm[1] = -delkmS[1]*fpair2 - fj[1];
	  fm[2] = -delkmS[2]*fpair2 - fj[2];

	  f[atomi][0] += fi[0];
	  f[atomi][1] += fi[1];
	  f[atomi][2] += fi[2];
	  f[atomj][0] += fj[0];
	  f[atomj][1] += fj[1];
	  f[atomj][2] += fj[2];
	  f[atomk][0] += fk[0];
	  f[atomk][1] += fk[1];
	  f[atomk][2] += fk[2];
	  f[atomm][0] += fm[0];
	  f[atomm][1] += fm[1];
	  f[atomm][2] += fm[2];

	  if (vflag_atom) {
	    delimS[0] = delikS[0] + delkmS[0];
	    delimS[1] = delikS[1] + delkmS[1];
	    delimS[2] = delikS[2] + delkmS[2];
	    v_tally4(atomi,atomj,atomk,atomm,fi,fj,fk,delimS,deljmS,delkmS);
	  }
	}
      } 
    }
  }
}

/* ----------------------------------------------------------------------
   Torsional forces and energy
------------------------------------------------------------------------- */

void PairAIREBO::TORSION(int eflag, int vflag)
{
  int i,j,k,l,ii,inum;
  double evdwl,fpair;
  double cos321;
  double w21,dw21,cos234,w34,dw34;
  double cross321[3],cross321mag,cross234[3],cross234mag;
  double w23,dw23,cw2,ekijl,Ec;
  double cw,cwnum,cwnom;
  double rij,rij2,rik,rjl,tspjik,dtsjik,tspijl,dtsijl,costmp,fcpc;
  double sin321,sin234,rjk2,rik2,ril2,rjl2;
  double rjk,ril;
  double Vtors;
  double dndij[3],tmpvec[3],dndik[3],dndjl[3];
  double dcidij,dcidik,dcidjk,dcjdji,dcjdjl,dcjdil;
  double dsidij,dsidik,dsidjk,dsjdji,dsjdjl,dsjdil;
  double dxidij,dxidik,dxidjk,dxjdji,dxjdjl,dxjdil;
  double ddndij,ddndik,ddndjk,ddndjl,ddndil,dcwddn,dcwdn,dvpdcw,Ftmp[3];
  double del32[3],rsq,r32,del23[3],del21[3],r21;
  double deljk[3],del34[3],delil[3],delkl[3],r23,r34;
  double fi[3],fj[3],fk[3],fl[3];
  int itype,jtype,ktype,ltype,kk,ll,jj;
  int *ilist,*REBO_neighs_i,*REBO_neighs_j;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int *tag = atom->tag;

  inum = list->inum;
  ilist = list->ilist;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = map[type[i]];
    if (itype != 0) continue;
    REBO_neighs_i = REBO_firstneigh[i];

    for (jj = 0; jj < REBO_numneigh[i]; jj++) {
      j = REBO_neighs_i[jj];
      if (tag[i] > tag[j]) continue;
      jtype = map[type[j]];
      if (jtype != 0) continue;
      del32[0] = x[j][0]-x[i][0];
      del32[1] = x[j][1]-x[i][1];
      del32[2] = x[j][2]-x[i][2];
      rsq = del32[0]*del32[0] + del32[1]*del32[1] + del32[2]*del32[2];
      r32 = sqrt(rsq);
      del23[0] = -del32[0];
      del23[1] = -del32[1];
      del23[2] = -del32[2];
      r23 = r32;
      w23 = Sp(r23,rcmin[itype][jtype],rcmax[itype][jtype],dw23);

      for (kk = 0; kk < REBO_numneigh[i]; kk++) {
	k = REBO_neighs_i[kk];
	ktype = map[type[k]];
	if (k == j) continue;
	del21[0] = x[i][0]-x[k][0];
	del21[1] = x[i][1]-x[k][1];
	del21[2] = x[i][2]-x[k][2];
	rsq = del21[0]*del21[0] + del21[1]*del21[1] + del21[2]*del21[2];
	r21 = sqrt(rsq);
	cos321 = - ((del21[0]*del32[0]) + (del21[1]*del32[1]) +
		    (del21[2]*del32[2])) / (r21*r32);
	cos321 = MIN(cos321,1.0);
	cos321 = MAX(cos321,-1.0);
	sin321 = sqrt(1.0 - cos321*cos321);
	if (sin321 < TOL) continue;
	
	deljk[0] = del21[0]-del23[0];
	deljk[1] = del21[1]-del23[1];
	deljk[2] = del21[2]-del23[2];
	rjk2 = deljk[0]*deljk[0] + deljk[1]*deljk[1] + deljk[2]*deljk[2];
	rjk=sqrt(rjk2);
	rik2 = r21*r21;
	w21 = Sp(r21,rcmin[itype][ktype],rcmax[itype][ktype],dw21);

	rij = r32;
	rik = r21;
	rij2 = r32*r32;
	rik2 = r21*r21;
	costmp = 0.5*(rij2+rik2-rjk2)/rij/rik;
	tspjik = Sp2(costmp,thmin,thmax,dtsjik);
	dtsjik = -dtsjik;

	REBO_neighs_j = REBO_firstneigh[j];
	for (ll = 0; ll < REBO_numneigh[j]; ll++) {
	  l = REBO_neighs_j[ll];
	  ltype = map[type[l]];
	  if (l == i || l == k) continue;
	  del34[0] = x[j][0]-x[l][0];
	  del34[1] = x[j][1]-x[l][1];
	  del34[2] = x[j][2]-x[l][2];
	  rsq = del34[0]*del34[0] + del34[1]*del34[1] + del34[2]*del34[2];
	  r34 = sqrt(rsq);
	  cos234 = (del32[0]*del34[0] + del32[1]*del34[1] +
		    del32[2]*del34[2]) / (r32*r34);
	  cos234 = MIN(cos234,1.0);
	  cos234 = MAX(cos234,-1.0);
	  sin234 = sqrt(1.0 - cos234*cos234);
	  if (sin234 < TOL) continue;
	  w34 = Sp(r34,rcmin[jtype][ltype],rcmax[jtype][ltype],dw34);
	  delil[0] = del23[0] + del34[0];
	  delil[1] = del23[1] + del34[1];
	  delil[2] = del23[2] + del34[2];
	  ril2 = delil[0]*delil[0] + delil[1]*delil[1] + delil[2]*delil[2];
	  ril=sqrt(ril2);
	  rjl2 = r34*r34;

	  rjl = r34;
	  rjl2 = r34*r34;
	  costmp = 0.5*(rij2+rjl2-ril2)/rij/rjl;
	  tspijl = Sp2(costmp,thmin,thmax,dtsijl);
	  dtsijl = -dtsijl; //need minus sign 
	  cross321[0] = (del32[1]*del21[2])-(del32[2]*del21[1]);
	  cross321[1] = (del32[2]*del21[0])-(del32[0]*del21[2]);
	  cross321[2] = (del32[0]*del21[1])-(del32[1]*del21[0]);
	  cross321mag = sqrt(cross321[0]*cross321[0]+
			     cross321[1]*cross321[1]+
			     cross321[2]*cross321[2]);
	  cross234[0] = (del23[1]*del34[2])-(del23[2]*del34[1]);
	  cross234[1] = (del23[2]*del34[0])-(del23[0]*del34[2]);
	  cross234[2] = (del23[0]*del34[1])-(del23[1]*del34[0]);
	  cross234mag = sqrt(cross234[0]*cross234[0]+
			     cross234[1]*cross234[1]+
			     cross234[2]*cross234[2]);
	  cwnum = (cross321[0]*cross234[0]) +
	    (cross321[1]*cross234[1])+(cross321[2]*cross234[2]);
	  cwnom = r21*r34*r32*r32*sin321*sin234;
	  cw = cwnum/cwnom;
	  
	  cw2 = (.5*(1.0-cw));
	  ekijl = epsilonT[ktype][ltype];
	  Ec = 256.0*ekijl/405.0;
	  Vtors = (Ec*(pow(cw2,5.0)))-(ekijl/10.0);

	  if (eflag) evdwl = Vtors*w21*w23*w34*(1.0-tspjik)*(1.0-tspijl);
	  
	  dndij[0] = (cross234[1]*del21[2])-(cross234[2]*del21[1]);
	  dndij[1] = (cross234[2]*del21[0])-(cross234[0]*del21[2]);
	  dndij[2] = (cross234[0]*del21[1])-(cross234[1]*del21[0]); 
	  
	  tmpvec[0] = (del34[1]*cross321[2])-(del34[2]*cross321[1]);
	  tmpvec[1] = (del34[2]*cross321[0])-(del34[0]*cross321[2]);
	  tmpvec[2] = (del34[0]*cross321[1])-(del34[1]*cross321[0]); 
	  
	  dndij[0] = dndij[0]+tmpvec[0];
	  dndij[1] = dndij[1]+tmpvec[1];
	  dndij[2] = dndij[2]+tmpvec[2];
	  
	  dndik[0] = (del23[1]*cross234[2])-(del23[2]*cross234[1]);
	  dndik[1] = (del23[2]*cross234[0])-(del23[0]*cross234[2]);
	  dndik[2] = (del23[0]*cross234[1])-(del23[1]*cross234[0]);
	  
	  dndjl[0] = (cross321[1]*del23[2])-(cross321[2]*del23[1]);
	  dndjl[1] = (cross321[2]*del23[0])-(cross321[0]*del23[2]);
	  dndjl[2] = (cross321[0]*del23[1])-(cross321[1]*del23[0]); 
	  
	  dcidij = ((r23*r23)-(r21*r21)+(rjk*rjk))/(2.0*r23*r23*r21);
	  dcidik = ((r21*r21)-(r23*r23)+(rjk*rjk))/(2.0*r23*r21*r21);
	  dcidjk = (-rjk)/(r23*r21);
	  dcjdji = ((r23*r23)-(r34*r34)+(ril*ril))/(2.0*r23*r23*r34);
	  dcjdjl = ((r34*r34)-(r23*r23)+(ril*ril))/(2.0*r23*r34*r34);
	  dcjdil = (-ril)/(r23*r34);
	  
	  dsidij = (-cos321/sin321)*dcidij;
	  dsidik = (-cos321/sin321)*dcidik;
	  dsidjk = (-cos321/sin321)*dcidjk;
	  
	  dsjdji = (-cos234/sin234)*dcjdji;
	  dsjdjl = (-cos234/sin234)*dcjdjl;
	  dsjdil = (-cos234/sin234)*dcjdil;
	  
	  dxidij = (r21*sin321)+(r23*r21*dsidij);
	  dxidik = (r23*sin321)+(r23*r21*dsidik);
	  dxidjk = (r23*r21*dsidjk);
	  
	  dxjdji = (r34*sin234)+(r23*r34*dsjdji);
	  dxjdjl = (r23*sin234)+(r23*r34*dsjdjl);
	  dxjdil = (r23*r34*dsjdil);
	  
	  ddndij = (dxidij*cross234mag)+(cross321mag*dxjdji);
	  ddndik = dxidik*cross234mag;
	  ddndjk = dxidjk*cross234mag;
	  ddndjl = cross321mag*dxjdjl;
	  ddndil = cross321mag*dxjdil;
	  dcwddn = -cwnum/(cwnom*cwnom);
	  dcwdn = 1.0/cwnom;
	  dvpdcw = (-1.0)*Ec*(-.5)*5.0*pow(cw2,4.0) * 
	    w23*w21*w34*(1.0-tspjik)*(1.0-tspijl);

	  Ftmp[0] = dvpdcw*((dcwdn*dndij[0])+(dcwddn*ddndij*del23[0]/r23));
	  Ftmp[1] = dvpdcw*((dcwdn*dndij[1])+(dcwddn*ddndij*del23[1]/r23));
	  Ftmp[2] = dvpdcw*((dcwdn*dndij[2])+(dcwddn*ddndij*del23[2]/r23));
	  fi[0] = Ftmp[0];
	  fi[1] = Ftmp[1];
	  fi[2] = Ftmp[2];
	  fj[0] = -Ftmp[0];
	  fj[1] = -Ftmp[1];
	  fj[2] = -Ftmp[2];
	  
	  Ftmp[0] = dvpdcw*((dcwdn*dndik[0])+(dcwddn*ddndik*del21[0]/r21));
	  Ftmp[1] = dvpdcw*((dcwdn*dndik[1])+(dcwddn*ddndik*del21[1]/r21));
	  Ftmp[2] = dvpdcw*((dcwdn*dndik[2])+(dcwddn*ddndik*del21[2]/r21));
	  fi[0] += Ftmp[0];
	  fi[1] += Ftmp[1];
	  fi[2] += Ftmp[2];
	  fk[0] = -Ftmp[0];
	  fk[1] = -Ftmp[1];
	  fk[2] = -Ftmp[2];
	  
	  Ftmp[0] = (dvpdcw*dcwddn*ddndjk*deljk[0])/rjk;
	  Ftmp[1] = (dvpdcw*dcwddn*ddndjk*deljk[1])/rjk;
	  Ftmp[2] = (dvpdcw*dcwddn*ddndjk*deljk[2])/rjk;
	  fj[0] += Ftmp[0];
	  fj[1] += Ftmp[1];
	  fj[2] += Ftmp[2];
	  fk[0] -= Ftmp[0];
	  fk[1] -= Ftmp[1];
	  fk[2] -= Ftmp[2];
	  
	  Ftmp[0] = dvpdcw*((dcwdn*dndjl[0])+(dcwddn*ddndjl*del34[0]/r34));
	  Ftmp[1] = dvpdcw*((dcwdn*dndjl[1])+(dcwddn*ddndjl*del34[1]/r34));
	  Ftmp[2] = dvpdcw*((dcwdn*dndjl[2])+(dcwddn*ddndjl*del34[2]/r34));
	  fj[0] += Ftmp[0];
	  fj[1] += Ftmp[1];
	  fj[2] += Ftmp[2];
	  fl[0] = -Ftmp[0];
	  fl[1] = -Ftmp[1];
	  fl[2] = -Ftmp[2];
	  
	  Ftmp[0] = (dvpdcw*dcwddn*ddndil*delil[0])/ril;
	  Ftmp[1] = (dvpdcw*dcwddn*ddndil*delil[1])/ril;
	  Ftmp[2] = (dvpdcw*dcwddn*ddndil*delil[2])/ril;
	  fi[0] += Ftmp[0];
	  fi[1] += Ftmp[1];
	  fi[2] += Ftmp[2];
	  fl[0] -= Ftmp[0];
	  fl[1] -= Ftmp[1];
	  fl[2] -= Ftmp[2];
	  
	  // coordination forces
	  
	  fpair = Vtors*dw21*w23*w34*(1.0-tspjik)*(1.0-tspijl) / r21;
	  fi[0] -= del21[0]*fpair;
	  fi[1] -= del21[1]*fpair;
	  fi[2] -= del21[2]*fpair;
	  fk[0] += del21[0]*fpair;
	  fk[1] += del21[1]*fpair;
	  fk[2] += del21[2]*fpair;
	  
	  fpair = Vtors*w21*dw23*w34*(1.0-tspjik)*(1.0-tspijl) / r23;
	  fi[0] -= del23[0]*fpair;
	  fi[1] -= del23[1]*fpair;
	  fi[2] -= del23[2]*fpair;
	  fj[0] += del23[0]*fpair;
	  fj[1] += del23[1]*fpair;
	  fj[2] += del23[2]*fpair;
	  
	  fpair = Vtors*w21*w23*dw34*(1.0-tspjik)*(1.0-tspijl) / r34;
	  fj[0] -= del34[0]*fpair;
	  fj[1] -= del34[1]*fpair;
	  fj[2] -= del34[2]*fpair;
	  fl[0] += del34[0]*fpair;
	  fl[1] += del34[1]*fpair;
	  fl[2] += del34[2]*fpair;

	  // additional cut off function forces

	  fcpc = -Vtors*w21*w23*w34*dtsjik*(1.0-tspijl);
	  fpair = fcpc*dcidij/rij;
	  fi[0] += fpair*del23[0];
	  fi[1] += fpair*del23[1];
	  fi[2] += fpair*del23[2];
	  fj[0] -= fpair*del23[0];
	  fj[1] -= fpair*del23[1];
	  fj[2] -= fpair*del23[2];

	  fpair = fcpc*dcidik/rik;
	  fi[0] += fpair*del21[0];
	  fi[1] += fpair*del21[1];
	  fi[2] += fpair*del21[2];
	  fk[0] -= fpair*del21[0];
	  fk[1] -= fpair*del21[1];
	  fk[2] -= fpair*del21[2];

	  fpair = fcpc*dcidjk/rjk;
	  fj[0] += fpair*deljk[0];
	  fj[1] += fpair*deljk[1];
	  fj[2] += fpair*deljk[2];
	  fk[0] -= fpair*deljk[0];
	  fk[1] -= fpair*deljk[1];
	  fk[2] -= fpair*deljk[2];
	  
	  fcpc = -Vtors*w21*w23*w34*(1.0-tspjik)*dtsijl;
	  fpair = fcpc*dcjdji/rij;
	  fi[0] += fpair*del23[0];
	  fi[1] += fpair*del23[1];
	  fi[2] += fpair*del23[2];
	  fj[0] -= fpair*del23[0];
	  fj[1] -= fpair*del23[1];
	  fj[2] -= fpair*del23[2];

	  fpair = fcpc*dcjdjl/rjl;
	  fj[0] += fpair*del34[0];
	  fj[1] += fpair*del34[1];
	  fj[2] += fpair*del34[2];
	  fl[0] -= fpair*del34[0];
	  fl[1] -= fpair*del34[1];
	  fl[2] -= fpair*del34[2];

	  fpair = fcpc*dcjdil/ril;
	  fi[0] += fpair*delil[0];
	  fi[1] += fpair*delil[1];
	  fi[2] += fpair*delil[2];
	  fl[0] -= fpair*delil[0];
	  fl[1] -= fpair*delil[1];
	  fl[2] -= fpair*delil[2];

	  // sum per-atom forces into atom force array

	  f[i][0] += fi[0]; f[i][1] += fi[1]; f[i][2] += fi[2];
	  f[j][0] += fj[0]; f[j][1] += fj[1]; f[j][2] += fj[2];
	  f[k][0] += fk[0]; f[k][1] += fk[1]; f[k][2] += fk[2];
	  f[l][0] += fl[0]; f[l][1] += fl[1]; f[l][2] += fl[2];

	  if (evflag) {
	    delkl[0] = delil[0] - del21[0];
	    delkl[1] = delil[1] - del21[1];
	    delkl[2] = delil[2] - del21[2];
	    ev_tally4(i,j,k,l,evdwl,fi,fj,fk,delil,del34,delkl);
	  }
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// S'(t) and S(t) cutoff functions
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   cutoff function Sprime
   return cutoff and dX = derivative
------------------------------------------------------------------------- */

double PairAIREBO::Sp(double Xij, double Xmin, double Xmax, double &dX)
{
  double cutoff;

  double t = (Xij-Xmin) / (Xmax-Xmin);
  if (t <= 0.0) {
    cutoff = 1.0;
    dX = 0.0;
  } 
  else if (t >= 1.0) {
    cutoff = 0.0;
    dX = 0.0;
  } 
  else {
    cutoff = 0.5 * (1.0+cos(PI*t));
    dX = (-0.5*PI*sin(PI*t)) / (Xmax-Xmin);
  }
  return cutoff;
}

/* ----------------------------------------------------------------------
   LJ cutoff function Sp2
   return cutoff and dX = derivative
------------------------------------------------------------------------- */

double PairAIREBO::Sp2(double Xij, double Xmin, double Xmax, double &dX)
{
  double cutoff;

  double t = (Xij-Xmin) / (Xmax-Xmin);
  if (t <= 0.0) {
    cutoff = 1.0;
    dX = 0.0;
  } 
  if (t >= 1.0) {
    cutoff = 0.0;
    dX = 0.0;
  } 
  if (t>0.0 && t<1.0) {
    cutoff = (1.0-(t*t*(3.0-2.0*t)));
    dX = 6.0*(t*t-t) / (Xmax-Xmin);
  }
  return cutoff;
}

/* ----------------------------------------------------------------------
   Bij function
------------------------------------------------------------------------- */

double PairAIREBO::bondorder(int i, int j, double rij[3],
			     double rijmag, double VA,
			     double **f, int vflag_atom)
{
  int atomi,atomj,k,n,l,atomk,atoml,atomn,atom1,atom2,atom3,atom4;
  int itype,jtype,ktype,ltype,ntype;
  double rik[3],rjl[3],rkn[3],rji[3],rki[3],rlj[3],rknmag,dNki,dwjl,bij;
  double NijC,NijH,NjiC,NjiH,wik,dwik,dwkn,wjl;
  double rikmag,rjlmag,cosjik,cosijl,g,tmp2,tmp3;
  double Etmp,pij,tmp,wij,dwij,NconjtmpI,NconjtmpJ,Nki,Nlj,dS;
  double lamdajik,lamdaijl,dgdc,dgdN,pji,Nijconj,piRC;
  double dcosjikdri[3],dcosijldri[3],dcosjikdrk[3];
  double dN2[2],dN3[3];
  double dcosjikdrj[3],dcosijldrj[3],dcosijldrl[3];
  double Tij;
  double r32[3],r32mag,cos321,r43[3],r13[3];
  double dNlj;
  double om1234,rln[3];
  double rlnmag,dwln,r23[3],r23mag,r21[3],r21mag;
  double w21,dw21,r34[3],r34mag,cos234,w34,dw34;
  double cross321[3],cross234[3],prefactor,SpN;
  double fcijpc,fcikpc,fcjlpc,fcjkpc,fcilpc;
  double dt2dik[3],dt2djl[3],dt2dij[3],aa,aaa1,aaa2,at2,cw,cwnum,cwnom;
  double sin321,sin234,rr,rijrik,rijrjl,rjk2,rik2,ril2,rjl2;
  double dctik,dctjk,dctjl,dctij,dctji,dctil,rik2i,rjl2i,sink2i,sinl2i;
  double rjk[3],ril[3],dt1dik,dt1djk,dt1djl,dt1dil,dt1dij;
  double F23[3],F12[3],F34[3],F31[3],F24[3],fi[3],fj[3],fk[3],fl[3];
  double f1[3],f2[3],f3[3],f4[4];
  double dcut321,PijS,PjiS;
  double rij2,tspjik,dtsjik,tspijl,dtsijl,costmp;
  int *REBO_neighs,*REBO_neighs_i,*REBO_neighs_j,*REBO_neighs_k,*REBO_neighs_l;

  double **x = atom->x;
  int *type = atom->type;
  
  atomi = i;
  atomj = j;
  itype = map[type[i]];
  jtype = map[type[j]];
  wij = Sp(rijmag,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
  NijC = nC[i]-(wij*kronecker(jtype,0));
  NijH = nH[i]-(wij*kronecker(jtype,1));
  NjiC = nC[j]-(wij*kronecker(itype,0));
  NjiH = nH[j]-(wij*kronecker(itype,1));
  bij = 0.0;
  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  dgdc = 0.0;
  dgdN = 0.0;
  NconjtmpI = 0.0;
  NconjtmpJ = 0.0;
  Etmp = 0.0;
  
  REBO_neighs = REBO_firstneigh[i];
  for (k = 0; k < REBO_numneigh[i]; k++) {
    atomk = REBO_neighs[k];
    if (atomk != atomj) {
      ktype = map[type[atomk]];
      rik[0] = x[atomi][0]-x[atomk][0];
      rik[1] = x[atomi][1]-x[atomk][1];
      rik[2] = x[atomi][2]-x[atomk][2];
      rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      lamdajik = 4.0*kronecker(itype,1) * 
	((rho[ktype][1]-rikmag)-(rho[jtype][1]-rijmag));
      wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dS);
      Nki = nC[atomk]-(wik*kronecker(itype,0))+nH[atomk] -
	(wik*kronecker(itype,1));
      cosjik = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2])) / 
	(rijmag*rikmag);
      cosjik = MIN(cosjik,1.0);
      cosjik = MAX(cosjik,-1.0);

      // evaluate splines g and derivatives dg

      g = gSpline(cosjik,(NijC+NijH),itype,&dgdc,&dgdN);
      Etmp = Etmp+(wik*g*exp(lamdajik)); 
      tmp3 = tmp3+(wik*dgdN*exp(lamdajik));
      NconjtmpI = NconjtmpI+(kronecker(ktype,0)*wik*Sp(Nki,Nmin,Nmax,dS));
    }
  }

  PijS = 0.0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  PijS = PijSpline(NijC,NijH,itype,jtype,dN2);
  pij = pow(1.0+Etmp+PijS,-0.5); 
  tmp = -0.5*pow(pij,3.0);

  // pij forces

  REBO_neighs = REBO_firstneigh[i];
  for (k = 0; k < REBO_numneigh[i]; k++) {
    atomk = REBO_neighs[k];
    if (atomk != atomj) {
      ktype = map[type[atomk]];
      rik[0] = x[atomi][0]-x[atomk][0];
      rik[1] = x[atomi][1]-x[atomk][1];
      rik[2] = x[atomi][2]-x[atomk][2];
      rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      lamdajik = 4.0*kronecker(itype,1) * 
	((rho[ktype][1]-rikmag)-(rho[jtype][1]-rijmag));
      wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
      cosjik = (rij[0]*rik[0] + rij[1]*rik[1] + rij[2]*rik[2]) / 
	(rijmag*rikmag);
      cosjik = MIN(cosjik,1.0);
      cosjik = MAX(cosjik,-1.0);

      dcosjikdri[0] = ((rij[0]+rik[0])/(rijmag*rikmag)) - 
	(cosjik*((rij[0]/(rijmag*rijmag))+(rik[0]/(rikmag*rikmag))));
      dcosjikdri[1] = ((rij[1]+rik[1])/(rijmag*rikmag)) - 
	(cosjik*((rij[1]/(rijmag*rijmag))+(rik[1]/(rikmag*rikmag))));
      dcosjikdri[2] = ((rij[2]+rik[2])/(rijmag*rikmag)) - 
	(cosjik*((rij[2]/(rijmag*rijmag))+(rik[2]/(rikmag*rikmag))));
      dcosjikdrk[0] = (-rij[0]/(rijmag*rikmag)) + 
	(cosjik*(rik[0]/(rikmag*rikmag)));
      dcosjikdrk[1] = (-rij[1]/(rijmag*rikmag)) + 
	(cosjik*(rik[1]/(rikmag*rikmag)));
      dcosjikdrk[2] = (-rij[2]/(rijmag*rikmag)) + 
	(cosjik*(rik[2]/(rikmag*rikmag)));
      dcosjikdrj[0] = (-rik[0]/(rijmag*rikmag)) + 
	(cosjik*(rij[0]/(rijmag*rijmag)));
      dcosjikdrj[1] = (-rik[1]/(rijmag*rikmag)) + 
	(cosjik*(rij[1]/(rijmag*rijmag)));
      dcosjikdrj[2] = (-rik[2]/(rijmag*rikmag)) + 
	(cosjik*(rij[2]/(rijmag*rijmag)));

      g = gSpline(cosjik,(NijC+NijH),itype,&dgdc,&dgdN);
      tmp2 = VA*.5*(tmp*wik*dgdc*exp(lamdajik));
      fj[0] = -tmp2*dcosjikdrj[0];
      fj[1] = -tmp2*dcosjikdrj[1];
      fj[2] = -tmp2*dcosjikdrj[2];
      fi[0] = -tmp2*dcosjikdri[0];
      fi[1] = -tmp2*dcosjikdri[1];
      fi[2] = -tmp2*dcosjikdri[2];
      fk[0] = -tmp2*dcosjikdrk[0];
      fk[1] = -tmp2*dcosjikdrk[1];
      fk[2] = -tmp2*dcosjikdrk[2];

      tmp2 = VA*.5*(tmp*wik*g*exp(lamdajik)*4.0*kronecker(itype,1));
      fj[0] -= tmp2*(-rij[0]/rijmag);
      fj[1] -= tmp2*(-rij[1]/rijmag);
      fj[2] -= tmp2*(-rij[2]/rijmag);
      fi[0] -= tmp2*((-rik[0]/rikmag)+(rij[0]/rijmag));
      fi[1] -= tmp2*((-rik[1]/rikmag)+(rij[1]/rijmag));
      fi[2] -= tmp2*((-rik[2]/rikmag)+(rij[2]/rijmag));
      fk[0] -= tmp2*(rik[0]/rikmag);
      fk[1] -= tmp2*(rik[1]/rikmag);
      fk[2] -= tmp2*(rik[2]/rikmag);

      // coordination forces

      // dwik forces

      tmp2 = VA*.5*(tmp*dwik*g*exp(lamdajik))/rikmag;
      fi[0] -= tmp2*rik[0];
      fi[1] -= tmp2*rik[1];
      fi[2] -= tmp2*rik[2];
      fk[0] += tmp2*rik[0];
      fk[1] += tmp2*rik[1];
      fk[2] += tmp2*rik[2];

      // PIJ forces

      tmp2 = VA*.5*(tmp*dN2[ktype]*dwik)/rikmag;
      fi[0] -= tmp2*rik[0];
      fi[1] -= tmp2*rik[1];
      fi[2] -= tmp2*rik[2];
      fk[0] += tmp2*rik[0];
      fk[1] += tmp2*rik[1];
      fk[2] += tmp2*rik[2];

      // dgdN forces

      tmp2 = VA*.5*(tmp*tmp3*dwik)/rikmag;
      fi[0] -= tmp2*rik[0];
      fi[1] -= tmp2*rik[1];
      fi[2] -= tmp2*rik[2];
      fk[0] += tmp2*rik[0];
      fk[1] += tmp2*rik[1];
      fk[2] += tmp2*rik[2]; 

      f[atomi][0] += fi[0]; f[atomi][1] += fi[1]; f[atomi][2] += fi[2];
      f[atomj][0] += fj[0]; f[atomj][1] += fj[1]; f[atomj][2] += fj[2];
      f[atomk][0] += fk[0]; f[atomk][1] += fk[1]; f[atomk][2] += fk[2];

      if (vflag_atom) {
	rji[0] = -rij[0]; rji[1] = -rij[1]; rji[2] = -rij[2];
	rki[0] = -rik[0]; rki[1] = -rik[1]; rki[2] = -rik[2];
	v_tally3(atomi,atomj,atomk,fj,fk,rji,rki);
      }
    }
  }

  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  Etmp = 0.0;

  REBO_neighs = REBO_firstneigh[j];
  for (l = 0; l < REBO_numneigh[j]; l++) {
    atoml = REBO_neighs[l];
    if (atoml != atomi) {
      ltype = map[type[atoml]];
      rjl[0] = x[atomj][0]-x[atoml][0];
      rjl[1] = x[atomj][1]-x[atoml][1];
      rjl[2] = x[atomj][2]-x[atoml][2];
      rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
      lamdaijl = 4.0*kronecker(jtype,1) * 
	((rho[ltype][1]-rjlmag)-(rho[itype][1]-rijmag));
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dS);
      Nlj = nC[atoml]-(wjl*kronecker(jtype,0)) + 
	nH[atoml]-(wjl*kronecker(jtype,1));
      cosijl = -1.0*((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2])) / 
	(rijmag*rjlmag);
      cosijl = MIN(cosijl,1.0);
      cosijl = MAX(cosijl,-1.0);

      // evaluate splines g and derivatives dg

      g = gSpline(cosijl,NjiC+NjiH,jtype,&dgdc,&dgdN);
      Etmp = Etmp+(wjl*g*exp(lamdaijl)); 
      tmp3 = tmp3+(wjl*dgdN*exp(lamdaijl)); 
      NconjtmpJ = NconjtmpJ+(kronecker(ltype,0)*wjl*Sp(Nlj,Nmin,Nmax,dS));
    }
  }

  PjiS = 0.0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  PjiS = PijSpline(NjiC,NjiH,jtype,itype,dN2);  
  pji = pow(1.0+Etmp+PjiS,-0.5);
  tmp = -0.5*pow(pji,3.0);

  REBO_neighs = REBO_firstneigh[j];
  for (l = 0; l < REBO_numneigh[j]; l++) {
    atoml = REBO_neighs[l];
    if (atoml != atomi) {
      ltype = map[type[atoml]];
      rjl[0] = x[atomj][0]-x[atoml][0];
      rjl[1] = x[atomj][1]-x[atoml][1];
      rjl[2] = x[atomj][2]-x[atoml][2];
      rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
      lamdaijl = 4.0*kronecker(jtype,1) * 
	((rho[ltype][1]-rjlmag)-(rho[itype][1]-rijmag));
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dwjl);
      cosijl = (-1.0*((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2]))) / 
	(rijmag*rjlmag);
      cosijl = MIN(cosijl,1.0);
      cosijl = MAX(cosijl,-1.0);

      dcosijldri[0] = (-rjl[0]/(rijmag*rjlmag)) - 
	(cosijl*rij[0]/(rijmag*rijmag));
      dcosijldri[1] = (-rjl[1]/(rijmag*rjlmag)) - 
	(cosijl*rij[1]/(rijmag*rijmag));
      dcosijldri[2] = (-rjl[2]/(rijmag*rjlmag)) - 
	(cosijl*rij[2]/(rijmag*rijmag));
      dcosijldrj[0] = ((-rij[0]+rjl[0])/(rijmag*rjlmag)) + 
	(cosijl*((rij[0]/pow(rijmag,2.0))-(rjl[0]/(rjlmag*rjlmag))));
      dcosijldrj[1] = ((-rij[1]+rjl[1])/(rijmag*rjlmag)) + 
	(cosijl*((rij[1]/pow(rijmag,2.0))-(rjl[1]/(rjlmag*rjlmag))));
      dcosijldrj[2] = ((-rij[2]+rjl[2])/(rijmag*rjlmag)) + 
	(cosijl*((rij[2]/pow(rijmag,2.0))-(rjl[2]/(rjlmag*rjlmag))));
      dcosijldrl[0] = (rij[0]/(rijmag*rjlmag))+(cosijl*rjl[0]/(rjlmag*rjlmag));
      dcosijldrl[1] = (rij[1]/(rijmag*rjlmag))+(cosijl*rjl[1]/(rjlmag*rjlmag));
      dcosijldrl[2] = (rij[2]/(rijmag*rjlmag))+(cosijl*rjl[2]/(rjlmag*rjlmag));

      // evaluate splines g and derivatives dg

      g = gSpline(cosijl,NjiC+NjiH,jtype,&dgdc,&dgdN);
      tmp2 = VA*.5*(tmp*wjl*dgdc*exp(lamdaijl));
      fi[0] = -tmp2*dcosijldri[0]; 
      fi[1] = -tmp2*dcosijldri[1]; 
      fi[2] = -tmp2*dcosijldri[2]; 
      fj[0] = -tmp2*dcosijldrj[0]; 
      fj[1] = -tmp2*dcosijldrj[1]; 
      fj[2] = -tmp2*dcosijldrj[2]; 
      fl[0] = -tmp2*dcosijldrl[0]; 
      fl[1] = -tmp2*dcosijldrl[1]; 
      fl[2] = -tmp2*dcosijldrl[2]; 

      tmp2 = VA*.5*(tmp*wjl*g*exp(lamdaijl)*4.0*kronecker(jtype,1));
      fi[0] -= tmp2*(rij[0]/rijmag);
      fi[1] -= tmp2*(rij[1]/rijmag);
      fi[2] -= tmp2*(rij[2]/rijmag);
      fj[0] -= tmp2*((-rjl[0]/rjlmag)-(rij[0]/rijmag));
      fj[1] -= tmp2*((-rjl[1]/rjlmag)-(rij[1]/rijmag));
      fj[2] -= tmp2*((-rjl[2]/rjlmag)-(rij[2]/rijmag));
      fl[0] -= tmp2*(rjl[0]/rjlmag);
      fl[1] -= tmp2*(rjl[1]/rjlmag);
      fl[2] -= tmp2*(rjl[2]/rjlmag);

      // coordination forces
      
      // dwik forces

      tmp2 = VA*.5*(tmp*dwjl*g*exp(lamdaijl))/rjlmag;
      fj[0] -= tmp2*rjl[0];
      fj[1] -= tmp2*rjl[1];
      fj[2] -= tmp2*rjl[2];
      fl[0] += tmp2*rjl[0];
      fl[1] += tmp2*rjl[1];
      fl[2] += tmp2*rjl[2];

      // PIJ forces 

      tmp2 = VA*.5*(tmp*dN2[ltype]*dwjl)/rjlmag;
      fj[0] -= tmp2*rjl[0];
      fj[1] -= tmp2*rjl[1];
      fj[2] -= tmp2*rjl[2];
      fl[0] += tmp2*rjl[0];
      fl[1] += tmp2*rjl[1];
      fl[2] += tmp2*rjl[2];

      // dgdN forces

      tmp2 = VA*.5*(tmp*tmp3*dwjl)/rjlmag;
      fj[0] -= tmp2*rjl[0];
      fj[1] -= tmp2*rjl[1];
      fj[2] -= tmp2*rjl[2];
      fl[0] += tmp2*rjl[0];
      fl[1] += tmp2*rjl[1];
      fl[2] += tmp2*rjl[2];  

      f[atomi][0] += fi[0]; f[atomi][1] += fi[1]; f[atomi][2] += fi[2];
      f[atomj][0] += fj[0]; f[atomj][1] += fj[1]; f[atomj][2] += fj[2];
      f[atoml][0] += fl[0]; f[atoml][1] += fl[1]; f[atoml][2] += fl[2];

      if (vflag_atom) {
	rlj[0] = -rjl[0]; rlj[1] = -rjl[1]; rlj[2] = -rjl[2];
	v_tally3(atomi,atomj,atoml,fi,fl,rij,rlj);
      }
    }
  }	
  
  // evaluate Nij conj

  Nijconj = 1.0+(NconjtmpI*NconjtmpI)+(NconjtmpJ*NconjtmpJ);
  piRC = piRCSpline(NijC+NijH,NjiC+NjiH,Nijconj,itype,jtype,dN3);
  
  // piRC forces

  REBO_neighs_i = REBO_firstneigh[i];
  for (k = 0; k < REBO_numneigh[i]; k++) {
    atomk = REBO_neighs_i[k];
    if (atomk !=atomj) {
      ktype = map[type[atomk]];
      rik[0] = x[atomi][0]-x[atomk][0];
      rik[1] = x[atomi][1]-x[atomk][1];
      rik[2] = x[atomi][2]-x[atomk][2];
      rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
      Nki = nC[atomk]-(wik*kronecker(itype,0))+nH[atomk] - 
	(wik*kronecker(itype,1));
      SpN = Sp(Nki,Nmin,Nmax,dNki);

      tmp2 = VA*dN3[0]*dwik/rikmag;
      f[atomi][0] -= tmp2*rik[0]; 
      f[atomi][1] -= tmp2*rik[1]; 
      f[atomi][2] -= tmp2*rik[2]; 
      f[atomk][0] += tmp2*rik[0]; 
      f[atomk][1] += tmp2*rik[1]; 
      f[atomk][2] += tmp2*rik[2]; 

      if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

      tmp2 = VA*dN3[2]*(2.0*NconjtmpI*dwik*SpN)/rikmag;
      f[atomi][0] -= tmp2*rik[0]; 
      f[atomi][1] -= tmp2*rik[1]; 
      f[atomi][2] -= tmp2*rik[2]; 
      f[atomk][0] += tmp2*rik[0]; 
      f[atomk][1] += tmp2*rik[1]; 
      f[atomk][2] += tmp2*rik[2];

      if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

      if (fabs(dNki) > TOL) {
	REBO_neighs_k = REBO_firstneigh[atomk];
	for (n = 0; n < REBO_numneigh[atomk]; n++) {
	  atomn = REBO_neighs_k[n];
	  if (atomn != atomi) {
	    ntype = map[type[atomn]];
	    rkn[0] = x[atomk][0]-x[atomn][0];
	    rkn[1] = x[atomk][1]-x[atomn][1];
	    rkn[2] = x[atomk][2]-x[atomn][2];
	    rknmag = sqrt((rkn[0]*rkn[0])+(rkn[1]*rkn[1])+(rkn[2]*rkn[2]));
	    Sp(rknmag,rcmin[ktype][ntype],rcmax[ktype][ntype],dwkn);

	    tmp2 = VA*dN3[2]*(2.0*NconjtmpI*wik*dNki*dwkn)/rknmag;
	    f[atomk][0] -= tmp2*rkn[0]; 
	    f[atomk][1] -= tmp2*rkn[1]; 
	    f[atomk][2] -= tmp2*rkn[2]; 
	    f[atomn][0] += tmp2*rkn[0]; 
	    f[atomn][1] += tmp2*rkn[1]; 
	    f[atomn][2] += tmp2*rkn[2];

	    if (vflag_atom) v_tally2(atomk,atomn,-tmp2,rkn);
	  }
	}
      } 
    }
  }  
  
  // piRC forces

  REBO_neighs = REBO_firstneigh[atomj];
  for (l = 0; l < REBO_numneigh[atomj]; l++) {
    atoml = REBO_neighs[l];
    if (atoml !=atomi) {
      ltype = map[type[atoml]];
      rjl[0] = x[atomj][0]-x[atoml][0];
      rjl[1] = x[atomj][1]-x[atoml][1];
      rjl[2] = x[atomj][2]-x[atoml][2];
      rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dwjl);
      Nlj = nC[atoml]-(wjl*kronecker(jtype,0))+nH[atoml] - 
	(wjl*kronecker(jtype,1));
      SpN = Sp(Nlj,Nmin,Nmax,dNlj);

      tmp2 = VA*dN3[1]*dwjl/rjlmag;
      f[atomj][0] -= tmp2*rjl[0]; 
      f[atomj][1] -= tmp2*rjl[1]; 
      f[atomj][2] -= tmp2*rjl[2];
      f[atoml][0] += tmp2*rjl[0]; 
      f[atoml][1] += tmp2*rjl[1]; 
      f[atoml][2] += tmp2*rjl[2];

      if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

      tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*dwjl*SpN)/rjlmag;
      f[atomj][0] -= tmp2*rjl[0]; 
      f[atomj][1] -= tmp2*rjl[1]; 
      f[atomj][2] -= tmp2*rjl[2];
      f[atoml][0] += tmp2*rjl[0]; 
      f[atoml][1] += tmp2*rjl[1]; 
      f[atoml][2] += tmp2*rjl[2];

      if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

      if (fabs(dNlj) > TOL) {
	REBO_neighs_l = REBO_firstneigh[atoml];
	for (n = 0; n < REBO_numneigh[atoml]; n++) {
	  atomn = REBO_neighs_l[n];
	  if (atomn != atomj) {
	    ntype = map[type[atomn]];
	    rln[0] = x[atoml][0]-x[atomn][0];
	    rln[1] = x[atoml][1]-x[atomn][1];
	    rln[2] = x[atoml][2]-x[atomn][2];
	    rlnmag = sqrt((rln[0]*rln[0])+(rln[1]*rln[1])+(rln[2]*rln[2]));
	    Sp(rlnmag,rcmin[ltype][ntype],rcmax[ltype][ntype],dwln);

	    tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*wjl*dNlj*dwln)/rlnmag;
	    f[atoml][0] -= tmp2*rln[0]; 
	    f[atoml][1] -= tmp2*rln[1]; 
	    f[atoml][2] -= tmp2*rln[2];
	    f[atomn][0] += tmp2*rln[0]; 
	    f[atomn][1] += tmp2*rln[1]; 
	    f[atomn][2] += tmp2*rln[2];
	    
	    if (vflag_atom) v_tally2(atoml,atomn,-tmp2,rln);
	  }
	}
      } 
    }
  }   
  
  Tij = 0.0;
  dN3[0] = 0.0;
  dN3[1] = 0.0;
  dN3[2] = 0.0;
  if (itype == 0 && jtype == 0)
    Tij=TijSpline((NijC+NijH),(NjiC+NjiH),Nijconj,dN3);
  Etmp = 0.0;

  if (fabs(Tij) > TOL) {
    atom2 = atomi;
    atom3 = atomj;
    r32[0] = x[atom3][0]-x[atom2][0];
    r32[1] = x[atom3][1]-x[atom2][1];
    r32[2] = x[atom3][2]-x[atom2][2];
    r32mag = sqrt((r32[0]*r32[0])+(r32[1]*r32[1])+(r32[2]*r32[2]));
    r23[0] = -r32[0];
    r23[1] = -r32[1];
    r23[2] = -r32[2];
    r23mag = r32mag;
    REBO_neighs_i = REBO_firstneigh[i];
    for (k = 0; k < REBO_numneigh[i]; k++) {
      atomk = REBO_neighs_i[k];
      atom1 = atomk;
      ktype = map[type[atomk]];
      if (atomk != atomj) {
	r21[0] = x[atom2][0]-x[atom1][0];
	r21[1] = x[atom2][1]-x[atom1][1];
	r21[2] = x[atom2][2]-x[atom1][2];
	r21mag = sqrt(r21[0]*r21[0] + r21[1]*r21[1] + r21[2]*r21[2]);
	cos321 = -1.0*((r21[0]*r32[0])+(r21[1]*r32[1])+(r21[2]*r32[2])) / 
	  (r21mag*r32mag);
	cos321 = MIN(cos321,1.0);
	cos321 = MAX(cos321,-1.0);
	Sp2(cos321,thmin,thmax,dcut321);
	sin321 = sqrt(1.0 - cos321*cos321);
	sink2i = 1.0/(sin321*sin321);
	rik2i = 1.0/(r21mag*r21mag);
	if (sin321 != 0.0) {
	  rr = (r23mag*r23mag)-(r21mag*r21mag);
	  rjk[0] = r21[0]-r23[0];
	  rjk[1] = r21[1]-r23[1];
	  rjk[2] = r21[2]-r23[2];
	  rjk2 = (rjk[0]*rjk[0])+(rjk[1]*rjk[1])+(rjk[2]*rjk[2]);
	  rijrik = 2.0*r23mag*r21mag;
	  rik2 = r21mag*r21mag;
	  dctik = (-rr+rjk2)/(rijrik*rik2);
	  dctij = (rr+rjk2)/(rijrik*r23mag*r23mag);
	  dctjk = -2.0/rijrik;
	  w21 = Sp(r21mag,rcmin[itype][ktype],rcmaxp[itype][ktype],dw21);
	  rijmag = r32mag;
	  rikmag = r21mag;
	  rij2 = r32mag*r32mag;
	  rik2 = r21mag*r21mag;
	  costmp = 0.5*(rij2+rik2-rjk2)/rijmag/rikmag;
	  tspjik = Sp2(costmp,thmin,thmax,dtsjik);
	  dtsjik = -dtsjik;

	  REBO_neighs_j = REBO_firstneigh[j];
	  for (l = 0; l < REBO_numneigh[j]; l++) {
	    atoml = REBO_neighs_j[l];
	    atom4 = atoml;
	    ltype = map[type[atoml]];
	    if (!(atoml == atomi || atoml == atomk)) {
	      r34[0] = x[atom3][0]-x[atom4][0];
	      r34[1] = x[atom3][1]-x[atom4][1];
	      r34[2] = x[atom3][2]-x[atom4][2];
	      r34mag = sqrt((r34[0]*r34[0])+(r34[1]*r34[1])+(r34[2]*r34[2]));
	      cos234 = (r32[0]*r34[0] + r32[1]*r34[1] + r32[2]*r34[2]) / 
		(r32mag*r34mag);
	      cos234 = MIN(cos234,1.0);
	      cos234 = MAX(cos234,-1.0);
	      sin234 = sqrt(1.0 - cos234*cos234);
	      sinl2i = 1.0/(sin234*sin234);
	      rjl2i = 1.0/(r34mag*r34mag);

	      if (sin234 != 0.0) {
		w34 = Sp(r34mag,rcmin[jtype][ltype],rcmaxp[jtype][ltype],dw34);
		rr = (r23mag*r23mag)-(r34mag*r34mag);
		ril[0] = r23[0]+r34[0];
		ril[1] = r23[1]+r34[1];
		ril[2] = r23[2]+r34[2];
		ril2 = (ril[0]*ril[0])+(ril[1]*ril[1])+(ril[2]*ril[2]);
		rijrjl = 2.0*r23mag*r34mag;
		rjl2 = r34mag*r34mag;
		dctjl = (-rr+ril2)/(rijrjl*rjl2);
		dctji = (rr+ril2)/(rijrjl*r23mag*r23mag);
		dctil = -2.0/rijrjl;
		rjlmag = r34mag;
		rjl2 = r34mag*r34mag;
		costmp = 0.5*(rij2+rjl2-ril2)/rijmag/rjlmag;
		tspijl = Sp2(costmp,thmin,thmax,dtsijl);
		dtsijl = -dtsijl;
		prefactor = VA*Tij;

		cross321[0] = (r32[1]*r21[2])-(r32[2]*r21[1]);
		cross321[1] = (r32[2]*r21[0])-(r32[0]*r21[2]);
		cross321[2] = (r32[0]*r21[1])-(r32[1]*r21[0]);
		cross234[0] = (r23[1]*r34[2])-(r23[2]*r34[1]);
		cross234[1] = (r23[2]*r34[0])-(r23[0]*r34[2]);
		cross234[2] = (r23[0]*r34[1])-(r23[1]*r34[0]);

		cwnum = (cross321[0]*cross234[0]) + 
		  (cross321[1]*cross234[1]) + (cross321[2]*cross234[2]);
		cwnom = r21mag*r34mag*r23mag*r23mag*sin321*sin234;
		om1234 = cwnum/cwnom;
		cw = om1234;
		Etmp += ((1.0-pow(om1234,2.0))*w21*w34) * 
		  (1.0-tspjik)*(1.0-tspijl);
		
		dt1dik = (rik2i)-(dctik*sink2i*cos321);
		dt1djk = (-dctjk*sink2i*cos321);
		dt1djl = (rjl2i)-(dctjl*sinl2i*cos234);
		dt1dil = (-dctil*sinl2i*cos234);
		dt1dij = (2.0/(r23mag*r23mag))-(dctij*sink2i*cos321) - 
		  (dctji*sinl2i*cos234);
		
		dt2dik[0] = (-r23[2]*cross234[1])+(r23[1]*cross234[2]);
		dt2dik[1] = (-r23[0]*cross234[2])+(r23[2]*cross234[0]);
		dt2dik[2] = (-r23[1]*cross234[0])+(r23[0]*cross234[1]);
		
		dt2djl[0] = (-r23[1]*cross321[2])+(r23[2]*cross321[1]);
		dt2djl[1] = (-r23[2]*cross321[0])+(r23[0]*cross321[2]);
		dt2djl[2] = (-r23[0]*cross321[1])+(r23[1]*cross321[0]);
		
		dt2dij[0] = (r21[2]*cross234[1])-(r34[2]*cross321[1]) - 
		  (r21[1]*cross234[2])+(r34[1]*cross321[2]);
		dt2dij[1] = (r21[0]*cross234[2])-(r34[0]*cross321[2]) - 
		  (r21[2]*cross234[0])+(r34[2]*cross321[0]);
		dt2dij[2] = (r21[1]*cross234[0])-(r34[1]*cross321[0]) - 
		  (r21[0]*cross234[1])+(r34[0]*cross321[1]);
		
		aa = (prefactor*2.0*cw/cwnom)*w21*w34 * 
		  (1.0-tspjik)*(1.0-tspijl);
		aaa1 = -prefactor*(1.0-pow(om1234,2.0)) * 
		  (1.0-tspjik)*(1.0-tspijl);
		aaa2 = aaa1*w21*w34;
		at2 = aa*cwnum;
		
		fcijpc = (-dt1dij*at2)+(aaa2*dtsjik*dctij*(1.0-tspijl)) + 
		  (aaa2*dtsijl*dctji*(1.0-tspjik));
		fcikpc = (-dt1dik*at2)+(aaa2*dtsjik*dctik*(1.0-tspijl));
		fcjlpc = (-dt1djl*at2)+(aaa2*dtsijl*dctjl*(1.0-tspjik));
		fcjkpc = (-dt1djk*at2)+(aaa2*dtsjik*dctjk*(1.0-tspijl));
		fcilpc = (-dt1dil*at2)+(aaa2*dtsijl*dctil*(1.0-tspjik));
		
		F23[0] = (fcijpc*r23[0])+(aa*dt2dij[0]);
		F23[1] = (fcijpc*r23[1])+(aa*dt2dij[1]);
		F23[2] = (fcijpc*r23[2])+(aa*dt2dij[2]);
		
		F12[0] = (fcikpc*r21[0])+(aa*dt2dik[0]);
		F12[1] = (fcikpc*r21[1])+(aa*dt2dik[1]);
		F12[2] = (fcikpc*r21[2])+(aa*dt2dik[2]);
		
		F34[0] = (fcjlpc*r34[0])+(aa*dt2djl[0]);
		F34[1] = (fcjlpc*r34[1])+(aa*dt2djl[1]);
		F34[2] = (fcjlpc*r34[2])+(aa*dt2djl[2]);
		
		F31[0] = (fcjkpc*rjk[0]);
		F31[1] = (fcjkpc*rjk[1]);
		F31[2] = (fcjkpc*rjk[2]);
		
		F24[0] = (fcilpc*ril[0]);
		F24[1] = (fcilpc*ril[1]);
		F24[2] = (fcilpc*ril[2]);
		
		f1[0] = -F12[0]-F31[0];
		f1[1] = -F12[1]-F31[1];
		f1[2] = -F12[2]-F31[2];
		f2[0] = F23[0]+F12[0]+F24[0];
		f2[1] = F23[1]+F12[1]+F24[1];
		f2[2] = F23[2]+F12[2]+F24[2];
		f3[0] = -F23[0]+F34[0]+F31[0];
		f3[1] = -F23[1]+F34[1]+F31[1];
		f3[2] = -F23[2]+F34[2]+F31[2];
		f4[0] = -F34[0]-F24[0];
		f4[1] = -F34[1]-F24[1];
		f4[2] = -F34[2]-F24[2];
		
		// coordination forces

		tmp2 = VA*Tij*((1.0-(om1234*om1234))) * 
		  (1.0-tspjik)*(1.0-tspijl)*dw21*w34/r21mag;
		f2[0] -= tmp2*r21[0];
		f2[1] -= tmp2*r21[1];
		f2[2] -= tmp2*r21[2];
		f1[0] += tmp2*r21[0];
		f1[1] += tmp2*r21[1];
		f1[2] += tmp2*r21[2];
		
		tmp2 = VA*Tij*((1.0-(om1234*om1234))) * 
		  (1.0-tspjik)*(1.0-tspijl)*w21*dw34/r34mag;
		f3[0] -= tmp2*r34[0];
		f3[1] -= tmp2*r34[1];
		f3[2] -= tmp2*r34[2];
		f4[0] += tmp2*r34[0];
		f4[1] += tmp2*r34[1];
		f4[2] += tmp2*r34[2];  

		f[atom1][0] += f1[0]; f[atom1][1] += f1[1];
		f[atom1][2] += f1[2];
		f[atom2][0] += f2[0]; f[atom2][1] += f2[1];
		f[atom2][2] += f2[2];
		f[atom3][0] += f3[0]; f[atom3][1] += f3[1];
		f[atom3][2] += f3[2];
		f[atom4][0] += f4[0]; f[atom4][1] += f4[1];
		f[atom4][2] += f4[2];
		
		if (vflag_atom) {
		  r13[0] = -rjk[0]; r13[1] = -rjk[1]; r13[2] = -rjk[2];
		  r43[0] = -r34[0]; r43[1] = -r34[1]; r43[2] = -r34[2];
		  v_tally4(atom1,atom2,atom3,atom4,f1,f2,f4,r13,r23,r43);
		}
	      }
	    }
	  }
	}
      }
    }
    
    // Tij forces now that we have Etmp

    REBO_neighs = REBO_firstneigh[i];
    for (k = 0; k < REBO_numneigh[i]; k++) {
      atomk = REBO_neighs[k];
      if (atomk != atomj) {
	ktype = map[type[atomk]];
	rik[0] = x[atomi][0]-x[atomk][0];
	rik[1] = x[atomi][1]-x[atomk][1];
	rik[2] = x[atomi][2]-x[atomk][2];
	rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
	wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
	Nki = nC[atomk]-(wik*kronecker(itype,0))+nH[atomk] - 
	  (wik*kronecker(itype,1));
	SpN = Sp(Nki,Nmin,Nmax,dNki);

	tmp2 = VA*dN3[0]*dwik*Etmp/rikmag;
	f[atomi][0] -= tmp2*rik[0]; 
	f[atomi][1] -= tmp2*rik[1]; 
	f[atomi][2] -= tmp2*rik[2]; 
	f[atomk][0] += tmp2*rik[0]; 
	f[atomk][1] += tmp2*rik[1]; 
	f[atomk][2] += tmp2*rik[2]; 

	if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

	tmp2 = VA*dN3[2]*(2.0*NconjtmpI*dwik*SpN)*Etmp/rikmag;
	f[atomi][0] -= tmp2*rik[0]; 
	f[atomi][1] -= tmp2*rik[1]; 
	f[atomi][2] -= tmp2*rik[2]; 
	f[atomk][0] += tmp2*rik[0]; 
	f[atomk][1] += tmp2*rik[1]; 
	f[atomk][2] += tmp2*rik[2]; 

	if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

	if (fabs(dNki) > TOL) {
	  REBO_neighs_k = REBO_firstneigh[atomk];
	  for (n = 0; n < REBO_numneigh[atomk]; n++) {
	    atomn = REBO_neighs_k[n];
	    ntype = map[type[atomn]];
	    if (atomn != atomi) {
	      rkn[0] = x[atomk][0]-x[atomn][0];
	      rkn[1] = x[atomk][1]-x[atomn][1];
	      rkn[2] = x[atomk][2]-x[atomn][2];
	      rknmag = sqrt((rkn[0]*rkn[0])+(rkn[1]*rkn[1])+(rkn[2]*rkn[2]));
	      Sp(rknmag,rcmin[ktype][ntype],rcmax[ktype][ntype],dwkn);

	      tmp2 = VA*dN3[2]*(2.0*NconjtmpI*wik*dNki*dwkn)*Etmp/rknmag;
	      f[atomk][0] -= tmp2*rkn[0]; 
	      f[atomk][1] -= tmp2*rkn[1]; 
	      f[atomk][2] -= tmp2*rkn[2]; 
	      f[atomn][0] += tmp2*rkn[0]; 
	      f[atomn][1] += tmp2*rkn[1]; 
	      f[atomn][2] += tmp2*rkn[2]; 

	      if (vflag_atom) v_tally2(atomk,atomn,-tmp2,rkn);
	    }
	  }
	}
      }
    }

    // Tij forces

    REBO_neighs = REBO_firstneigh[j];
    for (l = 0; l < REBO_numneigh[j]; l++) {
      atoml = REBO_neighs[l];
      if (atoml != atomi) {
	ltype = map[type[atoml]];
	rjl[0] = x[atomj][0]-x[atoml][0];
	rjl[1] = x[atomj][1]-x[atoml][1];
	rjl[2] = x[atomj][2]-x[atoml][2];
	rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
	wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dwjl);
	Nlj = nC[atoml]-(wjl*kronecker(jtype,0))+nH[atoml] - 
	  (wjl*kronecker(jtype,1));
	SpN = Sp(Nlj,Nmin,Nmax,dNlj);

	tmp2 = VA*dN3[1]*dwjl*Etmp/rjlmag;
	f[atomj][0] -= tmp2*rjl[0]; 
	f[atomj][1] -= tmp2*rjl[1]; 
	f[atomj][2] -= tmp2*rjl[2]; 
	f[atoml][0] += tmp2*rjl[0]; 
	f[atoml][1] += tmp2*rjl[1]; 
	f[atoml][2] += tmp2*rjl[2]; 

	if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

	tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*dwjl*SpN)*Etmp/rjlmag;
	f[atomj][0] -= tmp2*rjl[0]; 
	f[atomj][1] -= tmp2*rjl[1]; 
	f[atomj][2] -= tmp2*rjl[2]; 
	f[atoml][0] += tmp2*rjl[0]; 
	f[atoml][1] += tmp2*rjl[1]; 
	f[atoml][2] += tmp2*rjl[2]; 

	if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

	if (fabs(dNlj) > TOL) {
	  REBO_neighs_l = REBO_firstneigh[atoml];
	  for (n = 0; n < REBO_numneigh[atoml]; n++) {
	    atomn = REBO_neighs_l[n];
	    ntype = map[type[atomn]];
	    if (atomn !=atomj) {
	      rln[0] = x[atoml][0]-x[atomn][0];
	      rln[1] = x[atoml][1]-x[atomn][1];
	      rln[2] = x[atoml][2]-x[atomn][2];
	      rlnmag = sqrt((rln[0]*rln[0])+(rln[1]*rln[1])+(rln[2]*rln[2]));
	      Sp(rlnmag,rcmin[ltype][ntype],rcmax[ltype][ntype],dwln);

	      tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*wjl*dNlj*dwln)*Etmp/rlnmag;
	      f[atoml][0] -= tmp2*rln[0];
	      f[atoml][1] -= tmp2*rln[1];
	      f[atoml][2] -= tmp2*rln[2];
	      f[atomn][0] += tmp2*rln[0];
	      f[atomn][1] += tmp2*rln[1];
	      f[atomn][2] += tmp2*rln[2];

	      if (vflag_atom) v_tally2(atoml,atomn,-tmp2,rln);
	    }
	  }
	}
      }
    }
  }

  bij = (0.5*(pij+pji))+piRC+(Tij*Etmp);
  return bij;
}

/* ----------------------------------------------------------------------
   Bij* function
------------------------------------------------------------------------- */

double PairAIREBO::bondorderLJ(int i, int j, double rij[3], double rijmag,
			       double VA, double rij0[3], double rij0mag,
			       double **f, int vflag_atom)
{
  int k,n,l,atomk,atoml,atomn,atom1,atom2,atom3,atom4;
  int atomi,atomj,itype,jtype,ktype,ltype,ntype;
  double rik[3], rjl[3], rkn[3],rknmag,dNki;
  double NijC,NijH,NjiC,NjiH,wik,dwik,dwkn,wjl;
  double rikmag,rjlmag,cosjik,cosijl,g,tmp2,tmp3;
  double Etmp,pij,tmp,wij,dwij,NconjtmpI,NconjtmpJ;
  double Nki,Nlj,dS,lamdajik,lamdaijl,dgdc,dgdN,pji,Nijconj,piRC;
  double dcosjikdri[3],dcosijldri[3],dcosjikdrk[3];
  double dN2[2],dN3[3];
  double dcosijldrj[3],dcosijldrl[3],dcosjikdrj[3],dwjl;
  double Tij,crosskij[3],crosskijmag;
  double crossijl[3],crossijlmag,omkijl;
  double tmppij,tmppji,dN2PIJ[2],dN2PJI[2],dN3piRC[3],dN3Tij[3];
  double bij,tmp3pij,tmp3pji,Stb,dStb;
  double r32[3],r32mag,cos321;
  double om1234,rln[3];
  double rlnmag,dwln,r23[3],r23mag,r21[3],r21mag;
  double w21,dw21,r34[3],r34mag,cos234,w34,dw34;
  double cross321[3],cross234[3],prefactor,SpN;
  double fcijpc,fcikpc,fcjlpc,fcjkpc,fcilpc;
  double dt2dik[3],dt2djl[3],dt2dij[3],aa,aaa1,aaa2,at2,cw,cwnum,cwnom;
  double sin321,sin234,rr,rijrik,rijrjl,rjk2,rik2,ril2,rjl2;
  double dctik,dctjk,dctjl,dctij,dctji,dctil,rik2i,rjl2i,sink2i,sinl2i;
  double rjk[3],ril[3],dt1dik,dt1djk,dt1djl,dt1dil,dt1dij;
  double dNlj;
  double PijS,PjiS;
  double rij2,tspjik,dtsjik,tspijl,dtsijl,costmp; 
  int *REBO_neighs,*REBO_neighs_i,*REBO_neighs_j,*REBO_neighs_k,*REBO_neighs_l;
  double F12[3],F23[3],F34[3],F31[3],F24[3];
  double fi[3],fj[3],fk[3],fl[3],f1[3],f2[3],f3[3],f4[4];
  double rji[3],rki[3],rlj[3],r13[3],r43[3];

  double **x = atom->x;
  int *type = atom->type;
  
  atomi = i;
  atomj = j;
  itype = map[type[atomi]];
  jtype = map[type[atomj]];
  wij = Sp(rij0mag,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
  NijC = nC[atomi]-(wij*kronecker(jtype,0));
  NijH = nH[atomi]-(wij*kronecker(jtype,1));
  NjiC = nC[atomj]-(wij*kronecker(itype,0));
  NjiH = nH[atomj]-(wij*kronecker(itype,1));
  
  rij[0] = rij0[0];
  rij[1] = rij0[1];
  rij[2] = rij0[2];
  rijmag = rij0mag; 
  
  bij = 0.0;
  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  dgdc = 0.0;
  dgdN = 0.0;
  NconjtmpI = 0.0;
  NconjtmpJ = 0.0;
  Etmp = 0.0;

  REBO_neighs = REBO_firstneigh[i];
  for (k = 0; k < REBO_numneigh[i]; k++) {
    atomk = REBO_neighs[k];
    if (atomk != atomj) {
      ktype = map[type[atomk]];
      rik[0] = x[atomi][0]-x[atomk][0];
      rik[1] = x[atomi][1]-x[atomk][1];
      rik[2] = x[atomi][2]-x[atomk][2];
      rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      lamdajik = 4.0*kronecker(itype,1) * 
	((rho[ktype][1]-rikmag)-(rho[jtype][1]-rijmag));		
      wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dS);
      Nki = nC[atomk]-(wik*kronecker(itype,0)) + 
	nH[atomk]-(wik*kronecker(itype,1));
      cosjik = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2])) / 
	(rijmag*rikmag);
      cosjik = MIN(cosjik,1.0);
      cosjik = MAX(cosjik,-1.0);

      // evaluate splines g and derivatives dg

      g = gSpline(cosjik,(NijC+NijH),itype,&dgdc,&dgdN);
      Etmp = Etmp+(wik*g*exp(lamdajik)); 
      tmp3 = tmp3+(wik*dgdN*exp(lamdajik));
      NconjtmpI = NconjtmpI+(kronecker(ktype,0)*wik*Sp(Nki,Nmin,Nmax,dS));
    }
  }

  PijS = 0.0;
  dN2PIJ[0] = 0.0;
  dN2PIJ[1] = 0.0;
  PijS = PijSpline(NijC,NijH,itype,jtype,dN2PIJ);
  pij = pow(1.0+Etmp+PijS,-0.5); 
  tmppij = -.5*pow(pij,3.0);
  tmp3pij = tmp3;
  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  Etmp = 0.0;

  REBO_neighs = REBO_firstneigh[j];
  for (l = 0; l < REBO_numneigh[j]; l++) {
    atoml = REBO_neighs[l];
    if (atoml != atomi) {
      ltype = map[type[atoml]];
      rjl[0] = x[atomj][0]-x[atoml][0];
      rjl[1] = x[atomj][1]-x[atoml][1];
      rjl[2] = x[atomj][2]-x[atoml][2];
      rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
      lamdaijl = 4.0*kronecker(jtype,1) * 
	((rho[ltype][1]-rjlmag)-(rho[itype][1]-rijmag));
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dS);
      Nlj = nC[atoml]-(wjl*kronecker(jtype,0))+nH[atoml] - 
	(wjl*kronecker(jtype,1));
      cosijl = -1.0*((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2])) / 
	(rijmag*rjlmag);
      cosijl = MIN(cosijl,1.0);
      cosijl = MAX(cosijl,-1.0);

      // evaluate splines g and derivatives dg

      g = gSpline(cosijl,NjiC+NjiH,jtype,&dgdc,&dgdN);
      Etmp = Etmp+(wjl*g*exp(lamdaijl)); 
      tmp3 = tmp3+(wjl*dgdN*exp(lamdaijl)); 
      NconjtmpJ = NconjtmpJ+(kronecker(ltype,0)*wjl*Sp(Nlj,Nmin,Nmax,dS));
    }
  }

  PjiS = 0.0;
  dN2PJI[0] = 0.0;
  dN2PJI[1] = 0.0;
  PjiS = PijSpline(NjiC,NjiH,jtype,itype,dN2PJI);
  pji = pow(1.0+Etmp+PjiS,-0.5);
  tmppji = -.5*pow(pji,3.0);
  tmp3pji = tmp3;

  // evaluate Nij conj

  Nijconj = 1.0+(NconjtmpI*NconjtmpI)+(NconjtmpJ*NconjtmpJ);
  piRC = piRCSpline(NijC+NijH,NjiC+NjiH,Nijconj,itype,jtype,dN3piRC);
  Tij = 0.0;
  dN3Tij[0] = 0.0;
  dN3Tij[1] = 0.0;
  dN3Tij[2] = 0.0;
  if (itype == 0 && jtype == 0)
    Tij=TijSpline((NijC+NijH),(NjiC+NjiH),Nijconj,dN3Tij);

  Etmp = 0.0;
  if (fabs(Tij) > TOL) {
    REBO_neighs_i = REBO_firstneigh[i];
    for (k = 0; k < REBO_numneigh[i]; k++) {
      atomk = REBO_neighs_i[k];
      ktype = map[type[atomk]];
      if (atomk != atomj) {
	rik[0] = x[atomi][0]-x[atomk][0];
	rik[1] = x[atomi][1]-x[atomk][1];
	rik[2] = x[atomi][2]-x[atomk][2];
	rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
	cos321 = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2])) / 
	  (rijmag*rikmag);
	cos321 = MIN(cos321,1.0);
	cos321 = MAX(cos321,-1.0);

	rjk[0] = rik[0]-rij[0];
	rjk[1] = rik[1]-rij[1];
	rjk[2] = rik[2]-rij[2];
	rjk2 = (rjk[0]*rjk[0])+(rjk[1]*rjk[1])+(rjk[2]*rjk[2]);
	rij2 = rijmag*rijmag;
	rik2 = rikmag*rikmag;
	costmp = 0.5*(rij2+rik2-rjk2)/rijmag/rikmag;
	tspjik = Sp2(costmp,thmin,thmax,dtsjik);

	if (sqrt(1.0 - cos321*cos321) != 0.0) {
	  wik = Sp(rikmag,rcmin[itype][ktype],rcmaxp[itype][ktype],dwik);
	  REBO_neighs_j = REBO_firstneigh[j];
	  for (l = 0; l < REBO_numneigh[j]; l++) {
	    atoml = REBO_neighs_j[l];
	    ltype = map[type[atoml]];
	    if (!(atoml == atomi || atoml == atomk)) {
	      rjl[0] = x[atomj][0]-x[atoml][0];
	      rjl[1] = x[atomj][1]-x[atoml][1];
	      rjl[2] = x[atomj][2]-x[atoml][2];
	      rjlmag = sqrt(rjl[0]*rjl[0] + rjl[1]*rjl[1] + rjl[2]*rjl[2]);
	      cos234 = -((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2])) / 
		(rijmag*rjlmag);
	      cos234 = MIN(cos234,1.0);
	      cos234 = MAX(cos234,-1.0);

	      ril[0] = rij[0]+rjl[0];
	      ril[1] = rij[1]+rjl[1];
	      ril[2] = rij[2]+rjl[2];
	      ril2 = (ril[0]*ril[0])+(ril[1]*ril[1])+(ril[2]*ril[2]);
	      rijrjl = 2.0*rijmag*rjlmag;
	      rjl2 = rjlmag*rjlmag;
	      costmp = 0.5*(rij2+rjl2-ril2)/rijmag/rjlmag;
	      tspijl = Sp2(costmp,thmin,thmax,dtsijl);

	      if (sqrt(1.0 - cos234*cos234) != 0.0) {
		wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmaxp[jtype][ltype],dS);
		crosskij[0] = (rij[1]*rik[2]-rij[2]*rik[1]);
		crosskij[1] = (rij[2]*rik[0]-rij[0]*rik[2]);
		crosskij[2] = (rij[0]*rik[1]-rij[1]*rik[0]);
		crosskijmag = sqrt(crosskij[0]*crosskij[0] + 
				   crosskij[1]*crosskij[1] + 
				   crosskij[2]*crosskij[2]);
		crossijl[0] = (rij[1]*rjl[2]-rij[2]*rjl[1]);
		crossijl[1] = (rij[2]*rjl[0]-rij[0]*rjl[2]);
		crossijl[2] = (rij[0]*rjl[1]-rij[1]*rjl[0]);
		crossijlmag = sqrt(crossijl[0]*crossijl[0] + 
				   crossijl[1]*crossijl[1] + 
				   crossijl[2]*crossijl[2]);
		omkijl = -1.0*(((crosskij[0]*crossijl[0]) + 
				(crosskij[1]*crossijl[1]) + 
				(crosskij[2]*crossijl[2])) / 
			       (crosskijmag*crossijlmag));
		Etmp += ((1.0-pow(omkijl,2.0))*wik*wjl) * 
		  (1.0-tspjik)*(1.0-tspijl);
	      }	
	    }
	  }
	}
      }
    }
  }
  
  bij = (.5*(pij+pji))+piRC+(Tij*Etmp);
  Stb = Sp2(bij,bLJmin[itype][jtype],bLJmax[itype][jtype],dStb);
  VA = VA*dStb;

  if (dStb != 0.0) {
    tmp = tmppij;
    dN2[0] = dN2PIJ[0];
    dN2[1] = dN2PIJ[1];
    tmp3 = tmp3pij;

    // pij forces

    REBO_neighs_i = REBO_firstneigh[i];
    for (k = 0; k < REBO_numneigh[i]; k++) {
      atomk = REBO_neighs_i[k];
      if (atomk != atomj) {
	lamdajik = 0.0; 
	rik[0] = x[atomi][0]-x[atomk][0];
	rik[1] = x[atomi][1]-x[atomk][1];
	rik[2] = x[atomi][2]-x[atomk][2];
	rikmag = sqrt(rik[0]*rik[0] + rik[1]*rik[1] + rik[2]*rik[2]);
	lamdajik = 4.0*kronecker(itype,1) * 
	  ((rho[ktype][1]-rikmag)-(rho[jtype][1]-rijmag));
	wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
	cosjik = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2])) / 
	  (rijmag*rikmag);
	cosjik = MIN(cosjik,1.0);
	cosjik = MAX(cosjik,-1.0);

	dcosjikdri[0] = ((rij[0]+rik[0])/(rijmag*rikmag)) - 
	  (cosjik*((rij[0]/(rijmag*rijmag))+(rik[0]/(rikmag*rikmag))));
	dcosjikdri[1] = ((rij[1]+rik[1])/(rijmag*rikmag)) - 
	  (cosjik*((rij[1]/(rijmag*rijmag))+(rik[1]/(rikmag*rikmag))));
	dcosjikdri[2] = ((rij[2]+rik[2])/(rijmag*rikmag)) - 
	  (cosjik*((rij[2]/(rijmag*rijmag))+(rik[2]/(rikmag*rikmag))));
	dcosjikdrk[0] = (-rij[0]/(rijmag*rikmag)) + 
	  (cosjik*(rik[0]/(rikmag*rikmag)));
	dcosjikdrk[1] = (-rij[1]/(rijmag*rikmag)) + 
	  (cosjik*(rik[1]/(rikmag*rikmag)));
	dcosjikdrk[2] = (-rij[2]/(rijmag*rikmag)) + 
	  (cosjik*(rik[2]/(rikmag*rikmag)));
	dcosjikdrj[0] = (-rik[0]/(rijmag*rikmag)) + 
	  (cosjik*(rij[0]/(rijmag*rijmag)));
	dcosjikdrj[1] = (-rik[1]/(rijmag*rikmag)) + 
	  (cosjik*(rij[1]/(rijmag*rijmag)));
	dcosjikdrj[2] = (-rik[2]/(rijmag*rikmag)) + 
	  (cosjik*(rij[2]/(rijmag*rijmag)));

	g = gSpline(cosjik,(NijC+NijH),itype,&dgdc,&dgdN);

	tmp2 = VA*.5*(tmp*wik*dgdc*exp(lamdajik));
	fj[0] = -tmp2*dcosjikdrj[0]; 
	fj[1] = -tmp2*dcosjikdrj[1]; 
	fj[2] = -tmp2*dcosjikdrj[2]; 
	fi[0] = -tmp2*dcosjikdri[0]; 
	fi[1] = -tmp2*dcosjikdri[1]; 
	fi[2] = -tmp2*dcosjikdri[2]; 
	fk[0] = -tmp2*dcosjikdrk[0]; 
	fk[1] = -tmp2*dcosjikdrk[1]; 
	fk[2] = -tmp2*dcosjikdrk[2]; 

	tmp2 = VA*.5*(tmp*wik*g*exp(lamdajik)*4.0*kronecker(itype,1));
	fj[0] -= tmp2*(-rij[0]/rijmag);
	fj[1] -= tmp2*(-rij[1]/rijmag);
	fj[2] -= tmp2*(-rij[2]/rijmag);
	fi[0] -= tmp2*((-rik[0]/rikmag)+(rij[0]/rijmag));
	fi[1] -= tmp2*((-rik[1]/rikmag)+(rij[1]/rijmag));
	fi[2] -= tmp2*((-rik[2]/rikmag)+(rij[2]/rijmag));
	fk[0] -= tmp2*(rik[0]/rikmag);
	fk[1] -= tmp2*(rik[1]/rikmag);
	fk[2] -= tmp2*(rik[2]/rikmag);

	// coordination forces
	
	// dwik forces

	tmp2 = VA*.5*(tmp*dwik*g*exp(lamdajik))/rikmag;
	fi[0] -= tmp2*rik[0]; 
	fi[1] -= tmp2*rik[1]; 
	fi[2] -= tmp2*rik[2]; 
	fk[0] += tmp2*rik[0]; 
	fk[1] += tmp2*rik[1]; 
	fk[2] += tmp2*rik[2]; 

	// PIJ forces

	tmp2 = VA*.5*(tmp*dN2[ktype]*dwik)/rikmag;
	fi[0] -= tmp2*rik[0]; 
	fi[1] -= tmp2*rik[1]; 
	fi[2] -= tmp2*rik[2]; 
	fk[0] += tmp2*rik[0]; 
	fk[1] += tmp2*rik[1]; 
	fk[2] += tmp2*rik[2]; 

	// dgdN forces

	tmp2 = VA*.5*(tmp*tmp3*dwik)/rikmag;
	fi[0] -= tmp2*rik[0]; 
	fi[1] -= tmp2*rik[1]; 
	fi[2] -= tmp2*rik[2]; 
	fk[0] += tmp2*rik[0]; 
	fk[1] += tmp2*rik[1]; 
	fk[2] += tmp2*rik[2];

	f[atomi][0] += fi[0]; f[atomi][1] += fi[1]; f[atomi][2] += fi[2];
	f[atomj][0] += fj[0]; f[atomj][1] += fj[1]; f[atomj][2] += fj[2];
	f[atomk][0] += fk[0]; f[atomk][1] += fk[1]; f[atomk][2] += fk[2];
	
	if (vflag_atom) {
	  rji[0] = -rij[0]; rji[1] = -rij[1]; rji[2] = -rij[2];
	  rki[0] = -rik[0]; rki[1] = -rik[1]; rki[2] = -rik[2];
	  v_tally3(atomi,atomj,atomk,fj,fk,rji,rki);
	}
      }
    }

    tmp = tmppji;
    tmp3 = tmp3pji;
    dN2[0] = dN2PJI[0];
    dN2[1] = dN2PJI[1];
    REBO_neighs  =  REBO_firstneigh[j];
    for (l = 0; l < REBO_numneigh[j]; l++) {
      atoml = REBO_neighs[l];
      if (atoml !=atomi) {
	ltype = map[type[atoml]];
	rjl[0] = x[atomj][0]-x[atoml][0];
	rjl[1] = x[atomj][1]-x[atoml][1];
	rjl[2] = x[atomj][2]-x[atoml][2];
	rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
	lamdaijl = 4.0*kronecker(jtype,1) * 
	  ((rho[ltype][1]-rjlmag)-(rho[itype][1]-rijmag));
	wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dwjl);
	cosijl = (-1.0*((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2]))) / 
	  (rijmag*rjlmag);
	cosijl = MIN(cosijl,1.0);
	cosijl = MAX(cosijl,-1.0);

	dcosijldri[0] = (-rjl[0]/(rijmag*rjlmag)) - 
	  (cosijl*rij[0]/(rijmag*rijmag));
	dcosijldri[1] = (-rjl[1]/(rijmag*rjlmag)) - 
	  (cosijl*rij[1]/(rijmag*rijmag));
	dcosijldri[2] = (-rjl[2]/(rijmag*rjlmag)) - 
	  (cosijl*rij[2]/(rijmag*rijmag));
	dcosijldrj[0] = ((-rij[0]+rjl[0])/(rijmag*rjlmag)) + 
	  (cosijl*((rij[0]/pow(rijmag,2.0))-(rjl[0]/(rjlmag*rjlmag))));
	dcosijldrj[1] = ((-rij[1]+rjl[1])/(rijmag*rjlmag)) + 
	  (cosijl*((rij[1]/pow(rijmag,2.0))-(rjl[1]/(rjlmag*rjlmag))));
	dcosijldrj[2] = ((-rij[2]+rjl[2])/(rijmag*rjlmag)) + 
	  (cosijl*((rij[2]/pow(rijmag,2.0))-(rjl[2]/(rjlmag*rjlmag))));
	dcosijldrl[0] = (rij[0]/(rijmag*rjlmag)) + 
	  (cosijl*rjl[0]/(rjlmag*rjlmag));
	dcosijldrl[1] = (rij[1]/(rijmag*rjlmag)) + 
	  (cosijl*rjl[1]/(rjlmag*rjlmag));
	dcosijldrl[2] = (rij[2]/(rijmag*rjlmag)) + 
	  (cosijl*rjl[2]/(rjlmag*rjlmag));

	// evaluate splines g and derivatives dg

	g = gSpline(cosijl,NjiC+NjiH,jtype,&dgdc,&dgdN);
	tmp2 = VA*.5*(tmp*wjl*dgdc*exp(lamdaijl));
	fi[0] = -tmp2*dcosijldri[0]; 
	fi[1] = -tmp2*dcosijldri[1]; 
	fi[2] = -tmp2*dcosijldri[2]; 
	fj[0] = -tmp2*dcosijldrj[0]; 
	fj[1] = -tmp2*dcosijldrj[1]; 
	fj[2] = -tmp2*dcosijldrj[2]; 
	fl[0] = -tmp2*dcosijldrl[0]; 
	fl[1] = -tmp2*dcosijldrl[1]; 
	fl[2] = -tmp2*dcosijldrl[2]; 

	tmp2 = VA*.5*(tmp*wjl*g*exp(lamdaijl)*4.0*kronecker(jtype,1));
	fi[0] -= tmp2*(rij[0]/rijmag);
	fi[1] -= tmp2*(rij[1]/rijmag);
	fi[2] -= tmp2*(rij[2]/rijmag);
	fj[0] -= tmp2*((-rjl[0]/rjlmag)-(rij[0]/rijmag));
	fj[1] -= tmp2*((-rjl[1]/rjlmag)-(rij[1]/rijmag));
	fj[2] -= tmp2*((-rjl[2]/rjlmag)-(rij[2]/rijmag));
	fl[0] -= tmp2*(rjl[0]/rjlmag);
	fl[1] -= tmp2*(rjl[1]/rjlmag);
	fl[2] -= tmp2*(rjl[2]/rjlmag); 

 	// coordination forces
	// dwik forces

	tmp2 = VA*.5*(tmp*dwjl*g*exp(lamdaijl))/rjlmag;
	fj[0] -= tmp2*rjl[0];
	fj[1] -= tmp2*rjl[1];
	fj[2] -= tmp2*rjl[2];
	fl[0] += tmp2*rjl[0];
	fl[1] += tmp2*rjl[1];
	fl[2] += tmp2*rjl[2];

	// PIJ forces

	tmp2 = VA*.5*(tmp*dN2[ltype]*dwjl)/rjlmag;
	fj[0] -= tmp2*rjl[0];
	fj[1] -= tmp2*rjl[1];
	fj[2] -= tmp2*rjl[2];
	fl[0] += tmp2*rjl[0];
	fl[1] += tmp2*rjl[1];
	fl[2] += tmp2*rjl[2];

	// dgdN forces

	tmp2=VA*.5*(tmp*tmp3*dwjl)/rjlmag;
	fj[0] -= tmp2*rjl[0];
	fj[1] -= tmp2*rjl[1];
	fj[2] -= tmp2*rjl[2];
	fl[0] += tmp2*rjl[0];
	fl[1] += tmp2*rjl[1];
	fl[2] += tmp2*rjl[2];

	f[atomi][0] += fi[0]; f[atomi][1] += fi[1]; f[atomi][2] += fi[2];
	f[atomj][0] += fj[0]; f[atomj][1] += fj[1]; f[atomj][2] += fj[2];
	f[atoml][0] += fl[0]; f[atoml][1] += fl[1]; f[atoml][2] += fl[2];
	
	if (vflag_atom) {
	  rlj[0] = -rjl[0]; rlj[1] = -rjl[1]; rlj[2] = -rjl[2];
	  v_tally3(atomi,atomj,atoml,fi,fl,rij,rlj);
	}
      }
    }	
    
    // piRC forces

    dN3[0] = dN3piRC[0];
    dN3[1] = dN3piRC[1];
    dN3[2] = dN3piRC[2];

    REBO_neighs_i = REBO_firstneigh[i];
    for (k = 0; k < REBO_numneigh[i]; k++) {
      atomk = REBO_neighs_i[k];
      if (atomk != atomj) {
	ktype = map[type[atomk]];
	rik[0] = x[atomi][0]-x[atomk][0];
	rik[1] = x[atomi][1]-x[atomk][1];
	rik[2] = x[atomi][2]-x[atomk][2];
	rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
	wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
	Nki = nC[atomk]-(wik*kronecker(itype,0))+nH[atomk] -
	  (wik*kronecker(itype,1));
	SpN = Sp(Nki,Nmin,Nmax,dNki);

	tmp2 = VA*dN3[0]*dwik/rikmag;
	f[atomi][0] -= tmp2*rik[0]; 
	f[atomi][1] -= tmp2*rik[1]; 
	f[atomi][2] -= tmp2*rik[2]; 
	f[atomk][0] += tmp2*rik[0]; 
	f[atomk][1] += tmp2*rik[1]; 
	f[atomk][2] += tmp2*rik[2]; 

	if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

	tmp2 = VA*dN3[2]*(2.0*NconjtmpI*dwik*SpN)/rikmag;
	f[atomi][0] -= tmp2*rik[0]; 
	f[atomi][1] -= tmp2*rik[1]; 
	f[atomi][2] -= tmp2*rik[2]; 
	f[atomk][0] += tmp2*rik[0]; 
	f[atomk][1] += tmp2*rik[1]; 
	f[atomk][2] += tmp2*rik[2];

	if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

	if (fabs(dNki) > TOL) {
	  REBO_neighs_k = REBO_firstneigh[atomk];
	  for (n = 0; n < REBO_numneigh[atomk]; n++) {
	    atomn = REBO_neighs_k[n];
	    if (atomn != atomi) {
	      ntype = map[type[atomn]];
	      rkn[0] = x[atomk][0]-x[atomn][0];
	      rkn[1] = x[atomk][1]-x[atomn][1];
	      rkn[2] = x[atomk][2]-x[atomn][2];
	      rknmag = sqrt((rkn[0]*rkn[0])+(rkn[1]*rkn[1])+(rkn[2]*rkn[2]));
	      Sp(rknmag,rcmin[ktype][ntype],rcmax[ktype][ntype],dwkn);

	      tmp2 = VA*dN3[2]*(2.0*NconjtmpI*wik*dNki*dwkn)/rknmag;
	      f[atomk][0] -= tmp2*rkn[0]; 
	      f[atomk][1] -= tmp2*rkn[1]; 
	      f[atomk][2] -= tmp2*rkn[2]; 
	      f[atomn][0] += tmp2*rkn[0]; 
	      f[atomn][1] += tmp2*rkn[1]; 
	      f[atomn][2] += tmp2*rkn[2];

	      if (vflag_atom) v_tally2(atomk,atomn,-tmp2,rkn);
	    }
	  }
	}
      }
    }

    // piRC forces to J side

    REBO_neighs = REBO_firstneigh[j];
    for (l = 0; l < REBO_numneigh[j]; l++) {
      atoml = REBO_neighs[l];
      if (atoml != atomi) {
	ltype = map[type[atoml]];
	rjl[0] = x[atomj][0]-x[atoml][0];
	rjl[1] = x[atomj][1]-x[atoml][1];
	rjl[2] = x[atomj][2]-x[atoml][2];
	rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
	wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dwjl);
	Nlj = nC[atoml]-(wjl*kronecker(jtype,0))+nH[atoml] - 
	  (wjl*kronecker(jtype,1));
	SpN = Sp(Nlj,Nmin,Nmax,dNlj);

	tmp2 = VA*dN3[1]*dwjl/rjlmag;
	f[atomj][0] -= tmp2*rjl[0]; 
	f[atomj][1] -= tmp2*rjl[1]; 
	f[atomj][2] -= tmp2*rjl[2];
	f[atoml][0] += tmp2*rjl[0]; 
	f[atoml][1] += tmp2*rjl[1]; 
	f[atoml][2] += tmp2*rjl[2];

	if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

	tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*dwjl*SpN)/rjlmag;
	f[atomj][0] -= tmp2*rjl[0]; 
	f[atomj][1] -= tmp2*rjl[1]; 
	f[atomj][2] -= tmp2*rjl[2];
	f[atoml][0] += tmp2*rjl[0]; 
	f[atoml][1] += tmp2*rjl[1]; 
	f[atoml][2] += tmp2*rjl[2];

	if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

	if (fabs(dNlj) > TOL) {
	  REBO_neighs_l = REBO_firstneigh[atoml];
	  for (n = 0; n < REBO_numneigh[atoml]; n++) {
	    atomn = REBO_neighs_l[n];
	    if (atomn != atomj) {
	      ntype = map[type[atomn]];
	      rln[0] = x[atoml][0]-x[atomn][0];
	      rln[1] = x[atoml][1]-x[atomn][1];
	      rln[2] = x[atoml][2]-x[atomn][2];
	      rlnmag = sqrt((rln[0]*rln[0])+(rln[1]*rln[1])+(rln[2]*rln[2]));
	      Sp(rlnmag,rcmin[ltype][ntype],rcmax[ltype][ntype],dwln);

	      tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*wjl*dNlj*dwln)/rlnmag;
	      f[atoml][0] -= tmp2*rln[0]; 
	      f[atoml][1] -= tmp2*rln[1]; 
	      f[atoml][2] -= tmp2*rln[2];
	      f[atomn][0] += tmp2*rln[0]; 
	      f[atomn][1] += tmp2*rln[1]; 
	      f[atomn][2] += tmp2*rln[2];

	      if (vflag_atom) v_tally2(atoml,atomn,-tmp2,rln);
	    }
	  }
	}
      }
    }

    if (fabs(Tij) > TOL) {
      dN3[0] = dN3Tij[0];
      dN3[1] = dN3Tij[1];
      dN3[2] = dN3Tij[2];
      atom2 = atomi;
      atom3 = atomj;
      r32[0] = x[atom3][0]-x[atom2][0];
      r32[1] = x[atom3][1]-x[atom2][1];
      r32[2] = x[atom3][2]-x[atom2][2];
      r32mag = sqrt((r32[0]*r32[0])+(r32[1]*r32[1])+(r32[2]*r32[2]));
      r23[0] = -r32[0];
      r23[1] = -r32[1];
      r23[2] = -r32[2];
      r23mag = r32mag;

      REBO_neighs_i = REBO_firstneigh[i];
      for (k = 0; k < REBO_numneigh[i]; k++) {
	atomk = REBO_neighs_i[k];
	atom1 = atomk;
	ktype = map[type[atomk]];
	if (atomk != atomj) {
	  r21[0] = x[atom2][0]-x[atom1][0];
	  r21[1] = x[atom2][1]-x[atom1][1];
	  r21[2] = x[atom2][2]-x[atom1][2];
	  r21mag = sqrt(r21[0]*r21[0] + r21[1]*r21[1] + r21[2]*r21[2]);
	  cos321 = ((r21[0]*rij[0])+(r21[1]*rij[1])+(r21[2]*rij[2])) / 
	    (r21mag*rijmag);
	  cos321 = MIN(cos321,1.0);
	  cos321 = MAX(cos321,-1.0);
	  sin321 = sqrt(1.0 - cos321*cos321);
	  sink2i = 1.0/(sin321*sin321);
	  rik2i = 1.0/(r21mag*r21mag);

	  if (sin321 != 0.0) { 
	    rr = (rijmag*rijmag)-(r21mag*r21mag);
	    rjk[0] = r21[0]-rij[0];
	    rjk[1] = r21[1]-rij[1];
	    rjk[2] = r21[2]-rij[2];
	    rjk2 = (rjk[0]*rjk[0])+(rjk[1]*rjk[1])+(rjk[2]*rjk[2]);
	    rijrik = 2.0*rijmag*r21mag;
	    rik2 = r21mag*r21mag;
	    dctik = (-rr+rjk2)/(rijrik*rik2);
	    dctij = (rr+rjk2)/(rijrik*rijmag*rijmag);
	    dctjk = -2.0/rijrik;
	    w21 = Sp(r21mag,rcmin[itype][ktype],rcmaxp[itype][ktype],dw21);
	    rikmag = r21mag; 
	    rij2 = r32mag*r32mag;
	    rik2 = r21mag*r21mag;
	    costmp = 0.5*(rij2+rik2-rjk2)/rijmag/rikmag;
	    tspjik = Sp2(costmp,thmin,thmax,dtsjik);
	    dtsjik = -dtsjik;
 
	    REBO_neighs_j = REBO_firstneigh[j];
	    for (l = 0; l < REBO_numneigh[j]; l++) {
	      atoml = REBO_neighs_j[l];
	      atom4 = atoml;
	      ltype = map[type[atoml]];
	      if (!(atoml == atomi || atoml == atomk)) {
		r34[0] = x[atom3][0]-x[atom4][0];
		r34[1] = x[atom3][1]-x[atom4][1];
		r34[2] = x[atom3][2]-x[atom4][2];
		r34mag = sqrt(r34[0]*r34[0] + r34[1]*r34[1] + r34[2]*r34[2]);
		cos234 = -1.0*((rij[0]*r34[0])+(rij[1]*r34[1]) + 
			       (rij[2]*r34[2]))/(rijmag*r34mag);
		cos234 = MIN(cos234,1.0);
		cos234 = MAX(cos234,-1.0);
		sin234 = sqrt(1.0 - cos234*cos234);
		sinl2i = 1.0/(sin234*sin234);
		rjl2i = 1.0/(r34mag*r34mag);

		if (sin234 != 0.0) {
		  w34 = Sp(r34mag,rcmin[jtype][ltype],
			   rcmaxp[jtype][ltype],dw34);
		  rr = (r23mag*r23mag)-(r34mag*r34mag);
		  ril[0] = r23[0]+r34[0];
		  ril[1] = r23[1]+r34[1];
		  ril[2] = r23[2]+r34[2];
		  ril2 = (ril[0]*ril[0])+(ril[1]*ril[1])+(ril[2]*ril[2]);
		  rijrjl = 2.0*r23mag*r34mag;
		  rjl2 = r34mag*r34mag;
		  dctjl = (-rr+ril2)/(rijrjl*rjl2);
		  dctji = (rr+ril2)/(rijrjl*r23mag*r23mag);
		  dctil = -2.0/rijrjl;
		  rjlmag = r34mag;
		  rjl2 = r34mag*r34mag;
		  costmp = 0.5*(rij2+rjl2-ril2)/rijmag/rjlmag;
		  tspijl = Sp2(costmp,thmin,thmax,dtsijl);
		  dtsijl = -dtsijl; //need minus sign
		  prefactor = VA*Tij;

		  cross321[0] = (r32[1]*r21[2])-(r32[2]*r21[1]);
		  cross321[1] = (r32[2]*r21[0])-(r32[0]*r21[2]);
		  cross321[2] = (r32[0]*r21[1])-(r32[1]*r21[0]);
		  cross234[0] = (r23[1]*r34[2])-(r23[2]*r34[1]);
		  cross234[1] = (r23[2]*r34[0])-(r23[0]*r34[2]);
		  cross234[2] = (r23[0]*r34[1])-(r23[1]*r34[0]);

		  cwnum = (cross321[0]*cross234[0]) + 
		    (cross321[1]*cross234[1])+(cross321[2]*cross234[2]);
		  cwnom = r21mag*r34mag*r23mag*r23mag*sin321*sin234;
		  om1234 = cwnum/cwnom;
		  cw = om1234;
		  Etmp += ((1.0-pow(om1234,2.0))*w21*w34) * 
		    (1.0-tspjik)*(1.0-tspijl);
		  
		  dt1dik = (rik2i)-(dctik*sink2i*cos321);
		  dt1djk = (-dctjk*sink2i*cos321);
		  dt1djl = (rjl2i)-(dctjl*sinl2i*cos234);
		  dt1dil = (-dctil*sinl2i*cos234);
		  dt1dij = (2.0/(r23mag*r23mag)) - 
		    (dctij*sink2i*cos321)-(dctji*sinl2i*cos234);
		  
		  dt2dik[0] = (-r23[2]*cross234[1])+(r23[1]*cross234[2]);
		  dt2dik[1] = (-r23[0]*cross234[2])+(r23[2]*cross234[0]);
		  dt2dik[2] = (-r23[1]*cross234[0])+(r23[0]*cross234[1]);
		  
		  dt2djl[0] = (-r23[1]*cross321[2])+(r23[2]*cross321[1]);
		  dt2djl[1] = (-r23[2]*cross321[0])+(r23[0]*cross321[2]);
		  dt2djl[2] = (-r23[0]*cross321[1])+(r23[1]*cross321[0]);
		  
		  dt2dij[0] = (r21[2]*cross234[1]) - 
		    (r34[2]*cross321[1])-(r21[1]*cross234[2]) + 
		    (r34[1]*cross321[2]);
		  dt2dij[1] = (r21[0]*cross234[2]) - 
		    (r34[0]*cross321[2])-(r21[2]*cross234[0]) + 
		    (r34[2]*cross321[0]);
		  dt2dij[2] = (r21[1]*cross234[0]) -
		    (r34[1]*cross321[0])-(r21[0]*cross234[1]) + 
		    (r34[0]*cross321[1]);
		  
		  aa = (prefactor*2.0*cw/cwnom)*w21*w34 * 
		    (1.0-tspjik)*(1.0-tspijl);
		  aaa1 = -prefactor*(1.0-pow(om1234,2.0)) * 
		    (1.0-tspjik)*(1.0-tspijl);
		  aaa2 = aaa1*w21*w34;
		  at2 = aa*cwnum;
		  
		  fcijpc = (-dt1dij*at2)+(aaa2*dtsjik*dctij*(1.0-tspijl)) + 
		    (aaa2*dtsijl*dctji*(1.0-tspjik));
		  fcikpc = (-dt1dik*at2)+(aaa2*dtsjik*dctik*(1.0-tspijl));
		  fcjlpc = (-dt1djl*at2)+(aaa2*dtsijl*dctjl*(1.0-tspjik));
		  fcjkpc = (-dt1djk*at2)+(aaa2*dtsjik*dctjk*(1.0-tspijl));
		  fcilpc = (-dt1dil*at2)+(aaa2*dtsijl*dctil*(1.0-tspjik));
		  
		  F23[0] = (fcijpc*r23[0])+(aa*dt2dij[0]);
		  F23[1] = (fcijpc*r23[1])+(aa*dt2dij[1]);
		  F23[2] = (fcijpc*r23[2])+(aa*dt2dij[2]);
		  
		  F12[0] = (fcikpc*r21[0])+(aa*dt2dik[0]);
		  F12[1] = (fcikpc*r21[1])+(aa*dt2dik[1]);
		  F12[2] = (fcikpc*r21[2])+(aa*dt2dik[2]);
		  
		  F34[0] = (fcjlpc*r34[0])+(aa*dt2djl[0]);
		  F34[1] = (fcjlpc*r34[1])+(aa*dt2djl[1]);
		  F34[2] = (fcjlpc*r34[2])+(aa*dt2djl[2]);
		  
		  F31[0] = (fcjkpc*rjk[0]);
		  F31[1] = (fcjkpc*rjk[1]);
		  F31[2] = (fcjkpc*rjk[2]);
		  
		  F24[0] = (fcilpc*ril[0]);
		  F24[1] = (fcilpc*ril[1]);
		  F24[2] = (fcilpc*ril[2]);
		  
		  f1[0] = -F12[0]-F31[0];
		  f1[1] = -F12[1]-F31[1];
		  f1[2] = -F12[2]-F31[2];
		  f2[0] = F23[0]+F12[0]+F24[0];
		  f2[1] = F23[1]+F12[1]+F24[1];
		  f2[2] = F23[2]+F12[2]+F24[2];
		  f3[0] = -F23[0]+F34[0]+F31[0];
		  f3[1] = -F23[1]+F34[1]+F31[1];
		  f3[2] = -F23[2]+F34[2]+F31[2];
		  f4[0] = -F34[0]-F24[0];
		  f4[1] = -F34[1]-F24[1];
		  f4[2] = -F34[2]-F24[2];
		  
		  // coordination forces

		  tmp2 = VA*Tij*((1.0-(om1234*om1234))) * 
		    (1.0-tspjik)*(1.0-tspijl)*dw21*w34/r21mag;
		  f2[0] -= tmp2*r21[0];
		  f2[1] -= tmp2*r21[1];
		  f2[2] -= tmp2*r21[2];
		  f1[0] += tmp2*r21[0];
		  f1[1] += tmp2*r21[1];
		  f1[2] += tmp2*r21[2];
		  
		  tmp2 = VA*Tij*((1.0-(om1234*om1234))) * 
		    (1.0-tspjik)*(1.0-tspijl)*w21*dw34/r34mag;
		  f3[0] -= tmp2*r34[0];
		  f3[1] -= tmp2*r34[1];
		  f3[2] -= tmp2*r34[2];
		  f4[0] += tmp2*r34[0];
		  f4[1] += tmp2*r34[1];  
		  f4[2] += tmp2*r34[2]; 

		  f[atom1][0] += f1[0]; f[atom1][1] += f1[1];
		  f[atom1][2] += f1[2];
		  f[atom2][0] += f2[0]; f[atom2][1] += f2[1];
		  f[atom2][2] += f2[2];
		  f[atom3][0] += f3[0]; f[atom3][1] += f3[1];
		  f[atom3][2] += f3[2];
		  f[atom4][0] += f4[0]; f[atom4][1] += f4[1];
		  f[atom4][2] += f4[2];
		  
		  if (vflag_atom) {
		    r13[0] = -rjk[0]; r13[1] = -rjk[1]; r13[2] = -rjk[2];
		    r43[0] = -r34[0]; r43[1] = -r34[1]; r43[2] = -r34[2];
		    v_tally4(atom1,atom2,atom3,atom4,f1,f2,f4,r13,r23,r43);
		  }
		}
	      }
	    }
	  }
	}
      }
      
      REBO_neighs = REBO_firstneigh[i];
      for (k = 0; k < REBO_numneigh[i]; k++) {
	atomk = REBO_neighs[k];
	if (atomk != atomj) {
	  ktype = map[type[atomk]];
	  rik[0] = x[atomi][0]-x[atomk][0];
	  rik[1] = x[atomi][1]-x[atomk][1];
	  rik[2] = x[atomi][2]-x[atomk][2];
	  rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
	  wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
	  Nki = nC[atomk]-(wik*kronecker(itype,0))+nH[atomk] - 
	    (wik*kronecker(itype,1));
	  SpN = Sp(Nki,Nmin,Nmax,dNki);

	  tmp2 = VA*dN3[0]*dwik*Etmp/rikmag;
	  f[atomi][0] -= tmp2*rik[0]; 
	  f[atomi][1] -= tmp2*rik[1]; 
	  f[atomi][2] -= tmp2*rik[2]; 
	  f[atomk][0] += tmp2*rik[0]; 
	  f[atomk][1] += tmp2*rik[1]; 
	  f[atomk][2] += tmp2*rik[2]; 

	  if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

	  tmp2 = VA*dN3[2]*(2.0*NconjtmpI*dwik*SpN)*Etmp/rikmag;
	  f[atomi][0] -= tmp2*rik[0]; 
	  f[atomi][1] -= tmp2*rik[1]; 
	  f[atomi][2] -= tmp2*rik[2]; 
	  f[atomk][0] += tmp2*rik[0]; 
	  f[atomk][1] += tmp2*rik[1]; 
	  f[atomk][2] += tmp2*rik[2]; 

	  if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

	  if (fabs(dNki)  >TOL) {
	    REBO_neighs_k = REBO_firstneigh[atomk];
	    for (n = 0; n < REBO_numneigh[atomk]; n++) {
	      atomn = REBO_neighs_k[n];
	      ntype = map[type[atomn]];
	      if (atomn !=atomi) {
		rkn[0] = x[atomk][0]-x[atomn][0];
		rkn[1] = x[atomk][1]-x[atomn][1];
		rkn[2] = x[atomk][2]-x[atomn][2];
		rknmag = sqrt((rkn[0]*rkn[0])+(rkn[1]*rkn[1])+(rkn[2]*rkn[2]));
		Sp(rknmag,rcmin[ktype][ntype],rcmax[ktype][ntype],dwkn);

		tmp2 = VA*dN3[2]*(2.0*NconjtmpI*wik*dNki*dwkn)*Etmp/rknmag;
		f[atomk][0] -= tmp2*rkn[0]; 
		f[atomk][1] -= tmp2*rkn[1]; 
		f[atomk][2] -= tmp2*rkn[2]; 
		f[atomn][0] += tmp2*rkn[0]; 
		f[atomn][1] += tmp2*rkn[1]; 
		f[atomn][2] += tmp2*rkn[2]; 

		if (vflag_atom) v_tally2(atomk,atomn,-tmp2,rkn);
	      }
	    }
	  }
	}
      }

      // Tij forces

      REBO_neighs = REBO_firstneigh[j];
      for (l = 0; l < REBO_numneigh[j]; l++) {
	atoml = REBO_neighs[l];
	if (atoml != atomi) {
	  ltype = map[type[atoml]];
	  rjl[0] = x[atomj][0]-x[atoml][0];
	  rjl[1] = x[atomj][1]-x[atoml][1];
	  rjl[2] = x[atomj][2]-x[atoml][2];
	  rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
	  wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dwjl);
	  Nlj = nC[atoml]-(wjl*kronecker(jtype,0))+nH[atoml] - 
	    (wjl*kronecker(jtype,1));
	  SpN = Sp(Nlj,Nmin,Nmax,dNlj);

	  tmp2 = VA*dN3[1]*dwjl*Etmp/rjlmag;
	  f[atomj][0] -= tmp2*rjl[0]; 
	  f[atomj][1] -= tmp2*rjl[1]; 
	  f[atomj][2] -= tmp2*rjl[2]; 
	  f[atoml][0] += tmp2*rjl[0]; 
	  f[atoml][1] += tmp2*rjl[1]; 
	  f[atoml][2] += tmp2*rjl[2]; 

	  if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

	  tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*dwjl*SpN)*Etmp/rjlmag;
	  f[atomj][0] -= tmp2*rjl[0]; 
	  f[atomj][1] -= tmp2*rjl[1]; 
	  f[atomj][2] -= tmp2*rjl[2]; 
	  f[atoml][0] += tmp2*rjl[0]; 
	  f[atoml][1] += tmp2*rjl[1]; 
	  f[atoml][2] += tmp2*rjl[2]; 

	  if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

	  if (fabs(dNlj) > TOL) {
	    REBO_neighs_l = REBO_firstneigh[atoml];
	    for (n = 0; n < REBO_numneigh[atoml]; n++) {
	      atomn = REBO_neighs_l[n];
	      ntype = map[type[atomn]];
	      if (atomn != atomj) {
		rln[0] = x[atoml][0]-x[atomn][0];
		rln[1] = x[atoml][1]-x[atomn][1];
		rln[2] = x[atoml][2]-x[atomn][2];
		rlnmag = sqrt((rln[0]*rln[0])+(rln[1]*rln[1])+(rln[2]*rln[2]));
		Sp(rlnmag,rcmin[ltype][ntype],rcmax[ltype][ntype],dwln);

		tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*wjl*dNlj*dwln)*Etmp/rlnmag;
		f[atoml][0] -= tmp2*rln[0];
		f[atoml][1] -= tmp2*rln[1];
		f[atoml][2] -= tmp2*rln[2];
		f[atomn][0] += tmp2*rln[0];
		f[atomn][1] += tmp2*rln[1];
		f[atomn][2] += tmp2*rln[2];

		if (vflag_atom) v_tally2(atoml,atomn,-tmp2,rln);
	      }
	    }
	  }
	}
      }
    }
  }

  return Stb;
}

/* ----------------------------------------------------------------------
   G spline
------------------------------------------------------------------------- */

double PairAIREBO::gSpline(double costh, double Nij, int typei,
			   double *dgdc, double *dgdN)
{
  double coeffs[6],dS,g1,g2,dg1,dg2,cut,g;
  int i,j;

  i = 0;
  j = 0;
  g = 0.0;
  cut = 0.0;
  dS = 0.0;
  dg1 = 0.0;
  dg2 = 0.0;
  *dgdc = 0.0;
  *dgdN = 0.0;

  // central atom is Carbon

  if (typei == 0) {
    if (costh < gCdom[0]) costh = gCdom[0];
    if (costh > gCdom[4]) costh = gCdom[4];
    if (Nij >= NCmax) {
      for (i = 0; i < 4; i++) {
	if (costh >= gCdom[i] && costh <= gCdom[i+1]) {
	  for (j = 0; j < 6; j++) coeffs[j] = gC2[i][j];
	}
      }
      g2 = Sp5th(costh,coeffs,&dg2);
      g = g2;
      *dgdc = dg2;
      *dgdN = 0.0;
    }
    if (Nij <= NCmin) {
      for (i = 0; i < 4; i++) {
	if (costh >= gCdom[i] && costh <= gCdom[i+1]) {
	  for (j = 0; j < 6; j++) coeffs[j] = gC1[i][j];
	}
      }
      g1 = Sp5th(costh,coeffs,&dg1);
      g = g1;
      *dgdc = dg1;
      *dgdN = 0.0;
    }
    if (Nij > NCmin && Nij < NCmax) {
      for (i = 0; i < 4; i++) {
	if (costh >= gCdom[i] && costh <= gCdom[i+1]) {
	  for (j = 0; j < 6; j++) coeffs[j] = gC1[i][j];
	}
      }
      g1 = Sp5th(costh,coeffs,&dg1);
      for (i = 0; i < 4; i++) {
	if (costh >= gCdom[i] && costh <= gCdom[i+1]) {
	  for (j = 0; j < 6; j++) coeffs[j] = gC2[i][j];
	}
      }
      g2 = Sp5th(costh,coeffs,&dg2);
      cut = Sp(Nij,NCmin,NCmax,dS);
      g = g2+cut*(g1-g2);
      *dgdc = dg2+(cut*(dg1-dg2));
      *dgdN = dS*(g1-g2);
    }
  }
  
  // central atom is Hydrogen

  if (typei == 1) {
    if (costh < gHdom[0]) costh = gHdom[0];
    if (costh > gHdom[3]) costh = gHdom[3];
    for (i = 0; i < 3; i++) {
      if (costh >= gHdom[i] && costh <= gHdom[i+1]) {
	for (j = 0; j < 6; j++) coeffs[j] = gH[i][j];
      }
    }
    g = Sp5th(costh,coeffs,&dg1);
    *dgdN = 0.0;
    *dgdc = dg1;
  }

  return g;
}

/* ----------------------------------------------------------------------
   Pij spline
------------------------------------------------------------------------- */

double PairAIREBO::PijSpline(double NijC, double NijH, int typei, int typej,
			     double dN2[2])
{
  int x,y,i,done;
  double Pij,coeffs[16];
  
  for (i = 0; i < 16; i++) coeffs[i]=0.0;

  x = 0;
  y = 0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  done = 0;

  // if inputs are out of bounds set them back to a point in bounds

  if (typei == 0 && typej == 0) {
    if (NijC < pCCdom[0][0] || NijC > pCCdom[0][1] || 
	NijH < pCCdom[1][0] || NijH > pCCdom[1][1]) {
      NijC=0.0;
      NijH=0.0;
    }
    if (fabs(NijC-floor(NijC)) < TOL && fabs(NijH-floor(NijH)) < TOL) {
      Pij = PCCf[(int) NijC][(int) NijH];
      dN2[0] = PCCdfdx[(int) NijC][(int) NijH];
      dN2[1] = PCCdfdy[(int) NijC][(int) NijH];
      done = 1;
    }
    if (done == 0) {
      x = (int) (floor(NijC));
      y = (int) (floor(NijH));
      for (i = 0; i<16; i++) coeffs[i] = pCC[x][y][i];
      Pij = Spbicubic(NijC,NijH,coeffs,dN2);
    }
  }

  // if inputs are out of bounds set them back to a point in bounds
  
  if (!(typei == typej)) {
    if (NijC < pCHdom[0][0] || NijC > pCHdom[0][1] || 
	NijH < pCHdom[1][0] || NijH > pCHdom[1][1]) {
      NijC = 0.0;
      NijH = 0.0;
    }	
    if (fabs(NijC-floor(NijC)) < TOL && fabs(NijH-floor(NijH)) < TOL) {
      Pij = PCHf[(int) NijC][(int) NijH];
      dN2[0] = PCHdfdx[(int) NijC][(int) NijH];
      dN2[1] = PCHdfdy[(int) NijC][(int) NijH];
      done = 1;
    }
    if (done == 0) {
      x = (int) (floor(NijC));
      y = (int) (floor(NijH));
      for (i = 0; i<16; i++) coeffs[i] = pCH[x][y][i];
      Pij = Spbicubic(NijC,NijH,coeffs,dN2);
    }
  }
  
  if (typei == 1 && typej == 1) {
    Pij = 0.0;
    dN2[0] = 0.0;
    dN2[1] = 0.0;
  }
  return Pij;
}

/* ----------------------------------------------------------------------
   PiRC spline
------------------------------------------------------------------------- */

double PairAIREBO::piRCSpline(double Nij, double Nji, double Nijconj,
			      int typei, int typej, double dN3[3])
{
  int x,y,z,i,done;
  double piRC,coeffs[64];
  x=0;
  y=0;
  z=0;
  i=0;

  done=0;

  for (i=0; i<64; i++) coeffs[i]=0.0;

  if (typei==0 && typej==0) {
    //if the inputs are out of bounds set them back to a point in bounds
    if (Nij<piCCdom[0][0]) Nij=piCCdom[0][0];
    if (Nji<piCCdom[1][0]) Nji=piCCdom[1][0];
    if (Nji>piCCdom[1][1]) Nji=0.0;
    if (Nijconj<piCCdom[2][0]) Nijconj=piCCdom[2][0];
    if (Nijconj>piCCdom[2][1]) Nijconj=piCCdom[2][1];

    if (fabs(Nij-floor(Nij))<TOL && fabs(Nji-floor(Nji))<TOL && 
	fabs(Nijconj-floor(Nijconj))<TOL) {
      piRC=piCCf[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[0]=piCCdfdx[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[1]=piCCdfdy[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[2]=piCCdfdz[(int) Nij][(int) Nji][(int) Nijconj];
      done=1;
    }

    if (done==0) {
      for (i=0; i<piCCdom[0][1]; i++)
	if (Nij>=(double) i && Nij<=(double) i+1 || Nij==(double) i) x=i;
      for (i=0; i<piCCdom[1][1]; i++)
	if (Nji>=(double) i && Nji<=(double) i+1 || Nji==(double) i) y=i;
      for (i=0; i<piCCdom[2][1]; i++)
	if (Nijconj>=(double) i && Nijconj<=(double) i+1 || 
	    Nijconj==(double) i) z=i;

      for (i=0; i<64; i++) coeffs[i]=piCC[x][y][z][i];
      piRC=Sptricubic(Nij,Nji,Nijconj,coeffs,dN3);
    }
  }
  
  
  // CH interaction

  if (typei==0 && typej==1 || typei==1 && typej==0) {
    // if the inputs are out of bounds set them back to a point in bounds

    if (Nij<piCHdom[0][0] || Nij>piCHdom[0][1] || 
	Nji<piCHdom[1][0] || Nji>piCHdom[1][1] || 
	Nijconj<piCHdom[2][0] || Nijconj>piCHdom[2][1]) {
      if (Nij<piCHdom[0][0]) Nij=piCHdom[0][0];
      if (Nij>piCHdom[0][1]) Nij=0.0;
      if (Nji<piCHdom[1][0]) Nji=piCHdom[1][0];
      if (Nji>piCHdom[1][1]) Nji=0.0;
      if (Nijconj<piCHdom[2][0]) Nijconj=piCHdom[2][0];
      if (Nijconj>piCHdom[2][1]) Nijconj=piCHdom[2][1];
    }

    if (fabs(Nij-floor(Nij))<TOL && fabs(Nji-floor(Nji))<TOL && 
	fabs(Nijconj-floor(Nijconj))<TOL) {
      piRC=piCHf[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[0]=piCHdfdx[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[1]=piCHdfdy[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[2]=piCHdfdz[(int) Nij][(int) Nji][(int) Nijconj];
      done=1;
    }

    if (done==0) {
      for (i=0; i<piCHdom[0][1]; i++)
	if (Nij>=i && Nij<=i+1) x=i;
      for (i=0; i<piCHdom[1][1]; i++)
	if (Nji>=i && Nji<=i+1) y=i;
      for (i=0; i<piCHdom[2][1]; i++)
	if (Nijconj>=i && Nijconj<=i+1) z=i;

      for (i=0; i<64; i++) coeffs[i]=piCH[x][y][z][i];
      piRC=Sptricubic(Nij,Nji,Nijconj,coeffs,dN3);
    }
  }
  
  if (typei==1 && typej==1) {
    if (Nij<piHHdom[0][0] || Nij>piHHdom[0][1] || 
	Nji<piHHdom[1][0] || Nji>piHHdom[1][1] || 
	Nijconj<piHHdom[2][0] || Nijconj>piHHdom[2][1]) {
      Nij=0.0;
      Nji=0.0;
      Nijconj=0.0;
    }
    if (fabs(Nij-floor(Nij))<TOL && fabs(Nji-floor(Nji))<TOL && 
	fabs(Nijconj-floor(Nijconj))<TOL) {
      piRC=piHHf[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[0]=piHHdfdx[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[1]=piHHdfdy[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[2]=piHHdfdz[(int) Nij][(int) Nji][(int) Nijconj];
      done=1;
    } 
    if (done==0) {
      for (i=0; i<piHHdom[0][1]; i++)
	if (Nij>=i && Nij<=i+1) x=i;
      for (i=0; i<piHHdom[1][1]; i++)
	if (Nji>=i && Nji<=i+1) y=i;
      for (i=0; i<piHHdom[2][1]; i++)
	if (Nijconj>=i && Nijconj<=i+1) z=i;

      for (i=0; i<64; i++) coeffs[i]=piHH[x][y][z][i];
      piRC=Sptricubic(Nij,Nji,Nijconj,coeffs,dN3);
    }
  }

  return piRC;
}

/* ----------------------------------------------------------------------
   Tij spline
------------------------------------------------------------------------- */

double PairAIREBO::TijSpline(double Nij, double Nji,
			     double Nijconj, double dN3[3])
{
  int x,y,z,i,done;
  double Tijf,coeffs[64];

  x=0;
  y=0;
  z=0;
  i=0;
  Tijf=0.0;
  done=0;
  for (i=0; i<64; i++) coeffs[i]=0.0;

  //if the inputs are out of bounds set them back to a point in bounds

  if (Nij<Tijdom[0][0]) Nij=Tijdom[0][0];
  if (Nij>Tijdom[0][1]) Nij=Tijdom[0][1];
  if (Nji<Tijdom[1][0]) Nji=Tijdom[1][0];
  if (Nji>Tijdom[1][1]) Nji=Tijdom[1][1];
  if (Nijconj<Tijdom[2][0]) Nijconj=Tijdom[2][0];
  if (Nijconj>Tijdom[2][1]) Nijconj=Tijdom[2][1];

  if (fabs(Nij-floor(Nij))<TOL && fabs(Nji-floor(Nji))<TOL && 
      fabs(Nijconj-floor(Nijconj))<TOL) {
    Tijf=Tf[(int) Nij][(int) Nji][(int) Nijconj];
    dN3[0]=Tdfdx[(int) Nij][(int) Nji][(int) Nijconj];
    dN3[1]=Tdfdy[(int) Nij][(int) Nji][(int) Nijconj];
    dN3[2]=Tdfdz[(int) Nij][(int) Nji][(int) Nijconj];
    done=1;
  }

  if (done==0) {
    for (i=0; i<Tijdom[0][1]; i++) 
      if (Nij>=i && Nij<=i+1) x=i;
    for (i=0; i<Tijdom[1][1]; i++)
      if (Nji>=i && Nji<=i+1) y=i;
    for (i=0; i<Tijdom[2][1]; i++) 
      if (Nijconj>=i && Nijconj<=i+1) z=i;

    for (i=0; i<64; i++) coeffs[i]=Tijc[x][y][z][i];
    Tijf=Sptricubic(Nij,Nji,Nijconj,coeffs,dN3);
  }

  return Tijf;
}

/* ----------------------------------------------------------------------
   Kronecker delta function
------------------------------------------------------------------------- */

double PairAIREBO::kronecker(int a, int b)
{
  double kd;
  if (a == b) kd = 1.0;
  else kd = 0.0;
  return kd;
}

/* ----------------------------------------------------------------------
   add pages to REBO neighbor list, starting at npage
------------------------------------------------------------------------- */

void PairAIREBO::add_pages(int npage)
{
  maxpage += PGDELTA;
  pages = (int **) 
    memory->srealloc(pages,maxpage*sizeof(int *),"AIREBO:pages");
  for (int i = npage; i < maxpage; i++)
    pages[i] = (int *) memory->smalloc(pgsize*sizeof(int),"AIREBO:pages[i]");
}

/* ----------------------------------------------------------------------
   read AIREBO potential file
------------------------------------------------------------------------- */

void PairAIREBO::read_file(char *filename)
{
  int i,j,k,l,limit;
  char s[MAXLINE];
  
  // REBO Parameters (AIREBO)
  
  double rcmin_CC,rcmin_CH,rcmin_HH,rcmax_CC,rcmax_CH,
    rcmax_HH,rcmaxp_CC,rcmaxp_CH,rcmaxp_HH;
  double Q_CC,Q_CH,Q_HH,alpha_CC,alpha_CH,alpha_HH,A_CC,A_CH,A_HH;
  double BIJc_CC1,BIJc_CC2,BIJc_CC3,BIJc_CH1,BIJc_CH2,BIJc_CH3,
    BIJc_HH1,BIJc_HH2,BIJc_HH3;
  double Beta_CC1,Beta_CC2,Beta_CC3,Beta_CH1,Beta_CH2,Beta_CH3,
    Beta_HH1,Beta_HH2,Beta_HH3;
  double rho_CC,rho_CH,rho_HH;
  
  // LJ Parameters (AIREBO)
  
  double rcLJmin_CC,rcLJmin_CH,rcLJmin_HH,rcLJmax_CC,rcLJmax_CH,
    rcLJmax_HH,bLJmin_CC;
  double bLJmin_CH,bLJmin_HH,bLJmax_CC,bLJmax_CH,bLJmax_HH,
    epsilon_CC,epsilon_CH,epsilon_HH;
  double sigma_CC,sigma_CH,sigma_HH,epsilonT_CCCC,epsilonT_CCCH,epsilonT_HCCH;
  
  MPI_Comm_rank(world,&me);
  
  // read file on proc 0
  
  if (me == 0) {
    FILE *fp = fopen(filename,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open AIREBO potential file %s",filename);
      error->one(str);
    }

    // skip initial comment lines

    while (1) {
      fgets(s,MAXLINE,fp);
      if (s[0] != '#') break;
    }
    
    // read parameters
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmin_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmin_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmin_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmax_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmax_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmax_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmaxp_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmaxp_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmaxp_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&smin);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Nmin);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Nmax);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&NCmin);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&NCmax);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Q_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Q_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Q_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&alpha_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&alpha_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&alpha_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&A_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&A_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&A_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CC1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CC2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CC3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CH1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CH2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CH3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_HH1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_HH2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_HH3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CC1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CC2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CC3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CH1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CH2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CH3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_HH1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_HH2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_HH3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rho_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rho_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rho_HH);
    
    // LJ parameters
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmin_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmin_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmin_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmax_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmax_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmax_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmin_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmin_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmin_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmax_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmax_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmax_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilon_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilon_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilon_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&sigma_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&sigma_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&sigma_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilonT_CCCC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilonT_CCCH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilonT_HCCH);
    
    // gC spline
    
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    
    // number-1 = # of domains for the spline 	
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);
    
    for (i = 0; i < limit; i++) {
      fgets(s,MAXLINE,fp);
      sscanf(s,"%lg",&gCdom[i]);
    }
    fgets(s,MAXLINE,fp);
    for (i = 0; i < limit-1; i++) {
      for (j = 0; j < 6; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&gC1[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);
    for (i = 0; i < limit-1; i++) {
      for (j = 0; j < 6; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&gC2[i][j]);
      }
    }
    
    // gH spline
    
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);	
    
    for (i = 0; i < limit; i++) {
      fgets(s,MAXLINE,fp);
      sscanf(s,"%lg",&gHdom[i]);
    }
    
    fgets(s,MAXLINE,fp);
    
    for (i = 0; i < limit-1; i++) {
      for (j = 0; j < 6; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&gH[i][j]);
      }
    }
    
    // pCC spline
    
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);	
    
    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/2; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&pCCdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);
    
    for (i = 0; i < (int) pCCdom[0][1]; i++) {
      for (j = 0; j < (int) pCCdom[1][1]; j++) {
	for (k = 0; k < 16; k++) {
	  fgets(s,MAXLINE,fp);
	  sscanf(s,"%lg",&pCC[i][j][k]);
	}
      }
    }
    
    // pCH spline
    
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);	
    
    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/2; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&pCHdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);
    
    for (i = 0; i < (int) pCCdom[0][1]; i++) {
      for (j = 0; j < (int) pCCdom[1][1]; j++) {
	for (k = 0; k < 16; k++) {
	  fgets(s,MAXLINE,fp);
	  sscanf(s,"%lg",&pCH[i][j][k]);
	}
      }
    }
    
    // piCC cpline
    
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);
    
    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/3; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&piCCdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);
    
    for (i = 0; i < (int) piCCdom[0][1]; i++) {
      for (j = 0; j < (int) piCCdom[1][1]; j++) {
	for (k = 0; k < (int) piCCdom[2][1]; k++) {
	  for (l = 0; l < 64; l = l+1) {
	    fgets(s,MAXLINE,fp);
	    sscanf(s,"%lg",&piCC[i][j][k][l]);
	  }
	}
      }
    }
    
    // piCH spline
    
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);
    
    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/3; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&piCHdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);
    
    for (i = 0; i < (int) piCHdom[0][1]; i++) {
      for (j = 0; j < (int) piCHdom[1][1]; j++) {
	for (k = 0; k < (int) piCHdom[2][1]; k++) {
	  for (l = 0; l < 64; l = l+1) {
	    fgets(s,MAXLINE,fp);
	    sscanf(s,"%lg",&piCH[i][j][k][l]);
	  }
	}
      }
    }
    
    // piHH spline
    
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);
    
    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/3; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&piHHdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);
    
    for (i = 0; i < (int) piHHdom[0][1]; i++) {
      for (j = 0; j < (int) piHHdom[1][1]; j++) {
	for (k = 0; k < (int) piHHdom[2][1]; k++) {
	  for (l = 0; l < 64; l = l+1) {
	    fgets(s,MAXLINE,fp);
	    sscanf(s,"%lg",&piHH[i][j][k][l]);
	  }
	}
      }
    }
    
    // Tij spline
    
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);
    
    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/3; j++) {
	fgets(s,MAXLINE,fp);
	sscanf(s,"%lg",&Tijdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);
    
    for (i = 0; i < (int) Tijdom[0][1]; i++) {
      for (j = 0; j < (int) Tijdom[1][1]; j++) {
	for (k = 0; k < (int) Tijdom[2][1]; k++) {
	  for (l = 0; l < 64; l = l+1) {
	    fgets(s,MAXLINE,fp);
	    sscanf(s,"%lg",&Tijc[i][j][k][l]);
	  }
	}
      }
    }
    
    fclose(fp);
  }
  
  // store read-in values in arrays
  
  if (me == 0) {
    
    // REBO
    
    rcmin[0][0] = rcmin_CC;
    rcmin[0][1] = rcmin_CH;
    rcmin[1][0] = rcmin[0][1];
    rcmin[1][1] = rcmin_HH;
    
    rcmax[0][0] = rcmax_CC;
    rcmax[0][1] = rcmax_CH;
    rcmax[1][0] = rcmax[0][1];
    rcmax[1][1] = rcmax_HH;
    
    rcmaxsq[0][0] = rcmax[0][0]*rcmax[0][0];
    rcmaxsq[1][0] = rcmax[1][0]*rcmax[1][0];
    rcmaxsq[0][1] = rcmax[0][1]*rcmax[0][1];
    rcmaxsq[1][1] = rcmax[1][1]*rcmax[1][1];
    
    rcmaxp[0][0] = rcmaxp_CC;
    rcmaxp[0][1] = rcmaxp_CH;
    rcmaxp[1][0] = rcmaxp[0][1];
    rcmaxp[1][1] = rcmaxp_HH;
    
    Q[0][0] = Q_CC;
    Q[0][1] = Q_CH;
    Q[1][0] = Q[0][1];
    Q[1][1] = Q_HH;
    
    alpha[0][0] = alpha_CC;
    alpha[0][1] = alpha_CH;
    alpha[1][0] = alpha[0][1];
    alpha[1][1] = alpha_HH;
    
    A[0][0] = A_CC;
    A[0][1] = A_CH;
    A[1][0] = A[0][1];
    A[1][1] = A_HH;
    
    rho[0][0] = rho_CC;
    rho[0][1] = rho_CH;
    rho[1][0] = rho[0][1];
    rho[1][1] = rho_HH;
    
    BIJc[0][0][0] = BIJc_CC1;
    BIJc[0][0][1] = BIJc_CC2;
    BIJc[0][0][2] = BIJc_CC3;
    BIJc[0][1][0] = BIJc_CH1;
    BIJc[0][1][1] = BIJc_CH2;
    BIJc[0][1][2] = BIJc_CH3;
    BIJc[1][0][0] = BIJc_CH1;
    BIJc[1][0][1] = BIJc_CH2;
    BIJc[1][0][2] = BIJc_CH3;
    BIJc[1][1][0] = BIJc_HH1;
    BIJc[1][1][1] = BIJc_HH2;
    BIJc[1][1][2] = BIJc_HH3;
    
    Beta[0][0][0] = Beta_CC1;
    Beta[0][0][1] = Beta_CC2;
    Beta[0][0][2] = Beta_CC3;
    Beta[0][1][0] = Beta_CH1;
    Beta[0][1][1] = Beta_CH2;
    Beta[0][1][2] = Beta_CH3;
    Beta[1][0][0] = Beta_CH1;
    Beta[1][0][1] = Beta_CH2;
    Beta[1][0][2] = Beta_CH3;
    Beta[1][1][0] = Beta_HH1;
    Beta[1][1][1] = Beta_HH2;
    Beta[1][1][2] = Beta_HH3;
    
    // LJ
    
    rcLJmin[0][0] = rcLJmin_CC;
    rcLJmin[0][1] = rcLJmin_CH;
    rcLJmin[1][0] = rcLJmin[0][1];
    rcLJmin[1][1] = rcLJmin_HH;
    
    rcLJmax[0][0] = rcLJmax_CC;
    rcLJmax[0][1] = rcLJmax_CH;
    rcLJmax[1][0] = rcLJmax[0][1];
    rcLJmax[1][1] = rcLJmax_HH;
    
    rcLJmaxsq[0][0] = rcLJmax[0][0]*rcLJmax[0][0];
    rcLJmaxsq[1][0] = rcLJmax[1][0]*rcLJmax[1][0];
    rcLJmaxsq[0][1] = rcLJmax[0][1]*rcLJmax[0][1];
    rcLJmaxsq[1][1] = rcLJmax[1][1]*rcLJmax[1][1];

    bLJmin[0][0] = bLJmin_CC;
    bLJmin[0][1] = bLJmin_CH;
    bLJmin[1][0] = bLJmin[0][1];
    bLJmin[1][1] = bLJmin_HH;
    
    bLJmax[0][0] = bLJmax_CC;
    bLJmax[0][1] = bLJmax_CH;
    bLJmax[1][0] = bLJmax[0][1];
    bLJmax[1][1] = bLJmax_HH;
    
    epsilon[0][0] = epsilon_CC;
    epsilon[0][1] = epsilon_CH;
    epsilon[1][0] = epsilon[0][1];
    epsilon[1][1] = epsilon_HH;
    
    sigma[0][0] = sigma_CC;
    sigma[0][1] = sigma_CH;
    sigma[1][0] = sigma[0][1];
    sigma[1][1] = sigma_HH;
    
    // torsional

    thmin = -1.0;
    thmax = -0.995; 
    epsilonT[0][0] = epsilonT_CCCC;
    epsilonT[0][1] = epsilonT_CCCH;
    epsilonT[1][0] = epsilonT[0][1];
    epsilonT[1][1] = epsilonT_HCCH;
  }
  
  // broadcast read-in and setup values

  MPI_Bcast(&thmin,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&thmax,1,MPI_DOUBLE,0,world);
  
  MPI_Bcast(&smin,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&Nmin,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&Nmax,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&NCmin,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&NCmax,1,MPI_DOUBLE,0,world);
  
  
  MPI_Bcast(&rcmin[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&rcmax[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&rcmaxsq[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&rcmaxp[0][0],4,MPI_DOUBLE,0,world);
  
  MPI_Bcast(&Q[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&alpha[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&A[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&rho[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&BIJc[0][0][0],12,MPI_DOUBLE,0,world);
  MPI_Bcast(&Beta[0][0][0],12,MPI_DOUBLE,0,world);
  
  MPI_Bcast(&rcLJmin[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&rcLJmax[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&rcLJmaxsq[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&rcLJmin[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&rcLJmin[0][0],4,MPI_DOUBLE,0,world);
  
  MPI_Bcast(&rcLJmin[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&rcLJmax[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&bLJmin[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&bLJmax[0][0],4,MPI_DOUBLE,0,world);
  
  MPI_Bcast(&epsilon[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&epsilonT[0][0],4,MPI_DOUBLE,0,world);
  
  MPI_Bcast(&gCdom[0],5,MPI_DOUBLE,0,world);
  MPI_Bcast(&gC1[0][0],24,MPI_DOUBLE,0,world);
  MPI_Bcast(&gC2[0][0],24,MPI_DOUBLE,0,world);
  MPI_Bcast(&gHdom[0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&gH[0][0],18,MPI_DOUBLE,0,world);
  
  MPI_Bcast(&pCCdom[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&pCHdom[0][0],4,MPI_DOUBLE,0,world);
  MPI_Bcast(&pCC[0][0][0],256,MPI_DOUBLE,0,world);
  MPI_Bcast(&pCH[0][0][0],256,MPI_DOUBLE,0,world);
  
  MPI_Bcast(&piCCdom[0][0],6,MPI_DOUBLE,0,world);
  MPI_Bcast(&piCHdom[0][0],6,MPI_DOUBLE,0,world);
  MPI_Bcast(&piHHdom[0][0],6,MPI_DOUBLE,0,world);
  MPI_Bcast(&piCC[0][0][0][0],9216,MPI_DOUBLE,0,world);
  MPI_Bcast(&piCH[0][0][0][0],9216,MPI_DOUBLE,0,world);
  MPI_Bcast(&piHH[0][0][0][0],9216,MPI_DOUBLE,0,world);
  
  MPI_Bcast(&Tijdom[0][0],6,MPI_DOUBLE,0,world);
  MPI_Bcast(&Tijc[0][0][0][0],9216,MPI_DOUBLE,0,world);
}

// ----------------------------------------------------------------------
// generic Spline functions
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   fifth order spline evaluation
------------------------------------------------------------------------- */

double PairAIREBO::Sp5th(double x, double coeffs[6], double *df)
{
  double f;
  int i;
  i = 0;
  f = 0.0;
  *df = 0.0;

  for (i = 0; i<6; i++) {
    f += coeffs[i]*pow(x,((double) i));
    if (i > 0) *df += coeffs[i]*((double) i)*pow(x,((double) i-1.0));
  }

  return f;
}

/* ----------------------------------------------------------------------
   bicubic spline evaluation
------------------------------------------------------------------------- */

double PairAIREBO::Spbicubic(double x, double y,
			     double coeffs[16], double df[2])
{
  double f;
  int i,j,cn;
  
  f = 0.0;
  df[0] = 0.0;
  df[1] = 0.0;
  cn = 0;

  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      f += coeffs[cn]*pow(x,((double) i))*pow(y,((double) j));
      if (i > 0) df[0] +=
		   (coeffs[cn]*((double) i)*pow(x,((double) i-1.0)) *
		    pow(y,((double) j)));
      if (j > 0) df[1] += 
		   (coeffs[cn]*((double) j)*pow(x,((double) i)) * 
		    pow(y,((double) j-1.0)));
      cn++;
    }
  }

  return f;
}

/* ----------------------------------------------------------------------
   tricubic spline evaluation
------------------------------------------------------------------------- */

double PairAIREBO::Sptricubic(double x, double y, double z,
			      double coeffs[64], double df[3])
{
  double f,ir,jr,kr;
  int i,j,k,cn;

  f = 0.0;
  df[0] = 0.0;
  df[1] = 0.0;
  df[2] = 0.0;
  cn = 0;

  for (i = 0; i < 4; i++) {
    ir = (double) i;
    for (j = 0; j < 4; j++) {
      jr = (double) j;
      for (k = 0; k < 4; k++) {
	kr = (double) k;
	f += (coeffs[cn]*pow(x,ir)*pow(y,jr)*pow(z,kr));
	if (i > 0) df[0] +=
		     (coeffs[cn]*ir*pow(x,ir-1.0)*pow(y,jr)*pow(z,kr));
	if (j > 0) df[1] +=
		     (coeffs[cn]*jr*pow(x,ir)*pow(y,jr-1.0)*pow(z,kr));
	if (k > 0) df[2] +=
		     (coeffs[cn]*kr*pow(x,ir)*pow(y,jr)*pow(z,kr-1.0));
	cn++;
      }
    }
  }

  return f;
}

/* ----------------------------------------------------------------------
   initialize spline knot values
------------------------------------------------------------------------- */

void PairAIREBO::spline_init()
{
  int i,j,k;
  
  for (i = 0; i < 5; i++) {
    for (j = 0; j < 5; j++) {
      PCCf[i][j] = 0.0;
      PCCdfdx[i][j] = 0.0;
      PCCdfdy[i][j] = 0.0;
      PCHf[i][j] = 0.0;
      PCHdfdx[i][j] = 0.0;
      PCHdfdy[i][j] = 0.0;
    }
  }
  
  PCCf[0][2] = -0.00050;
  PCCf[0][3] = 0.0161250;
  PCCf[1][1] = -0.010960;
  PCCf[1][2] = 0.0063260;
  PCCf[2][0] = -0.0276030;
  PCCf[2][1] = 0.003180;
  
  PCHf[0][1] = 0.2093370;
  PCHf[0][2] = -0.064450;
  PCHf[0][3] = -0.3039280;
  PCHf[1][0] = 0.010;
  PCHf[1][1] = -0.1251230;
  PCHf[1][2] = -0.2989050;
  PCHf[2][0] = -0.1220420;
  PCHf[2][1] = -0.3005290;
  PCHf[3][0] = -0.3075850;
  
  for (i = 0; i < 5; i++) {
    for (j = 0; j < 5; j++) {
      for (k = 0; k < 10; k++) {
	piCCf[i][j][k] = 0.0;
	piCCdfdx[i][j][k] = 0.0;
	piCCdfdy[i][j][k] = 0.0;
	piCCdfdz[i][j][k] = 0.0;
	piCHf[i][j][k] = 0.0;
	piCHdfdx[i][j][k] = 0.0;
	piCHdfdy[i][j][k] = 0.0;
	piCHdfdz[i][j][k] = 0.0;
	piHHf[i][j][k] = 0.0;
	piHHdfdx[i][j][k] = 0.0;
	piHHdfdy[i][j][k] = 0.0;
	piHHdfdz[i][j][k] = 0.0;
	Tf[i][j][k] = 0.0;
	Tdfdx[i][j][k] = 0.0;
	Tdfdy[i][j][k] = 0.0;
	Tdfdz[i][j][k] = 0.0;
      }
    }
  }
  
  for (i = 3; i < 10; i++) piCCf[0][0][i] = 0.0049590;
  piCCf[1][0][1] = 0.0216940;
  piCCf[0][1][1] = 0.0216940;
  for (i = 2; i < 10; i++) piCCf[1][0][i] = 0.0049590;
  for (i = 2; i < 10; i++) piCCf[0][1][i] = 0.0049590;
  piCCf[1][1][1] = 0.05250;
  piCCf[1][1][2] = -0.0020890;
  for (i = 3; i < 10; i++) piCCf[1][1][i] = -0.0080430;
  piCCf[2][0][1] = 0.0246990;
  piCCf[0][2][1] = 0.0246990;
  piCCf[2][0][2] = -0.0059710;
  piCCf[0][2][2] = -0.0059710;
  for (i = 3; i < 10; i++) piCCf[2][0][i] = 0.0049590;
  for (i = 3; i < 10; i++) piCCf[0][2][i] = 0.0049590;
  piCCf[2][1][1] = 0.0048250;
  piCCf[1][2][1] = 0.0048250;
  piCCf[2][1][2] = 0.0150;
  piCCf[1][2][2] = 0.0150;
  piCCf[2][1][3] = -0.010;
  piCCf[1][2][3] = -0.010;
  piCCf[2][1][4] = -0.0116890;
  piCCf[1][2][4] = -0.0116890;
  piCCf[2][1][5] = -0.0133780;
  piCCf[1][2][5] = -0.0133780;
  piCCf[2][1][6] = -0.0150670;
  piCCf[1][2][6] = -0.0150670;
  for (i = 7; i < 10; i++) piCCf[2][1][i] = -0.0150670;
  for (i = 7; i < 10; i++) piCCf[1][2][i] = -0.0150670;
  piCCf[2][2][1] = 0.0472250;
  piCCf[2][2][2] = 0.0110;
  piCCf[2][2][3] = 0.0198530;
  piCCf[2][2][4] = 0.0165440;
  piCCf[2][2][5] = 0.0132350;
  piCCf[2][2][6] = 0.0099260;
  piCCf[2][2][7] = 0.0066180;
  piCCf[2][2][8] = 0.0033090;
  piCCf[3][0][1] = -0.0998990;
  piCCf[0][3][1] = -0.0998990;
  piCCf[3][0][2] = -0.0998990;
  piCCf[0][3][2] = -0.0998990;
  for (i = 3; i < 10; i++) piCCf[3][0][i] = 0.0049590;
  for (i = 3; i < 10; i++) piCCf[0][3][i] = 0.0049590;
  piCCf[3][1][2] = -0.0624180;
  piCCf[1][3][2] = -0.0624180;
  for (i = 3; i < 10; i++) piCCf[3][1][i] = -0.0624180;
  for (i = 3; i < 10; i++) piCCf[1][3][i] = -0.0624180;
  piCCf[3][2][1] = -0.0223550;
  piCCf[2][3][1] = -0.0223550;
  for (i = 2; i < 10; i++) piCCf[3][2][i] = -0.0223550;
  for (i = 2; i < 10; i++) piCCf[2][3][i] = -0.0223550;
  
  piCCdfdx[2][1][1] = -0.026250;
  piCCdfdx[2][1][5] = -0.0271880;
  piCCdfdx[2][1][6] = -0.0271880;
  for (i = 7; i < 10; i++) piCCdfdx[2][1][i] = -0.0271880;
  piCCdfdx[1][3][2] = 0.0375450;
  for (i = 2; i < 10; i++) piCCdfdx[2][3][i] = 0.0624180;
  
  piCCdfdy[1][2][1] = -0.026250;
  piCCdfdy[1][2][5] = -0.0271880;
  piCCdfdy[1][2][6] = -0.0271880;
  for (i = 7; i < 10; i++) piCCdfdy[1][2][i] = -0.0271880;
  piCCdfdy[3][1][2] = 0.0375450;
  for (i = 2; i < 10; i++) piCCdfdy[3][2][i] = 0.0624180;
  
  piCCdfdz[2][1][4] = -0.0100220;
  piCCdfdz[1][2][4] = -0.0100220;
  piCCdfdz[2][1][5] = -0.0100220;
  piCCdfdz[1][2][5] = -0.0100220;
  for (i = 4; i < 10; i++) piCCdfdz[2][2][i] = -0.0033090;

  piCHf[1][1][1] = -0.050;
  piCHf[1][1][2] = -0.050;
  piCHf[1][1][3] = -0.30;
  for (i = 4; i < 10; i++) piCHf[1][1][i] = -0.050;
  for (i = 5; i < 10; i++) piCHf[2][0][i] = -0.0045240;
  for (i = 5; i < 10; i++) piCHf[0][2][i] = -0.0045240;
  piCHf[2][1][2] = -0.250;
  piCHf[1][2][2] = -0.250;
  piCHf[2][1][3] = -0.250;
  piCHf[1][2][3] = -0.250;
  piCHf[3][1][1] = -0.10;
  piCHf[1][3][1] = -0.10;
  piCHf[3][1][2] = -0.1250;
  piCHf[1][3][2] = -0.1250;
  piCHf[3][1][3] = -0.1250;
  piCHf[1][3][3] = -0.1250;
  for (i = 4; i < 10; i++) piCHf[3][1][i] = -0.10;
  for (i = 4; i < 10; i++) piCHf[1][3][i] = -0.10;
  
  piHHf[1][1][1] = 0.1249160;
  
  Tf[2][2][1] = -0.035140;
  for (i = 2; i < 10; i++) Tf[2][2][i] = -0.0040480;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

double PairAIREBO::memory_usage()
{
  double bytes = 0.0;
  bytes += maxlocal * sizeof(int);
  bytes += maxlocal * sizeof(int *);
  bytes += maxpage * neighbor->pgsize * sizeof(int);
  bytes += 3 * maxlocal * sizeof(double);
  return bytes;
}
