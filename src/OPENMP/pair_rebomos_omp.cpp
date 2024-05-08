// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   References:

   This code:
   Stewart J A and Spearot D E (2013) Atomistic simulations of nanoindentation on the basal plane of crystalline molybdenum disulfide. Modelling Simul. Mater. Sci. Eng. 21.

   Based on:
   Liang T, Phillpot S R and Sinnott S B (2009) Parameterization of a reactive many-body potential for Mo2S systems. Phys. Rev. B79 245110.
   Liang T, Phillpot S R and Sinnott S B (2012) Erratum: Parameterization of a reactive many-body potential for Mo-S systems. (Phys. Rev. B79 245110 (2009)) Phys. Rev. B85 199903(E).

   LAMMPS file contributing authors: James Stewart, Khanh Dang and Douglas Spearot (University of Arkansas)
------------------------------------------------------------------------- */

// clang-format on

#include "pair_rebomos_omp.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "math_special.h"
#include "memory.h"
#include "my_page.h"
#include "neigh_list.h"

#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;
using MathSpecial::cube;
using MathSpecial::powint;
using MathSpecial::square;

static constexpr double TOL = 1.0e-9;

/* ---------------------------------------------------------------------- */

PairREBOMoSOMP::PairREBOMoSOMP(LAMMPS *lmp) : PairREBOMoS(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

// clang-format off

/* ---------------------------------------------------------------------- */

void PairREBOMoSOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  REBO_neigh_thr();

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    FREBO_thr(ifrom,ito,eflag,thr);
    FLJ_thr(ifrom,ito,eflag,thr);

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ----------------------------------------------------------------------
   create REBO neighbor list from main neighbor list
   REBO neighbor list stores neighbors of ghost atoms
------------------------------------------------------------------------- */

void PairREBOMoSOMP::REBO_neigh_thr()
{
  const int nthreads = comm->nthreads;

  if (atom->nmax > maxlocal) {
    maxlocal = atom->nmax;
    memory->destroy(REBO_numneigh);
    memory->sfree(REBO_firstneigh);
    memory->destroy(nM);
    memory->destroy(nS);
    memory->create(REBO_numneigh,maxlocal,"REBOMoS:numneigh");
    REBO_firstneigh = (int **) memory->smalloc(maxlocal*sizeof(int *),
                                               "REBOMoS:firstneigh");
    memory->create(nM,maxlocal,"REBOMoS:nM");
    memory->create(nS,maxlocal,"REBOMoS:nS");
  }

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE
#endif
  {
    int i,j,ii,jj,n,jnum,itype,jtype;
    double xtmp,ytmp,ztmp,delx,dely,delz,rsq,dS;
    int *ilist,*jlist,*numneigh,**firstneigh;
    int *neighptr;

    double **x = atom->x;
    int *type = atom->type;

    const int allnum = list->inum + list->gnum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif

    const int iidelta = 1 + allnum/nthreads;
    const int iifrom = tid*iidelta;
    const int iito = ((iifrom+iidelta)>allnum) ? allnum : (iifrom+iidelta);

    // store all REBO neighs of owned and ghost atoms
    // scan full neighbor list of I

    // each thread has its own page allocator
    MyPage<int> &ipg = ipage[tid];
    ipg.reset();

    for (ii = iifrom; ii < iito; ii++) {
      i = ilist[ii];

      n = 0;
      neighptr = ipg.vget();

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      itype = map[type[i]];
      nM[i] = nS[i] = 0.0;
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

        if (rsq < rcmaxsq[itype][jtype]) {
          neighptr[n++] = j;
          if (jtype == 0)
            nM[i] += Sp(sqrt(rsq),rcmin[itype][jtype],rcmax[itype][jtype],dS);
          else
            nS[i] += Sp(sqrt(rsq),rcmin[itype][jtype],rcmax[itype][jtype],dS);
        }
      }

      REBO_firstneigh[i] = neighptr;
      REBO_numneigh[i] = n;
      ipg.vgot(n);
      if (ipg.status())
        error->one(FLERR,"REBO list overflow, boost neigh_modify one");
    }
  }
}

/* ----------------------------------------------------------------------
   REBO forces and energy
------------------------------------------------------------------------- */

void PairREBOMoSOMP::FREBO_thr(int ifrom, int ito, int eflag, ThrData * const thr)
{
  int i,j,k,ii,itype,jtype;
  tagint itag, jtag;
  double delx,dely,delz,evdwl,fpair,xtmp,ytmp,ztmp;
  double rsq,rij,wij;
  double Qij,Aij,alphaij,VR,pre,dVRdi,VA,bij,dVAdi,dVA;
  double dwij,del[3];
  int *ilist,*REBO_neighs;

  evdwl = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const int * const type = atom->type;
  const tagint * const tag = atom->tag;
  const int nlocal = atom->nlocal;

  ilist = list->ilist;

  // two-body interactions from REBO neighbor list, skip half of them

  for (ii = ifrom; ii < ito; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    REBO_neighs = REBO_firstneigh[i];

    for (k = 0; k < REBO_numneigh[i]; k++) {
      j = REBO_neighs[k];
      jtag = tag[j];

      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

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
      VA = -wij * BIJc[itype][jtype] * exp(-Beta[itype][jtype]*rij);

      dVA = -Beta[itype][jtype] * VA;
      dVA += VA/wij * dwij;

      del[0] = delx;
      del[1] = dely;
      del[2] = delz;
      bij = bondorder_thr(i,j,del,rij,VA,thr);
      dVAdi = bij*dVA;

      fpair = -(dVRdi+dVAdi) / rij;
      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (eflag) evdwl = VR + bij*VA;
      if (evflag) ev_tally_thr(this,i,j,nlocal,/* newton_pair */1,evdwl,0.0,fpair,delx,dely,delz,thr);
    }
  }
}

/* ----------------------------------------------------------------------
   compute LJ forces and energy
------------------------------------------------------------------------- */

void PairREBOMoSOMP::FLJ_thr(int ifrom, int ito, int eflag, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  tagint itag,jtag;
  double evdwl,fpair,xtmp,ytmp,ztmp;
  double rij,delij[3],rijsq;
  double VLJ,dVLJ;
  double vdw,dvdw;
  double r2inv,r6inv;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double c2,c3,dr,drp,r6;

  // I-J interaction from full neighbor list
  // skip 1/2 of interactions since only consider each pair once

  evdwl = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const tagint * const tag = atom->tag;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = ifrom; ii < ito; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtag = tag[j];

      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }
      jtype = map[type[j]];

      delij[0] = xtmp - x[j][0];
      delij[1] = ytmp - x[j][1];
      delij[2] = ztmp - x[j][2];
      rijsq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];
      rij = sqrt(rijsq);

      // compute LJ forces and energy

      // Outside Rmax
      if (rij > rcLJmax[itype][jtype] || rij < rcLJmin[itype][jtype]){
          VLJ = 0;
          dVLJ = 0;
      }

      // Inside Rmax and above 0.95*sigma
      else if (rij <= rcLJmax[itype][jtype] && rij >= 0.95*sigma[itype][jtype]){
              r2inv = 1.0/rijsq;
              r6inv = r2inv*r2inv*r2inv;
              VLJ = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);
              dVLJ = -r6inv*(lj1[itype][jtype]*r6inv - lj2[itype][jtype])/rij;
      }

      // Below 0.95*sigma
      else if (rij < 0.95*sigma[itype][jtype] && rij >= rcLJmin[itype][jtype]){
              dr = 0.95*sigma[itype][jtype] - rcLJmin[itype][jtype];
              r6 = powint((sigma[itype][jtype]/(0.95*sigma[itype][jtype])),6);
              vdw = 4*epsilon[itype][jtype]*r6*(r6 - 1.0);
              dvdw = (-4*epsilon[itype][jtype]/(0.95*sigma[itype][jtype]))*r6*(12.0*r6 - 6.0);
              c2 = ((3.0/dr)*vdw - dvdw)/dr;
              c3 = (vdw/(dr*dr) - c2)/dr;

              drp = rij - rcLJmin[itype][jtype];
              VLJ = drp*drp*(drp*c3 + c2);
              dVLJ = drp*(3.0*drp*c3 + 2.0*c2);
      }

      fpair = -dVLJ/rij;
      f[i][0] += delij[0]*fpair;
      f[i][1] += delij[1]*fpair;
      f[i][2] += delij[2]*fpair;
      f[j][0] -= delij[0]*fpair;
      f[j][1] -= delij[1]*fpair;
      f[j][2] -= delij[2]*fpair;

      if (eflag) evdwl = VLJ;
      if (evflag) ev_tally_thr(this,i,j,nlocal,/*newton_pair*/1,evdwl,0.0,fpair,delij[0],delij[1],delij[2],thr);

    }
  }
}

/* ----------------------------------------------------------------------
   Bij function

   The bond order term modified the attractive portion of the REBO
   potential based on the number of atoms around a specific pair
   and the bond angle between sets of three atoms.

   The functions G(cos(theta)) and P(N) are evaluated and their
   derivatives are also computed for use in the force calculation.
------------------------------------------------------------------------- */

double PairREBOMoSOMP::bondorder_thr(int i, int j, double rij[3], double rijmag, double VA, ThrData *thr)
{
  int atomi,atomj,atomk,atoml;
  int k,l;
  int itype, jtype, ktype, ltype;
  double rik[3], rjl[3], rji[3], rki[3],rlj[3], dwjl, bij;
  double NijM,NijS,NjiM,NjiS,wik,dwik,wjl;
  double rikmag,rjlmag,cosjik,cosijl,g,tmp2;
  double Etmp,pij,tmp,dwij,dS;
  double dgdc,pji;
  double dcosjikdri[3],dcosijldri[3],dcosjikdrk[3];
  double dp;
  double dcosjikdrj[3],dcosijldrj[3],dcosijldrl[3];
  double fi[3],fj[3],fk[3],fl[3];
  double PijS, PjiS;
  int *REBO_neighs;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const int * const type = atom->type;

  atomi = i;
  atomj = j;
  itype = map[type[i]];
  jtype = map[type[j]];
  Sp(rijmag,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
  NijM = nM[i];
  NijS = nS[i];
  NjiM = nM[j];
  NjiS = nS[j];
  bij = 0.0;
  tmp = 0.0;
  tmp2 = 0.0;
  dgdc = 0.0;
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
      wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dS);
      cosjik = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2])) / (rijmag*rikmag);
      cosjik = MIN(cosjik,1.0);
      cosjik = MAX(cosjik,-1.0);

      // evaluate g and derivative dg

      g = gSpline(cosjik,itype,dgdc);
      Etmp = Etmp+(wik*g);
    }
  }

  dp = 0.0;
  PijS = PijSpline(NijM,NijS,itype,dp);
  pij = 1.0/sqrt(1.0+Etmp+PijS);
  tmp = -0.5*cube(pij);

  // derivative calculations

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
      cosjik = (rij[0]*rik[0] + rij[1]*rik[1] + rij[2]*rik[2]) / (rijmag*rikmag);
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

      g = gSpline(cosjik,itype,dgdc);
      tmp2 = VA*0.5*(tmp*wik*dgdc);
      fj[0] = -tmp2*dcosjikdrj[0];
      fj[1] = -tmp2*dcosjikdrj[1];
      fj[2] = -tmp2*dcosjikdrj[2];
      fi[0] = -tmp2*dcosjikdri[0];
      fi[1] = -tmp2*dcosjikdri[1];
      fi[2] = -tmp2*dcosjikdri[2];
      fk[0] = -tmp2*dcosjikdrk[0];
      fk[1] = -tmp2*dcosjikdrk[1];
      fk[2] = -tmp2*dcosjikdrk[2];

      // coordination forces

      // dwik forces (from partial derivative)

      tmp2 = VA*0.5*(tmp*dwik*g)/rikmag;
      fi[0] -= tmp2*rik[0];
      fi[1] -= tmp2*rik[1];
      fi[2] -= tmp2*rik[2];
      fk[0] += tmp2*rik[0];
      fk[1] += tmp2*rik[1];
      fk[2] += tmp2*rik[2];

      // PIJ forces (from coordination P(N) term)

      tmp2 = VA*0.5*(tmp*dp*dwik)/rikmag;
      fi[0] -= tmp2*rik[0];
      fi[1] -= tmp2*rik[1];
      fi[2] -= tmp2*rik[2];
      fk[0] += tmp2*rik[0];
      fk[1] += tmp2*rik[1];
      fk[2] += tmp2*rik[2];

      // dgdN forces are removed

      f[atomi][0] += fi[0]; f[atomi][1] += fi[1]; f[atomi][2] += fi[2];
      f[atomj][0] += fj[0]; f[atomj][1] += fj[1]; f[atomj][2] += fj[2];
      f[atomk][0] += fk[0]; f[atomk][1] += fk[1]; f[atomk][2] += fk[2];

      if (vflag_either) {
        rji[0] = -rij[0]; rji[1] = -rij[1]; rji[2] = -rij[2];
        rki[0] = -rik[0]; rki[1] = -rik[1]; rki[2] = -rik[2];
        v_tally3_thr(this,atomi,atomj,atomk,fj,fk,rji,rki,thr);
      }
    }
  }

  // PIJ force contribution additional term
  tmp2 = -VA*0.5*(tmp*dp*dwij)/rijmag;

  f[atomi][0] += rij[0]*tmp2;
  f[atomi][1] += rij[1]*tmp2;
  f[atomi][2] += rij[2]*tmp2;
  f[atomj][0] -= rij[0]*tmp2;
  f[atomj][1] -= rij[1]*tmp2;
  f[atomj][2] -= rij[2]*tmp2;

  if (vflag_either) v_tally2_thr(this,atomi,atomj,tmp2,rij,thr);

  tmp = 0.0;
  tmp2 = 0.0;
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
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dS);
      cosijl = -1.0*((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2])) / (rijmag*rjlmag);
      cosijl = MIN(cosijl,1.0);
      cosijl = MAX(cosijl,-1.0);

      // evaluate g and derivative dg

      g = gSpline(cosijl,jtype,dgdc);
      Etmp = Etmp+(wjl*g);
    }
  }

  dp = 0.0;
  PjiS = PijSpline(NjiM,NjiS,jtype,dp);
  pji = 1.0/sqrt(1.0+Etmp+PjiS);
  tmp = -0.5*cube(pji);

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
      cosijl = (-1.0*((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2]))) / (rijmag*rjlmag);
      cosijl = MIN(cosijl,1.0);
      cosijl = MAX(cosijl,-1.0);

      dcosijldri[0] = (-rjl[0]/(rijmag*rjlmag)) -
        (cosijl*rij[0]/(rijmag*rijmag));
      dcosijldri[1] = (-rjl[1]/(rijmag*rjlmag)) -
        (cosijl*rij[1]/(rijmag*rijmag));
      dcosijldri[2] = (-rjl[2]/(rijmag*rjlmag)) -
        (cosijl*rij[2]/(rijmag*rijmag));
      dcosijldrj[0] = ((-rij[0]+rjl[0])/(rijmag*rjlmag)) +
        (cosijl*((rij[0]/square(rijmag))-(rjl[0]/(rjlmag*rjlmag))));
      dcosijldrj[1] = ((-rij[1]+rjl[1])/(rijmag*rjlmag)) +
        (cosijl*((rij[1]/square(rijmag))-(rjl[1]/(rjlmag*rjlmag))));
      dcosijldrj[2] = ((-rij[2]+rjl[2])/(rijmag*rjlmag)) +
        (cosijl*((rij[2]/square(rijmag))-(rjl[2]/(rjlmag*rjlmag))));
      dcosijldrl[0] = (rij[0]/(rijmag*rjlmag))+(cosijl*rjl[0]/(rjlmag*rjlmag));
      dcosijldrl[1] = (rij[1]/(rijmag*rjlmag))+(cosijl*rjl[1]/(rjlmag*rjlmag));
      dcosijldrl[2] = (rij[2]/(rijmag*rjlmag))+(cosijl*rjl[2]/(rjlmag*rjlmag));

      // evaluate g and derivatives dg

      g = gSpline(cosijl,jtype,dgdc);
      tmp2 = VA*0.5*(tmp*wjl*dgdc);
      fi[0] = -tmp2*dcosijldri[0];
      fi[1] = -tmp2*dcosijldri[1];
      fi[2] = -tmp2*dcosijldri[2];
      fj[0] = -tmp2*dcosijldrj[0];
      fj[1] = -tmp2*dcosijldrj[1];
      fj[2] = -tmp2*dcosijldrj[2];
      fl[0] = -tmp2*dcosijldrl[0];
      fl[1] = -tmp2*dcosijldrl[1];
      fl[2] = -tmp2*dcosijldrl[2];

      // coordination forces

      // dwik forces (from partial derivative)

      tmp2 = VA*0.5*(tmp*dwjl*g)/rjlmag;
      fj[0] -= tmp2*rjl[0];
      fj[1] -= tmp2*rjl[1];
      fj[2] -= tmp2*rjl[2];
      fl[0] += tmp2*rjl[0];
      fl[1] += tmp2*rjl[1];
      fl[2] += tmp2*rjl[2];

      // PIJ forces (coordination)

      tmp2 = VA*0.5*(tmp*dp*dwjl)/rjlmag;
      fj[0] -= tmp2*rjl[0];
      fj[1] -= tmp2*rjl[1];
      fj[2] -= tmp2*rjl[2];
      fl[0] += tmp2*rjl[0];
      fl[1] += tmp2*rjl[1];
      fl[2] += tmp2*rjl[2];

      // dgdN forces are removed

      f[atomi][0] += fi[0]; f[atomi][1] += fi[1]; f[atomi][2] += fi[2];
      f[atomj][0] += fj[0]; f[atomj][1] += fj[1]; f[atomj][2] += fj[2];
      f[atoml][0] += fl[0]; f[atoml][1] += fl[1]; f[atoml][2] += fl[2];

      if (vflag_either) {
        rlj[0] = -rjl[0]; rlj[1] = -rjl[1]; rlj[2] = -rjl[2];
        v_tally3_thr(this,atomi,atomj,atoml,fi,fl,rij,rlj,thr);
      }
    }
  }

  // PIJ force contribution additional term

  tmp2 = -VA*0.5*(tmp*dp*dwij)/rijmag;
  f[atomi][0] += rij[0]*tmp2;
  f[atomi][1] += rij[1]*tmp2;
  f[atomi][2] += rij[2]*tmp2;
  f[atomj][0] -= rij[0]*tmp2;
  f[atomj][1] -= rij[1]*tmp2;
  f[atomj][2] -= rij[2]*tmp2;

  if (vflag_either) v_tally2_thr(this,atomi,atomj,tmp2,rij,thr);

  bij = (0.5*(pij+pji));
  return bij;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairREBOMoSOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairREBOMoS::memory_usage();

  return bytes;
}
