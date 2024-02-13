// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_airebo_omp.h"

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
using namespace MathSpecial;

static constexpr double TOL = 1.0e-9;

/* ---------------------------------------------------------------------- */

PairAIREBOOMP::PairAIREBOOMP(LAMMPS *lmp) :
  PairAIREBO(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairAIREBOOMP::compute(int eflag, int vflag)
{
  double pv0=0.0,pv1=0.0,pv2=0.0;

  ev_init(eflag,vflag);

  REBO_neigh_thr();

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag) reduction(+:pv0,pv1,pv2)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    FREBO_thr(ifrom,ito,eflag,&pv0,thr);
    if (ljflag) FLJ_thr(ifrom,ito,eflag,&pv1,thr);
    if (torflag) TORSION_thr(ifrom,ito,eflag,&pv2,thr);

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region

  pvector[0] = pv0;
  pvector[1] = pv1;
  pvector[2] = pv2;
}

/* ----------------------------------------------------------------------
   create REBO neighbor list from main neighbor list
   REBO neighbor list stores neighbors of ghost atoms
------------------------------------------------------------------------- */

void PairAIREBOOMP::REBO_neigh_thr()
{
  const int nthreads = comm->nthreads;

  if (atom->nmax > maxlocal) {
    maxlocal = atom->nmax;
    memory->destroy(REBO_numneigh);
    memory->sfree(REBO_firstneigh);
    memory->destroy(nC);
    memory->destroy(nH);
    memory->create(REBO_numneigh,maxlocal,"AIREBO:numneigh");
    REBO_firstneigh = (int **) memory->smalloc(maxlocal*sizeof(int *),
                                               "AIREBO:firstneigh");
    memory->create(nC,maxlocal,"AIREBO:nC");
    memory->create(nH,maxlocal,"AIREBO:nH");
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
      nC[i] = nH[i] = 0.0;
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
            nC[i] += Sp(sqrt(rsq),rcmin[itype][jtype],rcmax[itype][jtype],dS);
          else
            nH[i] += Sp(sqrt(rsq),rcmin[itype][jtype],rcmax[itype][jtype],dS);
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

void PairAIREBOOMP::FREBO_thr(int ifrom, int ito, int eflag, double *pv0, ThrData * const thr)
{
  int i,j,k,m,ii,itype,jtype;
  tagint itag,jtag;
  double delx,dely,delz,evdwl,fpair,xtmp,ytmp,ztmp;
  double rsq,rij,wij;
  double Qij,Aij,alphaij,VR,pre,dVRdi,VA,term,bij,dVAdi,dVA;
  double dwij,del[3];
  int *ilist,*REBO_neighs;

  evdwl = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  int *type = atom->type;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

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
      for (m = 0; m < 3; m++) {
        term = -wij * BIJc[itype][jtype][m] * exp(-Beta[itype][jtype][m]*rij);
        VA += term;
        dVA += -Beta[itype][jtype][m] * term;
      }
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

      if (eflag) *pv0 += evdwl = VR + bij*VA;
      if (evflag) ev_tally_thr(this,i,j,nlocal,/* newton_pair */ 1,
                               evdwl,0.0,fpair,delx,dely,delz,thr);
    }
  }
}

/* ----------------------------------------------------------------------
   compute LJ forces and energy
   find 3- and 4-step paths between atoms I,J via REBO neighbor lists
------------------------------------------------------------------------- */

void PairAIREBOOMP::FLJ_thr(int ifrom, int ito, int eflag, double *pv1, ThrData * const thr)
{
  int i,j,k,m,ii,jj,kk,mm,jnum,itype,jtype,ktype,mtype;
  tagint itag,jtag;
  int atomi,atomj,atomk,atomm;
  int testpath,npath,done;
  double evdwl,fpair,xtmp,ytmp,ztmp;
  double rsq,best,wik,wkm,cij,rij,dwij,dwik,dwkj,dwkm,dwmj;
  double delij[3],rijsq,delik[3],rik,deljk[3];
  double rkj,wkj,dC,VLJ,dVLJ,VA,Str,dStr,Stb;
  double vdw,slw,dvdw,dslw,drij,swidth,tee,tee2;
  double rljmin,rljmax;
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
  rljmin = 0.0;
  rljmax = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = ifrom; ii < ito; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    atomi = i;
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
      atomj = j;

      delij[0] = xtmp - x[j][0];
      delij[1] = ytmp - x[j][1];
      delij[2] = ztmp - x[j][2];
      rijsq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];

      // if outside of LJ cutoff, skip
      // if outside of 4-path cutoff, best = 0.0, no need to test paths
      // if outside of 2-path cutoff but inside 4-path cutoff,
      //   best = 0.0, test 3-,4-paths
      // if inside 2-path cutoff, best = wij, only test 3-,4-paths if best < 1
      npath = testpath = done = 0;
      best = 0.0;

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
          } else { dwik = wik = 0.0; rikS = rik = 1.0; }

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
              } else { dwkm = wkm = 0.0; rkmS = rkm = 1.0; }

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

      rljmin = sigma[itype][jtype];
      rljmax = sigcut * rljmin;
      rljmin = sigmin * rljmin;

      if (rij > rljmax) {
        slw = 0.0;
        dslw = 0.0;
      } else if (rij > rljmin) {
        drij = rij - rljmin;
        swidth = rljmax - rljmin;
        tee = drij / swidth;
        tee2 = tee*tee;
        slw = 1.0 - tee2 * (3.0 - 2.0 * tee);
        dslw = -6.0 * tee * (1.0 - tee) / swidth;
      } else {
        slw = 1.0;
        dslw = 0.0;
      }

      if (morseflag) {

        const double exr = exp(-rij*lj4[itype][jtype]);
        vdw = lj1[itype][jtype]*exr*(lj2[itype][jtype]*exr - 2);
        dvdw = lj3[itype][jtype]*exr*(1-lj2[itype][jtype]*exr);

      } else {

        r2inv = 1.0/rijsq;
        r6inv = r2inv*r2inv*r2inv;

        vdw = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);
        dvdw = -r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]) / rij;
      }

      // VLJ now becomes vdw * slw, derivaties, etc.

      VLJ = vdw * slw;
      dVLJ = dvdw * slw + vdw * dslw;

      Str = Sp2(rij,rcLJmin[itype][jtype],rcLJmax[itype][jtype],dStr);
      VA = Str*cij*VLJ;
      if (Str > 0.0) {
        scale = rcmin[itype][jtype] / rij;
        delscale[0] = scale * delij[0];
        delscale[1] = scale * delij[1];
        delscale[2] = scale * delij[2];
        Stb = bondorderLJ_thr(i,j,delscale,rcmin[itype][jtype],VA,delij,rij,thr);
      } else Stb = 0.0;

      fpair = -(dStr * (Stb*cij*VLJ - cij*VLJ) +
                dVLJ * (Str*Stb*cij + cij - Str*cij)) / rij;

      f[i][0] += delij[0]*fpair;
      f[i][1] += delij[1]*fpair;
      f[i][2] += delij[2]*fpair;
      f[j][0] -= delij[0]*fpair;
      f[j][1] -= delij[1]*fpair;
      f[j][2] -= delij[2]*fpair;

      if (eflag) *pv1 += evdwl = VA*Stb + (1.0-Str)*cij*VLJ;
      if (evflag) ev_tally_thr(this,i,j,nlocal,/* newton_pair */ 1,
                               evdwl,0.0,fpair,delij[0],delij[1],delij[2],thr);

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

          if (vflag_either) v_tally2_thr(this,atomi,atomj,fpair,delij,thr);

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

          if (vflag_either)
            v_tally3_thr(this,atomi,atomj,atomk,fi,fj,delikS,deljkS,thr);

        } else if (npath == 4) {
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

          if (vflag_either) {
            delimS[0] = delikS[0] + delkmS[0];
            delimS[1] = delikS[1] + delkmS[1];
            delimS[2] = delikS[2] + delkmS[2];
            v_tally4_thr(this,atomi,atomj,atomk,atomm,fi,fj,fk,delimS,deljmS,delkmS,thr);
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   torsional forces and energy
------------------------------------------------------------------------- */

void PairAIREBOOMP::TORSION_thr(int ifrom, int ito, int eflag, double *pv2, ThrData * const thr)
{
  int i,j,k,l,ii;
  tagint itag,jtag;
  double evdwl,fpair,xtmp,ytmp,ztmp;
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

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  int *type = atom->type;
  tagint *tag = atom->tag;

  ilist = list->ilist;

  for (ii = ifrom; ii < ito; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    if (itype != 0) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    REBO_neighs_i = REBO_firstneigh[i];

    for (jj = 0; jj < REBO_numneigh[i]; jj++) {
      j = REBO_neighs_i[jj];
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
          Vtors = (Ec*(powint(cw2,5)))-(ekijl/10.0);

          if (eflag) *pv2 += evdwl = Vtors*w21*w23*w34*(1.0-tspjik)*(1.0-tspijl);

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
          dvpdcw = (-1.0)*Ec*(-.5)*5.0*powint(cw2,4) *
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
            ev_tally4_thr(this,i,j,k,l,evdwl,fi,fj,fk,delil,del34,delkl,thr);
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Bij function
------------------------------------------------------------------------- */

double PairAIREBOOMP::bondorder_thr(int i, int j, double rij[3], double rijmag,
                                    double VA, ThrData * const thr)
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
  double dt2dik[3],dt2djl[3],dt2dij[3],aa,aaa2,at2,cw,cwnum,cwnom;
  double sin321,sin234,rr,rijrik,rijrjl,rjk2,rik2,ril2,rjl2;
  double dctik,dctjk,dctjl,dctij,dctji,dctil,rik2i,rjl2i,sink2i,sinl2i;
  double rjk[3],ril[3],dt1dik,dt1djk,dt1djl,dt1dil,dt1dij;
  double F23[3],F12[3],F34[3],F31[3],F24[3],fi[3],fj[3],fk[3],fl[3];
  double f1[3],f2[3],f3[3],f4[4];
  double dcut321,PijS,PjiS;
  double rij2,tspjik,dtsjik,tspijl,dtsijl,costmp;
  int *REBO_neighs,*REBO_neighs_i,*REBO_neighs_j,*REBO_neighs_k,*REBO_neighs_l;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
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
  pij = 1.0/sqrt(1.0+Etmp+PijS);
  tmp = -0.5*pij*pij*pij;

  const double invrijm = 1.0/rijmag;
  const double invrijm2 = invrijm*invrijm;

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

      const double invrikm = 1.0/rikmag;
      const double invrijkm = invrijm*invrikm;
      const double invrikm2 = invrikm*invrikm;

      cosjik = (rij[0]*rik[0] + rij[1]*rik[1] + rij[2]*rik[2]) * invrijkm;
      cosjik = MIN(cosjik,1.0);
      cosjik = MAX(cosjik,-1.0);

      dcosjikdri[0] = ((rij[0]+rik[0])*invrijkm) -
        (cosjik*((rij[0]*invrijm2)+(rik[0]*invrikm2)));
      dcosjikdri[1] = ((rij[1]+rik[1])*invrijkm) -
        (cosjik*((rij[1]*invrijm2)+(rik[1]*invrikm2)));
      dcosjikdri[2] = ((rij[2]+rik[2])*invrijkm) -
        (cosjik*((rij[2]*invrijm2)+(rik[2]*invrikm2)));
      dcosjikdrk[0] = (-rij[0]*invrijkm) + (cosjik*(rik[0]*invrikm2));
      dcosjikdrk[1] = (-rij[1]*invrijkm) + (cosjik*(rik[1]*invrikm2));
      dcosjikdrk[2] = (-rij[2]*invrijkm) + (cosjik*(rik[2]*invrikm2));
      dcosjikdrj[0] = (-rik[0]*invrijkm) + (cosjik*(rij[0]*invrijm2));
      dcosjikdrj[1] = (-rik[1]*invrijkm) + (cosjik*(rij[1]*invrijm2));
      dcosjikdrj[2] = (-rik[2]*invrijkm) + (cosjik*(rij[2]*invrijm2));

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
      fj[0] -= tmp2*(-rij[0]*invrijm);
      fj[1] -= tmp2*(-rij[1]*invrijm);
      fj[2] -= tmp2*(-rij[2]*invrijm);
      fi[0] -= tmp2*((-rik[0]*invrikm)+(rij[0]*invrijm));
      fi[1] -= tmp2*((-rik[1]*invrikm)+(rij[1]*invrijm));
      fi[2] -= tmp2*((-rik[2]*invrikm)+(rij[2]*invrijm));
      fk[0] -= tmp2*(rik[0]*invrikm);
      fk[1] -= tmp2*(rik[1]*invrikm);
      fk[2] -= tmp2*(rik[2]*invrikm);

      // coordination forces

      // dwik forces

      tmp2 = VA*.5*(tmp*dwik*g*exp(lamdajik))*invrikm;
      fi[0] -= tmp2*rik[0];
      fi[1] -= tmp2*rik[1];
      fi[2] -= tmp2*rik[2];
      fk[0] += tmp2*rik[0];
      fk[1] += tmp2*rik[1];
      fk[2] += tmp2*rik[2];

      // PIJ forces

      tmp2 = VA*.5*(tmp*dN2[ktype]*dwik)*invrikm;
      fi[0] -= tmp2*rik[0];
      fi[1] -= tmp2*rik[1];
      fi[2] -= tmp2*rik[2];
      fk[0] += tmp2*rik[0];
      fk[1] += tmp2*rik[1];
      fk[2] += tmp2*rik[2];

      // dgdN forces

      tmp2 = VA*.5*(tmp*tmp3*dwik)*invrikm;
      fi[0] -= tmp2*rik[0];
      fi[1] -= tmp2*rik[1];
      fi[2] -= tmp2*rik[2];
      fk[0] += tmp2*rik[0];
      fk[1] += tmp2*rik[1];
      fk[2] += tmp2*rik[2];

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
  pji = 1.0/sqrt(1.0+Etmp+PjiS);
  tmp = -0.5*pji*pji*pji;

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

      const double invrjlm = 1.0/rjlmag;
      const double invrijlm = invrijm*invrjlm;
      const double invrjlm2 = invrjlm*invrjlm;

      cosijl = (-1.0*((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2])))
        * invrijlm;

      cosijl = MIN(cosijl,1.0);
      cosijl = MAX(cosijl,-1.0);

      dcosijldri[0] = (-rjl[0]*invrijlm) - (cosijl*rij[0]*invrijm2);
      dcosijldri[1] = (-rjl[1]*invrijlm) - (cosijl*rij[1]*invrijm2);
      dcosijldri[2] = (-rjl[2]*invrijlm) - (cosijl*rij[2]*invrijm2);
      dcosijldrj[0] = ((-rij[0]+rjl[0])*invrijlm) +
        (cosijl*((rij[0]*invrijm2)-(rjl[0]*invrjlm2)));
      dcosijldrj[1] = ((-rij[1]+rjl[1])*invrijlm) +
        (cosijl*((rij[1]*invrijm2)-(rjl[1]*invrjlm2)));
      dcosijldrj[2] = ((-rij[2]+rjl[2])*invrijlm) +
        (cosijl*((rij[2]*invrijm2)-(rjl[2]*invrjlm2)));
      dcosijldrl[0] = (rij[0]*invrijlm)+(cosijl*rjl[0]*invrjlm2);
      dcosijldrl[1] = (rij[1]*invrijlm)+(cosijl*rjl[1]*invrjlm2);
      dcosijldrl[2] = (rij[2]*invrijlm)+(cosijl*rjl[2]*invrjlm2);

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
      fi[0] -= tmp2*(rij[0]*invrijm);
      fi[1] -= tmp2*(rij[1]*invrijm);
      fi[2] -= tmp2*(rij[2]*invrijm);
      fj[0] -= tmp2*((-rjl[0]*invrjlm)-(rij[0]*invrijm));
      fj[1] -= tmp2*((-rjl[1]*invrjlm)-(rij[1]*invrijm));
      fj[2] -= tmp2*((-rjl[2]*invrjlm)-(rij[2]*invrijm));
      fl[0] -= tmp2*(rjl[0]*invrjlm);
      fl[1] -= tmp2*(rjl[1]*invrjlm);
      fl[2] -= tmp2*(rjl[2]*invrjlm);

      // coordination forces

      // dwik forces

      tmp2 = VA*.5*(tmp*dwjl*g*exp(lamdaijl))*invrjlm;
      fj[0] -= tmp2*rjl[0];
      fj[1] -= tmp2*rjl[1];
      fj[2] -= tmp2*rjl[2];
      fl[0] += tmp2*rjl[0];
      fl[1] += tmp2*rjl[1];
      fl[2] += tmp2*rjl[2];

      // PIJ forces

      tmp2 = VA*.5*(tmp*dN2[ltype]*dwjl)*invrjlm;
      fj[0] -= tmp2*rjl[0];
      fj[1] -= tmp2*rjl[1];
      fj[2] -= tmp2*rjl[2];
      fl[0] += tmp2*rjl[0];
      fl[1] += tmp2*rjl[1];
      fl[2] += tmp2*rjl[2];

      // dgdN forces

      tmp2 = VA*.5*(tmp*tmp3*dwjl)*invrjlm;
      fj[0] -= tmp2*rjl[0];
      fj[1] -= tmp2*rjl[1];
      fj[2] -= tmp2*rjl[2];
      fl[0] += tmp2*rjl[0];
      fl[1] += tmp2*rjl[1];
      fl[2] += tmp2*rjl[2];

      f[atomi][0] += fi[0]; f[atomi][1] += fi[1]; f[atomi][2] += fi[2];
      f[atomj][0] += fj[0]; f[atomj][1] += fj[1]; f[atomj][2] += fj[2];
      f[atoml][0] += fl[0]; f[atoml][1] += fl[1]; f[atoml][2] += fl[2];

      if (vflag_either) {
        rlj[0] = -rjl[0]; rlj[1] = -rjl[1]; rlj[2] = -rjl[2];
        v_tally3_thr(this,atomi,atomj,atoml,fi,fl,rij,rlj,thr);
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

      if (vflag_either) v_tally2_thr(this,atomi,atomk,-tmp2,rik,thr);

      // due to kronecker(ktype, 0) term in contribution
      // to NconjtmpI and later Nijconj
      if (ktype != 0) continue;

      tmp2 = VA*dN3[2]*(2.0*NconjtmpI*dwik*SpN)/rikmag;
      f[atomi][0] -= tmp2*rik[0];
      f[atomi][1] -= tmp2*rik[1];
      f[atomi][2] -= tmp2*rik[2];
      f[atomk][0] += tmp2*rik[0];
      f[atomk][1] += tmp2*rik[1];
      f[atomk][2] += tmp2*rik[2];

      if (vflag_either) v_tally2_thr(this,atomi,atomk,-tmp2,rik,thr);

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

            if (vflag_either) v_tally2_thr(this,atomk,atomn,-tmp2,rkn,thr);
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

      if (vflag_either) v_tally2_thr(this,atomj,atoml,-tmp2,rjl,thr);

      // due to kronecker(ltype, 0) term in contribution
      // to NconjtmpJ and later Nijconj
      if (ltype != 0) continue;

      tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*dwjl*SpN)/rjlmag;
      f[atomj][0] -= tmp2*rjl[0];
      f[atomj][1] -= tmp2*rjl[1];
      f[atomj][2] -= tmp2*rjl[2];
      f[atoml][0] += tmp2*rjl[0];
      f[atoml][1] += tmp2*rjl[1];
      f[atoml][2] += tmp2*rjl[2];

      if (vflag_either) v_tally2_thr(this,atomj,atoml,-tmp2,rjl,thr);

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

            if (vflag_either) v_tally2_thr(this,atoml,atomn,-tmp2,rln,thr);
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
        if ((sin321 > TOL) && (r21mag > TOL)) { // XXX was sin321 != 0.0
          sink2i = 1.0/(sin321*sin321);
          rik2i = 1.0/(r21mag*r21mag);
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
            if (atoml != atomi && atoml != atomk) {
              r34[0] = x[atom3][0]-x[atom4][0];
              r34[1] = x[atom3][1]-x[atom4][1];
              r34[2] = x[atom3][2]-x[atom4][2];
              r34mag = sqrt((r34[0]*r34[0])+(r34[1]*r34[1])+(r34[2]*r34[2]));
              cos234 = (r32[0]*r34[0] + r32[1]*r34[1] + r32[2]*r34[2]) /
                (r32mag*r34mag);
              cos234 = MIN(cos234,1.0);
              cos234 = MAX(cos234,-1.0);
              sin234 = sqrt(1.0 - cos234*cos234);

              if ((sin234 > TOL) && (r34mag > TOL)) {  // XXX was sin234 != 0.0
                sinl2i = 1.0/(sin234*sin234);
                rjl2i = 1.0/(r34mag*r34mag);
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
                Etmp += ((1.0-square(om1234))*w21*w34) *
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
                aaa2 = -prefactor*(1.0-square(om1234)) * w21*w34;
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

                if (vflag_either) {
                  r13[0] = -rjk[0]; r13[1] = -rjk[1]; r13[2] = -rjk[2];
                  r43[0] = -r34[0]; r43[1] = -r34[1]; r43[2] = -r34[2];
                  v_tally4_thr(this,atom1,atom2,atom3,atom4,f1,f2,f4,r13,r23,r43,thr);
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

        if (vflag_either) v_tally2_thr(this,atomi,atomk,-tmp2,rik,thr);

        // due to kronecker(ktype, 0) term in contribution
        // to NconjtmpI and later Nijconj
        if (ktype != 0) continue;

        tmp2 = VA*dN3[2]*(2.0*NconjtmpI*dwik*SpN)*Etmp/rikmag;
        f[atomi][0] -= tmp2*rik[0];
        f[atomi][1] -= tmp2*rik[1];
        f[atomi][2] -= tmp2*rik[2];
        f[atomk][0] += tmp2*rik[0];
        f[atomk][1] += tmp2*rik[1];
        f[atomk][2] += tmp2*rik[2];

        if (vflag_either) v_tally2_thr(this,atomi,atomk,-tmp2,rik,thr);

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

              if (vflag_either) v_tally2_thr(this,atomk,atomn,-tmp2,rkn,thr);
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

        if (vflag_either) v_tally2_thr(this,atomj,atoml,-tmp2,rjl,thr);

        // due to kronecker(ltype, 0) term in contribution
        // to NconjtmpJ and later Nijconj
        if (ltype != 0) continue;

        tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*dwjl*SpN)*Etmp/rjlmag;
        f[atomj][0] -= tmp2*rjl[0];
        f[atomj][1] -= tmp2*rjl[1];
        f[atomj][2] -= tmp2*rjl[2];
        f[atoml][0] += tmp2*rjl[0];
        f[atoml][1] += tmp2*rjl[1];
        f[atoml][2] += tmp2*rjl[2];

        if (vflag_either) v_tally2_thr(this,atomj,atoml,-tmp2,rjl,thr);

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

              if (vflag_either) v_tally2_thr(this,atoml,atomn,-tmp2,rln,thr);
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
-------------------------------------------------------------------------

This function calculates S(t_b(b_ij*)) as specified in the AIREBO paper.
To do so, it needs to compute b_ij*, i.e. the bondorder given that the
atoms i and j are placed a fictitious distance rijmag_mod apart.
Now there are two approaches to calculate the resulting forces:
1. Carry through the fictitious distance and corresponding vector
   rij_mod, correcting afterwards using the derivative of r/|r|.
2. Perform all the calculations using the real distance, and do not
   use a correction, only using rijmag_mod where necessary.
This code opts for (2). Mathematically, the approaches are equivalent
if implemented correctly, since in terms where only the normalized
vector is used, both calculations necessarily lead to the same result
since if f(x) = g(x/|x|) then for x = y/|y| f(x) = g(y/|y|/1).
The actual modified distance is only used in the lamda terms.
Note that these do not contribute to the forces between i and j, since
rijmag_mod is a constant and the corresponding derivatives are
accordingly zero.
This function should be kept in sync with bondorder(), i.e. changes
there probably also need to be performed here.

*/

double PairAIREBOOMP::bondorderLJ_thr(int i, int j, double /* rij_mod */[3], double rijmag_mod,
                                      double VA, double rij[3], double rijmag, ThrData * const thr)
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
  double fcikpc,fcjlpc,fcjkpc,fcilpc,fcijpc;
  double dt2dik[3],dt2djl[3],dt2dij[3],aa,aaa2,at2,cw,cwnum,cwnom;
  double sin321,sin234,rr,rijrik,rijrjl,rjk2,rik2,ril2,rjl2;
  double dctik,dctjk,dctjl,dctij,dctji,dctil,rik2i,rjl2i,sink2i,sinl2i;
  double rjk[3],ril[3],dt1dik,dt1djk,dt1djl,dt1dil,dt1dij;
  double F23[3],F12[3],F34[3],F31[3],F24[3],fi[3],fj[3],fk[3],fl[3];
  double f1[3],f2[3],f3[3],f4[4];
  double PijS,PjiS;
  double rij2,tspjik,dtsjik,tspijl,dtsijl,costmp;
  int *REBO_neighs,*REBO_neighs_i,*REBO_neighs_j,*REBO_neighs_k,*REBO_neighs_l;
  double tmppij,tmppji,dN2PIJ[2],dN2PJI[2],dN3piRC[3],dN3Tij[3];
  double tmp3pij,tmp3pji,Stb,dStb;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const int * const type = atom->type;

  atomi = i;
  atomj = j;
  itype = map[type[atomi]];
  jtype = map[type[atomj]];
  wij = Sp(rijmag,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
  NijC = nC[atomi]-(wij*kronecker(jtype,0));
  NijH = nH[atomi]-(wij*kronecker(jtype,1));
  NjiC = nC[atomj]-(wij*kronecker(itype,0));
  NjiH = nH[atomj]-(wij*kronecker(itype,1));
  bij = 0.0;
  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  dgdc = 0.0;
  dgdN = 0.0;
  NconjtmpI = 0.0;
  NconjtmpJ = 0.0;
  Etmp = 0.0;
  Stb = 0.0;
  dStb = 0.0;

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
        ((rho[ktype][1]-rikmag)-(rho[jtype][1]-rijmag_mod));
      wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dS);
      Nki = nC[atomk]-(wik*kronecker(itype,0))+
        nH[atomk]-(wik*kronecker(itype,1));
      cosjik = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2])) /
        (rijmag*rikmag);
      cosjik = MIN(cosjik,1.0);
      cosjik = MAX(cosjik,-1.0);

      // evaluate splines g and derivatives dg

      g = gSpline(cosjik,(NijC+NijH),itype,&dgdc,&dgdN);
      Etmp += (wik*g*exp(lamdajik));
      tmp3 += (wik*dgdN*exp(lamdajik));
      NconjtmpI = NconjtmpI+(kronecker(ktype,0)*wik*Sp(Nki,Nmin,Nmax,dS));
    }
  }

  PijS = 0.0;
  dN2PIJ[0] = 0.0;
  dN2PIJ[1] = 0.0;
  PijS = PijSpline(NijC,NijH,itype,jtype,dN2PIJ);
  pij = 1.0/sqrt(1.0+Etmp+PijS);
  tmppij = -.5*pij*pij*pij;
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
        ((rho[ltype][1]-rjlmag)-(rho[itype][1]-rijmag_mod));
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dS);
      Nlj = nC[atoml]-(wjl*kronecker(jtype,0))+nH[atoml] -
        (wjl*kronecker(jtype,1));
      cosijl = -1.0*((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2])) /
        (rijmag*rjlmag);
      cosijl = MIN(cosijl,1.0);
      cosijl = MAX(cosijl,-1.0);

      // evaluate splines g and derivatives dg

      g = gSpline(cosijl,NjiC+NjiH,jtype,&dgdc,&dgdN);
      Etmp += (wjl*g*exp(lamdaijl));
      tmp3 += (wjl*dgdN*exp(lamdaijl));
      NconjtmpJ = NconjtmpJ+(kronecker(ltype,0)*wjl*Sp(Nlj,Nmin,Nmax,dS));
    }
  }

  PjiS = 0.0;
  dN2PJI[0] = 0.0;
  dN2PJI[1] = 0.0;
  PjiS = PijSpline(NjiC,NjiH,jtype,itype,dN2PJI);
  pji = 1.0/sqrt(1.0+Etmp+PjiS);
  tmppji = -.5*pji*pji*pji;
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
        sin321 = sqrt(1.0 - cos321*cos321);
        if ((sin321 > TOL) && (r21mag > TOL)) { // XXX was sin321 != 0.0
          w21 = Sp(r21mag,rcmin[itype][ktype],rcmaxp[itype][ktype],dw21);
          tspjik = Sp2(cos321,thmin,thmax,dtsjik);

          REBO_neighs_j = REBO_firstneigh[j];
          for (l = 0; l < REBO_numneigh[j]; l++) {
            atoml = REBO_neighs_j[l];
            atom4 = atoml;
            ltype = map[type[atoml]];
            if (atoml != atomi && atoml != atomk) {
              r34[0] = x[atom3][0]-x[atom4][0];
              r34[1] = x[atom3][1]-x[atom4][1];
              r34[2] = x[atom3][2]-x[atom4][2];
              r34mag = sqrt((r34[0]*r34[0])+(r34[1]*r34[1])+(r34[2]*r34[2]));
              cos234 = (r32[0]*r34[0] + r32[1]*r34[1] + r32[2]*r34[2]) /
                (r32mag*r34mag);
              cos234 = MIN(cos234,1.0);
              cos234 = MAX(cos234,-1.0);
              sin234 = sqrt(1.0 - cos234*cos234);

              if ((sin234 > TOL) && (r34mag > TOL)) {  // XXX was sin234 != 0.0
                w34 = Sp(r34mag,rcmin[jtype][ltype],rcmaxp[jtype][ltype],dw34);
                tspijl = Sp2(cos234,thmin,thmax,dtsijl);

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
                Etmp += ((1.0-square(om1234))*w21*w34) *
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
      ktype = map[type[atomk]];
      if (atomk != atomj) {
        lamdajik = 0.0;
        rik[0] = x[atomi][0]-x[atomk][0];
        rik[1] = x[atomi][1]-x[atomk][1];
        rik[2] = x[atomi][2]-x[atomk][2];
        rikmag = sqrt(rik[0]*rik[0] + rik[1]*rik[1] + rik[2]*rik[2]);
        lamdajik = 4.0*kronecker(itype,1) *
          ((rho[ktype][1]-rikmag)-(rho[jtype][1]-rijmag_mod));
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
        fi[0] += tmp2*(rik[0]/rikmag);
        fi[1] += tmp2*(rik[1]/rikmag);
        fi[2] += tmp2*(rik[2]/rikmag);
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

        if (vflag_either) {
          rji[0] = -rij[0]; rji[1] = -rij[1]; rji[2] = -rij[2];
          rki[0] = -rik[0]; rki[1] = -rik[1]; rki[2] = -rik[2];
          v_tally3_thr(this,atomi,atomj,atomk,fj,fk,rji,rki,thr);
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
          ((rho[ltype][1]-rjlmag)-(rho[itype][1]-rijmag_mod));
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
          (cosijl*((rij[0]/square(rijmag))-(rjl[0]/(rjlmag*rjlmag))));
        dcosijldrj[1] = ((-rij[1]+rjl[1])/(rijmag*rjlmag)) +
          (cosijl*((rij[1]/square(rijmag))-(rjl[1]/(rjlmag*rjlmag))));
        dcosijldrj[2] = ((-rij[2]+rjl[2])/(rijmag*rjlmag)) +
          (cosijl*((rij[2]/square(rijmag))-(rjl[2]/(rjlmag*rjlmag))));
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
        fj[0] += tmp2*(rjl[0]/rjlmag);
        fj[1] += tmp2*(rjl[1]/rjlmag);
        fj[2] += tmp2*(rjl[2]/rjlmag);
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

        if (vflag_either) {
          rlj[0] = -rjl[0]; rlj[1] = -rjl[1]; rlj[2] = -rjl[2];
          v_tally3_thr(this,atomi,atomj,atoml,fi,fl,rij,rlj,thr);
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

        if (vflag_either) v_tally2_thr(this,atomi,atomk,-tmp2,rik,thr);

        // due to kronecker(ktype, 0) term in contribution
        // to NconjtmpI and later Nijconj
        if (ktype != 0) continue;

        tmp2 = VA*dN3[2]*(2.0*NconjtmpI*dwik*SpN)/rikmag;
        f[atomi][0] -= tmp2*rik[0];
        f[atomi][1] -= tmp2*rik[1];
        f[atomi][2] -= tmp2*rik[2];
        f[atomk][0] += tmp2*rik[0];
        f[atomk][1] += tmp2*rik[1];
        f[atomk][2] += tmp2*rik[2];

        if (vflag_either) v_tally2_thr(this,atomi,atomk,-tmp2,rik,thr);

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

              if (vflag_either) v_tally2_thr(this,atomk,atomn,-tmp2,rkn,thr);
            }
          }
        }
      }
    }

    // piRC forces to J side

    REBO_neighs = REBO_firstneigh[atomj];
    for (l = 0; l < REBO_numneigh[atomj]; l++) {
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

        if (vflag_either) v_tally2_thr(this,atomj,atoml,-tmp2,rjl,thr);

        // due to kronecker(ltype, 0) term in contribution
        // to NconjtmpJ and later Nijconj
        if (ltype != 0) continue;

        tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*dwjl*SpN)/rjlmag;
        f[atomj][0] -= tmp2*rjl[0];
        f[atomj][1] -= tmp2*rjl[1];
        f[atomj][2] -= tmp2*rjl[2];
        f[atoml][0] += tmp2*rjl[0];
        f[atoml][1] += tmp2*rjl[1];
        f[atoml][2] += tmp2*rjl[2];

        if (vflag_either) v_tally2_thr(this,atomj,atoml,-tmp2,rjl,thr);

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

              if (vflag_either) v_tally2_thr(this,atoml,atomn,-tmp2,rln,thr);
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
          if ((sin321 > TOL) && (r21mag > TOL)) { // XXX was sin321 != 0.0
            sink2i = 1.0/(sin321*sin321);
            rik2i = 1.0/(r21mag*r21mag);
            rr = (rijmag*rijmag)-(r21mag*r21mag);
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
              if (atoml != atomi && atoml != atomk) {
                r34[0] = x[atom3][0]-x[atom4][0];
                r34[1] = x[atom3][1]-x[atom4][1];
                r34[2] = x[atom3][2]-x[atom4][2];
                r34mag = sqrt(r34[0]*r34[0] + r34[1]*r34[1] + r34[2]*r34[2]);
                cos234 = (r32[0]*r34[0] + r32[1]*r34[1] + r32[2]*r34[2]) /
                  (r32mag*r34mag);
                cos234 = MIN(cos234,1.0);
                cos234 = MAX(cos234,-1.0);
                sin234 = sqrt(1.0 - cos234*cos234);

                if ((sin234 > TOL) && (r34mag > TOL)) { // XXX was sin234 != 0.0
                  sinl2i = 1.0/(sin234*sin234);
                  rjl2i = 1.0/(r34mag*r34mag);
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
                  aaa2 = -prefactor*(1.0-square(om1234)) * w21*w34;
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

                  if (vflag_either) {
                    r13[0] = -rjk[0]; r13[1] = -rjk[1]; r13[2] = -rjk[2];
                    r43[0] = -r34[0]; r43[1] = -r34[1]; r43[2] = -r34[2];
                    v_tally4_thr(this,atom1,atom2,atom3,atom4,f1,f2,f4,r13,r23,r43,thr);
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

          if (vflag_either) v_tally2_thr(this,atomi,atomk,-tmp2,rik,thr);

          // due to kronecker(ktype, 0) term in contribution
          // to NconjtmpI and later Nijconj
          if (ktype != 0) continue;

          tmp2 = VA*dN3[2]*(2.0*NconjtmpI*dwik*SpN)*Etmp/rikmag;
          f[atomi][0] -= tmp2*rik[0];
          f[atomi][1] -= tmp2*rik[1];
          f[atomi][2] -= tmp2*rik[2];
          f[atomk][0] += tmp2*rik[0];
          f[atomk][1] += tmp2*rik[1];
          f[atomk][2] += tmp2*rik[2];

          if (vflag_either) v_tally2_thr(this,atomi,atomk,-tmp2,rik,thr);

          if (fabs(dNki) > TOL) {
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

                if (vflag_either) v_tally2_thr(this,atomk,atomn,-tmp2,rkn,thr);
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

          if (vflag_either) v_tally2_thr(this,atomj,atoml,-tmp2,rjl,thr);

          // due to kronecker(ltype, 0) term in contribution
          // to NconjtmpJ and later Nijconj
          if (ltype != 0) continue;

          tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*dwjl*SpN)*Etmp/rjlmag;
          f[atomj][0] -= tmp2*rjl[0];
          f[atomj][1] -= tmp2*rjl[1];
          f[atomj][2] -= tmp2*rjl[2];
          f[atoml][0] += tmp2*rjl[0];
          f[atoml][1] += tmp2*rjl[1];
          f[atoml][2] += tmp2*rjl[2];

          if (vflag_either) v_tally2_thr(this,atomj,atoml,-tmp2,rjl,thr);

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

                if (vflag_either) v_tally2_thr(this,atoml,atomn,-tmp2,rln,thr);
              }
            }
          }
        }
      }
    }
  }

  return Stb;
}

/* ---------------------------------------------------------------------- */

double PairAIREBOOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairAIREBO::memory_usage();

  return bytes;
}
