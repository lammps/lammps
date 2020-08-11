/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "omp_compat.h"
#include <cmath>
#include "pair_comb_omp.h"
#include "atom.h"
#include "comm.h"
#include "group.h"
#include "force.h"
#include "memory.h"
#include "my_page.h"
#include "neighbor.h"
#include "neigh_list.h"

#include "suffix.h"
using namespace LAMMPS_NS;

#define MAXNEIGH 24

/* ---------------------------------------------------------------------- */

PairCombOMP::PairCombOMP(LAMMPS *lmp) :
  PairComb(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairCombOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  // Build short range neighbor list

  Short_neigh_thr();

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, NULL, thr);

    if (evflag) {
      if (eflag) {
        if (vflag_atom) eval<1,1,1>(ifrom, ito, thr);
        else eval<1,1,0>(ifrom, ito, thr);
      } else {
        if (vflag_atom) eval<1,0,1>(ifrom, ito, thr);
        else eval<1,0,0>(ifrom, ito, thr);
      }
    } else eval<0,0,0>(ifrom, ito, thr);

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int VFLAG_ATOM>
void PairCombOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,k,ii,jj,kk,jnum,iparam_i;
  tagint itag,jtag;
  int itype,jtype,ktype,iparam_ij,iparam_ijk;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fi[3],fj[3],fk[3];
  double zeta_ij,prefactor;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int mr1,mr2,mr3;
  int rsc,inty;
  double elp_ij,filp[3],fjlp[3],fklp[3];
  double iq,jq;
  double yaself;
  double potal,fac11,fac11e;
  double vionij,fvionij,sr1,sr2,sr3,Eov,Fov;
  int sht_jnum, *sht_jlist, nj;

  evdwl = ecoul = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const double * const q = atom->q;
  const tagint * const tag = atom->tag;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  yaself = vionij = fvionij = Eov = Fov = 0.0;

  double fxtmp,fytmp,fztmp;
  double fjxtmp,fjytmp,fjztmp;

  // self energy correction term: potal

  potal_calc(potal,fac11,fac11e);

  // loop over full neighbor list of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp = fytmp = fztmp = 0.0;

    iq = q[i];
    NCo[i] = 0;
    nj = 0;
    iparam_i = elem2param[itype][itype][itype];

    // self energy, only on i atom

    yaself = self(&params[iparam_i],iq,potal);

    if (EVFLAG) ev_tally_thr(this,i,i,nlocal,0,yaself,
                             0.0,0.0,0.0,0.0,0.0,thr);

    // two-body interactions (long and short repulsive)

    jlist = firstneigh[i];
    jnum = numneigh[i];
    sht_jlist = sht_first[i];
    sht_jnum = sht_num[i];

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

      // Qj calculates 2-body Coulombic

      jtype = map[type[j]];
      jq = q[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      iparam_ij = elem2param[itype][jtype][jtype];

      // long range q-dependent

      if (rsq > params[iparam_ij].lcutsq) continue;

      inty = intype[itype][jtype];

      // polynomial three-point interpolation

      tri_point(rsq, mr1, mr2, mr3, sr1, sr2, sr3, itype);

      // 1/r energy and forces

      direct(inty,mr1,mr2,mr3,rsq,sr1,sr2,sr3,iq,jq,
             potal,fac11,fac11e,vionij,fvionij);

      // field correction to self energy

      field(&params[iparam_ij],rsq,iq,jq,vionij,fvionij);

      // polarization field
      // sums up long range forces

      fxtmp += delx*fvionij;
      fytmp += dely*fvionij;
      fztmp += delz*fvionij;
      f[j][0] -= delx*fvionij;
      f[j][1] -= dely*fvionij;
      f[j][2] -= delz*fvionij;

      if (EVFLAG)
        ev_tally_thr(this,i,j,nlocal,/* newton_pair */ 1,
                     0.0,vionij,fvionij,delx,dely,delz,thr);

      // short range q-independent

      if (rsq > params[iparam_ij].cutsq) continue;

      repulsive(&params[iparam_ij],rsq,fpair,EFLAG,evdwl,iq,jq);

      // repulsion is pure two-body, sums up pair repulsive forces

      fxtmp += delx*fpair;
      fytmp += dely*fpair;
      fztmp += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (EVFLAG)
        ev_tally_thr(this,i,j,nlocal,/* newton_pair */ 1,
                     evdwl,0.0,fpair,delx,dely,delz,thr);
    }

    // accumulate coordination number information

    if (cor_flag) {
      int numcoor = 0;
      for (jj = 0; jj < sht_jnum; jj++) {
        j = sht_jlist[jj];
        jtype = map[type[j]];
        iparam_ij = elem2param[itype][jtype][jtype];

        if(params[iparam_ij].hfocor > 0.0 ) {
          delr1[0] = x[j][0] - xtmp;
          delr1[1] = x[j][1] - ytmp;
          delr1[2] = x[j][2] - ztmp;
          rsq1 = vec3_dot(delr1,delr1);

          if (rsq1 > params[iparam_ij].cutsq) continue;
          ++numcoor;
        }
      }
      NCo[i] = numcoor;
    }

    // three-body interactions
    // half i-j loop

    for (jj = 0; jj < sht_jnum; jj++) {
      j = sht_jlist[jj];

      jtype = map[type[j]];
      iparam_ij = elem2param[itype][jtype][jtype];

      // this Qj for q-dependent BSi

      jq = q[j];

      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = vec3_dot(delr1,delr1);

      if (rsq1 > params[iparam_ij].cutsq) continue;
      nj ++;

      // accumulate bondorder zeta for each i-j interaction via loop over k

      fjxtmp = fjytmp = fjztmp = 0.0;
      zeta_ij = 0.0;
      cuo_flag1 = 0; cuo_flag2 = 0;

      for (kk = 0; kk < sht_jnum; kk++) {
        k = sht_jlist[kk];
        if (j == k) continue;
        ktype = map[type[k]];
        iparam_ijk = elem2param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = vec3_dot(delr2,delr2);

        if (rsq2 > params[iparam_ijk].cutsq) continue;

        zeta_ij += zeta(&params[iparam_ijk],rsq1,rsq2,delr1,delr2);

        if (params[iparam_ijk].hfocor == -2.0) cuo_flag1 = 1;
        if (params[iparam_ijk].hfocor == -1.0) cuo_flag2 = 1;
      }

      if (cuo_flag1 && cuo_flag2) cuo_flag = 1;
      else cuo_flag = 0;

      force_zeta(&params[iparam_ij],EFLAG,i,nj,rsq1,zeta_ij,
                 iq,jq,fpair,prefactor,evdwl);

      // over-coordination correction for HfO2

      if (cor_flag && NCo[i] != 0)
        Over_cor(&params[iparam_ij],rsq1,NCo[i],Eov, Fov);
      evdwl +=  Eov;
      fpair +=  Fov;

      fxtmp += delr1[0]*fpair;
      fytmp += delr1[1]*fpair;
      fztmp += delr1[2]*fpair;
      fjxtmp -= delr1[0]*fpair;
      fjytmp -= delr1[1]*fpair;
      fjztmp -= delr1[2]*fpair;

      if (EVFLAG) ev_tally_thr(this,i,j,nlocal,/* newton_pair */ 1,evdwl,0.0,
                               -fpair,-delr1[0],-delr1[1],-delr1[2],thr);

      // attractive term via loop over k (3-body forces)

      for (kk = 0; kk < sht_jnum; kk++) {
        k = sht_jlist[kk];
        if (j == k) continue;
        ktype = map[type[k]];
        iparam_ijk = elem2param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = vec3_dot(delr2,delr2);
        if (rsq2 > params[iparam_ijk].cutsq) continue;

        for (rsc = 0; rsc < 3; rsc++)
          fi[rsc] = fj[rsc] = fk[rsc] = 0.0;

        attractive(&params[iparam_ijk],prefactor,
                   rsq1,rsq2,delr1,delr2,fi,fj,fk);

        // 3-body LP and BB correction and forces

        elp_ij = elp(&params[iparam_ijk],rsq1,rsq2,delr1,delr2);
        flp(&params[iparam_ijk],rsq1,rsq2,delr1,delr2,filp,fjlp,fklp);

        fxtmp += fi[0] + filp[0];
        fytmp += fi[1] + filp[1];
        fztmp += fi[2] + filp[2];
        fjxtmp += fj[0] + fjlp[0];
        fjytmp += fj[1] + fjlp[1];
        fjztmp += fj[2] + fjlp[2];
        f[k][0] += fk[0] + fklp[0];
        f[k][1] += fk[1] + fklp[1];
        f[k][2] += fk[2] + fklp[2];

        if (EVFLAG)
          ev_tally_thr(this,i,j,nlocal,/* newton_pair */ 1,
                       elp_ij,0.0,0.0,0.0,0.0,0.0, thr);
        if (VFLAG_ATOM) v_tally3_thr(i,j,k,fj,fk,delr1,delr2,thr);
      }
      f[j][0] += fjxtmp;
      f[j][1] += fjytmp;
      f[j][2] += fjztmp;
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;

    if (cuo_flag) params[iparam_i].cutsq *= 0.65;
  }
  cuo_flag = 0;
}

/* ---------------------------------------------------------------------- */

double PairCombOMP::yasu_char(double *qf_fix, int &igroup)
{
  int ii;
  double potal,fac11,fac11e;

  const double * const * const x = atom->x;
  const double * const q = atom->q;
  const int * const type = atom->type;
  const tagint * const tag = atom->tag;

  const int inum = list->inum;
  const int * const ilist = list->ilist;
  const int * const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;

  const int * const mask = atom->mask;
  const int groupbit = group->bitmask[igroup];

  qf = qf_fix;
  for (ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (mask[i] & groupbit)
      qf[i] = 0.0;
  }

  // communicating charge force to all nodes, first forward then reverse

  comm->forward_comm_pair(this);

  // self energy correction term: potal

  potal_calc(potal,fac11,fac11e);

  // loop over full neighbor list of my atoms
#if defined(_OPENMP)
#pragma omp parallel for private(ii) LMP_DEFAULT_NONE LMP_SHARED(potal,fac11e)
#endif
  for (ii = 0; ii < inum; ii ++) {
    double fqi,fqj,fqij,fqji,fqjj,delr1[3];
    double sr1,sr2,sr3;
    int mr1,mr2,mr3;

    const int i = ilist[ii];
    const tagint itag = tag[i];
    int nj = 0;

    if (mask[i] & groupbit) {
      fqi = fqj = fqij = fqji = fqjj = 0.0; // should not be needed.
      int itype = map[type[i]];
      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      const double iq = q[i];
      const int iparam_i = elem2param[itype][itype][itype];

      // charge force from self energy

      fqi = qfo_self(&params[iparam_i],iq,potal);

      // two-body interactions

      const int * const jlist = firstneigh[i];
      const int jnum = numneigh[i];

      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;
        const tagint jtag = tag[j];

        if (itag > jtag) {
          if ((itag+jtag) % 2 == 0) continue;
        } else if (itag < jtag) {
          if ((itag+jtag) % 2 == 1) continue;
        } else {
          if (x[j][2] < ytmp) continue;
          if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
          if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
        }

        const int jtype = map[type[j]];
        double jq = q[j];

        delr1[0] = x[j][0] - xtmp;
        delr1[1] = x[j][1] - ytmp;
        delr1[2] = x[j][2] - ztmp;
        double rsq1 = vec3_dot(delr1,delr1);

        const int iparam_ij = elem2param[itype][jtype][jtype];

        // long range q-dependent

        if (rsq1 > params[iparam_ij].lcutsq) continue;

        const int inty = intype[itype][jtype];

        // polynomial three-point interpolation

        tri_point(rsq1,mr1,mr2,mr3,sr1,sr2,sr3,itype);

        // 1/r charge forces

        qfo_direct(inty,mr1,mr2,mr3,rsq1,sr1,sr2,sr3,fac11e,fqij);

        // field correction to self energy and charge force

        qfo_field(&params[iparam_ij],rsq1,iq,jq,fqji,fqjj);
        fqi   += jq * fqij + fqji;
#if defined(_OPENMP) && !defined(__NVCC__)
#pragma omp atomic
#endif
        qf[j] += (iq * fqij + fqjj);
      }

        // three-body interactions

      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;
        const int jtype = map[type[j]];
        const double jq = q[j];

        delr1[0] = x[j][0] - xtmp;
        delr1[1] = x[j][1] - ytmp;
        delr1[2] = x[j][2] - ztmp;
        double rsq1 = vec3_dot(delr1,delr1);

        const int iparam_ij = elem2param[itype][jtype][jtype];

        if (rsq1 > params[iparam_ij].cutsq) continue;
        nj ++;

        // charge force in Aij and Bij

        qfo_short(&params[iparam_ij],i,nj,rsq1,iq,jq,fqij,fqjj);
        fqi += fqij;
#if defined(_OPENMP) && !defined(__NVCC__)
#pragma omp atomic
#endif
        qf[j] += fqjj;
      }

#if defined(_OPENMP) && !defined(__NVCC__)
#pragma omp atomic
#endif
      qf[i] += fqi;
    }
  }

  comm->reverse_comm_pair(this);

  // sum charge force on each node and return it

  double eneg = 0.0;
  for (ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (mask[i] & groupbit)
      eneg += qf[i];
  }
  double enegtot;
  MPI_Allreduce(&eneg,&enegtot,1,MPI_DOUBLE,MPI_SUM,world);
  return enegtot;
}

/* ---------------------------------------------------------------------- */

double PairCombOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairComb::memory_usage();

  return bytes;
}
/* ---------------------------------------------------------------------- */

void PairCombOMP::Short_neigh_thr()
{

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->sfree(sht_first);
    sht_first = (int **) memory->smalloc(nmax*sizeof(int *),
                                         "pair:sht_first");
    memory->grow(sht_num,nmax,"pair:sht_num");
    memory->grow(NCo,nmax,"pair:NCo");
    memory->grow(bbij,nmax,MAXNEIGH,"pair:bbij");
  }

  const int nthreads = comm->nthreads;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE
#endif
  {
    int nj,*neighptrj;
    int *ilist,*jlist,*numneigh,**firstneigh;
    int jnum,i,j,ii,jj;
    double xtmp,ytmp,ztmp,rsq,delrj[3];
    double **x = atom->x;

    const int inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;


#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif

    const int iidelta = 1 + inum/nthreads;
    const int iifrom = tid*iidelta;
    const int iito   = ((iifrom + iidelta) > inum) ? inum : (iifrom+iidelta);

    // each thread has its own page allocator
    MyPage<int> &ipg = ipage[tid];
    ipg.reset();

    // create Comb neighbor list

    for (ii = iifrom; ii < iito; ii++) {
      i = ilist[ii];

      nj = 0;
      neighptrj = ipg.vget();

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delrj[0] = xtmp - x[j][0];
        delrj[1] = ytmp - x[j][1];
        delrj[2] = ztmp - x[j][2];
        rsq = vec3_dot(delrj,delrj);

        if (rsq > cutmin) continue;
        neighptrj[nj++] = j;
      }
      sht_first[i] = neighptrj;
      sht_num[i] = nj;
      ipg.vgot(nj);
      if (ipg.status())
        error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    }
  }
}
