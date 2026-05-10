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

#include "pair_lcbop_omp.h"

#include "atom.h"
#include "comm.h"
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

static constexpr double TOL = 1.0e-9;

/* ---------------------------------------------------------------------- */

PairLCBOPOMP::PairLCBOPOMP(LAMMPS *lmp) : PairLCBOP(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairLCBOPOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  SR_neigh_thr();

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

    FSR_thr(ifrom, ito, eflag, thr);
    FLR_thr(ifrom, ito, eflag, thr);

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

double PairLCBOPOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLCBOP::memory_usage();
  return bytes;
}

/* ----------------------------------------------------------------------
   create SR neighbor list from main neighbor list (parallel version)
   SR neighbor list stores neighbors of ghost atoms
------------------------------------------------------------------------- */

void PairLCBOPOMP::SR_neigh_thr()
{
  const int nthreads = comm->nthreads;

  if (atom->nmax > maxlocal) {
    maxlocal = atom->nmax;
    memory->destroy(SR_numneigh);
    memory->sfree(SR_firstneigh);
    memory->destroy(N);
    memory->destroy(M);
    memory->create(SR_numneigh,maxlocal,"LCBOP:numneigh");
    SR_firstneigh = (int **) memory->smalloc(maxlocal*sizeof(int *),
                           "LCBOP:firstneigh");
    memory->create(N,maxlocal,"LCBOP:N");
    memory->create(M,maxlocal,"LCBOP:M");
  }

  const int allnum = list->inum + list->gnum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // First pass: store all SR neighs of owned and ghost atoms,
  // compute N[i] = sum of cutoff functions with SR neighbors.
  // Each thread handles a contiguous block of atoms.

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(ilist,numneigh,firstneigh,allnum,nthreads)
#endif
  {
    int i,j,ii,jj,n,jnum;
    double xtmp,ytmp,ztmp,delx,dely,delz,rsq,dS;
    int *jlist,*neighptr;
    double **x = atom->x;

#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif

    const int iidelta = 1 + allnum/nthreads;
    const int iifrom = tid*iidelta;
    const int iito = ((iifrom+iidelta) > allnum) ? allnum : (iifrom+iidelta);

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
      N[i] = 0.0;
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        if (rsq < r_2_sq) {
          neighptr[n++] = j;
          N[i] += f_c(sqrt(rsq),r_1,r_2,&dS);
        }
      }

      SR_firstneigh[i] = neighptr;
      SR_numneigh[i] = n;
      ipg.vgot(n);
      if (ipg.status())
        error->one(FLERR, Error::NOLASTLINE,
                   "Neighbor list overflow, boost neigh_modify one" + utils::errorurl(36));
    }
  }
  // implicit OpenMP barrier after first parallel region

  // Second pass: compute M[i] = sum_j f_C_ij*F(N_j - f_C_ij).
  // Reads N[] (now complete) and writes M[i] (no overlap between threads).

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(ilist,allnum,nthreads)
#endif
  {
    int i,j,ii,jj,jnum;
    double xtmp,ytmp,ztmp,delx,dely,delz,rsq,dS;
    int *jlist;
    double **x = atom->x;

#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif

    const int iidelta = 1 + allnum/nthreads;
    const int iifrom = tid*iidelta;
    const int iito = ((iifrom+iidelta) > allnum) ? allnum : (iifrom+iidelta);

    for (ii = iifrom; ii < iito; ii++) {
      i = ilist[ii];

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      M[i] = 0.0;

      jlist = SR_firstneigh[i];
      jnum = SR_numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        if (rsq < r_2_sq) {
          double f_c_ij = f_c(sqrt(rsq),r_1,r_2,&dS);
          double Nji = N[j]-f_c_ij;
          // F(xij) = 1-f_c_LR(Nji, 2,3,&dummy)
          M[i] += f_c_ij * ( 1-f_c_LR(Nji, 2,3,&dS) );
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLCBOPOMP::FNij_thr(int i, int j, double factor, double **f, ThrData *const thr)
{
  int atomi = i;
  int atomj = j;
  int *SR_neighs = SR_firstneigh[i];
  double **x = atom->x;
  for (int k = 0; k < SR_numneigh[i]; k++) {
    int atomk = SR_neighs[k];
    if (atomk != atomj) {
      double rik[3];
      rik[0] = x[atomi][0]-x[atomk][0];
      rik[1] = x[atomi][1]-x[atomk][1];
      rik[2] = x[atomi][2]-x[atomk][2];
      double riksq = (rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]);
      if (riksq > r_1*r_1) {
        double rikmag = sqrt(riksq);
        double df_c_ik;
        f_c( rikmag, r_1, r_2, &df_c_ik );

        double fpair = -factor*df_c_ik / rikmag;
        f[atomi][0] += rik[0]*fpair;
        f[atomi][1] += rik[1]*fpair;
        f[atomi][2] += rik[2]*fpair;
        f[atomk][0] -= rik[0]*fpair;
        f[atomk][1] -= rik[1]*fpair;
        f[atomk][2] -= rik[2]*fpair;

        if (vflag_either) v_tally2_thr(this,atomi,atomk,fpair,rik,thr);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLCBOPOMP::FMij_thr(int i, int j, double factor, double **f, ThrData *const thr)
{
  int atomi = i;
  int atomj = j;
  int *SR_neighs = SR_firstneigh[i];
  double **x = atom->x;
  for (int k = 0; k < SR_numneigh[i]; k++) {
    int atomk = SR_neighs[k];
    if (atomk != atomj) {
      double rik[3];
      rik[0] = x[atomi][0]-x[atomk][0];
      rik[1] = x[atomi][1]-x[atomk][1];
      rik[2] = x[atomi][2]-x[atomk][2];
      double rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      double df_c_ik;
      double f_c_ik = f_c( rikmag, r_1, r_2, &df_c_ik );
      double Nki = N[k]-f_c_ik;
      double dF=0;
      double Fx = 1-f_c_LR(Nki, 2,3,&dF);
      dF = -dF;

      if (df_c_ik > TOL) {
        double factor2 = factor*df_c_ik*Fx;
        double fpair = -factor2 / rikmag;
        f[atomi][0] += rik[0]*fpair;
        f[atomi][1] += rik[1]*fpair;
        f[atomi][2] += rik[2]*fpair;
        f[atomk][0] -= rik[0]*fpair;
        f[atomk][1] -= rik[1]*fpair;
        f[atomk][2] -= rik[2]*fpair;
        if (vflag_either) v_tally2_thr(this,atomi,atomk,fpair,rik,thr);
      }

      if (dF > TOL) {
        double factor2 = factor*f_c_ik*dF;
        FNij_thr(atomk, atomi, factor2, f, thr);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairLCBOPOMP::b_thr(int i, int j, double rij[3], double rijmag,
                            double VA, double **f, ThrData *const thr)
{
  int *SR_neighs = SR_firstneigh[i];
  double **x = atom->x;
  int atomi = i;
  int atomj = j;

  // calculate bij magnitude
  double bij = 1.0;
  for (int k = 0; k < SR_numneigh[i]; k++) {
    int atomk = SR_neighs[k];
    if (atomk != atomj) {
      double rik[3];
      rik[0] = x[atomi][0]-x[atomk][0];
      rik[1] = x[atomi][1]-x[atomk][1];
      rik[2] = x[atomi][2]-x[atomk][2];
      double rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      double delta_ijk = rijmag-rikmag;
      double dummy;
      double f_c_ik = f_c( rikmag, r_1, r_2, &dummy );
      double cos_ijk = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2]))
                / (rijmag*rikmag);
      cos_ijk = MIN(cos_ijk,1.0);
      cos_ijk = MAX(cos_ijk,-1.0);

      double G = gSpline(cos_ijk,   &dummy);
      double H = hSpline(delta_ijk, &dummy);
      bij += (f_c_ik*G*H);
    }
  }
  bij = pow( bij, -delta );

  // bij forces

  for (int k = 0; k < SR_numneigh[i]; k++) {
    int atomk = SR_neighs[k];
    if (atomk != atomj) {
      double rik[3];
      rik[0] = x[atomi][0]-x[atomk][0];
      rik[1] = x[atomi][1]-x[atomk][1];
      rik[2] = x[atomi][2]-x[atomk][2];
      double rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      double delta_ijk = rijmag-rikmag;
      double df_c_ik;
      double f_c_ik = f_c( rikmag, r_1, r_2, &df_c_ik );
      double cos_ijk = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2]))
                / (rijmag*rikmag);
      cos_ijk = MIN(cos_ijk,1.0);
      cos_ijk = MAX(cos_ijk,-1.0);

      double dcos_ijk_dri[3],dcos_ijk_drj[3],dcos_ijk_drk[3];
      dcos_ijk_drj[0] = -rik[0] / (rijmag*rikmag)
             + cos_ijk * rij[0] / (rijmag*rijmag);
      dcos_ijk_drj[1] = -rik[1] / (rijmag*rikmag)
             + cos_ijk * rij[1] / (rijmag*rijmag);
      dcos_ijk_drj[2] = -rik[2] / (rijmag*rikmag)
             + cos_ijk * rij[2] / (rijmag*rijmag);

      dcos_ijk_drk[0] = -rij[0] / (rijmag*rikmag)
             + cos_ijk * rik[0] / (rikmag*rikmag);
      dcos_ijk_drk[1] = -rij[1] / (rijmag*rikmag)
             + cos_ijk * rik[1] / (rikmag*rikmag);
      dcos_ijk_drk[2] = -rij[2] / (rijmag*rikmag)
             + cos_ijk * rik[2] / (rikmag*rikmag);

      dcos_ijk_dri[0] = -dcos_ijk_drk[0] - dcos_ijk_drj[0];
      dcos_ijk_dri[1] = -dcos_ijk_drk[1] - dcos_ijk_drj[1];
      dcos_ijk_dri[2] = -dcos_ijk_drk[2] - dcos_ijk_drj[2];

      double dG, dH;
      double G = gSpline( cos_ijk,   &dG );
      double H = hSpline( delta_ijk, &dH );
      double tmp = -VA*0.5*(-0.5*bij*bij*bij);

      double fi[3], fj[3], fk[3];

      double tmp2 = -tmp*df_c_ik*G*H/rikmag;
      fi[0] =  tmp2*rik[0];
      fi[1] =  tmp2*rik[1];
      fi[2] =  tmp2*rik[2];
      fk[0] = -tmp2*rik[0];
      fk[1] = -tmp2*rik[1];
      fk[2] = -tmp2*rik[2];

      tmp2 = -tmp*f_c_ik*dG*H;
      fi[0] += tmp2*dcos_ijk_dri[0];
      fi[1] += tmp2*dcos_ijk_dri[1];
      fi[2] += tmp2*dcos_ijk_dri[2];
      fj[0] =  tmp2*dcos_ijk_drj[0];
      fj[1] =  tmp2*dcos_ijk_drj[1];
      fj[2] =  tmp2*dcos_ijk_drj[2];
      fk[0] += tmp2*dcos_ijk_drk[0];
      fk[1] += tmp2*dcos_ijk_drk[1];
      fk[2] += tmp2*dcos_ijk_drk[2];

      tmp2 = -tmp*f_c_ik*G*dH;
      fi[0] += tmp2*( rij[0]/rijmag - rik[0]/rikmag );
      fi[1] += tmp2*( rij[1]/rijmag - rik[1]/rikmag );
      fi[2] += tmp2*( rij[2]/rijmag - rik[2]/rikmag );
      fj[0] += tmp2*( -rij[0]/rijmag );
      fj[1] += tmp2*( -rij[1]/rijmag );
      fj[2] += tmp2*( -rij[2]/rijmag );
      fk[0] += tmp2*( rik[0]/rikmag );
      fk[1] += tmp2*( rik[1]/rikmag );
      fk[2] += tmp2*( rik[2]/rikmag );

      f[atomi][0] += fi[0]; f[atomi][1] += fi[1]; f[atomi][2] += fi[2];
      f[atomj][0] += fj[0]; f[atomj][1] += fj[1]; f[atomj][2] += fj[2];
      f[atomk][0] += fk[0]; f[atomk][1] += fk[1]; f[atomk][2] += fk[2];

      if (vflag_either) {
        double rji[3], rki[3];
        rji[0] = -rij[0]; rji[1] = -rij[1]; rji[2] = -rij[2];
        rki[0] = -rik[0]; rki[1] = -rik[1]; rki[2] = -rik[2];
        v_tally3_thr(this,atomi,atomj,atomk,fj,fk,rji,rki,thr);
      }
    }
  }

  return bij;
}

/* ---------------------------------------------------------------------- */

double PairLCBOPOMP::bondorder_thr(int i, int j, double rij[3], double rijmag,
                                   double VA, double **f, ThrData *const thr)
{
  double bij, bji;
  {
    double rji[3];
    rji[0] = -rij[0]; rji[1] = -rij[1]; rji[2] = -rij[2];
    bij = b_thr(i,j,rij,rijmag,VA,f,thr);
    bji = b_thr(j,i,rji,rijmag,VA,f,thr);
  }

  double Fij_conj;
  {
    double dummy;

    double df_c_ij;
    double f_c_ij = f_c(rijmag, r_1, r_2, &df_c_ij);
    double Nij = MIN(3,N[i]-f_c_ij);
    double Nji = MIN(3,N[j]-f_c_ij);

    double Mij = M[i] - f_c_ij*(1-f_c(Nji,2,3,&dummy));
    double Mji = M[j] - f_c_ij*(1-f_c(Nij,2,3,&dummy));
    Mij = MIN(Mij,3);
    Mji = MIN(Mji,3);

    double Nij_el, dNij_el_dNij, dNij_el_dMij;
    double Nji_el, dNji_el_dNji, dNji_el_dMji;
    {
      double num_Nij_el = 4 - Mij;
      double num_Nji_el = 4 - Mji;
      double den_Nij_el = Nij + 1 - Mij;
      double den_Nji_el = Nji + 1 - Mji;
      Nij_el = num_Nij_el / den_Nij_el;
      Nji_el = num_Nji_el / den_Nji_el;
      dNij_el_dNij = -Nij_el/den_Nij_el;
      dNji_el_dNji = -Nji_el/den_Nji_el;
      dNij_el_dMij = (-1 + Nij_el) /den_Nij_el;
      dNji_el_dMji = (-1 + Nji_el) /den_Nji_el;
    }

    double Nconj;
    double dNconj_dNij;
    double dNconj_dNji;
    double dNconj_dNel;
    {
      double num_Nconj = ( Nij+1 )*( Nji+1 )*( Nij_el+Nji_el ) - 4*( Nij+Nji+2);
      double den_Nconj = Nij*( 3-Nij )*( Nji+1 ) + Nji*( 3-Nji )*( Nij+1 ) + eps;
      Nconj = num_Nconj / den_Nconj;
      if (Nconj <= 0) {
        Nconj = 0;
        dNconj_dNij = 0;
        dNconj_dNji = 0;
        dNconj_dNel = 0;
      } else if (Nconj >= 1) {
        Nconj = 1;
        dNconj_dNij = 0;
        dNconj_dNji = 0;
        dNconj_dNel = 0;
      } else {
        dNconj_dNij = (
            ( (Nji+1)*(Nij_el + Nji_el)-4)
            - Nconj*( (Nji+1)*(3-2*Nij) + Nji*(3-Nji) )
          ) /den_Nconj;
        dNconj_dNji = (
            ( (Nij+1)*(Nji_el + Nij_el)-4)
            - Nconj*( (Nij+1)*(3-2*Nji) + Nij*(3-Nij) )
          ) /den_Nconj;
        dNconj_dNel = (Nij+1)*(Nji+1) / den_Nconj;
      }
    }

    double dF_dNij, dF_dNji, dF_dNconj;
    Fij_conj = F_conj(Nij, Nji, Nconj, &dF_dNij, &dF_dNji, &dF_dNconj);

    if (3-Nij > TOL) {
      double factor = -VA*0.5*(dF_dNij + dF_dNconj*(dNconj_dNij + dNconj_dNel*dNij_el_dNij));
      FNij_thr(i, j, factor, f, thr);
    }
    if (3-Nji > TOL) {
      double factor = -VA*0.5*(dF_dNji + dF_dNconj*(dNconj_dNji + dNconj_dNel*dNji_el_dNji));
      FNij_thr(j, i, factor, f, thr);
    }
    if (3-Mij > TOL) {
      double factor = -VA*0.5*(dF_dNconj*dNconj_dNel*dNij_el_dMij);
      FMij_thr(i, j, factor, f, thr);
    }
    if (3-Mji > TOL) {
      double factor = -VA*0.5*(dF_dNconj*dNconj_dNel*dNji_el_dMji);
      FMij_thr(j, i, factor, f, thr);
    }
  }

  return 0.5*(bij + bji + Fij_conj);
}

/* ---------------------------------------------------------------------- */

void PairLCBOPOMP::FSR_thr(int ifrom, int ito, int eflag, ThrData *const thr)
{
  int i,j,jj,ii;
  tagint itag,jtag;
  double delx,dely,delz,fpair,xtmp,ytmp,ztmp;
  double r_sq,rijmag,f_c_ij,df_c_ij;
  double VR,dVRdi,VA,Bij,dVAdi,dVA;
  double del[3];
  int *ilist,*SR_neighs;

  double **x = atom->x;
  double **f = thr->get_f();
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  ilist = list->ilist;

  for (ii = ifrom; ii < ito; ii++) {
    i = ilist[ii];
    itag = tag[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    SR_neighs = SR_firstneigh[i];

    for (jj = 0; jj < SR_numneigh[i]; jj++) {
      j = SR_neighs[jj];
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

      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];
      r_sq = delx*delx + dely*dely + delz*delz;
      rijmag = sqrt(r_sq);
      f_c_ij = f_c( rijmag,r_1,r_2,&df_c_ij );
      if (f_c_ij <= TOL) continue;

      VR = A*exp(-alpha*rijmag);
      dVRdi = -alpha*VR;
      dVRdi = dVRdi*f_c_ij + df_c_ij*VR;
      VR *= f_c_ij;

      VA = dVA = 0.0;
      {
        double term = B_1 * exp(-beta_1*rijmag);
        VA += term;
        dVA += -beta_1 * term;
        term = B_2 * exp(-beta_2*rijmag);
        VA += term;
        dVA += -beta_2 * term;
      }
      dVA = dVA*f_c_ij + df_c_ij*VA;
      VA *= f_c_ij;
      del[0] = delx;
      del[1] = dely;
      del[2] = delz;
      Bij = bondorder_thr(i,j,del,rijmag,VA,f,thr);
      dVAdi = Bij*dVA;

      fpair = -(dVRdi-dVAdi) / rijmag;
      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      double evdwl=0.0;
      if (eflag) evdwl = VR - Bij*VA;
      if (evflag) ev_tally_thr(this,i,j,nlocal,/* newton_pair */ 1,
                               evdwl,0.0,fpair,delx,dely,delz,thr);
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLCBOPOMP::FLR_thr(int ifrom, int ito, int eflag, ThrData *const thr)
{
  int i,j,jj,ii;
  tagint itag,jtag;
  double delx,dely,delz,fpair,xtmp,ytmp,ztmp;
  double r_sq,rijmag,f_c_ij,df_c_ij;
  double V,dVdi;

  double **x = atom->x;
  double **f = thr->get_f();
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  for (ii = ifrom; ii < ito; ii++) {
    i = ilist[ii];
    itag = tag[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    int *neighs = firstneigh[i];

    for (jj = 0; jj < numneigh[i]; jj++) {
      j = neighs[jj];
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

      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];
      r_sq = delx*delx + dely*dely + delz*delz;
      rijmag = sqrt(r_sq);
      f_c_ij = 1-f_c( rijmag,r_1,r_2,&df_c_ij );
      df_c_ij = -df_c_ij;
      f_c_ij *= f_c_LR( rijmag, r_1_LR, r_2_LR, &df_c_ij );
      if (f_c_ij <= TOL) continue;

      V = dVdi = 0;
      if (rijmag<r_0) {
        double exp_part = exp( -lambda_1*(rijmag-r_0) );
        V = eps_1*( exp_part*exp_part - 2*exp_part) + v_1;
        dVdi = 2*eps_1*lambda_1*exp_part*( 1-exp_part );
      } else {
        double exp_part = exp( -lambda_2*(rijmag-r_0) );
        V = eps_2*( exp_part*exp_part - 2*exp_part) + v_2;
        dVdi = 2*eps_2*lambda_2*exp_part*( 1-exp_part );
      }
      dVdi = dVdi*f_c_ij + df_c_ij*V;
      V *= f_c_ij;

      fpair = -dVdi / rijmag;
      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      double evdwl=0.0;
      if (eflag) evdwl = V;
      if (evflag) ev_tally_thr(this,i,j,nlocal,/* newton_pair */ 1,
                               evdwl,0.0,fpair,delx,dely,delz,thr);
    }
  }
}
