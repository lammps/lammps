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

#include "pair_eam_omp.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairEAMOMP::PairEAMOMP(LAMMPS *lmp) :
  PairEAM(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairEAMOMP::~PairEAMOMP()
{
  for (auto *p : forcelist_thr) memory->destroy(p);
}

/* ---------------------------------------------------------------------- */

void PairEAMOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;
  int beyond_rhomax = 0;

  // grow energy and fp arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(fp);
    memory->destroy(numforce);
    nmax = atom->nmax;
    memory->create(rho,nthreads*nmax,"pair:rho");
    memory->create(fp,nmax,"pair:fp");
    memory->create(numforce,nmax,"pair:numforce");
  }

  if (nthreads > static_cast<int>(forcelist_thr.size())) {
    forcelist_thr.resize(nthreads, nullptr);
    maxforce_thr.resize(nthreads, 0);
  }

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag) reduction(+:beyond_rhomax)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    if (force->newton_pair)
      thr->init_eam(nall, rho);
    else
      thr->init_eam(atom->nlocal, rho);

    if (he_flag) {
      if (evflag) {
        if (eflag) {
          if (force->newton_pair) eval<1,1,1,1>(ifrom, ito, &beyond_rhomax, thr);
          else eval<1,1,0,1>(ifrom, ito, &beyond_rhomax, thr);
        } else {
          if (force->newton_pair) eval<1,0,1,1>(ifrom, ito, &beyond_rhomax, thr);
          else eval<1,0,0,1>(ifrom, ito, &beyond_rhomax, thr);
        }
      } else {
        if (force->newton_pair) eval<0,0,1,1>(ifrom, ito, &beyond_rhomax, thr);
        else eval<0,0,0,1>(ifrom, ito, &beyond_rhomax, thr);
      }
    } else {
      if (evflag) {
        if (eflag) {
          if (force->newton_pair) eval<1,1,1,0>(ifrom, ito, &beyond_rhomax, thr);
          else eval<1,1,0,0>(ifrom, ito, &beyond_rhomax, thr);
        } else {
          if (force->newton_pair) eval<1,0,1,0>(ifrom, ito, &beyond_rhomax, thr);
          else eval<1,0,0,0>(ifrom, ito, &beyond_rhomax, thr);
        }
      } else {
        if (force->newton_pair) eval<0,0,1,0>(ifrom, ito, &beyond_rhomax, thr);
        else eval<0,0,0,0>(ifrom, ito, &beyond_rhomax, thr);
      }
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region

  if (eflag && (!exceeded_rhomax)) {
    MPI_Allreduce(&beyond_rhomax, &exceeded_rhomax, 1, MPI_INT, MPI_SUM, world);
    if (exceeded_rhomax) {
      if (comm->me == 0)
        error->warning(FLERR,
                       "A per-atom density exceeded rhomax of EAM potential table - "
                       "a linear extrapolation to the energy was made");
    }
  }
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR, int HE>
void PairEAMOMP::eval(int iifrom, int iito, int *beyond_rhomax, ThrData * const thr)
{
  int i,j,ii,jj,m,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,p,rhoip,rhojp,z2,z2p,recip,phip,psip,phi;
  double *coeff;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  double * _noalias const rho_t = thr->get_rho();
  const int tid = thr->get_tid();
  const int nthreads = comm->nthreads;

  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  const double cutforcesq = this->cutforcesq;
  const double rdr = this->rdr;
  const int nr = this->nr;

  double fxtmp,fytmp,fztmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  const bigint init_cap = static_cast<bigint>(iito - iifrom) * 20 + 64;
  if (maxforce_thr[tid] < init_cap) {
    maxforce_thr[tid] = init_cap;
    memory->grow(forcelist_thr[tid],maxforce_thr[tid],"pair:forcelist_thr");
  }
  int * _noalias flist = forcelist_thr[tid];
  bigint maxf = maxforce_thr[tid];
  bigint nlist = 0;

  // rho = density at each atom
  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ii++) {
    i = ilist[ii];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (nlist + jnum > maxf) {
      const bigint needed = nlist + jnum;
      maxf = needed + (needed >> 2) + 64;
      maxforce_thr[tid] = maxf;
      memory->grow(forcelist_thr[tid],maxf,"pair:forcelist_thr");
      flist = forcelist_thr[tid];
    }

    double rhotmp = rho_t[i];
    int nforce = 0;
    const int * _noalias const type2rhor_i = type2rhor[itype];
    double * const * _noalias const rhor_spline_ii = rhor_spline[type2rhor_i[itype]];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
        flist[nlist++] = j;
        ++nforce;
        jtype = type[j];
        p = sqrt(rsq)*rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);
        if (jtype == itype) {
          coeff = rhor_spline_ii[m];
          const double rho_val = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
          rhotmp += rho_val;
          if (NEWTON_PAIR || j < nlocal) rho_t[j] += rho_val;
        } else {
          const int rhor_ji = type2rhor[jtype][itype];
          coeff = rhor_spline[rhor_ji][m];
          const double rho_val = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
          rhotmp += rho_val;
          if (NEWTON_PAIR || j < nlocal) {
            const int rhor_ij = type2rhor_i[jtype];
            if (rhor_ij == rhor_ji) rho_t[j] += rho_val;
            else {
              coeff = rhor_spline[rhor_ij][m];
              rho_t[j] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
            }
          }
        }
      }
    }
    numforce[i] = nforce;
    rho_t[i] = rhotmp;
  }

  // wait until all threads are done with computation
  sync_threads();

  // communicate and sum densities

  if (NEWTON_PAIR) {
    // reduce per thread density
    thr->timer(Timer::PAIR);
    data_reduce_thr(rho, nall, nthreads, 1, tid);

    // wait until reduction is complete
    sync_threads();

#if defined(_OPENMP)
#if _OPENMP >= 202011
#pragma omp masked
#else
#pragma omp master
#endif
#endif
    { comm->reverse_comm(this); }

    // wait until master thread is done with communication
    sync_threads();

  } else {
    thr->timer(Timer::PAIR);
    data_reduce_thr(rho, nlocal, nthreads, 1, tid);

    // wait until reduction is complete
    sync_threads();
  }

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // if rho > rhomax (e.g. due to close approach of two atoms) the table is
  //   exceeded, so add linear term to conserve energy; for eam/he the table
  //   starts at rhomin and may be exceeded on either side

  for (ii = iifrom; ii < iito; ii++) {
    i = ilist[ii];
    embedding_index<HE>(rho[i],m,p);
    coeff = frho_spline[type2frho[type[i]]][m];
    fp[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
    if (EFLAG) {
      phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
      if (HE && (rho[i] < rhomin)) {
        phi += fp[i] * (rho[i]-rhomin);
      } else if (rho[i] > rhomax) {
        phi += fp[i] * (rho[i]-rhomax);
        if (!HE) *beyond_rhomax = 1;
      }
      e_tally_thr(this, i, i, nlocal, NEWTON_PAIR, scale[type[i]][type[i]]*phi, 0.0, thr);
    }
  }

  // wait until all theads are done with computation
  sync_threads();

  // communicate derivative of embedding function
  // MPI communication only on master thread
#if defined(_OPENMP)
#if _OPENMP >= 202011
#pragma omp masked
#else
#pragma omp master
#endif
#endif
  { comm->forward_comm(this); }

  // wait until master thread is done with communication
  sync_threads();

  // compute forces on each atom
  // loop over neighbors of my atoms

  const double * _noalias const fp = this->fp;
  bigint cur = 0;

  for (ii = iifrom; ii < iito; ii++) {
    i = ilist[ii];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    fxtmp = fytmp = fztmp = 0.0;
    const int * _noalias const type2rhor_i = type2rhor[itype];
    const int * _noalias const type2z2r_i = type2z2r[itype];
    const double * _noalias const scale_i = scale[itype];
    double * const * _noalias const rhor_spline_ii = rhor_spline[type2rhor_i[itype]];
    double * const * _noalias const z2r_spline_ii = z2r_spline[type2z2r_i[itype]];
    const double fptmp = fp[i];
    const int nforce = numforce[i];
    const int * _noalias const jforce = flist + cur;
    cur += nforce;

    for (jj = 0; jj < nforce; jj++) {
      j = jforce[jj];

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;

      jtype = type[j];
      r = sqrt(rsq);
      p = r*rdr + 1.0;
      m = static_cast<int> (p);
      m = MIN(m,nr-1);
      p -= m;
      p = MIN(p,1.0);

      // rhoip = derivative of (density at atom j due to atom i)
      // rhojp = derivative of (density at atom i due to atom j)
      // phi = pair potential energy
      // phip = phi'
      // z2 = phi * r
      // z2p = (phi * r)' = (phi' r) + phi
      // psip needs both fp[i] and fp[j] terms since r_ij appears in two
      //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
      //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip

      if (jtype == itype) {
        coeff = rhor_spline_ii[m];
        rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
        rhojp = rhoip;
        coeff = z2r_spline_ii[m];
      } else {
        const int rhor_ij = type2rhor_i[jtype];
        const int rhor_ji = type2rhor[jtype][itype];
        coeff = rhor_spline[rhor_ij][m];
        rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
        if (rhor_ji == rhor_ij) rhojp = rhoip;
        else {
          coeff = rhor_spline[rhor_ji][m];
          rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
        }
        coeff = z2r_spline[type2z2r_i[jtype]][m];
      }
      z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
      z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

      recip = 1.0/r;
      phi = z2*recip;
      phip = z2p*recip - phi*recip;
      psip = fptmp*rhojp + fp[j]*rhoip + phip;
      fpair = -scale_i[jtype]*psip*recip;

      fxtmp += delx*fpair;
      fytmp += dely*fpair;
      fztmp += delz*fpair;
      if (NEWTON_PAIR || j < nlocal) {
        f[j].x -= delx*fpair;
        f[j].y -= dely*fpair;
        f[j].z -= delz*fpair;
      }

      if (EFLAG) evdwl = scale_i[jtype]*phi;
      if (EVFLAG) ev_tally_thr(this, i,j,nlocal,NEWTON_PAIR,
                               evdwl,0.0,fpair,delx,dely,delz,thr);
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairEAMOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairEAM::memory_usage();
  for (auto cap : maxforce_thr) bytes += (double) cap * sizeof(int);

  return bytes;
}
