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

#include <cmath>
#include <cstring>

#include "pair_eam_cd_omp.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"

#include "suffix.h"
using namespace LAMMPS_NS;

// This is for debugging purposes. The ASSERT() macro is used in the code to check
// if everything runs as expected. Change this to #if 0 if you don't need the checking.
#if 0
        #define ASSERT(cond) ((!(cond)) ? my_failure(error,__FILE__,__LINE__) : my_noop())

        inline void my_noop() {}
        inline void my_failure(Error* error, const char* file, int line) {
                char str[1024];
                sprintf(str,"Assertion failure: File %s, line %i", file, line);
                error->one(FLERR,str);
        }
#else
        #define ASSERT(cond)
#endif

/* ---------------------------------------------------------------------- */

PairEAMCDOMP::PairEAMCDOMP(LAMMPS *lmp, int _cdeamVersion) :
  PairEAM(lmp), PairEAMCD(lmp,_cdeamVersion), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairEAMCDOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  // grow energy and fp arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(rhoB);
    memory->destroy(D_values);
    memory->destroy(fp);
    nmax = atom->nmax;
    memory->create(rho,nthreads*nmax,"pair:rho");
    memory->create(rhoB,nthreads*nmax,"pair:mu");
    memory->create(D_values,nthreads*nmax,"pair:D_values");
    memory->create(fp,nmax,"pair:fp");
  }

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    if (force->newton_pair)
      thr->init_cdeam(nall, rho, rhoB, D_values);
    else
      thr->init_cdeam(atom->nlocal, rho, rhoB, D_values);

    switch (cdeamVersion) {

    case 1:

      if (evflag) {
        if (eflag) {
          if (force->newton_pair) eval<1,1,1,1>(ifrom, ito, thr);
          else eval<1,1,0,1>(ifrom, ito, thr);
        } else {
          if (force->newton_pair) eval<1,0,1,1>(ifrom, ito, thr);
          else eval<1,0,0,1>(ifrom, ito, thr);
        }
      } else {
        if (force->newton_pair) eval<0,0,1,1>(ifrom, ito, thr);
        else eval<0,0,0,1>(ifrom, ito, thr);
      }
      break;

    case 2:

      if (evflag) {
        if (eflag) {
          if (force->newton_pair) eval<1,1,1,2>(ifrom, ito, thr);
          else eval<1,1,0,2>(ifrom, ito, thr);
        } else {
          if (force->newton_pair) eval<1,0,1,2>(ifrom, ito, thr);
          else eval<1,0,0,2>(ifrom, ito, thr);
        }
      } else {
        if (force->newton_pair) eval<0,0,1,2>(ifrom, ito, thr);
        else eval<0,0,0,2>(ifrom, ito, thr);
      }
      break;

    default:
      {
#if defined(_OPENMP)
#pragma omp master
#endif
        error->all(FLERR,"unsupported eam/cd pair style variant");
      }
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR, int CDEAMVERSION>
void PairEAMCDOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rhoip,rhojp,recip,phi;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
  double * const rho_t = thr->get_rho();
  double * const rhoB_t = thr->get_rhoB();
  double * const D_values_t = thr->get_D_values();
  const int tid = thr->get_tid();
  const int nthreads = comm->nthreads;

  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  double fxtmp,fytmp,fztmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Stage I

  // Compute rho and rhoB at each local atom site.
  // Additionally calculate the D_i values here if we are using the one-site formulation.
  // For the two-site formulation we have to calculate the D values in an extra loop (Stage II).

  for (ii = iifrom; ii < iito; ii++) {
    i = ilist[ii];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;

      if(rsq < cutforcesq) {
        jtype = type[j];
        double r = sqrt(rsq);
        const EAMTableIndex index = radiusToTableIndex(r);
        double localrho = RhoOfR(index, jtype, itype);
        rho_t[i] += localrho;
        if(jtype == speciesB) rhoB_t[i] += localrho;
        if(NEWTON_PAIR || j < nlocal) {
          localrho = RhoOfR(index, itype, jtype);
          rho_t[j] += localrho;
          if(itype == speciesB) rhoB_t[j] += localrho;
        }

        if(CDEAMVERSION == 1 && itype != jtype) {
          // Note: if the i-j interaction is not concentration dependent (because either
          // i or j are not species A or B) then its contribution to D_i and D_j should
          // be ignored.
          // This if-clause is only required for a ternary.
          if((itype == speciesA && jtype == speciesB)
             || (jtype == speciesA && itype == speciesB)) {
            double Phi_AB = PhiOfR(index, itype, jtype, 1.0 / r);
            D_values_t[i] += Phi_AB;
            if(NEWTON_PAIR || j < nlocal)
              D_values_t[j] += Phi_AB;
          }
        }
      }
    }
  }

  // wait until all threads are done with computation
  sync_threads();

  // communicate and sum densities

  if (NEWTON_PAIR) {
    // reduce per thread density
    thr->timer(Timer::PAIR);
    data_reduce_thr(rho, nall, nthreads, 1, tid);
    data_reduce_thr(rhoB, nall, nthreads, 1, tid);
    if (CDEAMVERSION==1)
      data_reduce_thr(D_values, nall, nthreads, 1, tid);

    // wait until reduction is complete
    sync_threads();

#if defined(_OPENMP)
#pragma omp master
#endif
    { communicationStage = 1;
      comm->reverse_comm_pair(this); }

    // wait until master thread is done with communication
    sync_threads();

  } else {
    // reduce per thread density
    thr->timer(Timer::PAIR);
    data_reduce_thr(rho, nlocal, nthreads, 1, tid);
    data_reduce_thr(rhoB, nlocal, nthreads, 1, tid);
    if (CDEAMVERSION==1)
      data_reduce_thr(D_values, nlocal, nthreads, 1, tid);

    // wait until reduction is complete
    sync_threads();
  }

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom

  for (ii = iifrom; ii < iito; ii++) {
    i = ilist[ii];
    EAMTableIndex index = rhoToTableIndex(rho[i]);
    fp[i] = FPrimeOfRho(index, type[i]);
    if(EFLAG) {
      phi = FofRho(index, type[i]);
      e_tally_thr(this, i, i, nlocal, NEWTON_PAIR, phi, 0.0, thr);
    }
  }

  // wait until all theads are done with computation
  sync_threads();

  // Communicate derivative of embedding function and densities
  // and D_values (this for one-site formulation only).
#if defined(_OPENMP)
#pragma omp master
#endif
  { communicationStage = 2;
    comm->forward_comm_pair(this); }

  // wait until master thread is done with communication
  sync_threads();


  // The electron densities may not drop to zero because then the concentration would no longer be defined.
  // But the concentration is not needed anyway if there is no interaction with another atom, which is the case
  // if the electron density is exactly zero. That's why the following lines have been commented out.
  //
  //for(i = 0; i < nlocal + atom->nghost; i++) {
  //        if(rho[i] == 0 && (type[i] == speciesA || type[i] == speciesB))
  //                error->one(FLERR,"CD-EAM potential routine: Detected atom with zero electron density.");
  //}

  // Stage II
  // This is only required for the original two-site formulation of the CD-EAM potential.

  if(CDEAMVERSION == 2) {
    // Compute intermediate value D_i for each atom.
    for (ii = iifrom; ii < iito; ii++) {
      i = ilist[ii];
      xtmp = x[i].x;
      ytmp = x[i].y;
      ztmp = x[i].z;
      itype = type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // This code line is required for ternary alloys.
      if(itype != speciesA && itype != speciesB) continue;

      double x_i = rhoB[i] / rho[i];        // Concentration at atom i.

      for(jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        jtype = type[j];
        if(itype == jtype) continue;

        // This code line is required for ternary alloys.
        if(jtype != speciesA && jtype != speciesB) continue;

        delx = xtmp - x[j].x;
        dely = ytmp - x[j].y;
        delz = ztmp - x[j].z;
        rsq = delx*delx + dely*dely + delz*delz;

        if(rsq < cutforcesq) {
          double r = sqrt(rsq);
          const EAMTableIndex index = radiusToTableIndex(r);

          // The concentration independent part of the cross pair potential.
          double Phi_AB = PhiOfR(index, itype, jtype, 1.0 / r);

          // Average concentration of two sites
          double x_ij = 0.5 * (x_i + rhoB[j]/rho[j]);

          // Calculate derivative of h(x_ij) polynomial function.
          double h_prime = evalHprime(x_ij);

          D_values_t[i] += h_prime * Phi_AB / (2.0 * rho[i] * rho[i]);
          if(NEWTON_PAIR || j < nlocal)
            D_values_t[j] += h_prime * Phi_AB / (2.0 * rho[j] * rho[j]);
        }
      }
    }

    if (NEWTON_PAIR) {
    thr->timer(Timer::PAIR);
      data_reduce_thr(D_values, nall, nthreads, 1, tid);

      // wait until reduction is complete
      sync_threads();

#if defined(_OPENMP)
#pragma omp master
#endif
      { communicationStage = 3;
        comm->reverse_comm_pair(this); }

      // wait until master thread is done with communication
      sync_threads();

  } else {
    thr->timer(Timer::PAIR);
      data_reduce_thr(D_values, nlocal, nthreads, 1, tid);

    // wait until reduction is complete
    sync_threads();
  }

#if defined(_OPENMP)
#pragma omp master
#endif
    { communicationStage = 4;
      comm->forward_comm_pair(this); }

    // wait until master thread is done with communication
    sync_threads();
  }

  // Stage III

  // Compute force acting on each atom.
  for (ii = iifrom; ii < iito; ii++) {
    i = ilist[ii];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];
    fxtmp = fytmp = fztmp = 0.0;

    jlist = firstneigh[i];
    jnum = numneigh[i];

    // Concentration at site i
    double x_i = -1.0;                // The value -1 indicates: no concentration dependence for all interactions of atom i.
    // It will be replaced by the concentration at site i if atom i is either A or B.

    double D_i, h_prime_i;

    // This if-clause is only required for ternary alloys.
    if((itype == speciesA || itype == speciesB) && rho[i] != 0.0) {

      // Compute local concentration at site i.
      x_i = rhoB[i]/rho[i];
      ASSERT(x_i >= 0 && x_i<=1.0);

      if(CDEAMVERSION == 1) {
        // Calculate derivative of h(x_i) polynomial function.
        h_prime_i = evalHprime(x_i);
        D_i = D_values[i] * h_prime_i / (2.0 * rho[i] * rho[i]);
      } else if(CDEAMVERSION == 2) {
        D_i = D_values[i];
      } else {
        ASSERT(false);
      }
    }

    for(jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;

      if(rsq < cutforcesq) {
        jtype = type[j];
        double r = sqrt(rsq);
        const EAMTableIndex index = radiusToTableIndex(r);

        // rhoip = derivative of (density at atom j due to atom i)
        // rhojp = derivative of (density at atom i due to atom j)
        // psip needs both fp[i] and fp[j] terms since r_ij appears in two
        //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
        //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip
        rhoip = RhoPrimeOfR(index, itype, jtype);
        rhojp = RhoPrimeOfR(index, jtype, itype);
        fpair = fp[i]*rhojp + fp[j]*rhoip;
        recip = 1.0/r;

        double x_j = -1;  // The value -1 indicates: no concentration dependence for this i-j pair
        // because atom j is not of species A nor B.

        // This code line is required for ternary alloy.
        if(jtype == speciesA || jtype == speciesB) {
          ASSERT(rho[i] != 0.0);
          ASSERT(rho[j] != 0.0);

          // Compute local concentration at site j.
          x_j = rhoB[j]/rho[j];
          ASSERT(x_j >= 0 && x_j<=1.0);

          double D_j;
          if(CDEAMVERSION == 1) {
            // Calculate derivative of h(x_j) polynomial function.
            double h_prime_j = evalHprime(x_j);
            D_j = D_values[j] * h_prime_j / (2.0 * rho[j] * rho[j]);
          } else if(CDEAMVERSION == 2) {
            D_j = D_values[j];
          } else {
            ASSERT(false);
          }
          double t2 = -rhoB[j];
          if(itype == speciesB) t2 += rho[j];
          fpair += D_j * rhoip * t2;
        }

        // This if-clause is only required for a ternary alloy.
        // Actually we don't need it at all because D_i should be zero anyway if
        // atom i has no concentration dependent interactions (because it is not species A or B).
        if(x_i != -1.0) {
          double t1 = -rhoB[i];
          if(jtype == speciesB) t1 += rho[i];
          fpair += D_i * rhojp * t1;
        }

        double phip;
        double phi = PhiOfR(index, itype, jtype, recip, phip);
        if(itype == jtype || x_i == -1.0 || x_j == -1.0) {
          // Case of no concentration dependence.
          fpair += phip;
        } else {
          // We have a concentration dependence for the i-j interaction.
          double h;
          if(CDEAMVERSION == 1) {
            // Calculate h(x_i) polynomial function.
            double h_i = evalH(x_i);
            // Calculate h(x_j) polynomial function.
            double h_j = evalH(x_j);
            h = 0.5 * (h_i + h_j);
          } else if(CDEAMVERSION == 2) {
            // Average concentration.
            double x_ij = 0.5 * (x_i + x_j);
            // Calculate h(x_ij) polynomial function.
            h = evalH(x_ij);
          } else {
            ASSERT(false);
          }
          fpair += h * phip;
          phi *= h;
        }

        // Divide by r_ij and negate to get forces from gradient.
        fpair /= -r;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if(NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx*fpair;
          f[j].y -= dely*fpair;
          f[j].z -= delz*fpair;
        }

        if(EFLAG) evdwl = phi;
        if(EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,evdwl,0.0,
                                fpair,delx,dely,delz,thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairEAMCDOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairEAMCD::memory_usage();

  return bytes;
}
