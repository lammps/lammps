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
   Contributing authors: Roy Pollock (LLNL), Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "ewald_omp.h"
#include <mpi.h>
#include <cmath>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "timer.h"
#include "math_const.h"

#include "suffix.h"
using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.00001

/* ---------------------------------------------------------------------- */

EwaldOMP::EwaldOMP(LAMMPS *lmp) : Ewald(lmp), ThrOMP(lmp, THR_KSPACE)
{
  triclinic_support = 0;
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */
void EwaldOMP::allocate()
{
  Ewald::allocate();

  // always re-allocate for simplicity.
  delete[] sfacrl;
  delete[] sfacim;

  sfacrl = new double[kmax3d*comm->nthreads];
  sfacim = new double[kmax3d*comm->nthreads];
}

/* ----------------------------------------------------------------------
   compute the Ewald long-range force, energy, virial
------------------------------------------------------------------------- */

void EwaldOMP::compute(int eflag, int vflag)
{
  // set energy/virial flags

  ev_init(eflag,vflag);

  // extend size of per-atom arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(ek);
    memory->destroy3d_offset(cs,-kmax_created);
    memory->destroy3d_offset(sn,-kmax_created);
    nmax = atom->nmax;
    memory->create(ek,nmax,3,"ewald:ek");
    memory->create3d_offset(cs,-kmax,kmax,3,nmax,"ewald:cs");
    memory->create3d_offset(sn,-kmax,kmax,3,nmax,"ewald:sn");
    kmax_created = kmax;
  }

  // partial structure factors on each processor
  // total structure factor by summing over procs

  eik_dot_r();
  MPI_Allreduce(sfacrl,sfacrl_all,kcount,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(sfacim,sfacim_all,kcount,MPI_DOUBLE,MPI_SUM,world);

  // update qsum and qsqsum, if atom count has changed and energy needed
  // (n.b. needs to be done outside of the multi-threaded region)

  if ((eflag_global || eflag_atom) && atom->natoms != natoms_original) {
    qsum_qsq();
    natoms_original = atom->natoms;
  }

  // K-space portion of electric field
  // double loop over K-vectors and local atoms

  double * const * const f = atom->f;
  const double * const q = atom->q;
  const int nthreads = comm->nthreads;
  const int nlocal = atom->nlocal;
  const double qscale = force->qqrd2e * scale;

  double eng_tmp = 0.0;
  double v0,v1,v2,v3,v4,v5;
  v0=v1=v2=v3=v4=v5=0.0;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag) reduction(+:eng_tmp,v0,v1,v2,v3,v4,v5)
#endif
  {

    int i,j,k,ifrom,ito,tid;
    int kx,ky,kz;
    double cypz,sypz,exprl,expim,partial;

    loop_setup_thr(ifrom, ito, tid, nlocal, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, 0, NULL, NULL, thr);

    for (i = ifrom; i < ito; i++) {
      ek[i][0] = 0.0;
      ek[i][1] = 0.0;
      ek[i][2] = 0.0;
    }

    for (k = 0; k < kcount; k++) {
      kx = kxvecs[k];
      ky = kyvecs[k];
      kz = kzvecs[k];

      for (i = ifrom; i < ito; i++) {
        cypz = cs[ky][1][i]*cs[kz][2][i] - sn[ky][1][i]*sn[kz][2][i];
        sypz = sn[ky][1][i]*cs[kz][2][i] + cs[ky][1][i]*sn[kz][2][i];
        exprl = cs[kx][0][i]*cypz - sn[kx][0][i]*sypz;
        expim = sn[kx][0][i]*cypz + cs[kx][0][i]*sypz;
        partial = expim*sfacrl_all[k] - exprl*sfacim_all[k];
        ek[i][0] += partial*eg[k][0];
        ek[i][1] += partial*eg[k][1];
        ek[i][2] += partial*eg[k][2];

        if (evflag_atom) {
          const double partial_peratom = exprl*sfacrl_all[k] + expim*sfacim_all[k];
          if (eflag_atom) eatom[i] += q[i]*ug[k]*partial_peratom;
          if (vflag_atom)
            for (j = 0; j < 6; j++)
              vatom[i][j] += ug[k]*vg[k][j]*partial_peratom;
        }
      }
    }

    // convert E-field to force

    for (i = ifrom; i < ito; i++) {
      const double fac = qscale*q[i];
      f[i][0] += fac*ek[i][0];
      f[i][1] += fac*ek[i][1];
      f[i][2] += fac*ek[i][2];
    }

    // global energy

    if (eflag_global) {
#if defined(_OPENMP)
#pragma omp for private(k)
#endif
      for (k = 0; k < kcount; k++)
        eng_tmp += ug[k] * (sfacrl_all[k]*sfacrl_all[k] +
                            sfacim_all[k]*sfacim_all[k]);
    }

    // global virial

    if (vflag_global) {
#if defined(_OPENMP)
#pragma omp for private(k)
#endif
      for (k = 0; k < kcount; k++) {
        double uk = ug[k] * (sfacrl_all[k]*sfacrl_all[k] + sfacim_all[k]*sfacim_all[k]);
        v0 += uk*vg[k][0];
        v1 += uk*vg[k][1];
        v2 += uk*vg[k][2];
        v3 += uk*vg[k][3];
        v4 += uk*vg[k][4];
        v5 += uk*vg[k][5];
      }
    }

    // per-atom energy/virial
    // energy includes self-energy correction

    if (evflag_atom) {
      if (eflag_atom) {
        for (i = ifrom; i < ito; i++) {
          eatom[i] -= g_ewald*q[i]*q[i]/MY_PIS + MY_PI2*q[i]*qsum /
            (g_ewald*g_ewald*volume);
          eatom[i] *= qscale;
        }
      }

      if (vflag_atom)
        for (i = ifrom; i < ito; i++)
          for (j = 0; j < 6; j++) vatom[i][j] *= q[i]*qscale;
    }

    thr->timer(Timer::KSPACE);
    reduce_thr(this, eflag,vflag,thr);
  } // end of omp parallel region

  if (eflag_global) {
    eng_tmp -= g_ewald*qsqsum/MY_PIS +
      MY_PI2*qsum*qsum / (g_ewald*g_ewald*volume);
    energy = eng_tmp * qscale;
  }

  if (vflag_global) {
    virial[0] = v0 * qscale;
    virial[1] = v1 * qscale;
    virial[2] = v2 * qscale;
    virial[3] = v3 * qscale;
    virial[4] = v4 * qscale;
    virial[5] = v5 * qscale;
  }

  if (slabflag) slabcorr();
}

/* ---------------------------------------------------------------------- */

void EwaldOMP::eik_dot_r()
{
  const double * const * const x = atom->x;
  const double * const q = atom->q;
  const int nlocal = atom->nlocal;
  const int nthreads = comm->nthreads;

#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {
    int i,ifrom,ito,k,l,m,n,ic,tid;
    double cstr1,sstr1,cstr2,sstr2,cstr3,sstr3,cstr4,sstr4;
    double sqk,clpm,slpm;

    loop_setup_thr(ifrom, ito, tid, nlocal, nthreads);

    double * const sfacrl_thr = sfacrl + tid*kmax3d;
    double * const sfacim_thr = sfacim + tid*kmax3d;

    n = 0;

    // (k,0,0), (0,l,0), (0,0,m)

    for (ic = 0; ic < 3; ic++) {
      sqk = unitk[ic]*unitk[ic];
      if (sqk <= gsqmx) {
        cstr1 = 0.0;
        sstr1 = 0.0;
        for (i = ifrom; i < ito; i++) {
          cs[0][ic][i] = 1.0;
          sn[0][ic][i] = 0.0;
          cs[1][ic][i] = cos(unitk[ic]*x[i][ic]);
          sn[1][ic][i] = sin(unitk[ic]*x[i][ic]);
          cs[-1][ic][i] = cs[1][ic][i];
          sn[-1][ic][i] = -sn[1][ic][i];
          cstr1 += q[i]*cs[1][ic][i];
          sstr1 += q[i]*sn[1][ic][i];
        }
        sfacrl_thr[n] = cstr1;
        sfacim_thr[n++] = sstr1;
      }
    }

    for (m = 2; m <= kmax; m++) {
      for (ic = 0; ic < 3; ic++) {
        sqk = m*unitk[ic] * m*unitk[ic];
        if (sqk <= gsqmx) {
          cstr1 = 0.0;
          sstr1 = 0.0;
          for (i = ifrom; i < ito; i++) {
            cs[m][ic][i] = cs[m-1][ic][i]*cs[1][ic][i] -
              sn[m-1][ic][i]*sn[1][ic][i];
            sn[m][ic][i] = sn[m-1][ic][i]*cs[1][ic][i] +
              cs[m-1][ic][i]*sn[1][ic][i];
            cs[-m][ic][i] = cs[m][ic][i];
            sn[-m][ic][i] = -sn[m][ic][i];
            cstr1 += q[i]*cs[m][ic][i];
            sstr1 += q[i]*sn[m][ic][i];
          }
          sfacrl_thr[n] = cstr1;
          sfacim_thr[n++] = sstr1;
        }
      }
    }

    // 1 = (k,l,0), 2 = (k,-l,0)

    for (k = 1; k <= kxmax; k++) {
      for (l = 1; l <= kymax; l++) {
        sqk = (k*unitk[0] * k*unitk[0]) + (l*unitk[1] * l*unitk[1]);
        if (sqk <= gsqmx) {
          cstr1 = 0.0;
          sstr1 = 0.0;
          cstr2 = 0.0;
          sstr2 = 0.0;
          for (i = ifrom; i < ito; i++) {
            cstr1 += q[i]*(cs[k][0][i]*cs[l][1][i] - sn[k][0][i]*sn[l][1][i]);
            sstr1 += q[i]*(sn[k][0][i]*cs[l][1][i] + cs[k][0][i]*sn[l][1][i]);
            cstr2 += q[i]*(cs[k][0][i]*cs[l][1][i] + sn[k][0][i]*sn[l][1][i]);
            sstr2 += q[i]*(sn[k][0][i]*cs[l][1][i] - cs[k][0][i]*sn[l][1][i]);
          }
          sfacrl_thr[n] = cstr1;
          sfacim_thr[n++] = sstr1;
          sfacrl_thr[n] = cstr2;
          sfacim_thr[n++] = sstr2;
        }
      }
    }

    // 1 = (0,l,m), 2 = (0,l,-m)

    for (l = 1; l <= kymax; l++) {
      for (m = 1; m <= kzmax; m++) {
        sqk = (l*unitk[1] * l*unitk[1]) + (m*unitk[2] * m*unitk[2]);
        if (sqk <= gsqmx) {
          cstr1 = 0.0;
          sstr1 = 0.0;
          cstr2 = 0.0;
          sstr2 = 0.0;
          for (i = ifrom; i < ito; i++) {
            cstr1 += q[i]*(cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i]);
            sstr1 += q[i]*(sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i]);
            cstr2 += q[i]*(cs[l][1][i]*cs[m][2][i] + sn[l][1][i]*sn[m][2][i]);
            sstr2 += q[i]*(sn[l][1][i]*cs[m][2][i] - cs[l][1][i]*sn[m][2][i]);
          }
          sfacrl_thr[n] = cstr1;
          sfacim_thr[n++] = sstr1;
          sfacrl_thr[n] = cstr2;
          sfacim_thr[n++] = sstr2;
        }
      }
    }

    // 1 = (k,0,m), 2 = (k,0,-m)

    for (k = 1; k <= kxmax; k++) {
      for (m = 1; m <= kzmax; m++) {
        sqk = (k*unitk[0] * k*unitk[0]) + (m*unitk[2] * m*unitk[2]);
        if (sqk <= gsqmx) {
          cstr1 = 0.0;
          sstr1 = 0.0;
          cstr2 = 0.0;
          sstr2 = 0.0;
          for (i = ifrom; i < ito; i++) {
            cstr1 += q[i]*(cs[k][0][i]*cs[m][2][i] - sn[k][0][i]*sn[m][2][i]);
            sstr1 += q[i]*(sn[k][0][i]*cs[m][2][i] + cs[k][0][i]*sn[m][2][i]);
            cstr2 += q[i]*(cs[k][0][i]*cs[m][2][i] + sn[k][0][i]*sn[m][2][i]);
            sstr2 += q[i]*(sn[k][0][i]*cs[m][2][i] - cs[k][0][i]*sn[m][2][i]);
          }
          sfacrl_thr[n] = cstr1;
          sfacim_thr[n++] = sstr1;
          sfacrl_thr[n] = cstr2;
          sfacim_thr[n++] = sstr2;
        }
      }
    }

    // 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)

    for (k = 1; k <= kxmax; k++) {
      for (l = 1; l <= kymax; l++) {
        for (m = 1; m <= kzmax; m++) {
          sqk = (k*unitk[0] * k*unitk[0]) + (l*unitk[1] * l*unitk[1]) +
            (m*unitk[2] * m*unitk[2]);
          if (sqk <= gsqmx) {
            cstr1 = 0.0;
            sstr1 = 0.0;
            cstr2 = 0.0;
            sstr2 = 0.0;
            cstr3 = 0.0;
            sstr3 = 0.0;
            cstr4 = 0.0;
            sstr4 = 0.0;
            for (i = ifrom; i < ito; i++) {
              clpm = cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i];
              slpm = sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i];
              cstr1 += q[i]*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
              sstr1 += q[i]*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);

              clpm = cs[l][1][i]*cs[m][2][i] + sn[l][1][i]*sn[m][2][i];
              slpm = -sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i];
              cstr2 += q[i]*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
              sstr2 += q[i]*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);

              clpm = cs[l][1][i]*cs[m][2][i] + sn[l][1][i]*sn[m][2][i];
              slpm = sn[l][1][i]*cs[m][2][i] - cs[l][1][i]*sn[m][2][i];
              cstr3 += q[i]*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
              sstr3 += q[i]*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);

              clpm = cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i];
              slpm = -sn[l][1][i]*cs[m][2][i] - cs[l][1][i]*sn[m][2][i];
              cstr4 += q[i]*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
              sstr4 += q[i]*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);
            }
            sfacrl_thr[n] = cstr1;
            sfacim_thr[n++] = sstr1;
            sfacrl_thr[n] = cstr2;
            sfacim_thr[n++] = sstr2;
            sfacrl_thr[n] = cstr3;
            sfacim_thr[n++] = sstr3;
            sfacrl_thr[n] = cstr4;
            sfacim_thr[n++] = sstr4;
          }
        }
      }
    }

    sync_threads();
    data_reduce_thr(sfacrl,kmax3d,nthreads,1,tid);
    data_reduce_thr(sfacim,kmax3d,nthreads,1,tid);

  } // end of parallel region
}
