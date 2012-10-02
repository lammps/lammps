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

#include "math.h"
#include "pair_lj_cut_coul_long_proxy_tip4p_omp.h"
#include "pppm_tip4p_proxy.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"

#include <string.h>

#include "suffix.h"
using namespace LAMMPS_NS;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairLJCutCoulLongProxyTIP4POMP::PairLJCutCoulLongProxyTIP4POMP(LAMMPS *lmp) :
  PairLJCutCoulLongTIP4P(lmp), ThrOMP(lmp, THR_PAIR|THR_PROXY)
{
  proxyflag = 1;
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  nproxy=1;

  kspace = NULL;

  // TIP4P cannot compute virial as F dot r
  // due to finding bonded H atoms which are not near O atom

  no_virial_fdotr_compute = 1;
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulLongProxyTIP4POMP::init_style()
{
  if (comm->nthreads < 2)
    error->all(FLERR,"need at least two threads per MPI task for this pair style");

  kspace = static_cast<PPPMTIP4PProxy *>(force->kspace);

  PairLJCutCoulLongTIP4P::init_style();
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulLongProxyTIP4POMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  // reallocate hneigh & newsite if necessary
  // initialize hneigh[0] to -1 on steps when reneighboring occurred
  // initialize hneigh[2] to 0 every step

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(hneigh);
    memory->create(hneigh,nmax,3,"pair:hneigh");
    memory->destroy(newsite);
    memory->create(newsite,nmax,3,"pair:newsite");
  }

  // XXX: this could be threaded, too.
  int i;
  if (neighbor->ago == 0)
    for (i = 0; i < nall; i++) hneigh[i][0] = -1;
  for (i = 0; i < nall; i++) hneigh[i][2] = 0;

  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads, nproxy);
    ThrData *thr = fix->get_thr(tid);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, thr);

    // thread id 0 runs pppm, the rest the pair style
    if (tid < nproxy) {
      kspace->compute_proxy(eflag,vflag);
    } else {
      if (!ncoultablebits) {
        if (evflag) {
          if (eflag) {
            if (vflag) eval<1,1,1,1>(ifrom, ito, thr);
            else eval<1,1,1,0>(ifrom, ito, thr);
          } else {
            if (vflag) eval<1,1,0,1>(ifrom, ito, thr);
            else eval<1,1,0,0>(ifrom, ito, thr);
          }
        } else eval<1,0,0,0>(ifrom, ito, thr);
      } else {
        if (evflag) {
          if (eflag) {
            if (vflag) eval<0,1,1,1>(ifrom, ito, thr);
            else eval<0,1,1,0>(ifrom, ito, thr);
          } else {
            if (vflag) eval<0,1,0,1>(ifrom, ito, thr);
            else eval<0,1,0,0>(ifrom, ito, thr);
          }
        } else eval<0,0,0,0>(ifrom, ito, thr);
      }
    }

    sync_threads();
    reduce_thr(this, eflag, vflag, thr, nproxy);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int CTABLE, int EVFLAG, int EFLAG, int VFLAG>
void PairLJCutCoulLongProxyTIP4POMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i,j,ii,jj,jnum,itype,jtype,itable;
  int n,vlist[6];
  int iH1,iH2,jH1,jH2;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul;
  double fraction,table;
  double delxOM,delyOM,delzOM;
  double r,rsq,r2inv,r6inv,forcecoul,forcelj,cforce;
  double factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc,ddotf;
  double v[6],xH1[3],xH2[3];
  double fdx,fdy,fdz,f1x,f1y,f1z,fOx,fOy,fOz,fHx,fHy,fHz;
  const double *x1,*x2;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const double * const q = atom->q;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const int tid = thr->get_tid();
  const double * const special_coul = force->special_coul;
  const double * const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  const double cut_coulsqplus = (cut_coul+2.0*qdist) * (cut_coul+2.0*qdist);

  double fxtmp,fytmp,fztmp;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {

    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    // if atom I = water O, set x1 = offset charge site
    // else x1 = x of atom I
    // NOTE: to make this part thread safe, we need to
    // make sure that the hneigh[][] entries only get
    // updated, when all data is in place. worst case,
    // some calculation is repeated, but since the results
    // will be the same, there is no race condition.
    if (itype == typeO) {
      if (hneigh[i][0] < 0) {
        iH1 = atom->map(atom->tag[i] + 1);
        iH2 = atom->map(atom->tag[i] + 2);
        if (iH1 == -1 || iH2 == -1)
          error->one(FLERR,"TIP4P hydrogen is missing");
        if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
          error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
        compute_newsite_thr(x[i],x[iH1],x[iH2],newsite[i]);
        hneigh[i][2] = 1;
        hneigh[i][1] = iH2;
        hneigh[i][0] = iH1;
      } else {
        iH1 = hneigh[i][0];
        iH2 = hneigh[i][1];
        if (hneigh[i][2] == 0) {
          compute_newsite_thr(x[i],x[iH1],x[iH2],newsite[i]);
          hneigh[i][2] = 1;
        }
      }
      x1 = newsite[i];
    } else x1 = x[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      // LJ interaction based on true rsq

      if (rsq < cut_ljsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
        forcelj *= factor_lj * r2inv;

        fxtmp += delx*forcelj;
        fytmp += dely*forcelj;
        fztmp += delz*forcelj;
        f[j][0] -= delx*forcelj;
        f[j][1] -= dely*forcelj;
        f[j][2] -= delz*forcelj;

        if (EFLAG) {
          evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
            offset[itype][jtype];
          evdwl *= factor_lj;
        } else evdwl = 0.0;

        if (EVFLAG) ev_tally_thr(this,i,j,nlocal, /* newton_pair = */ 1,
                                 evdwl,0.0,forcelj,delx,dely,delz,thr);
      }

      // adjust rsq and delxyz for off-site O charge(s) if necessary
      // but only if they are within reach
      // NOTE: to make this part thread safe, we need to
      // make sure that the hneigh[][] entries only get
      // updated, when all data is in place. worst case,
      // some calculation is repeated, but since the results
      // will be the same, there is no race condition.
      if (rsq < cut_coulsqplus) {
        if (itype == typeO || jtype == typeO) {

          // if atom J = water O, set x2 = offset charge site
          // else x2 = x of atom J

          if (jtype == typeO) {
            if (hneigh[j][0] < 0) {
              jH1 = atom->map(atom->tag[j] + 1);
              jH2 = atom->map(atom->tag[j] + 2);
              if (jH1 == -1 || jH2 == -1)
                error->one(FLERR,"TIP4P hydrogen is missing");
              if (atom->type[jH1] != typeH || atom->type[jH2] != typeH)
                error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
              compute_newsite_thr(x[j],x[jH1],x[jH2],newsite[j]);
              hneigh[j][2] = 1;
              hneigh[j][1] = jH2;
              hneigh[j][0] = jH1;
            } else {
              jH1 = hneigh[j][0];
              jH2 = hneigh[j][1];
              if (hneigh[j][2] == 0) {
                compute_newsite_thr(x[j],x[jH1],x[jH2],newsite[j]);
                hneigh[j][2] = 1;
              }
            }
            x2 = newsite[j];
          } else x2 = x[j];

          delx = x1[0] - x2[0];
          dely = x1[1] - x2[1];
          delz = x1[2] - x2[2];
          rsq = delx*delx + dely*dely + delz*delz;
        }

        // Coulombic interaction based on modified rsq

        if (rsq < cut_coulsq) {
          r2inv = 1 / rsq;
          if (CTABLE || rsq <= tabinnersq) {
            r = sqrt(rsq);
            grij = g_ewald * r;
            expm2 = exp(-grij*grij);
            t = 1.0 / (1.0 + EWALD_P*grij);
            erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
            prefactor = qqrd2e * qtmp*q[j]/r;
            forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
            if (factor_coul < 1.0) {
              forcecoul -= (1.0-factor_coul)*prefactor;
            }
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
            table = ftable[itable] + fraction*dftable[itable];
            forcecoul = qtmp*q[j] * table;
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction*dctable[itable];
              prefactor = qtmp*q[j] * table;
              forcecoul -= (1.0-factor_coul)*prefactor;
            }
          }

          cforce = forcecoul * r2inv;

          // if i,j are not O atoms, force is applied directly
          // if i or j are O atoms, force is on fictitious atom & partitioned
          // force partitioning due to Feenstra, J Comp Chem, 20, 786 (1999)
          // f_f = fictitious force, fO = f_f (1 - 2 alpha), fH = alpha f_f
          // preserves total force and torque on water molecule
          // virial = sum(r x F) where each water's atoms are near xi and xj
          // vlist stores 2,4,6 atoms whose forces contribute to virial

          n = 0;

          if (itype != typeO) {
            fxtmp += delx * cforce;
            fytmp += dely * cforce;
            fztmp += delz * cforce;

            if (VFLAG) {
              v[0] = x[i][0] * delx * cforce;
              v[1] = x[i][1] * dely * cforce;
              v[2] = x[i][2] * delz * cforce;
              v[3] = x[i][0] * dely * cforce;
              v[4] = x[i][0] * delz * cforce;
              v[5] = x[i][1] * delz * cforce;
              vlist[n++] = i;
            }

          } else {

            fdx = delx*cforce;
            fdy = dely*cforce;
            fdz = delz*cforce;

            delxOM = x[i][0] - x1[0];
            delyOM = x[i][1] - x1[1];
            delzOM = x[i][2] - x1[2];

            ddotf = (delxOM * fdx + delyOM * fdy + delzOM * fdz) /
              (qdist*qdist);

            f1x = alpha * (fdx - ddotf * delxOM);
            f1y = alpha * (fdy - ddotf * delyOM);
            f1z = alpha * (fdz - ddotf * delzOM);

            fOx = fdx - f1x;
            fOy = fdy - f1y;
            fOz = fdz - f1z;

            fHx = 0.5 * f1x;
            fHy = 0.5 * f1y;
            fHz = 0.5 * f1z;

            fxtmp += fOx;
            fytmp += fOy;
            fztmp += fOz;

            f[iH1][0] += fHx;
            f[iH1][1] += fHy;
            f[iH1][2] += fHz;

            f[iH2][0] += fHx;
            f[iH2][1] += fHy;
            f[iH2][2] += fHz;

            if (VFLAG) {
              domain->closest_image(x[i],x[iH1],xH1);
              domain->closest_image(x[i],x[iH2],xH2);

              v[0] = x[i][0]*fOx + xH1[0]*fHx + xH2[0]*fHx;
              v[1] = x[i][1]*fOy + xH1[1]*fHy + xH2[1]*fHy;
              v[2] = x[i][2]*fOz + xH1[2]*fHz + xH2[2]*fHz;
              v[3] = x[i][0]*fOy + xH1[0]*fHy + xH2[0]*fHy;
              v[4] = x[i][0]*fOz + xH1[0]*fHz + xH2[0]*fHz;
              v[5] = x[i][1]*fOz + xH1[1]*fHz + xH2[1]*fHz;

              vlist[n++] = i;
              vlist[n++] = iH1;
              vlist[n++] = iH2;
            }
          }

          if (jtype != typeO) {
            f[j][0] -= delx * cforce;
            f[j][1] -= dely * cforce;
            f[j][2] -= delz * cforce;

            if (VFLAG) {
              v[0] -= x[j][0] * delx * cforce;
              v[1] -= x[j][1] * dely * cforce;
              v[2] -= x[j][2] * delz * cforce;
              v[3] -= x[j][0] * dely * cforce;
              v[4] -= x[j][0] * delz * cforce;
              v[5] -= x[j][1] * delz * cforce;
              vlist[n++] = j;
            }

          } else {

            fdx = -delx*cforce;
            fdy = -dely*cforce;
            fdz = -delz*cforce;

            delxOM = x[j][0] - x2[0];
            delyOM = x[j][1] - x2[1];
            delzOM = x[j][2] - x2[2];

            ddotf = (delxOM * fdx + delyOM * fdy + delzOM * fdz) /
              (qdist*qdist);

            f1x = alpha * (fdx - ddotf * delxOM);
            f1y = alpha * (fdy - ddotf * delyOM);
            f1z = alpha * (fdz - ddotf * delzOM);

            fOx = fdx - f1x;
            fOy = fdy - f1y;
            fOz = fdz - f1z;

            fHx = 0.5 * f1x;
            fHy = 0.5 * f1y;
            fHz = 0.5 * f1z;

            f[j][0] += fOx;
            f[j][1] += fOy;
            f[j][2] += fOz;

            f[jH1][0] += fHx;
            f[jH1][1] += fHy;
            f[jH1][2] += fHz;

            f[jH2][0] += fHx;
            f[jH2][1] += fHy;
            f[jH2][2] += fHz;

            if (VFLAG) {
              domain->closest_image(x[j],x[jH1],xH1);
              domain->closest_image(x[j],x[jH2],xH2);

              v[0] += x[j][0]*fOx + xH1[0]*fHx + xH2[0]*fHx;
              v[1] += x[j][1]*fOy + xH1[1]*fHy + xH2[1]*fHy;
              v[2] += x[j][2]*fOz + xH1[2]*fHz + xH2[2]*fHz;
              v[3] += x[j][0]*fOy + xH1[0]*fHy + xH2[0]*fHy;
              v[4] += x[j][0]*fOz + xH1[0]*fHz + xH2[0]*fHz;
              v[5] += x[j][1]*fOz + xH1[1]*fHz + xH2[1]*fHz;

              vlist[n++] = j;
              vlist[n++] = jH1;
              vlist[n++] = jH2;
            }
          }

          if (EFLAG) {
            if (CTABLE || rsq <= tabinnersq)
              ecoul = prefactor*erfc;
            else {
              table = etable[itable] + fraction*detable[itable];
              ecoul = qtmp*q[j] * table;
            }
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;

          if (EVFLAG) ev_tally_list_thr(this,n,vlist,ecoul,v,thr);
        }
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ---------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
  compute position xM of fictitious charge site for O atom and 2 H atoms
  return it as xM
------------------------------------------------------------------------- */

void PairLJCutCoulLongProxyTIP4POMP::compute_newsite_thr(const double * xO,
                                                    const double * xH1,
                                                    const double * xH2,
                                                    double * xM) const
{
  double delx1 = xH1[0] - xO[0];
  double dely1 = xH1[1] - xO[1];
  double delz1 = xH1[2] - xO[2];
  domain->minimum_image(delx1,dely1,delz1);

  double delx2 = xH2[0] - xO[0];
  double dely2 = xH2[1] - xO[1];
  double delz2 = xH2[2] - xO[2];
  domain->minimum_image(delx2,dely2,delz2);

  const double prefac = alpha * 0.5;
  xM[0] = xO[0] + prefac * (delx1 + delx2);
  xM[1] = xO[1] + prefac * (dely1 + dely2);
  xM[2] = xO[2] + prefac * (delz1 + delz2);
}

/* ---------------------------------------------------------------------- */

double PairLJCutCoulLongProxyTIP4POMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJCutCoulLongTIP4P::memory_usage();
  return bytes;
}
