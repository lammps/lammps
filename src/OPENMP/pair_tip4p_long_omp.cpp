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

#include "omp_compat.h"
#include <cmath>
#include "pair_tip4p_long_omp.h"
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "error.h"
#include "memory.h"
#include "neigh_list.h"

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

PairTIP4PLongOMP::PairTIP4PLongOMP(LAMMPS *lmp) :
  PairTIP4PLong(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  newsite_thr = nullptr;
  hneigh_thr = nullptr;

  // TIP4P cannot compute virial as F dot r
  // due to finding bonded H atoms which are not near O atom

  no_virial_fdotr_compute = 1;
}

/* ---------------------------------------------------------------------- */

PairTIP4PLongOMP::~PairTIP4PLongOMP()
{
  memory->destroy(hneigh_thr);
  memory->destroy(newsite_thr);
}

/* ---------------------------------------------------------------------- */

void PairTIP4PLongOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  // reallocate hneigh_thr & newsite_thr if necessary
  // initialize hneigh_thr[0] to -1 on steps when reneighboring occurred
  // initialize hneigh_thr[2] to 0 every step

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(hneigh_thr);
    memory->create(hneigh_thr,nmax,"pair:hneigh_thr");
    memory->destroy(newsite_thr);
    memory->create(newsite_thr,nmax,"pair:newsite_thr");
  }

  int i;
  // tag entire list as completely invalid after a neighbor
  // list update, since that can change the order of atoms.
  if (neighbor->ago == 0)
    for (i = 0; i < nall; i++) hneigh_thr[i].a = -1;

  // indicate that the coordinates for the M point need to
  // be updated. this needs to be done in every step.
  for (i = 0; i < nall; i++) hneigh_thr[i].t = 0;

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

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int CTABLE, int EVFLAG, int EFLAG, int VFLAG>
void PairTIP4PLongOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul;
  double fraction,table;
  double r,rsq,r2inv,forcecoul,cforce;
  double factor_coul;
  double grij,expm2,prefactor,t,erfc;
  double v[6];
  double fdx,fdy,fdz,fOx,fOy,fOz,fHx,fHy,fHz;
  dbl3_t x1,x2,xH1,xH2;

  int *ilist,*jlist,*numneigh,**firstneigh;
  int i,j,ii,jj,jnum,itype,jtype,itable, key;
  int n,vlist[6];
  int iH1,iH2,jH1,jH2;

  ecoul = 0.0;

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const double * _noalias const q = atom->q;
  const int * _noalias const type = atom->type;
  const double * _noalias const special_coul = force->special_coul;
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
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    itype = type[i];

    // if atom I = water O, set x1 = offset charge site
    // else x1 = x of atom I
    // NOTE: to make this part thread safe, we need to
    // make sure that the hneigh_thr[][] entries only get
    // updated, when all data is in place. worst case,
    // some calculation is repeated, but since the results
    // will be the same, there is no race condition.
    if (itype == typeO) {
      if (hneigh_thr[i].a < 0) {
        iH1 = atom->map(atom->tag[i] + 1);
        iH2 = atom->map(atom->tag[i] + 2);
        if (iH1 == -1 || iH2 == -1)
          error->one(FLERR,"TIP4P hydrogen is missing");
        if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
          error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
        // set iH1,iH2 to index of closest image to O
        iH1 = domain->closest_image(i,iH1);
        iH2 = domain->closest_image(i,iH2);
        compute_newsite_thr(x[i],x[iH1],x[iH2],newsite_thr[i]);
        hneigh_thr[i].t = 1;
        hneigh_thr[i].b = iH2;
        hneigh_thr[i].a = iH1;

      } else {
        iH1 = hneigh_thr[i].a;
        iH2 = hneigh_thr[i].b;
        if (hneigh_thr[i].t == 0) {
          compute_newsite_thr(x[i],x[iH1],x[iH2],newsite_thr[i]);
          hneigh_thr[i].t = 1;
        }
      }
      x1 = newsite_thr[i];
    } else x1 = x[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp=fytmp=fztmp=0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      // adjust rsq and delxyz for off-site O charge(s) if necessary
      // but only if they are within reach
      // NOTE: to make this part thread safe, we need to
      // make sure that the hneigh_thr[][] entries only get
      // updated, when all data is in place. worst case,
      // some calculation is repeated, but since the results
      // will be the same, there is no race condition.
      if (rsq < cut_coulsqplus) {
        if (itype == typeO || jtype == typeO) {

          // if atom J = water O, set x2 = offset charge site
          // else x2 = x of atom J

          if (jtype == typeO) {
            if (hneigh_thr[j].a < 0) {
              jH1 = atom->map(atom->tag[j] + 1);
              jH2 = atom->map(atom->tag[j] + 2);
              if (jH1 == -1 || jH2 == -1)
                error->one(FLERR,"TIP4P hydrogen is missing");
              if (atom->type[jH1] != typeH || atom->type[jH2] != typeH)
                error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
              // set jH1,jH2 to closest image to O
              jH1 = domain->closest_image(j,jH1);
              jH2 = domain->closest_image(j,jH2);
              compute_newsite_thr(x[j],x[jH1],x[jH2],newsite_thr[j]);
              hneigh_thr[j].t = 1;
              hneigh_thr[j].b = jH2;
              hneigh_thr[j].a = jH1;

            } else {
              jH1 = hneigh_thr[j].a;
              jH2 = hneigh_thr[j].b;
              if (hneigh_thr[j].t == 0) {
                compute_newsite_thr(x[j],x[jH1],x[jH2],newsite_thr[j]);
                hneigh_thr[j].t = 1;
              }
            }
            x2 = newsite_thr[j];
          } else x2 = x[j];

          delx = x1.x - x2.x;
          dely = x1.y - x2.y;
          delz = x1.z - x2.z;
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

          if (EVFLAG) {
            n = 0;
            key = 0;
          }

          if (itype != typeO) {
            fxtmp += delx * cforce;
            fytmp += dely * cforce;
            fztmp += delz * cforce;

            if (VFLAG) {
              v[0] = x[i].x * delx * cforce;
              v[1] = x[i].y * dely * cforce;
              v[2] = x[i].z * delz * cforce;
              v[3] = x[i].x * dely * cforce;
              v[4] = x[i].x * delz * cforce;
              v[5] = x[i].y * delz * cforce;
            }
            if (EVFLAG) vlist[n++] = i;

          } else {
            if (EVFLAG) key++;

            fdx = delx*cforce;
            fdy = dely*cforce;
            fdz = delz*cforce;

            fOx = fdx*(1 - alpha);
            fOy = fdy*(1 - alpha);
            fOz = fdz*(1 - alpha);

            fHx = 0.5*alpha * fdx;
            fHy = 0.5*alpha * fdy;
            fHz = 0.5*alpha * fdz;

            fxtmp += fOx;
            fytmp += fOy;
            fztmp += fOz;

            f[iH1].x += fHx;
            f[iH1].y += fHy;
            f[iH1].z += fHz;

            f[iH2].x += fHx;
            f[iH2].y += fHy;
            f[iH2].z += fHz;

            if (VFLAG) {
              xH1 = x[iH1];
              xH2 = x[iH2];
              v[0] = x[i].x*fOx + xH1.x*fHx + xH2.x*fHx;
              v[1] = x[i].y*fOy + xH1.y*fHy + xH2.y*fHy;
              v[2] = x[i].z*fOz + xH1.z*fHz + xH2.z*fHz;
              v[3] = x[i].x*fOy + xH1.x*fHy + xH2.x*fHy;
              v[4] = x[i].x*fOz + xH1.x*fHz + xH2.x*fHz;
              v[5] = x[i].y*fOz + xH1.y*fHz + xH2.y*fHz;
            }
            if (EVFLAG) {
              vlist[n++] = i;
              vlist[n++] = iH1;
              vlist[n++] = iH2;
            }
          }

          if (jtype != typeO) {
            f[j].x -= delx * cforce;
            f[j].y -= dely * cforce;
            f[j].z -= delz * cforce;

            if (VFLAG) {
              v[0] -= x[j].x * delx * cforce;
              v[1] -= x[j].y * dely * cforce;
              v[2] -= x[j].z * delz * cforce;
              v[3] -= x[j].x * dely * cforce;
              v[4] -= x[j].x * delz * cforce;
              v[5] -= x[j].y * delz * cforce;
            }
            if (EVFLAG) vlist[n++] = j;

          } else {
            if (EVFLAG) key += 2;

            fdx = -delx*cforce;
            fdy = -dely*cforce;
            fdz = -delz*cforce;

            fOx = fdx*(1 - alpha);
            fOy = fdy*(1 - alpha);
            fOz = fdz*(1 - alpha);

            fHx = 0.5*alpha * fdx;
            fHy = 0.5*alpha * fdy;
            fHz = 0.5*alpha * fdz;

            f[j].x += fOx;
            f[j].y += fOy;
            f[j].z += fOz;

            f[jH1].x += fHx;
            f[jH1].y += fHy;
            f[jH1].z += fHz;

            f[jH2].x += fHx;
            f[jH2].y += fHy;
            f[jH2].z += fHz;

            if (VFLAG) {
              xH1 = x[jH1];
              xH2 = x[jH2];
              v[0] += x[j].x*fOx + xH1.x*fHx + xH2.x*fHx;
              v[1] += x[j].y*fOy + xH1.y*fHy + xH2.y*fHy;
              v[2] += x[j].z*fOz + xH1.z*fHz + xH2.z*fHz;
              v[3] += x[j].x*fOy + xH1.x*fHy + xH2.x*fHy;
              v[4] += x[j].x*fOz + xH1.x*fHz + xH2.x*fHz;
              v[5] += x[j].y*fOz + xH1.y*fHz + xH2.y*fHz;
            }
            if (EVFLAG) {
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

          if (EVFLAG) ev_tally_list_thr(this,key,vlist,v,ecoul,alpha,thr);
        }
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ----------------------------------------------------------------------
  compute position xM of fictitious charge site for O atom and 2 H atoms
  return it as xM
------------------------------------------------------------------------- */

void PairTIP4PLongOMP::compute_newsite_thr(const dbl3_t &xO,
                                                const dbl3_t &xH1,
                                                const dbl3_t &xH2,
                                                dbl3_t &xM) const
{
  double delx1 = xH1.x - xO.x;
  double dely1 = xH1.y - xO.y;
  double delz1 = xH1.z - xO.z;

  double delx2 = xH2.x - xO.x;
  double dely2 = xH2.y - xO.y;
  double delz2 = xH2.z - xO.z;

  const double prefac = alpha * 0.5;
  xM.x = xO.x + prefac * (delx1 + delx2);
  xM.y = xO.y + prefac * (dely1 + dely2);
  xM.z = xO.z + prefac * (delz1 + delz2);
}

/* ---------------------------------------------------------------------- */

double PairTIP4PLongOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairTIP4PLong::memory_usage();
  return bytes;
}
