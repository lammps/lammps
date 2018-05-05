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
   OPT version: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include <cmath>
#include "pair_lj_cut_tip4p_long_opt.h"
#include "atom.h"
#include "domain.h"
#include "force.h"
#include "error.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairLJCutTIP4PLongOpt::PairLJCutTIP4PLongOpt(LAMMPS *lmp) :
  PairLJCutTIP4PLong(lmp)
{
  single_enable = 0;
  respa_enable = 0;

  // TIP4P cannot compute virial as F dot r
  // due to finding bonded H atoms which are not near O atom

  no_virial_fdotr_compute = 1;
}

/* ---------------------------------------------------------------------- */

void PairLJCutTIP4PLongOpt::compute(int eflag, int vflag)
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

  int i;
  if (neighbor->ago == 0)
    for (i = 0; i < nall; i++) hneigh[i][0] = -1;
  for (i = 0; i < nall; i++) hneigh[i][2] = 0;


  if (!ncoultablebits) {
    if (evflag) {
      if (eflag) {
        if (vflag) return eval<1,1,1,1>();
        else return eval<1,1,1,0>();
      } else {
        if (vflag) return eval<1,1,0,1>();
        else return eval<1,1,0,0>();
      }
    } else return eval<1,0,0,0>();
  } else {
    if (evflag) {
      if (eflag) {
        if (vflag) return eval<0,1,1,1>();
        else return eval<0,1,1,0>();
      } else {
        if (vflag) return eval<0,1,0,1>();
        else return eval<0,1,0,0>();
      }
    } else return eval<0,0,0,0>();
  }
}

/* ---------------------------------------------------------------------- */

template < const int CTABLE, const int EVFLAG,
           const int EFLAG, const int VFLAG>
void PairLJCutTIP4PLongOpt::eval()
{
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul;
  double fraction,table;
  double r,rsq,r2inv,r6inv,forcecoul,forcelj,cforce;
  double factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  double v[6];
  double fdx,fdy,fdz,fOx,fOy,fOz,fHx,fHy,fHz;
  const double *x1,*x2,*xH1,*xH2;

  int *ilist,*jlist,*numneigh,**firstneigh;
  int i,j,ii,jj,inum,jnum,itype,jtype,itable,key;
  int n,vlist[6];
  int iH1,iH2,jH1,jH2;

  evdwl = ecoul = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = atom->f;
  const double * const q = atom->q;
  const tagint * const tag = atom->tag;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_coul = force->special_coul;
  const double * const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  const double cut_coulsqplus = (cut_coul+2.0*qdist) * (cut_coul+2.0*qdist);

  double fxtmp,fytmp,fztmp;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    // if atom I = water O, set x1 = offset charge site
    // else x1 = x of atom I

    if (itype == typeO) {
      if (hneigh[i][0] < 0) {
        iH1 = atom->map(tag[i] + 1);
        iH2 = atom->map(tag[i] + 2);
        hneigh[i][2] = 1;
        if (iH1 == -1 || iH2 == -1)
          error->one(FLERR,"TIP4P hydrogen is missing");
        if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
          error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
        // set iH1,iH2 to closest image to O
        iH1 = domain->closest_image(i,iH1);
        iH2 = domain->closest_image(i,iH2);
        compute_newsite_opt(x[i],x[iH1],x[iH2],newsite[i]);
        hneigh[i][0] = iH1;
        hneigh[i][1] = iH2;
        hneigh[i][2] = 1;

      } else {
        iH1 = hneigh[i][0];
        iH2 = hneigh[i][1];
        if (hneigh[i][2] == 0) {
          hneigh[i][2] = 1;
          compute_newsite_opt(x[i],x[iH1],x[iH2],newsite[i]);
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

        if (EVFLAG) ev_tally(i,j,nlocal,/* newton_pair = */ 1,
                             evdwl,0.0,forcelj,delx,dely,delz);
      }

      // adjust rsq and delxyz for off-site O charge(s) if necessary
      // but only if they are within reach

      if (rsq < cut_coulsqplus) {
        if (itype == typeO || jtype == typeO) {

          // if atom J = water O, set x2 = offset charge site
          // else x2 = x of atom J

          if (jtype == typeO) {
            if (hneigh[j][0] < 0) {
              jH1 = atom->map(tag[j] + 1);
              jH2 = atom->map(tag[j] + 2);
              hneigh[j][2] = 1;
              if (jH1 == -1 || jH2 == -1)
                error->one(FLERR,"TIP4P hydrogen is missing");
              if (atom->type[jH1] != typeH || atom->type[jH2] != typeH)
                error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
              // set jH1,jH2 to closest image to O
              jH1 = domain->closest_image(j,jH1);
              jH2 = domain->closest_image(j,jH2);
              compute_newsite_opt(x[j],x[jH1],x[jH2],newsite[j]);
              hneigh[j][0] = jH1;
              hneigh[j][1] = jH2;
              hneigh[j][2] = 1;

            } else {
              jH1 = hneigh[j][0];
              jH2 = hneigh[j][1];
              if (hneigh[j][2] == 0) {
                hneigh[j][2] = 1;
                compute_newsite_opt(x[j],x[jH1],x[jH2],newsite[j]);
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

          if (EVFLAG) {
            n = 0;
            key = 0;
          }

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
            }
            if (EVFLAG) vlist[n++] = i;

          } else {
            if (EVFLAG) key += 1;

            fdx = delx*cforce;
            fdy = dely*cforce;
            fdz = delz*cforce;

            fOx = fdx*(1-alpha);
            fOy = fdy*(1-alpha);
            fOz = fdz*(1-alpha);

            fHx = 0.5 * alpha * fdx;
            fHy = 0.5 * alpha * fdy;
            fHz = 0.5 * alpha * fdz;

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
              xH1 = x[iH1];
              xH2 = x[iH2];
              v[0] = x[i][0]*fOx + xH1[0]*fHx + xH2[0]*fHx;
              v[1] = x[i][1]*fOy + xH1[1]*fHy + xH2[1]*fHy;
              v[2] = x[i][2]*fOz + xH1[2]*fHz + xH2[2]*fHz;
              v[3] = x[i][0]*fOy + xH1[0]*fHy + xH2[0]*fHy;
              v[4] = x[i][0]*fOz + xH1[0]*fHz + xH2[0]*fHz;
              v[5] = x[i][1]*fOz + xH1[1]*fHz + xH2[1]*fHz;
            }
            if (EVFLAG) {
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
            }
            if (EVFLAG) vlist[n++] = j;

          } else {
            if (EVFLAG) key += 2;

            fdx = -delx*cforce;
            fdy = -dely*cforce;
            fdz = -delz*cforce;

            fOx = fdx*(1-alpha);
            fOy = fdy*(1-alpha);
            fOz = fdz*(1-alpha);

            fHx = 0.5 * alpha * fdx;
            fHy = 0.5 * alpha * fdy;
            fHz = 0.5 * alpha * fdz;

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
              xH1 = x[jH1];
              xH2 = x[jH2];
              v[0] += x[j][0]*fOx + xH1[0]*fHx + xH2[0]*fHx;
              v[1] += x[j][1]*fOy + xH1[1]*fHy + xH2[1]*fHy;
              v[2] += x[j][2]*fOz + xH1[2]*fHz + xH2[2]*fHz;
              v[3] += x[j][0]*fOy + xH1[0]*fHy + xH2[0]*fHy;
              v[4] += x[j][0]*fOz + xH1[0]*fHz + xH2[0]*fHz;
              v[5] += x[j][1]*fOz + xH1[1]*fHz + xH2[1]*fHz;
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
          if (EVFLAG) ev_tally_tip4p(key,vlist,v,ecoul,alpha);
        }
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ----------------------------------------------------------------------
  compute position xM of fictitious charge site for O atom and 2 H atoms
  return it as xM
------------------------------------------------------------------------- */

void PairLJCutTIP4PLongOpt::compute_newsite_opt(const double * xO,
                                                    const double * xH1,
                                                    const double * xH2,
                                                    double * xM) const
{
  double delx1 = xH1[0] - xO[0];
  double dely1 = xH1[1] - xO[1];
  double delz1 = xH1[2] - xO[2];

  double delx2 = xH2[0] - xO[0];
  double dely2 = xH2[1] - xO[1];
  double delz2 = xH2[2] - xO[2];

  const double prefac = alpha * 0.5;
  xM[0] = xO[0] + prefac * (delx1 + delx2);
  xM[1] = xO[1] + prefac * (dely1 + dely2);
  xM[2] = xO[2] + prefac * (delz1 + delz2);
}

/* ---------------------------------------------------------------------- */

double PairLJCutTIP4PLongOpt::memory_usage()
{
  double bytes = PairLJCutTIP4PLong::memory_usage();

  return bytes;
}
