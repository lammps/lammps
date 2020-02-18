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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include "pair_spin_neel.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSpinNeel::~PairSpinNeel()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cut_spin_neel);
    memory->destroy(g1);
    memory->destroy(g1_mech);
    memory->destroy(g2);
    memory->destroy(g3);
    memory->destroy(q1);
    memory->destroy(q1_mech);
    memory->destroy(q2);
    memory->destroy(q3);
    memory->destroy(cutsq); // to be deleted
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSpinNeel::settings(int narg, char **arg)
{
  PairSpin::settings(narg,arg);

  cut_spin_neel_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = i+1; j <= atom->ntypes; j++) {
        if (setflag[i][j]) {
          cut_spin_neel[i][j] = cut_spin_neel_global;
        }
      }
    }
  }

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type spin pairs (only one for now)
------------------------------------------------------------------------- */

void PairSpinNeel::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  // check if args correct

  if (strcmp(arg[2],"neel") != 0)
    error->all(FLERR,"Incorrect args in pair_style command");
  if (narg != 10)
    error->all(FLERR,"Incorrect args in pair_style command");

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  const double rij = force->numeric(FLERR,arg[3]);
  const double k1 = force->numeric(FLERR,arg[4]);
  const double k2 = force->numeric(FLERR,arg[5]);
  const double k3 = force->numeric(FLERR,arg[6]);
  const double l1 = force->numeric(FLERR,arg[7]);
  const double l2 = force->numeric(FLERR,arg[8]);
  const double l3 = force->numeric(FLERR,arg[9]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut_spin_neel[i][j] = rij;
      g1[i][j] = k1/hbar;
      q1[i][j] = l1/hbar;
      g1_mech[i][j] = k1;
      q1_mech[i][j] = l1;
      g2[i][j] = k2;
      g3[i][j] = k3;
      q2[i][j] = l2;
      q3[i][j] = l3;
      setflag[i][j] = 1;
      count++;
    }
  }
  if (count == 0)
    error->all(FLERR,"Incorrect args in pair_style command");

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSpinNeel::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  g1[j][i] = g1[i][j];
  g1_mech[j][i] = g1_mech[i][j];
  g2[j][i] = g2[i][j];
  g3[j][i] = g3[i][j];
  q1[j][i] = q1[i][j];
  q1_mech[j][i] = q1_mech[i][j];
  q2[j][i] = q2[i][j];
  q3[j][i] = q3[i][j];

  return cut_spin_neel_global;
}

/* ----------------------------------------------------------------------
   extract the larger cutoff
------------------------------------------------------------------------- */

void *PairSpinNeel::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut") == 0) return (void *) &cut_spin_neel_global;
  return NULL;
}

/* ---------------------------------------------------------------------- */

void PairSpinNeel::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double evdwl,ecoul;
  double xi[3], rij[3], eij[3];
  double spi[3], spj[3];
  double fi[3], fmi[3];
  double local_cut2;
  double rsq, inorm;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double **fm = atom->fm;
  double **sp = atom->sp;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // computation of the neel interaction
  // loop over atoms and their neighbors

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    xi[0] = x[i][0];
    xi[1] = x[i][1];
    xi[2] = x[i][2];
    spi[0] = sp[i][0];
    spi[1] = sp[i][1];
    spi[2] = sp[i][2];

    // loop on neighbors

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      spj[0] = sp[j][0];
      spj[1] = sp[j][1];
      spj[2] = sp[j][2];

      evdwl = 0.0;

      fi[0] = fi[1] = fi[2] = 0.0;
      fmi[0] = fmi[1] = fmi[2] = 0.0;
      rij[0] = rij[1] = rij[2] = 0.0;

      rij[0] = x[j][0] - xi[0];
      rij[1] = x[j][1] - xi[1];
      rij[2] = x[j][2] - xi[2];
      rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
      inorm = 1.0/sqrt(rsq);
      eij[0] = rij[0]*inorm;
      eij[1] = rij[1]*inorm;
      eij[2] = rij[2]*inorm;

      itype = type[i];
      jtype = type[j];

      local_cut2 = cut_spin_neel[itype][jtype]*cut_spin_neel[itype][jtype];

      // compute neel interaction

      if (rsq <= local_cut2) {
        compute_neel(i,j,rsq,eij,fmi,spi,spj);
        if (lattice_flag) {
          compute_neel_mech(i,j,rsq,eij,fi,spi,spj);
        }
      }

      f[i][0] += fi[0];
      f[i][1] += fi[1];
      f[i][2] += fi[2];
      fm[i][0] += fmi[0];
      fm[i][1] += fmi[1];
      fm[i][2] += fmi[2];

      if (newton_pair || j < nlocal) {
        f[j][0] -= fi[0];
        f[j][1] -= fi[1];
        f[j][2] -= fi[2];
      }

      if (eflag) {
        // evdwl = (spi[0]*fmi[0] + spi[1]*fmi[1] + spi[2]*fmi[2]);
        evdwl = compute_neel_energy(i,j,rsq,eij,spi,spj);
        evdwl *= 0.5*hbar;
      } else evdwl = 0.0;

      if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
          evdwl,ecoul,fi[0],fi[1],fi[2],rij[0],rij[1],rij[2]);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();

}

/* ----------------------------------------------------------------------
   update the pair interactions fmi acting on the spin ii
------------------------------------------------------------------------- */

void PairSpinNeel::compute_single_pair(int ii, double fmi[3])
{
  int *type = atom->type;
  double **x = atom->x;
  double **sp = atom->sp;
  double local_cut2;

  double xi[3], rij[3], eij[3];
  double spi[3], spj[3];

  int j,jnum,itype,jtype,ntypes;
  int k,locflag;
  int *jlist,*numneigh,**firstneigh;

  double rsq, inorm;

  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // check if interaction applies to type of ii

  itype = type[ii];
  ntypes = atom->ntypes;
  locflag = 0;
  k = 1;
  while (k <= ntypes) {
    if (k <= itype) {
      if (setflag[k][itype] == 1) {
        locflag =1;
        break;
      }
      k++;
    } else if (k > itype) {
      if (setflag[itype][k] == 1) {
        locflag =1;
        break;
      }
      k++;
    } else error->all(FLERR,"Wrong type number");
  }

  // if interaction applies to type ii,
  // locflag = 1 and compute pair interaction

  if (locflag == 1) {

    spi[0] = sp[ii][0];
    spi[1] = sp[ii][1];
    spi[2] = sp[ii][2];

    xi[0] = x[ii][0];
    xi[1] = x[ii][1];
    xi[2] = x[ii][2];

    eij[0] = eij[1] = eij[2] = 0.0;

    jlist = firstneigh[ii];
    jnum = numneigh[ii];

    for (int jj = 0; jj < jnum; jj++) {

      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      local_cut2 = cut_spin_neel[itype][jtype]*cut_spin_neel[itype][jtype];

      spj[0] = sp[j][0];
      spj[1] = sp[j][1];
      spj[2] = sp[j][2];

      rij[0] = x[j][0] - xi[0];
      rij[1] = x[j][1] - xi[1];
      rij[2] = x[j][2] - xi[2];
      rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
      inorm = 1.0/sqrt(rsq);
      eij[0] = inorm*rij[0];
      eij[1] = inorm*rij[1];
      eij[2] = inorm*rij[2];

      if (rsq <= local_cut2) {
        compute_neel(ii,j,rsq,eij,fmi,spi,spj);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairSpinNeel::compute_neel(int i, int j, double rsq, double eij[3], double fmi[3], double spi[3], double spj[3])
{
  int *type = atom->type;
  int itype, jtype;
  itype = type[i];
  jtype = type[j];

  double qr,gr,g1r,q1r,q2r,ra;
  double pdx, pdy, pdz;
  double pq1x, pq1y, pq1z;
  double pq2x, pq2y, pq2z;
  double eij_si,eij_sj,si_sj,eij_si_2,eij_sj_3,coeff1;

  // compute Neel's functions

  ra = rsq/g3[itype][jtype]/g3[itype][jtype];
  gr = 4.0*g1[itype][jtype]*ra;
  gr *= (1.0-g2[itype][jtype]*ra);
  gr *= exp(-ra);

  ra = rsq/q3[itype][jtype]/q3[itype][jtype];
  qr = 4.0*q1[itype][jtype]*ra;
  qr *= (1.0-q2[itype][jtype]*ra);
  qr *= exp(-ra);

  g1r = (gr + 12.0*qr/35.0);
  q1r = 9.0*qr/5.0;
  q2r = -2.0*qr/5.0;

  // pseudo-dipolar component
  
  eij_si = eij[0]*spi[0] + eij[1]*spi[1] + eij[2]*spi[2];
  eij_sj = eij[0]*spj[0] + eij[1]*spj[1] + eij[2]*spj[2];
  si_sj = spi[0]*spj[0] + spi[1]*spj[1] + spi[2]*spj[2];

  pdx = g1r*(eij_sj*eij[0] - spj[0]/3.0);
  pdy = g1r*(eij_sj*eij[1] - spj[1]/3.0);
  pdz = g1r*(eij_sj*eij[2] - spj[2]/3.0);

  // pseudo-quadrupolar components

  eij_si_2 = eij_si*eij_si;
  pq1x = -(eij_si_2 - si_sj/3.0)*spj[0]/3.0;
  pq1y = -(eij_si_2 - si_sj/3.0)*spj[1]/3.0;
  pq1z = -(eij_si_2 - si_sj/3.0)*spj[2]/3.0;

  coeff1 = (eij_sj*eij_sj-si_sj/3.0);
  pq1x += coeff1*(2.0*eij_si*eij[0] - spj[0]/3.0);
  pq1y += coeff1*(2.0*eij_si*eij[1] - spj[1]/3.0);
  pq1z += coeff1*(2.0*eij_si*eij[2] - spj[2]/3.0);

  pq1x *= q1r;
  pq1y *= q1r;
  pq1z *= q1r;

  eij_sj_3 = eij_sj*eij_sj*eij_sj;
  pq2x = 3.0*eij_si_2*eij_sj*eij[0] + eij_sj_3*eij[0];
  pq2y = 3.0*eij_si_2*eij_sj*eij[1] + eij_sj_3*eij[1];
  pq2z = 3.0*eij_si_2*eij_sj*eij[2] + eij_sj_3*eij[2];

  pq2x *= q2r;
  pq2y *= q2r;
  pq2z *= q2r;

  // adding three contributions

  fmi[0] += (pdx + pq1x + pq2x);
  fmi[1] += (pdy + pq1y + pq2y);
  fmi[2] += (pdz + pq1z + pq2z);
}

/* ---------------------------------------------------------------------- */

void PairSpinNeel::compute_neel_mech(int i, int j, double rsq, double eij[3], double fi[3],  double spi[3], double spj[3])
{
  int *type = atom->type;
  int itype, jtype;
  itype = type[i];
  jtype = type[j];

  double g_mech, gij, dgij;
  double q_mech, q1ij, dq1ij;
  double q2ij, dq2ij;
  double pdx, pdy, pdz;
  double pq1x, pq1y, pq1z;
  double pq2x, pq2y, pq2z;
  double ra, rr, drij, ig3, iq3;

  drij = sqrt(rsq);

  double scalar_si_sj = spi[0]*spj[0]+spi[1]*spj[1]+spi[2]*spj[2];
  double scalar_eij_si = eij[0]*spi[0]+eij[1]*spi[1]+eij[2]*spi[2];
  double scalar_eij_sj = eij[0]*spj[0]+eij[1]*spj[1]+eij[2]*spj[2];

  // pseudo-dipolar component

  g_mech = g1_mech[itype][jtype];
  ig3 = 1.0/(g3[itype][jtype]*g3[itype][jtype]);

  ra = rsq*ig3;
  rr = drij*ig3;

  gij = 4.0*g_mech*ra;
  gij *= (1.0-g2[itype][jtype]*ra);
  gij *= exp(-ra);

  dgij = 1.0-ra-g2[itype][jtype]*ra*(2.0-ra);
  dgij *= 8.0*g_mech*rr*exp(-ra);

  double pdt1 = (dgij-2.0*gij/drij)*scalar_eij_si*scalar_eij_sj;
  pdt1 -= scalar_si_sj*dgij/3.0;
  double pdt2 = scalar_eij_sj*gij/drij;
  double pdt3 = scalar_eij_si*gij/drij;
  pdx = -(pdt1*eij[0] + pdt2*spi[0] + pdt3*spj[0]);
  pdy = -(pdt1*eij[1] + pdt2*spi[1] + pdt3*spj[1]);
  pdz = -(pdt1*eij[2] + pdt2*spi[2] + pdt3*spj[2]);

  // pseudo-quadrupolar component

  q_mech = q1_mech[itype][jtype];
  iq3 = 1.0/(q3[itype][jtype]*q3[itype][jtype]);

  ra = rsq*iq3;
  rr = drij*iq3;

  q1ij = 4.0*q_mech*ra;
  q1ij *= (1.0-q2[itype][jtype]*ra);
  q1ij *= exp(-ra);
  q2ij = -2.0*q1ij/9.0;

  dq1ij = 1.0-ra-q2[itype][jtype]*ra*(2.0-ra);
  dq1ij *= 8.0*q_mech*rr*exp(-ra);
  dq2ij = -2.0*dq1ij/9.0;

  double scalar_eij_si_2 = scalar_eij_si*scalar_eij_si;
  double scalar_eij_sj_2 = scalar_eij_sj*scalar_eij_sj;
  double pqt1 = scalar_eij_si_2 - scalar_si_sj/3.0;
  double pqt2 = scalar_eij_sj_2 - scalar_si_sj/3.0;
  pq1x = dq1ij * pqt1 * pqt2 * eij[0];
  pq1y = dq1ij * pqt1 * pqt2 * eij[1];
  pq1z = dq1ij * pqt1 * pqt2 * eij[2];

  double scalar_eij_si_3 = scalar_eij_si*scalar_eij_si*scalar_eij_si;
  double scalar_eij_sj_3 = scalar_eij_sj*scalar_eij_sj*scalar_eij_sj;
  double scalar_si_sj_2 = scalar_si_sj*scalar_si_sj;
/*  double pqt3 = 2.0*scalar_eij_si*scalar_eij_sj_3/drij;
  double pqt4 = 2.0*scalar_eij_sj*scalar_eij_si_3/drij;
  double pqt5 = -2.0*scalar_si_sj*scalar_eij_si/(3.0*drij);
  double pqt6 = -2.0*scalar_si_sj*scalar_eij_sj/(3.0*drij);
//  pq1x += q1ij*(spi[0]*(pqt3+pqt6) + spj[0]*(pqt4+));
  pq1x += q1ij*(pqt3*spi[0]+pqt4*spj[0]+pqt5*spi[0]+pqt6*spi[0]);
  pq1y += q1ij*(pqt3*spi[1]+pqt4*spj[1]+pqt5*spi[1]+pqt6*spj[1]);
  pq1z += q1ij*(pqt3*spi[2]+pqt4*spj[2]+pqt5*spi[2]+pqt6*spj[2]);
*/
  double pqt3 = 2.0*scalar_eij_si*(scalar_eij_sj_2-scalar_si_sj/3.0)/drij;
  double pqt4 = 2.0*scalar_eij_sj*(scalar_eij_si_2-scalar_si_sj/3.0)/drij;
  pq1x += q1ij*(pqt3*spi[0] + pqt4*spj[0]);
  pq1y += q1ij*(pqt3*spi[1] + pqt4*spj[1]);
  pq1z += q1ij*(pqt3*spi[2] + pqt4*spj[2]);
  double pqt7 = 4.0*scalar_eij_si_2*scalar_eij_sj_2/drij;
  double pqt8 = 2.0*scalar_si_sj_2*scalar_eij_sj/(3.0*drij);
  double pqt9 = 2.0*scalar_si_sj_2*scalar_eij_si/(3.0*drij);
  pq1x -= q1ij*(pqt7 + pqt8 + pqt9)*eij[0];
  pq1y -= q1ij*(pqt7 + pqt8 + pqt9)*eij[1];
  pq1z -= q1ij*(pqt7 + pqt8 + pqt9)*eij[2];

/*
  double pqt3 = 2.0*scalar_eij_si*(scalar_eij_sj_2-scalar_si_sj/3.0)/drij;
  double pqt4 = 2.0*scalar_eij_sj*(scalar_eij_si_2-scalar_si_sj/3.0)/drij;
  pq1x += q1ij*(pqt3*spi[0] + pqt4*spj[0]);
  pq1y += q1ij*(pqt3*spi[1] + pqt4*spj[1]);
  pq1z += q1ij*(pqt3*spi[2] + pqt4*spj[2]);
*/

  //double scalar_eij_si_3 = scalar_eij_si*scalar_eij_si*scalar_eij_si;
  //double scalar_eij_sj_3 = scalar_eij_sj*scalar_eij_sj*scalar_eij_sj;
  double pqt10 = scalar_eij_sj*scalar_eij_si_3;
  double pqt11 = scalar_eij_si*scalar_eij_sj_3;
  pq2x = dq2ij*(pqt10 + pqt11)*eij[0];
  pq2y = dq2ij*(pqt10 + pqt11)*eij[1];
  pq2z = dq2ij*(pqt10 + pqt11)*eij[2];

  double pqt12 = scalar_eij_si_3/drij;
  double pqt13 = scalar_eij_sj_3/drij;
  double pqt14 = 3.0*scalar_eij_sj*scalar_eij_si_2/drij;
  double pqt15 = 3.0*scalar_eij_si*scalar_eij_sj_2/drij;
  pq2x += q2ij*((pqt12+pqt15)*spj[0]+(pqt13+pqt14)*spi[0]);
  pq2y += q2ij*((pqt12+pqt15)*spj[1]+(pqt13+pqt14)*spi[1]);
  pq2z += q2ij*((pqt12+pqt15)*spj[2]+(pqt13+pqt14)*spi[2]);
  double pqt16 = 4.0*scalar_eij_sj*scalar_eij_si_3/drij;
  double pqt17 = 4.0*scalar_eij_si*scalar_eij_sj_3/drij;
  pq2x -= q2ij*(pqt16 + pqt17)*eij[0];
  pq2y -= q2ij*(pqt16 + pqt17)*eij[1];
  pq2z -= q2ij*(pqt16 + pqt17)*eij[2];

  // adding three contributions

  fi[0] = pdx + pq1x + pq2x;
  fi[1] = pdy + pq1y + pq2y;
  fi[2] = pdz + pq1z + pq2z;
}

/* ---------------------------------------------------------------------- */

double PairSpinNeel::compute_neel_energy(int i, int j, double rsq, double eij[3], double spi[3], double spj[3]) 
{
  int *type = atom->type;
  int itype, jtype;
  itype = type[i];
  jtype = type[j];

  double qr,gr,g1r,q1r,q2r,ra;
  double epd,epq1,epq2;
  double eij_si,eij_sj,si_sj;
  double eij_si_2,eij_sj_2,eij_si_3,eij_sj_3;

  // compute Neel's functions

  ra = rsq/g3[itype][jtype]/g3[itype][jtype];
  gr = 4.0*g1[itype][jtype]*ra;
  gr *= (1.0-g2[itype][jtype]*ra);
  gr *= exp(-ra);

  ra = rsq/q3[itype][jtype]/q3[itype][jtype];
  qr = 4.0*q1[itype][jtype]*ra;
  qr *= (1.0-q2[itype][jtype]*ra);
  qr *= exp(-ra);

  g1r = (gr + 12.0*qr/35.0);
  q1r = 9.0*qr/5.0;
  q2r = -2.0*qr/5.0;

  eij_si = eij[0]*spi[0] + eij[1]*spi[1] + eij[2]*spi[2];
  eij_sj = eij[0]*spj[0] + eij[1]*spj[1] + eij[2]*spj[2];
  si_sj = spi[0]*spj[0] + spi[1]*spj[1] + spi[2]*spj[2];
  epd = g1r*(eij_si*eij_sj-si_sj/3.0);
  eij_si_2 = eij_si*eij_si;
  eij_sj_2 = eij_sj*eij_sj;
  epq1 = q1r*(eij_si_2-si_sj/3.0)*(eij_sj_2-si_sj/3.0);
  eij_si_3 = eij_si*eij_si_2;
  eij_sj_3 = eij_sj*eij_sj_2;
  epq2 = q2r*(eij_si*eij_sj_3+eij_sj*eij_si_3);

  return (epd+epq1+epq2);
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSpinNeel::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cut_spin_neel,n+1,n+1,"pair/spin/soc/neel:cut_spin_neel");

  memory->create(g1,n+1,n+1,"pair/spin/soc/neel:g1");
  memory->create(g1_mech,n+1,n+1,"pair/spin/soc/neel:g1_mech");
  memory->create(g2,n+1,n+1,"pair/spin/soc/neel:g2");
  memory->create(g3,n+1,n+1,"pair/spin/soc/neel:g3");

  memory->create(q1,n+1,n+1,"pair/spin/soc/neel:q1");
  memory->create(q1_mech,n+1,n+1,"pair/spin/soc/neel:q1_mech");
  memory->create(q2,n+1,n+1,"pair/spin/soc/neel:q2");
  memory->create(q3,n+1,n+1,"pair/spin/soc/neel:q3");

  memory->create(cutsq,n+1,n+1,"pair/spin/soc/neel:cutsq");
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinNeel::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&g1[i][j],sizeof(double),1,fp);
        fwrite(&g1_mech[i][j],sizeof(double),1,fp);
        fwrite(&g2[i][j],sizeof(double),1,fp);
        fwrite(&g3[i][j],sizeof(double),1,fp);
        fwrite(&q1[i][j],sizeof(double),1,fp);
        fwrite(&q1_mech[i][j],sizeof(double),1,fp);
        fwrite(&q2[i][j],sizeof(double),1,fp);
        fwrite(&q3[i][j],sizeof(double),1,fp);
        fwrite(&cut_spin_neel[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinNeel::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,NULL,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&g1[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&g1_mech[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&g2[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&g3[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&q1[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&q1_mech[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&q2[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&q3[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&cut_spin_neel[i][j],sizeof(double),1,fp,NULL,error);
        }
        MPI_Bcast(&g1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&g1_mech[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&g2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&g3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&q1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&q1_mech[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&q2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&q3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_spin_neel[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinNeel::write_restart_settings(FILE *fp)
{
  fwrite(&cut_spin_neel_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinNeel::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_spin_neel_global,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,NULL,error);
  }
  MPI_Bcast(&cut_spin_neel_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}
