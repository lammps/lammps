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

#include <cmath>
#include <cstdlib>
#include <cstring>

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "fix.h"
#include "fix_nve_spin.h"
#include "force.h"
#include "pair_hybrid.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "pair_spin_neel.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairSpinNeel::PairSpinNeel(LAMMPS *lmp) : PairSpin(lmp),
lockfixnvespin(NULL)
{
  single_enable = 0;
  no_virial_fdotr_compute = 1;
  lattice_flag = 0;
}

/* ---------------------------------------------------------------------- */

PairSpinNeel::~PairSpinNeel()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(rs);
    memory->destroy(rl);
    memory->destroy(dr);
    memory->destroy(r0);

    memory->destroy(g1);
    memory->destroy(g2);
    memory->destroy(g3);
    memory->destroy(g4);
    memory->destroy(g5);
    memory->destroy(gm1);
    memory->destroy(gm2);
    memory->destroy(gm3);
    memory->destroy(gm4);
    memory->destroy(gm5);
    
    memory->destroy(q1);
    memory->destroy(q2);
    memory->destroy(q3);
    memory->destroy(q4);
    memory->destroy(q5);
    memory->destroy(qm1);
    memory->destroy(qm2);
    memory->destroy(qm3);
    memory->destroy(qm4);
    memory->destroy(qm5);
    memory->destroy(cutsq); // to be deleted
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSpinNeel::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2)
    error->all(FLERR,"Incorrect number of args in pair_style pair/spin command");

  if (strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"Spin simulations require metal unit style");

  cut_spin_neel_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = i+1; j <= atom->ntypes; j++) {
        if (setflag[i][j]) {
          rs[i][j] = cut_spin_neel_global;
          rl[i][j] = cut_spin_neel_global;
          dr[i][j] = 0.0;
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
  if (narg != 17)
    error->all(FLERR,"Incorrect args in pair_style command");

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  const double r1 = force->numeric(FLERR,arg[3]);
  const double r2 = force->numeric(FLERR,arg[4]);
  const double r3 = force->numeric(FLERR,arg[5]);
  const double ro = force->numeric(FLERR,arg[6]);
  const double k1 = force->numeric(FLERR,arg[7]);
  const double k2 = force->numeric(FLERR,arg[8]);
  const double k3 = force->numeric(FLERR,arg[9]);
  const double k4 = force->numeric(FLERR,arg[10]);
  const double k5 = force->numeric(FLERR,arg[11]);
  const double l1 = force->numeric(FLERR,arg[12]);
  const double l2 = force->numeric(FLERR,arg[13]);
  const double l3 = force->numeric(FLERR,arg[14]);
  const double l4 = force->numeric(FLERR,arg[15]);
  const double l5 = force->numeric(FLERR,arg[16]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      rs[i][j] = r1;
      rl[i][j] = r2;
      dr[i][j] = r3;
      r0[i][j] = ro;

      g1[i][j] = k1/hbar;
      g2[i][j] = k2/hbar;
      g3[i][j] = k3/hbar;
      g4[i][j] = k4/hbar;
      g5[i][j] = k5/hbar;
      gm1[i][j] = k1;
      gm2[i][j] = k2;
      gm3[i][j] = k3;
      gm4[i][j] = k4;
      gm5[i][j] = k5;
      
      q1[i][j] = l1/hbar;
      q2[i][j] = l2/hbar;
      q3[i][j] = l3/hbar;
      q4[i][j] = l4/hbar;
      q5[i][j] = l5/hbar;
      qm1[i][j] = l1;
      qm2[i][j] = l2;
      qm3[i][j] = l3;
      qm4[i][j] = l4;
      qm5[i][j] = l5;
      setflag[i][j] = 1;
      count++;
    }
  }
  if (count == 0)
    error->all(FLERR,"Incorrect args in pair_style command");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSpinNeel::init_style()
{
  if (!atom->sp_flag)
    error->all(FLERR,"Pair spin requires atom/spin style");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  // checking if nve/spin is a listed fix

  int ifix = 0;
  while (ifix < modify->nfix) {
    if (strcmp(modify->fix[ifix]->style,"nve/spin") == 0) break;
    if (strcmp(modify->fix[ifix]->style,"neb/spin") == 0) break;
    ifix++;
  }
  if ((ifix == modify->nfix) && (comm->me == 0))
    error->warning(FLERR,"Using pair/spin style without nve/spin or neb/spin");

  // get the lattice_flag from nve/spin

  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"nve/spin") == 0) {
      lockfixnvespin = (FixNVESpin *) modify->fix[i];
      lattice_flag = lockfixnvespin->lattice_flag;
    }
  }

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSpinNeel::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  rs[j][i] = rs[i][j];
  rl[j][i] = rl[i][j];
  dr[j][i] = dr[i][j];
  r0[j][i] = r0[i][j];

  g1[j][i] = g1[i][j];
  g2[j][i] = g2[i][j];
  g3[j][i] = g3[i][j];
  g4[j][i] = g4[i][j];
  g5[j][i] = g5[i][j];
  gm1[j][i] = gm1[i][j];
  gm2[j][i] = gm2[i][j];
  gm3[j][i] = gm3[i][j];
  gm4[j][i] = gm4[i][j];
  gm5[j][i] = gm5[i][j];

  q1[j][i] = q1[i][j];
  q2[j][i] = q2[i][j];
  q3[j][i] = q3[i][j];
  q4[j][i] = q4[i][j];
  q5[j][i] = q5[i][j];
  qm1[j][i] = qm1[i][j];
  qm2[j][i] = qm2[i][j];
  qm3[j][i] = qm3[i][j];
  qm4[j][i] = qm4[i][j];
  qm5[j][i] = qm5[i][j];

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
  double r_s,r_l,d_r,r1,r2,r_ij;
  double pi,sm;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;
  pi = MY_PI;

  double **x = atom->x;
  double **f = atom->f;
  double *emag = atom->emag;
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

      //itype = type[i];
      jtype = type[j];

      r_s = rs[itype][jtype];
      r_l = rl[itype][jtype];
      d_r = dr[itype][jtype];

      local_cut2 = (r_l + d_r) * (r_l + d_r);
     
      // define smoothing factor
      // compute neel interaction and energy

      if (rsq <= local_cut2) {
        r_ij = sqrt(rsq);
        r1 = fabs(r_s - r_ij);
        r2 = fabs(r_l - r_ij);

        sm = 0.0;
        if (r1 <= d_r) {
          sm = (1.0 - sin(pi*(r_s-r_ij)/(2.0*d_r))) / 2.0;
        } else if (r_ij > (r_s+d_r) && r_ij < (r_l-d_r)) {
          sm = 1.0;
        } else if (r2 <= d_r) {
          sm = (1.0 - sin(pi*(r_ij-r_l)/(2.0*d_r))) / 2.0;
        } else sm = 0.0;

        compute_neel(i,j,r_ij,eij,fmi,spi,spj);
        if (lattice_flag) {
          compute_neel_mech(i,j,r_ij,eij,fi,spi,spj);
        }
        
	if (eflag) {
          evdwl = sm * compute_neel_energy(i,j,r_ij,eij,spi,spj);
          emag[i] += evdwl;
        } else evdwl = 0.0;

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
            evdwl,ecoul,fi[0],fi[1],fi[2],rij[0],rij[1],rij[2]);
      }

      f[i][0] += sm * fi[0];
      f[i][1] += sm * fi[1];
      f[i][2] += sm * fi[2];
      fm[i][0] += sm * fmi[0];
      fm[i][1] += sm * fmi[1];
      fm[i][2] += sm * fmi[2];

      if (newton_pair || j < nlocal) {
        f[j][0] -= fi[0];
        f[j][1] -= fi[1];
        f[j][2] -= fi[2];
      }

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
  double local_cut2,rsq, inorm;
  double r_s,r_l,d_r,r1,r2,r_ij;
  double pi,sm;
  double xi[3],rij[3],eij[3];
  double spi[3],spj[3],fmij[3];

  int i,j,jnum,itype,jtype,ntypes;
  int k,locflag;
  int *jlist,*numneigh,**firstneigh;

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

  // if interaction applies to type of ii,
  // locflag = 1 and compute pair interaction

  if (locflag == 1) {

    spi[0] = sp[ii][0];
    spi[1] = sp[ii][1];
    spi[2] = sp[ii][2];
    
    xi[0] = x[ii][0];
    xi[1] = x[ii][1];
    xi[2] = x[ii][2];
    eij[0] = eij[1] = eij[2] = 0.0;

    pi = MY_PI;
    jlist = firstneigh[ii];
    jnum = numneigh[ii];

    for (int jj = 0; jj < jnum; jj++) {

      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

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
      
      fmij[0] = fmij[1] = fmij[2] = 0.0;

      r_s = rs[itype][jtype];
      r_l = rl[itype][jtype];
      d_r = dr[itype][jtype];

      local_cut2 = (r_l + d_r) * (r_l + d_r);
      
      // define smoothing factor
      // compute neel interaction

      if (rsq <= local_cut2) {
        r_ij = sqrt(rsq);
        r1 = fabs(r_s - r_ij);
        r2 = fabs(r_l - r_ij);

        // to be corrected (error in sw func)

        sm = 0.0;
        if (r1 <= d_r) {
          sm = (1.0 - sin(pi*(r_s-r_ij)/(2.0*d_r))) / 2.0;
        } else if (r_ij > (r_s+d_r) && r_ij < (r_l-d_r)) {
          sm = 1.0;
        } else if (r2 <= d_r) {
          sm = (1.0 - sin(pi*(r_ij-r_l)/(2.0*d_r))) / 2.0;
        } else sm = 0.0;

        compute_neel(i,j,r_ij,eij,fmij,spi,spj);
      }
      fmi[0] += sm * fmij[0]; 
      fmi[1] += sm * fmij[1]; 
      fmi[2] += sm * fmij[2]; 
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairSpinNeel::compute_neel(int i, int j, double rij, double eij[3], 
    double fmi[3], double spi[3], double spj[3])
{
  int *type = atom->type;
  int itype, jtype;
  itype = type[i];
  jtype = type[j];

  double ro,ri3;
  double r1,r2,r3,r4;
  double qr,gr,g1r,q1r,q2r;
  double ga,gb,gc,gd,ge;
  double qa,qb,qc,qd,qe;
  double eij_si,eij_sj,si_sj;
  double coeff0,coeff1,coeff2,coeff3,coeff4;
  double fmix,fmiy,fmiz;

  ro = r0[itype][jtype];

  ga = g1[itype][jtype];
  gb = g2[itype][jtype];
  gc = g3[itype][jtype];
  gd = g4[itype][jtype];
  ge = g5[itype][jtype];

  qa = q1[itype][jtype];
  qb = q2[itype][jtype];
  qc = q3[itype][jtype];
  qd = q4[itype][jtype];
  qe = q5[itype][jtype];

  ri3 = 1.0 / 3.0;
  r1 = (rij-ro);
  r2 = r1*r1;
  r3 = r2*r1;
  r4 = r2*r2;

  gr = ga + gb*r1 + gc*r2 + gd*r3 + ge*r4;
  qr = qa + qb*r1 + qc*r2 + qd*r3 + qe*r4;

  g1r = gr + 12.0*qr/35.0;
  q1r = 9.0*qr/5.0;
  //q1r = 0.0;
  q2r = -2.0*qr/5.0;
  //q2r = 0.0;

  //printf("test g1r: %g, q1r: %g \n",g1r,q1r);
  //printf("test %g %g %g | %g \n",g1r,q1r,q2r,rij);
  
  eij_si = eij[0]*spi[0] + eij[1]*spi[1] + eij[2]*spi[2];
  eij_sj = eij[0]*spj[0] + eij[1]*spj[1] + eij[2]*spj[2];
  si_sj = spi[0]*spj[0] + spi[1]*spj[1] + spi[2]*spj[2];
  
  fmix = g1r * (eij[0]*eij_sj - spj[0]*ri3);
  fmiy = g1r * (eij[1]*eij_sj - spj[1]*ri3);
  fmiz = g1r * (eij[2]*eij_sj - spj[2]*ri3);
 
  coeff0 = q1r * (eij_sj*eij_sj - si_sj*ri3);
  coeff1 = coeff0 * 2.0 * eij_si;
  coeff2 = q1r * (eij_si*eij_si - si_sj*ri3);
  fmix += coeff1*eij[0] - spj[0]*ri3*(coeff0 + coeff2);
  fmiy += coeff1*eij[1] - spj[1]*ri3*(coeff0 + coeff2);
  fmiz += coeff1*eij[2] - spj[2]*ri3*(coeff0 + coeff2);

  coeff3 = q2r*eij_sj*eij_sj*eij_sj;
  coeff4 = 3.0*q2r*eij_sj*eij_si*eij_si;
  fmix += (coeff3 + coeff4) * eij[0];
  fmiy += (coeff3 + coeff4) * eij[1];
  fmiz += (coeff3 + coeff4) * eij[2];

  fmi[0] += fmix;
  fmi[1] += fmiy;
  fmi[2] += fmiz;
}

/* ---------------------------------------------------------------------- */

void PairSpinNeel::compute_neel_mech(int i, int j, double rij, double eij[3], 
    double fi[3], double spi[3], double spj[3])
{
  int *type = atom->type;
  int itype, jtype;
  itype = type[i];
  jtype = type[j];
  
  double ro,ri3;
  double r1,r2,r3,r4;
  double qr,gr,g1r,q1r,q2r;
  double dqr,dgr,dg1r,dq1r,dq2r;
  double ga,gb,gc,gd,ge;
  double qa,qb,qc,qd,qe;
  double dga,dgb,dgc,dgd,dge;
  double dqa,dqb,dqc,dqd,dqe;
  double eij_si,eij_sj,si_sj;
  double coeff1,coeff2,coeff3,coeff4;
  double coeff5,coeff6,coeff7,coeff8;
  double coeff9,coeff10,coeff11,coeff12;
  double fix,fiy,fiz;

  ro = r0[itype][jtype];

  ga = gm1[itype][jtype];
  gb = gm2[itype][jtype];
  gc = gm3[itype][jtype];
  gd = gm4[itype][jtype];
  ge = gm5[itype][jtype];

  qa = qm1[itype][jtype];
  qb = qm2[itype][jtype];
  qc = qm3[itype][jtype];
  qd = qm4[itype][jtype];
  qe = qm5[itype][jtype];

  ri3 = 1.0 / 3.0; 
  r1 = (rij-ro);
  r2 = r1*r1;
  r3 = r2*r1;
  r4 = r2*r2;

  gr = ga + gb*r1 + gc*r2 + gd*r3 + ge*r4;
  qr = qa + qb*r1 + qc*r2 + qd*r3 + qe*r4;
  
  dgr = gb + 2.0*gc*r1 + 3.0*gd*r2 + 4.0*ge*r3;
  dqr = qb + 2.0*qc*r1 + 3.0*qd*r2 + 4.0*qe*r3;
  
  g1r = gr + 12.0*qr/35.0;
  q1r = 9.0*qr/5.0;
  q2r = -2.0*qr/5.0;
  
  dg1r = dgr + 12.0*dqr/35.0;
  dq1r = 9.0*dqr/5.0;
  dq2r = -2.0*dqr/5.0;

  eij_si = eij[0]*spi[0] + eij[1]*spi[1] + eij[2]*spi[2];
  eij_sj = eij[0]*spj[0] + eij[1]*spj[1] + eij[2]*spj[2];
  si_sj = spi[0]*spj[0] + spi[1]*spj[1] + spi[2]*spj[2];

  coeff1 = dg1r * (eij_si*eij_sj - si_sj*ri3);
  coeff2 = dq1r * (eij_si*eij_si - si_sj*ri3) * (eij_sj*eij_sj - si_sj*ri3); 
  coeff3 = dq2r * (eij_sj*eij_si*eij_si*eij_si + eij_si*eij_sj*eij_sj*eij_sj);

  fix = (coeff1 + coeff2 + coeff3) * eij[0];
  fiy = (coeff1 + coeff2 + coeff3) * eij[1];
  fiz = (coeff1 + coeff2 + coeff3) * eij[2];

  coeff4 = g1r * eij_sj / rij;
  coeff5 = g1r * eij_si / rij;
  coeff6 = -2.0 * g1r * eij_si * eij_sj / rij;

  fix += coeff4*spi[0] + coeff5*spj[0] + coeff6*eij[0];
  fiy += coeff4*spi[1] + coeff5*spj[1] + coeff6*eij[1];
  fiz += coeff4*spi[2] + coeff5*spj[2] + coeff6*eij[2];

  coeff7 = q1r*2.0*eij_si*(eij_sj*eij_sj - si_sj*ri3)/rij;
  coeff8 = q1r*2.0*eij_sj*(eij_si*eij_si - si_sj*ri3)/rij;
  coeff9 = q1r*2.0*(2.0*eij_si*eij_si*eij_sj*eij_sj 
      + si_sj*si_sj*ri3*(eij_si*eij_si + eij_sj*eij_sj))/rij;

  fix += coeff7*spi[0] + coeff8*spj[0] - coeff9*eij[0];
  fiy += coeff7*spi[1] + coeff8*spj[1] - coeff9*eij[1];
  fiz += coeff7*spi[2] + coeff8*spj[2] - coeff9*eij[2];

  coeff10 = q2r*eij_sj*(eij_sj*eij_sj + 3.0*eij_si*eij_si)/rij;
  coeff11 = q2r*eij_si*(eij_si*eij_si + 3.0*eij_sj*eij_sj)/rij;
  coeff12 = -q2r*4.0*eij_si*eij_sj*(eij_si*eij_si + eij_sj*eij_sj)/rij;

  fix += coeff10*spi[0] + coeff11*spj[0] + coeff12*eij[0];
  fix += coeff10*spi[1] + coeff11*spj[1] + coeff12*eij[1];
  fix += coeff10*spi[2] + coeff11*spj[2] + coeff12*eij[2];

  fi[0] -= fix;
  fi[1] -= fiy;
  fi[2] -= fiz;
}

/* ---------------------------------------------------------------------- */

double PairSpinNeel::compute_neel_energy(int i, int j, double rij, double eij[3], 
    double spi[3], double spj[3])
{
  int *type = atom->type;
  int itype, jtype;
  itype = type[i];
  jtype = type[j];

  double ro,ri3;
  double r1,r2,r3,r4;
  double gr,qr,g1r,q1r,q2r;
  double ga,gb,gc,gd,ge;
  double qa,qb,qc,qd,qe;
  double eij_si,eij_sj,si_sj;
  double energy = 0.0;

  ro = r0[itype][jtype];

  ga = gm1[itype][jtype];
  gb = gm2[itype][jtype];
  gc = gm3[itype][jtype];
  gd = gm4[itype][jtype];
  ge = gm5[itype][jtype];

  qa = qm1[itype][jtype];
  qb = qm2[itype][jtype];
  qc = qm3[itype][jtype];
  qd = qm4[itype][jtype];
  qe = qm5[itype][jtype];

  ri3 = 1.0 / 3.0;
  r1 = (rij-ro);
  r2 = r1*r1;
  r3 = r2*r1;
  r4 = r2*r2;

  gr = ga + gb*r1 + gc*r2 + gd*r3 + ge*r4;
  qr = qa + qb*r1 + qc*r2 + qd*r3 + qe*r4;

  g1r = gr + 12.0*qr/35.0;
  q1r = 9.0*qr/5.0;
  //q1r = 0.0;
  q2r = -2.0*qr/5.0;
  //q2r = 0.0;

  //printf("test energy g1r: %g, q1r: %g, q2r: %g\n",g1r,q1r,q2r);
  //printf("test gi %g %g %g %g %g \n",ga,gb,gc,gd,ge);
  //printf("test %g %g %g | %g \n",g1r,q1r,q2r,rij);
 
  eij_si = eij[0]*spi[0] + eij[1]*spi[1] + eij[2]*spi[2];
  eij_sj = eij[0]*spj[0] + eij[1]*spj[1] + eij[2]*spj[2];
  si_sj = spi[0]*spj[0] + spi[1]*spj[1] + spi[2]*spj[2];
  
  energy = g1r*(si_sj*ri3 - eij_si*eij_sj);
  //printf("test energy1: %g\n",energy);
  energy -= q1r*(eij_si*eij_si - si_sj*ri3)*(eij_sj*eij_sj - si_sj*ri3);
  //printf("test energy2: %g\n",energy);
  energy -= q2r*eij_sj*eij_si*(eij_si*eij_si + eij_sj*eij_sj);
  //printf("test energy3: %g\n",energy);
  
  return energy;
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

  memory->create(rs,n+1,n+1,"pair/spin/soc/neel:rs");
  memory->create(rl,n+1,n+1,"pair/spin/soc/neel:rl");
  memory->create(dr,n+1,n+1,"pair/spin/soc/neel:dr");
  memory->create(r0,n+1,n+1,"pair/spin/soc/neel:r0");

  memory->create(g1,n+1,n+1,"pair/spin/neel:g1");
  memory->create(g2,n+1,n+1,"pair/spin/neel:g2");
  memory->create(g3,n+1,n+1,"pair/spin/neel:g3");
  memory->create(g4,n+1,n+1,"pair/spin/neel:g4");
  memory->create(g5,n+1,n+1,"pair/spin/neel:g5");
  memory->create(gm1,n+1,n+1,"pair/spin/neel:gm1");
  memory->create(gm2,n+1,n+1,"pair/spin/neel:gm2");
  memory->create(gm3,n+1,n+1,"pair/spin/neel:gm3");
  memory->create(gm4,n+1,n+1,"pair/spin/neel:gm4");
  memory->create(gm5,n+1,n+1,"pair/spin/neel:gm5");

  memory->create(q1,n+1,n+1,"pair/spin/neel:q1");
  memory->create(q2,n+1,n+1,"pair/spin/neel:q2");
  memory->create(q3,n+1,n+1,"pair/spin/neel:q3");
  memory->create(q4,n+1,n+1,"pair/spin/neel:q4");
  memory->create(q5,n+1,n+1,"pair/spin/neel:q5");
  memory->create(qm1,n+1,n+1,"pair/spin/neel:qm1");
  memory->create(qm2,n+1,n+1,"pair/spin/neel:qm2");
  memory->create(qm3,n+1,n+1,"pair/spin/neel:qm3");
  memory->create(qm4,n+1,n+1,"pair/spin/neel:qm4");
  memory->create(qm5,n+1,n+1,"pair/spin/neel:qm5");

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
        fwrite(&q1[i][j],sizeof(double),1,fp);
        fwrite(&q2[i][j],sizeof(double),1,fp);
        fwrite(&q3[i][j],sizeof(double),1,fp);
        fwrite(&q4[i][j],sizeof(double),1,fp);
        fwrite(&q5[i][j],sizeof(double),1,fp);
        fwrite(&qm1[i][j],sizeof(double),1,fp);
        fwrite(&qm2[i][j],sizeof(double),1,fp);
        fwrite(&qm3[i][j],sizeof(double),1,fp);
        fwrite(&qm4[i][j],sizeof(double),1,fp);
        fwrite(&qm5[i][j],sizeof(double),1,fp);

        fwrite(&q1[i][j],sizeof(double),1,fp);
        fwrite(&q2[i][j],sizeof(double),1,fp);
        fwrite(&q3[i][j],sizeof(double),1,fp);
        fwrite(&q4[i][j],sizeof(double),1,fp);
        fwrite(&q5[i][j],sizeof(double),1,fp);
        fwrite(&qm1[i][j],sizeof(double),1,fp);
        fwrite(&qm2[i][j],sizeof(double),1,fp);
        fwrite(&qm3[i][j],sizeof(double),1,fp);
        fwrite(&qm4[i][j],sizeof(double),1,fp);
        fwrite(&qm5[i][j],sizeof(double),1,fp);
        fwrite(&rs[i][j],sizeof(double),1,fp);
        fwrite(&rl[i][j],sizeof(double),1,fp);
        fwrite(&dr[i][j],sizeof(double),1,fp);
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
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&r0[i][j],sizeof(double),1,fp);

          fread(&g1[i][j],sizeof(double),1,fp);
          fread(&g2[i][j],sizeof(double),1,fp);
          fread(&g3[i][j],sizeof(double),1,fp);
          fread(&g4[i][j],sizeof(double),1,fp);
          fread(&g5[i][j],sizeof(double),1,fp);
          fread(&gm1[i][j],sizeof(double),1,fp);
          fread(&gm2[i][j],sizeof(double),1,fp);
          fread(&gm3[i][j],sizeof(double),1,fp);
          fread(&gm4[i][j],sizeof(double),1,fp);
          fread(&gm5[i][j],sizeof(double),1,fp);
          
          fread(&q1[i][j],sizeof(double),1,fp);
          fread(&q2[i][j],sizeof(double),1,fp);
          fread(&q3[i][j],sizeof(double),1,fp);
          fread(&q4[i][j],sizeof(double),1,fp);
          fread(&q5[i][j],sizeof(double),1,fp);
          fread(&qm1[i][j],sizeof(double),1,fp);
          fread(&qm2[i][j],sizeof(double),1,fp);
          fread(&qm3[i][j],sizeof(double),1,fp);
          fread(&qm4[i][j],sizeof(double),1,fp);
          fread(&qm5[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&r0[i][j],1,MPI_DOUBLE,0,world);
        
	MPI_Bcast(&g1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&g2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&g3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&g4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&g5[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gm1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gm2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gm3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gm4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gm5[i][j],1,MPI_DOUBLE,0,world);
        
	MPI_Bcast(&q1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&q2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&q3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&q4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&q5[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&qm1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&qm2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&qm3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&qm4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&qm5[i][j],1,MPI_DOUBLE,0,world);
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
    fread(&cut_spin_neel_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_spin_neel_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}
