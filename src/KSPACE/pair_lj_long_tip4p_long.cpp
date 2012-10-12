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
   Contributing authors: Amalie Frischknecht and Ahmed Ismail (SNL)
                         Rolf Isele-Holder (Aachen University)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_lj_long_tip4p_long.h"
#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "kspace.h"
#include "update.h"
#include "respa.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairLJLongTIP4PLong::PairLJLongTIP4PLong(LAMMPS *lmp) : 
  PairLJLongCoulLong(lmp)
{
  tip4pflag = 1;
  single_enable = 0;
  respa_enable = 0;

  nmax = 0;
  hneigh = NULL;
  newsite = NULL;

  // TIP4P cannot compute virial as F dot r
  // due to find_M() finding bonded H atoms which are not near O atom

  no_virial_fdotr_compute = 1;
}

/* ---------------------------------------------------------------------- */

PairLJLongTIP4PLong::~PairLJLongTIP4PLong()
{
  memory->destroy(hneigh);
  memory->destroy(newsite);
}

/* ---------------------------------------------------------------------- */

void PairLJLongTIP4PLong::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,itable;
  int n,vlist[6];
  int key;
  int iH1,iH2,jH1,jH2;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul;
  double fraction,table;
  double r,r2inv,forcecoul,forcelj,cforce;
  double factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  double xiM[3],xjM[3],fO[3],fH[3],fd[3],v[6],xH1[3],xH2[3];// f1[3];
  double *x1,*x2;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq;
 
  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  // reallocate hneigh & newsite if necessary
  // initialize hneigh[0] to -1 on steps when reneighboring occurred
  // initialize hneigh[2] to 0 every step

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(hneigh);
    memory->create(hneigh,nmax,3,"pair:hneigh");
    memory->destroy(newsite);
    memory->create(newsite,nmax,3,"pair:newsite");
  }
  if (neighbor->ago == 0)
    for (i = 0; i < nall; i++) hneigh[i][0] = -1;
  for (i = 0; i < nall; i++) hneigh[i][2] = 0;

  double **f = atom->f;
  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;
  double cut_coulsqplus = (cut_coul+2.0*qdist)*(cut_coul+2.0*qdist);

  int order1 = ewald_order&(1<<1), order6 = ewald_order&(1<<6);
  int ni;
  double  *cut_ljsqi, *lj1i, *lj2i, *lj3i, *lj4i, *offseti;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2;

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
    if (itype == typeO) {
      if (hneigh[i][0] < 0) {
        hneigh[i][0] = iH1 = atom->map(atom->tag[i] + 1);
        hneigh[i][1] = iH2 = atom->map(atom->tag[i] + 2);
        hneigh[i][2] = 1;
        if (iH1 == -1 || iH2 == -1)
          error->one(FLERR,"TIP4P hydrogen is missing");
        if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
          error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
        compute_newsite(x[i],x[iH1],x[iH2],newsite[i]);
      } else {
        iH1 = hneigh[i][0];
        iH2 = hneigh[i][1];
        if (hneigh[i][2] == 0) {
          hneigh[i][2] = 1;
          compute_newsite(x[i],x[iH1],x[iH2],newsite[i]);
        }
      }
      x1 = newsite[i];
    } else x1 = x[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    offseti = offset[itype];
    lj1i = lj1[itype]; lj2i = lj2[itype]; lj3i = lj3[itype]; lj4i = lj4[itype];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      ni = sbmask(j);
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
        
      if (rsq < cut_ljsq[itype][jtype]) {			// lj
        r2inv = 1.0/rsq;
       	if (order6) {					// long-range lj
	  register double rn = r2inv*r2inv*r2inv;
	  register double x2 = g2*rsq, a2 = 1.0/x2;
	  x2 = a2*exp(-x2)*lj4i[jtype];
	  if (ni == 0) {
	    forcelj =
	      (rn*=rn)*lj1i[jtype]-g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq;
	    if (eflag)
	      evdwl = rn*lj3i[jtype]-g6*((a2+1.0)*a2+0.5)*x2;
	  }
	  else {					// special case
	    register double f = special_lj[ni], t = rn*(1.0-f);
	    forcelj = f*(rn *= rn)*lj1i[jtype]-
	      g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*lj2i[jtype];
	    if (eflag) 
	      evdwl = f*rn*lj3i[jtype]-g6*((a2+1.0)*a2+0.5)*x2+t*lj4i[jtype];
	  }
	}
	else {						// cut lj
	  register double rn = r2inv*r2inv*r2inv;
	  if (ni == 0) {
	    forcelj = rn*(rn*lj1i[jtype]-lj2i[jtype]);
	    if (eflag) evdwl = rn*(rn*lj3i[jtype]-lj4i[jtype])-offseti[jtype];
	  }
	  else {					// special case
	    register double f = special_lj[ni];
	    forcelj = f*rn*(rn*lj1i[jtype]-lj2i[jtype]);
	    if (eflag)
	      evdwl = f * (rn*(rn*lj3i[jtype]-lj4i[jtype])-offseti[jtype]);
	  }
        }

        forcelj *= r2inv;
	f[i][0] += delx*forcelj;
	f[i][1] += dely*forcelj;
	f[i][2] += delz*forcelj;
	f[j][0] -= delx*forcelj;
	f[j][1] -= dely*forcelj;
	f[j][2] -= delz*forcelj;

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
	      	             evdwl,0.0,forcelj,delx,dely,delz);
      }

      
      // adjust rsq and delxyz for off-site O charge(s)
      // ADDITIONAL REQEUST REQUIRED HERE!!!!!

      if (rsq < cut_coulsqplus) {
        if (itype == typeO || jtype == typeO) { 
	  if (jtype == typeO) {
            if (hneigh[j][0] < 0) {
              hneigh[j][0] = jH1 = atom->map(atom->tag[j] + 1);
              hneigh[j][1] = jH2 = atom->map(atom->tag[j] + 2);
              hneigh[j][2] = 1;
              if (jH1 == -1 || jH2 == -1)
                error->one(FLERR,"TIP4P hydrogen is missing");
              if (atom->type[jH1] != typeH || atom->type[jH2] != typeH)
                error->one(FLERR,"TIP4P hydrogen has incorrect atom type");
              compute_newsite(x[j],x[jH1],x[jH2],newsite[j]);
            } else {
              jH1 = hneigh[j][0];
              jH2 = hneigh[j][1];
              if (hneigh[j][2] == 0) {
                hneigh[j][2] = 1;
                compute_newsite(x[j],x[jH1],x[jH2],newsite[j]);
              }
            }
            x2 = newsite[j];
	  } else x2 = x[j];
	  delx = x1[0] - x2[0];
	  dely = x1[1] - x2[1];
	  delz = x1[2] - x2[2];
	  rsq = delx*delx + dely*dely + delz*delz;
        }

	// test current rsq against cutoff and compute Coulombic force
      
        if (rsq < cut_coulsq && order1) {
	  r2inv = 1.0 / rsq;
	  if (!ncoultablebits || rsq <= tabinnersq) {
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

	  //if (evflag) ev_tally(i,j,nlocal,newton_pair,
          //		       evdwl,0.0,cforce,delx,dely,delz);

	  // if i,j are not O atoms, force is applied directly
	  // if i or j are O atoms, force is on fictitious atom & partitioned
	  // force partitioning due to Feenstra, J Comp Chem, 20, 786 (1999)
	  // f_f = fictitious force, fO = f_f (1 - 2 alpha), fH = alpha f_f
	  // preserves total force and torque on water molecule
	  // virial = sum(r x F) where each water's atoms are near xi and xj
	  // vlist stores 2,4,6 atoms whose forces contribute to virial

	  n = 0;
          key = 0;

	  if (itype != typeO) {
	    f[i][0] += delx * cforce;
	    f[i][1] += dely * cforce;
	    f[i][2] += delz * cforce;

            if (vflag) {
              v[0] = x[i][0] * delx * cforce;
              v[1] = x[i][1] * dely * cforce;
              v[2] = x[i][2] * delz * cforce;
              v[3] = x[i][0] * dely * cforce;
              v[4] = x[i][0] * delz * cforce;
              v[5] = x[i][1] * delz * cforce;
            }
          vlist[n++] = i;

	  } else {
            key += 1;
            fd[0] = delx*cforce;
            fd[1] = dely*cforce;
            fd[2] = delz*cforce;

            fO[0] = fd[0]*(1 - alpha);
            fO[1] = fd[1]*(1 - alpha);
            fO[2] = fd[2]*(1 - alpha);

            fH[0] = 0.5 * alpha * fd[0];
            fH[1] = 0.5 * alpha * fd[1];
            fH[2] = 0.5 * alpha * fd[2];

            f[i][0] += fO[0];
            f[i][1] += fO[1];
            f[i][2] += fO[2];

            f[iH1][0] += fH[0];
            f[iH1][1] += fH[1];
            f[iH1][2] += fH[2];

            f[iH2][0] += fH[0];
            f[iH2][1] += fH[1];
            f[iH2][2] += fH[2];

	    if (vflag) {
	      domain->closest_image(x[i],x[iH1],xH1);
	      domain->closest_image(x[i],x[iH2],xH2);

	      v[0] = x[i][0]*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
	      v[1] = x[i][1]*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
	      v[2] = x[i][2]*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
	      v[3] = x[i][0]*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
	      v[4] = x[i][0]*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
	      v[5] = x[i][1]*fO[2] + xH1[1]*fH[2] + xH2[1]*fH[2];
	    }
	    vlist[n++] = i;
	    vlist[n++] = iH1;
	    vlist[n++] = iH2;
  	  }

	  if (jtype != typeO) {
	    f[j][0] -= delx * cforce;
	    f[j][1] -= dely * cforce;
	    f[j][2] -= delz * cforce;

	    if (vflag) {
	      v[0] -= x[j][0] * delx * cforce;
	      v[1] -= x[j][1] * dely * cforce;
	      v[2] -= x[j][2] * delz * cforce;
	      v[3] -= x[j][0] * dely * cforce;
	      v[4] -= x[j][0] * delz * cforce;
	      v[5] -= x[j][1] * delz * cforce;
            }
	    vlist[n++] = j;

	  } else {
            key += 2;

	    fd[0] = -delx*cforce;
	    fd[1] = -dely*cforce;
	    fd[2] = -delz*cforce;

            fO[0] = fd[0]*(1 - alpha);
            fO[1] = fd[1]*(1 - alpha);
            fO[2] = fd[2]*(1 - alpha);

            fH[0] = 0.5 * alpha * fd[0];
            fH[1] = 0.5 * alpha * fd[1];
            fH[2] = 0.5 * alpha * fd[2]; 

	    f[j][0] += fO[0];
	    f[j][1] += fO[1];
	    f[j][2] += fO[2];

	    f[jH1][0] += fH[0];
	    f[jH1][1] += fH[1];
	    f[jH1][2] += fH[2];

	    f[jH2][0] += fH[0];
	    f[jH2][1] += fH[1];
	    f[jH2][2] += fH[2];

	    if (vflag) {
	      domain->closest_image(x[j],x[jH1],xH1);
	      domain->closest_image(x[j],x[jH2],xH2);

	      v[0] += x[j][0]*fO[0] + xH1[0]*fH[0] + xH2[0]*fH[0];
	      v[1] += x[j][1]*fO[1] + xH1[1]*fH[1] + xH2[1]*fH[1];
	      v[2] += x[j][2]*fO[2] + xH1[2]*fH[2] + xH2[2]*fH[2];
	      v[3] += x[j][0]*fO[1] + xH1[0]*fH[1] + xH2[0]*fH[1];
	      v[4] += x[j][0]*fO[2] + xH1[0]*fH[2] + xH2[0]*fH[2];
	      v[5] += x[j][1]*fO[2] + xH1[1]*fH[2] + xH2[1]*fH[2];
            }
      	    vlist[n++] = j;
	    vlist[n++] = jH1;
	    vlist[n++] = jH2;
	  }

	  if (eflag) {
	    if (!ncoultablebits || rsq <= tabinnersq)
	      ecoul = prefactor*erfc;
	    else {
	      table = etable[itable] + fraction*detable[itable];
	      ecoul = qtmp*q[j] * table;
	    }
	    if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;
 
          if (evflag) ev_tally_tip4p(key,vlist,v,ecoul,alpha);
	}
      }
    }
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */
#define PAIR_CUTOFF	"Only one cut-off allowed when requesting all long"
#define PAIR_COUL_CUT	"Coulombic cut not supported in pair_style lj/coul"
#define PAIR_LARGEST	"Using largest cut-off for lj/coul long long"
#define PAIR_MIX	"Mixing forced for lj coefficients"

void PairLJLongTIP4PLong::settings(int narg, char **arg)
{
  if (narg < 8 || narg > 9) error->all(FLERR,"Illegal pair_style command");

  ewald_off = 0;
  ewald_order = 0;
  options(arg, 6);
  options(++arg, 1);
  if (!comm->me && ewald_order&(1<<6)) error->warning(FLERR,PAIR_MIX);
  if (!comm->me && ewald_order==((1<<1)|(1<<6))) error->warning(FLERR,PAIR_LARGEST);
  if (!((ewald_order^ewald_off)&(1<<1))) error->all(FLERR,PAIR_COUL_CUT);
  typeO = force->inumeric(arg[1]);
  typeH = force->inumeric(arg[2]);
  typeB = force->inumeric(arg[3]);
  typeA = force->inumeric(arg[4]);
  qdist = force->numeric(arg[5]);


  cut_lj_global = force->numeric(arg[6]);
  if (narg == 8) cut_coul = cut_lj_global;
  else cut_coul = force->numeric(arg[7]);
 
  
  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j]) cut_lj[i][j] = cut_lj_global;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJLongTIP4PLong::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style lj/coul/tip4p requires atom IDs");
  if (!force->newton_pair) 
    error->all(FLERR,"Pair style lj/coul/tip4p requires newton pair on");
  if (!atom->q_flag)
    error->all(FLERR,"Pair style lj/coul/tip4p requires atom attribute q");
  if (force->bond == NULL)
    error->all(FLERR,"Must use a bond style with TIP4P potential");
  if (force->angle == NULL)
    error->all(FLERR,"Must use an angle style with TIP4P potential");

  PairLJLongCoulLong::init_style();

  // set alpha parameter

  double theta = force->angle->equilibrium_angle(typeA);
  double blen = force->bond->equilibrium_distance(typeB);
  alpha = qdist / (cos(0.5*theta) * blen);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJLongTIP4PLong::init_one(int i, int j)
{
  double cut = PairLJLongCoulLong::init_one(i,j);

  // check that LJ epsilon = 0.0 for water H
  // set LJ cutoff to 0.0 for any interaction involving water H
  // so LJ term isn't calculated in compute()

  if ((i == typeH && epsilon[i][i] != 0.0))
    error->all(FLERR,"Water H epsilon must be 0.0 for "
               "pair style lj/coul/tip4p");

  if (i == typeH || j == typeH)
    cut_ljsq[j][i] = cut_ljsq[i][j] = 0.0;

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJLongTIP4PLong::write_restart_settings(FILE *fp)
{
  fwrite(&typeO,sizeof(int),1,fp);
  fwrite(&typeH,sizeof(int),1,fp);
  fwrite(&typeB,sizeof(int),1,fp);
  fwrite(&typeA,sizeof(int),1,fp);
  fwrite(&qdist,sizeof(double),1,fp);

  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&ncoultablebits,sizeof(int),1,fp);
  fwrite(&tabinner,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJLongTIP4PLong::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&typeO,sizeof(int),1,fp);
    fread(&typeH,sizeof(int),1,fp);
    fread(&typeB,sizeof(int),1,fp);
    fread(&typeA,sizeof(int),1,fp);
    fread(&qdist,sizeof(double),1,fp);

    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&ncoultablebits,sizeof(int),1,fp);
    fread(&tabinner,sizeof(double),1,fp);
  }

  MPI_Bcast(&typeO,1,MPI_INT,0,world);
  MPI_Bcast(&typeH,1,MPI_INT,0,world);
  MPI_Bcast(&typeB,1,MPI_INT,0,world);
  MPI_Bcast(&typeA,1,MPI_INT,0,world);
  MPI_Bcast(&qdist,1,MPI_DOUBLE,0,world);

  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&ncoultablebits,1,MPI_INT,0,world);
  MPI_Bcast(&tabinner,1,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
  compute position xM of fictitious charge site for O atom and 2 H atoms
  return it as xM
------------------------------------------------------------------------- */

void PairLJLongTIP4PLong::compute_newsite(double *xO, double *xH1,
                                             double *xH2, double *xM)
{
  double delx1 = xH1[0] - xO[0];
  double dely1 = xH1[1] - xO[1];
  double delz1 = xH1[2] - xO[2];
  domain->minimum_image(delx1,dely1,delz1);

  double delx2 = xH2[0] - xO[0];
  double dely2 = xH2[1] - xO[1];
  double delz2 = xH2[2] - xO[2];
  domain->minimum_image(delx2,dely2,delz2);

  xM[0] = xO[0] + alpha * 0.5 * (delx1 + delx2);
  xM[1] = xO[1] + alpha * 0.5 * (dely1 + dely2);
  xM[2] = xO[2] + alpha * 0.5 * (delz1 + delz2);
}

/* ---------------------------------------------------------------------- */

void *PairLJLongTIP4PLong::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"qdist") == 0) return (void *) &qdist;
  if (strcmp(str,"typeO") == 0) return (void *) &typeO;
  if (strcmp(str,"typeH") == 0) return (void *) &typeH;
  if (strcmp(str,"typeA") == 0) return (void *) &typeA;
  if (strcmp(str,"typeB") == 0) return (void *) &typeB;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;

  const char *ids[] = {
    "B", "sigma", "epsilon", "ewald_order", "ewald_cut", "cut_coul", 
    "ewald_mix", "cut_LJ", NULL};
  void *ptrs[] = {
    lj4, sigma, epsilon, &ewald_order, &cut_coul, &cut_coul, 
    &mix_flag, &cut_lj_global, NULL};
  int i;

  for (i=0; ids[i]&&strcmp(ids[i], str); ++i);
  if (i <= 2) dim = 2;
  else dim = 0;
  return ptrs[i];
  return NULL;
}

/* ----------------------------------------------------------------------
   memory usage of hneigh
------------------------------------------------------------------------- */

double PairLJLongTIP4PLong::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  bytes += 2 * nmax * sizeof(double);
  return bytes;
}
