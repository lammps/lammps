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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "pair_table_rx.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "modify.h"
#include "fix.h"

using namespace LAMMPS_NS;

enum{NONE,RLINEAR,RSQ,BMP};

#ifdef DBL_EPSILON
  #define MY_EPSILON (10.0*DBL_EPSILON)
#else
  #define MY_EPSILON (10.0*2.220446049250313e-16)
#endif

#define OneFluidValue (-1)
#define isOneFluid(_site_) ( (_site_) == OneFluidValue )

/* ---------------------------------------------------------------------- */

PairTableRX::PairTableRX(LAMMPS *lmp) : PairTable(lmp)
{
  fractionalWeighting = true;
  site1 = NULL;
  site2 = NULL;
}

/* ---------------------------------------------------------------------- */

PairTableRX::~PairTableRX()
{
  delete [] site1;
  delete [] site2;
}

/* ---------------------------------------------------------------------- */

void PairTableRX::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,itable;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,evdwlOld,fpair;
  double rsq,factor_lj,fraction,value,a,b;
  int *ilist,*jlist,*numneigh,**firstneigh;
  Table *tb;

  union_int_float_t rsq_lookup;
  int tlm1 = tablength - 1;

  fraction = 0.0;
  a = 0.0;
  b = 0.0;

  evdwlOld = 0.0;
  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  double mixWtSite1old_i, mixWtSite1old_j;
  double mixWtSite2old_i, mixWtSite2old_j;
  double mixWtSite1_i, mixWtSite1_j;
  double mixWtSite2_i, mixWtSite2_j;
  double *uCG = atom->uCG;
  double *uCGnew = atom->uCGnew;

  double *mixWtSite1old = NULL;
  double *mixWtSite2old = NULL;
  double *mixWtSite1 = NULL;
  double *mixWtSite2 = NULL;

  {
    const int ntotal = atom->nlocal + atom->nghost;
    memory->create(mixWtSite1old, ntotal, "PairTableRx::compute::mixWtSite1old");
    memory->create(mixWtSite2old, ntotal, "PairTableRx::compute::mixWtSite2old");
    memory->create(mixWtSite1, ntotal, "PairTableRx::compute::mixWtSite1");
    memory->create(mixWtSite2, ntotal, "PairTableRx::compute::mixWtSite2");

    for (int i = 0; i < ntotal; ++i)
      getMixingWeights(i, mixWtSite1old[i], mixWtSite2old[i], mixWtSite1[i], mixWtSite2[i]);
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    double uCG_i = 0.0;
    double uCGnew_i = 0.0;
    double fx_i = 0.0, fy_i = 0.0, fz_i = 0.0;

    mixWtSite1old_i = mixWtSite1old[i];
    mixWtSite2old_i = mixWtSite2old[i];
    mixWtSite1_i = mixWtSite1[i];
    mixWtSite2_i = mixWtSite2[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        mixWtSite1old_j = mixWtSite1old[j];
        mixWtSite2old_j = mixWtSite2old[j];
        mixWtSite1_j = mixWtSite1[j];
        mixWtSite2_j = mixWtSite2[j];

        tb = &tables[tabindex[itype][jtype]];
        if (rsq < tb->innersq)
          error->one(FLERR,"Pair distance < table inner cutoff");

        if (tabstyle == LOOKUP) {
          itable = static_cast<int> ((rsq - tb->innersq) * tb->invdelta);
          if (itable >= tlm1)
            error->one(FLERR,"Pair distance > table outer cutoff");
          fpair = factor_lj * tb->f[itable];
        } else if (tabstyle == LINEAR) {
          itable = static_cast<int> ((rsq - tb->innersq) * tb->invdelta);
          if (itable >= tlm1)
            error->one(FLERR,"Pair distance > table outer cutoff");
          fraction = (rsq - tb->rsq[itable]) * tb->invdelta;
          value = tb->f[itable] + fraction*tb->df[itable];
          fpair = factor_lj * value;
        } else if (tabstyle == SPLINE) {
          itable = static_cast<int> ((rsq - tb->innersq) * tb->invdelta);
          if (itable >= tlm1)
            error->one(FLERR,"Pair distance > table outer cutoff");
          b = (rsq - tb->rsq[itable]) * tb->invdelta;
          a = 1.0 - b;
          value = a * tb->f[itable] + b * tb->f[itable+1] +
            ((a*a*a-a)*tb->f2[itable] + (b*b*b-b)*tb->f2[itable+1]) *
            tb->deltasq6;
          fpair = factor_lj * value;
        } else {
          rsq_lookup.f = rsq;
          itable = rsq_lookup.i & tb->nmask;
          itable >>= tb->nshiftbits;
          fraction = (rsq_lookup.f - tb->rsq[itable]) * tb->drsq[itable];
          value = tb->f[itable] + fraction*tb->df[itable];
          fpair = factor_lj * value;
        }
        if (isite1 == isite2) fpair = sqrt(mixWtSite1old_i*mixWtSite2old_j)*fpair;
        else fpair = (sqrt(mixWtSite1old_i*mixWtSite2old_j) + sqrt(mixWtSite2old_i*mixWtSite1old_j))*fpair;

        fx_i += delx*fpair;
        fy_i += dely*fpair;
        fz_i += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (tabstyle == LOOKUP)
          evdwl = tb->e[itable];
        else if (tabstyle == LINEAR || tabstyle == BITMAP){
          evdwl = tb->e[itable] + fraction*tb->de[itable];
        }
        else
          evdwl = a * tb->e[itable] + b * tb->e[itable+1] +
            ((a*a*a-a)*tb->e2[itable] + (b*b*b-b)*tb->e2[itable+1]) *
            tb->deltasq6;
        if (isite1 == isite2){
          evdwlOld = sqrt(mixWtSite1old_i*mixWtSite2old_j)*evdwl;
          evdwl = sqrt(mixWtSite1_i*mixWtSite2_j)*evdwl;
        } else {
          evdwlOld = (sqrt(mixWtSite1old_i*mixWtSite2old_j) + sqrt(mixWtSite2old_i*mixWtSite1old_j))*evdwl;
          evdwl = (sqrt(mixWtSite1_i*mixWtSite2_j) + sqrt(mixWtSite2_i*mixWtSite1_j))*evdwl;
        }
        evdwlOld *= factor_lj;
        evdwl *= factor_lj;

        uCG_i += 0.5*evdwlOld;
        uCG[j] += 0.5*evdwlOld;

        uCGnew_i += 0.5*evdwl;
        uCGnew[j] += 0.5*evdwl;
        evdwl = evdwlOld;

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }

    uCG[i] += uCG_i;
    uCGnew[i] += uCGnew_i;

    f[i][0] += fx_i;
    f[i][1] += fy_i;
    f[i][2] += fz_i;
  }
  if (vflag_fdotr) virial_fdotr_compute();

  memory->destroy(mixWtSite1old);
  memory->destroy(mixWtSite2old);
  memory->destroy(mixWtSite1);
  memory->destroy(mixWtSite2);
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairTableRX::settings(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal pair_style command");

  // new settings

  if (strcmp(arg[0],"lookup") == 0) tabstyle = LOOKUP;
  else if (strcmp(arg[0],"linear") == 0) tabstyle = LINEAR;
  else if (strcmp(arg[0],"spline") == 0) tabstyle = SPLINE;
  else if (strcmp(arg[0],"bitmap") == 0) tabstyle = BITMAP;
  else error->all(FLERR,"Unknown table style in pair_style command");

  tablength = force->inumeric(FLERR,arg[1]);
  if (tablength < 2) error->all(FLERR,"Illegal number of pair table entries");

  // optional keywords
  // assert the tabulation is compatible with a specific long-range solver

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"ewald") == 0) ewaldflag = 1;
    else if (strcmp(arg[iarg],"pppm") == 0) pppmflag = 1;
    else if (strcmp(arg[iarg],"msm") == 0) msmflag = 1;
    else if (strcmp(arg[iarg],"dispersion") == 0) dispersionflag = 1;
    else if (strcmp(arg[iarg],"tip4p") == 0) tip4pflag = 1;
    else if (strcmp(arg[iarg],"fractional") == 0) fractionalWeighting = true;
    else if (strcmp(arg[iarg],"molecular") == 0) fractionalWeighting = false;
    else error->all(FLERR,"Illegal pair_style command");
    iarg++;
  }

  // delete old tables, since cannot just change settings

  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(tabindex);
  }
  allocated = 0;

  ntables = 0;
  tables = NULL;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairTableRX::coeff(int narg, char **arg)
{
  if (narg != 6 && narg != 7) error->all(FLERR,"Illegal pair_coeff command");
  if (!allocated) allocate();

  bool rx_flag = false;
  for (int i = 0; i < modify->nfix; i++)
    if (strncmp(modify->fix[i]->style,"rx",2) == 0) rx_flag = true;
  if (!rx_flag) error->all(FLERR,"PairTableRX requires a fix rx command.");

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  int me;
  MPI_Comm_rank(world,&me);
  tables = (Table *)
    memory->srealloc(tables,(ntables+1)*sizeof(Table),"pair:tables");
  Table *tb = &tables[ntables];
  null_table(tb);
  if (me == 0) read_table(tb,arg[2],arg[3]);
  bcast_table(tb);

  nspecies = atom->nspecies_dpd;
  if(nspecies==0) error->all(FLERR,"There are no rx species specified.");
  int n;
  n = strlen(arg[4]) + 1;
  site1 = new char[n];
  strcpy(site1,arg[4]);

  int ispecies;
  for (ispecies = 0; ispecies < nspecies; ispecies++){
    if (strcmp(site1,&atom->dname[ispecies][0]) == 0) break;
  }
  if (ispecies == nspecies && strcmp(site1,"1fluid") != 0)
    error->all(FLERR,"Site1 name not recognized in pair coefficients");

  n = strlen(arg[5]) + 1;
  site2 = new char[n];
  strcpy(site2,arg[5]);

  for (ispecies = 0; ispecies < nspecies; ispecies++){
    if (strcmp(site2,&atom->dname[ispecies][0]) == 0) break;
  }
  if (ispecies == nspecies && strcmp(site2,"1fluid") != 0)
    error->all(FLERR,"Site2 name not recognized in pair coefficients");

  // set table cutoff

  if (narg == 7) tb->cut = force->numeric(FLERR,arg[6]);
  else if (tb->rflag) tb->cut = tb->rhi;
  else tb->cut = tb->rfile[tb->ninput-1];

  // error check on table parameters
  // insure cutoff is within table
  // for BITMAP tables, file values can be in non-ascending order

  if (tb->ninput <= 1) error->one(FLERR,"Invalid pair table length");
  double rlo,rhi;
  if (tb->rflag == 0) {
    rlo = tb->rfile[0];
    rhi = tb->rfile[tb->ninput-1];
  } else {
    rlo = tb->rlo;
    rhi = tb->rhi;
  }
  if (tb->cut <= rlo || tb->cut > rhi)
    error->all(FLERR,"Invalid pair table cutoff");
  if (rlo <= 0.0) error->all(FLERR,"Invalid pair table cutoff");

  // match = 1 if don't need to spline read-in tables
  // this is only the case if r values needed by final tables
  //   exactly match r values read from file
  // for tabstyle SPLINE, always need to build spline tables

  tb->match = 0;
  if (tabstyle == LINEAR && tb->ninput == tablength &&
      tb->rflag == RSQ && tb->rhi == tb->cut) tb->match = 1;
  if (tabstyle == BITMAP && tb->ninput == 1 << tablength &&
      tb->rflag == BMP && tb->rhi == tb->cut) tb->match = 1;
  if (tb->rflag == BMP && tb->match == 0)
    error->all(FLERR,"Bitmapped table in file does not match requested table");

  // spline read-in values and compute r,e,f vectors within table

  if (tb->match == 0) spline_table(tb);
  compute_table(tb);

  // store ptr to table in tabindex

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      tabindex[i][j] = ntables;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Illegal pair_coeff command");
  ntables++;

  {
     if ( strcmp(site1,"1fluid") == 0 )
       isite1 = OneFluidValue;
     else {
       isite1 = nspecies;

       for (int k = 0; k < nspecies; k++){
         if (strcmp(site1, atom->dname[k]) == 0){
           isite1 = k;
           break;
         }
       }

       if (isite1 == nspecies) error->all(FLERR,"isite1 == nspecies");
     }

     if ( strcmp(site2,"1fluid") == 0 )
       isite2 = OneFluidValue;
     else {
       isite2 = nspecies;

       for (int k = 0; k < nspecies; k++){
         if (strcmp(site2, atom->dname[k]) == 0){
           isite2 = ispecies;
           break;
         }
       }

       if (isite2 == nspecies)
         error->all(FLERR,"isite2 == nspecies");
     }
  }

}

/* ---------------------------------------------------------------------- */

double PairTableRX::single(int i, int j, int itype, int jtype, double rsq,
                         double /*factor_coul*/, double factor_lj,
                         double &fforce)
{
  int itable;
  double fraction,value,a,b,phi;
  int tlm1 = tablength - 1;

  Table *tb = &tables[tabindex[itype][jtype]];
  double mixWtSite1_i, mixWtSite1_j;
  double mixWtSite2_i, mixWtSite2_j;
  double mixWtSite1old_i, mixWtSite1old_j;
  double mixWtSite2old_i, mixWtSite2old_j;

  fraction = 0.0;
  a = 0.0;
  b = 0.0;

  getMixingWeights(i,mixWtSite1old_i,mixWtSite2old_i,mixWtSite1_i,mixWtSite2_i);
  getMixingWeights(j,mixWtSite1old_j,mixWtSite2old_j,mixWtSite1_j,mixWtSite2_j);

  if (rsq < tb->innersq) error->one(FLERR,"Pair distance < table inner cutoff");

  if (tabstyle == LOOKUP) {
    itable = static_cast<int> ((rsq-tb->innersq) * tb->invdelta);
    if (itable >= tlm1) error->one(FLERR,"Pair distance > table outer cutoff");
    fforce = factor_lj * tb->f[itable];
  } else if (tabstyle == LINEAR) {
    itable = static_cast<int> ((rsq-tb->innersq) * tb->invdelta);
    if (itable >= tlm1) error->one(FLERR,"Pair distance > table outer cutoff");
    fraction = (rsq - tb->rsq[itable]) * tb->invdelta;
    value = tb->f[itable] + fraction*tb->df[itable];
    fforce = factor_lj * value;
  } else if (tabstyle == SPLINE) {
    itable = static_cast<int> ((rsq-tb->innersq) * tb->invdelta);
    if (itable >= tlm1) error->one(FLERR,"Pair distance > table outer cutoff");
    b = (rsq - tb->rsq[itable]) * tb->invdelta;
    a = 1.0 - b;
    value = a * tb->f[itable] + b * tb->f[itable+1] +
      ((a*a*a-a)*tb->f2[itable] + (b*b*b-b)*tb->f2[itable+1]) *
      tb->deltasq6;
    fforce = factor_lj * value;
  } else {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    itable = rsq_lookup.i & tb->nmask;
    itable >>= tb->nshiftbits;
    fraction = (rsq_lookup.f - tb->rsq[itable]) * tb->drsq[itable];
    value = tb->f[itable] + fraction*tb->df[itable];
    fforce = factor_lj * value;
  }

  if (isite1 == isite2) fforce = sqrt(mixWtSite1_i*mixWtSite2_j)*fforce;
  else fforce = (sqrt(mixWtSite1_i*mixWtSite2_j) + sqrt(mixWtSite2_i*mixWtSite1_j))*fforce;

  if (tabstyle == LOOKUP)
    phi = tb->e[itable];
  else if (tabstyle == LINEAR || tabstyle == BITMAP)
    phi = tb->e[itable] + fraction*tb->de[itable];
  else
    phi = a * tb->e[itable] + b * tb->e[itable+1] +
      ((a*a*a-a)*tb->e2[itable] + (b*b*b-b)*tb->e2[itable+1]) * tb->deltasq6;

  if (isite1 == isite2) phi = sqrt(mixWtSite1_i*mixWtSite2_j)*phi;
  else phi = (sqrt(mixWtSite1_i*mixWtSite2_j) + sqrt(mixWtSite2_i*mixWtSite1_j))*phi;

  return factor_lj*phi;
}

/* ---------------------------------------------------------------------- */

void PairTableRX::getMixingWeights(int id, double &mixWtSite1old, double &mixWtSite2old, double &mixWtSite1, double &mixWtSite2)
{
  double fractionOFAold, fractionOFA;
  double fractionOld1, fraction1;
  double fractionOld2, fraction2;
  double nMoleculesOFAold, nMoleculesOFA;
  double nMoleculesOld1, nMolecules1;
  double nMoleculesOld2, nMolecules2;
  double nTotal, nTotalOld;

  nTotal = 0.0;
  nTotalOld = 0.0;
  for (int ispecies = 0; ispecies < nspecies; ++ispecies){
    nTotal += atom->dvector[ispecies][id];
    nTotalOld += atom->dvector[ispecies+nspecies][id];
  }
  if(nTotal < MY_EPSILON || nTotalOld < MY_EPSILON)
    error->all(FLERR,"The number of molecules in CG particle is less than 10*DBL_EPSILON.");

  if (isOneFluid(isite1) == false){
    nMoleculesOld1 = atom->dvector[isite1+nspecies][id];
    nMolecules1 = atom->dvector[isite1][id];
    fractionOld1 = nMoleculesOld1/nTotalOld;
    fraction1 = nMolecules1/nTotal;
  }
  if (isOneFluid(isite2) == false){
    nMoleculesOld2 = atom->dvector[isite2+nspecies][id];
    nMolecules2 = atom->dvector[isite2][id];
    fractionOld2 = nMoleculesOld2/nTotalOld;
    fraction2 = nMolecules2/nTotal;
  }

  if (isOneFluid(isite1) || isOneFluid(isite2)){
    nMoleculesOFAold  = 0.0;
    nMoleculesOFA  = 0.0;
    fractionOFAold  = 0.0;
    fractionOFA  = 0.0;

    for (int ispecies = 0; ispecies < nspecies; ispecies++){
      if (isite1 == ispecies || isite2 == ispecies) continue;
      nMoleculesOFAold += atom->dvector[ispecies+nspecies][id];
      nMoleculesOFA += atom->dvector[ispecies][id];
      fractionOFAold += atom->dvector[ispecies+nspecies][id]/nTotalOld;
      fractionOFA += atom->dvector[ispecies][id]/nTotal;
    }
    if(isOneFluid(isite1)){
      nMoleculesOld1 = 1.0-(nTotalOld-nMoleculesOFAold);
      nMolecules1 = 1.0-(nTotal-nMoleculesOFA);
      fractionOld1 = fractionOFAold;
      fraction1 = fractionOFA;
    }
    if(isOneFluid(isite2)){
      nMoleculesOld2 = 1.0-(nTotalOld-nMoleculesOFAold);
      nMolecules2 = 1.0-(nTotal-nMoleculesOFA);
      fractionOld2 = fractionOFAold;
      fraction2 = fractionOFA;
    }
  }

  if(fractionalWeighting){
    mixWtSite1old = fractionOld1;
    mixWtSite1 = fraction1;
    mixWtSite2old = fractionOld2;
    mixWtSite2 = fraction2;
  } else {
    mixWtSite1old = nMoleculesOld1;
    mixWtSite1 = nMolecules1;
    mixWtSite2old = nMoleculesOld2;
    mixWtSite2 = nMolecules2;
  }
}
