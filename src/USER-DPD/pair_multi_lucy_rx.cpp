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

/* ----------------------------------------------------------------------------------------
   Contributing authors:
   James Larentzos and Joshua Moore (U.S. Army Research Laboratory)

   Please cite the related publications:
   J.D. Moore, B.C. Barnes, S. Izvekov, M. Lisal, M.S. Sellers, D.E. Taylor & J.K. Brennan
   "A coarse-grain force field for RDX: Density dependent and energy conserving"
   The Journal of Chemical Physics, 2016, 144, 104501.
------------------------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include "math_const.h"
#include <cstdlib>
#include <cstring>
#include "pair_multi_lucy_rx.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "utils.h"
#include "citeme.h"
#include "modify.h"
#include "fix.h"
#include "utils.h"

using namespace LAMMPS_NS;

enum{NONE,RLINEAR,RSQ};

#define MAXLINE 1024

#ifdef DBL_EPSILON
  #define MY_EPSILON (10.0*DBL_EPSILON)
#else
  #define MY_EPSILON (10.0*2.220446049250313e-16)
#endif

#define oneFluidParameter (-1)
#define isOneFluid(_site) ( (_site) == oneFluidParameter )

static const char cite_pair_multi_lucy_rx[] =
  "pair_style multi/lucy/rx command:\n\n"
  "@Article{Moore16,\n"
  " author = {J.D. Moore, B.C. Barnes, S. Izvekov, M. Lisal, M.S. Sellers, D.E. Taylor and J. K. Brennan},\n"
  " title = {A coarse-grain force field for RDX:  Density dependent and energy conserving},\n"
  " journal = {J. Chem. Phys.},\n"
  " year =    2016,\n"
  " volume =  144\n"
  " pages =   {104501}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

PairMultiLucyRX::PairMultiLucyRX(LAMMPS *lmp) : Pair(lmp),
  ntables(0), tables(NULL), tabindex(NULL), site1(NULL), site2(NULL)
{
  if (lmp->citeme) lmp->citeme->add(cite_pair_multi_lucy_rx);

  if (atom->rho_flag != 1) error->all(FLERR,"Pair multi/lucy/rx command requires atom_style with density (e.g. dpd, meso)");

  ntables = 0;
  tables = NULL;

  comm_forward = 1;
  comm_reverse = 1;

  fractionalWeighting = true;
}

/* ---------------------------------------------------------------------- */

PairMultiLucyRX::~PairMultiLucyRX()
{
  if (copymode) return;

  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(tabindex);
  }
}

/* ---------------------------------------------------------------------- */

void PairMultiLucyRX::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,itable;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,evdwlOld,fpair;
  double rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  Table *tb;

  int tlm1 = tablength - 1;

  evdwlOld = 0.0;
  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int newton_pair = force->newton_pair;

  double mixWtSite1old_i,mixWtSite1old_j;
  double mixWtSite2old_i,mixWtSite2old_j;
  double mixWtSite1_i;
  double *uCG = atom->uCG;
  double *uCGnew = atom->uCGnew;

  double pi = MathConst::MY_PI;
  double A_i, A_j;
  double fraction_i,fraction_j;
  int jtable;
  double *rho = atom->rho;

  double *mixWtSite1old = NULL;
  double *mixWtSite2old = NULL;
  double *mixWtSite1 = NULL;
  double *mixWtSite2 = NULL;

  {
    const int ntotal = nlocal + nghost;
    memory->create(mixWtSite1old, ntotal, "PairMultiLucyRX::mixWtSite1old");
    memory->create(mixWtSite2old, ntotal, "PairMultiLucyRX::mixWtSite2old");
    memory->create(mixWtSite1, ntotal, "PairMultiLucyRX::mixWtSite1");
    memory->create(mixWtSite2, ntotal, "PairMultiLucyRX::mixWtSite2");

    for (int i = 0; i < ntotal; ++i)
       getMixingWeights(i, mixWtSite1old[i], mixWtSite2old[i], mixWtSite1[i], mixWtSite2[i]);
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  computeLocalDensity();

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    double fx_i = 0.0;
    double fy_i = 0.0;
    double fz_i = 0.0;

    mixWtSite1old_i = mixWtSite1old[i];
    mixWtSite2old_i = mixWtSite2old[i];
    mixWtSite1_i = mixWtSite1[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        fpair = 0.0;

        mixWtSite1old_j = mixWtSite1old[j];
        mixWtSite2old_j = mixWtSite2old[j];

        tb = &tables[tabindex[itype][jtype]];
        if (rho[i]*rho[i] < tb->innersq || rho[j]*rho[j] < tb->innersq){
          printf("Table inner cutoff = %lf\n",sqrt(tb->innersq));
          printf("rho[%d]=%lf\n",i,rho[i]);
          printf("rho[%d]=%lf\n",j,rho[j]);
          error->one(FLERR,"Density < table inner cutoff");
        }
        if (tabstyle == LOOKUP) {
          itable = static_cast<int> (((rho[i]*rho[i]) - tb->innersq) * tb->invdelta);
          jtable = static_cast<int> (((rho[j]*rho[j]) - tb->innersq) * tb->invdelta);
          if (itable >= tlm1 || jtable >= tlm1){
            printf("Table outer index = %d\n",tlm1);
            printf("itableIndex=%d rho[%d]=%lf\n",itable,i,rho[i]);
            printf("jtableIndex=%d rho[%d]=%lf\n",jtable,j,rho[j]);
            error->one(FLERR,"Density > table outer cutoff");
          }
          A_i = tb->f[itable];
          A_j = tb->f[jtable];

          const double rfactor = 1.0-sqrt(rsq/cutsq[itype][jtype]);
          fpair = 0.5*(A_i + A_j)*(4.0-3.0*rfactor)*rfactor*rfactor*rfactor;
          fpair /= sqrt(rsq);

        } else if (tabstyle == LINEAR) {
          itable = static_cast<int> ((rho[i]*rho[i] - tb->innersq) * tb->invdelta);
          jtable = static_cast<int> (((rho[j]*rho[j]) - tb->innersq) * tb->invdelta);
          if (itable >= tlm1 || jtable >= tlm1){
            printf("Table outer index = %d\n",tlm1);
            printf("itableIndex=%d rho[%d]=%lf\n",itable,i,rho[i]);
            printf("jtableIndex=%d rho[%d]=%lf\n",jtable,j,rho[j]);
            error->one(FLERR,"Density > table outer cutoff");
          }
          if(itable<0) itable=0;
          if(itable>=tlm1) itable=tlm1;
          if(jtable<0) jtable=0;
          if(jtable>=tlm1)jtable=tlm1;

          fraction_i = (((rho[i]*rho[i]) - tb->rsq[itable]) * tb->invdelta);
          fraction_j = (((rho[j]*rho[j]) - tb->rsq[jtable]) * tb->invdelta);
          if(itable==0) fraction_i=0.0;
          if(itable==tlm1) fraction_i=0.0;
          if(jtable==0) fraction_j=0.0;
          if(jtable==tlm1) fraction_j=0.0;

          A_i = tb->f[itable] + fraction_i*tb->df[itable];
          A_j = tb->f[jtable] + fraction_j*tb->df[jtable];

          const double rfactor = 1.0-sqrt(rsq/cutsq[itype][jtype]);
          fpair = 0.5*(A_i + A_j)*(4.0-3.0*rfactor)*rfactor*rfactor*rfactor;
          fpair /= sqrt(rsq);

        } else error->one(FLERR,"Only LOOKUP and LINEAR table styles have been implemented for pair multi/lucy/rx");

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
        if (evflag) ev_tally(i,j,nlocal,newton_pair,0.0,0.0,fpair,delx,dely,delz);
      }
    }

    f[i][0] += fx_i;
    f[i][1] += fy_i;
    f[i][2] += fz_i;

    tb = &tables[tabindex[itype][itype]];
    itable = static_cast<int> (((rho[i]*rho[i]) - tb->innersq) * tb->invdelta);
    if (tabstyle == LOOKUP) evdwl = tb->e[itable];
    else if (tabstyle == LINEAR){
      if (itable >= tlm1){
        printf("itableIndex=%d rho[%d]=%lf\n",itable,i,rho[i]);
        error->one(FLERR,"Density > table outer cutoff");
      }
      if(itable==0) fraction_i=0.0;
      else fraction_i = (((rho[i]*rho[i]) - tb->rsq[itable]) * tb->invdelta);
      evdwl = tb->e[itable] + fraction_i*tb->de[itable];
    } else error->one(FLERR,"Only LOOKUP and LINEAR table styles have been implemented for pair multi/lucy/rx");

    evdwl *=(pi*cutsq[itype][itype]*cutsq[itype][itype])/84.0;
    evdwlOld = mixWtSite1old_i*evdwl;
    evdwl = mixWtSite1_i*evdwl;

    uCG[i] += evdwlOld;
    uCGnew[i] += evdwl;

    evdwl = evdwlOld;

    if (evflag) ev_tally(0,0,nlocal,newton_pair,evdwl,0.0,0.0,0.0,0.0,0.0);
  }

  if (vflag_fdotr) virial_fdotr_compute();

  memory->destroy(mixWtSite1old);
  memory->destroy(mixWtSite2old);
  memory->destroy(mixWtSite1);
  memory->destroy(mixWtSite2);
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMultiLucyRX::allocate()
{
  allocated = 1;
  const int nt = atom->ntypes + 1;

  memory->create(setflag,nt,nt,"pair:setflag");
  memory->create(cutsq,nt,nt,"pair:cutsq");
  memory->create(tabindex,nt,nt,"pair:tabindex");

  memset(&setflag[0][0],0,nt*nt*sizeof(int));
  memset(&cutsq[0][0],0,nt*nt*sizeof(double));
  memset(&tabindex[0][0],0,nt*nt*sizeof(int));
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMultiLucyRX::settings(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal pair_style command");

  // new settings

  if (strcmp(arg[0],"lookup") == 0) tabstyle = LOOKUP;
  else if (strcmp(arg[0],"linear") == 0) tabstyle = LINEAR;
  else error->all(FLERR,"Unknown table style in pair_style command");

  tablength = force->inumeric(FLERR,arg[1]);
  if (tablength < 2) error->all(FLERR,"Illegal number of pair table entries");

  // optional keywords

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"fractional") == 0)   fractionalWeighting = true;
    else if (strcmp(arg[iarg],"molecular") == 0)   fractionalWeighting = false;
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

void PairMultiLucyRX::coeff(int narg, char **arg)
{
  if (narg != 6 && narg != 7) error->all(FLERR,"Illegal pair_coeff command");

  bool rx_flag = false;
  for (int i = 0; i < modify->nfix; i++)
    if (strncmp(modify->fix[i]->style,"rx",2) == 0) rx_flag = true;
  if (!rx_flag) error->all(FLERR,"PairMultiLucyRX requires a fix rx command.");

  if (!allocated) allocate();

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
  int n;
  n = strlen(arg[4]) + 1;
  site1 = new char[n];
  strcpy(site1,arg[4]);

  n = strlen(arg[5]) + 1;
  site2 = new char[n];
  strcpy(site2,arg[5]);

  // set table cutoff

  if (narg == 7) tb->cut = force->numeric(FLERR,arg[6]);
  else if (tb->rflag) tb->cut = tb->rhi;
  else tb->cut = tb->rfile[tb->ninput-1];

  // error check on table parameters
  // insure cutoff is within table

  if (tb->ninput <= 1) error->one(FLERR,"Invalid pair table length");
  if (tb->rflag == 0) {
    rho_0 = tb->rfile[0];
  } else {
    rho_0 = tb->rlo;
  }

  tb->match = 0;
  if (tabstyle == LINEAR && tb->ninput == tablength &&
      tb->rflag == RSQ) tb->match = 1;

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

  // Match site* to isite values.

  if (strcmp(site1, "1fluid") == 0)
     isite1 = oneFluidParameter;
  else {
     isite1 = nspecies;
     for (int ispecies = 0; ispecies < nspecies; ++ispecies)
        if (strcmp(site1, atom->dname[ispecies]) == 0){
           isite1 = ispecies;
           break;
        }

     if (isite1 == nspecies)
        error->all(FLERR,"Pair_multi_lucy_rx site1 is invalid.");
  }

  if (strcmp(site2, "1fluid") == 0)
     isite2 = oneFluidParameter;
  else {
     isite2 = nspecies;
     for (int ispecies = 0; ispecies < nspecies; ++ispecies)
        if (strcmp(site2, atom->dname[ispecies]) == 0){
           isite2 = ispecies;
           break;
        }

     if (isite2 == nspecies)
        error->all(FLERR,"Pair_multi_lucy_rx site2 is invalid.");
  }

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMultiLucyRX::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  tabindex[j][i] = tabindex[i][j];

  return tables[tabindex[i][j]].cut;
}

/* ----------------------------------------------------------------------
   read a table section from a tabulated potential file
   only called by proc 0
   this function sets these values in Table:
     ninput,rfile,efile,ffile,rflag,rlo,rhi,fpflag,fplo,fphi
------------------------------------------------------------------------- */

void PairMultiLucyRX::read_table(Table *tb, char *file, char *keyword)
{
  char line[MAXLINE];
  char *r_token;

  // open file

  FILE *fp = force->open_potential(file);
  if (fp == NULL) {
    char str[128];
    snprintf(str,128,"Cannot open file %s",file);
    error->one(FLERR,str);
  }

  // loop until section found with matching keyword

  while (1) {
    if (fgets(line,MAXLINE,fp) == NULL)
      error->one(FLERR,"Did not find keyword in table file");
    if (strspn(line," \t\n\r") == strlen(line)) continue;  // blank line
    if (line[0] == '#') continue;                          // comment
    r_token = line;
    char *word = utils::strtok_r(r_token," \t\n\r",&r_token);
    if (strcmp(word,keyword) == 0) break;           // matching keyword
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);                         // no match, skip section
    param_extract(tb,line);
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
    for (int i = 0; i < tb->ninput; i++) utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
  }

  // read args on 2nd line of section
  // allocate table arrays for file values

  utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
  param_extract(tb,line);
  memory->create(tb->rfile,tb->ninput,"pair:rfile");
  memory->create(tb->efile,tb->ninput,"pair:efile");
  memory->create(tb->ffile,tb->ninput,"pair:ffile");

  // read r,e,f table values from file
  // if rflag set, compute r
  // if rflag not set, use r from file

  int itmp;
  double rtmp;

  utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
  for (int i = 0; i < tb->ninput; i++) {
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
    sscanf(line,"%d %lg %lg %lg",&itmp,&rtmp,&tb->efile[i],&tb->ffile[i]);

    if (tb->rflag == RLINEAR)
      rtmp = tb->rlo + (tb->rhi - tb->rlo)*i/(tb->ninput-1);
    else if (tb->rflag == RSQ) {
      rtmp = tb->rlo*tb->rlo +
        (tb->rhi*tb->rhi - tb->rlo*tb->rlo)*i/(tb->ninput-1);
      rtmp = sqrt(rtmp);
    }

    tb->rfile[i] = rtmp;
  }

  // close file

  fclose(fp);
}

/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,rfile,efile,ffile,rflag,rlo,rhi,fpflag,fplo,fphi
------------------------------------------------------------------------- */

void PairMultiLucyRX::bcast_table(Table *tb)
{
  MPI_Bcast(&tb->ninput,1,MPI_INT,0,world);

  int me;
  MPI_Comm_rank(world,&me);
  if (me > 0) {
    memory->create(tb->rfile,tb->ninput,"pair:rfile");
    memory->create(tb->efile,tb->ninput,"pair:efile");
    memory->create(tb->ffile,tb->ninput,"pair:ffile");
  }

  MPI_Bcast(tb->rfile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->efile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->ffile,tb->ninput,MPI_DOUBLE,0,world);

  MPI_Bcast(&tb->rflag,1,MPI_INT,0,world);
  if (tb->rflag) {
    MPI_Bcast(&tb->rlo,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&tb->rhi,1,MPI_DOUBLE,0,world);
  }
  MPI_Bcast(&tb->fpflag,1,MPI_INT,0,world);
  if (tb->fpflag) {
    MPI_Bcast(&tb->fplo,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&tb->fphi,1,MPI_DOUBLE,0,world);
  }
}

/* ----------------------------------------------------------------------
   build spline representation of e,f over entire range of read-in table
   this function sets these values in Table: e2file,f2file
------------------------------------------------------------------------- */

void PairMultiLucyRX::spline_table(Table *tb)
{
  memory->create(tb->e2file,tb->ninput,"pair:e2file");
  memory->create(tb->f2file,tb->ninput,"pair:f2file");

  double ep0 = - tb->ffile[0];
  double epn = - tb->ffile[tb->ninput-1];
  spline(tb->rfile,tb->efile,tb->ninput,ep0,epn,tb->e2file);

  if (tb->fpflag == 0) {
    tb->fplo = (tb->ffile[1] - tb->ffile[0]) / (tb->rfile[1] - tb->rfile[0]);
    tb->fphi = (tb->ffile[tb->ninput-1] - tb->ffile[tb->ninput-2]) /
      (tb->rfile[tb->ninput-1] - tb->rfile[tb->ninput-2]);
  }

  double fp0 = tb->fplo;
  double fpn = tb->fphi;
  spline(tb->rfile,tb->ffile,tb->ninput,fp0,fpn,tb->f2file);
}

/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: N value R/RSQ lo hi FP fplo fphi
   N is required, other params are optional
------------------------------------------------------------------------- */

void PairMultiLucyRX::param_extract(Table *tb, char *line)
{
  char *r_token;
  tb->ninput = 0;
  tb->rflag = NONE;
  tb->fpflag = 0;
  r_token = line;

  char *word = utils::strtok_r(r_token," \t\n\r\f",&r_token);
  while (word) {
    if (strcmp(word,"N") == 0) {
      word = utils::strtok_r(NULL," \t\n\r\f",&r_token);
      tb->ninput = atoi(word);
    } else if (strcmp(word,"R") == 0 || strcmp(word,"RSQ") == 0) {
      if (strcmp(word,"R") == 0) tb->rflag = RLINEAR;
      else if (strcmp(word,"RSQ") == 0) tb->rflag = RSQ;
      word = utils::strtok_r(NULL," \t\n\r\f",&r_token);
      tb->rlo = atof(word);
      word = utils::strtok_r(NULL," \t\n\r\f",&r_token);
      tb->rhi = atof(word);
    } else if (strcmp(word,"FP") == 0) {
      tb->fpflag = 1;
      word = utils::strtok_r(NULL," \t\n\r\f",&r_token);
      tb->fplo = atof(word);
      word = utils::strtok_r(NULL," \t\n\r\f",&r_token);
      tb->fphi = atof(word);
    } else {
      printf("WORD: %s\n",word);
      error->one(FLERR,"Invalid keyword in pair table parameters");
    }
    word = utils::strtok_r(NULL," \t\n\r\f",&r_token);
  }

  if (tb->ninput == 0) error->one(FLERR,"Pair table parameters did not set N");
}

/* ----------------------------------------------------------------------
   compute r,e,f vectors from splined values
------------------------------------------------------------------------- */

void PairMultiLucyRX::compute_table(Table *tb)
{
  int tlm1 = tablength-1;

  // inner = inner table bound
  // cut = outer table bound
  // delta = table spacing in rsq for N-1 bins

  double inner;
  if (tb->rflag) inner = tb->rlo;
  else inner = tb->rfile[0];
  tb->innersq = inner*inner;
  tb->delta = (tb->rhi*tb->rhi - tb->innersq) / tlm1;
  tb->invdelta = 1.0/tb->delta;

  // direct lookup tables
  // N-1 evenly spaced bins in rsq from inner to cut
  // e,f = value at midpt of bin
  // e,f are N-1 in length since store 1 value at bin midpt
  // f is converted to f/r when stored in f[i]
  // e,f are never a match to read-in values, always computed via spline interp

  if (tabstyle == LOOKUP) {
    memory->create(tb->e,tlm1,"pair:e");
    memory->create(tb->f,tlm1,"pair:f");

    double r,rsq;
    for (int i = 0; i < tlm1; i++) {
      rsq = tb->innersq + (i+0.5)*tb->delta;
      r = sqrt(rsq);
      tb->e[i] = splint(tb->rfile,tb->efile,tb->e2file,tb->ninput,r);
      tb->f[i] = splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,r);
    }
  }

  // linear tables
  // N-1 evenly spaced bins in rsq from inner to cut
  // rsq,e,f = value at lower edge of bin
  // de,df values = delta from lower edge to upper edge of bin
  // rsq,e,f are N in length so de,df arrays can compute difference
  // f is converted to f/r when stored in f[i]
  // e,f can match read-in values, else compute via spline interp

  if (tabstyle == LINEAR) {
    memory->create(tb->rsq,tablength,"pair:rsq");
    memory->create(tb->e,tablength,"pair:e");
    memory->create(tb->f,tablength,"pair:f");
    memory->create(tb->de,tlm1,"pair:de");
    memory->create(tb->df,tlm1,"pair:df");

    double r,rsq;
    for (int i = 0; i < tablength; i++) {
      rsq = tb->innersq + i*tb->delta;
      r = sqrt(rsq);
      tb->rsq[i] = rsq;
      if (tb->match) {
        tb->e[i] = tb->efile[i];
        tb->f[i] = tb->ffile[i];
      } else {
        tb->e[i] = splint(tb->rfile,tb->efile,tb->e2file,tb->ninput,r);
        tb->f[i] = splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,r);
      }
    }

    for (int i = 0; i < tlm1; i++) {
      tb->de[i] = tb->e[i+1] - tb->e[i];
      tb->df[i] = tb->f[i+1] - tb->f[i];
    }
  }
}

/* ----------------------------------------------------------------------
   set all ptrs in a table to NULL, so can be freed safely
------------------------------------------------------------------------- */

void PairMultiLucyRX::null_table(Table *tb)
{
  tb->rfile = tb->efile = tb->ffile = NULL;
  tb->e2file = tb->f2file = NULL;
  tb->rsq = tb->drsq = tb->e = tb->de = NULL;
  tb->f = tb->df = tb->e2 = tb->f2 = NULL;
}

/* ----------------------------------------------------------------------
   free all arrays in a table
------------------------------------------------------------------------- */

void PairMultiLucyRX::free_table(Table *tb)
{
  memory->destroy(tb->rfile);
  memory->destroy(tb->efile);
  memory->destroy(tb->ffile);
  memory->destroy(tb->e2file);
  memory->destroy(tb->f2file);

  memory->destroy(tb->rsq);
  memory->destroy(tb->drsq);
  memory->destroy(tb->e);
  memory->destroy(tb->de);
  memory->destroy(tb->f);
  memory->destroy(tb->df);
  memory->destroy(tb->e2);
  memory->destroy(tb->f2);
}

/* ----------------------------------------------------------------------
   spline and splint routines modified from Numerical Recipes
------------------------------------------------------------------------- */

void PairMultiLucyRX::spline(double *x, double *y, int n,
                       double yp1, double ypn, double *y2)
{
  int i,k;
  double p,qn,sig,un;
  double *u = new double[n];

  if (yp1 > 0.99e30) y2[0] = u[0] = 0.0;
  else {
    y2[0] = -0.5;
    u[0] = (3.0/(x[1]-x[0])) * ((y[1]-y[0]) / (x[1]-x[0]) - yp1);
  }
  for (i = 1; i < n-1; i++) {
    sig = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
    p = sig*y2[i-1] + 2.0;
    y2[i] = (sig-1.0) / p;
    u[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1]);
    u[i] = (6.0*u[i] / (x[i+1]-x[i-1]) - sig*u[i-1]) / p;
  }
  if (ypn > 0.99e30) qn = un = 0.0;
  else {
    qn = 0.5;
    un = (3.0/(x[n-1]-x[n-2])) * (ypn - (y[n-1]-y[n-2]) / (x[n-1]-x[n-2]));
  }
  y2[n-1] = (un-qn*u[n-2]) / (qn*y2[n-2] + 1.0);
  for (k = n-2; k >= 0; k--) y2[k] = y2[k]*y2[k+1] + u[k];

  delete [] u;
}

/* ---------------------------------------------------------------------- */

double PairMultiLucyRX::splint(double *xa, double *ya, double *y2a, int n, double x)
{
  int klo,khi,k;
  double h,b,a,y;

  klo = 0;
  khi = n-1;
  while (khi-klo > 1) {
    k = (khi+klo) >> 1;
    if (xa[k] > x) khi = k;
    else klo = k;
  }
  h = xa[khi]-xa[klo];
  a = (xa[khi]-x) / h;
  b = (x-xa[klo]) / h;
  y = a*ya[klo] + b*ya[khi] +
    ((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi]) * (h*h)/6.0;
  return y;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMultiLucyRX::write_restart(FILE *fp)
{
  write_restart_settings(fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMultiLucyRX::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMultiLucyRX::write_restart_settings(FILE *fp)
{
  fwrite(&tabstyle,sizeof(int),1,fp);
  fwrite(&tablength,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMultiLucyRX::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&tabstyle,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&tablength,sizeof(int),1,fp,NULL,error);
  }
  MPI_Bcast(&tabstyle,1,MPI_INT,0,world);
  MPI_Bcast(&tablength,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

void PairMultiLucyRX::computeLocalDensity()
{
  double **x = atom->x;
  const int *type = atom->type;
  const int nlocal = atom->nlocal;

  const int inum = list->inum;
  const int *ilist = list->ilist;
  const int *numneigh = list->numneigh;
        int **firstneigh = list->firstneigh;

  const double pi = MathConst::MY_PI;

  const bool newton_pair = force->newton_pair;
  const bool one_type = (atom->ntypes == 1);

  // Special cut-off values for when there's only one type.
  const double cutsq_type11 = cutsq[1][1];
  const double rcut_type11 = sqrt(cutsq_type11);
  const double factor_type11 = 84.0/(5.0*pi*rcut_type11*rcut_type11*rcut_type11);

  double *rho = atom->rho;

  // zero out density
  if (newton_pair) {
    const int m = nlocal + atom->nghost;
    for (int i = 0; i < m; i++) rho[i] = 0.0;
  }
  else
    for (int i = 0; i < nlocal; i++) rho[i] = 0.0;

// rho = density at each atom
// loop over neighbors of my atoms
  for (int ii = 0; ii < inum; ii++){
    const int i = ilist[ii];

    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];

    double rho_i = rho[i];

    const int itype = type[i];
    const int *jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++){
      const int j = (jlist[jj] & NEIGHMASK);
      const int jtype = type[j];

      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx*delx + dely*dely + delz*delz;

      if (one_type) {
        if (rsq < cutsq_type11) {
          const double rcut = rcut_type11;
          const double r_over_rcut = sqrt(rsq) / rcut;
          const double tmpFactor = 1.0 - r_over_rcut;
          const double tmpFactor4 = tmpFactor*tmpFactor*tmpFactor*tmpFactor;
          const double factor = factor_type11*(1.0 + 1.5*r_over_rcut)*tmpFactor4;
          rho_i += factor;
          if (newton_pair || j < nlocal)
            rho[j] += factor;
        }
      } else if (rsq < cutsq[itype][jtype]) {
        const double rcut = sqrt(cutsq[itype][jtype]);
        const double tmpFactor = 1.0-sqrt(rsq)/rcut;
        const double tmpFactor4 = tmpFactor*tmpFactor*tmpFactor*tmpFactor;
        const double factor = (84.0/(5.0*pi*rcut*rcut*rcut))*(1.0+3.0*sqrt(rsq)/(2.0*rcut))*tmpFactor4;
        rho_i += factor;
        if (newton_pair || j < nlocal)
          rho[j] += factor;
      }
    }

    rho[i] = rho_i;
  }
  if (newton_pair) comm->reverse_comm_pair(this);

  comm->forward_comm_pair(this);

}

/* ---------------------------------------------------------------------- */

void PairMultiLucyRX::getMixingWeights(int id, double &mixWtSite1old, double &mixWtSite2old, double &mixWtSite1, double &mixWtSite2)
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
  for (int ispecies = 0; ispecies < nspecies; ispecies++){
    nTotal += atom->dvector[ispecies][id];
    nTotalOld += atom->dvector[ispecies+nspecies][id];
  }

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
      fractionOFAold += atom->dvector[ispecies+nspecies][id] / nTotalOld;
      fractionOFA += atom->dvector[ispecies][id] / nTotal;
    }
    if (isOneFluid(isite1)){
      nMoleculesOld1 = 1.0-(nTotalOld-nMoleculesOFAold);
      nMolecules1 = 1.0-(nTotal-nMoleculesOFA);
      fractionOld1 = fractionOFAold;
      fraction1 = fractionOFA;
    }
    if (isOneFluid(isite2)){
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

/* ---------------------------------------------------------------------- */

int PairMultiLucyRX::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;
  double *rho = atom->rho;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = rho[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairMultiLucyRX::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;
  double *rho = atom->rho;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) rho[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

int PairMultiLucyRX::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  double *rho = atom->rho;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = rho[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void PairMultiLucyRX::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  double *rho = atom->rho;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    rho[j] += buf[m++];
  }
}
