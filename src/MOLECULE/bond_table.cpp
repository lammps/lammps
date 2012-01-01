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
   Contributing author: Chuanfu Luo (luochuanfu@gmail.com)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "bond_table.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{LINEAR,SPLINE};

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

BondTable::BondTable(LAMMPS *lmp) : Bond(lmp) 
{
  ntables = 0;
  tables = NULL;
}

/* ---------------------------------------------------------------------- */

BondTable::~BondTable()
{
  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(r0);
    memory->destroy(tabindex);
  }
}

/* ---------------------------------------------------------------------- */

void BondTable::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r;
  double u,mdu;

  ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    domain->minimum_image(delx,dely,delz);

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);

    // force & energy

    uf_lookup(type,r,u,mdu);
    fbond = mdu/r;
    ebond = u;

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx*fbond;
      f[i1][1] += dely*fbond;
      f[i1][2] += delz*fbond;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx*fbond;
      f[i2][1] -= dely*fbond;
      f[i2][2] -= delz*fbond;
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondTable::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(tabindex,n+1,"bond:tabindex");
  memory->create(r0,n+1,"bond:r0");
  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}


/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void BondTable::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal bond_style command");

  if (strcmp(arg[0],"linear") == 0) tabstyle = LINEAR;
  else if (strcmp(arg[0],"spline") == 0) tabstyle = SPLINE;
  else error->all(FLERR,"Unknown table style in bond style table");

  tablength = force->inumeric(arg[1]);
  if (tablength < 2) error->all(FLERR,"Illegal number of bond table entries");

  // delete old tables, since cannot just change settings

  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
     memory->destroy(setflag);
     memory->destroy(tabindex);
  }
  allocated = 0;

  ntables = 0;
  tables = NULL;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void BondTable::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal bond_coeff command");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);
  
  int me;
  MPI_Comm_rank(world,&me);
  tables = (Table *) 
    memory->srealloc(tables,(ntables+1)*sizeof(Table),"bond:tables");
  Table *tb = &tables[ntables];
  null_table(tb);
  if (me == 0) read_table(tb,arg[1],arg[2]);
  bcast_table(tb);

  // error check on table parameters

  if (tb->ninput <= 1) error->one(FLERR,"Invalid bond table length");

  tb->lo = tb->rfile[0];
  tb->hi = tb->rfile[tb->ninput-1];
  if (tb->lo >= tb->hi) error->all(FLERR,"Bond table values are not increasing");

  // spline read-in and compute r,e,f vectors within table

  spline_table(tb);
  compute_table(tb);

  // store ptr to table in tabindex

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    tabindex[i] = ntables;
    r0[i] = tb->r0;
    setflag[i] = 1;
    count++;
  }
  ntables++;

  if (count == 0) error->all(FLERR,"Illegal bond_coeff command");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
   should not be used, since don't know minimum of tabulated function
------------------------------------------------------------------------- */

double BondTable::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void BondTable::write_restart(FILE *fp)
{
  fwrite(&tabstyle,sizeof(int),1,fp);
  fwrite(&tablength,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void BondTable::read_restart(FILE *fp)
{
  if (comm->me == 0) {
    fread(&tabstyle,sizeof(int),1,fp);
    fread(&tablength,sizeof(int),1,fp);
  }
  MPI_Bcast(&tabstyle,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&tablength,1,MPI_INT,0,world);

  allocate();
}

/* ---------------------------------------------------------------------- */

double BondTable::single(int type, double rsq, int i, int j)
{
  double r = sqrt(rsq);
  double u;
  u_lookup(type,r,u);
  return u;
}

/* ---------------------------------------------------------------------- */

void BondTable::null_table(Table *tb)
{
  tb->rfile = tb->efile = tb->ffile = NULL;
  tb->e2file = tb->f2file = NULL;
  tb->r = tb->e = tb->de = NULL;
  tb->f = tb->df = tb->e2 = tb->f2 = NULL;
}

/* ---------------------------------------------------------------------- */

void BondTable::free_table(Table *tb)
{
  memory->destroy(tb->rfile);
  memory->destroy(tb->efile);
  memory->destroy(tb->ffile);
  memory->destroy(tb->e2file);
  memory->destroy(tb->f2file);
  
  memory->destroy(tb->r);
  memory->destroy(tb->e);
  memory->destroy(tb->de);
  memory->destroy(tb->f);
  memory->destroy(tb->df);
  memory->destroy(tb->e2);
  memory->destroy(tb->f2);
}

/* ----------------------------------------------------------------------
   read table file, only called by proc 0
------------------------------------------------------------------------- */

void BondTable::read_table(Table *tb, char *file, char *keyword)
{
  char line[MAXLINE];

  // open file

  FILE *fp = fopen(file,"r");
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",file);
    error->one(FLERR,str);
  }

  // loop until section found with matching keyword

  while (1) {
    if (fgets(line,MAXLINE,fp) == NULL)
      error->one(FLERR,"Did not find keyword in table file");
    if (strspn(line," \t\n") == strlen(line)) continue;    // blank line
    if (line[0] == '#') continue;                          // comment
    if (strstr(line,keyword) == line) break;               // matching keyword
    fgets(line,MAXLINE,fp);                         // no match, skip section
    param_extract(tb,line);
    fgets(line,MAXLINE,fp);
    for (int i = 0; i < tb->ninput; i++) fgets(line,MAXLINE,fp);
  }

  // read args on 2nd line of section
  // allocate table arrays for file values

  fgets(line,MAXLINE,fp);
  param_extract(tb,line);
  memory->create(tb->rfile,tb->ninput,"bond:rfile");
  memory->create(tb->efile,tb->ninput,"bond:efile");
  memory->create(tb->ffile,tb->ninput,"bond:ffile");

  // read r,e,f table values from file

  int itmp;
  fgets(line,MAXLINE,fp);
  for (int i = 0; i < tb->ninput; i++) {
    fgets(line,MAXLINE,fp);
    sscanf(line,"%d %lg %lg %lg",
      &itmp,&tb->rfile[i],&tb->efile[i],&tb->ffile[i]);
  }

  fclose(fp);
}

/* ----------------------------------------------------------------------
   build spline representation of e,f over entire range of read-in table
   this function sets these values in e2file,f2file
------------------------------------------------------------------------- */

void BondTable::spline_table(Table *tb)
{
  memory->create(tb->e2file,tb->ninput,"bond:e2file");
  memory->create(tb->f2file,tb->ninput,"bond:f2file");

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
   compute r,e,f vectors from splined values
------------------------------------------------------------------------- */

void BondTable::compute_table(Table *tb)
{
  // delta = table spacing for N-1 bins
  int tlm1 = tablength-1;

  tb->delta = (tb->hi - tb->lo)/ tlm1;
  tb->invdelta = 1.0/tb->delta;
  tb->deltasq6 = tb->delta*tb->delta / 6.0;
  
  // N-1 evenly spaced bins in r from min to max
  // r,e,f = value at lower edge of bin
  // de,df values = delta values of e,f
  // r,e,f are N in length so de,df arrays can compute difference

  memory->create(tb->r,tablength,"bond:r");
  memory->create(tb->e,tablength,"bond:e");
  memory->create(tb->de,tlm1,"bond:de");
  memory->create(tb->f,tablength,"bond:f");
  memory->create(tb->df,tlm1,"bond:df");
  memory->create(tb->e2,tablength,"bond:e2");
  memory->create(tb->f2,tablength,"bond:f2");

  double a;
  for (int i = 0; i < tablength; i++) {
    a = tb->lo + i*tb->delta;
    tb->r[i] = a;
    tb->e[i] = splint(tb->rfile,tb->efile,tb->e2file,tb->ninput,a);
    tb->f[i] = splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,a);
  }
        
  for (int i = 0; i < tlm1; i++) {
    tb->de[i] = tb->e[i+1] - tb->e[i];
    tb->df[i] = tb->f[i+1] - tb->f[i];
  }
     
  double ep0 = - tb->f[0];
  double epn = - tb->f[tlm1];
  spline(tb->r,tb->e,tablength,ep0,epn,tb->e2);  
  spline(tb->r,tb->f,tablength,tb->fplo,tb->fphi,tb->f2);
}

/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: N value FP fplo fphi EQ r0
   N is required, other params are optional
------------------------------------------------------------------------- */

void BondTable::param_extract(Table *tb, char *line)
{
  tb->ninput = 0;
  tb->fpflag = 0;
  tb->r0 = 0.0;
  
  char *word = strtok(line," \t\n\r\f");
  while (word) {
    if (strcmp(word,"N") == 0) {
      word = strtok(NULL," \t\n\r\f");
      tb->ninput = atoi(word);
    } else if (strcmp(word,"FP") == 0) {
      tb->fpflag = 1;
      word = strtok(NULL," \t\n\r\f");
      tb->fplo = atof(word);
      word = strtok(NULL," \t\n\r\f");
      tb->fphi = atof(word);
    } else if (strcmp(word,"EQ") == 0) {
      word = strtok(NULL," \t\n\r\f");
      tb->r0 = atof(word);
    } else {
      error->one(FLERR,"Invalid keyword in bond table parameters");
    }
    word = strtok(NULL," \t\n\r\f");
  }

  if (tb->ninput == 0) error->one(FLERR,"Bond table parameters did not set N");
}

/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,rfile,efile,ffile,fpflag,fplo,fphi,r0
------------------------------------------------------------------------- */

void BondTable::bcast_table(Table *tb)
{
  MPI_Bcast(&tb->ninput,1,MPI_INT,0,world);
  MPI_Bcast(&tb->r0,1,MPI_INT,0,world);

  int me;
  MPI_Comm_rank(world,&me);
  if (me > 0) {
    memory->create(tb->rfile,tb->ninput,"angle:rfile");
    memory->create(tb->efile,tb->ninput,"angle:efile");
    memory->create(tb->ffile,tb->ninput,"angle:ffile");
  }

  MPI_Bcast(tb->rfile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->efile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->ffile,tb->ninput,MPI_DOUBLE,0,world);

  MPI_Bcast(&tb->fpflag,1,MPI_INT,0,world);
  if (tb->fpflag) {
    MPI_Bcast(&tb->fplo,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&tb->fphi,1,MPI_DOUBLE,0,world);
  }
  MPI_Bcast(&tb->r0,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   spline and splint routines modified from Numerical Recipes
------------------------------------------------------------------------- */

void BondTable::spline(double *x, double *y, int n,
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

double BondTable::splint(double *xa, double *ya, double *y2a, int n, double x)
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
   calculate potential u and force f at distance x
   insure x is between bond min/max
------------------------------------------------------------------------- */

void BondTable::uf_lookup(int type, double x, double &u, double &f)
{
  int itable;
  double fraction,a,b;

  Table *tb = &tables[tabindex[type]];
  x = MAX(x,tb->lo);
  x = MIN(x,tb->hi);

  if (tabstyle == LINEAR) {
    itable = static_cast<int> ((x - tb->lo) * tb->invdelta);
    fraction = (x - tb->r[itable]) * tb->invdelta;
    u = tb->e[itable] + fraction*tb->de[itable];
    f = tb->f[itable] + fraction*tb->df[itable];
  } else if (tabstyle == SPLINE) {
    itable = static_cast<int> ((x - tb->lo) * tb->invdelta);
    fraction = (x - tb->r[itable]) * tb->invdelta;
    
    b = (x - tb->r[itable]) * tb->invdelta;
    a = 1.0 - b;
    u = a * tb->e[itable] + b * tb->e[itable+1] + 
      ((a*a*a-a)*tb->e2[itable] + (b*b*b-b)*tb->e2[itable+1]) * 
      tb->deltasq6;
    f = a * tb->f[itable] + b * tb->f[itable+1] + 
      ((a*a*a-a)*tb->f2[itable] + (b*b*b-b)*tb->f2[itable+1]) * 
      tb->deltasq6;
  } 
}

/* ----------------------------------------------------------------------
   calculate potential u at distance x
   insure x is between bond min/max
------------------------------------------------------------------------- */

void BondTable::u_lookup(int type, double x, double &u)
{
  int itable;
  double fraction,a,b;

  Table *tb = &tables[tabindex[type]];
  x = MAX(x,tb->lo);
  x = MIN(x,tb->hi);
  
  if (tabstyle == LINEAR) {
    itable = static_cast<int> ((x - tb->lo) * tb->invdelta);
    fraction = (x - tb->r[itable]) * tb->invdelta;
    u = tb->e[itable] + fraction*tb->de[itable];
  } else if (tabstyle == SPLINE) {
    itable = static_cast<int> ((x - tb->lo) * tb->invdelta);
    fraction = (x - tb->r[itable]) * tb->invdelta;
    
    b = (x - tb->r[itable]) * tb->invdelta;
    a = 1.0 - b;
    u = a * tb->e[itable] + b * tb->e[itable+1] + 
      ((a*a*a-a)*tb->e2[itable] + (b*b*b-b)*tb->e2[itable+1]) * 
      tb->deltasq6;
  } 
}
