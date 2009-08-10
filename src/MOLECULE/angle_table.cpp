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
#include "angle_table.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

enum{LINEAR,SPLINE};

#define MAXLINE 1024
#define SMALL 0.001
#define TINY  1.E-10

/* ---------------------------------------------------------------------- */

AngleTable::AngleTable(LAMMPS *lmp) : Angle(lmp) 
{
  ntables = 0;
  tables = NULL;
}

/* ---------------------------------------------------------------------- */

AngleTable::~AngleTable()
{
  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);
  
  if (allocated) {
    memory->sfree(setflag);
    memory->sfree(theta0);
    memory->sfree(tabindex);
  }
}

/* ---------------------------------------------------------------------- */

void AngleTable::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type,factor;
  double eangle,f1[3],f3[3];
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22,vx1,vx2,vy1,vy2,vz1,vz2;
  double theta,u,mdu; //mdu: minus du, -du/dx=f

  eangle = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];

    // 1st bond

    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];
    domain->minimum_image(delx1,dely1,delz1);

    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);

    // 2nd bond

    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];
    domain->minimum_image(delx2,dely2,delz2);

    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);

    // angle (cos and sin)

    c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    c /= r1*r2;
        
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
             
    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;
    s = 1.0/s;

    // tabulated force & energy

    theta = acos(c);
    uf_lookup(type,theta,u,mdu);
    
    if (eflag) eangle = u;

    a = mdu * s;               
    a11 = a*c / rsq1;
    a12 = -a / (r1*r2);
    a22 = a*c / rsq2;

    f1[0] = a11*delx1 + a12*delx2;
    f1[1] = a11*dely1 + a12*dely2;
    f1[2] = a11*delz1 + a12*delz2;
    f3[0] = a22*delx2 + a12*delx1;
    f3[1] = a22*dely2 + a12*dely1;
    f3[2] = a22*delz2 + a12*delz1;

    // apply force to each of 3 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
			 delx1,dely1,delz1,delx2,dely2,delz2);
  }
}

/* ---------------------------------------------------------------------- */

void AngleTable::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  theta0 = (double *) memory->smalloc((n+1)*sizeof(double),"angle:theta0");
  tabindex = (int *) memory->smalloc((n+1)*sizeof(int),"angle:tabindex");

  setflag = (int *) memory->smalloc((n+1)*sizeof(int),"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void AngleTable::settings(int narg, char **arg)
{
  if (narg != 2) error->all("Illegal angle_style command");

  if (strcmp(arg[0],"linear") == 0) tabstyle = LINEAR;
  else if (strcmp(arg[0],"spline") == 0) tabstyle = SPLINE;
  else error->all("Unknown table style in angle style table");

  n = atoi(arg[1]);
  nm1 = n - 1;

  // delete old tables, since cannot just change settings

  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
     memory->sfree(setflag);
     memory->sfree(tabindex);
  }
  allocated = 0;

  ntables = 0;
  tables = NULL;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void AngleTable::coeff(int which, int narg, char **arg)
{
  if (which > 0) return;
  if (narg != 3) error->all("Illegal angle_coeff command");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nangletypes,ilo,ihi);
  
  int me;
  MPI_Comm_rank(world,&me);
  tables = (Table *) 
    memory->srealloc(tables,(ntables+1)*sizeof(Table),"angle:tables");
  Table *tb = &tables[ntables];
  null_table(tb);
  if (me == 0) read_table(tb,arg[1],arg[2]);
  bcast_table(tb);

  // error check on table parameters

  if (tb->ninput <= 1) error->one("Invalid angle table length");
  double alo,ahi;
  alo = tb->afile[0];
  ahi = tb->afile[tb->ninput-1];

  if (fabs(alo-0.0) > TINY || fabs(ahi-180.0) > TINY)
    error->all("Angle table must range from 0 to 180 degrees");
    
  // convert theta from degrees to radians

  for (int i = 0; i < tb->ninput; i++){
    tb->afile[i] *= PI/180.0;
    tb->ffile[i] *= 180.0/PI; 
  }

  // spline read-in and compute a,e,f vectors within table

  spline_table(tb);
  compute_table(tb);

  // store ptr to table in tabindex

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    tabindex[i] = ntables;
    setflag[i] = 1;
    theta0[i] = tb->theta0;
    count++;
  }
  ntables++;

  if (count == 0) error->all("Illegal angle_coeff command");
}

/* ----------------------------------------------------------------------
   return an equilbrium angle length
   should not be used, since don't know minimum of tabulated function
------------------------------------------------------------------------- */

double AngleTable::equilibrium_angle(int i)
{
  return theta0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void AngleTable::write_restart(FILE *fp)
{
  fwrite(&tabstyle,sizeof(int),1,fp);
  fwrite(&n,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void AngleTable::read_restart(FILE *fp)
{
  if (comm->me == 0) {
    fread(&tabstyle,sizeof(int),1,fp);
    fread(&n,sizeof(int),1,fp);
  }
  MPI_Bcast(&tabstyle,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&n,1,MPI_INT,0,world);
  nm1 = n - 1;

  allocate();
}

/* ---------------------------------------------------------------------- */

double AngleTable::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;

  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(delx1,dely1,delz1);
  double r1 = sqrt(delx1*delx1 + dely1*dely1 + delz1*delz1);
  
  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2,dely2,delz2);
  double r2 = sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2);

  double c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  double theta = acos(c);
  double u;
  u_lookup(type,theta,u);
  return u;
}

/* ---------------------------------------------------------------------- */

void AngleTable::null_table(Table *tb)
{
  tb->afile = tb->efile = tb->ffile = NULL;
  tb->e2file = tb->f2file = NULL;
  tb->ang = tb->e = tb->de = NULL;
  tb->f = tb->df = tb->e2 = tb->f2 = NULL;
}

/* ---------------------------------------------------------------------- */

void AngleTable::free_table(Table *tb)
{
  memory->sfree(tb->afile);
  memory->sfree(tb->efile);
  memory->sfree(tb->ffile);
  memory->sfree(tb->e2file);
  memory->sfree(tb->f2file);
  
  memory->sfree(tb->ang);
  memory->sfree(tb->e);
  memory->sfree(tb->de);
  memory->sfree(tb->f);
  memory->sfree(tb->df);
  memory->sfree(tb->e2);
  memory->sfree(tb->f2);
}

/* ----------------------------------------------------------------------
   read table file, only called by proc 0
------------------------------------------------------------------------- */

void AngleTable::read_table(Table *tb, char *file, char *keyword)
{
  char line[MAXLINE];

  // open file

  FILE *fp = fopen(file,"r");
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",file);
    error->one(str);
  }

  // loop until section found with matching keyword

  while (1) {
    if (fgets(line,MAXLINE,fp) == NULL)
      error->one("Did not find keyword in table file");
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
  tb->afile = (double *) 
    memory->smalloc(tb->ninput*sizeof(double),"angle:afile");
  tb->efile = (double *) 
    memory->smalloc(tb->ninput*sizeof(double),"angle:efile");
  tb->ffile = (double *) 
    memory->smalloc(tb->ninput*sizeof(double),"angle:ffile");

  // read a,e,f table values from file

  int itmp;
  fgets(line,MAXLINE,fp);
  for (int i = 0; i < tb->ninput; i++) {
    fgets(line,MAXLINE,fp);
    sscanf(line,"%d %lg %lg %lg",
      &itmp,&tb->afile[i],&tb->efile[i],&tb->ffile[i]);
  }

  fclose(fp);
}

/* ----------------------------------------------------------------------
   build spline representation of e,f over entire range of read-in table
   this function sets these values in e2file,f2file
------------------------------------------------------------------------- */

void AngleTable::spline_table(Table *tb)
{
  tb->e2file = (double *) 
    memory->smalloc(tb->ninput*sizeof(double),"angle:e2file");
  tb->f2file = (double *) 
    memory->smalloc(tb->ninput*sizeof(double),"angle:f2file");

  double ep0 = - tb->ffile[0];
  double epn = - tb->ffile[tb->ninput-1];
  spline(tb->afile,tb->efile,tb->ninput,ep0,epn,tb->e2file);

  if (tb->fpflag == 0) {
    tb->fplo = (tb->ffile[1] - tb->ffile[0]) / (tb->afile[1] - tb->afile[0]);
    tb->fphi = (tb->ffile[tb->ninput-1] - tb->ffile[tb->ninput-2]) / 
      (tb->afile[tb->ninput-1] - tb->afile[tb->ninput-2]);
  }

  double fp0 = tb->fplo;
  double fpn = tb->fphi;
  spline(tb->afile,tb->ffile,tb->ninput,fp0,fpn,tb->f2file);
}

/* ----------------------------------------------------------------------
   compute a,e,f vectors from splined values
------------------------------------------------------------------------- */

void AngleTable::compute_table(Table *tb)
{
  // delta = table spacing in angle for N-1 bins

  tb->delta = PI/ nm1;
  tb->invdelta = 1.0/tb->delta;
  tb->deltasq6 = tb->delta*tb->delta / 6.0;
  
  // N-1 evenly spaced bins in angle from 0 to PI
  // ang,e,f = value at lower edge of bin
  // de,df values = delta values of e,f
  // ang,e,f are N in length so de,df arrays can compute difference

  tb->ang = (double *) memory->smalloc(n*sizeof(double),"angle:ang");
  tb->e = (double *) memory->smalloc(n*sizeof(double),"angle:e");
  tb->de = (double *) memory->smalloc(nm1*sizeof(double),"angle:de");
  tb->f = (double *) memory->smalloc(n*sizeof(double),"angle:f");
  tb->df = (double *) memory->smalloc(nm1*sizeof(double),"angle:df");
  tb->e2 = (double *) memory->smalloc(n*sizeof(double),"angle:e2");
  tb->f2 = (double *) memory->smalloc(n*sizeof(double),"angle:f2");

  double a;
  for (int i = 0; i < n; i++) {
    a = i*tb->delta;
    tb->ang[i] = a;
	  tb->e[i] = splint(tb->afile,tb->efile,tb->e2file,tb->ninput,a);
	  tb->f[i] = splint(tb->afile,tb->ffile,tb->f2file,tb->ninput,a);
  }
        
  for (int i = 0; i < nm1; i++) {
    tb->de[i] = tb->e[i+1] - tb->e[i];
    tb->df[i] = tb->f[i+1] - tb->f[i];
  }
     
  double ep0 = - tb->f[0];
  double epn = - tb->f[nm1];
  spline(tb->ang,tb->e,n,ep0,epn,tb->e2);  
  spline(tb->ang,tb->f,n,tb->fplo,tb->fphi,tb->f2);
}

/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: N value FP fplo fphi EQ theta0
   N is required, other params are optional
------------------------------------------------------------------------- */

void AngleTable::param_extract(Table *tb, char *line)
{
  tb->ninput = 0;
  tb->fpflag = 0;
  tb->theta0 = 180.0; 
  
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
      tb->fplo *= (180.0/PI)*(180.0/PI);
      tb->fphi *= (180.0/PI)*(180.0/PI);
    } else if (strcmp(word,"EQ") == 0) {
      word = strtok(NULL," \t\n\r\f");
      tb->theta0 = atof(word);
    } else {
      error->one("Invalid keyword in angle table parameters");
    }
    word = strtok(NULL," \t\n\r\f");
  }

  if (tb->ninput == 0) error->one("Angle table parameters did not set N");
}

/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,afile,efile,ffile,fpflag,fplo,fphi,theta0
------------------------------------------------------------------------- */

void AngleTable::bcast_table(Table *tb)
{
  MPI_Bcast(&tb->ninput,1,MPI_INT,0,world);

  int me;
  MPI_Comm_rank(world,&me);
  if (me > 0) {
    tb->afile = (double *) 
      memory->smalloc(tb->ninput*sizeof(double),"angle:afile");
    tb->efile = (double *) 
      memory->smalloc(tb->ninput*sizeof(double),"angle:efile");
    tb->ffile = (double *) 
      memory->smalloc(tb->ninput*sizeof(double),"angle:ffile");
  }

  MPI_Bcast(tb->afile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->efile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->ffile,tb->ninput,MPI_DOUBLE,0,world);

  MPI_Bcast(&tb->fpflag,1,MPI_INT,0,world);
  if (tb->fpflag) {
    MPI_Bcast(&tb->fplo,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&tb->fphi,1,MPI_DOUBLE,0,world);
  }
  MPI_Bcast(&tb->theta0,1,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   spline and splint routines modified from Numerical Recipes
------------------------------------------------------------------------- */

void AngleTable::spline(double *x, double *y, int n,
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

double AngleTable::splint(double *xa, double *ya, double *y2a, int n, double x)
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
   calculate potential u and force f at angle x
------------------------------------------------------------------------- */

void AngleTable::uf_lookup(int type, double x, double &u, double &f)
{
  int itable;
  double fraction,value,a,b;

  Table *tb = &tables[tabindex[type]];
  
  if (tabstyle == LINEAR) {
    itable = static_cast<int> ( x * tb->invdelta);
    fraction = (x - tb->ang[itable]) * tb->invdelta;
    u = tb->e[itable] + fraction*tb->de[itable];
    f = tb->f[itable] + fraction*tb->df[itable];
  } else if (tabstyle == SPLINE) {
    itable = static_cast<int> ( x * tb->invdelta);
    fraction = (x - tb->ang[itable]) * tb->invdelta;
    
    b = (x - tb->ang[itable]) * tb->invdelta;
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
   calculate potential u at angle x
------------------------------------------------------------------------- */

void AngleTable::u_lookup(int type, double x, double &u)
{
  int itable;
  double fraction,value,a,b;

  Table *tb = &tables[tabindex[type]];
  
  if (tabstyle == LINEAR) {
    itable = static_cast<int> ( x * tb->invdelta);
    fraction = (x - tb->ang[itable]) * tb->invdelta;
    u = tb->e[itable] + fraction*tb->de[itable];
  } else if (tabstyle == SPLINE) {
    itable = static_cast<int> ( x * tb->invdelta);
    fraction = (x - tb->ang[itable]) * tb->invdelta;
    
    b = (x - tb->ang[itable]) * tb->invdelta;
    a = 1.0 - b;
    u = a * tb->e[itable] + b * tb->e[itable+1] + 
      ((a*a*a-a)*tb->e2[itable] + (b*b*b-b)*tb->e2[itable+1]) * 
      tb->deltasq6;
  } 
}
