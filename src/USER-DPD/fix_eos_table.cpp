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
   Contributing author: James Larentzos (U.S. Army Research Laboratory)
------------------------------------------------------------------------- */

#include "fix_eos_table.h"
#include <mpi.h>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "utils.h"

#define MAXLINE 1024

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEOStable::FixEOStable(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), ntables(0), tables(NULL)
{
  if (narg != 7) error->all(FLERR,"Illegal fix eos/table command");
  nevery = 1;

  if (strcmp(arg[3],"linear") == 0) tabstyle = LINEAR;
  else error->all(FLERR,"Unknown table style in fix eos/table");

  tablength = force->inumeric(FLERR,arg[5]);
  if (tablength < 2) error->all(FLERR,"Illegal number of eos/table entries");

  ntables = 0;
  tables = NULL;
  int me;
  MPI_Comm_rank(world,&me);
  tables = (Table *)
    memory->srealloc(tables,(ntables+2)*sizeof(Table),"eos:tables");
  Table *tb = &tables[ntables];
  Table *tb2 = &tables[ntables+1];
  null_table(tb);
  null_table(tb2);
  if (me == 0) read_table(tb,tb2,arg[4],arg[6]);
  bcast_table(tb);
  bcast_table(tb2);

  // error check on table parameters

  if (tb->ninput <= 1) error->one(FLERR,"Invalid eos/table length");

  tb->lo = tb->rfile[0];
  tb->hi = tb->rfile[tb->ninput-1];
  if (tb->lo >= tb->hi) error->all(FLERR,"eos/table values are not increasing");

  if (tb2->ninput <= 1) error->one(FLERR,"Invalid eos/table length");

  tb2->lo = tb2->rfile[0];
  tb2->hi = tb2->rfile[tb2->ninput-1];
  if (tb2->lo >= tb2->hi) error->all(FLERR,"eos/table values are not increasing");

  // spline read-in and compute r,e,f vectors within table

  spline_table(tb);
  compute_table(tb);
  spline_table(tb2);
  compute_table(tb2);
  ntables++;

  if (atom->dpd_flag != 1)
    error->all(FLERR,"FixEOStable requires atom_style with internal temperature and energies (e.g. dpd)");
}

/* ---------------------------------------------------------------------- */

FixEOStable::~FixEOStable()
{
  for (int m = 0; m < 2*ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);
}

/* ---------------------------------------------------------------------- */

int FixEOStable::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEOStable::init()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *dpdTheta = atom->dpdTheta;
  double tmp;

  if(this->restart_reset){
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        temperature_lookup(uCond[i]+uMech[i],dpdTheta[i]);
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if(dpdTheta[i] <= 0.0)
          error->one(FLERR,"Internal temperature <= zero");
        energy_lookup(dpdTheta[i],tmp);
        uCond[i] = 0.0;
        uMech[i] = tmp;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixEOStable::post_integrate()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *dpdTheta = atom->dpdTheta;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit){
      temperature_lookup(uCond[i]+uMech[i],dpdTheta[i]);
      if(dpdTheta[i] <= 0.0)
        error->one(FLERR,"Internal temperature <= zero");
    }
}

/* ---------------------------------------------------------------------- */

void FixEOStable::end_of_step()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *dpdTheta = atom->dpdTheta;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit){
      temperature_lookup(uCond[i]+uMech[i],dpdTheta[i]);
      if(dpdTheta[i] <= 0.0)
        error->one(FLERR,"Internal temperature <= zero");
    }
}

/* ---------------------------------------------------------------------- */

void FixEOStable::null_table(Table *tb)
{
  tb->rfile = tb->efile = NULL;
  tb->e2file = NULL;
  tb->r = tb->e = tb->de = NULL;
  tb->e2 = NULL;
}

/* ---------------------------------------------------------------------- */

void FixEOStable::free_table(Table *tb)
{
  memory->destroy(tb->rfile);
  memory->destroy(tb->efile);
  memory->destroy(tb->e2file);

  memory->destroy(tb->r);
  memory->destroy(tb->e);
  memory->destroy(tb->de);
  memory->destroy(tb->e2);
}

/* ----------------------------------------------------------------------
   read table file, only called by proc 0
------------------------------------------------------------------------- */

void FixEOStable::read_table(Table *tb, Table *tb2, char *file, char *keyword)
{
  char line[MAXLINE];

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
    if (strspn(line," \t\n\r") == strlen(line)) continue;    // blank line
    if (line[0] == '#') continue;                          // comment
    char *word = strtok(line," \t\n\r");
    if (strcmp(word,keyword) == 0) break;           // matching keyword
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);                         // no match, skip section
    param_extract(tb,tb2,line);
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
    for (int i = 0; i < tb->ninput; i++) utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
  }

  // read args on 2nd line of section
  // allocate table arrays for file values

  utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
  param_extract(tb,tb2,line);
  memory->create(tb->rfile,tb->ninput,"eos:rfile");
  memory->create(tb->efile,tb->ninput,"eos:efile");
  memory->create(tb2->rfile,tb2->ninput,"eos:rfile2");
  memory->create(tb2->efile,tb2->ninput,"eos:efile2");

  // read r,e table values from file

  int itmp;
  utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
  for (int i = 0; i < tb->ninput; i++) {
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
    sscanf(line,"%d %lg %lg",&itmp,&tb->rfile[i],&tb->efile[i]);
    sscanf(line,"%d %lg %lg",&itmp,&tb2->efile[i],&tb2->rfile[i]);
  }

  fclose(fp);
}

/* ----------------------------------------------------------------------
   build spline representation of e,f over entire range of read-in table
   this function sets these values in e2file
------------------------------------------------------------------------- */

void FixEOStable::spline_table(Table *tb)
{
  memory->create(tb->e2file,tb->ninput,"eos:e2file");

  double ep0 = 0.0;
  double epn = 0.0;
  spline(tb->rfile,tb->efile,tb->ninput,ep0,epn,tb->e2file);
}

/* ----------------------------------------------------------------------
   compute r,e,f vectors from splined values
------------------------------------------------------------------------- */

void FixEOStable::compute_table(Table *tb)
{
  // delta = table spacing for N-1 bins
  int tlm1 = tablength-1;

  tb->delta = (tb->hi - tb->lo)/ tlm1;
  tb->invdelta = 1.0/tb->delta;
  tb->deltasq6 = tb->delta*tb->delta / 6.0;

  // N-1 evenly spaced bins in r from min to max
  // r,e = value at lower edge of bin
  // de values = delta values of e,f
  // r,e are N in length so de arrays can compute difference

  memory->create(tb->r,tablength,"eos:r");
  memory->create(tb->e,tablength,"eos:e");
  memory->create(tb->de,tlm1,"eos:de");
  memory->create(tb->e2,tablength,"eos:e2");

  double a;
  for (int i = 0; i < tablength; i++) {
    a = tb->lo + i*tb->delta;
    tb->r[i] = a;
    tb->e[i] = splint(tb->rfile,tb->efile,tb->e2file,tb->ninput,a);
  }

  for (int i = 0; i < tlm1; i++) {
    tb->de[i] = tb->e[i+1] - tb->e[i];
  }
}

/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: N value
   N is required, other params are optional
------------------------------------------------------------------------- */

void FixEOStable::param_extract(Table *tb, Table *tb2, char *line)
{
  tb->ninput = 0;
  tb2->ninput = 0;

  char *word = strtok(line," \t\n\r\f");
  while (word) {
    if (strcmp(word,"N") == 0) {
      word = strtok(NULL," \t\n\r\f");
      tb->ninput = atoi(word);
      tb2->ninput = atoi(word);
    } else {
      error->one(FLERR,"Invalid keyword in fix eos/table parameters");
    }
    word = strtok(NULL," \t\n\r\f");
  }

  if (tb->ninput == 0) error->one(FLERR,"fix eos/table parameters did not set N");
  if (tb2->ninput == 0) error->one(FLERR,"fix eos/table parameters did not set N");
}

/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,rfile,efile
------------------------------------------------------------------------- */

void FixEOStable::bcast_table(Table *tb)
{
  MPI_Bcast(&tb->ninput,1,MPI_INT,0,world);

  int me;
  MPI_Comm_rank(world,&me);
  if (me > 0) {
    memory->create(tb->rfile,tb->ninput,"eos:rfile");
    memory->create(tb->efile,tb->ninput,"eos:efile");
  }

  MPI_Bcast(tb->rfile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->efile,tb->ninput,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   spline and splint routines modified from Numerical Recipes
------------------------------------------------------------------------- */

void FixEOStable::spline(double *x, double *y, int n,
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

double FixEOStable::splint(double *xa, double *ya, double *y2a, int n, double x)
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
   calculate internal energy u at temperature t
   insure t is between min/max
------------------------------------------------------------------------- */

void FixEOStable::energy_lookup(double t, double &u)
{
  int itable;
  double fraction;

  Table *tb = &tables[0];
  if(t < tb->lo || t > tb->hi){
    printf("Temperature=%lf TableMin=%lf TableMax=%lf\n",t,tb->lo,tb->hi);
    error->one(FLERR,"Temperature is not within table cutoffs");
  }

  if (tabstyle == LINEAR) {
    itable = static_cast<int> ((t - tb->lo) * tb->invdelta);
    fraction = (t - tb->r[itable]) * tb->invdelta;
    u = tb->e[itable] + fraction*tb->de[itable];
  }
}
/* ----------------------------------------------------------------------
   calculate temperature t at energy u
   insure u is between min/max
------------------------------------------------------------------------- */

void FixEOStable::temperature_lookup(double u, double &t)
{
  int itable;
  double fraction;

  Table *tb = &tables[1];
  if(u < tb->lo || u > tb->hi){
    printf("Energy=%lf TableMin=%lf TableMax=%lf\n",u,tb->lo,tb->hi);
    error->one(FLERR,"Energy is not within table cutoffs");
  }

  if (tabstyle == LINEAR) {
    itable = static_cast<int> ((u - tb->lo) * tb->invdelta);
    fraction = (u - tb->r[itable]) * tb->invdelta;
    t = tb->e[itable] + fraction*tb->de[itable];
  }
}
