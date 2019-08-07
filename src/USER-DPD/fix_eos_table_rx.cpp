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

#include "fix_eos_table_rx.h"
#include <mpi.h>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "atom.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "comm.h"
#include "modify.h"

#define MAXLINE 1024

#ifdef DBL_EPSILON
  #define MY_EPSILON (10.0*DBL_EPSILON)
#else
  #define MY_EPSILON (10.0*2.220446049250313e-16)
#endif

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEOStableRX::FixEOStableRX(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), ntables(0), tables(NULL),
  tables2(NULL), dHf(NULL), eosSpecies(NULL)
{
  if (narg != 8 && narg != 10) error->all(FLERR,"Illegal fix eos/table/rx command");
  nevery = 1;

  rx_flag = false;
  nspecies = 1;
  for (int i = 0; i < modify->nfix; i++)
    if (strncmp(modify->fix[i]->style,"rx",2) == 0){
      rx_flag = true;
      nspecies = atom->nspecies_dpd;
      if(nspecies==0) error->all(FLERR,"There are no rx species specified.");
    }

  if (strcmp(arg[3],"linear") == 0) tabstyle = LINEAR;
  else error->all(FLERR,"Unknown table style in fix eos/table/rx");

  tablength = force->inumeric(FLERR,arg[5]);
  if (tablength < 2) error->all(FLERR,"Illegal number of eos/table/rx entries");

  ntables = 0;
  tables = NULL;
  tables2 = NULL;
  eosSpecies = NULL;

  int me;
  MPI_Comm_rank(world,&me);
  for (int ii=0;ii<nspecies;ii++){
    tables = (Table *)
      memory->srealloc(tables,(ntables+1)*sizeof(Table),"eos:table/rx");
    tables2 = (Table *)
      memory->srealloc(tables2,(ntables+1)*sizeof(Table),"eos:table/rx");

    Table *tb = &tables[ntables];
    Table *tb2 = &tables2[ntables];

    null_table(tb);
    null_table(tb2);

    ntables++;
  }

  ntables = 0;
  Table *tb = &tables[ntables];
  Table *tb2 = &tables2[ntables];

  if (me == 0) read_table(tb,tb2,arg[4],arg[6]);

  for (int ii=0;ii<nspecies;ii++){
    Table *tb = &tables[ntables];
    Table *tb2 = &tables2[ntables];

    bcast_table(tb);
    bcast_table(tb2);

    // error check on table parameters

    if (tb->ninput <= 1) error->one(FLERR,"Invalid eos/table/rx length");

    tb->lo = tb->rfile[0];
    tb->hi = tb->rfile[tb->ninput-1];
    if (tb->lo >= tb->hi) error->all(FLERR,"eos/table/rx values are not increasing");

    if (tb2->ninput <= 1) error->one(FLERR,"Invalid eos/table/rx length");

    tb2->lo = tb2->rfile[0];
    tb2->hi = tb2->rfile[tb2->ninput-1];
    if (tb2->lo >= tb2->hi) error->all(FLERR,"eos/table/rx values are not increasing");

    // spline read-in and compute r,e,f vectors within table

    spline_table(tb);
    compute_table(tb);
    spline_table(tb2);
    compute_table(tb2);
    ntables++;
  }

  // Read the Formation Enthalpies and Correction Coefficients
  dHf = new double[nspecies];
  energyCorr = new double[nspecies];
  tempCorrCoeff = new double[nspecies];
  moleculeCorrCoeff= new double[nspecies];
  for (int ii=0; ii<nspecies; ii++){
    dHf[ii] = 0.0;
    energyCorr[ii] = 0.0;
    tempCorrCoeff[ii] = 0.0;
    moleculeCorrCoeff[ii] = 0.0;
  }

  if(rx_flag) read_file(arg[7]);
  else dHf[0] = atof(arg[7]);

  if(narg==10){
    energyCorr[0] = atof(arg[8]);
    tempCorrCoeff[0] = atof(arg[9]);
  }

  comm_forward = 3;
  comm_reverse = 2;

  if (atom->dpd_flag != 1)
    error->all(FLERR,"FixEOStableRX requires atom_style with internal temperature and energies (e.g. dpd)");
}

/* ---------------------------------------------------------------------- */

FixEOStableRX::~FixEOStableRX()
{
  if (copymode) return;

  for (int m = 0; m < ntables; m++) {
    free_table(&tables[m]);
    free_table(&tables2[m]);
  }
  memory->sfree(tables);
  memory->sfree(tables2);

  delete [] dHf;
  delete [] eosSpecies;
  delete [] energyCorr;
  delete [] tempCorrCoeff;
  delete [] moleculeCorrCoeff;
}

/* ---------------------------------------------------------------------- */

int FixEOStableRX::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEOStableRX::setup(int /*vflag*/)
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *uChem = atom->uChem;
  double *dpdTheta = atom->dpdTheta;
  double duChem;
  double *uCG   = atom->uCG;
  double *uCGnew = atom->uCGnew;

  if(!this->restart_reset){
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit){
        duChem = uCG[i] - uCGnew[i];
        uChem[i] += duChem;
        uCG[i] = 0.0;
        uCGnew[i] = 0.0;
      }
  }

  // Communicate the updated momenta and velocities to all nodes
  comm->forward_comm_fix(this);

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      temperature_lookup(i,uCond[i]+uMech[i]+uChem[i],dpdTheta[i]);
}

/* ---------------------------------------------------------------------- */

void FixEOStableRX::init()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *uChem = atom->uChem;
  double *dpdTheta = atom->dpdTheta;
  double tmp;

  if(this->restart_reset){
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        temperature_lookup(i,uCond[i]+uMech[i]+uChem[i],dpdTheta[i]);
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if(dpdTheta[i] <= 0.0)
          error->one(FLERR,"Internal temperature <= zero");
        energy_lookup(i,dpdTheta[i],tmp);
        uCond[i] = 0.0;
        uMech[i] = tmp;
        uChem[i] = 0.0;
      }
  }
}


/* ---------------------------------------------------------------------- */

void FixEOStableRX::post_integrate()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *uChem = atom->uChem;
  double *dpdTheta = atom->dpdTheta;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit){
      temperature_lookup(i,uCond[i]+uMech[i]+uChem[i],dpdTheta[i]);
      if(dpdTheta[i] <= 0.0)
        error->one(FLERR,"Internal temperature <= zero");
    }
}

/* ---------------------------------------------------------------------- */

void FixEOStableRX::end_of_step()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *uChem = atom->uChem;
  double *dpdTheta = atom->dpdTheta;
  double duChem;
  double *uCG   = atom->uCG;
  double *uCGnew = atom->uCGnew;

  // Communicate the ghost uCGnew
  comm->reverse_comm_fix(this);

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit){
      duChem = uCG[i] - uCGnew[i];
      uChem[i] += duChem;
      uCG[i] = 0.0;
      uCGnew[i] = 0.0;
    }

  // Communicate the updated momenta and velocities to all nodes
  comm->forward_comm_fix(this);

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit){
      temperature_lookup(i,uCond[i]+uMech[i]+uChem[i],dpdTheta[i]);
      if(dpdTheta[i] <= 0.0)
        error->one(FLERR,"Internal temperature <= zero");
    }
}

/* ---------------------------------------------------------------------- */

void FixEOStableRX::read_file(char *file)
{
  int min_params_per_line = 2;
  int max_params_per_line = 5;
  char **words = new char*[max_params_per_line+1];

  // open file on proc 0

  FILE *fp;
  fp = NULL;
  if (comm->me == 0) {
    fp = fopen(file,"r");
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open eos table/rx potential file %s",file);
      error->one(FLERR,str);
    }
  }

  // one set of params can span multiple lines
  int n,nwords,ispecies;
  char line[MAXLINE],*ptr;
  int eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < min_params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != min_params_per_line && nwords != max_params_per_line)
      error->all(FLERR,"Incorrect format in eos table/rx potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    for (ispecies = 0; ispecies < nspecies; ispecies++)
      if (strcmp(words[0],&atom->dname[ispecies][0]) == 0) break;

    if (ispecies < nspecies){
      dHf[ispecies] = atof(words[1]);
      if(nwords > min_params_per_line+1){
        energyCorr[ispecies] = atof(words[2]);
        tempCorrCoeff[ispecies] = atof(words[3]);
        moleculeCorrCoeff[ispecies] = atof(words[4]);
      }
    }
  }

  delete [] words;
}

/* ---------------------------------------------------------------------- */

void FixEOStableRX::null_table(Table *tb)
{
  tb->rfile = tb->efile = NULL;
  tb->e2file = NULL;
  tb->r = tb->e = tb->de = NULL;
  tb->e2 = NULL;
}

/* ---------------------------------------------------------------------- */

void FixEOStableRX::free_table(Table *tb)
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

void FixEOStableRX::read_table(Table *tb, Table *tb2, char *file, char *keyword)
{
  char line[MAXLINE];

  // open file

  FILE *fp = fopen(file,"r");
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
    fgets(line,MAXLINE,fp);                         // no match, skip section
    param_extract(tb,line);
    fgets(line,MAXLINE,fp);
    for (int i = 0; i < tb->ninput; i++) fgets(line,MAXLINE,fp);
  }

  // read args on 2nd line of section
  // allocate table arrays for file values

  fgets(line,MAXLINE,fp);
  param_extract(tb,line);
  tb2->ninput = tb->ninput;
  memory->create(tb->rfile,tb->ninput,"eos:rfile");
  memory->create(tb->efile,tb->ninput,"eos:efile");
  memory->create(tb2->rfile,tb2->ninput,"eos:rfile");
  memory->create(tb2->efile,tb2->ninput,"eos:efile");

  for (int ispecies=1;ispecies<nspecies;ispecies++){
    Table *tbl = &tables[ispecies];
    Table *tbl2 = &tables2[ispecies];
    tbl->ninput = tb->ninput;
    tbl2->ninput = tb2->ninput;

    memory->create(tbl->rfile,tbl->ninput,"eos:rfile");
    memory->create(tbl->efile,tbl->ninput,"eos:efile");
    memory->create(tbl2->rfile,tbl2->ninput,"eos:rfile");
    memory->create(tbl2->efile,tbl2->ninput,"eos:efile");
  }

  // read r,e table values from file

  double rtmp, tmpE;
  int nwords;
  char * word;
  int ispecies;
  int ninputs = tb->ninput;

  fgets(line,MAXLINE,fp);
  for (int i = 0; i < ninputs; i++) {
    fgets(line,MAXLINE,fp);

    nwords = atom->count_words(line);
    if(nwords != nspecies+2){
      printf("nwords=%d  nspecies=%d\n",nwords,nspecies);
      error->all(FLERR,"Illegal fix eos/table/rx command");
    }
    nwords = 0;
    word = strtok(line," \t\n\r\f");
    word = strtok(NULL," \t\n\r\f");
    rtmp = atof(word);

    for (int icolumn=0;icolumn<ncolumn;icolumn++){
      ispecies = eosSpecies[icolumn];

      Table *tbl = &tables[ispecies];
      Table *tbl2 = &tables2[ispecies];

      word = strtok(NULL," \t\n\r\f");
      tmpE = atof(word);

      tbl->rfile[i] = rtmp;
      tbl->efile[i] = tmpE;

      tbl2->rfile[i] = tmpE;
      tbl2->efile[i] = rtmp;
    }
  }
  fclose(fp);
}

/* ----------------------------------------------------------------------
   build spline representation of e,f over entire range of read-in table
   this function sets these values in e2file
------------------------------------------------------------------------- */

void FixEOStableRX::spline_table(Table *tb)
{
  memory->create(tb->e2file,tb->ninput,"eos:e2file");

  double ep0 = 0.0;
  double epn = 0.0;
  spline(tb->rfile,tb->efile,tb->ninput,ep0,epn,tb->e2file);

}

/* ----------------------------------------------------------------------
   compute r,e,f vectors from splined values
------------------------------------------------------------------------- */

void FixEOStableRX::compute_table(Table *tb)
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

void FixEOStableRX::param_extract(Table *tb, char *line)
{
  int ispecies;
  ncolumn = 0;

  if (!eosSpecies)
    eosSpecies = new int[nspecies];
  for (ispecies = 0; ispecies < nspecies; ispecies++)
    eosSpecies[ispecies] = -1;

  tb->ninput = 0;

  char *word = strtok(line," \t\n\r\f");
  if (strcmp(word,"N") == 0) {
    word = strtok(NULL," \t\n\r\f");
    tb->ninput = atoi(word);
  } else
    error->one(FLERR,"Invalid keyword in fix eos/table/rx parameters");
  word = strtok(NULL," \t\n\r\f");

  if(rx_flag){
    while (word) {
      for (ispecies = 0; ispecies < nspecies; ispecies++)
        if (strcmp(word,&atom->dname[ispecies][0]) == 0){
          eosSpecies[ncolumn] =  ispecies;
          ncolumn++;
          break;
        }
      if (ispecies == nspecies){
        printf("name=%s not found in species list\n",word);
        error->one(FLERR,"Invalid keyword in fix eos/table/rx parameters");
      }
      word = strtok(NULL," \t\n\r\f");
    }

    for (int icolumn = 0; icolumn < ncolumn; icolumn++)
      if(eosSpecies[icolumn]==-1)
        error->one(FLERR,"EOS data is missing from fix eos/table/rx tabe");
    if(ncolumn != nspecies){
      printf("ncolumns=%d nspecies=%d\n",ncolumn,nspecies);
      error->one(FLERR,"The number of columns in fix eos/table/rx does not match the number of species");
    }
  } else {
    eosSpecies[0] = 0;
    ncolumn++;
  }

  if (tb->ninput == 0) error->one(FLERR,"fix eos/table/rx parameters did not set N");

}

/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,rfile,efile
------------------------------------------------------------------------- */

void FixEOStableRX::bcast_table(Table *tb)
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

void FixEOStableRX::spline(double *x, double *y, int n,
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

double FixEOStableRX::splint(double *xa, double *ya, double *y2a, int n, double x)
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
   calculate potential ui at temperature thetai
------------------------------------------------------------------------- */

void FixEOStableRX::energy_lookup(int id, double thetai, double &ui)
{
  int itable, nPG;
  double fraction, uTmp, nMolecules, nTotal, nTotalPG;
  double tolerance = 1.0e-10;

  ui = 0.0;
  nTotal = 0.0;
  nTotalPG = 0.0;
  nPG = 0;

  if(rx_flag){
    for(int ispecies=0;ispecies<nspecies;ispecies++){
      nTotal += atom->dvector[ispecies][id];
      if(fabs(moleculeCorrCoeff[ispecies]) > tolerance){
        nPG++;
        nTotalPG += atom->dvector[ispecies][id];
      }
    }
  } else {
    nTotal = 1.0;
  }

  for(int ispecies=0;ispecies<nspecies;ispecies++){
    Table *tb = &tables[ispecies];
    thetai = MAX(thetai,tb->lo);
    thetai = MIN(thetai,tb->hi);

    if (tabstyle == LINEAR) {
      itable = static_cast<int> ((thetai - tb->lo) * tb->invdelta);
      fraction = (thetai - tb->r[itable]) * tb->invdelta;
      uTmp = tb->e[itable] + fraction*tb->de[itable];

      uTmp += dHf[ispecies];
      uTmp += tempCorrCoeff[ispecies]*thetai; // temperature correction
      uTmp += energyCorr[ispecies]; // energy correction
      if(nPG > 0) ui += moleculeCorrCoeff[ispecies]*nTotalPG/double(nPG); // molecule correction

      if(rx_flag) nMolecules = atom->dvector[ispecies][id];
      else nMolecules = 1.0;
      ui += nMolecules*uTmp;
    }
  }
  ui = ui - double(nTotal+1.5)*force->boltz*thetai;
}

/* ----------------------------------------------------------------------
   calculate temperature thetai at energy ui
------------------------------------------------------------------------- */

void FixEOStableRX::temperature_lookup(int id, double ui, double &thetai)
{
  Table *tb = &tables[0];

  int it;
  double t1,t2,u1,u2,f1,f2;
  double maxit = 100;
  double temp;
  double delta = 0.001;
  double tolerance = 1.0e-10;

  // Store the current thetai in t1
  t1 = MAX(thetai,tb->lo);
  t1 = MIN(t1,tb->hi);
  if(t1==tb->hi) delta = -delta;

  // Compute u1 at thetai
  energy_lookup(id,t1,u1);

  // Compute f1
  f1 = u1 - ui;

  // Compute guess of t2
  t2 = (1.0 + delta)*t1;

  // Compute u2 at t2
  energy_lookup(id,t2,u2);

  // Compute f1
  f2 = u2 - ui;

  // Apply the Secant Method
  for(it=0; it<maxit; it++){
    if(fabs(f2-f1) < MY_EPSILON){
      if(std::isnan(f1) || std::isnan(f2)) error->one(FLERR,"NaN detected in secant solver.");
      temp = t1;
      temp = MAX(temp,tb->lo);
      temp = MIN(temp,tb->hi);
      char str[256];
      sprintf(str,"Secant solver did not converge because table bounds were exceeded:  it=%d id=%d ui=%lf thetai=%lf t1=%lf t2=%lf f1=%lf f2=%lf dpdTheta=%lf\n",it,id,ui,thetai,t1,t2,f1,f2,temp);
      error->warning(FLERR,str);
      break;
    }
    temp = t2 - f2*(t2-t1)/(f2-f1);
    if(fabs(temp-t2) < tolerance) break;
    f1 = f2;
    t1 = t2;
    t2 = temp;
    energy_lookup(id,t2,u2);
    f2 = u2 - ui;
  }
  if(it==maxit){
    char str[256];
    sprintf(str,"Maxit exceeded in secant solver:  id=%d ui=%lf thetai=%lf t1=%lf t2=%lf f1=%lf f2=%lf\n",id,ui,thetai,t1,t2,f1,f2);
    if(std::isnan(f1) || std::isnan(f2) || std::isnan(ui) || std::isnan(thetai) || std::isnan(t1) || std::isnan(t2))
      error->one(FLERR,"NaN detected in secant solver.");
    error->one(FLERR,str);
  }
  thetai = temp;
}

/* ---------------------------------------------------------------------- */

int FixEOStableRX::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int ii,jj,m;
  double *uChem = atom->uChem;
  double *uCG = atom->uCG;
  double *uCGnew = atom->uCGnew;

  m = 0;
  for (ii = 0; ii < n; ii++) {
    jj = list[ii];
    buf[m++] = uChem[jj];
    buf[m++] = uCG[jj];
    buf[m++] = uCGnew[jj];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixEOStableRX::unpack_forward_comm(int n, int first, double *buf)
{
  int ii,m,last;
  double *uChem = atom->uChem;
  double *uCG = atom->uCG;
  double *uCGnew = atom->uCGnew;

  m = 0;
  last = first + n ;
  for (ii = first; ii < last; ii++){
    uChem[ii]  = buf[m++];
    uCG[ii]    = buf[m++];
    uCGnew[ii] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int FixEOStableRX::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  double *uCG = atom->uCG;
  double *uCGnew = atom->uCGnew;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = uCG[i];
    buf[m++] = uCGnew[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixEOStableRX::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  double *uCG = atom->uCG;
  double *uCGnew = atom->uCGnew;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];

    uCG[j] += buf[m++];
    uCGnew[j] += buf[m++];
  }
}
