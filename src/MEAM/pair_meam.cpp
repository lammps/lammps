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
   Contributing author: Greg Wagner (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_meam.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define MAXLINE 1024

enum{FCC,BCC,HCP,DIM,DIAMOND,B1,C11,L12};
int nkeywords = 16;
char *keywords[] = {"Ec","alpha","rho0","delta","lattce",
		    "attrac","repuls","nn2","Cmin","Cmax","rc","delr",
		    "augt1","gsmooth_factor","re","ialloy"};

/* ---------------------------------------------------------------------- */

PairMEAM::PairMEAM(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  one_coeff = 1;

  nmax = 0;
  rho = rho0 = rho1 = rho2 = rho3 = frhop = NULL;
  gamma = dgamma1 = dgamma2 = dgamma3 = arho2b = NULL;
  arho1 = arho2 = arho3 = arho3b = t_ave = tsq_ave = NULL;

  maxneigh = 0;
  scrfcn = dscrfcn = fcpair = NULL;

  nelements = 0;
  elements = NULL;
  mass = NULL;

  // set comm size needed by this Pair

  comm_forward = 35;
  comm_reverse = 27;
}

/* ----------------------------------------------------------------------
   free all arrays
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairMEAM::~PairMEAM()
{
  meam_cleanup_();

  memory->sfree(rho);
  memory->sfree(rho0);
  memory->sfree(rho1);
  memory->sfree(rho2);
  memory->sfree(rho3);
  memory->sfree(frhop);
  memory->sfree(gamma);
  memory->sfree(dgamma1);
  memory->sfree(dgamma2);
  memory->sfree(dgamma3);
  memory->sfree(arho2b);

  memory->destroy_2d_double_array(arho1);
  memory->destroy_2d_double_array(arho2);
  memory->destroy_2d_double_array(arho3);
  memory->destroy_2d_double_array(arho3b);
  memory->destroy_2d_double_array(t_ave);
  memory->destroy_2d_double_array(tsq_ave);

  memory->sfree(scrfcn);
  memory->sfree(dscrfcn);
  memory->sfree(fcpair);
  
  for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  delete [] mass;

  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);
    delete [] map;
    delete [] fmap;
  }
}

/* ---------------------------------------------------------------------- */

void PairMEAM::compute(int eflag, int vflag)
{
  int i,j,ii,n,inum_half,itype,jtype,errorflag;
  double evdwl;
  int *ilist_half,*jlist_half,*numneigh_half,**firstneigh_half;
  int *numneigh_full,**firstneigh_full;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int newton_pair = force->newton_pair;

  // grow local arrays if necessary

  if (atom->nmax > nmax) {
    memory->sfree(rho);
    memory->sfree(rho0);
    memory->sfree(rho1);
    memory->sfree(rho2);
    memory->sfree(rho3);
    memory->sfree(frhop);
    memory->sfree(gamma);
    memory->sfree(dgamma1);
    memory->sfree(dgamma2);
    memory->sfree(dgamma3);
    memory->sfree(arho2b);
    memory->destroy_2d_double_array(arho1);
    memory->destroy_2d_double_array(arho2);
    memory->destroy_2d_double_array(arho3);
    memory->destroy_2d_double_array(arho3b);
    memory->destroy_2d_double_array(t_ave);
    memory->destroy_2d_double_array(tsq_ave);

    nmax = atom->nmax;

    rho = (double *) memory->smalloc(nmax*sizeof(double),"pair:rho");
    rho0 = (double *) memory->smalloc(nmax*sizeof(double),"pair:rho0");
    rho1 = (double *) memory->smalloc(nmax*sizeof(double),"pair:rho1");
    rho2 = (double *) memory->smalloc(nmax*sizeof(double),"pair:rho2");
    rho3 = (double *) memory->smalloc(nmax*sizeof(double),"pair:rho3");
    frhop = (double *) memory->smalloc(nmax*sizeof(double),"pair:frhop");
    gamma = (double *) memory->smalloc(nmax*sizeof(double),"pair:gamma");
    dgamma1 = (double *) memory->smalloc(nmax*sizeof(double),"pair:dgamma1");
    dgamma2 = (double *) memory->smalloc(nmax*sizeof(double),"pair:dgamma2");
    dgamma3 = (double *) memory->smalloc(nmax*sizeof(double),"pair:dgamma3");
    arho2b = (double *) memory->smalloc(nmax*sizeof(double),"pair:arho2b");
    arho1 = memory->create_2d_double_array(nmax,3,"pair:arho1");
    arho2 = memory->create_2d_double_array(nmax,6,"pair:arho2");
    arho3 = memory->create_2d_double_array(nmax,10,"pair:arho3");
    arho3b = memory->create_2d_double_array(nmax,3,"pair:arho3b");
    t_ave = memory->create_2d_double_array(nmax,3,"pair:t_ave");
    tsq_ave = memory->create_2d_double_array(nmax,3,"pair:tsq_ave");
  }

  // neighbor list info

  inum_half = listhalf->inum;
  ilist_half = listhalf->ilist;
  numneigh_half = listhalf->numneigh;
  firstneigh_half = listhalf->firstneigh;
  numneigh_full = listfull->numneigh;
  firstneigh_full = listfull->firstneigh;

  // check size of scrfcn based on half neighbor list

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  n = 0;
  for (ii = 0; ii < inum_half; ii++) n += numneigh_half[ilist_half[ii]];

  if (n > maxneigh) {
    memory->sfree(scrfcn);
    memory->sfree(dscrfcn);
    memory->sfree(fcpair);
    maxneigh = n;
    scrfcn =
      (double *) memory->smalloc(maxneigh*sizeof(double),"pair:scrfcn");
    dscrfcn = 
      (double *) memory->smalloc(maxneigh*sizeof(double),"pair:dscrfcn");
    fcpair = 
      (double *) memory->smalloc(maxneigh*sizeof(double),"pair:fcpair");
  }

  // zero out local arrays

  for (i = 0; i < nall; i++) {
    rho0[i] = 0.0;
    arho2b[i] = 0.0;
    arho1[i][0] = arho1[i][1] = arho1[i][2] = 0.0;
    for (j = 0; j < 6; j++) arho2[i][j] = 0.0;
    for (j = 0; j < 10; j++) arho3[i][j] = 0.0;
    arho3b[i][0] = arho3b[i][1] = arho3b[i][2] = 0.0;
    t_ave[i][0] = t_ave[i][1] = t_ave[i][2] = 0.0;
    tsq_ave[i][0] = tsq_ave[i][1] = tsq_ave[i][2] = 0.0;
  }

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int ntype = atom->ntypes;

  // change neighbor list indices to Fortran indexing

  neigh_c2f(inum_half,ilist_half,numneigh_half,firstneigh_half);
  neigh_c2f(inum_half,ilist_half,numneigh_full,firstneigh_full);

  // 3 stages of MEAM calculation
  // loop over my atoms followed by communication

  int ifort;
  int offset = 0;
  errorflag = 0;

  for (ii = 0; ii < inum_half; ii++) {
    i = ilist_half[ii];
    ifort = i+1;
    meam_dens_init_(&ifort,&nmax,&ntype,type,fmap,&x[0][0],
		    &numneigh_half[i],firstneigh_half[i],
		    &numneigh_full[i],firstneigh_full[i],
		    &scrfcn[offset],&dscrfcn[offset],&fcpair[offset],
		    rho0,&arho1[0][0],&arho2[0][0],arho2b,
		    &arho3[0][0],&arho3b[0][0],&t_ave[0][0],&tsq_ave[0][0],
                    &errorflag);
    if (errorflag) {
      char str[128];
      sprintf(str,"MEAM library error %d",errorflag);
      error->one(str);
    }
    offset += numneigh_half[i];
  }

  comm->reverse_comm_pair(this);

  meam_dens_final_(&nlocal,&nmax,&eflag_either,&eflag_global,&eflag_atom,
		   &eng_vdwl,eatom,&ntype,type,fmap,
		   &arho1[0][0],&arho2[0][0],arho2b,&arho3[0][0],
		   &arho3b[0][0],&t_ave[0][0],&tsq_ave[0][0],gamma,dgamma1,
		   dgamma2,dgamma3,rho,rho0,rho1,rho2,rho3,frhop,&errorflag);
  if (errorflag) {
    char str[128];
    sprintf(str,"MEAM library error %d",errorflag);
    error->one(str);
  }

  comm->comm_pair(this);

  offset = 0;

  // vptr is first value in vatom if it will be used by meam_force()
  // else vatom may not exist, so pass dummy ptr

  double *vptr;
  if (vflag_atom) vptr = &vatom[0][0];
  else vptr = &cutmax;

  for (ii = 0; ii < inum_half; ii++) {
    i = ilist_half[ii];
    ifort = i+1;
    meam_force_(&ifort,&nmax,&eflag_either,&eflag_global,&eflag_atom,
		&vflag_atom,&eng_vdwl,eatom,&ntype,type,fmap,&x[0][0],
		&numneigh_half[i],firstneigh_half[i],
		&numneigh_full[i],firstneigh_full[i],
		&scrfcn[offset],&dscrfcn[offset],&fcpair[offset],
		dgamma1,dgamma2,dgamma3,rho0,rho1,rho2,rho3,frhop,
		&arho1[0][0],&arho2[0][0],arho2b,&arho3[0][0],&arho3b[0][0],
		&t_ave[0][0],&tsq_ave[0][0],&f[0][0],vptr,&errorflag);
    if (errorflag) {
      char str[128];
      sprintf(str,"MEAM library error %d",errorflag);
      error->one(str);
    }
    offset += numneigh_half[i];
  }

  // change neighbor list indices back to C indexing

  neigh_f2c(inum_half,ilist_half,numneigh_half,firstneigh_half);
  neigh_f2c(inum_half,ilist_half,numneigh_full,firstneigh_full);

  if (vflag_fdotr) virial_compute();
}

/* ---------------------------------------------------------------------- */

void PairMEAM::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  map = new int[n+1];
  fmap = new int[n];
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairMEAM::settings(int narg, char **arg)
{
  if (narg != 0) error->all("Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMEAM::coeff(int narg, char **arg)
{
  int i,j,m,n;

  if (!allocated) allocate();

  if (narg < 6) error->all("Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all("Incorrect args for pair coefficients");

  // read MEAM element names between 2 filenames
  // nelements = # of MEAM elements
  // elements = list of unique element names

  if (nelements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
    delete [] mass;
  }
  nelements = narg - 4 - atom->ntypes;
  if (nelements < 1) error->all("Incorrect args for pair coefficients");
  elements = new char*[nelements];
  mass = new double[nelements];
  
  for (i = 0; i < nelements; i++) {
    n = strlen(arg[i+3]) + 1;
    elements[i] = new char[n];
    strcpy(elements[i],arg[i+3]);
  }

  // read MEAM library and parameter files
  // pass all parameters to MEAM package
  // tell MEAM package that setup is done

  read_files(arg[2],arg[2+nelements+1]);
  meam_setup_done_(&cutmax);
  
  // read args that map atom types to MEAM elements
  // map[i] = which element the Ith atom type is, -1 if not mapped

  for (i = 4 + nelements; i < narg; i++) {
    m = i - (4+nelements) + 1;
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    if (j < nelements) map[m] = j;
    else if (strcmp(arg[i],"NULL") == 0) map[m] = -1;
    else error->all("Incorrect args for pair coefficients");
  }

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  // set mass for i,i in atom class

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
	setflag[i][j] = 1;
	if (i == j) atom->set_mass(i,mass[map[i]]);
	count++;
      }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairMEAM::init_style()
{
  if (force->newton_pair == 0)
    error->all("Pair style MEAM requires newton pair on");

  // need full and half neighbor list

  int irequest_full = neighbor->request(this);
  neighbor->requests[irequest_full]->id = 1;
  neighbor->requests[irequest_full]->half = 0;
  neighbor->requests[irequest_full]->full = 1;
  int irequest_half = neighbor->request(this);
  neighbor->requests[irequest_half]->id = 2;
  neighbor->requests[irequest_half]->half = 0;
  neighbor->requests[irequest_half]->half_from_full = 1;
  neighbor->requests[irequest_half]->otherlist = irequest_full;

  // setup Fortran-style mapping array needed by MEAM package
  // fmap is indexed from 1:ntypes by Fortran and stores a Fortran index
  // if type I is not a MEAM atom, fmap stores a 0

  for (int i = 1; i <= atom->ntypes; i++) fmap[i-1] = map[i] + 1;
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   half or full
------------------------------------------------------------------------- */

void PairMEAM::init_list(int id, NeighList *ptr)
{
  if (id == 1) listfull = ptr;
  else if (id == 2) listhalf = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMEAM::init_one(int i, int j)
{
  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairMEAM::read_files(char *globalfile, char *userfile)
{
  // open global meamf file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = fopen(globalfile,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open MEAM potential file %s",globalfile);
      error->one(str);
    }
  }

  // allocate parameter arrays

  int params_per_line = 19;

  int *lat = new int[nelements];
  int *ielement = new int[nelements];
  int *ibar = new int[nelements];
  double *z = new double[nelements];
  double *atwt = new double[nelements];
  double *alpha = new double[nelements];
  double *b0 = new double[nelements];
  double *b1 = new double[nelements];
  double *b2 = new double[nelements];
  double *b3 = new double[nelements];
  double *alat = new double[nelements];
  double *esub = new double[nelements];
  double *asub = new double[nelements];
  double *t0 = new double[nelements];
  double *t1 = new double[nelements];
  double *t2 = new double[nelements];
  double *t3 = new double[nelements];
  double *rozero = new double[nelements];

  bool *found = new bool[nelements];
  for (int i = 0; i < nelements; i++) found[i] = false;

  // read each set of params from global MEAM file
  // one set of params can span multiple lines
  // store params if element name is in element list
  // if element name appears multiple times, only store 1st entry

  int i,n,nwords;
  char **words = new char*[params_per_line+1];
  char line[MAXLINE],*ptr;
  int eof = 0;

  int nset = 0;
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

    if (ptr = strchr(line,'#')) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
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
      if (ptr = strchr(line,'#')) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all("Incorrect format in MEAM potential file");

    // words = ptrs to all words in line
    // strip single and double quotes from words

    nwords = 0;
    words[nwords++] = strtok(line,"' \t\n\r\f");
    while (words[nwords++] = strtok(NULL,"' \t\n\r\f")) continue;

    // skip if element name isn't in element list

    for (i = 0; i < nelements; i++)
      if (strcmp(words[0],elements[i]) == 0) break;
    if (i == nelements) continue;

    // skip if element already appeared

    if (found[i] == true) continue;
    found[i] = true;

    // map lat string to an integer

    if (strcmp(words[1],"fcc") == 0) lat[i] = FCC;
    else if (strcmp(words[1],"bcc") == 0) lat[i] = BCC;
    else if (strcmp(words[1],"hcp") == 0) lat[i] = HCP;
    else if (strcmp(words[1],"dim") == 0) lat[i] = DIM;
    else if (strcmp(words[1],"dia") == 0) lat[i] = DIAMOND;
    else error->all("Unrecognized lattice type in MEAM file 1");

    // store parameters

    z[i] = atof(words[2]);
    ielement[i] = atoi(words[3]);
    atwt[i] = atof(words[4]);
    alpha[i] = atof(words[5]);
    b0[i] = atof(words[6]);
    b1[i] = atof(words[7]);
    b2[i] = atof(words[8]);
    b3[i] = atof(words[9]);
    alat[i] = atof(words[10]);
    esub[i] = atof(words[11]);
    asub[i] = atof(words[12]);
    t0[i] = atof(words[13]);
    t1[i] = atof(words[14]);
    t2[i] = atof(words[15]);
    t3[i] = atof(words[16]);
    rozero[i] = atof(words[17]);
    ibar[i] = atoi(words[18]);

    nset++;
  }

  // error if didn't find all elements in file

  if (nset != nelements)
    error->all("Did not find all elements in MEAM library file");

  // pass element parameters to MEAM package

  meam_setup_global_(&nelements,lat,z,ielement,atwt,alpha,b0,b1,b2,b3,
  		     alat,esub,asub,t0,t1,t2,t3,rozero,ibar);

  // set element masses

  for (i = 0; i < nelements; i++) mass[i] = atwt[i];

  // clean-up memory

  delete [] words;

  delete [] lat;
  delete [] ielement;
  delete [] ibar;
  delete [] z;
  delete [] atwt;
  delete [] alpha;
  delete [] b0;
  delete [] b1;
  delete [] b2;
  delete [] b3;
  delete [] alat;
  delete [] esub;
  delete [] asub;
  delete [] t0;
  delete [] t1;
  delete [] t2;
  delete [] t3;
  delete [] rozero;
  delete [] found;

  // done if user param file is NULL

  if (strcmp(userfile,"NULL") == 0) return;

  // open user param file on proc 0

  if (comm->me == 0) {
    fp = fopen(userfile,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open MEAM potential file %s",userfile);
      error->one(str);
    }
  }

  // read settings
  // pass them one at a time to MEAM package
  // match strings to list of corresponding ints

  int which;
  double value;
  int nindex,index[3];
  int maxparams = 6;
  char **params = new char*[maxparams];
  int nparams;

  eof = 0;
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

    if (ptr = strchr(line,'#')) *ptr = '\0';
    nparams = atom->count_words(line);
    if (nparams == 0) continue;

    // words = ptrs to all words in line

    nparams = 0;
    params[nparams++] = strtok(line,"=(), '\t\n\r\f");
    while (nparams < maxparams && 
	   (params[nparams++] = strtok(NULL,"=(), '\t\n\r\f")))
      continue;
    nparams--;

    for (which = 0; which < nkeywords; which++)
      if (strcmp(params[0],keywords[which]) == 0) break;
    if (which == nkeywords) {
      char str[128];
      sprintf(str,"Keyword %s in MEAM parameter file not recognized",
	      params[0]);
      error->all(str);
    }
    nindex = nparams - 2;
    for (i = 0; i < nindex; i++) index[i] = atoi(params[i+1]);

    // map lattce_meam value to an integer

    if (which == 4) {
      if (strcmp(params[nparams-1],"fcc") == 0) value = FCC;
      else if (strcmp(params[nparams-1],"bcc") == 0) value = BCC;
      else if (strcmp(params[nparams-1],"hcp") == 0) value = HCP;
      else if (strcmp(params[nparams-1],"dim") == 0) value = DIM;
      else if (strcmp(params[nparams-1],"dia") == 0) value = DIAMOND;
      else if (strcmp(params[nparams-1],"b1")  == 0) value = B1;
      else if (strcmp(params[nparams-1],"c11") == 0) value = C11;
      else if (strcmp(params[nparams-1],"l12") == 0) value = L12;
      else error->all("Unrecognized lattice type in MEAM file 2");
    }
    else value = atof(params[nparams-1]);

    // pass single setting to MEAM package

    int errorflag = 0;
    meam_setup_param_(&which,&value,&nindex,index,&errorflag);
    if (errorflag) {
      char str[128];
      sprintf(str,"MEAM library error %d",errorflag);
      error->all(str);
    }
  }

  delete [] params;
}

/* ---------------------------------------------------------------------- */

int PairMEAM::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,k,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = rho0[j];
    buf[m++] = rho1[j];
    buf[m++] = rho2[j];
    buf[m++] = rho3[j];
    buf[m++] = frhop[j];
    buf[m++] = gamma[j];
    buf[m++] = dgamma1[j];
    buf[m++] = dgamma2[j];
    buf[m++] = dgamma3[j];
    buf[m++] = arho2b[j];
    buf[m++] = arho1[j][0];
    buf[m++] = arho1[j][1];
    buf[m++] = arho1[j][2];
    buf[m++] = arho2[j][0];
    buf[m++] = arho2[j][1];
    buf[m++] = arho2[j][2];
    buf[m++] = arho2[j][3];
    buf[m++] = arho2[j][4];
    buf[m++] = arho2[j][5];
    for (k = 0; k < 10; k++) buf[m++] = arho3[j][k];
    buf[m++] = arho3b[j][0];
    buf[m++] = arho3b[j][1];
    buf[m++] = arho3b[j][2];
    buf[m++] = t_ave[j][0];
    buf[m++] = t_ave[j][1];
    buf[m++] = t_ave[j][2];
    buf[m++] = tsq_ave[j][0];
    buf[m++] = tsq_ave[j][1];
    buf[m++] = tsq_ave[j][2];
  }
  return 38;
}

/* ---------------------------------------------------------------------- */

void PairMEAM::unpack_comm(int n, int first, double *buf)
{
  int i,k,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    rho0[i] = buf[m++];
    rho1[i] = buf[m++];
    rho2[i] = buf[m++];
    rho3[i] = buf[m++];
    frhop[i] = buf[m++];
    gamma[i] = buf[m++];
    dgamma1[i] = buf[m++];
    dgamma2[i] = buf[m++];
    dgamma3[i] = buf[m++];
    arho2b[i] = buf[m++];
    arho1[i][0] = buf[m++];
    arho1[i][1] = buf[m++];
    arho1[i][2] = buf[m++];
    arho2[i][0] = buf[m++];
    arho2[i][1] = buf[m++];
    arho2[i][2] = buf[m++];
    arho2[i][3] = buf[m++];
    arho2[i][4] = buf[m++];
    arho2[i][5] = buf[m++];
    for (k = 0; k < 10; k++) arho3[i][k] = buf[m++];
    arho3b[i][0] = buf[m++];
    arho3b[i][1] = buf[m++];
    arho3b[i][2] = buf[m++];
    t_ave[i][0] = buf[m++];
    t_ave[i][1] = buf[m++];
    t_ave[i][2] = buf[m++];
    tsq_ave[i][0] = buf[m++];
    tsq_ave[i][1] = buf[m++];
    tsq_ave[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int PairMEAM::pack_reverse_comm(int n, int first, double *buf)
{
  int i,k,m,last,size;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = rho0[i];
    buf[m++] = arho2b[i];
    buf[m++] = arho1[i][0];
    buf[m++] = arho1[i][1];
    buf[m++] = arho1[i][2];
    buf[m++] = arho2[i][0];
    buf[m++] = arho2[i][1];
    buf[m++] = arho2[i][2];
    buf[m++] = arho2[i][3];
    buf[m++] = arho2[i][4];
    buf[m++] = arho2[i][5];
    for (k = 0; k < 10; k++) buf[m++] = arho3[i][k];
    buf[m++] = arho3b[i][0];
    buf[m++] = arho3b[i][1];
    buf[m++] = arho3b[i][2];
    buf[m++] = t_ave[i][0];
    buf[m++] = t_ave[i][1];
    buf[m++] = t_ave[i][2];
    buf[m++] = tsq_ave[i][0];
    buf[m++] = tsq_ave[i][1];
    buf[m++] = tsq_ave[i][2];
  }

  return 30;
}

/* ---------------------------------------------------------------------- */

void PairMEAM::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,k,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    rho0[j] += buf[m++];
    arho2b[j] += buf[m++];
    arho1[j][0] += buf[m++];
    arho1[j][1] += buf[m++];
    arho1[j][2] += buf[m++];
    arho2[j][0] += buf[m++];
    arho2[j][1] += buf[m++];
    arho2[j][2] += buf[m++];
    arho2[j][3] += buf[m++];
    arho2[j][4] += buf[m++];
    arho2[j][5] += buf[m++];
    for (k = 0; k < 10; k++) arho3[j][k] += buf[m++];
    arho3b[j][0] += buf[m++];
    arho3b[j][1] += buf[m++];
    arho3b[j][2] += buf[m++];
    t_ave[j][0] += buf[m++];
    t_ave[j][1] += buf[m++];
    t_ave[j][2] += buf[m++];
    tsq_ave[j][0] += buf[m++];
    tsq_ave[j][1] += buf[m++];
    tsq_ave[j][2] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

double PairMEAM::memory_usage()
{
  double bytes = 11 * nmax * sizeof(double);
  bytes += (3 + 6 + 10 + 3 + 3 + 3) * nmax * sizeof(double);
  bytes += 3 * maxneigh * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   toggle neighbor list indices between zero- and one-based values
   needed for access by MEAM Fortran library
------------------------------------------------------------------------- */

void PairMEAM::neigh_f2c(int inum, int *ilist, int *numneigh, int **firstneigh)
{
  int i,j,ii,jnum;
  int *jlist;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (j = 0; j < jnum; j++) jlist[j]--;
  }
}

void PairMEAM::neigh_c2f(int inum, int *ilist, int *numneigh, int **firstneigh)
{
  int i,j,ii,jnum;
  int *jlist;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (j = 0; j < jnum; j++) jlist[j]++;
  }
}
