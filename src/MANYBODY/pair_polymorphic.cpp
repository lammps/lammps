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
   Contributing authors: Xiaowang Zhou, Reese Jones (SNL)
   Based on pair_tersoff by Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_polymorphic.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 1024
#define DELTA 4

/* ====================================================================== */

PairPolymorphic::PairPolymorphic(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  one_coeff = 1;

  nelements = 0;
  elements = NULL;
  pairParameters = NULL;
  tripletParameters = NULL;
  elem2param = NULL;
  elem3param = NULL;
  type_map = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairPolymorphic::~PairPolymorphic()
{
  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  memory->destroy(pairParameters);
  memory->destroy(tripletParameters);
  memory->destroy(elem2param);
  memory->destroy(elem3param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] type_map;
  }
}

/* ---------------------------------------------------------------------- */

void PairPolymorphic::compute(int eflag, int vflag)
{
  tagint itag,jtag;
  int i,j,k,ii,jj,kk,inum,jnum;
  int iel,jel,kel,iparam_ij,iparam_ik,iparam_ijk;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2,r0,r1,r2;
  double delr1[3],delr2[3],fi[3],fj[3],fk[3];
  double zeta_ij,prefactor,wfac,pfac,gfac,fa,fa_d,bij,bij_d;
  double costheta;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = vflag_atom = 0;

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    iel = type_map[type[i]];
    if (iel < 0) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // two-body interactions, skip half of them

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtag = tag[j];

      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < x[i][2]) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }


      jel = type_map[type[j]];
      if (jel < 0) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      iparam_ij = elem2param[iel][jel];
      PairParameters & p = pairParameters[iparam_ij];
      if (rsq > p.cutsq) continue;
      r0 = sqrt(rsq);
      if (eflag) evdwl = (p.U)->value(r0);
      fpair = (p.U)->derivative(r0);
      fpair = -fpair/r0;

      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
        evdwl,0.0,fpair,delx,dely,delz);
    }

    if (eta) {

      iparam_ij = elem2param[iel][iel];
      PairParameters & p = pairParameters[iparam_ij];

      // accumulate bondorder zeta for each i-j interaction via loop over k

      zeta_ij = 0.0;

      for (kk = 0; kk < jnum; kk++) {
        k = jlist[kk];
        k &= NEIGHMASK;
        kel = type_map[type[k]];
        if (kel < 0) continue;

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

        iparam_ik = elem2param[kel][kel];
        PairParameters & q = pairParameters[iparam_ik];
        if (rsq2 > q.cutsq) continue;
        r2 = sqrt(rsq2);

        wfac = (q.W)->value(r2);

        zeta_ij += wfac;
      }

      // pairwise force due to zeta

      bij   = (p.F)->value(zeta_ij);
      bij_d = (p.F)->derivative(zeta_ij);

      prefactor = 0.5* bij_d;
      if (eflag) evdwl = -0.5*bij;

      if (evflag) ev_tally(i,i,nlocal,newton_pair,
                           evdwl,0.0,0.0,-delr1[0],-delr1[1],-delr1[2]);

      // attractive term via loop over k

      for (kk = 0; kk < jnum; kk++) {
        k = jlist[kk];
        k &= NEIGHMASK;
        kel = type_map[type[k]];
        if (kel < 0) continue;

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

        iparam_ik = elem2param[kel][kel];
        PairParameters & q = pairParameters[iparam_ik];
        if (rsq2 > q.cutsq) continue;
        r2 = sqrt(rsq2);

        fpair = (q.W)->derivative(r2);
        fpair = -prefactor*fpair/r2;

        f[i][0] += delr2[0]*fpair;
        f[i][1] += delr2[1]*fpair;
        f[i][2] += delr2[2]*fpair;
        f[k][0] -= delr2[0]*fpair;
        f[k][1] -= delr2[1]*fpair;
        f[k][2] -= delr2[2]*fpair;

        if (vflag_atom) v_tally2(i, k, -fpair, delr2);
      }

    } else {

      // three-body interactions
      // skip immediately if I-J is not within cutoff

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        jel = type_map[type[j]];
        if (jel < 0) continue;

        delr1[0] = x[j][0] - xtmp;
        delr1[1] = x[j][1] - ytmp;
        delr1[2] = x[j][2] - ztmp;
        rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

        iparam_ij = elem2param[iel][jel];
        PairParameters & p = pairParameters[iparam_ij];
        if (rsq1 > p.cutsq) continue;
        r1 = sqrt(rsq1);

        // accumulate bondorder zeta for each i-j interaction via loop over k

        zeta_ij = 0.0;

        for (kk = 0; kk < jnum; kk++) {
          if (jj == kk) continue;
          k = jlist[kk];
          k &= NEIGHMASK;
          kel = type_map[type[k]];
          if (kel < 0) continue;

          delr2[0] = x[k][0] - xtmp;
          delr2[1] = x[k][1] - ytmp;
          delr2[2] = x[k][2] - ztmp;
          rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

          iparam_ik = elem2param[iel][kel];
          PairParameters & q = pairParameters[iparam_ik];
          if (rsq2 > q.cutsq) continue;
          r2 = sqrt(rsq2);

          costheta = (delr1[0]*delr2[0] + delr1[1]*delr2[1] +
                      delr1[2]*delr2[2]) / (r1*r2);
          iparam_ijk = elem3param[jel][iel][kel]; 
          TripletParameters & trip = tripletParameters[iparam_ijk];

          wfac= (q.W)->value(r2);
          pfac= (q.P)->value(r1-(p.xi)*r2);
          gfac= (trip.G)->value(costheta);

          zeta_ij += wfac*pfac*gfac;
        }

        // pairwise force due to zeta

        fa   = (p.V)->value(r1);
        fa_d = (p.V)->derivative(r1);
        bij   = (p.F)->value(zeta_ij);
        bij_d = (p.F)->derivative(zeta_ij);
        fpair = -0.5*bij*fa_d / r1;
        prefactor = 0.5* fa * bij_d;
        if (eflag) evdwl = -0.5*bij*fa;

        f[i][0] += delr1[0]*fpair;
        f[i][1] += delr1[1]*fpair;
        f[i][2] += delr1[2]*fpair;
        f[j][0] -= delr1[0]*fpair;
        f[j][1] -= delr1[1]*fpair;
        f[j][2] -= delr1[2]*fpair;

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,-fpair,-delr1[0],-delr1[1],-delr1[2]);

        // attractive term via loop over k

        for (kk = 0; kk < jnum; kk++) {
          if (jj == kk) continue;
          k = jlist[kk];
          k &= NEIGHMASK;
          kel = type_map[type[k]];
          if (kel < 0) continue;

          delr2[0] = x[k][0] - xtmp;
          delr2[1] = x[k][1] - ytmp;
          delr2[2] = x[k][2] - ztmp;
          rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

          iparam_ik = elem2param[iel][kel];
          PairParameters & q = pairParameters[iparam_ik];
          if (rsq2 > q.cutsq) continue;
          r2 = sqrt(rsq2);

          iparam_ijk = elem3param[jel][iel][kel];
          TripletParameters & trip = tripletParameters[iparam_ijk];
 
          attractive(&q,&trip,prefactor,r1,r2,delr1,delr2,fi,fj,fk);

          f[i][0] += fi[0];
          f[i][1] += fi[1];
          f[i][2] += fi[2];
          f[j][0] += fj[0];
          f[j][1] += fj[1];
          f[j][2] += fj[2];
          f[k][0] += fk[0];
          f[k][1] += fk[1];
          f[k][2] += fk[2];

          if (vflag_atom) v_tally3(i,j,k,fj,fk,delr1,delr2);
        }
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairPolymorphic::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  type_map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairPolymorphic::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairPolymorphic::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style polymorphic requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style polymorphic requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPolymorphic::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}


/* ---------------------------------------------------------------------- */

void PairPolymorphic::setup()
{
  int i,j,k,n;

  memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,"pair:elem2param");
  memory->destroy(elem3param);
  memory->create(elem3param,nelements,nelements,nelements,"pair:elem3param");

  // map atom pair to parameter index, as read from potential file
  n = 0;
  for (i = 0; i < nelements; i++) { // note self first
    elem2param[i][i] = n;
    n++;
  }
  for (i = 0; i < nelements; i++)
  for (j = i+1; j < nelements; j++) {
    elem2param[i][j] = n;
    elem2param[j][i] = n;
    n++;
  }

  // map atom triplet to parameter index
  n = 0;
  for (i = 0; i < nelements; i++)
  for (j = 0; j < nelements; j++)
  for (k = 0; k < nelements; k++) {
    elem3param[i][j][k] = n;
    n++;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairPolymorphic::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *
  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that type_map atom types to elements in potential file
  int ntypes = atom->ntypes;
  // type_map = atom type to element in potential file
  if (type_map) { delete [] type_map; }
  type_map = new int[ntypes+1];
  for (int i = 0; i < ntypes+1; i++) { 
    type_map[i] = -1; 
  }
  // elements = list of requested element names (ntypes long)
  char** elements = new char*[ntypes];
  for (int i = 0; i < ntypes; i++) { elements[i] = NULL; }
  // parse and store
  for (int i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") != 0) {
      int n = strlen(arg[i]) + 1;
      elements[i-3] = new char[n];
      strcpy(elements[i-3],arg[i]);
    }
  }

  // read potential file and initialize potential parameters
  read_file(arg[2],elements);
  setup();

  if (elements)
    for (int i = 0; i < ntypes; i++) 
      if (elements[i]) delete [] elements[i];
  delete [] elements;

  // clear setflag since coeff() called once with I,J = * *
  int n = atom->ntypes;
  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) {
      setflag[i][j] = 0;
    }
  }

  // set setflag i,j for type pairs where both are type_mapped to elements
  int count = 0;
  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) {
      if ((type_map[i] > -1)  && (type_map[j] > -1)) {
        setflag[i][j] = 1;
        count++;
      }
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}
/* ---------------------------------------------------------------------- */

void PairPolymorphic::read_file(char *file, char** elements)
{
  char line[MAXLINE],*ptr, ftype[MAXLINE];
  int n;
  // open file on proc 0
  FILE *fp=NULL;
  if (comm->me == 0) {
    fp = force->open_potential(file);
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open polymorphic potential file %s",file);
      error->one(FLERR,str);
    }
    // move past comments to first data line
    fgets(line,MAXLINE,fp);
    while (line == strchr(line,'#')) fgets(line,MAXLINE,fp);
    n = strlen(line) + 1;
  }
  MPI_Bcast(&n,1,MPI_INT,0,world);
  MPI_Bcast(line,n,MPI_CHAR,0,world);
  ptr = strtok(line," \t\n\r\f"); // 1st line, 1st token : nelements
  nelements = atoi(ptr); // number of elements in potential file
  ptr = strtok(NULL," \t\n\r\f"); // 1st line, 2nd token : indicator eta
  eta = (atoi(ptr)>0) ? true:false;
  if (comm->me == 0) { printf("%d elements in: %s,",nelements,file); }

  // type_map the elements in the potential file to LAMMPS atom types
  for (int i = 0; i < nelements; i++) {
    read_line(fp,line);
    ptr = strtok(line," \t\n\r\f"); // 1st token atomic number
    ptr = strtok(NULL," \t\n\r\f"); // 2st token atomic mass
    ptr = strtok(NULL," \t\n\r\f"); // 3st token atomic symbol
    if (comm->me == 0) { printf(" %s",ptr); }
    int j = 0;
    for (j = 0; j < atom->ntypes; j++) {
      if (elements[j] && strcmp(ptr,elements[j]) == 0) {
        type_map[j+1] = i;
        if (comm->me == 0) { printf("=%d ",j+1); }
        break;
      }
    }
    if (j == nelements) 
      error->all(FLERR,"Element not defined in potential file");
  }
  if (comm->me == 0) { printf("\n"); }

  // size
  npair = nelements*(nelements+1)/2;
  ntriple = nelements*nelements*nelements;
  pairParameters = (PairParameters*)
    memory->srealloc(pairParameters,npair*sizeof(PairParameters),
    "pair:pairParameters");
  tripletParameters = (TripletParameters*)
    memory->srealloc(tripletParameters,ntriple*sizeof(TripletParameters),
    "pair:tripletParameters");

  // pairwise cutoffs
  for (int i = 0; i < npair; i++) {
    PairParameters & p = pairParameters[i];
    read_line(fp,line);
    ptr = strtok(line," \t\n\r\f"); // 1st token: cutoff
    p.cut = atof(ptr);
    p.cutsq = p.cut*p.cut;
    ptr = strtok(NULL," \t\n\r\f"); // 2nd token: indicator xi
    p.xi = (atoi(ptr)>0) ? true:false;
  }

  // set cutmax to max of all params
  cutmax = 0.0;
  for (int i = 0; i < npair; i++) {
    PairParameters & p = pairParameters[i];
    if (p.cut > cutmax) cutmax = p.cut;
  }

  // start reading functions
  for (int i = 0; i < npair; i++) { // U
    PairParameters & p = pairParameters[i];
    read_line(fp,line);
    ptr = strtok(line," \t\n\r\f"); // 1st token
    strcpy(ftype,ptr);
    p.U = create_function(ftype,fp); 
  }
  for (int i = 0; i < npair; i++) { // V
    PairParameters & p = pairParameters[i];
    read_line(fp,line);
    ptr = strtok(line," \t\n\r\f"); // 1st token
    strcpy(ftype,ptr);
    p.V = create_function(ftype,fp); 
  }
  for (int i = 0; i < npair; i++) { // W
    PairParameters & p = pairParameters[i];
    read_line(fp,line);
    ptr = strtok(line," \t\n\r\f"); // 1st token
    strcpy(ftype,ptr);
    p.W = create_function(ftype,fp); 
  }
  for (int i = 0; i < npair; i++) { // P
    PairParameters & p = pairParameters[i];
    read_line(fp,line);
    ptr = strtok(line," \t\n\r\f"); // 1st token
    strcpy(ftype,ptr);
    p.P = create_function(ftype,fp); 
  }
  for (int i = 0; i < ntriple; i++) { // G
    TripletParameters & p = tripletParameters[i];
    read_line(fp,line);
    ptr = strtok(line," \t\n\r\f"); // 1st token
    strcpy(ftype,ptr);
    p.G = create_function(ftype,fp); 
  }
  for (int i = 0; i < npair; i++) { // F
    PairParameters & p = pairParameters[i];
    read_line(fp,line);
    ptr = strtok(line," \t\n\r\f"); // 1st token
    strcpy(ftype,ptr);
    p.F = create_function(ftype,fp); 
  }
  if (comm->me == 0) { fclose(fp); }
}

/* ---------------------------------------------------------------------- */

C1function * PairPolymorphic::create_function(char* ftype, FILE* fp)
{
  char * ptr;
  if (strcmp(ftype,"spline")==0) { // N, min, max, values
    C1tabularFunction * f = new C1tabularFunction();
    ptr = strtok(NULL," \t\n\r\f"); 
    int n = atof(ptr);
    ptr = strtok(NULL," \t\n\r\f"); 
    double xmin = atof(ptr);
    ptr = strtok(NULL," \t\n\r\f"); 
    double xmax = atof(ptr);
    double * table = new double[n];
    read_array(fp,n,table);
    f->set_values(n,xmin,xmax,table);
    delete [] table;
    return f;
  }
  else if (strcmp(ftype,"constant") == 0) {
    ptr = strtok(NULL," \t\n\r\f"); 
    double c = atof(ptr);
    return new C1constant(c);
  }
  else if (strcmp(ftype,"exponential") == 0) {
    ptr = strtok(NULL," \t\n\r\f"); 
    double c = atof(ptr);
    ptr = strtok(NULL," \t\n\r\f"); 
    double lambda = atof(ptr);
    return new C1exponential(c,lambda);
  }
  else if (strcmp(ftype,"sine") == 0) {
    ptr = strtok(NULL," \t\n\r\f"); 
    double c = atof(ptr);
    ptr = strtok(NULL," \t\n\r\f"); 
    double w = atof(ptr);
    return new C1sine(c,w);
  }
  else if (strcmp(ftype,"cosine") == 0) {
    ptr = strtok(NULL," \t\n\r\f"); 
    double c = atof(ptr);
    ptr = strtok(NULL," \t\n\r\f"); 
    double w = atof(ptr);
    return new C1cosine(c,w);
  }
  else { error->all(FLERR,"unknown function type"); }
  return NULL;
}
/* ---------------------------------------------------------------------- */

void PairPolymorphic::read_line(FILE *fp, char *line)
{
  int n = 0;
  if (comm->me == 0) {
    fgets(line,MAXLINE,fp);
    n = strlen(line) + 1;
  }
  MPI_Bcast(&n,1,MPI_INT,0,world);
  MPI_Bcast(line,n,MPI_CHAR,0,world);
}
void PairPolymorphic::read_array(FILE *fp, int n, double *list)
{
  if (comm->me == 0) {
    char *ptr;
    char line[MAXLINE];
    int i = 0;
    while (i < n) {
      fgets(line,MAXLINE,fp);
      ptr = strtok(line," \t\n\r\f");
      list[i++] = atof(ptr);
      while ((ptr = strtok(NULL," \t\n\r\f"))) list[i++] = atof(ptr);
    }
  }
  MPI_Bcast(list,n,MPI_DOUBLE,0,world);
}


/* ----------------------------------------------------------------------
   attractive term
------------------------------------------------------------------------- */

void PairPolymorphic::attractive(PairParameters *p, TripletParameters *trip,
                            double prefactor, double rij, double rik,
                            double *delrij, double *delrik,
                            double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rijinv,rikinv;

  rijinv = 1.0/rij;
  vec3_scale(rijinv,delrij,rij_hat);
  
  rikinv = 1.0/rik;
  vec3_scale(rikinv,delrik,rik_hat);

  ters_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk,p,trip);
}

/* ---------------------------------------------------------------------- */

void PairPolymorphic::ters_zetaterm_d(double prefactor,
                                 double *rij_hat, double rij,
                                 double *rik_hat, double rik,
                                 double *dri, double *drj, double *drk,
                                 PairParameters *p, TripletParameters *trip)
{
  double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

  cos_theta = vec3_dot(rij_hat,rik_hat);

  fc  = (p->W)->value(rik);
  dfc = (p->W)->derivative(rik);
  ex_delr   = (p->P)->value(rij-(p->xi)*rik);
  ex_delr_d = (p->P)->derivative(rij-(p->xi)*rik);
  gijk   = (trip->G)->value(cos_theta);
  gijk_d = (trip->G)->derivative(cos_theta);

  costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

  // compute the derivative wrt Ri
  // dri = -dfc*gijk*ex_delr*rik_hat;
  // dri += fc*gijk_d*ex_delr*dcosdri;
  // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);

  vec3_scale(-dfc*gijk*ex_delr,rik_hat,dri);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdri,dri,dri);
  vec3_scaleadd(fc*gijk*ex_delr_d,rik_hat,dri,dri);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rij_hat,dri,dri);
  vec3_scale(prefactor,dri,dri);

  // compute the derivative wrt Rj
  // drj = fc*gijk_d*ex_delr*dcosdrj;
  // drj += fc*gijk*ex_delr_d*rij_hat;

  vec3_scale(fc*gijk_d*ex_delr,dcosdrj,drj);
  vec3_scaleadd(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  vec3_scale(prefactor,drj,drj);

  // compute the derivative wrt Rk
  // drk = dfc*gijk*ex_delr*rik_hat;
  // drk += fc*gijk_d*ex_delr*dcosdrk;
  // drk += -fc*gijk*ex_delr_d*rik_hat;

  vec3_scale(dfc*gijk*ex_delr,rik_hat,drk);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdrk,drk,drk);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rik_hat,drk,drk);
  vec3_scale(prefactor,drk,drk);
}

/* ---------------------------------------------------------------------- */

void PairPolymorphic::costheta_d(double *rij_hat, double rij,
        double *rik_hat, double rik,
        double *dri, double *drj, double *drk)
{
  // first element is devative wrt Ri, second wrt Rj, third wrt Rk

  double cos_theta = vec3_dot(rij_hat,rik_hat);

  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,drj);
  vec3_scale(1.0/rij,drj,drj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,drk);
  vec3_scale(1.0/rik,drk,drk);
  vec3_add(drj,drk,dri);
  vec3_scale(-1.0,dri,dri);
}

