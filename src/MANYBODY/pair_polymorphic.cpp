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
   Contributing author: Reese Jones, Xiaowang Zhou (SNL)
   This modifies from pair_tersoff.cpp by Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "pair_polymorphic.h"
#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

PairPolymorphic::PairPolymorphic(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;

  nelements = 0;
  elements = NULL;
  pairParameters = NULL;
  tripletParameters = NULL;
  elem2param = NULL;
  elem3param = NULL;
  map = NULL;
  epsilon = 0.0;
  neighsize = 0;
  firstneighV = NULL;
  firstneighW = NULL;
  firstneighW1 = NULL;
  delxV = NULL;
  delyV = NULL;
  delzV = NULL;
  drV = NULL;
  delxW = NULL;
  delyW = NULL;
  delzW = NULL;
  drW = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairPolymorphic::~PairPolymorphic()
{
  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  delete [] match;
  memory->destroy(pairParameters);
  memory->destroy(tripletParameters);
  memory->destroy(elem2param);
  memory->destroy(elem3param);
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
    delete [] firstneighV;
    delete [] firstneighW;
    delete [] firstneighW1;
    delete [] delxV;
    delete [] delyV;
    delete [] delzV;
    delete [] drV;
    delete [] delxW;
    delete [] delyW;
    delete [] delzW;
    delete [] drW;
  }
}

/* ---------------------------------------------------------------------- */

void PairPolymorphic::compute(int eflag, int vflag)
{
  tagint itag,jtag;
  int i,j,k,ii,jj,kk,kk1,inum,jnum;
  int itype,jtype,ktype;
  int iparam_ii,iparam_jj,iparam_kk,iparam_ij,iparam_ik,iparam_ijk;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r0,r1,r2;
  double delr1[3],delr2[3],fi[3],fj[3],fk[3];
  double zeta_ij,prefactor,wfac,pfac,gfac,fa,fa_d,bij,bij_d;
  double costheta;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double emb;

  evdwl = 0.0;
  ev_init(eflag,vflag);

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
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (neighsize < jnum) {
      delete [] firstneighV;
      delete [] delxV;
      delete [] delyV;
      delete [] delzV;
      delete [] drV;
      delete [] firstneighW;
      delete [] delxW;
      delete [] delyW;
      delete [] delzW;
      delete [] drW;
      delete [] firstneighW1;
      neighsize = jnum + 20;
      firstneighV = new int[neighsize];
      delxV = new double[neighsize];
      delyV = new double[neighsize];
      delzV = new double[neighsize];
      drV = new double[neighsize];
      firstneighW = new int[neighsize];
      delxW = new double[neighsize];
      delyW = new double[neighsize];
      delzW = new double[neighsize];
      drW = new double[neighsize];
      firstneighW1 = new int[neighsize];
    }

    if (eta) {
      iparam_ii = elem2param[itype][itype];
      PairParameters & p = pairParameters[iparam_ii];
      emb = (p.F)->get_vmax();
    }

    numneighV = -1;
    numneighW = -1;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq >= cutmaxsq) continue;
      r0 = sqrt(rsq);

      iparam_ij = elem2param[itype][jtype];
      PairParameters & p = pairParameters[iparam_ij];

// do not include the neighbor if get_vmax() <= epsilon because the function is near zero
      if (eta) {
        if (emb > epsilon) {
          iparam_jj = elem2param[jtype][jtype];
          PairParameters & q = pairParameters[iparam_jj];
          if (rsq < (q.W)->get_xmaxsq() && (q.W)->get_vmax() > epsilon) {
            numneighW = numneighW + 1;
            firstneighW[numneighW] = j;
            delxW[numneighW] = delx;
            delyW[numneighW] = dely;
            delzW[numneighW] = delz;
            drW[numneighW] =  r0;
          }
        }
      } else {
        if ((p.F)->get_vmax() > epsilon) {
          if (rsq < (p.V)->get_xmaxsq() && (p.V)->get_vmax() > epsilon) {
            numneighV = numneighV + 1;
            firstneighV[numneighV] = j;
            delxV[numneighV] = delx;
            delyV[numneighV] = dely;
            delzV[numneighV] = delz;
            drV[numneighV] =  r0;
          }
          if (rsq < (p.W)->get_xmaxsq() && (p.W)->get_vmax() > epsilon) {
            numneighW = numneighW + 1;
            firstneighW[numneighW] = j;
            delxW[numneighW] = delx;
            delyW[numneighW] = dely;
            delzW[numneighW] = delz;
            drW[numneighW] =  r0;
          }
        }
      }

    // two-body interactions, skip half of them

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

      if (rsq >= (p.U)->get_xmaxsq() || (p.U)->get_vmax() <= epsilon) continue;
      (p.U)->value(r0,evdwl,eflag,fpair,1);
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

      if (emb > epsilon) {

        iparam_ii = elem2param[itype][itype];
        PairParameters & p = pairParameters[iparam_ii];

        // accumulate bondorder zeta for each i-j interaction via loop over k

        zeta_ij = 0.0;

        for (kk = 0; kk <= numneighW; kk++) {
          k = firstneighW[kk];
          ktype = map[type[k]];

          iparam_kk = elem2param[ktype][ktype];
          PairParameters & q = pairParameters[iparam_kk];

          (q.W)->value(drW[kk],wfac,1,fpair,0);

          zeta_ij += wfac;
        }

        // pairwise force due to zeta

        (p.F)->value(zeta_ij,bij,1,bij_d,1);

        prefactor = 0.5* bij_d;
        if (eflag) evdwl = -0.5*bij;

        if (evflag) ev_tally(i,i,nlocal,newton_pair,evdwl,0.0,0.0,delx,dely,delz);

        // attractive term via loop over k

        for (kk = 0; kk <= numneighW; kk++) {
          k = firstneighW[kk];
          ktype = map[type[k]];

          delr2[0] = -delxW[kk];
          delr2[1] = -delyW[kk];
          delr2[2] = -delzW[kk];

          iparam_kk = elem2param[ktype][ktype];
          PairParameters & q = pairParameters[iparam_kk];

          (q.W)->value(drW[kk],wfac,0,fpair,1);
          fpair = -prefactor*fpair/drW[kk];

          f[i][0] += delr2[0]*fpair;
          f[i][1] += delr2[1]*fpair;
          f[i][2] += delr2[2]*fpair;
          f[k][0] -= delr2[0]*fpair;
          f[k][1] -= delr2[1]*fpair;
          f[k][2] -= delr2[2]*fpair;

          if (vflag_atom) v_tally2(i, k, -fpair, delr2);
        }
      }

    } else {

      for (jj = 0; jj <= numneighV; jj++) {
        j = firstneighV[jj];
        jtype = map[type[j]];

        iparam_ij = elem2param[itype][jtype];
        PairParameters & p = pairParameters[iparam_ij];

        delr1[0] = -delxV[jj];
        delr1[1] = -delyV[jj];
        delr1[2] = -delzV[jj];
        r1 = drV[jj];

        // accumulate bondorder zeta for each i-j interaction via loop over k

        zeta_ij = 0.0;

        numneighW1 = -1;
        for (kk = 0; kk <= numneighW; kk++) {
          k = firstneighW[kk];
          if (j == k) continue;
          ktype = map[type[k]];
          iparam_ijk = elem3param[jtype][itype][ktype];
          TripletParameters & trip = tripletParameters[iparam_ijk];
          if ((trip.G)->get_vmax() <= epsilon) continue;

          numneighW1 = numneighW1 + 1;
          firstneighW1[numneighW1] = kk;

          delr2[0] = -delxW[kk];
          delr2[1] = -delyW[kk];
          delr2[2] = -delzW[kk];
          r2 = drW[kk];

          costheta = (delr1[0]*delr2[0] + delr1[1]*delr2[1] +
                      delr1[2]*delr2[2]) / (r1*r2);

          iparam_ik = elem2param[itype][ktype];
          PairParameters & q = pairParameters[iparam_ik];

          (q.W)->value(r2,wfac,1,fpair,0);
          (q.P)->value(r1-(p.xi)*r2,pfac,1,fpair,0);
          (trip.G)->value(costheta,gfac,1,fpair,0);

          zeta_ij += wfac*pfac*gfac;
        }

        // pairwise force due to zeta

        (p.V)->value(r1,fa,1,fa_d,1);
        (p.F)->value(zeta_ij,bij,1,bij_d,1);
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

        for (kk1 = 0; kk1 <= numneighW1; kk1++) {
          kk = firstneighW1[kk1];
          k = firstneighW[kk];
          ktype = map[type[k]];
          iparam_ijk = elem3param[jtype][itype][ktype];
          TripletParameters & trip = tripletParameters[iparam_ijk];

          delr2[0] = -delxW[kk];
          delr2[1] = -delyW[kk];
          delr2[2] = -delzW[kk];
          r2 = drW[kk];

          iparam_ik = elem2param[itype][ktype];
          PairParameters & q = pairParameters[iparam_ik];

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

  map = new int[n+1];

  neighsize = 40;
  firstneighV = new int[neighsize];
  delxV = new double[neighsize];
  delyV = new double[neighsize];
  delzV = new double[neighsize];
  drV = new double[neighsize];
  firstneighW = new int[neighsize];
  delxW = new double[neighsize];
  delyW = new double[neighsize];
  delzW = new double[neighsize];
  drW = new double[neighsize];
  firstneighW1 = new int[neighsize];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairPolymorphic::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairPolymorphic::coeff(int narg, char **arg)
{
  int i,j,n;

  if (!allocated) allocate();

  if (narg == 4 + atom->ntypes) {
     narg--;
     epsilon = atof(arg[narg]);
  } else if (narg != 3 + atom->ntypes) {
    error->all(FLERR,"Incorrect args for pair coefficients");
  }

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) elements[i] = NULL;

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup_params();

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
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

void PairPolymorphic::read_file(char *file)
{
  char line[MAXLINE],*ptr;
  int n;
  char *r_token;

  // open file on proc 0
  FILE *fp=NULL;
  if (comm->me == 0) {
    fp = force->open_potential(file);
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open polymorphic potential file %s",file);
      error->one(FLERR,str);
    }
    // move past comments to first data line
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
    while (line == strchr(line,'#')) utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
    n = strlen(line) + 1;
  }
  MPI_Bcast(&n,1,MPI_INT,0,world);
  MPI_Bcast(line,n,MPI_CHAR,0,world);
  r_token = line;
  ptr = utils::strtok_r(r_token," \t\n\r\f",&r_token); // 1st line, 1st token
  int ntypes = atoi(ptr);
  if (ntypes != nelements)
    error->all(FLERR,"Incorrect number of elements in potential file");
  match = new int[nelements];
  ptr = utils::strtok_r(NULL," \t\n\r\f",&r_token); // 1st line, 2nd token
  eta = (atoi(ptr)>0) ? true:false;

  // map the elements in the potential file to LAMMPS atom types
  for (int i = 0; i < nelements; i++) {
    if (comm->me == 0) {
      utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
      n = strlen(line) + 1;
    }
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);
    r_token = line;
    ptr = utils::strtok_r(r_token," \t\n\r\f",&r_token); // 1st token
    ptr = utils::strtok_r(NULL," \t\n\r\f",&r_token); // 2st token
    ptr = utils::strtok_r(NULL," \t\n\r\f",&r_token); // 3st token
    int j;
    for (j = 0; j < nelements; j++) {
      if (strcmp(ptr,elements[j]) == 0) break;
    }
    if (j == nelements)
      error->all(FLERR,"Element not defined in potential file");
    match[i] = j;
  }
  // sizes
  if (comm->me == 0) {
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
    n = strlen(line) + 1;
  }

  // Note: the format of this line has changed between the
  // 2015-06-06 and 2015-12-09 versions of the pair style.

  MPI_Bcast(&n,1,MPI_INT,0,world);
  MPI_Bcast(line,n,MPI_CHAR,0,world);
  r_token = line;
  nr = ng = nx = 0;
  ptr = utils::strtok_r(r_token," \t\n\r\f",&r_token); // 1st token
  if (ptr) nr = atoi(ptr);
  ptr = utils::strtok_r(NULL," \t\n\r\f",&r_token); // 2nd token
  if (ptr) ng = atoi(ptr);
  ptr = utils::strtok_r(NULL," \t\n\r\f",&r_token); // 3rd token
  if (ptr) nx = atoi(ptr);
  ptr = utils::strtok_r(NULL," \t\n\r\f",&r_token); // 4th token
  if (ptr) maxX = atof(ptr);
  if (ptr == NULL)
    error->all(FLERR,"Potential file incompatible with this pair style version");
  if ((ng == 0) || (nr == 0) || (nx == 0))
    error->all(FLERR,"Error reading potential file header");

  npair = nelements*(nelements+1)/2;
  ntriple = nelements*nelements*nelements;
  pairParameters = (PairParameters*)
    memory->srealloc(pairParameters,npair*sizeof(PairParameters),
    "pair:pairParameters");
  tripletParameters = (TripletParameters*)
    memory->srealloc(tripletParameters,ntriple*sizeof(TripletParameters),
    "pair:tripletParameters");

  // cutoffs
  for (int i = 0; i < npair; i++) {
    PairParameters & p = pairParameters[i];
    if (comm->me == 0) {
      utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
      n = strlen(line) + 1;
    }
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);
    r_token = line;
    ptr = utils::strtok_r(r_token," \t\n\r\f",&r_token); // 1st token
    p.cut = atof(ptr);
    p.cutsq = p.cut*p.cut;
    ptr = utils::strtok_r(NULL," \t\n\r\f",&r_token); // 2nd token
    p.xi = (atoi(ptr)>0) ? true:false;
  }

  // set cutmax to max of all params
  cutmax = 0.0;
  for (int i = 0; i < npair; i++) {
    PairParameters & p = pairParameters[i];
    if (p.cut > cutmax) cutmax = p.cut;
  }
  cutmaxsq = cutmax*cutmax;

  // start reading tabular functions
  double * singletable = new double[nr];
  for (int i = 0; i < npair; i++) { // U
    PairParameters & p = pairParameters[i];
    if (comm->me == 0) {
      grab(fp,nr,singletable);
    }
    MPI_Bcast(singletable,nr,MPI_DOUBLE,0,world);
    p.U = new tabularFunction(nr,0.0,p.cut);
    (p.U)->set_values(nr,0.0,p.cut,singletable,epsilon);
  }
  for (int i = 0; i < npair; i++) { // V
    PairParameters & p = pairParameters[i];
    if (comm->me == 0) {
      grab(fp,nr,singletable);
    }
    MPI_Bcast(singletable,nr,MPI_DOUBLE,0,world);
    p.V = new tabularFunction(nr,0.0,p.cut);
    (p.V)->set_values(nr,0.0,p.cut,singletable,epsilon);
  }
  for (int i = 0; i < npair; i++) { // W
    PairParameters & p = pairParameters[i];
    if (comm->me == 0) {
      grab(fp,nr,singletable);
    }
    MPI_Bcast(singletable,nr,MPI_DOUBLE,0,world);
    p.W = new tabularFunction(nr,0.0,p.cut);
    (p.W)->set_values(nr,0.0,p.cut,singletable,epsilon);
  }
  for (int i = 0; i < npair; i++) { // P
    PairParameters & p = pairParameters[i];
    if (comm->me == 0) {
      grab(fp,nr,singletable);
    }
    MPI_Bcast(singletable,nr,MPI_DOUBLE,0,world);
    p.P = new tabularFunction(nr,-cutmax,cutmax);
    (p.P)->set_values(nr,-cutmax,cutmax,singletable,epsilon);
  }
  delete[] singletable;
  singletable = new double[ng];
  for (int i = 0; i < ntriple; i++) { // G
    TripletParameters & p = tripletParameters[i];
    if (comm->me == 0) {
      grab(fp,ng,singletable);
    }
    MPI_Bcast(singletable,ng,MPI_DOUBLE,0,world);
    p.G = new tabularFunction(ng,-1.0,1.0);
    (p.G)->set_values(ng,-1.0,1.0,singletable,epsilon);
  }
  delete[] singletable;
  singletable = new double[nx];
  for (int i = 0; i < npair; i++) { // F
    PairParameters & p = pairParameters[i];
    if (comm->me == 0) {
      grab(fp,nx,singletable);
    }
    MPI_Bcast(singletable,nx,MPI_DOUBLE,0,world);
    p.F = new tabularFunction(nx,0.0,maxX);
    (p.F)->set_values(nx,0.0,maxX,singletable,epsilon);
  }
  delete[] singletable;
  if (comm->me == 0) {
    fclose(fp);
  }

}

/* ---------------------------------------------------------------------- */

void PairPolymorphic::setup_params()
{
  int i,j,k,n;

  memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,"pair:elem2param");
  memory->destroy(elem3param);
  memory->create(elem3param,nelements,nelements,nelements,"pair:elem3param");

  // map atom pair to parameter index

  n = 0;
  for (i = 0; i < nelements; i++) {
    elem2param[match[i]][match[i]] = n;
    n++;
  }
  for (i = 0; i < nelements-1; i++) {
  for (j = i+1; j < nelements; j++) {
    elem2param[match[i]][match[j]] = n;
    elem2param[match[j]][match[i]] = n;
    n++;
  }
  }

  // map atom triplet to parameter index

  n = 0;
  for (i = 0; i < nelements; i++)
  for (j = 0; j < nelements; j++)
  for (k = 0; k < nelements; k++) {
    elem3param[match[i]][match[j]][match[k]] = n;
    n++;
  }

//   for debugging, call write_tables() to check the tabular functions
//   if (comm->me == 0) {
//     write_tables(51);
//     errorX->all(FLERR,"Test potential tables");
//   }
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

  (p->W)->value(rik,fc,1,dfc,1);
  (p->P)->value(rij-(p->xi)*rik,ex_delr,1,ex_delr_d,1);
  (trip->G)->value(cos_theta,gijk,1,gijk_d,1);

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

/* ----------------------------------------------------------------------
 *    grab n values from file fp and put them in list
 *       values can be several to a line
 *          only called by proc 0
 *          ------------------------------------------------------------------------- */

void PairPolymorphic::grab(FILE *fp, int n, double *list)
{
  char *ptr;
  char line[MAXLINE];
  char *r_token;

  int i = 0;
  while (i < n) {
    utils::sfgets(FLERR,line,MAXLINE,fp,NULL,error);
    r_token = line;
    ptr = utils::strtok_r(r_token," \t\n\r\f",&r_token);
    list[i++] = atof(ptr);
    while ((ptr = utils::strtok_r(NULL," \t\n\r\f",&r_token)))
      list[i++] = atof(ptr);
  }
}

/* ---------------------------------------------------------------------- */

void PairPolymorphic::write_tables(int npts)
{
  char tag[6] = "";
  if (comm->me != 0) sprintf(tag,"%d",comm->me);
  FILE* fp =  NULL;
  double  xmin,xmax,x,uf,vf,wf,pf,gf,ff,ufp,vfp,wfp,pfp,gfp,ffp;
  char line[MAXLINE];
  for (int i = 0; i < nelements; i++) {
  for (int j = 0; j < nelements; j++) {
    strcpy(line,elements[i]);
    strcat(line,elements[j]);
    strcat(line,"_UVW");
    strcat(line,tag);
    fp = fopen(line, "w");
    int iparam_ij = elem2param[i][j];
    PairParameters & pair = pairParameters[iparam_ij];
    xmin = (pair.U)->get_xmin();
    xmax = (pair.U)->get_xmax();
    double xl = xmax - xmin;
    xmin = xmin - 0.5*xl;
    xmax = xmax + 0.5*xl;
    for (int k = 0; k < npts; k++) {
      x = xmin + (xmax-xmin) * k / (npts-1);
      (pair.U)->value(x, uf, 1, ufp, 1);
      (pair.V)->value(x, vf, 1, vfp, 1);
      (pair.W)->value(x, wf, 1, wfp, 1);
      fprintf(fp,"%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f \n",x,uf,vf,wf,ufp,vfp,wfp);
    }
    fclose(fp);
  }
  }
  for (int i = 0; i < nelements; i++) {
  for (int j = 0; j < nelements; j++) {
    strcpy(line,elements[i]);
    strcat(line,elements[j]);
    strcat(line,"_P");
    strcat(line,tag);
    fp = fopen(line, "w");
    int iparam_ij = elem2param[i][j];
    PairParameters & pair = pairParameters[iparam_ij];
    xmin = (pair.P)->get_xmin();
    xmax = (pair.P)->get_xmax();
    double xl = xmax - xmin;
    xmin = xmin - 0.5*xl;
    xmax = xmax + 0.5*xl;
    for (int k = 0; k < npts; k++) {
      x = xmin + (xmax-xmin) * k / (npts-1);
      (pair.P)->value(x, pf, 1, pfp, 1);
      fprintf(fp,"%12.4f %12.4f %12.4f \n",x,pf,pfp);
    }
    fclose(fp);
  }
  }
  for (int i = 0; i < nelements; i++) {
  for (int j = 0; j < nelements; j++) {
  for (int k = 0; k < nelements; k++) {
    strcpy(line,elements[i]);
    strcat(line,elements[j]);
    strcat(line,elements[k]);
    strcat(line,"_G");
    strcat(line,tag);
    fp = fopen(line, "w");
    int iparam_ij = elem3param[i][j][k];
    TripletParameters & pair = tripletParameters[iparam_ij];
    xmin = (pair.G)->get_xmin();
    xmax = (pair.G)->get_xmax();
    for (int n = 0; n < npts; n++) {
      x = xmin + (xmax-xmin) * n / (npts-1);
      (pair.G)->value(x, gf, 1, gfp, 1);
      fprintf(fp,"%12.4f %12.4f %12.4f \n",x,gf,gfp);
    }
    fclose(fp);
  }
  }
  }
  for (int i = 0; i < nelements; i++) {
  for (int j = 0; j < nelements; j++) {
    strcpy(line,elements[i]);
    strcat(line,elements[j]);
    strcat(line,"_F");
    strcat(line,tag);
    fp = fopen(line, "w");
    int iparam_ij = elem2param[i][j];
    PairParameters & pair = pairParameters[iparam_ij];
    xmin = (pair.F)->get_xmin();
    xmax = (pair.F)->get_xmax();
    double xl = xmax - xmin;
    xmin = xmin - 0.5*xl;
    xmax = xmax + 0.5*xl;
    for (int k = 0; k < npts; k++) {
      x = xmin + (xmax-xmin) * k / (npts-1);
      (pair.F)->value(x, ff, 1, ffp, 1);
      fprintf(fp,"%12.4f %12.4f %12.4f \n",x,ff,ffp);
    }
    fclose(fp);
  }
  }

}

