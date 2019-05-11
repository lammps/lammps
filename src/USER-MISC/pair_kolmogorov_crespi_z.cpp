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
   Contributing author: Jaap Kroes (Radboud Universiteit)
   e-mail: jaapkroes at gmail dot com
   based on previous versions by Merel van Wijk and by Marco Raguzzoni

   This is a simplified version of the potential described in
   [Kolmogorov & Crespi, Phys. Rev. B 71, 235415 (2005)]
   The simplification is that all normals are taken along the z-direction
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_kolmogorov_crespi_z.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

PairKolmogorovCrespiZ::PairKolmogorovCrespiZ(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;

  // initialize element to parameter maps
  nelements = 0;
  elements = NULL;
  nparams = maxparam = 0;
  params = NULL;
  elem2param = NULL;
  map = NULL;

  // always compute energy offset
  offset_flag = 1;
}

/* ---------------------------------------------------------------------- */

PairKolmogorovCrespiZ::~PairKolmogorovCrespiZ()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(offset);
  }

  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  memory->destroy(params);
  memory->destroy(elem2param);
  if (allocated) delete [] map;
}

/* ---------------------------------------------------------------------- */

void PairKolmogorovCrespiZ::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair, fpair1;
  double rsq,r,rhosq,exp1,exp2,r6,r8;
  double frho,sumC,sumC2,sumCff,fsum,rdsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

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
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      // rho^2 = r^2 - (n,r) = r^2 - z^2
      rhosq = delx*delx + dely*dely;
      rsq = rhosq + delz*delz;

      if (rsq < cutsq[itype][jtype]) {

        int iparam_ij = elem2param[map[itype]][map[jtype]];
        Param& p = params[iparam_ij];

        r = sqrt(rsq);
        r6 = rsq*rsq*rsq;
        r8 = r6*rsq;
        rdsq = rhosq*p.delta2inv; // (rho/delta)^2

        // store exponents
        exp1 = exp(-p.lambda*(r-p.z0));
        exp2 = exp(-rdsq);

        // note that f(rho_ij) equals f(rho_ji) as normals are all along z
        sumC = p.C0+p.C2*rdsq+p.C4*rdsq*rdsq;
        sumC2 = (2*p.C2+4*p.C4*rdsq)*p.delta2inv;
        frho = exp2*sumC;
        sumCff = p.C + 2*frho;

        // derivatives
        fpair = -6.0*p.A*p.z06/r8+p.lambda*exp1/r*sumCff;
        fpair1 = exp1*exp2*(4.0*p.delta2inv*sumC-2.0*sumC2);
        fsum = fpair + fpair1;

        f[i][0] += delx*fsum;
        f[i][1] += dely*fsum;
        // fi_z does not contain contributions from df/dr
        // because rho_ij does not depend on z_i or z_j
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fsum;
          f[j][1] -= dely*fsum;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = -p.A*p.z06/r6+ exp1*sumCff - offset[itype][jtype];
        }

        if (evflag){
          ev_tally_xyz(i,j,nlocal,newton_pair,evdwl,0,
                       fsum,fsum,fpair,delx,dely,delz);
        }
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairKolmogorovCrespiZ::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(offset,n+1,n+1,"pair:offset");
  map = new int[atom->ntypes+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairKolmogorovCrespiZ::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");
  if (strcmp(force->pair_style,"hybrid/overlay")!=0)
    error->all(FLERR,"ERROR: requires hybrid/overlay pair_style");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairKolmogorovCrespiZ::coeff(int narg, char **arg)
{
  int i,j,n;

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

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

  read_file(arg[2]);

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
        cut[i][j] = cut_global;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairKolmogorovCrespiZ::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  if (!offset_flag)
    error->all(FLERR,"Must use 'pair_modify shift yes' with this pair style");

  if (offset_flag && (cut[i][j] > 0.0)) {
    int iparam_ij = elem2param[map[i]][map[j]];
    Param& p = params[iparam_ij];
    offset[i][j] = -p.A*pow(p.z0/cut[i][j],6);
  } else offset[i][j] = 0.0;
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   read Kolmogorov-Crespi potential file
------------------------------------------------------------------------- */

void PairKolmogorovCrespiZ::read_file(char *filename)
{
  int params_per_line = 11;
  char **words = new char*[params_per_line+1];
  memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential(filename);
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open KC potential file %s",filename);
      error->one(FLERR,str);
    }
  }

  // read each line out of file, skipping blank lines or leading '#'
  // store line of params if all 3 element tags are in element list

  int i,j,n,m,nwords,ielement,jelement;
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
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR,"Insufficient format in KC potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    // ielement,jelement = 1st args
    // if these 2 args are in element list, then parse this line
    // else skip to next line (continue)

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                          "pair:params");
    }

    params[nparams].ielement = ielement;
    params[nparams].jelement = jelement;
    params[nparams].z0       = atof(words[2]);
    params[nparams].C0       = atof(words[3]);
    params[nparams].C2       = atof(words[4]);
    params[nparams].C4       = atof(words[5]);
    params[nparams].C        = atof(words[6]);
    params[nparams].delta    = atof(words[7]);
    params[nparams].lambda   = atof(words[8]);
    params[nparams].A        = atof(words[9]);
    // S provides a convenient scaling of all energies
    params[nparams].S        = atof(words[10]);

    // energies in meV further scaled by S
    double meV = 1.0e-3*params[nparams].S;
    params[nparams].C *= meV;
    params[nparams].A *= meV;
    params[nparams].C0 *= meV;
    params[nparams].C2 *= meV;
    params[nparams].C4 *= meV;

    // precompute some quantities
    params[nparams].delta2inv = pow(params[nparams].delta,-2);
    params[nparams].z06 = pow(params[nparams].z0,6);

    nparams++;
    if(nparams >= pow(atom->ntypes,3)) break;
  }
  memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,"pair:elem2param");
  for (i = 0; i < nelements; i++) {
    for (j = 0; j < nelements; j++) {
      n = -1;
      for (m = 0; m < nparams; m++) {
        if (i == params[m].ielement && j == params[m].jelement) {
          if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
          n = m;
        }
      }
      if (n < 0) error->all(FLERR,"Potential file is missing an entry");
      elem2param[i][j] = n;
    }
  }
  delete [] words;
}

/* ---------------------------------------------------------------------- */
