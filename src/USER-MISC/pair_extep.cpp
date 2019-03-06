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
   Contributing author: Jan Los
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_extep.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "my_page.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 1024
#define DELTA 4
#define PGDELTA 1

/* ---------------------------------------------------------------------- */

PairExTeP::PairExTeP(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  ghostneigh = 1;

  nelements = 0;
  elements = NULL;
  nparams = maxparam = 0;
  params = NULL;
  elem2param = NULL;

  maxlocal = 0;
  SR_numneigh = NULL;
  SR_firstneigh = NULL;
  ipage = NULL;
  pgsize = oneatom = 0;
  map = NULL;

  Nt = NULL;
  Nd = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairExTeP::~PairExTeP()
{
  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  memory->destroy(params);
  memory->destroy(elem2param);

  memory->destroy(SR_numneigh);
  memory->sfree(SR_firstneigh);
  delete [] ipage;
  memory->destroy(Nt);
  memory->destroy(Nd);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cutghost);
    delete [] map;
  }
}

/* ----------------------------------------------------------------------
   create SR neighbor list from main neighbor list
   SR neighbor list stores neighbors of ghost atoms
------------------------------------------------------------------------- */

void PairExTeP::SR_neigh()
{
  int i,j,ii,jj,n,allnum,jnum,itype,jtype,iparam_ij;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *neighptr;

  double **x = atom->x;
  int *type = atom->type;

  if (atom->nmax > maxlocal) {  // ensure there is enough space
    maxlocal = atom->nmax;      // for atoms and ghosts allocated
    memory->destroy(SR_numneigh);
    memory->sfree(SR_firstneigh);
    memory->destroy(Nt);
    memory->destroy(Nd);
    memory->create(SR_numneigh,maxlocal,"ExTeP:numneigh");
    SR_firstneigh = (int **) memory->smalloc(maxlocal*sizeof(int *),
                           "ExTeP:firstneigh");
    memory->create(Nt,maxlocal,"ExTeP:Nt");
    memory->create(Nd,maxlocal,"ExTeP:Nd");
  }

  allnum = list->inum + list->gnum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // store all SR neighs of owned and ghost atoms
  // scan full neighbor list of I

  ipage->reset();

  for (ii = 0; ii < allnum; ii++) {
    i = ilist[ii];
    itype=map[type[i]];

    n = 0;
    neighptr = ipage->vget();

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    Nt[i] = 0.0;
    Nd[i] = 0.0;

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      jtype=map[type[j]];
      iparam_ij = elem2param[itype][jtype][jtype];

      if (rsq < params[iparam_ij].cutsq) {
        neighptr[n++] = j;
        double tmp_fc = ters_fc(sqrt(rsq),&params[iparam_ij]);
        Nt[i] += tmp_fc;
        if (itype!=jtype) {
          Nd[i] += tmp_fc;
        }
      }
    }
  //printf("SR_neigh : N[%d] = %f\n",i,N[i]);

    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }
}


/* ---------------------------------------------------------------------- */

void PairExTeP::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum;
  int itype,jtype,ktype,iparam_ij,iparam_ijk;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2,r2;
  double delr1[3],delr2[3],fi[3],fj[3],fk[3];
  double zeta_ij,prefactor;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = vflag_atom = 0;

  SR_neigh();

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

      jtype = map[type[j]];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      iparam_ij = elem2param[itype][jtype][jtype];
      if (rsq > params[iparam_ij].cutsq) continue;

      repulsive(&params[iparam_ij],rsq,fpair,eflag,evdwl);

      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,fpair,delx,dely,delz);
    }

    // three-body interactions      -(bij + Fcorrection) * fA
    // skip immediately if I-J is not within cutoff

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtag = tag[j];
      jtype = map[type[j]];
      iparam_ij = elem2param[itype][jtype][jtype];

      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 > params[iparam_ij].cutsq) continue;

      // accumulate bondorder zeta for each i-j interaction via loop over k

      zeta_ij = 0.0;

      /* F_IJ (1) */
      // compute correction to energy and forces
      // dE/dr = -Fij(Zi,Zj) dV/dr
      //         - dFij/dZi dZi/dr V
      //         (conjugate term is computed when j is a central atom)

      double FXY, dFXY_dNdij, dFXY_dNdji, fa, fa_d, deng, fpair;
      double Ntij = Nt[i];
      double Ndij = Nd[i];
      double Ntji = Nt[j];
      double Ndji = Nd[j];
      double r = sqrt(rsq1);
      double fc_ij = ters_fc(r,&params[iparam_ij]);

      Ntij -= fc_ij;
      Ntji -= fc_ij;
      if (jtype!=itype) {
        Ndij -= fc_ij;
        Ndji -= fc_ij;
      }
      if (Ntij<0) { Ntij=0.; }
      if (Ndij<0) { Ndij=0.; }
      if (Ntji<0) { Ntji=0.; }
      if (Ndji<0) { Ndji=0.; }
      FXY = F_corr(itype, jtype, Ndij, Ndji, &dFXY_dNdij, &dFXY_dNdji);

      // envelop functions
      double fenv, dfenv_ij;
      fenv = envelop_function(Ntij, Ntji, &dfenv_ij);
      //
      double Fc = fenv * FXY;
      double dFc_dNtij = dfenv_ij * FXY;
      double dFc_dNdij = fenv * dFXY_dNdij;

      fa = ters_fa(r,&params[iparam_ij]);
      fa_d = ters_fa_d(r,&params[iparam_ij]);
      deng = 0.5 * fa * Fc;
      fpair = 0.5 * fa_d * Fc / r;

      f[i][0] += delr1[0]*fpair;
      f[i][1] += delr1[1]*fpair;
      f[i][2] += delr1[2]*fpair;
      f[j][0] -= delr1[0]*fpair;
      f[j][1] -= delr1[1]*fpair;
      f[j][2] -= delr1[2]*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           deng,0.0,-fpair,-delr1[0],-delr1[1],-delr1[2]);
      /* END F_IJ (1) */

      for (kk = 0; kk < jnum; kk++) {
        if (jj == kk) continue;
        k = jlist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        iparam_ijk = elem2param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 > params[iparam_ijk].cutsq) continue;

        r2 = sqrt(rsq2);

        zeta_ij += zeta(&params[iparam_ijk],r,r2,delr1,delr2);

        /* F_IJ (2) */
        // compute force components due to spline derivatives
        // uses only the part with FXY_x (FXY_y is done when i and j are inversed)
        int iparam_ik = elem2param[itype][ktype][0];
        double fc_ik_d = ters_fc_d(r2,&params[iparam_ik]);
        double fc_prefac_ik_0 = 1.0 * fc_ik_d * fa / r2;
        double fc_prefac_ik = dFc_dNtij * fc_prefac_ik_0;
        f[i][0] += fc_prefac_ik * delr2[0];
        f[i][1] += fc_prefac_ik * delr2[1];
        f[i][2] += fc_prefac_ik * delr2[2];
        f[k][0] -= fc_prefac_ik * delr2[0];
        f[k][1] -= fc_prefac_ik * delr2[1];
        f[k][2] -= fc_prefac_ik * delr2[2];
        if ( itype != ktype ) {
          fc_prefac_ik = dFc_dNdij * fc_prefac_ik_0;
          f[i][0] += fc_prefac_ik * delr2[0];
          f[i][1] += fc_prefac_ik * delr2[1];
          f[i][2] += fc_prefac_ik * delr2[2];
          f[k][0] -= fc_prefac_ik * delr2[0];
          f[k][1] -= fc_prefac_ik * delr2[1];
          f[k][2] -= fc_prefac_ik * delr2[2];
        }
        /* END F_IJ (2) */

      }

      // pairwise force due to zeta

      force_zeta(&params[iparam_ij],r,zeta_ij,fpair,prefactor,eflag,evdwl);

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
        ktype = map[type[k]];
        iparam_ijk = elem2param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 > params[iparam_ijk].cutsq) continue;

        attractive(&params[iparam_ijk],prefactor,
                   rsq1,rsq2,delr1,delr2,fi,fj,fk);


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

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairExTeP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cutghost,n+1,n+1,"pair:cutghost");

  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairExTeP::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairExTeP::coeff(int narg, char **arg)
{
  int i,j,n;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

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
  spline_init();
  setup();

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

void PairExTeP::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style ExTeP requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style ExTeP requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  // including neighbors of ghosts
  neighbor->requests[irequest]->ghost = 1;

  // create pages if first time or if neighbor pgsize/oneatom has changed

  int create = 0;
  if (ipage == NULL) create = 1;
  if (pgsize != neighbor->pgsize) create = 1;
  if (oneatom != neighbor->oneatom) create = 1;

  if (create) {
    delete [] ipage;
    pgsize = neighbor->pgsize;
    oneatom = neighbor->oneatom;

    int nmypage= comm->nthreads;
    ipage = new MyPage<int>[nmypage];
    for (int i = 0; i < nmypage; i++)
      ipage[i].init(oneatom,pgsize,PGDELTA);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairExTeP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  cutghost[i][j] = cutmax ;
  cutghost[j][i] = cutghost[i][j];

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairExTeP::read_file(char *file)
{
  int params_per_line = 17;
  char **words = new char*[params_per_line+1];

  memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential(file);
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open ExTeP potential file %s",file);
      error->one(FLERR,str);
    }
  }

  // read each line out of file, skipping blank lines or leading '#'
  // store line of params if all 3 element tags are in element list

  int n,nwords,ielement,jelement,kelement;
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
      error->all(FLERR,"Insufficient spline parameters in potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    // ielement,jelement,kelement = 1st args
    // if all 3 args are in element list, then parse this line
    // else skip to next line

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;
    for (kelement = 0; kelement < nelements; kelement++)
      if (strcmp(words[2],elements[kelement]) == 0) break;
    if (kelement == nelements) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                          "pair:params");
    }

    params[nparams].ielement = ielement;
    params[nparams].jelement = jelement;
    params[nparams].kelement = kelement;
    params[nparams].powerm = atof(words[3]);
    params[nparams].gamma = atof(words[4]);
    params[nparams].lam3 = atof(words[5]);
    params[nparams].c = atof(words[6]);
    params[nparams].d = atof(words[7]);
    params[nparams].h = atof(words[8]);
    params[nparams].powern = atof(words[9]);
    params[nparams].beta = atof(words[10]);
    params[nparams].lam2 = atof(words[11]);
    params[nparams].bigb = atof(words[12]);
    params[nparams].bigr = atof(words[13]);
    params[nparams].bigd = atof(words[14]);
    params[nparams].lam1 = atof(words[15]);
    params[nparams].biga = atof(words[16]);

    // currently only allow m exponent of 1 or 3

    params[nparams].powermint = int(params[nparams].powerm);

    if (params[nparams].c < 0.0 || params[nparams].d < 0.0 ||
        params[nparams].powern < 0.0 || params[nparams].beta < 0.0 ||
        params[nparams].lam2 < 0.0 || params[nparams].bigb < 0.0 ||
        params[nparams].bigr < 0.0 ||params[nparams].bigd < 0.0 ||
        params[nparams].bigd > params[nparams].bigr ||
        params[nparams].lam1 < 0.0 || params[nparams].biga < 0.0 ||
        params[nparams].powerm - params[nparams].powermint != 0.0 ||
        (params[nparams].powermint != 3 && params[nparams].powermint != 1) ||
        params[nparams].gamma < 0.0)
      error->all(FLERR,"Illegal ExTeP parameter");

    nparams++;
    if (nparams >= pow(atom->ntypes,3)) break;
  }

  // deallocate words array
  delete [] words;

  /* F_IJ (3) */
  // read the spline coefficients
  params_per_line = 8;
  // reallocate with new size
  words = new char*[params_per_line+1];

  // intialize F_corr_data to all zeros
  for (int iel=0;iel<atom->ntypes;iel++)
    for (int jel=0;jel<atom->ntypes;jel++)
      for (int in=0;in<4;in++)
        for (int jn=0;jn<4;jn++)
          for (int ivar=0;ivar<3;ivar++)
            F_corr_data[iel][jel][in][jn][ivar]=0;

  // loop until EOF
  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      //fputs(line,stdout);
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

    if (nwords != params_per_line)
      error->all(FLERR,"Incorrect format in ExTeP potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    // ielement,jelement = 1st args
    // if all 3 args are in element list, then parse this line
    // else skip to next line
    // these lines set ielement and jelement to the
    // integers matching the strings from the input

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;

    int Ni  = atoi(words[2]);
    int Nj  = atoi(words[3]);
    double spline_val = atof(words[4]);
    double spline_derx = atof(words[5]);
    double spline_dery = atof(words[6]);

    // Set value for all pairs of ielement,jelement  (any kelement)
    for (int iparam = 0; iparam < nparams; iparam++) {
      if ( ielement == params[iparam].ielement
           && jelement == params[iparam].jelement) {
        F_corr_data[ielement][jelement][Ni][Nj][0] = spline_val;
        F_corr_data[ielement][jelement][Ni][Nj][1] = spline_derx;
        F_corr_data[ielement][jelement][Ni][Nj][2] = spline_dery;

        F_corr_data[jelement][ielement][Nj][Ni][0] = spline_val;
        F_corr_data[jelement][ielement][Nj][Ni][1] = spline_dery;
        F_corr_data[jelement][ielement][Nj][Ni][2] = spline_derx;
      }
    }
  }

  delete [] words;
  /* END F_IJ (3) */

}

/* ---------------------------------------------------------------------- */

void PairExTeP::setup()
{
  int i,j,k,m,n;

  // set elem2param for all element triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,nelements,"pair:elem2param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement &&
              k == params[m].kelement) {
            if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
            n = m;
          }
        }
        if (n < 0) error->all(FLERR,"Potential file is missing an entry");
        elem2param[i][j][k] = n;
      }

  // compute parameter values derived from inputs

  for (m = 0; m < nparams; m++) {
    params[m].cut = params[m].bigr + params[m].bigd;
    params[m].cutsq = params[m].cut*params[m].cut;

    params[m].c1 = pow(2.0*params[m].powern*1.0e-16,-1.0/params[m].powern);
    params[m].c2 = pow(2.0*params[m].powern*1.0e-8,-1.0/params[m].powern);
    params[m].c3 = 1.0/params[m].c2;
    params[m].c4 = 1.0/params[m].c1;
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++)
    if (params[m].cut > cutmax) cutmax = params[m].cut;
}

/* ---------------------------------------------------------------------- */

void PairExTeP::repulsive(Param *param, double rsq, double &fforce,
                            int eflag, double &eng)
{
  double r,tmp_fc,tmp_fc_d,tmp_exp;

  r = sqrt(rsq);
  tmp_fc = ters_fc(r,param);
  tmp_fc_d = ters_fc_d(r,param);
  tmp_exp = exp(-param->lam1 * r);
  fforce = -param->biga * tmp_exp * (tmp_fc_d - tmp_fc*param->lam1) / r;
  if (eflag) eng = tmp_fc * param->biga * tmp_exp;
}

/* ---------------------------------------------------------------------- */

double PairExTeP::zeta(Param *param, double rij, double rik,
                         double *delrij, double *delrik)
{
  double costheta,arg,ex_delr;

  costheta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] +
              delrij[2]*delrik[2]) / (rij*rik);

  if (param->powermint == 3) arg = pow(param->lam3 * (rij-rik),3.0);
  else arg = param->lam3 * (rij-rik);

  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(arg);

  return ters_fc(rik,param) * ters_gijk(costheta,param) * ex_delr;
}

/* ---------------------------------------------------------------------- */

void PairExTeP::force_zeta(Param *param, double r, double zeta_ij,
                             double &fforce, double &prefactor,
                             int eflag, double &eng)
{
  double fa,fa_d,bij;

  fa = ters_fa(r,param);
  fa_d = ters_fa_d(r,param);
  bij = ters_bij(zeta_ij,param);
  fforce = 0.5*bij*fa_d / r;
  prefactor = -0.5*fa * ( ters_bij_d(zeta_ij,param) );
  if (eflag) eng = 0.5*bij*fa;
}

/* ----------------------------------------------------------------------
   attractive term
   use param_ij cutoff for rij test
   use param_ijk cutoff for rik test
------------------------------------------------------------------------- */

void PairExTeP::attractive(Param *param, double prefactor,
                             double rsqij, double rsqik,
                             double *delrij, double *delrik,
                             double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;

  rij = sqrt(rsqij);
  rijinv = 1.0/rij;
  vec3_scale(rijinv,delrij,rij_hat);

  rik = sqrt(rsqik);
  rikinv = 1.0/rik;
  vec3_scale(rikinv,delrik,rik_hat);

  ters_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk,param);
}

/* ---------------------------------------------------------------------- */

double PairExTeP::ters_fc(double r, Param *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;

  if (r < ters_R-ters_D) return 1.0;
  if (r > ters_R+ters_D) return 0.0;
  return 0.5*(1.0 - sin(MY_PI2*(r - ters_R)/ters_D));
}

/* ---------------------------------------------------------------------- */

double PairExTeP::ters_fc_d(double r, Param *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;

  if (r < ters_R-ters_D) return 0.0;
  if (r > ters_R+ters_D) return 0.0;
  return -(MY_PI4/ters_D) * cos(MY_PI2*(r - ters_R)/ters_D);
}

/* ---------------------------------------------------------------------- */

double PairExTeP::ters_fa(double r, Param *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return -param->bigb * exp(-param->lam2 * r) * ters_fc(r,param);
}

/* ---------------------------------------------------------------------- */

double PairExTeP::ters_fa_d(double r, Param *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return param->bigb * exp(-param->lam2 * r) *
    (param->lam2 * ters_fc(r,param) - ters_fc_d(r,param));
}

/* ---------------------------------------------------------------------- */

double PairExTeP::ters_bij(double zeta, Param *param)
{
  double tmp = param->beta * zeta;
  if (tmp > param->c1) return 1.0/sqrt(tmp);
  if (tmp > param->c2)
    return (1.0 - pow(tmp,-param->powern) / (2.0*param->powern))/sqrt(tmp);
  if (tmp < param->c4) return 1.0;
  if (tmp < param->c3)
    return 1.0 - pow(tmp,param->powern)/(2.0*param->powern);
  return pow(1.0 + pow(tmp,param->powern), -1.0/(2.0*param->powern));
}

/* ---------------------------------------------------------------------- */

double PairExTeP::ters_bij_d(double zeta, Param *param)
{
  double tmp = param->beta * zeta;
  if (tmp > param->c1) return param->beta * -0.5*pow(tmp,-1.5);
  if (tmp > param->c2)
    return param->beta * (-0.5*pow(tmp,-1.5) *
                          (1.0 - 0.5*(1.0 +  1.0/(2.0*param->powern)) *
                           pow(tmp,-param->powern)));
  if (tmp < param->c4) return 0.0;
  if (tmp < param->c3)
    return -0.5*param->beta * pow(tmp,param->powern-1.0);

  double tmp_n = pow(tmp,param->powern);
  return -0.5 * pow(1.0+tmp_n, -1.0-(1.0/(2.0*param->powern)))*tmp_n / zeta;
}

/* ---------------------------------------------------------------------- */

void PairExTeP::ters_zetaterm_d(double prefactor,
                                double *rij_hat, double rij,
                                double *rik_hat, double rik,
                                double *dri, double *drj, double *drk,
                                Param *param)
{
  double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc(rik,param);
  dfc = ters_fc_d(rik,param);
  if (param->powermint == 3) tmp = pow(param->lam3 * (rij-rik),3.0);
  else tmp = param->lam3 * (rij-rik);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (param->powermint == 3)
    ex_delr_d = 3.0*pow(param->lam3,3.0) * pow(rij-rik,2.0)*ex_delr;
  else ex_delr_d = param->lam3 * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk(cos_theta,param);
  gijk_d = ters_gijk_d(cos_theta,param);
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

void PairExTeP::costheta_d(double *rij_hat, double rij,
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


/* ---------------------------------------------------------------------- */

/* F_IJ (4) */
// initialize spline for F_corr (based on PairLCBOP::F_conj)

void PairExTeP::spline_init() {
  for ( int iel=0; iel<atom->ntypes; iel++) {
    for ( int jel=0; jel<atom->ntypes; jel++) {
      for ( int N_ij=0; N_ij<4; N_ij++ ) {
        for ( int N_ji=0; N_ji<4; N_ji++ ) {
          TF_corr_param &f = F_corr_param[iel][jel][N_ij][N_ji];

          // corner points for each spline function
          f.f_00 = F_corr_data[iel][jel][N_ij  ][N_ji  ][0];
          f.f_01 = F_corr_data[iel][jel][N_ij  ][N_ji+1][0];
          f.f_10 = F_corr_data[iel][jel][N_ij+1][N_ji  ][0];
          f.f_11 = F_corr_data[iel][jel][N_ij+1][N_ji+1][0];

          // compute f tilde according to [Los & Fasolino, PRB 68, 024107 2003]
          f.f_x_00 =   F_corr_data[iel][jel][N_ij  ][N_ji  ][1] - f.f_10 + f.f_00;
          f.f_x_01 =   F_corr_data[iel][jel][N_ij  ][N_ji+1][1] - f.f_11 + f.f_01;
          f.f_x_10 = -(F_corr_data[iel][jel][N_ij+1][N_ji  ][1] - f.f_10 + f.f_00);
          f.f_x_11 = -(F_corr_data[iel][jel][N_ij+1][N_ji+1][1] - f.f_11 + f.f_01);

          f.f_y_00 =   F_corr_data[iel][jel][N_ij  ][N_ji  ][2] - f.f_01 + f.f_00;
          f.f_y_01 = -(F_corr_data[iel][jel][N_ij  ][N_ji+1][2] - f.f_01 + f.f_00);
          f.f_y_10 =   F_corr_data[iel][jel][N_ij+1][N_ji  ][2] - f.f_11 + f.f_10;
          f.f_y_11 = -(F_corr_data[iel][jel][N_ij+1][N_ji+1][2] - f.f_11 + f.f_10);
        }
      }
    }
  }
}

double PairExTeP::envelop_function(double x, double y, double *func_der) {
  double fx,fy,fxy,dfx,dfxydx;
  double del, delsq;

  fxy = 1.0;
  dfxydx = 0.0;

  if (x <= 3.0) {
    fx = 1.0;
    dfx = 0.0;
    if (x < 1.0 && y < 1.0) {
      double gx=(1.0-x);
      double gy=(1.0-y);
      double gxsq=gx*gx;
      double gysq=gy*gy;
      fxy = 1.0 - gxsq*gysq;
      dfxydx = 2.0*gx*gysq;
    }
  } else if (x < 4.0) {
    del = 4.0-x;
    delsq = del*del;
    fx = (3.0-2.0*del)*delsq;
    dfx = - 6.0*del*(1.0-del);
  } else {
    fx = 0.0;
    dfx = 0.0;
  }
  if (y <= 3.0) {
    fy = 1.0;
  } else if (y < 4.0) {
    del = 4.0-y;
    delsq = del*del;
    fy = (3.0-2.0*del)*delsq;
  } else {
    fy = 0.0;
  }

  double func_val = fxy*fx*fy;
  *func_der = dfxydx*fx*fy+fxy*dfx*fy;

  return func_val;
}

double PairExTeP::F_corr(int iel, int jel, double Ndij, double Ndji, double *dFN_x, double *dFN_y ) {

  // compute F_XY

  int Ndij_int         = static_cast<int>( floor( Ndij ) );
  int Ndji_int         = static_cast<int>( floor( Ndji ) );
  double x                = Ndij - Ndij_int;
  double y                = Ndji - Ndji_int;
  TF_corr_param &f  = F_corr_param[iel][jel][Ndij_int][Ndji_int];
  double F   = 0;
  double dF_dx = 0, dF_dy = 0;
  double l, r;
  if (Ndij_int < 4 && Ndji_int < 4) {
    l = (1-y)* (1-x);
    r = ( f.f_00 + x*x* f.f_x_10 + y*y* f.f_y_01 );
    F += l*r;
    dF_dx += -(1-y)*r +l*2*x* f.f_x_10;
    dF_dy += -(1-x)*r +l*2*y* f.f_y_01;
    l = (1-y)*x;
    r = ( f.f_10 + (1-x)*(1-x)*f.f_x_00 + y*  y* f.f_y_11 );
    F += l*r;
    dF_dx += (1-y)*r -l*2*(1-x)*f.f_x_00;
    dF_dy += -x*r +l*2*y* f.f_y_11;
    l = y*  (1-x);
    r = ( f.f_01 + x*x* f.f_x_11 + (1-y)*(1-y)*f.f_y_00 );
    F += l*r;
    dF_dx += -y*r +l*2*x* f.f_x_11;
    dF_dy += (1-x)*r -l*2*(1-y)*f.f_y_00;
    l = y*  x;
    r = ( f.f_11 + (1-x)*(1-x)*f.f_x_01 + (1-y)*(1-y)*f.f_y_10 );
    F += l*r;
    dF_dx += y*r -l*2*(1-x)*f.f_x_01;
    dF_dy += x*r -l*2*(1-y)*f.f_y_10;
  }
  double result = F;
  *dFN_x = dF_dx;
  *dFN_y = dF_dy;

  return result;
}
/* F_IJ (4) */
