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
   Contributing author: Tzu-Ray Shan (U Florida, present: tnshan@sandia.gov)
   LAMMPS implementation of the Charge-optimized many-body (COMB) potential
   based on the HELL MD program (Prof Simon Phillpot, UF, sphil@mse.ufl.edu)
   and Aidan Thompson's Tersoff code in LAMMPS
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_comb.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "group.h"
#include "update.h"
#include "my_page.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 1024
#define DELTA 4
#define PGDELTA 1
#define MAXNEIGH 24

/* ---------------------------------------------------------------------- */

PairComb::PairComb(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  nmax = 0;
  NCo = NULL;
  bbij = NULL;

  nelements = 0;
  elements = NULL;
  nparams = 0;
  maxparam = 0;
  params = NULL;
  elem2param = NULL;

  intype = NULL;
  fafb = NULL;
  dfafb = NULL;
  ddfafb = NULL;
  phin = NULL;
  dphin = NULL;
  erpaw = NULL;
  
  sht_num = NULL;
  sht_first = NULL;

  ipage = NULL;
  pgsize = oneatom = 0;

  // set comm size needed by this Pair

  comm_forward = 1;
  comm_reverse = 1;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairComb::~PairComb()
{
  memory->destroy(NCo);

  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];

  delete [] elements;
  memory->sfree(params);
  memory->destroy(elem2param);

  memory->destroy(intype);
  memory->destroy(fafb);
  memory->destroy(dfafb);
  memory->destroy(ddfafb);
  memory->destroy(phin);
  memory->destroy(dphin);
  memory->destroy(erpaw);
  memory->destroy(bbij);
  memory->destroy(sht_num);
  memory->sfree(sht_first);

  delete [] ipage;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
    delete [] esm;
  }
}

/* ---------------------------------------------------------------------- */

void PairComb::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum,iparam_i;
  int itype,jtype,ktype,iparam_ij,iparam_ijk;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fi[3],fj[3],fk[3];
  double zeta_ij,prefactor;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int mr1,mr2,mr3;
  int rsc,inty;
  double elp_ij,filp[3],fjlp[3],fklp[3];
  double iq,jq;
  double yaself;
  double potal,fac11,fac11e;
  double vionij,fvionij,sr1,sr2,sr3,Eov,Fov;
  int sht_jnum, *sht_jlist, nj;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = vflag_atom = 0;

  // Build short range neighbor list

  Short_neigh();

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  yaself = vionij = fvionij = Eov = Fov = 0.0;

  // self energy correction term: potal

  potal_calc(potal,fac11,fac11e);

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    iq = q[i];
    NCo[i] = 0;
    nj = 0;
    iparam_i = elem2param[itype][itype][itype];

    // self energy, only on i atom

    yaself = self(&params[iparam_i],iq,potal);

    if (evflag) ev_tally(i,i,nlocal,0,yaself,0.0,0.0,0.0,0.0,0.0);

    // two-body interactions (long and short repulsive)

    jlist = firstneigh[i];
    jnum = numneigh[i];
    sht_jlist = sht_first[i];
    sht_jnum = sht_num[i];

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

      // Qj calculates 2-body Coulombic

      jtype = map[type[j]];
      jq = q[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      iparam_ij = elem2param[itype][jtype][jtype];

      // long range q-dependent

      if (rsq > params[iparam_ij].lcutsq) continue;

      inty = intype[itype][jtype];

      // polynomial three-point interpolation

      tri_point(rsq, mr1, mr2, mr3, sr1, sr2, sr3, itype);

      // 1/r energy and forces

      direct(inty,mr1,mr2,mr3,rsq,sr1,sr2,sr3,iq,jq,
             potal,fac11,fac11e,vionij,fvionij);

      // field correction to self energy

      field(&params[iparam_ij],rsq,iq,jq,vionij,fvionij);

      // polarization field
      // sums up long range forces

      f[i][0] += delx*fvionij;
      f[i][1] += dely*fvionij;
      f[i][2] += delz*fvionij;
      f[j][0] -= delx*fvionij;
      f[j][1] -= dely*fvionij;
      f[j][2] -= delz*fvionij;

      if (evflag)
        ev_tally(i,j,nlocal,newton_pair,0.0,vionij,fvionij,delx,dely,delz);

      // short range q-independent

      if (rsq > params[iparam_ij].cutsq) continue;

      repulsive(&params[iparam_ij],rsq,fpair,eflag,evdwl,iq,jq);

      // repulsion is pure two-body, sums up pair repulsive forces

      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag)
        ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
    }

    // accumulate coordination number information

    if (cor_flag) {
      for (jj = 0; jj < sht_jnum; jj++) {
        j = sht_jlist[jj];
        jtype = map[type[j]];
        iparam_ij = elem2param[itype][jtype][jtype];

        if(params[iparam_ij].hfocor > 0.0 ) {
          delr1[0] = x[j][0] - xtmp;
          delr1[1] = x[j][1] - ytmp;
          delr1[2] = x[j][2] - ztmp;
          rsq1 = vec3_dot(delr1,delr1);

          if (rsq1 > params[iparam_ij].cutsq) continue;
          NCo[i] += 1;
        }
      }
    }

    // three-body interactions
    // half i-j loop

    for (jj = 0; jj < sht_jnum; jj++) {
      j = sht_jlist[jj];

      jtype = map[type[j]];
      iparam_ij = elem2param[itype][jtype][jtype];

      // this Qj for q-dependent BSi

      jq = q[j];

      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = vec3_dot(delr1,delr1);

      if (rsq1 > params[iparam_ij].cutsq) continue;
      nj ++;

      // accumulate bondorder zeta for each i-j interaction via loop over k

      zeta_ij = 0.0;
      cuo_flag1 = 0; cuo_flag2 = 0;

      for (kk = 0; kk < sht_jnum; kk++) {
        k = sht_jlist[kk];
        if (j == k) continue;
        ktype = map[type[k]];
        iparam_ijk = elem2param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = vec3_dot(delr2,delr2);

        if (rsq2 > params[iparam_ijk].cutsq) continue;

        zeta_ij += zeta(&params[iparam_ijk],rsq1,rsq2,delr1,delr2);

        if (params[iparam_ijk].hfocor == -2.0) cuo_flag1 = 1;
        if (params[iparam_ijk].hfocor == -1.0) cuo_flag2 = 1;
      }

      if (cuo_flag1 && cuo_flag2) cuo_flag = 1;
      else cuo_flag = 0;

      force_zeta(&params[iparam_ij],eflag,i,nj,rsq1,zeta_ij,
                 iq,jq,fpair,prefactor,evdwl);

      // over-coordination correction for HfO2

      if (cor_flag && NCo[i] != 0)
        Over_cor(&params[iparam_ij],rsq1,NCo[i],Eov, Fov);
      evdwl +=  Eov;
      fpair +=  Fov;

      f[i][0] += delr1[0]*fpair;
      f[i][1] += delr1[1]*fpair;
      f[i][2] += delr1[2]*fpair;
      f[j][0] -= delr1[0]*fpair;
      f[j][1] -= delr1[1]*fpair;
      f[j][2] -= delr1[2]*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,-fpair,-delr1[0],-delr1[1],-delr1[2]);

      // attractive term via loop over k (3-body forces)

      for (kk = 0; kk < sht_jnum; kk++) {
        k = sht_jlist[kk];
        if (j == k) continue;
        ktype = map[type[k]];
        iparam_ijk = elem2param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = vec3_dot(delr2,delr2);
        if (rsq2 > params[iparam_ijk].cutsq) continue;

        for (rsc = 0; rsc < 3; rsc++)
          fi[rsc] = fj[rsc] = fk[rsc] = 0.0;

        attractive(&params[iparam_ijk],prefactor,
                   rsq1,rsq2,delr1,delr2,fi,fj,fk);

        // 3-body LP and BB correction and forces

        elp_ij = elp(&params[iparam_ijk],rsq1,rsq2,delr1,delr2);
        flp(&params[iparam_ijk],rsq1,rsq2,delr1,delr2,filp,fjlp,fklp);

        for (rsc = 0; rsc < 3; rsc++) {
          fi[rsc] += filp[rsc];
          fj[rsc] += fjlp[rsc];
          fk[rsc] += fklp[rsc];
        }

        for (rsc = 0; rsc < 3; rsc++) {
          f[i][rsc] += fi[rsc];
          f[j][rsc] += fj[rsc];
          f[k][rsc] += fk[rsc];
        }

        if (evflag)
          ev_tally(i,j,nlocal,newton_pair,elp_ij,0.0,0.0,0.0,0.0,0.0);
        if (vflag_atom) v_tally3(i,j,k,fj,fk,delr1,delr2);

      }
    }

    if (cuo_flag) params[iparam_i].cutsq *= 0.65;
  }

  cuo_flag = 0;

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairComb::allocate()
{
 allocated = 1;
 int n = atom->ntypes;

 memory->create(setflag,n+1,n+1,"pair:setflag");
 memory->create(cutsq,n+1,n+1,"pair:cutsq");

 map = new int[n+1];
 esm = new double[n];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairComb::settings(int narg, char **arg)
{
  if (narg > 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairComb::coeff(int narg, char **arg)
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
  setup();

  n = atom->ntypes;

  // generate streitz-mintmire direct 1/r energy look-up table

  if (comm->me == 0 && screen) fprintf(screen,"Pair COMB:\n");
  if (comm->me == 0 && screen)
    fprintf(screen,"  generating Coulomb integral lookup table ...\n");
  sm_table();

  if (cor_flag && comm->me == 0 && screen)
    fprintf(screen,"  will apply over-coordination correction ...\n");
  if (!cor_flag&& comm->me == 0 && screen)
    fprintf(screen,"  will not apply over-coordination correction ...\n");

  // clear setflag since coeff() called once with I,J = * *

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

void PairComb::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style COMB requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style COMB requires newton pair on");
  if (!atom->q_flag)
    error->all(FLERR,"Pair style COMB requires atom attribute q");

  // ptr to QEQ fix

  //for (i = 0; i < modify->nfix; i++)
  //  if (strcmp(modify->fix[i]->style,"qeq") == 0) break;
  //if (i < modify->nfix) fixqeq = (FixQEQ *) modify->fix[i];
  //else fixqeq = NULL;

  // need a full neighbor list

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  // local Comb neighbor list
  // create pages if first time or if neighbor pgsize/oneatom has changed

  int create = 0;
  if (ipage == NULL) create = 1;
  if (pgsize != neighbor->pgsize) create = 1;
  if (oneatom != neighbor->oneatom) create = 1;

  if (create) {
    delete [] ipage;
    pgsize = neighbor->pgsize;
    oneatom = neighbor->oneatom;

    int nmypage = comm->nthreads;
    ipage = new MyPage<int>[nmypage];
    for (int i = 0; i < nmypage; i++)
      ipage[i].init(oneatom,pgsize);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairComb::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairComb::read_file(char *file)
{
  int params_per_line = 49;
  char **words = new char*[params_per_line+1];

  memory->sfree(params);
  params = NULL;
  nparams = 0;
  maxparam = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = open_potential(file);
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open COMB potential file %s",file);
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
      error->all(FLERR,"Incorrect format in COMB potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while (words[nwords++] = strtok(NULL," \t\n\r\f")) continue;

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
    params[nparams].c = atof(words[4]);
    params[nparams].d = atof(words[5]);
    params[nparams].h = atof(words[6]);
    params[nparams].powern = atof(words[7]);
    params[nparams].beta = atof(words[8]);
    params[nparams].lam21 = atof(words[9]);
    params[nparams].lam22 = atof(words[10]);
    params[nparams].bigb1 = atof(words[11]);
    params[nparams].bigb2 = atof(words[12]);
    params[nparams].bigr = atof(words[13]);
    params[nparams].bigd = atof(words[14]);
    params[nparams].lam11 = atof(words[15]);
    params[nparams].lam12 = atof(words[16]);
    params[nparams].biga1 = atof(words[17]);
    params[nparams].biga2 = atof(words[18]);
    params[nparams].plp1 = atof(words[19]);
    params[nparams].plp3 = atof(words[20]);
    params[nparams].plp6 = atof(words[21]);
    params[nparams].a123 = atof(words[22]);
    params[nparams].aconf= atof(words[23]);
    params[nparams].addrep = atof(words[24]);
    params[nparams].romigb = atof(words[25]);
    params[nparams].romigc = atof(words[26]);
    params[nparams].romigd = atof(words[27]);
    params[nparams].romiga = atof(words[28]);
    params[nparams].QL1 = atof(words[29]);
    params[nparams].QU1 = atof(words[30]);
    params[nparams].DL1 = atof(words[31]);
    params[nparams].DU1 = atof(words[32]);
    params[nparams].QL2 = atof(words[33]);
    params[nparams].QU2 = atof(words[34]);
    params[nparams].DL2 = atof(words[35]);
    params[nparams].DU2 = atof(words[36]);
    params[nparams].chi = atof(words[37]);
    params[nparams].dj  = atof(words[38]);
    params[nparams].dk  = atof(words[39]);
    params[nparams].dl  = atof(words[40]);
    params[nparams].dm  = atof(words[41]);
    params[nparams].esm1 = atof(words[42]);
    params[nparams].cmn1 = atof(words[43]);
    params[nparams].cml1 = atof(words[44]);
    params[nparams].cmn2 = atof(words[45]);
    params[nparams].cml2 = atof(words[46]);
    params[nparams].coulcut = atof(words[47]);
    params[nparams].hfocor = atof(words[48]);

    params[nparams].powermint = int(params[nparams].powerm);

    // parameter sanity checks

    if (params[nparams].lam11 < 0.0 || params[nparams].lam12 < 0.0 ||
        params[nparams].c < 0.0 || params[nparams].d < 0.0 ||
        params[nparams].powern < 0.0 || params[nparams].beta < 0.0 ||
        params[nparams].lam21 < 0.0 || params[nparams].lam22 < 0.0 ||
        params[nparams].bigb1< 0.0 || params[nparams].bigb2< 0.0 ||
        params[nparams].biga1< 0.0 || params[nparams].biga2< 0.0 ||
        params[nparams].bigr < 0.0 || params[nparams].bigd < 0.0 ||
        params[nparams].bigd > params[nparams].bigr ||
        params[nparams].powerm - params[nparams].powermint != 0.0 ||
        (params[nparams].powermint != 3 && params[nparams].powermint != 1) ||
        params[nparams].plp1 < 0.0 || params[nparams].plp3 < 0.0 ||
        params[nparams].plp6 < 0.0  ||
        params[nparams].a123 > 360.0 || params[nparams].aconf < 0.0 ||
        params[nparams].addrep < 0.0 || params[nparams].romigb < 0.0 ||
        params[nparams].romigc < 0.0 || params[nparams].romigd < 0.0 ||
        params[nparams].romiga < 0.0 ||
        params[nparams].QL1 > 0.0 || params[nparams].QU1 < 0.0 ||
        params[nparams].DL1 < 0.0 || params[nparams].DU1 > 0.0 ||
        params[nparams].QL2 > 0.0 || params[nparams].QU2 < 0.0 ||
        params[nparams].DL2 < 0.0 || params[nparams].DU2 > 0.0 ||
        params[nparams].chi < 0.0 ||
//        params[nparams].dj < 0.0 || params[nparams].dk < 0.0 ||
//        params[nparams].dl < 0.0 || params[nparams].dm < 0.0 ||
        params[nparams].esm1 < 0.0)
      error->all(FLERR,"Illegal COMB parameter");

    if (params[nparams].lam11 < params[nparams].lam21 ||
        params[nparams].lam12 < params[nparams].lam22 ||
        params[nparams].biga1< params[nparams].bigb1 ||
        params[nparams].biga2< params[nparams].bigb2)
      error->all(FLERR,"Illegal COMB parameter");

    nparams++;
  }

  delete [] words;
}

/* ---------------------------------------------------------------------- */

void PairComb::setup()
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
    params[m].rlm1 = 0.5*(params[m].lam11+params[m].lam12)*params[m].romigc;
    params[m].rlm2 = 0.5*(params[m].lam21+params[m].lam22)*params[m].romigd;

    params[m].Qo1 = (params[m].QU1+params[m].QL1)/2.0; // (A22)
    params[m].dQ1 = (params[m].QU1-params[m].QL1)/2.0; // (A21)
    params[m].aB1 = 1.0 /
      (1.0-pow(fabs(params[m].Qo1/params[m].dQ1),10.0)); // (A20)
    params[m].bB1 = pow(fabs(params[m].aB1),0.1)/params[m].dQ1; // (A19)
    params[m].nD1 = log(params[m].DU1/(params[m].DU1-params[m].DL1))/
                    log(params[m].QU1/(params[m].QU1-params[m].QL1));
    params[m].bD1 = (pow((params[m].DL1-params[m].DU1),(1.0/params[m].nD1)))/
                    (params[m].QU1-params[m].QL1);

    params[m].Qo2 = (params[m].QU2+params[m].QL2)/2.0; // (A22)
    params[m].dQ2 = (params[m].QU2-params[m].QL2)/2.0; // (A21)
    params[m].aB2 = 1.0 /
      (1.0-pow(fabs(params[m].Qo2/params[m].dQ2),10.0)); // (A20)
    params[m].bB2 = pow(fabs(params[m].aB2),0.1)/params[m].dQ2; // (A19)
    params[m].nD2 = log(params[m].DU2/(params[m].DU2-params[m].DL2))/
                    log(params[m].QU2/(params[m].QU2-params[m].QL2));
    params[m].bD2 = (pow((params[m].DL2-params[m].DU2),(1.0/params[m].nD2)))/
                    (params[m].QU2-params[m].QL2);

    params[m].lcut = params[m].coulcut;
    params[m].lcutsq = params[m].lcut*params[m].lcut;

    params[m].gamma = 1.0;        // for the change in pair_comb.h
  }

  // set cutmax to max of all params

  cutmax = cutmin = 0.0;
  cor_flag = 0;
  for (m = 0; m < nparams; m++) {
    if (params[m].cut > cutmax) cutmax = params[m].cut;
    if (params[m].lcut > cutmax) cutmax = params[m].lcut;
    if (params[m].cutsq > cutmin) cutmin = params[m].cutsq+0.2;
    if (params[m].hfocor > 0.0001) cor_flag = 1;
  }
}

/* ---------------------------------------------------------------------- */

void PairComb::repulsive(Param *param, double rsq, double &fforce,
                    int eflag, double &eng, double iq, double jq)
{
  double r,tmp_fc,tmp_fc_d,tmp_exp,Di,Dj;
  double bigA,Asi,Asj,vrcs,fvrcs,fforce_tmp;
  double rslp,rslp2,rslp4,arr1,arr2,fc2j,fc3j,fcp2j,fcp3j;

  double romi = param->addrep;
  double rrcs = param->bigr + param->bigd;

  r = sqrt(rsq);
  if (r > rrcs) return ;

  tmp_fc = comb_fc(r,param);
  tmp_fc_d = comb_fc_d(r,param);
  tmp_exp = exp(-param->rlm1 * r);

  arr1 = 2.22850; arr2 = 1.89350;
  fc2j = comb_fc2(r);
  fc3j = comb_fc3(r);
  fcp2j = comb_fc2_d(r);
  fcp3j = comb_fc3_d(r);

  Di = param->DU1 + pow(fabs(param->bD1*(param->QU1-iq)),param->nD1);
  Dj = param->DU2 + pow(fabs(param->bD2*(param->QU2-jq)),param->nD2);
  Asi = param->biga1 * exp(param->lam11*Di);
  Asj = param->biga2 * exp(param->lam12*Dj);

  if ( Asi > 0.0 && Asj > 0.0 )
    bigA = sqrt(Asi*Asj)*param->romiga;
  else
    bigA = 0.0;

  fforce = -bigA * tmp_exp * (tmp_fc_d - tmp_fc*param->rlm1) / r;

  // additional repulsion for TiO2 and HfO2 (switch by cor_flag)

  vrcs = 0.0; fvrcs = 0.0;
  if (romi > 0.0) {
    if (!cor_flag) {
      vrcs = romi * pow((1.0-r/rrcs),2.0);
      fvrcs= romi * 2.0 * (r/rrcs-1.0)/rrcs; }
    else if (cor_flag) {
      rslp = ((arr1-r)/(arr1-arr2));
      rslp2 = rslp * rslp; rslp4 = rslp2 * rslp2;
      vrcs = fc2j * fc3j * romi * ((50.0*rslp4-30.0*rslp2+4.50))/8.0;
      fvrcs = fcp2j*fcp3j*romi*rslp*(-25.0*rslp2+7.50)/(arr1-arr2);
    }
    fforce_tmp = fforce*vrcs - (tmp_fc * bigA * tmp_exp * fvrcs);
    fforce += fforce_tmp;
  }

  // eng = repulsive energy

  if (eflag) eng = (tmp_fc * bigA * tmp_exp)*(1.0+vrcs);
}

/* ---------------------------------------------------------------------- */

double PairComb::zeta(Param *param, double rsqij, double rsqik,
                         double *delrij, double *delrik)
{
  double rij,rik,costheta,arg,ex_delr;

  rij = sqrt(rsqij);
  if (rij > param->bigr+param->bigd) return 0.0;
  rik = sqrt(rsqik);
  costheta = vec3_dot(delrij,delrik) / (rij*rik);

  if (param->powermint == 3) arg = pow(param->rlm2 * (rij-rik),3.0);
  else arg = param->rlm2 * (rij-rik);

  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(arg);

  return comb_fc(rik,param) * comb_gijk(costheta,param) * ex_delr;
}

/* ----------------------------------------------------------------------
   Legendre polynomial bond angle correction to energy
------------------------------------------------------------------------- */

double PairComb::elp(Param *param, double rsqij, double rsqik,
                     double *delrij, double *delrik)
{
  if (param->aconf > 1.0e-6 || param->plp1 > 1.0e-6 ||
      param->plp3 > 1.0e-6 || param->plp6 > 1.0e-6) {
    double rij,rik,costheta,lp1,lp3,lp6;
    double rmu,rmu2,comtt,fcj,fck;
    double pplp1 = param->plp1, pplp3 = param->plp3, pplp6 = param->plp6;
    double c123 = cos(param->a123*MY_PI/180.0);

    // cos(theta) of the i-j-k
    // cutoff function of rik

    rij = sqrt(rsqij);
    rik = sqrt(rsqik);
    costheta = vec3_dot(delrij,delrik) / (rij*rik);
    fcj = comb_fc(rij,param);
    fck = comb_fc(rik,param);
    rmu = costheta;

    // Legendre Polynomial functions

    if (param->plp1 > 1.0e-6 || param->plp3 > 1.0e-6 || param->plp6 > 1.0e-6) {
      rmu2 = rmu*rmu;
      lp1 = rmu; lp3 = 0.5*(5.0*rmu2*rmu-3.0*rmu);
      lp6 = (231.0*rmu2*rmu2*rmu2-315.0*rmu2*rmu2+105.0*rmu2-5.0)/16.0;
      comtt = pplp1*lp1 + pplp3*lp3 + pplp6*lp6;
    } else comtt = 0.0;

    // bond-bending terms

    if (param->aconf>1e-4) {
      if (param->hfocor >= 0.0)
        comtt += param->aconf *(rmu-c123)*(rmu-c123);
      else if (param->hfocor < 0.0)
        comtt += param->aconf *(4.0-(rmu-c123)*(rmu-c123));
    }

    return 0.5 * fcj * fck * comtt;
  }

  return 0.0;
}

/* ----------------------------------------------------------------------
   Legendre polynomial bond angle correction to forces
------------------------------------------------------------------------- */

void PairComb::flp(Param *param, double rsqij, double rsqik,
                   double *delrij, double *delrik, double *drilp,
                   double *drjlp, double *drklp)
{
  double ffj1,ffj2,ffk1,ffk2;
  ffj1 = 0.0; ffj2 = 0.0; ffk1 = 0.0; ffk2 = 0.0;

  if (param->aconf > 1.0e-4 || param->plp1 > 1.0e-6 ||
      param->plp3 > 1.0e-6 || param->plp6 > 1.0e-6) {
    double rij,rik,costheta,lp1,lp1_d,lp3,lp3_d,lp6,lp6_d;
    double rmu,rmu2,comtt,comtt_d,com4k,com5,fcj,fck,fck_d;

    double pplp1 = param->plp1;
    double pplp3 = param->plp3;
    double pplp6 = param->plp6;
    double c123 = cos(param->a123*MY_PI/180.0);

    // fck_d = derivative of cutoff function

    rij = sqrt(rsqij); rik = sqrt(rsqik);
    costheta = vec3_dot(delrij,delrik) / (rij*rik);
    fcj = comb_fc(rij,param);
    fck = comb_fc(rik,param);
    fck_d = comb_fc_d(rik,param);
    rmu = costheta;

    // Legendre Polynomial functions and derivatives

    if (param->plp1 > 1.0e-6 || param->plp3 > 1.0e-6 || param->plp6 > 1.0e-6) {
      rmu2 = rmu*rmu;
      lp1 = rmu; lp3 = (2.5*rmu2*rmu-1.5*rmu);
      lp6 = (231.0*rmu2*rmu2*rmu2-315.0*rmu2*rmu2+105.0*rmu2-5.0)/16.0;
      lp1_d = 1.0;lp3_d = (7.5*rmu2-1.5);
      lp6_d = (1386.0*rmu2*rmu2*rmu-1260.0*rmu2*rmu+210.0)/16.0;
      comtt   = pplp1*lp1   + pplp3*lp3   + pplp6*lp6;
      comtt_d = pplp1*lp1_d + pplp3*lp3_d + pplp6*lp6_d;
    } else {
      comtt = 0.0;
      comtt_d = 0.0;
    }

    // bond-bending terms derivatives

    if (param->aconf > 1.0e-4) {
      if (param->hfocor >= 0.0) {
        comtt += param->aconf *(rmu-c123)*(rmu-c123);
        comtt_d += 2.0*param->aconf*(rmu-c123);
      } else if (param->hfocor < 0.0) {
        comtt += param->aconf *(4.0-(rmu-c123)*(rmu-c123));
        comtt_d += -2.0*param->aconf*(rmu-c123);
      }
    }

    com4k = 2.0 * fcj * fck_d * comtt;
    com5 = fcj * fck * comtt_d;

    ffj1 =-0.5*(com5/(rij*rik));
    ffj2 = 0.5*(com5*rmu/rsqij);
    ffk1 = ffj1;
    ffk2 = 0.5*(-com4k/rik+com5*rmu/rsqik);

  } else {
    ffj1 = 0.0; ffj2 = 0.0;
    ffk1 = 0.0; ffk2 = 0.0;
  }

  // j-atom

  vec3_scale(ffj1,delrik,drjlp);             // (k,x[],y[]), y[]=k*x[]
  vec3_scaleadd(ffj2,delrij,drjlp,drjlp);   // (k,x[],y[],z[]), z[]=k*x[]+y[]

  // k-atom

  vec3_scale(ffk1,delrij,drklp);
  vec3_scaleadd(ffk2,delrik,drklp,drklp);

  // i-atom

  vec3_add(drjlp,drklp,drilp);                    // (x[],y[],z[]), z[]=x[]+y[]
  vec3_scale(-1.0,drilp,drilp);
}

/* ---------------------------------------------------------------------- */

void PairComb::force_zeta(Param *param, int eflag, int i, int j, double rsq,
                double zeta_ij, double iq, double jq, double &fforce,
                double &prefactor, double &eng)
{
  double r,fa,fa_d,bij;

  r = sqrt(rsq);
  if (r > param->bigr+param->bigd) return;
  fa = comb_fa(r,param,iq,jq);
  fa_d = comb_fa_d(r,param,iq,jq);
  bij = comb_bij(zeta_ij,param);
  bbij[i][j] = bij;

  // force
  fforce = 0.5*bij*fa_d / r;
  prefactor = -0.5*fa * comb_bij_d(zeta_ij,param);

  // eng = attractive energy
  if (eflag) eng = 0.5*bij*fa;
}

/* ---------------------------------------------------------------------- */

double PairComb::comb_fc(double r, Param *param)
{
  double comb_R = param->bigr;
  double comb_D = param->bigd;

  if (r < comb_R-comb_D) return 1.0;
  if (r > comb_R+comb_D) return 0.0;
  return 0.5*(1.0 - sin(MY_PI2*(r - comb_R)/comb_D));
}

/* ---------------------------------------------------------------------- */

double PairComb::comb_fc_d(double r, Param *param)
{
  double comb_R = param->bigr;
  double comb_D = param->bigd;

  if (r < comb_R-comb_D) return 0.0;
  if (r > comb_R+comb_D) return 0.0;
  return -(MY_PI4/comb_D) * cos(MY_PI2*(r - comb_R)/comb_D);
}

/* ---------------------------------------------------------------------- */

double PairComb::comb_fc2(double r)
{
  double comb_R = 1.89350;
  double comb_D = comb_R + 0.050;

  if (r < comb_R) return 0.0;
  if (r > comb_D) return 1.0;
  return 0.5*(1.0 + cos(MY_PI*(r - comb_R)/(comb_D-comb_R)));
}

/* ---------------------------------------------------------------------- */

double PairComb::comb_fc2_d(double r)
{
  double comb_R = 1.89350;
  double comb_D = comb_R + 0.050;

  if (r < comb_R) return 0.0;
  if (r > comb_D) return 0.0;
  return -(MY_PI2/(comb_D-comb_R)) * sin(MY_PI*(r - comb_R)/(comb_D-comb_R));
}

/* ---------------------------------------------------------------------- */

double PairComb::comb_fc3(double r)
{
  double comb_R = 2.51350;
  double comb_D = comb_R + 0.050;

  if (r < comb_R) return 1.0;
  if (r > comb_D) return 0.0;
  return 0.5*(1.0 + cos(MY_PI*(r - comb_R)/(comb_D-comb_R)));
}

/* ---------------------------------------------------------------------- */

double PairComb::comb_fc3_d(double r)
{
  double comb_R = 2.51350;
  double comb_D = comb_R + 0.050;

  if (r < comb_R) return 0.0;
  if (r > comb_D) return 0.0;
  return -(MY_PI2/(comb_D-comb_R)) * sin(MY_PI*(r - comb_R)/(comb_D-comb_R));
}

/* ---------------------------------------------------------------------- */

double PairComb::self(Param *param, double qi, double selfpot)
{
 double self_tmp, cmin, cmax, qmin, qmax;
 double s1=param->chi, s2=param->dj, s3=param->dk, s4=param->dl, s5=param->dm;

 self_tmp = 0.0;
 qmin = param->QL1*0.90;
 qmax = param->QU1*0.90;
 cmin = cmax = 1000.0;

 self_tmp = qi*(s1+qi*(s2+selfpot+qi*(s3+qi*(s4+qi*qi*s5))));

 if (qi < qmin) self_tmp += cmin * pow((qi-qmin),4.0);
 if (qi > qmax) self_tmp += cmax * pow((qi-qmax),4.0);

 return self_tmp;
}

/* ---------------------------------------------------------------------- */

double PairComb::comb_fa(double r, Param *param, double iq, double jq)
{
  double bigB,Bsi,Bsj;
  double qi,qj,Di,Dj;

  if (r > param->bigr + param->bigd) return 0.0;
  qi = iq; qj = jq;
  Di = Dj = Bsi = Bsj = bigB = 0.0;
  Di = param->DU1 + pow(fabs(param->bD1*(param->QU1-qi)),param->nD1);
  Dj = param->DU2 + pow(fabs(param->bD2*(param->QU2-qj)),param->nD2);
  Bsi = param->bigb1 * exp(param->lam21*Di)*
       (param->aB1-fabs(pow(param->bB1*(qi-param->Qo1),10.0)));
  Bsj = param->bigb2 * exp(param->lam22*Dj)*
       (param->aB2-fabs(pow(param->bB2*(qj-param->Qo2),10.0)));
  if (Bsi > 0.0 && Bsj > 0.0) bigB = sqrt(Bsi*Bsj)*param->romigb;
  else bigB = 0.0;

  return -bigB * exp(-param->rlm2 * r) * comb_fc(r,param);
}

/* ---------------------------------------------------------------------- */

double PairComb::comb_fa_d(double r, Param *param, double iq, double jq)
{
  double bigB,Bsi,Bsj;
  double qi,qj,Di,Dj;

  if (r > param->bigr + param->bigd) return 0.0;
  qi = iq; qj = jq;
  Di = Dj = Bsi = Bsj = bigB = 0.0;
  Di = param->DU1 + pow(fabs(param->bD1*(param->QU1-qi)),param->nD1);
  Dj = param->DU2 + pow(fabs(param->bD2*(param->QU2-qj)),param->nD2);
  Bsi = param->bigb1 * exp(param->lam21*Di)*
       (param->aB1-fabs(pow(param->bB1*(qi-param->Qo1),10.0)));
  Bsj = param->bigb2 * exp(param->lam22*Dj)*
       (param->aB2-fabs(pow(param->bB2*(qj-param->Qo2),10.0)));
  if (Bsi > 0.0 && Bsj > 0.0) bigB = sqrt(Bsi*Bsj)*param->romigb;
  else bigB = 0.0;

  return bigB * exp(-param->rlm2 * r) *
    (param->rlm2 * comb_fc(r,param) - comb_fc_d(r,param));
}

/* ---------------------------------------------------------------------- */

double PairComb::comb_bij(double zeta, Param *param)
{
  double tmp = param->beta * zeta;
  if (tmp > param->c1) return 1.0/sqrt(tmp);
  if (tmp > param->c2)
    return (1.0 - pow(tmp,-1.0*param->powern) / (2.0*param->powern))/sqrt(tmp);
  if (tmp < param->c4) return 1.0;
  if (tmp < param->c3)
    return 1.0 - pow(tmp,param->powern)/(2.0*param->powern);
  return pow(1.0 + pow(tmp,param->powern), -1.0/(2.0*param->powern));
}

/* ---------------------------------------------------------------------- */

double PairComb::comb_bij_d(double zeta, Param *param)
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

void PairComb::attractive(Param *param, double prefactor,
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

  comb_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk,param);
}

/* ---------------------------------------------------------------------- */

void PairComb::comb_zetaterm_d(double prefactor, double *rij_hat, double rij,
                               double *rik_hat, double rik, double *dri,
                               double *drj, double *drk, Param *param)
{
  double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = comb_fc(rik,param);
  dfc = comb_fc_d(rik,param);
  if (param->powermint == 3) tmp = pow(param->rlm2 * (rij-rik),3.0);
  else tmp = param->rlm2 * (rij-rik);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp); // ex_delr is Ygexp

  if (param->powermint == 3)
    ex_delr_d = 3.0*pow(param->rlm2,3.0) * pow(rij-rik,2.0)*ex_delr; // com3
  else ex_delr_d = param->rlm2 * ex_delr; // com3

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = comb_gijk(cos_theta,param);
  gijk_d = comb_gijk_d(cos_theta,param);
  costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

  // compute the derivative wrt Ri
  // dri = -dfc*gijk*ex_delr*rik_hat;
  // dri += fc*gijk_d*ex_delr*dcosdri;
  // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);
  // (k,x[],y[]), y[]=k*x[]
  // (k,x[],y[],z[]), z[]=k*x[]+y[]

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

void PairComb::costheta_d(double *rij_hat, double rij,
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

void PairComb::sm_table()
{
  int i,j,k,m,nntypes,ncoul;
  int inty, itype, jtype;
  int iparam_i, iparam_ij, iparam_ji;
  double r,dra,drin,rc,z,zr,zrc,ea,eb,ea3,eb3,alf;
  double exp2er,exp2ersh,fafash,dfafash,F1,dF1,ddF1,E1,E2,E3,E4;
  double exp2ear,exp2ebr,exp2earsh,exp2ebrsh,fafbsh,dfafbsh;

  int n = atom->ntypes;
  int nmax = atom->nmax;

  dra  = 0.001;  // lookup table step size
  drin = 0.1;    // starting distance of 1/r
  rc = cutmax;
  alf = 0.20;

  nntypes = int((n+1)*n/2); // interaction types
  ncoul = int((rc-drin)/dra)+1;

  // allocate arrays

  memory->create(intype,n,n,"pair:intype");
  memory->create(fafb,ncoul,nntypes,"pair:fafb");
  memory->create(dfafb,ncoul,nntypes,"pair:dfafb");
  memory->create(ddfafb,ncoul,nntypes,"pair:ddfafb");
  memory->create(phin,ncoul,nntypes,"pair:phin");
  memory->create(dphin,ncoul,nntypes,"pair:dphin");
  memory->create(erpaw,25000,2,"pair:erpaw");
  memory->create(NCo,nmax,"pair:NCo");
  memory->create(bbij,nmax,MAXNEIGH,"pair:bbij");
  memory->create(sht_num,nmax,"pair:sht_num");
  sht_first = (int **) memory->smalloc(nmax*sizeof(int *),"pair:sht_first");

  // set interaction number: 0-0=0, 1-1=1, 0-1=1-0=2

  m = 0; k = n;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (j == i) {
        intype[i][j] = m;
        m += 1;
      } else if (j != i && j > i) {
        intype[i][j] = k;
        k += 1;
      } else if (j != i && j < i) {
        intype[i][j] = intype[j][i];
      }
    }
  }

  // default arrays to zero

  for (i = 0; i < ncoul; i ++) {
    for (j = 0; j < nntypes; j ++) {
      fafb[i][j] = 0.0;
      dfafb[i][j] = 0.0;
      ddfafb[i][j] = 0.0;
      phin[i][j] = 0.0;
      dphin[i][j] = 0.0;
    }
  }

  // direct 1/r energy with Slater 1S orbital overlap

  for (i = 0; i < n; i++) {
    r = drin;
    itype = params[i].ielement;
    iparam_i = elem2param[itype][itype][itype];
    z = params[iparam_i].esm1;
    for (j = 0; j < ncoul; j++) {
      exp2er = exp(-2.0 * z * r);
      phin[j][i] = 1.0 - exp2er * (1.0 + 2.0 * z * r * (1.0 + z * r));
      dphin[j][i] = (4.0 * exp2er * z * z * z * r * r);
      r += dra;
    }
  }

  for (i = 0; i < n; i ++) {
    for (j = 0; j < n; j ++) {
      r = drin;
      if (j == i) {
        itype = params[i].ielement;
        inty = intype[itype][itype];
        iparam_i = elem2param[itype][itype][itype];
        z = params[iparam_i].esm1;
        zrc = z * rc;
        exp2ersh = exp(-2.0 * zrc);
        fafash = -exp2ersh * (1.0 / rc +
                              z * (11.0/8.0 + 3.0/4.0*zrc + zrc*zrc/6.0));
        dfafash = exp2ersh * (1.0/(rc*rc) + 2.0*z/rc +
                              z*z*(2.0 + 7.0/6.0*zrc + zrc*zrc/3.0));
        for (k = 0; k < ncoul; k ++) {
          zr = z * r;
          exp2er = exp(-2.0*zr);
          F1 = -exp2er * (1.0 / r +
                          z * (11.0/8.0 + 3.0/4.0*zr + zr*zr/6.0));
          dF1 = exp2er * (1.0/(r*r) + 2.0*z/r +
                          z*z*(2.0 + 7.0/6.0*zr + zr*zr/3.0));
          ddF1 = -exp2er * (2.0/(r*r*r) + 4.0*z/(r*r) -
                            z*z*z/3.0*(17.0/2.0 + 5.0*zr + 2.0*zr*zr));
          fafb[k][inty] = F1-fafash-(r-rc)*dfafash;
          dfafb[k][inty] = (dF1 - dfafash);
          ddfafb[k][inty] = ddF1;
          r += dra;
        }
      } else if (j != i) {
        itype = params[i].ielement;
        jtype = params[j].ielement;
        inty = intype[itype][jtype];
        iparam_ij = elem2param[itype][jtype][jtype];
        ea = params[iparam_ij].esm1;
        ea3 = ea*ea*ea;
        iparam_ji = elem2param[jtype][itype][itype];
        eb = params[iparam_ji].esm1;
        eb3 = eb*eb*eb;
        E1 = ea*eb3*eb/((ea+eb)*(ea+eb)*(ea-eb)*(ea-eb));
        E2 = eb*ea3*ea/((ea+eb)*(ea+eb)*(eb-ea)*(eb-ea));
        E3 = (3.0*ea*ea*eb3*eb-eb3*eb3) /
          ((ea+eb)*(ea+eb)*(ea+eb)*(ea-eb)*(ea-eb)*(ea-eb));
        E4 = (3.0*eb*eb*ea3*ea-ea3*ea3) /
          ((ea+eb)*(ea+eb)*(ea+eb)*(eb-ea)*(eb-ea)*(eb-ea));
        exp2earsh = exp(-2.0*ea*rc);
        exp2ebrsh = exp(-2.0*eb*rc);
        fafbsh = -exp2earsh*(E1 + E3/rc)-exp2ebrsh*(E2 + E4/rc);
        dfafbsh =
          exp2earsh*(2.0*ea*(E1+E3/rc)+E3/(rc*rc)) +
          exp2ebrsh*(2.0*eb*(E2+E4/rc)+E4/(rc*rc));
        for (k = 0; k < ncoul; k ++) {
          exp2ear = exp(-2.0*ea*r);
          exp2ebr = exp(-2.0*eb*r);
          fafb[k][inty] =
            - exp2ear*(E1+E3/r) - exp2ebr*(E2+E4/r)
            - fafbsh - (r-rc) * dfafbsh;
          dfafb[k][inty] = (exp2ear*(2.0*ea*(E1+E3/r) + E3/(r*r))
                            + exp2ebr*(2.0*eb*(E2+E4/r) + E4/(r*r))- dfafbsh);
          ddfafb[k][inty] = (- exp2ear*(E3/(r*r)*(1.0/r+2.0*ea/r+2.0/(r*r))
                                        + 2.0*ea*(E1+E3/r))-
                             exp2ebr*(E4/(r*r)
                                      *(1.0/r+2.0*eb/r+2.0/(r*r)) +
                                      2.0*eb*(E2+E4/r)));
          r += dra;
        }
      }
    }
  }

  for (i = 0; i < 25000; i ++) {
    r = dra * i + drin;
    erpaw[i][0] = erfc(r*alf);
    erpaw[i][1] = exp(-r*r*alf*alf);
  }
}

/* ---------------------------------------------------------------------- */

void PairComb::potal_calc(double &calc1, double &calc2, double &calc3)
{
  double alf,rcoul,esucon;
  int m;

  rcoul = 0.0;
  for (m = 0; m < nparams; m++)
    if (params[m].lcut > rcoul) rcoul = params[m].lcut;

  alf = 0.20;
  esucon = force->qqr2e;

  calc2 = (erfc(rcoul*alf)/rcoul/rcoul+2.0*alf/MY_PIS*
           exp(-alf*alf*rcoul*rcoul)/rcoul)*esucon/rcoul;
  calc3 = (erfc(rcoul*alf)/rcoul)*esucon;
  calc1 = -(alf/MY_PIS*esucon+calc3*0.5);
}

/* ---------------------------------------------------------------------- */

void PairComb::tri_point(double rsq, int &mr1, int &mr2,
                         int &mr3, double &sr1, double &sr2,
                         double &sr3, int &itype)
{
 double r, rin, dr, dd, rr1, rridr, rridr2;

 rin = 0.10; dr = 0.0010;
 r = sqrt(rsq);
 if (r < rin + 2.0*dr) r = rin + 2.0*dr;
 if (r > cutmax - 2.0*dr) r = cutmax - 2.0*dr;
 rridr = (r-rin)/dr;

 mr1 = int(rridr)-1;
 dd = rridr - float(mr1);
 if (dd > 0.5) mr1 += 1;
 mr2 = mr1 + 1;
 mr3 = mr2 + 1;

 rr1 = float(mr1)*dr;
 rridr = (r - rin - rr1)/dr;
 rridr2 = rridr * rridr;

 sr1 = (rridr2 - rridr) * 0.50;
 sr2 = 1.0 - rridr2;
 sr3 = (rridr2 + rridr) * 0.50;
}

/* ---------------------------------------------------------------------- */

void PairComb::direct(int inty, int mr1, int mr2, int mr3, double rsq,
                      double sr1, double sr2, double sr3,
                      double iq, double jq,
                      double potal, double fac11, double fac11e,
                      double &pot_tmp, double &pot_d)
{
 double r,erfcc,fafbn1,potij,sme2,esucon;
 double r3,erfcd,dfafbn1,smf2,dvdrr,alf,alfdpi;

 r = sqrt(rsq);
 r3 = r * rsq;
 alf = 0.20;
 alfdpi = 2.0*alf/MY_PIS;
 esucon = force->qqr2e;
 pot_tmp = 0.0;
 pot_d = 0.0;

 // 1/r energy

 erfcc = sr1*erpaw[mr1][0] + sr2*erpaw[mr2][0] + sr3*erpaw[mr3][0];
 fafbn1= sr1*fafb[mr1][inty] + sr2*fafb[mr2][inty] + sr3*fafb[mr3][inty];
 potij = (erfcc/r * esucon - fac11e);
 sme2 = potij + fafbn1 * esucon;
 pot_tmp = 1.0 * iq * jq *sme2;

 // 1/r force (wrt r)

 erfcd = sr1*erpaw[mr1][1] + sr2*erpaw[mr2][1] + sr3*erpaw[mr3][1];
 dfafbn1= sr1*dfafb[mr1][inty] + sr2*dfafb[mr2][inty] + sr3*dfafb[mr3][inty];
 dvdrr = (erfcc/r3+alfdpi*erfcd/rsq)*esucon-fac11;
 smf2 = dvdrr - dfafbn1 * esucon/r;
 pot_d =  1.0 * iq * jq * smf2;

}

/* ---------------------------------------------------------------------- */

void PairComb::field(Param *param, double rsq, double iq,double jq,
                     double &vionij,double &fvionij)
{
 double r,r5,r6,rc,rc5,rc6,rf5,drf6,smpn,smpl,rfx1,rfx2;
 double cmi1,cmi2,cmj1,cmj2;

 r = sqrt(rsq);
 r5 = r*r*r*r*r;
 r6 = r5 * r;
 rc = param->lcut;
 rc5 = rc*rc*rc*rc*rc;
 rc6 = rc5 * rc;
 cmi1 = param->cmn1;
 cmi2 = param->cmn2;
 cmj1 = param->cml1;
 cmj2 = param->cml2;
 rf5 = 1.0/r5 - 1.0/rc5 + 5.0*(r-rc)/rc6;
 drf6 = 5.0/rc6 - 5.0/r6;

 // field correction energy

 smpn = rf5*jq*(cmi1+jq*cmi2);
 smpl = rf5*iq*(cmj1+iq*cmj2);
 vionij += 1.0 * (smpn + smpl);

 // field correction force

 rfx1 = jq*drf6*(cmi1+jq*cmi2)/r;
 rfx2 = iq*drf6*(cmj1+iq*cmj2)/r;
 fvionij -= 1.0 * (rfx1 + rfx2);
}

/* ---------------------------------------------------------------------- */

double PairComb::yasu_char(double *qf_fix, int &igroup)
{
  int i,j,ii,jj,jnum;
  int itype,jtype,iparam_i,iparam_ij;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp;
  double rsq1,delr1[3];
  int *ilist,*jlist,*numneigh,**firstneigh;
  double iq,jq,fqi,fqj,fqij,fqjj;
  double potal,fac11,fac11e,sr1,sr2,sr3;
  int mr1,mr2,mr3,inty,nj;


  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  tagint *tag = atom->tag;

  int inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  int *mask = atom->mask;
  int groupbit = group->bitmask[igroup];

  qf = qf_fix;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit)
      qf[i] = 0.0;
  }

  // communicating charge force to all nodes, first forward then reverse

  comm->forward_comm_pair(this);

  // self energy correction term: potal

  potal_calc(potal,fac11,fac11e);

  // loop over full neighbor list of my atoms

  fqi = fqj = fqij = fqjj = 0.0;

  for (ii = 0; ii < inum; ii ++) {
    i = ilist[ii];
    itag = tag[i];
    nj = 0;
    if (mask[i] & groupbit) {
      itype = map[type[i]];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      iq = q[i];
      iparam_i = elem2param[itype][itype][itype];

      // charge force from self energy

      fqi = qfo_self(&params[iparam_i],iq,potal);

      // two-body interactions

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
        jq = q[j];

        delr1[0] = x[j][0] - xtmp;
        delr1[1] = x[j][1] - ytmp;
        delr1[2] = x[j][2] - ztmp;
        rsq1 = vec3_dot(delr1,delr1);

        iparam_ij = elem2param[itype][jtype][jtype];

        // long range q-dependent

        if (rsq1 > params[iparam_ij].lcutsq) continue;

        inty = intype[itype][jtype];

        // polynomial three-point interpolation

        tri_point(rsq1,mr1,mr2,mr3,sr1,sr2,sr3,itype);

        // 1/r charge forces

        qfo_direct(inty,mr1,mr2,mr3,rsq1,sr1,sr2,sr3,fac11e,fqij);
        fqi += jq * fqij;  qf[j] += iq * fqij;

        // field correction to self energy and charge force

        qfo_field(&params[iparam_ij],rsq1,iq,jq,fqij,fqjj);
        fqi += fqij;
        qf[j] += fqjj;
      }

        // three-body interactions

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        jtype = map[type[j]];
        jq = q[j];

        delr1[0] = x[j][0] - xtmp;
        delr1[1] = x[j][1] - ytmp;
        delr1[2] = x[j][2] - ztmp;
        rsq1 = vec3_dot(delr1,delr1);

        iparam_ij = elem2param[itype][jtype][jtype];

        if (rsq1 > params[iparam_ij].cutsq) continue;
        nj ++;

        // charge force in Aij and Bij

        qfo_short(&params[iparam_ij],i,nj,rsq1,iq,jq,fqij,fqjj);
        fqi += fqij;  qf[j] += fqjj;
      }
      qf[i] += fqi;
    }
  }

  comm->reverse_comm_pair(this);

  // sum charge force on each node and return it

  double eneg = 0.0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit)
      eneg += qf[i];
  }
  double enegtot;
  MPI_Allreduce(&eneg,&enegtot,1,MPI_DOUBLE,MPI_SUM,world);
  return enegtot;
}

/* ---------------------------------------------------------------------- */

double PairComb::qfo_self(Param *param, double qi, double selfpot)
{
 double self_d,cmin,cmax,qmin,qmax;
 double s1 = param->chi;
 double s2 = param->dj;
 double s3 = param->dk;
 double s4 = param->dl;
 double s5 = param->dm;

 self_d = 0.0;
 qmin = param->QL1*0.90;
 qmax = param->QU1*0.90;
 cmin = cmax = 1000.0;

 self_d = s1+qi*(2.0*(s2+selfpot)+qi*(3.0*s3+qi*(4.0*s4+qi*qi*6.0*s5)));

 if (qi < qmin) {
   // char str[128];
   // sprintf(str,"Pair COMB charge %.10f with force %.10f hit min barrier",
   // qi,self_d);
   // error->warning(FLERR,str,0);
   self_d += 4.0 * cmin * pow((qi-qmin),3.0);
 }
 if (qi > qmax) {
   // char str[128];
   // sprintf(str,"Pair COMB charge %.10f with force %.10f hit max barrier",
   //           qi,self_d);
   // error->warning(FLERR,str,0);
   self_d += 4.0 * cmax * pow((qi-qmax),3.0);
 }

 return self_d;
}

/* ---------------------------------------------------------------------- */

void PairComb::qfo_direct(int inty, int mr1, int mr2, int mr3,
                          double rsq, double sr1, double sr2,
                          double sr3, double fac11e, double &fqij)
{
 double r, erfcc, fafbn1, vm, esucon;

 r = sqrt(rsq);
 esucon=force->qqr2e;

 // 1/r force (wrt q)

 erfcc = sr1*erpaw[mr1][0]   + sr2*erpaw[mr2][0]   + sr3*erpaw[mr3][0];
 fafbn1= sr1*fafb[mr1][inty] + sr2*fafb[mr2][inty] + sr3*fafb[mr3][inty];
 vm = (erfcc/r * esucon - fac11e);
 fqij = 1.0 * (vm+esucon*fafbn1);
}

/* ---------------------------------------------------------------------- */

void PairComb::qfo_field(Param *param, double rsq,double iq,double jq,
                         double &fqij, double &fqjj)
{
 double r,r5,r6,rc,rc5,rc6;
 double cmi1,cmi2,cmj1,cmj2,rf5;

 fqij = fqjj = 0.0;
 r  = sqrt(rsq);
 r5 = r*r*r*r*r;
 r6 = r5 * r;
 rc = param->lcut;
 rc5 = rc*rc*rc*rc*rc;
 rc6 = rc5 * rc;
 cmi1 = param->cmn1;
 cmi2 = param->cmn2;
 cmj1 = param->cml1;
 cmj2 = param->cml2;
 rf5 = 1.0/r5 - 1.0/rc5 + 5.0*(r-rc)/rc6;

 // field correction charge force

 fqij = 1.0 * rf5 * (cmj1 + 2.0 * iq * cmj2);
 fqjj = 1.0 * rf5 * (cmi1 + 2.0 * jq * cmi2);
}

/* ---------------------------------------------------------------------- */

void PairComb::qfo_short(Param *param, int i, int j, double rsq,
                         double iq, double jq, double &fqij, double &fqjj)
{
  double r,tmp_fc,tmp_fc_d,tmp_exp1,tmp_exp2;
  double bigA,Asi,Asj,vrcs;
  double romi = param->addrep,rrcs = param->bigr + param->bigd;
  double qi,qj,Di,Dj,bigB,Bsi,Bsj;
  double QUchi,QOchi,QUchj,QOchj,YYDiqp,YYDjqp;
  double YYAsiqp,YYAsjqp,YYBsiqp,YYBsjqp;
  double caj,cbj,bij,cfqr,cfqs;
  double romie = param->romiga;
  double romib = param->romigb;
  double ca1,ca2,ca3,ca4;
  double rslp,rslp2,rslp4,arr1,arr2,fc2j,fc3j,fcp2j,fcp3j;

  qi = iq; qj = jq; r = sqrt(rsq);
  Di = Dj = Asi = Asj = bigA = Bsi = Bsj = bigB = 0.0;
  QUchi = QOchi = QUchj = QOchj = YYDiqp = YYDjqp =0.0;
  YYAsiqp = YYAsjqp = YYBsiqp = YYBsjqp = 0.0;
  caj = cbj = vrcs = cfqr = cfqs = 0.0;

  tmp_fc = comb_fc(r,param);
  tmp_fc_d = comb_fc_d(r,param);
  tmp_exp1 = exp(-param->rlm1 * r);
  tmp_exp2 = exp(-param->rlm2 * r);
  bij = bbij[i][j]; //comb_bij(zeta_ij,param);

  arr1 = 2.22850; arr2 = 1.89350;
  fc2j = comb_fc2(r);
  fc3j = comb_fc3(r);
  fcp2j = comb_fc2_d(r);
  fcp3j = comb_fc3_d(r);

  vrcs = 0.0;
  if (romi > 0.0) {
    if (!cor_flag) vrcs = romi * pow((1.0-r/rrcs),2.0);
    else if (cor_flag) {
      rslp = ((arr1-r)/(arr1-arr2));
      rslp2 = rslp * rslp; rslp4 = rslp2 * rslp2;
      vrcs = fc2j * fc3j * romi * ((50.0*rslp4-30.0*rslp2+4.50))/8.0;
    }
  }

  Di = param->DU1 + pow(fabs(param->bD1*(param->QU1-qi)),param->nD1);
  Dj = param->DU2 + pow(fabs(param->bD2*(param->QU2-qj)),param->nD2);

  Asi = param->biga1 * exp(param->lam11*Di);
  Asj = param->biga2 * exp(param->lam12*Dj);
  Bsi = param->bigb1 * exp(param->lam21*Di)*
    (param->aB1-fabs(pow(param->bB1*(qi-param->Qo1),10.0)));
  Bsj = param->bigb2 * exp(param->lam22*Dj)*
    (param->aB2-fabs(pow(param->bB2*(qj-param->Qo2),10.0)));

  QUchi = (param->QU1-qi)*param->bD1;
  QUchj = (param->QU2-qj)*param->bD2;
  QOchi = (qi-param->Qo1)*param->bB1;
  QOchj = (qj-param->Qo2)*param->bB2;

  if (QUchi == 0.0) YYDiqp = 0.0;
  else YYDiqp = -param->nD1 * QUchi * param->bD1 *
         pow(fabs(QUchi),(param->nD1-2.0));

  if (QUchj == 0.0) YYDjqp = 0.0;
  else YYDjqp = -param->nD2 * QUchj * param->bD2 *
         pow(fabs(QUchj),(param->nD2-2.0));

  YYAsiqp = Asi * param->lam11 * YYDiqp;
  YYAsjqp = Asj * param->lam12 * YYDjqp;

  if (QOchi == 0.0)
    YYBsiqp=Bsi*param->lam21*YYDiqp;
  else
    YYBsiqp=Bsi*param->lam21*YYDiqp-param->bigb1*exp(param->lam21*Di)*
      10.0*QOchi*param->bB1*pow(fabs(QOchi),(10.0-2.0));

  if (QOchj == 0.0)
    YYBsjqp=Bsj*param->lam22*YYDjqp;
  else
    YYBsjqp=Bsj*param->lam22*YYDjqp-param->bigb2*exp(param->lam22*Dj)*
      10.0*QOchj*param->bB2*pow(fabs(QOchj),(10.0-2.0));

  if (Asi > 0.0 && Asj > 0.0) caj = 1.0/(2.0*sqrt(Asi*Asj)) * romie;
  else caj = 0.0;

  if (Bsi > 0.0 && Bsj > 0.0) cbj = 1.0/(2.0*sqrt(Bsi*Bsj)) * romib ;
  else cbj = 0.0;

  cfqr =  0.50 * tmp_fc * (1.0 + vrcs); // 0.5 b/c full atom loop
  cfqs = -0.50 * tmp_fc *  bij;

  ca1 = Asj * caj * YYAsiqp;
  ca2 = Bsj * cbj * YYBsiqp;
  ca3 = Asi * caj * YYAsjqp;
  ca4 = Bsi * cbj * YYBsjqp;

  fqij  = cfqr * tmp_exp1 * ca1;
  fqij += cfqs * tmp_exp2 * ca2;
  fqjj  = cfqr * tmp_exp1 * ca3;
  fqjj += cfqs * tmp_exp2 * ca4;
}

/* ---------------------------------------------------------------------- */

void PairComb::Over_cor(Param *param, double rsq1, int NCoi,
                        double &Eov, double &Fov)
{
  double ECo,BCo,tmp_fc,tmp_fc_d;
  double r = sqrt(rsq1);
  int NCon = NCoi - 7;

  tmp_fc = comb_fc(r,param);
  tmp_fc_d = comb_fc(r,param);
  Eov = 0.0; Fov = 0.0;
  ECo = param->hfocor;
  BCo = 0.1;

  if (NCon >= 0.20) {
    Eov = tmp_fc * ECo * NCon/(1.0+exp(BCo*NCon));
    Fov = -(tmp_fc_d*Eov + tmp_fc*ECo/(1.0+exp(BCo*NCon)) -
            (tmp_fc*ECo*NCon*BCo*exp(BCo*NCon)) /
            ((1.0+exp(BCo*NCon))*(1.0+exp(BCo*NCon))));
    Fov /= r;
  }
}

/* ---------------------------------------------------------------------- */

int PairComb::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i ++) {
    j = list[i];
    buf[m++] = qf[j];
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

void PairComb::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n ;
  for (i = first; i < last; i++) qf[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

int PairComb::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = qf[i];
  return 1;
}

/* ---------------------------------------------------------------------- */

void PairComb::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    qf[j] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void PairComb::Short_neigh()
{
  int nj,itype,jtype,iparam_ij;
  int inum,jnum,i,j,ii,jj;
  int *neighptrj,*ilist,*jlist,*numneigh;
  int **firstneigh;
  double xtmp,ytmp,ztmp,rr,rsq,delrj[3];

  double **x = atom->x;
  int *type  = atom->type;
  int nlocal = atom->nlocal;
  int ntype = atom->ntypes;

  if (atom->nmax > nmax) {
    memory->sfree(sht_first);
    nmax = atom->nmax;
    sht_first = (int **) memory->smalloc(nmax*sizeof(int *),
                                         "pair:sht_first");
    memory->grow(sht_num,nmax,"pair:sht_num");
    memory->grow(NCo,nmax,"pair:NCo");
    memory->grow(bbij,nmax,MAXNEIGH,"pair:bbij");
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // create Comb neighbor list

  ipage->reset();

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    nj = 0;
    neighptrj = ipage->vget();

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      delrj[0] = xtmp - x[j][0];
      delrj[1] = ytmp - x[j][1];
      delrj[2] = ztmp - x[j][2];
      rsq = vec3_dot(delrj,delrj);

      if (rsq > cutmin) continue;
      neighptrj[nj++] = j;
    }

    sht_first[i] = neighptrj;
    sht_num[i] = nj;
    ipage->vgot(nj);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairComb::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  bytes += nmax * sizeof(int);
  bytes += nmax * sizeof(int *);

  for (int i = 0; i < comm->nthreads; i++)
    bytes += ipage[i].size();

  bytes += nmax * sizeof(int);
  bytes += MAXNEIGH*nmax * sizeof(double);
  return bytes;
}
