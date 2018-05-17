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
   Contributing author: Dave Humbird (DWH Process Consulting)
   Starting point provided by Ray Shan (Materials Design)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "pair_tersoffHG.h"

#include "atom.h"
#include "update.h"
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
#define TOL 1.0e-9
#define PGDELTA 1

//int BICUBIC_DIAG_=0;
//short hcc_opt_=0;
//short hcf_opt_=0;

//static  double (***cCC)[4][4];
//static  double (***cCF)[4][4];
//static  double (***cSiCl)[4][4];


/* ---------------------------------------------------------------------- */

PairTERSOFFHG::PairTERSOFFHG(LAMMPS *lmp) : PairTersoff(lmp)
{
  // See Eqs A7 and A8
  moliere_aB   = 0.5292; 
  moliere_c[0] = 0.35;
  moliere_c[1] = 0.55;
  moliere_c[2] = 0.10;
  moliere_d[0] = 0.30;
  moliere_d[1] = 1.20;
  moliere_d[2] = 6.00;

  firsov_const = pow(9.0*MY_PI*MY_PI/128.0, 1.0/3.0);
  NCl = NSi = NULL;
  nmax = 0;

  // set comm size needed by this Pair
  comm_forward = 1;
  comm_reverse = 1;

  maxlocal = 0;
  REBO_numneigh = NULL;
  REBO_firstneigh = NULL;
  ipage = NULL;
  pgsize = oneatom = 0;

  // read spline coefficients for coordination term, A12
  //bicubic_genCoef();
  //read_lib();
}

/* ---------------------------------------------------------------------- */

void PairTERSOFFHG::read_file(char *file)
{
  int params_per_line = 27;
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
      sprintf(str,"Cannot open TersoffHG potential file %s",file);
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
      error->all(FLERR,"Incorrect format in TERSOFFHG potential file");

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
    params[nparams].powerm = atof(words[3]);    // beta (A13)
    params[nparams].gamma = atof(words[4]);     // not used
    params[nparams].lam3 = atof(words[5]);      // alpha (A13)
    params[nparams].c = atof(words[6]);         // c (A14)
    params[nparams].d = atof(words[7]);         // d (A14)
    params[nparams].h = atof(words[8]);         // h (A14)
    params[nparams].powern = atof(words[9]);    // delta (A12)
    params[nparams].beta = atof(words[10]);     // not used
    params[nparams].lam2 = atof(words[11]);     // mu (A4)
    params[nparams].bigb = atof(words[12]);     // B (A4)
    params[nparams].bigr = atof(words[13]);     // rmin (A5)
    params[nparams].bigd = atof(words[14]);     // rmax - rmin (A5)
    params[nparams].lam1 = atof(words[15]);     // lambda (A3)
    params[nparams].biga = atof(words[16]);     // A (A3)
    params[nparams].powereta = atof(words[17]); // eta (A12)
    params[nparams].Z_i = atof(words[18]);      // Z_i (A7)
    params[nparams].Z_j = atof(words[19]);      // Z_j (A7)
    params[nparams].spl_ra = atof(words[20]);   // spl_ra (A10)
    params[nparams].spl_rb = atof(words[21]);   // spl_rb (A10)
    params[nparams].spl_a = atof(words[22]);    // spl_a (A10)
    params[nparams].spl_b = atof(words[23]);    // spl_b (A10)
    params[nparams].spl_c = atof(words[24]);    // spl_c (A10)
    params[nparams].spl_s = atof(words[25]);    // spl_s (A10)
    params[nparams].Re = atof(words[26]);       // Re (A13)

    // currently only allow m exponent of 1 or 3

    params[nparams].powermint = int(params[nparams].powerm);

    if (
        params[nparams].lam3 < 0.0 || params[nparams].c < 0.0 ||
        params[nparams].d < 0.0 || params[nparams].powern < 0.0 ||
        params[nparams].beta < 0.0 || params[nparams].lam2 < 0.0 ||
        params[nparams].bigb < 0.0 || params[nparams].bigr < 0.0 ||
        params[nparams].bigd < 0.0 ||
        params[nparams].bigd > params[nparams].bigr ||
        params[nparams].lam3 < 0.0 || params[nparams].biga < 0.0 ||
        params[nparams].powerm - params[nparams].powermint != 0.0 ||
        params[nparams].gamma < 0.0 ||
        params[nparams].Z_i < 1.0 || params[nparams].Z_j < 1.0 ||
        params[nparams].spl_ra < 0.0 || params[nparams].spl_rb < 0.0)
      error->all(FLERR,"Illegal REBO1 parameter");

    nparams++;
  }

  delete [] words;
}

/* ---------------------------------------------------------------------- */

void PairTERSOFFHG::read_lib()
{
  unsigned int maxlib = 1024;
  int i,j,k,l,nwords,m;
  int ii,jj,kk,ll,mm,iii;
  char s[maxlib];
  char **words = new char*[80];
  int nspl = 401;

  // open libraray file on proc 0

  FILE *fp;
  if (comm->me == 0) 
  {
    fp = force->open_potential("lib.rebo1");
    if (fp == NULL)
    {
      char str[128];
      sprintf(str,"Cannot open REBO1 lib.rebo1 file");
      error->one(FLERR,str);
    }

    for (i=0; i<4; i++) 
    {
      for (l=0; l<nspl; l++) 
      {
        fgets(s,maxlib,fp);
        nwords = 0;
        words[nwords++] = strtok(s," \t\n\r\f");
        while ((words[nwords++] = strtok(NULL," \t\n\r\f")))continue;
        coordnumber[i][l] = atof(words[0]);
        coordenergy[i][l] = atof(words[1]);
        coordforce[i][l]  = atof(words[2]);
      }
      for (l=0; l<nspl; l++) 
      {
        coordenergy[4][l] = 0.0;
        coordforce[4][l]  = 0.0;
      }
    }
  }

  MPI_Bcast(&coordnumber[0][0],1604,MPI_DOUBLE,0,world);
  MPI_Bcast(&coordenergy[0][0],1604,MPI_DOUBLE,0,world);
  MPI_Bcast(&coordforce[0][0],1604,MPI_DOUBLE,0,world);

  delete [] words;

}

/* ---------------------------------------------------------------------- */

void PairTERSOFFHG::count_neigh()
{
  int i, j, ii, jj, itype, jtype, n;
  int inum, jnum, param, *ilist, *jlist, *numneigh, **firstneigh;
  double r, rsq, delrij[3];
  const double cutshortsq = cutmax*cutmax;

  double **x = atom->x;
  int *type = atom->type;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->grow(NCl,nmax,"pair:NCl");
    memory->grow(NSi,nmax,"pair:NSi");
  }

  inum  = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = map[type[i]];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // skip immediately if center atom is not Si
    if (strcmp(elements[map[itype+1]],"Si") != 0) continue;

    NCl[i] = NSi[i] = 0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj] & NEIGHMASK;
      jtype = map[type[j]];


      delrij[0] = x[i][0] - x[j][0];
      delrij[1] = x[i][1] - x[j][1];
      delrij[2] = x[i][2] - x[j][2];
      rsq = vec3_dot(delrij,delrij);
      param = elem2param[itype][jtype][jtype];

      if (rsq > cutshortsq) continue;

      r = sqrt(rsq);
      if (strcmp(elements[map[jtype+1]],"Cl") == 0) NCl[i] += ters_fc(r,&params[param]);
      if (strcmp(elements[map[jtype+1]],"Si") == 0) NSi[i] += ters_fc(r,&params[param]);
    }
  }

  // communicating coordination number to all nodes
  pack_flag = 1;
  comm->forward_comm_pair(this);
  pack_flag = 2;
  comm->forward_comm_pair(this);

}

/* ---------------------------------------------------------------------- */

// void PairTERSOFFHG::compute(int eflag, int vflag)
// {
//   int i,j,k,ii,jj,kk,inum,jnum;
//   int itype,jtype,ktype,iparam_ij,iparam_ijk,iparam_ik;
//   tagint itag,jtag;
//   double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
//   double rsq,rsq1,rsq2;
//   double delr1[3],delr2[3],fi[3],fj[3],fk[3];
//   double zeta_ij,prefactor;
//   int nncl, nnsi;
//   int *ilist,*jlist,*numneigh,**firstneigh;

//   evdwl = 0.0;
//   if (eflag || vflag) ev_setup(eflag,vflag);
//   else evflag = vflag_fdotr = vflag_atom = 0;

//   double **x = atom->x;
//   double **f = atom->f;
//   tagint *tag = atom->tag;
//   int *type = atom->type;
//   int nlocal = atom->nlocal;
//   int newton_pair = force->newton_pair;
//   const double cutshortsq = cutmax*cutmax;

//   inum = list->inum;
//   ilist = list->ilist;
//   numneigh = list->numneigh;
//   firstneigh = list->firstneigh;

//   // count number of nearest neighbors; needs communications
//   count_neigh();

//   double fxtmp,fytmp,fztmp;

//   // loop over full neighbor list of my atoms

//   for (ii = 0; ii < inum; ii++) 
//   {
//     i = ilist[ii];
//     itag = tag[i];
//     itype = map[type[i]];
//     xtmp = x[i][0];
//     ytmp = x[i][1];
//     ztmp = x[i][2];
//     fxtmp = fytmp = fztmp = 0.0;

//     // two-body interactions, skip half of them

//     jlist = firstneigh[i];
//     jnum = numneigh[i];
//     int numshort = 0; 

//     nncl = int((NCl[i]-1.0)*100.0);
//     nnsi = int(NSi[i]);

//     for (jj = 0; jj < jnum; jj++) 
//     {
//       j = jlist[jj];
//       j &= NEIGHMASK;

//       delx = xtmp - x[j][0];
//       dely = ytmp - x[j][1];
//       delz = ztmp - x[j][2];
//       rsq = delx*delx + dely*dely + delz*delz;

//       if (rsq < cutshortsq) 
//       {
//         neighshort[numshort++] = j;
//         if (numshort >= maxshort) 
//         {
//           maxshort += maxshort/2;
//           memory->grow(neighshort,maxshort,"pair:neighshort");
//         }
//       }

//       jtag = tag[j];
//       if (itag > jtag) 
//       {
//         if ((itag+jtag) % 2 == 0) continue;
//       } 
//       else if (itag < jtag) 
//       {
//         if ((itag+jtag) % 2 == 1) continue;
//       } 
//       else 
//       {
//         if (x[j][2] < x[i][2]) continue;
//         if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
//         if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
//       }

//       jtype = map[type[j]];
//       iparam_ij = elem2param[itype][jtype][jtype];
//       if (rsq >= params[iparam_ij].cutsq) continue;

//       repulsive(&params[iparam_ij],rsq,fpair,eflag,evdwl);
//       //if (tag[j] > tag[i])
//       //  std::cerr<<tag[i]<<" "<<tag[j]<<" "<<evdwl<<std::endl;

//       fxtmp += delx*fpair;
//       fytmp += dely*fpair;
//       fztmp += delz*fpair;
//       f[j][0] -= delx*fpair;
//       f[j][1] -= dely*fpair;
//       f[j][2] -= delz*fpair;

//       if (evflag) ev_tally(i,j,nlocal,newton_pair,
//                            evdwl,0.0,fpair,delx,dely,delz);


//     //std::cerr<<tag[i]<<" "<<tag[j]<<" "<<fpair<<std::endl;
//     //" "<<delx*fpair<<" "<<dely*fpair<<" "<<dely*fpair<<" "
//     //<<forces<<std::endl;
//     }

//     // three-body interactions
//     // skip immediately if I-J is not within cutoff
//     double fjxtmp,fjytmp,fjztmp;

//     for (jj = 0; jj < numshort; jj++) 
//     {
//       j = neighshort[jj];
//       jtype = map[type[j]];
//       iparam_ij = elem2param[itype][jtype][jtype];
  
//       delr1[0] = x[j][0] - xtmp;
//       delr1[1] = x[j][1] - ytmp;
//       delr1[2] = x[j][2] - ztmp;
//       rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
//       if (rsq1 >= params[iparam_ij].cutsq) continue;

//       // accumulate bondorder zeta for each i-j interaction via loop over k

//       fjxtmp = fjytmp = fjztmp = 0.0;
//       zeta_ij = 0.0;

//       for (kk = 0; kk < numshort; kk++) 
//       {
//         if (jj == kk) continue;
//         k = neighshort[kk];
//         ktype = map[type[k]];
//         iparam_ijk = elem2param[itype][jtype][ktype];
//       	iparam_ik  = elem2param[itype][ktype][ktype];

//         delr2[0] = x[k][0] - xtmp;
//         delr2[1] = x[k][1] - ytmp;
//         delr2[2] = x[k][2] - ztmp;
//         rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
//         if (rsq2 >= params[iparam_ijk].cutsq) continue;

//         zeta_ij += zeta(&params[iparam_ijk],&params[iparam_ik],rsq1,rsq2,delr1,delr2);
//       }

//       // pairwise force due to zeta

//       force_zeta(&params[iparam_ij],rsq1,zeta_ij,fpair,prefactor,eflag,evdwl,nncl,nnsi);

//       fxtmp += delr1[0]*fpair;
//       fytmp += delr1[1]*fpair;
//       fztmp += delr1[2]*fpair;
//       fjxtmp -= delr1[0]*fpair;
//       fjytmp -= delr1[1]*fpair;
//       fjztmp -= delr1[2]*fpair;

//       if (evflag) ev_tally(i,j,nlocal,newton_pair,
//                            evdwl,0.0,-fpair,-delr1[0],-delr1[1],-delr1[2]);

//       // attractive term via loop over k

//       for (kk = 0; kk < numshort; kk++) {
//         if (jj == kk) continue;
//         k = neighshort[kk];
//         ktype = map[type[k]];
//         iparam_ijk = elem2param[itype][jtype][ktype];
//       	iparam_ik  = elem2param[itype][ktype][ktype];

//         delr2[0] = x[k][0] - xtmp;
//         delr2[1] = x[k][1] - ytmp;
//         delr2[2] = x[k][2] - ztmp;
//         rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
//         if (rsq2 >= params[iparam_ijk].cutsq) continue;

//         attractive(&params[iparam_ijk],&params[iparam_ik],prefactor,
//                    rsq1,rsq2,delr1,delr2,fi,fj,fk);

//         fxtmp += fi[0];
//         fytmp += fi[1];
//         fztmp += fi[2];
//         fjxtmp += fj[0];
//         fjytmp += fj[1];
//         fjztmp += fj[2];
//         f[k][0] += fk[0];
//         f[k][1] += fk[1];
//         f[k][2] += fk[2];

//         if (vflag_atom) v_tally3(i,j,k,fj,fk,delr1,delr2);

//       }
//       f[j][0] += fjxtmp;
//       f[j][1] += fjytmp;
//       f[j][2] += fjztmp;
//     }
//     f[i][0] += fxtmp;
//     f[i][1] += fytmp;
//     f[i][2] += fztmp;

//   }

//   if (vflag_fdotr) virial_fdotr_compute();

  
//   // for (ii = 0; ii < inum; ii++) 
//   // {
//   //   i = ilist[ii];
//   //   std::cerr<<tag[i]<<" "<<f[i][0]<<" "<<f[i][1]<<" "<<f[i][2]<<std::endl;
//   // }
// }

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
void PairTERSOFFHG::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum;
  int itype,jtype,ktype,iparam_ij,iparam_ijk,iparam_ik;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fi[3],fj[3],fk[3];
  double zeta_ij,prefactor;
  int nncl, nnsi;
  int *ilist,*jlist,*numneigh,**firstneigh;
  HGvector Rij, Rij1;
  double rij;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = vflag_atom = 0;

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  const double cutshortsq = cutmax*cutmax;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // count number of nearest neighbors; needs communications
  count_neigh();

  int n,allnum;
  double dS;
  int *neighptr;
  

  if (atom->nmax > maxlocal) {
    maxlocal = atom->nmax;
    memory->destroy(REBO_numneigh);
    memory->sfree(REBO_firstneigh);
//    memory->destroy(nC);
//    memory->destroy(nH);
    memory->create(REBO_numneigh,maxlocal,"AIREBO:numneigh");
    REBO_firstneigh = (int **) memory->smalloc(maxlocal*sizeof(int *),
                                               "AIREBO:firstneigh");
//    memory->create(nC,maxlocal,"AIREBO:nC");
//    memory->create(nH,maxlocal,"AIREBO:nH");
  }

  allnum = list->inum + list->gnum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // store all REBO neighs of owned and ghost atoms
  // scan full neighbor list of I

  ipage->reset();
  

  for (ii = 0; ii < allnum; ii++) {
    i = ilist[ii];

    n = 0;
    neighptr = ipage->vget();

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = map[type[i]];
    //nC[i] = nH[i] = 0.0;
    jlist = firstneigh[i];
    jnum = numneigh[i];

    double rcmaxsq = 9.0;
    double rcmin, rcmax;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < rcmaxsq){//[itype][jtype]) {
        rcmin = params[iparam_ij].bigr;
        rcmax = rcmin + params[iparam_ij].bigd;
        neighptr[n++] = j;
//        if (jtype == 0)
//          nC[i] += Sp(sqrt(rsq),rcmin,rcmax,dS);
          //nC[i] += Sp(sqrt(rsq),rcmin[itype][jtype],rcmax[itype][jtype],dS);
//        else
//          nH[i] += Sp(sqrt(rsq),rcmin,rcmax,dS);
          //nH[i] += Sp(sqrt(rsq),rcmin[itype][jtype],rcmax[itype][jtype],dS);
      }
    }

    REBO_firstneigh[i] = neighptr;
    REBO_numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  for (ii = 0; ii < inum; ii++) 
  {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    jlist = REBO_firstneigh[i];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

//     //nncl = int((NCl[i]-1.0)*100.0);
//     //nnsi = int(NSi[i]);
    for (jj = 0; jj < REBO_numneigh[i]; jj++) 
    {
      j = jlist[jj];
//      j &= NEIGHMASK;
      jtag = tag[j];
//        if (tag[j]==9)
//          std::cerr<<"A "<<f[j][0]<<" "<<f[j][1]<<" "<<f[j][2]<<"\n";

      if (itag > jtag) 
      {
        if ((itag+jtag) % 2 == 0) continue;
      } 
      else if (itag < jtag) 
      {
        if ((itag+jtag) % 2 == 1) continue;
      } 
      else 
      {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      Rij.set(x[i][0]-x[j][0], x[i][1]-x[j][1], x[i][2]-x[j][2]);
      rij=Rij.mag();
      rsq = Rij.sqmag();
      Rij1 = Rij/rij;
      jtype = map[type[j]];
      iparam_ij = elem2param[itype][jtype][jtype];

      //if (rsq <= params[iparam_ij].cutsq && 
      //if (jtag > itag)
      //{ 

        double ij_f, ij_fprime;
        double B = params[iparam_ij].bigb;
        double MORSE_MU = params[iparam_ij].lam2;
        double MORSE_LAM = params[iparam_ij].lam1;
        double bigr= params[iparam_ij].bigr;
        double bigd = params[iparam_ij].bigd;
        ij_f = Sp(rij, bigr, bigr+bigd, ij_fprime);
        double aa = B * SpExp(-MORSE_MU * rij);
        VA_ij = ij_f * aa;
        dVA_ij = aa * (ij_fprime - ij_f*MORSE_MU);
        if (rij < params[iparam_ij].spl_ra)
        {
          double moa = 0.388874*pow(params[iparam_ij].Z_i+params[iparam_ij].Z_j,-2./3.);
          double moeps=14.39965*params[iparam_ij].Z_i*params[iparam_ij].Z_j;
          double rija=rij/moa;
          VR_ij = moeps/rij*(( 0.35*SpExp(-0.3*rija)
                            +0.55*SpExp(-1.2*rija)
                            +0.1*SpExp(-6*rija))) + params[iparam_ij].spl_s;
          dVR_ij = -moeps/rij*(( 0.35*SpExp(-0.3*rija)*(1/rij+0.3/moa)
                              +0.55*SpExp(-1.2*rija)*(1/rij+1.2/moa)
                              +0.1*SpExp(-6*rija)*(1/rij+6/moa)));
        }
        else if (rij < params[iparam_ij].spl_rb)
        {
          VR_ij = params[iparam_ij].spl_c + SpExp(params[iparam_ij].spl_a * rij + params[iparam_ij].spl_b);
          dVR_ij = params[iparam_ij].spl_a * SpExp(params[iparam_ij].spl_a * rij + params[iparam_ij].spl_b);
        }
        else
        {
          aa=params[iparam_ij].biga * SpExp(-MORSE_LAM * rij);
          VR_ij = ij_f * aa;
          dVR_ij = aa * (ij_fprime - ij_f*MORSE_LAM);
        }
        double b_ij=BondOrder(i, j, VA_ij, eflag, vflag);

        evdwl = VR_ij - b_ij*VA_ij;
        fpair = (-dVR_ij+b_ij*dVA_ij)/rij;
        HGvector Fij=(-dVR_ij+b_ij*dVA_ij)*Rij1;

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                          evdwl, 0.0, fpair, Rij.x, Rij.y, Rij.z);

//        if (tag[i]==1)
//          std::cerr<<f[i][0]<<" "<<f[i][1]<<" "<<f[i][2]<<" ";
//        if (tag[j]==9)
//          std::cerr<<"B "<<f[j][0]<<" "<<f[j][1]<<" "<<f[j][2]<<"\n";

        f[i][0] += Fij.x;
        f[i][1] += Fij.y;
        f[i][2] += Fij.z;

        f[j][0] -= Fij.x;
        f[j][1] -= Fij.y;
        f[j][2] -= Fij.z;
//        if (tag[i]==1)
//          std::cerr<<f[i][0]<<" "<<f[i][1]<<" "<<f[i][2]<<"\n";
//        if (tag[j]==9)
//          std::cerr<<"C "<<f[j][0]<<" "<<f[j][1]<<" "<<f[j][2]<<"\n";

//        std::cerr<<itag<<" "<<jtag<<" "<<Fij.x<<" "<<Fij.y<<" "<<Fij.z<<" "<<std::endl;

//        std::cerr<<itag<<" "<<jtag<<" "<<" "<<b_ij<<" "<<" "<<std::endl;
      //}
    }
  }
//  for (ii = 0; ii < allnum; ii++) {
//    i = ilist[ii];
//    if (tag[i] < 10)
//      std::cerr<<tag[i]<<" "<<f[i][0]<<" "<<f[i][1]<<" "<<f[i][2]<<" "<<std::endl;
//  }

}

//******************************************************
void PairTERSOFFHG::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style AIREBO requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style AIREBO requires newton pair on");

  // need a full neighbor list, including neighbors of ghosts

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->ghost = 1;

  // local REBO neighbor list
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

//******************************************************
double PairTERSOFFHG::BondOrder(int i, int j, double Pre, int eflag, int vflag){
  HGvector Rij, Rik, Rjk, Rij1, Rik1, Rjk1, Fij;
  double rij, rik, rjk;
  double ik_f, ij_f, jk_f, ik_fprime, ij_fprime, jk_fprime;
  iFv_ij=0; iFv_ik.clear(); iFv_jk.clear(); iFv_ik2.clear();
  jFv_ij=0; jFv_ik.clear(); jFv_jk.clear(); jFv_jk2.clear();
  double cos_theta;
  scr_Rhat.clear();
  iNconj=jNconj=0;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *klist;
  int k,ii,jj,kk,inum,jnum,knum;
  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;  
  tagint itag, jtag, ktag;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  int itype,jtype,ktype,iparam_ij,iparam_ji,iparam_ijk,iparam_ik,iparam_jk;
  itype = map[type[i]];
  jtype = map[type[j]];
  iparam_ij = elem2param[itype][jtype][jtype];
  iparam_ji = elem2param[jtype][itype][itype];
    
  Rij.set(x[i][0]-x[j][0], x[i][1]-x[j][1], x[i][2]-x[j][2]);
  rij=Rij.mag();
  Rij1=Rij/rij;
 
  // double NF_ij, NC_ij, NSi_ij, Nt_ij, NF_ji, NC_ji, NSi_ji, Nt_ji;
  // double NCl_ij, NCl_ji;
  // double P_ij, P_ji, dFP_ij, dFP_ji, dCP_ij, dCP_ji, dlam;
  double el, elf, alpha, beta, Re, n, dlam, g, g1;//, xik, Yik, Y1ik, F_ij, dFi_ij, dFj_ij, dFc_ij;
  double eta_i, eta_j, delta_i, delta_j;

  // NF_ij  = i->Nmap[9]  - (j->id==9)  * bond_ij->f;
  // NC_ij  = i->Nmap[6]  - (j->id==6)  * bond_ij->f;
  // NSi_ij = i->Nmap[14] - (j->id==14) * bond_ij->f;
  // NCl_ij = i->Nmap[17] - (j->id==17) * bond_ij->f;
  // NF_ji  = j->Nmap[9]  - (i->id==9)  * bond_ij->f;
  // NC_ji  = j->Nmap[6]  - (i->id==6)  * bond_ij->f;
  // NSi_ji = j->Nmap[14] - (i->id==14) * bond_ij->f;
  // NCl_ji = j->Nmap[17] - (i->id==17) * bond_ij->f;

  // Nt_ij = NF_ij + NC_ij + NSi_ij + NCl_ij;
  // Nt_ji = NF_ji + NC_ji + NSi_ji + NCl_ji;

  // P_ij=P_ji=dFP_ij=dFP_ji=dCP_ij=dCP_ji=0;
  
  // if (bond_ij->type==12)
  // {
  //   Pcc_bicubicint(NF_ij, NC_ij+NSi_ij, &P_ij, &dFP_ij,&dCP_ij);
  //   Pcc_bicubicint(NF_ji, NC_ji+NSi_ji, &P_ji, &dFP_ji,&dCP_ji);
  // }
  // else if (bond_ij->type==15)
  // {
  //   if (i->id==6) Pcf_bicubicint(NF_ij, NC_ij+NSi_ij, &P_ij, &dFP_ij,&dCP_ij);
  //   else Pcf_bicubicint(NF_ji, NC_ji+NSi_ji, &P_ji, &dFP_ji,&dCP_ji);
  // }
  // else if (bond_ij->type==23)
  // {
  //   if (i->id==14) Psif_bicubicint(NF_ij, NC_ij+NSi_ij, &P_ij, &dFP_ij,&dCP_ij);
  //   else Psif_bicubicint(NF_ji,  NC_ji+NSi_ji, &P_ji, &dFP_ji,&dCP_ji);
  // }
  // else if (bond_ij->type==31)
  // {
  //   if (i->id==14) Psicl_bicubicint(NCl_ij, NC_ij+NSi_ij, &P_ij, &dFP_ij,&dCP_ij);
  //   else Psicl_bicubicint(NCl_ji,  NC_ji+NSi_ji, &P_ji, &dFP_ji,&dCP_ji);
  // }
  // bsp_ij=P_ij; bsp_ji=P_ji;
  bsp_ij=bsp_ji=0;

  // //***************for k neighbor of i*************************
  klist = REBO_firstneigh[i];
  for (kk = 0; kk < REBO_numneigh[i]; kk++) 
  {
    k = klist[kk];
    k &= NEIGHMASK;
    ktag = tag[k];
    ktype = map[type[k]];
    iparam_ik = elem2param[itype][ktype][ktype];
    iparam_ijk = elem2param[itype][jtype][ktype];
    Rik.set(x[i][0]-x[k][0], x[i][1]-x[k][1], x[i][2]-x[k][2]);
    rik=Rik.mag();
    Rik1=Rik/rik;
    if (tag[k]!=tag[j] && Rik.sqmag() <= params[iparam_ik].cutsq)
    {
      Rjk.set(x[j][0]-x[k][0], x[j][1]-x[k][1], x[j][2]-x[k][2]);
      rjk=Rjk.mag();
      Rjk1=Rjk/rjk;
      scr_Rhat.push_back(Rjk1);
      double bigr= params[iparam_ik].bigr;
      double bigd = params[iparam_ik].bigd;
      ik_f = Sp(rik, bigr, bigr+bigd, ik_fprime);
      // //*******b-sigma-pi calculations************
      alpha = params[iparam_ijk].lam3;
      if (alpha)
      {
        beta = params[iparam_ijk].powerm;
        //Re = params[iparam_ijk].Re;
        n = (rij - params[iparam_ij].Re) - (rik - params[iparam_ik].Re);
        el=SpExp(alpha*pow(n, beta));
        dlam=alpha;
        if (beta!=1){
          dlam*=beta*pow(n, beta-1);
        }
      }
      else{
        el = 1.0; dlam=0;
      }
      elf=el*ik_f;
      cos_theta = Rij1 ^ Rik1;
      g = params[iparam_ijk].c + params[iparam_ijk].d*pow(params[iparam_ijk].h-cos_theta,2);
      g1 = -2*params[iparam_ijk].d*(params[iparam_ijk].h-cos_theta);
      
      iFv_ij+=(elf*(g1*(1/rik - cos_theta/rij) + g*dlam));
      iFv_ik.push_back(ik_fprime*
           (g*el + 0)//((k->id==9||k->id==17) ? dFP_ij : dCP_ij))
           +elf*(g1*(1/rij - cos_theta/rik) - g*dlam));
      iFv_jk.push_back(elf*(-g1*rjk/rik/rij));
            
      bsp_ij += g*elf;
    }
  }

  // //***************for k neighbor of j*************************
  klist = REBO_firstneigh[j];
  for (kk = 0; kk < REBO_numneigh[j]; kk++) 
  {
    k = klist[kk];
    k &= NEIGHMASK;
    ktag = tag[k];
    ktype = map[type[k]];
    iparam_jk = elem2param[jtype][ktype][ktype];
    iparam_ijk = elem2param[jtype][itype][ktype];
    Rjk.set(x[j][0]-x[k][0], x[j][1]-x[k][1], x[j][2]-x[k][2]);
    rjk=Rjk.mag();
    Rjk1=Rjk/rjk;

    if (tag[k]!=tag[i] && Rjk.sqmag() <= params[iparam_jk].cutsq)
    {
      Rik.set(x[i][0]-x[k][0], x[i][1]-x[k][1], x[i][2]-x[k][2]);
      rik=Rik.mag();
      Rik1=Rik/rik;
      scr_Rhat.push_back(Rik1);
      double bigr= params[iparam_jk].bigr;
      double bigd = params[iparam_jk].bigd;
      jk_f = Sp(rjk, bigr, bigr+bigd, jk_fprime);
      // //*******b-sigma-pi calculations************
      alpha = params[iparam_ijk].lam3;
      if (alpha)
      {
        beta = params[iparam_ijk].powerm;
        //Re = params[iparam_ijk].Re;
        n = (rij - params[iparam_ji].Re) - (rjk - params[iparam_jk].Re);
        el=SpExp(alpha*pow(n, beta));
        dlam=alpha;
        if (beta!=1){
          dlam*=beta*pow(n, beta-1);
        }
      }
      else{
        el = 1.0; dlam=0;
      }
      elf=el*jk_f;
      cos_theta = - (Rij1 ^ Rjk1);
      g = params[iparam_ijk].c + params[iparam_ijk].d*pow(params[iparam_ijk].h-cos_theta,2);
      g1 = -2*params[iparam_ijk].d*(params[iparam_ijk].h-cos_theta);
      
      jFv_ij+=(elf*(g1*(1/rjk - cos_theta/rij) + g*dlam));
      jFv_jk.push_back(jk_fprime*
           (g*el + 0)//((k->id==9||k->id==17) ? dFP_ij : dCP_ij))
           +elf*(g1*(1/rij - cos_theta/rjk) - g*dlam));
      jFv_ik.push_back(elf*(-g1*rik/rjk/rij));
            
      bsp_ji += g*elf;
    }
  }
  //for (std::vector<double>::iterator v=jFv_jk.begin(); v!=jFv_jk.end(); v++)
  //  std::cerr<<*v<<" ";
  //std::cerr<<std::endl;

  //***********************************************

  Pre*=0.5;
  eta_i = params[iparam_ij].powereta;
  eta_j = params[iparam_ji].powereta;
  delta_i = params[iparam_ij].powern;
  delta_j = params[iparam_ji].powern;

  double bbar_ij = 0.5*(pow(1 + pow(bsp_ij, eta_i), -delta_i)
     + pow(1 + pow(bsp_ji, eta_j), -delta_j)); // +F_ij);

  if (eta_i!=1)
  {
    if (bsp_ij <= 0) bsp_ij=0; //unphysical, but happens due to roundoff
    else
      bsp_ij=pow(1 + pow(bsp_ij, eta_i), -delta_i-1)*
                 eta_i*pow(bsp_ij, eta_i-1);
  }
  else
    bsp_ij=pow(1 + bsp_ij, -delta_i-1);
  bsp_ij*=-Pre*delta_i;
  
  if (eta_j!=1)
  {
    if (bsp_ji <= 0) bsp_ji=0;
    else
      bsp_ji=pow(1 + pow(bsp_ji, eta_j), -delta_j-1)*
                 eta_j*pow(bsp_ji, eta_j-1);
  }
  else  
    bsp_ji=pow(1 + bsp_ji, -delta_j-1);
  
  bsp_ji*=-Pre*delta_j;
  Fij=(bsp_ij*iFv_ij+bsp_ji*jFv_ij)*Rij1;
  
  if (evflag) ev_tally(i,j,nlocal,newton_pair,
        0.0,0.0,(bsp_ij*iFv_ij+bsp_ji*jFv_ij)/rij,Rij.x,Rij.y,Rij.z);

  f[i][0] += Fij.x;
  f[i][1] += Fij.y;
  f[i][2] += Fij.z;

  f[j][0] -= Fij.x;
  f[j][1] -= Fij.y;
  f[j][2] -= Fij.z;

  cnt=0;
  scrcount=0;
  // //***************for k neighbor of i*************************
  klist = REBO_firstneigh[i];
  for (kk = 0; kk < REBO_numneigh[i]; kk++) 
  {
    k = klist[kk];
    k &= NEIGHMASK;
    ktag = tag[k];
    ktype = map[type[k]];
    iparam_ik = elem2param[itype][ktype][ktype];
    iparam_ijk = elem2param[itype][jtype][ktype];
    Rik.set(x[i][0]-x[k][0], x[i][1]-x[k][1], x[i][2]-x[k][2]);
    rik=Rik.mag();
    Rik1=Rik/rik;
    if (tag[k]!=tag[j] && Rik.sqmag() <= params[iparam_ik].cutsq)
    {
      Rjk.set(x[j][0]-x[k][0], x[j][1]-x[k][1], x[j][2]-x[k][2]);
      rjk=Rjk.mag();
      Rjk1=Rjk/rjk;
      force_ik=bsp_ij*iFv_ik[cnt];
      double bigr= params[iparam_ik].bigr;
      double bigd = params[iparam_ik].bigd;
      ik_f = Sp(rik, bigr, bigr+bigd, ik_fprime);
      Fij = force_ik*Rik1;
      if (evflag) ev_tally(i,k,nlocal,newton_pair,
            0.0,0.0,force_ik/rik,Rik.x,Rik.y,Rik.z);

      f[i][0] += Fij.x;
      f[i][1] += Fij.y;
      f[i][2] += Fij.z;

      f[k][0] -= Fij.x;
      f[k][1] -= Fij.y;
      f[k][2] -= Fij.z;
      
      force_jk = bsp_ij * iFv_jk[cnt];
      Fij =  force_jk * Rjk1; //scr_Rhat[scrcount++];

      if (evflag) ev_tally(j,k,nlocal,newton_pair,
            0.0,0.0,force_jk/rjk,Rjk.x,Rjk.y,Rjk.z);

      f[j][0] += Fij.x;
      f[j][1] += Fij.y;
      f[j][2] += Fij.z;

      f[k][0] -= Fij.x;
      f[k][1] -= Fij.y;
      f[k][2] -= Fij.z;

      cnt++;
    }
  }
  // //***************for k neighbor of j*************************
  cnt=0;
  klist = REBO_firstneigh[j];
  for (kk = 0; kk < REBO_numneigh[j]; kk++) 
  {
    k = klist[kk];
    k &= NEIGHMASK;
    ktag = tag[k];
    ktype = map[type[k]];
    iparam_jk = elem2param[jtype][ktype][ktype];
    iparam_ijk = elem2param[jtype][itype][ktype];
    Rjk.set(x[j][0]-x[k][0], x[j][1]-x[k][1], x[j][2]-x[k][2]);
    rjk=Rjk.mag();
    Rjk1=Rjk/rjk;
    if (tag[k]!=tag[i] && Rjk.sqmag() <= params[iparam_jk].cutsq)
    {
      Rik.set(x[i][0]-x[k][0], x[i][1]-x[k][1], x[i][2]-x[k][2]);
      rik=Rik.mag();
      Rik1=Rik/rik;
      force_jk=bsp_ji*jFv_jk[cnt];
      double bigr= params[iparam_jk].bigr;
      double bigd = params[iparam_jk].bigd;
      //jk_f = Sp(rjk, bigr, bigr+bigd, jk_fprime);
      Fij = force_jk*Rjk1;
//      std::cerr<<jFv_jk[cnt]<<std::endl;


      if (evflag) ev_tally(j,k,nlocal,newton_pair,
            0.0,0.0,force_jk/rjk,Rjk.x,Rjk.y,Rjk.z);

      f[j][0] += Fij.x;
      f[j][1] += Fij.y;
      f[j][2] += Fij.z;

      f[k][0] -= Fij.x;
      f[k][1] -= Fij.y;
      f[k][2] -= Fij.z;
      
      force_ik = bsp_ji * jFv_ik[cnt];
      Fij =  force_ik * Rik1; //scr_Rhat[scrcount++];
      //std::cerr<<Fij.x<<" "<<Fij.y<<" "<<Fij.z<<std::endl;

      if (evflag) ev_tally(i,k,nlocal,newton_pair,
            0.0,0.0,force_ik/rik,Rik.x,Rik.y,Rik.z);

      f[i][0] += Fij.x;
      f[i][1] += Fij.y;
      f[i][2] += Fij.z;

      f[k][0] -= Fij.x;
      f[k][1] -= Fij.y;
      f[k][2] -= Fij.z;
      cnt++;
    }
  }  
  return bbar_ij;
}
/* ---------------------------------------------------------------------- */

void PairTERSOFFHG::repulsive(Param *param, double rsq, double &fforce,
                               int eflag, double &eng)
{
  double r,rinv,tmp_fc,tmp_fc_d,tmp_exp;
  double energy, forces;
  r = sqrt(rsq);
  rinv = 1.0/r;

  double spl_ra = param->spl_ra;
  double spl_rb = param->spl_rb;

  if (r > spl_rb) {
  // Tersoff repulsive (of Morse form): when r > rb; Eq. A3
 
    tmp_fc = ters_fc(r,param);
    tmp_fc_d = ters_fc_d(r,param);
    tmp_exp = exp(-param->lam1 * r);
    forces = param->biga * tmp_exp * (tmp_fc_d - tmp_fc*param->lam1);
    energy = tmp_fc * param->biga * tmp_exp;

  } else if (r < spl_rb && r > spl_ra) {
  // Spline repulsive: when ra < r < rb; Eq. A9

    double spl_a = param->spl_a;
    double spl_b = param->spl_b;
    double spl_c = param->spl_c;
    energy = spl_c + exp(spl_a * r + spl_b);
    forces = spl_a * exp(spl_a * r + spl_b);

  } else if (r < spl_ra) {
  // Moliere repulsive: when r < ra; Eqa. A7-8

    double zi = param->Z_i;
    double zj = param->Z_j;
    double spl_s = param->spl_s;

    // Eq. A7, first term
    double zizj = zi * zj * force->qqr2e * force->qelectron * force->qelectron * rinv;

    // Eq. A8
    double zizj_half = pow(zi,0.5) + pow(zj,0.5);
    double screen_a = 0.83 * firsov_const * moliere_aB * pow(zizj_half,-2.0/3.0);

    // Eq. A7, second term
    double cexpdra = 0.0;
    for (int i = 0; i < 3; i++)
      cexpdra += moliere_c[i] * exp(-1.0 * moliere_d[i] * r / screen_a);

    energy = zizj * cexpdra + spl_s;

    // forces: derivative of first term, Eq. A7
    double r2inv = 1.0/rsq;
    double dzizj = -1.0 * zizj * r2inv;
  
    // forces: derivative of second term, Eq. A7
    double dcexpdra = 0.0;
    for (int i = 0; i < 3; i++)
      dcexpdra += -moliere_c[i] * moliere_d[i] / screen_a * exp(-moliere_d[i] * r / screen_a);

    forces = dzizj * cexpdra + zizj * dcexpdra;
  }
  fforce = -forces / r;
  if (eflag) eng = energy;
}

/* ---------------------------------------------------------------------- */

double PairTERSOFFHG::ters_fa(double r, Param *param)
{
  // Tersoff attraction term; Eq. A4

  if (r > param->bigr + param->bigd) return 0.0;
  return -param->bigb * exp(-param->lam2 * r) * ters_fc(r,param);
}

/* ---------------------------------------------------------------------- */

double PairTERSOFFHG::ters_fa_d(double r, Param *param)
{
  // derivative of Tersoff attraction; Eq. A4
  
  if (r > param->bigr + param->bigd) return 0.0;
  return param->bigb * exp(-param->lam2 * r) *
    (param->lam2 * ters_fc(r,param) - ters_fc_d(r,param));
}

/* ---------------------------------------------------------------------- */

double PairTERSOFFHG::ters_fc(double r, Param *param)
{
  // Teroff cutoff function; Eq. A6
  
  double ters_R = param->bigr;	// Rmin
  double ters_D = param->bigd;	

  if (r < ters_R) return 1.0;
  if (r > ters_R+ters_D) return 0.0;
  double t = (r - ters_R)/ters_D;
  double t2 = t*t;
  double t3 = t2*t;
  return 1.0 - t3*(6.0*t2-15.0*t+10.0);
}

/* ---------------------------------------------------------------------- */

double PairTERSOFFHG::ters_fc_d(double r, Param *param)
{
  // derivative of Tersoff cutoff function; Eq. A6
  
  double ters_R = param->bigr;
  double ters_D = param->bigd;

  if (r < ters_R) return 0.0;
  if (r > ters_R+ters_D) return 0.0;
  double t = (r - ters_R);
  double t2 = t*t;
  return -30.0*t2*(t-ters_D)*(t-ters_D)/pow(ters_D,5);

}

/* ---------------------------------------------------------------------- */

double PairTERSOFFHG::zeta(Param *paramij, Param *paramik, double rsqij, double rsqik,
                         double *delrij, double *delrik)
{
  // zeta term inside the bond order term; Eq. A13
  
  double rij,rik,cos_theta,tmp,ex_delr;
  double Re_ij = paramij->Re;
  double Re_ik = paramik->Re;

  rij = sqrt(rsqij);
  rik = sqrt(rsqik);
  cos_theta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] +
              delrij[2]*delrik[2]) / (rij*rik);

  tmp = paramij->lam3 * pow(((rij-Re_ij)-(rik-Re_ik)),paramij->powermint);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  return ters_fc(rik,paramij) * ters_gijk(cos_theta,paramij) * ex_delr;
}

/* ---------------------------------------------------------------------- */

void PairTERSOFFHG::force_zeta(Param *param, double rsq, double zeta_ij,
                             double &fforce, double &prefactor,
                             int eflag, double &eng, int nncl, int nnsi)
{
  double r,fa,fa_d,bij;

  r = sqrt(rsq);
  fa = ters_fa(r,param);
  fa_d = ters_fa_d(r,param);
  bij = ters_bij(zeta_ij,param,nncl,nnsi);
  //std::cerr<<bij<<std::endl;
  fforce = 0.5*bij*fa_d / r;
  prefactor = -0.5*fa * ters_bij_d(zeta_ij,param,nncl,nnsi);
  if (eflag) eng = 0.5*bij*fa;
}

/* ---------------------------------------------------------------------- */

void PairTERSOFFHG::ters_zetaterm_d(Param *paramij, Param *paramik, double prefactor,
                                  double *rij_hat, double rij,
                                  double *rik_hat, double rik,
                                  double *dri, double *drj, double *drk)
{
  // derivate of zeta term inside the bond order term; Eq. A13
  
  double gijk,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];
  double gijkdri[3],gijkdrj[3],gijkdrk[3];
  double Re_ij = paramij->Re;
  double Re_ik = paramik->Re;

  fc = ters_fc(rik,paramij);
  dfc = ters_fc_d(rik,paramij);
  tmp = paramij->lam3 * pow(((rij-Re_ij)-(rik-Re_ik)),paramij->powermint);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  ex_delr_d = paramij->lam3 * paramij->powermint * pow(((rij-Re_ij)-(rik-Re_ik)),paramij->powermint-1) * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk(cos_theta,paramij);
  ters_gijk_d(paramij,rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

  // compute the derivative wrt Ri
  // dri = -dfc*gijk*ex_delr*rik_hat;
  // dri += fc*gijk_d*ex_delr*dcosdri;
  // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);

  vec3_scale(-dfc*gijk*ex_delr,rik_hat,dri);
  vec3_scaleadd(fc*ex_delr,dcosdri,dri,dri);
  vec3_scaleadd(fc*gijk*ex_delr_d,rik_hat,dri,dri);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rij_hat,dri,dri);
  vec3_scale(prefactor,dri,dri);

  // compute the derivative wrt Rj
  // drj = fc*gijk_d*ex_delr*dcosdrj;
  // drj += fc*gijk*ex_delr_d*rij_hat;

  vec3_scale(fc*ex_delr,dcosdrj,drj);
  vec3_scaleadd(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  vec3_scale(prefactor,drj,drj);

  // compute the derivative wrt Rk
  // drk = dfc*gijk*ex_delr*rik_hat;
  // drk += fc*gijk_d*ex_delr*dcosdrk;
  // drk += -fc*gijk*ex_delr_d*rik_hat;

  vec3_scale(dfc*gijk*ex_delr,rik_hat,drk);
  vec3_scaleadd(fc*ex_delr,dcosdrk,drk,drk);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rik_hat,drk,drk);
  vec3_scale(prefactor,drk,drk);
}

/* ---------------------------------------------------------------------- */

void PairTERSOFFHG::ters_gijk_d(Param *param, double *rij_hat, double rij,
                             double *rik_hat, double rik,
                             double *dri, double *drj, double *drk)
{
  // derivative of cos(theta) in gijk term; Eq. A14
  // first element is devative wrt Ri, second wrt Rj, third wrt Rk

  double cos_theta = vec3_dot(rij_hat,rik_hat);

  const double ters_c = param->c;
  const double ters_d = param->d;
  const double hcth = param->h - cos_theta;
  const double tmp = 2.0 * ters_d * hcth;

  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,drj);
  vec3_scale(-tmp/rij,drj,drj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,drk);
  vec3_scale(-tmp/rik,drk,drk);
  vec3_add(drj,drk,dri);
  vec3_scale(-1.0,dri,dri);
}

/* ---------------------------------------------------------------------- */

double PairTERSOFFHG::ters_bij(double zeta, Param *param, int nncl, int nnsi)
{
  // REBO1 bond order term; Eq. A12

  double Hcoord;
  if (nnsi > 4) nnsi = 4;
  if (nnsi > 0 && nncl > 0) Hcoord = coordenergy[nnsi][nncl];
  else Hcoord = 0.0;

  double powern = param->powern;
  double powereta = param->powereta;
  double tmp = zeta + Hcoord;
  if (tmp < 0.0) tmp = 1.0e-6;

  return pow(1.0 + pow(tmp,powereta),-powern);
}

/* ---------------------------------------------------------------------- */

double PairTERSOFFHG::ters_bij_d(double zeta, Param *param, int nncl, int nnsi)
{
  // derivative of REBO1 bond order term; Eq. A12
  
  double Hcoord;
  if (nnsi > 4) nnsi = 4;
  if (nnsi > 0 && nncl > 0) Hcoord = coordenergy[nnsi][nncl];
  else Hcoord = 0.0;

  double powern = param->powern;
  double powereta = param->powereta;
  double tmp = zeta + Hcoord;
  if (tmp < 0.0) tmp = 1.0e-6;

  return -powereta * powern * pow(tmp,powereta-1) * 
  	pow((1+pow(tmp,powereta)),-powern-1);
}

/* ---------------------------------------------------------------------- */

void PairTERSOFFHG::attractive(Param *paramij, Param *paramik, double prefactor,
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

  ters_zetaterm_d(paramij, paramik, prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk);
}

/* ---------------------------------------------------------------------- */

int PairTERSOFFHG::pack_forward_comm(int n, int *list, double *buf,
                                 int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  if (pack_flag == 1) {
    for (i = 0; i < n; i ++) {
      j = list[i];
      buf[m++] = NCl[j];
    }
  } else if (pack_flag == 2) {
    for (i = 0; i < n; i ++) {
      j = list[i];
      buf[m++] = NSi[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairTERSOFFHG::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n ;
  if (pack_flag == 1) {
    for (i = first; i < last; i++)
      NCl[i] = buf[m++];
  } else if (pack_flag == 2) {
    for (i = first; i < last; i++)
      NSi[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int PairTERSOFFHG::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (pack_flag == 1) {
    for (i = first; i < last; i++)
      buf[m++] = NCl[i];
  } else if (pack_flag == 2) {
    for (i = first; i < last; i++)
      buf[m++] = NSi[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairTERSOFFHG::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  if (pack_flag == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      NCl[j] += buf[m++];
    }
  } else if (pack_flag == 2) {
    for (i = 0; i < n; i++) {
      j = list[i];
      NSi[j] += buf[m++];
    }
  }
}


void * PairTERSOFFHG::xmalloc(size_t n){
  void *p;
  if ((p=malloc(n)) == NULL){
    fprintf(stderr, "xmalloc:memory allocation\n");
    exit(1);
  }
  return p;
}

double PairTERSOFFHG::cSiF[4] = {};

//void PairTERSOFFHG::bicubic_genCoef (double y[X1_NGRIDPOINTS][X2_NGRIDPOINTS],
//          double y1[X1_NGRIDPOINTS][X2_NGRIDPOINTS], 
//          double y2[X1_NGRIDPOINTS][X2_NGRIDPOINTS],
//          double y12[X1_NGRIDPOINTS][X2_NGRIDPOINTS])//,double (****c)[4][4])
void PairTERSOFFHG::bicubic_genCoef()
{
  int i, j;
    for (i = 0; i < 4; i++)
    {
      cSiF[i] = 0;
      std::cerr<<cSiF[i]<<std::endl;
    }

  // double z[4], z1[4], z2[4], z12[4];
  // double (****c)[4][4] = &PairTERSOFFHG::cSiF;
  // /* Initialize *c as a matrix of matrices */
  // (*c) = (double (***)[4][4])xmalloc(ms*sizeof(double (**)[4][4]));
  // for (i = 0; i < ms; i++)
  // {
  //   (*c)[i] = (double (**)[4][4])xmalloc(ns*sizeof(double (*)[4][4]));
  //   for (j = 0; j < ns; j++)
  //     ((*c)[i])[j] = (double (*)[4][4])xmalloc(16*sizeof(double));
  //double z[4], z1[4], z2[4], z12[4];
  //double (****c)[4][4] = &PairTERSOFFHG::cSiF;
  /* Initialize *c as a matrix of matrices */
  //(*c) = (double (***)[4][4])xmalloc(ms*sizeof(double (**)[4][4]));
  //for (i = 0; i < ms; i++)
  //{
  //  (*c)[i] = (double (**)[4][4])xmalloc(ns*sizeof(double (*)[4][4]));
  //  for (j = 0; j < ns; j++)
  //    ((*c)[i])[j] = (double (*)[4][4])xmalloc(16*sizeof(double));
  //}
}
// void Psif_genCoef (void){
//   int i, j;
//   double pp;
//   double p[X1_NGRIDPOINTS][X2_NGRIDPOINTS]; 
//   double p1[X1_NGRIDPOINTS][X2_NGRIDPOINTS]; 
//   double p2[X1_NGRIDPOINTS][X2_NGRIDPOINTS]; 
//   double p12[X1_NGRIDPOINTS][X2_NGRIDPOINTS];
  
//   for (i = 0; i < X1_NGRIDPOINTS; i++)
//     for (j = 0; j < X2_NGRIDPOINTS; j++)
//     {
//       p[i][j] = 0;
//       p1[i][j] = 0;
//       p2[i][j] = 0;
//       p12[i][j] = 0.0;
//     }
  
//   // ifstream fin("/home/dhumbird/fcmd/src/bicubic/Si_F", ios::in);
//   // while (!fin.eof())
//   // {
//   //   fin>>i>>j>>pp;
//   //   p[i][j] = pp;
//   // }
//   // fin.close();

//   /* Set values of Pcf at the various integer grid points. */
  
//   //p[NF][NC]
//   //Fit to Walch's bond strengths
//   p[0][0] = 0; //SiF;
//   p[1][0] = -0.1755; //SiF2 
//   p[2][0] = .496; //SiF3 changed to 4.23 eV to make exothermic
//   p[3][0] = -.181; //SiF4
//   //****cluster matching
//   // change 4
//   p[0][3] = -0.12;
//   p[0][2] = 0.085;
//   p[2][1] = -.07;//-0.06;
//   p[1][2] = -0.08;//0.2;//-0.0845;
//   //these werent used
//   //p[1][3] = -.21;//-.25;
//   //p[2][2] = -.1;
//   /* Derivative values by centered finite difference */
//   p1[1][1] = 0.5*(p[2][1] - p[0][1]);
//   p1[1][0] = 0.5*(p[2][0] - p[0][0]); 
//   p1[2][0] = 0.5*(p[3][0] - p[1][0]);
//   p2[0][1] = 0.5*(p[0][2] - p[0][0]);
//   p2[0][2] = 0.5*(p[0][3] - p[0][1]);
//   p2[1][1] = 0.5*(p[1][2] - p[1][0]); 
  
//   /* Create the set of 4x4 coefficient matrices */
// //  bicubic_genCoef(p, p1, p2, p12, &cSiF);
// }

