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
#include <fstream>
#include <malloc.h>
#include <map>

#include "pair_tersoffHG.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "math_special.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "potential_file_reader.h"
#include "suffix.h"
#include "tokenizer.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;
using namespace MathExtra;

#define MAXLINE 1024
#define DELTA 4
#define TOL 1.0e-9
#define PGDELTA 1

/* ---------------------------------------------------------------------- */

PairTERSOFFHG::PairTERSOFFHG(LAMMPS *lmp) : PairTersoff(lmp)
{

  nmax = 0;

  // set comm size needed by this Pair
  comm_forward = 1;
  comm_reverse = 1;

  maxlocal = 0;
  REBO_numneigh = NULL;
  REBO_firstneigh = NULL;
  ipage = NULL;
  pgsize = oneatom = 0;
}

/* ---------------------------------------------------------------------- */

void PairTERSOFFHG::read_file(char *file)
{
//  int params_per_line = 28;
//  char **words = new char*[params_per_line+1];

  // memory->sfree(params);
  // params = NULL;
  // nparams = maxparam = 0;

  // // open file on proc 0

  // FILE *fp;
  // if (comm->me == 0) {
  //   fp = force->open_potential(file);
  //   if (fp == NULL) {
  //     char str[128];
  //     sprintf(str,"Cannot open TersoffHG potential file %s",file);
  //     error->one(FLERR,str);
  //   }
  // }

  // read each line out of file, skipping blank lines or leading '#'
  // store line of params if all 3 element tags are in element list

  // int n,nwords,ielement,jelement,kelement;
  // char line[MAXLINE],*ptr;
  // int eof = 0;

  // while (1) {
  //   if (comm->me == 0) {
  //     ptr = fgets(line,MAXLINE,fp);
  //     if (ptr == NULL) {
  //       eof = 1;
  //       fclose(fp);
  //     } else n = strlen(line) + 1;
  //   }
  //   MPI_Bcast(&eof,1,MPI_INT,0,world);
  //   if (eof) break;
  //   MPI_Bcast(&n,1,MPI_INT,0,world);
  //   MPI_Bcast(line,n,MPI_CHAR,0,world);

  //   // strip comment, skip line if blank

  //   if ((ptr = strchr(line,'#'))) *ptr = '\0';
  //   nwords = atom->count_words(line);
  //   if (nwords == 0) continue;

  //   // concatenate additional lines until have params_per_line words

  //   while (nwords < params_per_line) {
  //     n = strlen(line);
  //     if (comm->me == 0) {
  //       ptr = fgets(&line[n],MAXLINE-n,fp);
  //       if (ptr == NULL) {
  //         eof = 1;
  //         fclose(fp);
  //       } else n = strlen(line) + 1;
  //     }
  //     MPI_Bcast(&eof,1,MPI_INT,0,world);
  //     if (eof) break;
  //     MPI_Bcast(&n,1,MPI_INT,0,world);
  //     MPI_Bcast(line,n,MPI_CHAR,0,world);
  //     if ((ptr = strchr(line,'#'))) *ptr = '\0';
  //     nwords = atom->count_words(line);
  //   }

  //   if (nwords != params_per_line)
  //     error->all(FLERR,"Incorrect format in TERSOFFHG potential file");

  //   // words = ptrs to all words in line

  //   nwords = 0;
  //   words[nwords++] = strtok(line," \t\n\r\f");
  //   while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    // ielement,jelement,kelement = 1st args
    // if all 3 args are in element list, then parse this line
    // else skip to next line
  memory->sfree(params);
  params = nullptr;
  nparams = maxparam = 0;

  // open file on proc 0

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, file, "tersoffHG", unit_convert_flag);
    char *line;

    // transparently convert units for supported conversions

    int unit_convert = reader.get_unit_convert();
    double conversion_factor = utils::get_conversion_factor(utils::ENERGY,
                                                            unit_convert);
    while ((line = reader.next_line(NPARAMS_PER_LINE))) {
      try {
        ValueTokenizer values(line);

        std::string iname = values.next_string();
        std::string jname = values.next_string();
        std::string kname = values.next_string();

        // ielement,jelement,kelement = 1st args
        // if all 3 args are in element list, then parse this line
        // else skip to next entry in file
        int ielement, jelement, kelement;
        for (ielement = 0; ielement < nelements; ielement++)
          if (iname == elements[ielement]) break;
        if (ielement == nelements) continue;
        for (jelement = 0; jelement < nelements; jelement++)
          if (jname == elements[jelement]) break;
        if (jelement == nelements) continue;
        for (kelement = 0; kelement < nelements; kelement++)
          if (kname == elements[kelement]) break;
        if (kelement == nelements) continue;

    // load up parameter settings and error check their values
        if (nparams == maxparam) {
          maxparam += DELTA;
          params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                              "pair:params");
          memset(params + nparams, 0, DELTA*sizeof(Param));

        }
        params[nparams].ielement = ielement;
        params[nparams].jelement = jelement;
        params[nparams].kelement = kelement;
        params[nparams].powerm    = values.next_double(); // beta (A13)
        params[nparams].gamma = values.next_double(); //atof(words[4]);     // not used
        params[nparams].lam3 = values.next_double(); //atof(words[5]);      // alpha (A13)
        params[nparams].c = values.next_double(); //atof(words[6]);         // c (A14)
        params[nparams].d = values.next_double(); //atof(words[7]);         // d (A14)
        params[nparams].h = values.next_double(); //atof(words[8]);         // h (A14)
        params[nparams].powern = values.next_double(); //atof(words[9]);    // delta (A12)
        params[nparams].beta = values.next_double(); //atof(words[10]);     // not used
        params[nparams].lam2 = values.next_double(); //atof(words[11]);     // mu (A4)
        params[nparams].bigb = values.next_double(); //atof(words[12]);     // B (A4)
        params[nparams].bigr = values.next_double(); //atof(words[13]);     // rmin (A5)
        params[nparams].bigd = values.next_double(); //atof(words[14]);     // rmax - rmin (A5)
        params[nparams].lam1 = values.next_double(); //atof(words[15]);     // lambda (A3)
        params[nparams].biga = values.next_double(); //atof(words[16]);     // A (A3)
        params[nparams].powereta = values.next_double(); //atof(words[17]); // eta (A12)
        params[nparams].Z_i = values.next_double(); //atof(words[18]);      // Z_i (A7)
        params[nparams].Z_j = values.next_double(); //atof(words[19]);      // Z_j (A7)
        params[nparams].spl_ra = values.next_double(); //atof(words[20]);   // spl_ra (A10)
        params[nparams].spl_rb = values.next_double(); //atof(words[21]);   // spl_rb (A10)
        params[nparams].spl_a = values.next_double(); //atof(words[22]);    // spl_a (A10)
        params[nparams].spl_b = values.next_double(); //atof(words[23]);    // spl_b (A10)
        params[nparams].spl_c = values.next_double(); //atof(words[24]);    // spl_c (A10)
        params[nparams].spl_s = values.next_double(); //atof(words[25]);    // spl_s (A10)
        params[nparams].Re = values.next_double(); //atof(words[26]);       // Re (A13)
        std::string pfile=values.next_string();
        params[nparams].PxyFile = pfile; //words[27];       // File contaning P nodes
        params[nparams].powermint = int(params[nparams].powerm);

        if (unit_convert) {
          params[nparams].biga *= conversion_factor;
          params[nparams].bigb *= conversion_factor;
        }
      } catch (TokenizerException &e) {
        error->one(FLERR, e.what());
      }
        // currently only allow m exponent of 1 or 3


      //add a check to make sure the spline file exists
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
        error->one(FLERR,"Illegal REBO1 parameter");
      read_lib(&params[nparams]);
      nparams++;
    }
  }

  //delete [] words;
  MPI_Bcast(&nparams, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxparam, 1, MPI_INT, 0, world);

  if (comm->me != 0) {
    params = (Param *) memory->srealloc(params,maxparam*sizeof(Param), "pair:params");
  }

  MPI_Bcast(params, maxparam*sizeof(Param), MPI_BYTE, 0, world);
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
  double rcmaxsq;
  double rcmin, rcmax;


  // count number of nearest neighbors; needs communications
//  count_neigh();
  int n,allnum;
  double dS;
  int *neighptr;
  
  if (atom->nmax > maxlocal) {
    maxlocal = atom->nmax;
    memory->destroy(REBO_numneigh);
    memory->sfree(REBO_firstneigh);
    memory->create(REBO_numneigh,maxlocal,"AIREBO:numneigh");
    REBO_firstneigh = (int **) memory->smalloc(maxlocal*sizeof(int *),
                                               "AIREBO:firstneigh");
  }

  allnum = list->inum + list->gnum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // From pair_AIREBO: store all REBO neighs of owned and ghost atoms
  ipage->reset();
  Nmap.clear();
  
  for (ii = 0; ii < allnum; ii++) 
  {
    i = ilist[ii];

    n = 0;
    neighptr = ipage->vget();

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = map[type[i]];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) 
    {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];
      iparam_ij = elem3param[itype][jtype][jtype];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      rcmin = params[iparam_ij].bigr;
      rcmax = rcmin + params[iparam_ij].bigd;
      rcmaxsq = rcmax*rcmax;
      if (rsq < rcmaxsq){//[itype][jtype]) {
        neighptr[n++] = j;
        Nmap[i][jtype]+=Sp(sqrt(rsq),rcmin,rcmax,dS);
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

    for (jj = 0; jj < REBO_numneigh[i]; jj++) 
    {
      j = jlist[jj];
      jtag = tag[j];

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
        //std::cerr<<tag[i]<<" "<<tag[j]<<" "<<rij<<" "<<std::endl;

      Rij1 = Rij/rij;
      jtype = map[type[j]];
      iparam_ij = elem3param[itype][jtype][jtype];

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
//      std::cerr<<tag[i]<<" "<<tag[j]<<" "<<VA_ij<<"\n";

      double b_ij=BondOrder(i, j, ij_f, VA_ij, eflag, vflag);
//      std::cerr<<tag[i]<<" "<<tag[j]<<" "<<VA_ij<<"\n";

      evdwl = VR_ij - b_ij*VA_ij;
      fpair = (-dVR_ij+b_ij*dVA_ij)/rij;
      HGvector Fij=(-dVR_ij+b_ij*dVA_ij)*Rij1;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                        evdwl, 0.0, fpair, Rij.x, Rij.y, Rij.z);

      f[i][0] += Fij.x;
      f[i][1] += Fij.y;
      f[i][2] += Fij.z;

      f[j][0] -= Fij.x;
      f[j][1] -= Fij.y;
      f[j][2] -= Fij.z;
    }
  }
}


//******************************************************
double PairTERSOFFHG::BondOrder(int i, int j, double ij_f, double Pre, int eflag, int vflag){
  HGvector Rij, Rik, Rjk, Rij1, Rik1, Rjk1, Fij;
  double rij, rik, rjk, bbar_ij;
  double ik_f, jk_f, ik_fprime, ij_fprime, jk_fprime;
  iFv_ij=0; iFv_ik.clear(); iFv_jk.clear(); //iFv_ik2.clear();
  jFv_ij=0; jFv_ik.clear(); jFv_jk.clear(); //jFv_jk2.clear();
  double cos_theta, bigr, bigd;
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
  iparam_ij = elem3param[itype][jtype][jtype];
  iparam_ji = elem3param[jtype][itype][itype];
    
  Rij.set(x[i][0]-x[j][0], x[i][1]-x[j][1], x[i][2]-x[j][2]);
  rij=Rij.mag();
  Rij1=Rij/rij;
 
  //double NF_ij, NC_ij, NSi_ij, Nt_ij, NF_ji, NC_ji, NSi_ji, Nt_ji;
  //double NCl_ij, NCl_ji;
  double P_ij, P_ji, dFP_ij, dFP_ji, dCP_ij, dCP_ji;
  double el, elf, alpha, beta, Re, n, dlam, g, g1;//, xik, Yik, Y1ik, F_ij, dFi_ij, dFj_ij, dFc_ij;
  double eta_i, eta_j, delta_i, delta_j;

  // NC_ij  = i->Nmap[6]  - (j->id==6)  * bond_ij->f;
  // NSi_ij = i->Nmap[14] - (j->id==14) * bond_ij->f;
  // NCl_ij = i->Nmap[17] - (j->id==17) * bond_ij->f;
  // NF_ji  = j->Nmap[9]  - (i->id==9)  * bond_ij->f;
  // NC_ji  = j->Nmap[6]  - (i->id==6)  * bond_ij->f;
  // NSi_ji = j->Nmap[14] - (i->id==14) * bond_ij->f;
  // NCl_ji = j->Nmap[17] - (i->id==17) * bond_ij->f;

  // Nt_ij = NF_ij + NC_ij + NSi_ij + NCl_ij;
  // Nt_ji = NF_ji + NC_ji + NSi_ji + NCl_ji;

  P_ij=P_ji=dFP_ij=dFP_ji=dCP_ij=dCP_ji=0;
  //std::cerr<<"BONDORDER"<<std::endl;

  //these only work if 0 is Si, 1 is F/Cl. Later, find a way to virtualize
  bicubicint (Nmap[i][1]-(jtype==1)*ij_f, Nmap[i][0]-(jtype==0)*ij_f, 
              &P_ij, &dFP_ij, &dCP_ij, &params[iparam_ij]);
  //std::cerr<<" "<<Nmap[i][1]-(jtype==1)*ij_f<<" "<<Nmap[i][0]-(jtype==0)*ij_f<<" "<<
  //            P_ij<<" "<<dFP_ij<<" "<<dCP_ij<<std::endl;

  bicubicint (Nmap[j][1]-(itype==1)*ij_f, Nmap[j][0]-(itype==0)*ij_f,
              &P_ji, &dFP_ji, &dCP_ji, &params[iparam_ji]);
  
  //std::cerr<<" "<<Nmap[j][1]-(itype==1)*ij_f<<" "<<Nmap[j][0]-(itype==0)*ij_f<<" "<<
  //            P_ji<<" "<<dFP_ij<<" "<<dCP_ij<<std::endl;  
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
  bsp_ij=P_ij; bsp_ji=P_ji;
  
  // //***************for k neighbor of i*************************
  klist = REBO_firstneigh[i];
  for (kk = 0; kk < REBO_numneigh[i]; kk++) 
  {
    k = klist[kk];
    k &= NEIGHMASK;
    ktag = tag[k];
    ktype = map[type[k]];
    iparam_ik = elem3param[itype][ktype][ktype];
    iparam_ijk = elem3param[itype][jtype][ktype];
    Rik.set(x[i][0]-x[k][0], x[i][1]-x[k][1], x[i][2]-x[k][2]);
    rik=Rik.mag();
    Rik1=Rik/rik;
    if (tag[k]!=tag[j] && Rik.sqmag() <= params[iparam_ik].cutsq)
    {
      Rjk.set(x[j][0]-x[k][0], x[j][1]-x[k][1], x[j][2]-x[k][2]);
      rjk=Rjk.mag();
      Rjk1=Rjk/rjk;
      bigr= params[iparam_ik].bigr;
      bigd = params[iparam_ik].bigd;
      ik_f = Sp(rik, bigr, bigr+bigd, ik_fprime);
      // //*******b-sigma-pi calculations************
      alpha = params[iparam_ijk].lam3;
      if (alpha)
      {
        beta = params[iparam_ijk].powerm;
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
           (g*el + dFP_ij + dCP_ij)
//           (g*el + 0)//((k->id==9||k->id==17) ? dFP_ij : dCP_ij))
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
    iparam_jk = elem3param[jtype][ktype][ktype];
    iparam_ijk = elem3param[jtype][itype][ktype];
    Rjk.set(x[j][0]-x[k][0], x[j][1]-x[k][1], x[j][2]-x[k][2]);
    rjk=Rjk.mag();
    Rjk1=Rjk/rjk;

    if (tag[k]!=tag[i] && Rjk.sqmag() <= params[iparam_jk].cutsq)
    {
      Rik.set(x[i][0]-x[k][0], x[i][1]-x[k][1], x[i][2]-x[k][2]);
      rik=Rik.mag();
      Rik1=Rik/rik;
      bigr= params[iparam_jk].bigr;
      bigd = params[iparam_jk].bigd;
      jk_f = Sp(rjk, bigr, bigr+bigd, jk_fprime);
      // //*******b-sigma-pi calculations************
      alpha = params[iparam_ijk].lam3;
      if (alpha)
      {
        beta = params[iparam_ijk].powerm;
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
           (g*el + dFP_ji + dCP_ji)
//           (g*el + 0)//((k->id==9||k->id==17) ? dFP_ji : dCP_ji))
           +elf*(g1*(1/rij - cos_theta/rjk) - g*dlam));
      jFv_ik.push_back(elf*(-g1*rik/rjk/rij));
            
      bsp_ji += g*elf;
    }
  }

  //***********************************************

  Pre*=0.5;
  eta_i = params[iparam_ij].powereta;
  eta_j = params[iparam_ji].powereta;
  delta_i = params[iparam_ij].powern;
  delta_j = params[iparam_ji].powern;

  bbar_ij = 0.5*(pow(1 + pow(bsp_ij, eta_i), -delta_i)
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
    else bsp_ji=pow(1 + pow(bsp_ji, eta_j), -delta_j-1)
                *eta_j*pow(bsp_ji, eta_j-1);
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
    iparam_ik = elem3param[itype][ktype][ktype];
    iparam_ijk = elem3param[itype][jtype][ktype];
    Rik.set(x[i][0]-x[k][0], x[i][1]-x[k][1], x[i][2]-x[k][2]);
    rik=Rik.mag();
    Rik1=Rik/rik;
    if (tag[k]!=tag[j] && Rik.sqmag() <= params[iparam_ik].cutsq)
    {
      Rjk.set(x[j][0]-x[k][0], x[j][1]-x[k][1], x[j][2]-x[k][2]);
      rjk=Rjk.mag();
      Rjk1=Rjk/rjk;
      force_ik=bsp_ij*iFv_ik[cnt];
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
      Fij =  force_jk * Rjk1; 

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
    iparam_jk = elem3param[jtype][ktype][ktype];
    iparam_ijk = elem3param[jtype][itype][ktype];
    Rjk.set(x[j][0]-x[k][0], x[j][1]-x[k][1], x[j][2]-x[k][2]);
    rjk=Rjk.mag();
    Rjk1=Rjk/rjk;
    if (tag[k]!=tag[i] && Rjk.sqmag() <= params[iparam_jk].cutsq)
    {
      Rik.set(x[i][0]-x[k][0], x[i][1]-x[k][1], x[i][2]-x[k][2]);
      rik=Rik.mag();
      Rik1=Rik/rik;
      force_jk=bsp_ji*jFv_jk[cnt];
      Fij = force_jk*Rjk1;

      if (evflag) ev_tally(j,k,nlocal,newton_pair,
            0.0,0.0,force_jk/rjk,Rjk.x,Rjk.y,Rjk.z);

      f[j][0] += Fij.x;
      f[j][1] += Fij.y;
      f[j][2] += Fij.z;

      f[k][0] -= Fij.x;
      f[k][1] -= Fij.y;
      f[k][2] -= Fij.z;
      
      force_ik = bsp_ji * jFv_ik[cnt];
      Fij =  force_ik * Rik1; 

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

//************************************************************************
void PairTERSOFFHG::read_lib(Param *param)
{
  unsigned int maxlib = 1024;
  int i,j,k,l,nwords,m;
  double p[X1_NGRIDPOINTS][X2_NGRIDPOINTS]; 
  std::string splfile = param->PxyFile;
  double pp;

  for (i = 0; i < X1_NGRIDPOINTS; i++)
    for (j = 0; j < X2_NGRIDPOINTS; j++)
      p[i][j] = 0;

  if (splfile.compare("NULL")!=0 && comm->me == 0) 
  {
//    std::cerr<<splfile<<std::endl;
    std::ifstream fp(splfile, std::ios::in);
    if (!fp)
    {
//      char str[128];
//      sprintf(str,"Cannot open specified spline file");
      error->one(FLERR,"Cannot open specified spline file");
    }
    while (!fp.eof())
    {
      fp>>i>>j>>pp;
      p[i][j] = pp;
    }
    fp.close();
  }

  bicubic_genCoef(p, param);
  // need to add BCAST HERE e.g.  
  // MPI_Bcast(&cSiHal[0][0],1604,MPI_DOUBLE,0,world);

}

void PairTERSOFFHG::bicubic_genCoef (double y[X1_NGRIDPOINTS][X2_NGRIDPOINTS], Param *param)
{
  int i, j;
  double z[4], z1[4], z2[4], z12[4];
  double y1[X1_NGRIDPOINTS][X2_NGRIDPOINTS]; 
  double y2[X1_NGRIDPOINTS][X2_NGRIDPOINTS]; 
  double y12[X1_NGRIDPOINTS][X2_NGRIDPOINTS];
  for (i = 0; i < X1_NGRIDPOINTS; i++)
    for (j = 0; j < X2_NGRIDPOINTS; j++)
      y1[i][j] = y2[i][j] = y12[i][j] = 0.0;
  
  y1[1][0] = 0.5*(y[2][0] - y[0][0]); 
  y1[1][1] = 0.5*(y[2][1] - y[0][1]);
  y1[2][0] = 0.5*(y[3][0] - y[1][0]);
 
  y2[0][1] = 0.5*(y[0][2] - y[0][0]);
  y2[0][2] = 0.5*(y[0][3] - y[0][1]);
  y2[1][1] = 0.5*(y[1][2] - y[1][0]); 
  
  /* Initialize PxyInterp as a matrix of matrices */
  param->PxyInterp = (double (***)[4][4])xmalloc(X1_NGRIDSQUARES*sizeof(double (**)[4][4]));
  for (i = 0; i < X1_NGRIDSQUARES; i++){
    param->PxyInterp[i] = (double (**)[4][4])xmalloc(X2_NGRIDSQUARES*sizeof(double (*)[4][4]));
    for (j = 0; j < X2_NGRIDSQUARES; j++)
      (param->PxyInterp[i])[j] = (double (*)[4][4])xmalloc(16*sizeof(double));
  }
  /* For each grid square (i,j) compute the 4x4 matrix c_(ij).
   * 
   *          (4)------(3)
   *           |        |
   *           |        |
   *           |        |
   *          (1)------(2)
   */
  for (i = 0; i < X1_NGRIDSQUARES; i++){
    for (j = 0; j < X2_NGRIDSQUARES; j++){
      /* The array z[] holds the four function values at the
       * corners of this grid square. */
      z[0] = y[i][j];
      z[1] = y[i+1][j];
      z[2] = y[i+1][j+1];
      z[3] = y[i][j+1];
      /* The array z1[] holds the first derivative of the function
       * with respect to x1 at each of the corners of this grid square.*/
      z1[0] = y1[i][j];
      z1[1] = y1[i+1][j];
      z1[2] = y1[i+1][j+1];
      z1[3] = y1[i][j+1];
      /* The array z2[] holds the first derivative of the function
       * with respect to x2 at each of the corners of this grid square.*/
      z2[0] = y2[i][j];
      z2[1] = y2[i+1][j];
      z2[2] = y2[i+1][j+1];
      z2[3] = y2[i][j+1];
      /* The array z12[] holds the cross derivative of the function
       * at each of the corners of this grid square. */
      z12[0] = y12[i][j];
      z12[1] = y12[i+1][j];
      z12[2] = y12[i+1][j+1];
      z12[3] = y12[i][j+1];     
      /* Use bcucof() to compute the 4x4 coefficient matrix c_(ij) */
//      bcucof(z, z1, z2, z12, 1.0, 1.0, *((*c)[i][j]));
      bcucof(z, z1, z2, z12, 1.0, 1.0, *(param->PxyInterp[i][j]));
    }
  }
}

/* bcucof:  (Numerical Recipes, 2nd Ed., p 126)  Generates the
 * 16-member coefficient matrix c_(ij), i=1->4, j=1->4.  c_(ij)
 * is used in the "bcuint" function to interpolate a function
 * value at a specified off-lattice point. */
void PairTERSOFFHG::bcucof(double y[], double y1[], double y2[], double y12[],
            double d1, double d2, double c[4][4]){
  int wt[16][16] = /* Weighting factors; pretabulated long long ago */
  {
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    -3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0,
    2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1,
    0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1,
    -3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
    9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2,
    -6 ,6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2,
    2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0,
    -6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1,
    4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1
  };
  int i, j, k, l;
  double xx, d1d2, cl[16], x[16];
  
  d1d2 = d1*d2;
  for (i = 0; i < 4; i++){ /* Build a temporary vector */
    x[i] = y[i];
    x[i+4] = y1[i]*d1;
    x[i+8] = y2[i]*d2;
    x[i+12] = y12[i]*d1d2;
  }
  //Multiply the matrix wt[][] by vector x[] and store result in vector cl[]
  for (i = 0; i < 16; i++){ 
    xx = 0.0;
    for (k = 0; k < 16; k++) xx += wt[i][k]*x[k];
    cl[i] = xx;
  }
  //Place members of cl[] into their proper locations in c[][].
  l = 0;
  for (i = 0; i < 4; i++) 
    for (j = 0; j < 4; j++) {
      c[i][j] = cl[l++];
    }
}

void PairTERSOFFHG::bicubicint (double x1, double x2,
         double *y, double *y1, double *y2, Param *param)//double (****c)[4][4])
{
  int i, j;
  double x1max = X1_NGRIDSQUARES-1.e-10, x2max = X1_NGRIDSQUARES-1.e-10;
  x1 = (x1 > x1max ? x1max : x1);
  x2 = (x2 > x2max ? x2max : x2);
  i = (int)x1;
  j = (int)x2;
//  std::cerr<<param->PxyInterp[i][j]<<std::endl;
  bcuint(0.0, 1.0, 0.0, 1.0, x1-i, x2-j, y, y1, y2, *(param->PxyInterp[i][j]));
}

void PairTERSOFFHG::bcuint(double x1l, double x1u, double x2l, double x2u, 
      double x1, double x2, double *ansy, double *ansy1, double *ansy2,
      double c[4][4]){
  int i;
  double t, u, d1, d2;
  d1 = x1u - x1l;
  d2 = x2u - x2l;
  if (!d1 || !d2){
    fprintf(stderr, "2D Cubic Interpolation: Bad gridpoint coords:\n");
    fprintf(stderr, "\tx1u(%.5e) x1l(%.5e) x2u(%.5e) x2l(%.5e)\n", 
      x1u, x1l, x2u, x2l);
    fprintf(stderr, "Program exits\n");
    exit(0);
  }
  
  t = (x1-x1l)/d1;
  u = (x2-x2l)/d2;
  *ansy = (*ansy2) = (*ansy1) = 0.0;

  for (i = 3; i >= 0; i--){
    *ansy = t*(*ansy) + ((c[i][3]*u + c[i][2])*u + c[i][1])*u + c[i][0];
    *ansy2 = t*(*ansy2) + (3.0*c[i][3]*u + 2.0*c[i][2])*u + c[i][1];
    *ansy1 = u*(*ansy1) + (3.0*c[3][i]*t + 2.0*c[2][i])*t + c[1][i];
  }
  *ansy1 /= d1;
  *ansy2 /= d2;
}

/* ---------------------------------------------------------------------- */

// void PairTERSOFFHG::count_neigh()
// {
//   int i, j, ii, jj, itype, jtype, n;
//   int inum, jnum, param, *ilist, *jlist, *numneigh, **firstneigh;
//   double r, rsq, delrij[3];
//   const double cutshortsq = cutmax*cutmax;

//   double **x = atom->x;
//   int *type = atom->type;

//   if (atom->nmax > nmax) {
//     nmax = atom->nmax;
//     memory->grow(NCl,nmax,"pair:NCl");
//     memory->grow(NSi,nmax,"pair:NSi");
//   }

//   inum  = list->inum;
//   ilist = list->ilist;
//   numneigh = list->numneigh;
//   firstneigh = list->firstneigh;
//   for (ii = 0; ii < inum; ii++) {
//     i = ilist[ii];
//     itype = map[type[i]];
//     jlist = firstneigh[i];
//     jnum = numneigh[i];

//     // skip immediately if center atom is not Si
//     if (strcmp(elements[map[itype+1]],"Si") != 0) continue;

//     NCl[i] = NSi[i] = 0.0;

//     for (jj = 0; jj < jnum; jj++) {
//       j = jlist[jj] & NEIGHMASK;
//       jtype = map[type[j]];


//       delrij[0] = x[i][0] - x[j][0];
//       delrij[1] = x[i][1] - x[j][1];
//       delrij[2] = x[i][2] - x[j][2];
//       rsq = vec3_dot(delrij,delrij);
//       param = elem3param[itype][jtype][jtype];

//       if (rsq > cutshortsq) continue;

//       r = sqrt(rsq);
//       //if (strcmp(elements[map[jtype+1]],"Cl") == 0) NCl[i] += ters_fc(r,&params[param]);
//       //if (strcmp(elements[map[jtype+1]],"Si") == 0) NSi[i] += ters_fc(r,&params[param]);
//     }
//   }

//   // communicating coordination number to all nodes
//   pack_flag = 1;
//   comm->forward_comm_pair(this);
//   pack_flag = 2;
//   comm->forward_comm_pair(this);

// }

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
//       iparam_ij = elem3param[itype][jtype][jtype];
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
//       iparam_ij = elem3param[itype][jtype][jtype];
  
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
//         iparam_ijk = elem3param[itype][jtype][ktype];
//        iparam_ik  = elem3param[itype][ktype][ktype];

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
//         iparam_ijk = elem3param[itype][jtype][ktype];
//        iparam_ik  = elem3param[itype][ktype][ktype];

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
