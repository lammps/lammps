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
   Contributing authors: D.K. Ward (donward@sandia.gov) and X.W. Zhou (Sandia)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   The formulation for this work follows (a) D.G. Pettifor, et al., Mat.
   Sci. and Eng. A365, 2-13, (2004);(b) D.A. Murdick, et al., Phys.
   Rev. B 73, 045206 (2006);(c) D.G. Pettifor and I.I. Oleinik., Phys
   Rev. Lett. 84, 4124 (2000); (d) D.K. Ward, et al., Phys. Rev. B 85,
   115206 (2012).

   Copyright (2012) Sandia Corporation.  Under the terms of Contract DE-
   AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   rights in this software.

   pairbop v 1.0 comes with no warranty of any kind.  pairbop v 1.0 is a
   copyrighted code that is distributed free-of-charge, under the terms
   of the GNU Public License (GPL).  See "Open-Source
   Rules"_http://lammps.sandia.gov/open_source.html
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include "pair_bop.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "force.h"
#include "timer.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "math_special.h"

using namespace LAMMPS_NS;
using namespace MathSpecial;

#define MAXLINE 1024
#define EPSILON 1.0e-6

/* ---------------------------------------------------------------------- */

PairBOP::PairBOP(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  one_coeff = 1;
  manybody_flag = 1;

  map = NULL;
  pi_a = NULL;
  pro_delta = NULL;
  pi_delta = NULL;
  pi_p = NULL;
  pi_c = NULL;
  sigma_r0 = NULL;
  pi_r0 = NULL;
  phi_r0 = NULL;
  sigma_rc = NULL;
  pi_rc = NULL;
  phi_rc = NULL;
  r1 = NULL;
  sigma_beta0 = NULL;
  pi_beta0 = NULL;
  phi0 = NULL;
  sigma_n = NULL;
  pi_n = NULL;
  phi_m = NULL;
  sigma_nc = NULL;
  pi_nc = NULL;
  phi_nc = NULL;
  pro = NULL;
  sigma_delta = NULL;
  sigma_c = NULL;
  sigma_a = NULL;
  sigma_g0 = NULL;
  sigma_g1 = NULL;
  sigma_g2 = NULL;
  sigma_g3 = NULL;
  sigma_g4 = NULL;
  sigma_f = NULL;
  sigma_k = NULL;
  small3 = NULL;
  rcut = NULL;
  dr = NULL;
  rdr = NULL;
  disij = NULL;
  rij = NULL;
  cosAng = NULL;
  betaS = NULL;
  dBetaS = NULL;
  betaP = NULL;
  dBetaP = NULL;
  repul = NULL;
  dRepul = NULL;
  itypeSigBk = NULL;
  nSigBk = NULL;
  sigB = NULL;
  sigB1 = NULL;
  itypePiBk = NULL;
  nPiBk = NULL;
  piB = NULL;
  pBetaS = NULL;
  pBetaS1 = NULL;
  pBetaS2 = NULL;
  pBetaS3 = NULL;
  pBetaS4 = NULL;
  pBetaS5 = NULL;
  pBetaS6 = NULL;
  pBetaP = NULL;
  pBetaP1 = NULL;
  pBetaP2 = NULL;
  pBetaP3 = NULL;
  pBetaP4 = NULL;
  pBetaP5 = NULL;
  pBetaP6 = NULL;
  pRepul = NULL;
  pRepul1 = NULL;
  pRepul2 = NULL;
  pRepul3 = NULL;
  pRepul4 = NULL;
  pRepul5 = NULL;
  pRepul6 = NULL;
  FsigBO = NULL;
  FsigBO1 = NULL;
  FsigBO2 = NULL;
  FsigBO3 = NULL;
  FsigBO4 = NULL;
  FsigBO5 = NULL;
  FsigBO6 = NULL;
  rcmin = NULL;
  rcmax = NULL;
  rcmaxp = NULL;
  setflag = NULL;
  cutsq = NULL;
  cutghost = NULL;

  ghostneigh = 1;
  bt_sg=NULL;
  bt_pi=NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairBOP::~PairBOP()
{
  if(allocated) {
    memory_theta_destroy();
    if (otfly==0) memory->destroy(cos_index);
    delete [] map;

    memory->destroy(BOP_index);
    memory->destroy(rcut);
    memory->destroy(dr);
    memory->destroy(rdr);
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cutghost);
    memory->destroy(pBetaS);
    memory->destroy(pBetaS1);
    memory->destroy(pBetaS2);
    memory->destroy(pBetaS3);
    memory->destroy(pBetaS4);
    memory->destroy(pBetaS5);
    memory->destroy(pBetaS6);
    memory->destroy(pBetaP);
    memory->destroy(pBetaP1);
    memory->destroy(pBetaP2);
    memory->destroy(pBetaP3);
    memory->destroy(pBetaP4);
    memory->destroy(pBetaP5);
    memory->destroy(pBetaP6);
    memory->destroy(pRepul);
    memory->destroy(pRepul1);
    memory->destroy(pRepul2);
    memory->destroy(pRepul3);
    memory->destroy(pRepul4);
    memory->destroy(pRepul5);
    memory->destroy(pRepul6);
    memory->destroy(FsigBO);
    memory->destroy(FsigBO1);
    memory->destroy(FsigBO2);
    memory->destroy(FsigBO3);
    memory->destroy(FsigBO4);
    memory->destroy(FsigBO5);
    memory->destroy(FsigBO6);
    if(table==0) {
      memory->destroy(pi_a);
      memory->destroy(pro_delta);
      memory->destroy(pi_delta);
      memory->destroy(pi_p);
      memory->destroy(pi_c);
      memory->destroy(sigma_r0);
      memory->destroy(pi_r0);
      memory->destroy(phi_r0);
      memory->destroy(sigma_rc);
      memory->destroy(pi_rc);
      memory->destroy(phi_rc);
      memory->destroy(r1);
      memory->destroy(sigma_beta0);
      memory->destroy(pi_beta0);
      memory->destroy(phi0);
      memory->destroy(sigma_n);
      memory->destroy(pi_n);
      memory->destroy(phi_m);
      memory->destroy(sigma_nc);
      memory->destroy(pi_nc);
      memory->destroy(phi_nc);
      memory->destroy(pro);
      memory->destroy(sigma_delta);
      memory->destroy(sigma_c);
      memory->destroy(sigma_a);
      memory->destroy(sigma_g0);
      memory->destroy(sigma_g1);
      memory->destroy(sigma_g2);
      memory->destroy(sigma_g3);
      memory->destroy(sigma_g4);
      memory->destroy(sigma_f);
      memory->destroy(sigma_k);
      memory->destroy(small3);
    }
    else {
      memory->destroy(pi_a);
      memory->destroy(pro_delta);
      memory->destroy(pi_delta);
      memory->destroy(pi_p);
      memory->destroy(pi_c);
      memory->destroy(r1);
      memory->destroy(pro);
      memory->destroy(sigma_delta);
      memory->destroy(sigma_c);
      memory->destroy(sigma_a);
      memory->destroy(sigma_g0);
      memory->destroy(sigma_g1);
      memory->destroy(sigma_g2);
      memory->destroy(sigma_f);
      memory->destroy(sigma_k);
      memory->destroy(small3);
    }
  }
  if(allocate_sigma) {
    destroy_sigma();
  }
  if(allocate_pi) {
    destroy_pi();
  }
}

/* ---------------------------------------------------------------------- */

void PairBOP::compute(int eflag, int vflag)
{
  int ago,delay,every;
  int i,j,ii,jj,iij;
  int n,inum,temp_ij,ks;
  int itype,jtype;
  tagint i_tag,j_tag;
  int *ilist,*iilist,*numneigh;
  int **firstneigh;
  double dpr1,ps;
  double ftmp1,ftmp2,ftmp3,dE;
  double dis_ij[3],rsq_ij,r_ij;
  double betaS_ij,dBetaS_ij;
  double betaP_ij,dBetaP_ij;
  double repul_ij,dRepul_ij;
  double totE;

  double **f = atom->f;
  double **x = atom->x;
  int *type = atom->type;
  tagint *tag = atom->tag;
  int newton_pair = force->newton_pair;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  ago=neighbor->ago;
  delay=neighbor->delay;
  every=neighbor->every;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  // BOP Neighbor lists must be updated every time
  // atoms are moved between processors

  if ((ago ==0)||bop_step==0||(ago>=delay&&(ago%every)==0)||(nall>maxnall))
    gneigh();

  // For non on the fly calculations cos and derivatives
  // are calculated in advance and stored

  if(otfly==0) theta();
  else theta_mod();

  // Calculate Sigma Bond-Order

  if(a_flag==1) {
    if (otfly==0) sigmaBo_noa();
    else sigmaBo_noa_otf();
  }
  else {
    if (otfly==0) sigmaBo();
    else sigmaBo_otf();
  }

  // Calculate Pi Bond-Order

  if (otfly==0) PiBo();
  else PiBo_otf();

  n=0;
  totE=0;
  for (ii = 0; ii < inum; ii++) {
    i=ilist[ii];
    i_tag=tag[i];
    itype=map[type[i]]+1;
    iilist=firstneigh[i];
    for(jj=0;jj<numneigh[i];jj++) {
      temp_ij=BOP_index[i]+jj;
      j=iilist[jj];
      j_tag=tag[j];
      jtype=map[type[j]]+1;
      if(j_tag>=i_tag) {
        if(otfly==0) {
          if(neigh_flag[temp_ij]) {
            dpr1=(dRepul[temp_ij]-2.0*dBetaS[temp_ij]*sigB[n]
                -2.0*dBetaP[temp_ij]*piB[n])/rij[temp_ij];
            ftmp1=dpr1*disij[0][temp_ij];
            ftmp2=dpr1*disij[1][temp_ij];
            ftmp3=dpr1*disij[2][temp_ij];
            f[i][0]=f[i][0]+ftmp1;
            f[i][1]=f[i][1]+ftmp2;
            f[i][2]=f[i][2]+ftmp3;
            f[j][0]=f[j][0]-ftmp1;
            f[j][1]=f[j][1]-ftmp2;
            f[j][2]=f[j][2]-ftmp3;

            // add repulsive and bond order components to total energy
            // (d) Eq.1

            dE=-2.0*betaS[temp_ij]*sigB[n]-2.0*betaP[temp_ij]*piB[n];
            totE+=dE+repul[temp_ij];
            if(evflag) {
              ev_tally_full(i,repul[temp_ij],dE,0.0,0.0,0.0,0.0);
              ev_tally_full(j,repul[temp_ij],dE,0.0,0.0,0.0,0.0);
              ev_tally_xyz(i,j,nlocal,newton_pair,0.0,0.0,-ftmp1,-ftmp2,-ftmp3,
                  disij[0][temp_ij],disij[1][temp_ij],disij[2][temp_ij]);
            }
            n++;
          }
        }
        else {
          if(itype==jtype)
            iij=itype-1;
          else if(itype<jtype)
            iij=itype*bop_types-itype*(itype+1)/2+jtype-1;
          else
            iij=jtype*bop_types-jtype*(jtype+1)/2+itype-1;
          dis_ij[0]=x[j][0]-x[i][0];
          dis_ij[1]=x[j][1]-x[i][1];
          dis_ij[2]=x[j][2]-x[i][2];
          rsq_ij=dis_ij[0]*dis_ij[0]
              +dis_ij[1]*dis_ij[1]
              +dis_ij[2]*dis_ij[2];
          r_ij=sqrt(rsq_ij);
          if(r_ij<=rcut[iij]) {
            ps=r_ij*rdr[iij]+1.0;
            ks=(int)ps;
            if(nr-1<ks)
              ks=nr-1;
            ps=ps-ks;
            if(ps>1.0)
              ps=1.0;
            betaS_ij=((pBetaS3[iij][ks-1]*ps+pBetaS2[iij][ks-1])*ps
                +pBetaS1[iij][ks-1])*ps+pBetaS[iij][ks-1];
            dBetaS_ij=(pBetaS6[iij][ks-1]*ps+pBetaS5[iij][ks-1])*ps
                +pBetaS4[iij][ks-1];
            betaP_ij=((pBetaP3[iij][ks-1]*ps+pBetaP2[iij][ks-1])*ps
                +pBetaP1[iij][ks-1])*ps+pBetaP[iij][ks-1];
            dBetaP_ij=(pBetaP6[iij][ks-1]*ps+pBetaP5[iij][ks-1])*ps
                +pBetaP4[iij][ks-1];
            repul_ij=((pRepul3[iij][ks-1]*ps+pRepul2[iij][ks-1])*ps
                +pRepul1[iij][ks-1])*ps+pRepul[iij][ks-1];
            dRepul_ij=(pRepul6[iij][ks-1]*ps+pRepul5[iij][ks-1])*ps
                +pRepul4[iij][ks-1];
            dpr1=(dRepul_ij-2.0*dBetaS_ij*sigB[n]
                -2.0*dBetaP_ij*piB[n])/r_ij;
            ftmp1=dpr1*dis_ij[0];
            ftmp2=dpr1*dis_ij[1];
            ftmp3=dpr1*dis_ij[2];
            f[i][0]=f[i][0]+ftmp1;
            f[i][1]=f[i][1]+ftmp2;
            f[i][2]=f[i][2]+ftmp3;
            f[j][0]=f[j][0]-ftmp1;
            f[j][1]=f[j][1]-ftmp2;
            f[j][2]=f[j][2]-ftmp3;

            // add repulsive and bond order components to total energy
            // (d) Eq. 1

            dE=-2.0*betaS_ij*sigB[n]-2.0*betaP_ij*piB[n];
            totE+=dE+repul_ij;
            if(evflag) {
              ev_tally_full(i,repul_ij,dE,0.0,0.0,0.0,0.0);
              ev_tally_full(j,repul_ij,dE,0.0,0.0,0.0,0.0);
              ev_tally_xyz(i,j,nlocal,newton_pair,0.0,0.0,-ftmp1,-ftmp2,-ftmp3,
                  dis_ij[0],dis_ij[1],dis_ij[2]);
            }
            n++;
          }
        }
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
  bop_step = 1;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBOP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(rcut,npairs,"BOP:rcut");
  memory->create(dr,npairs,"BOP:dr");
  memory->create(rdr,npairs,"BOP:dr");
  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cutghost,n+1,n+1,"pair:cutghost");
  memory->create(pBetaS,npairs,nr,"BOP:pBetaS");
  memory->create(pBetaS1,npairs,nr,"BOP:pBetaS1");
  memory->create(pBetaS2,npairs,nr,"BOP:pBetaS2");
  memory->create(pBetaS3,npairs,nr,"BOP:pBetaS3");
  memory->create(pBetaS4,npairs,nr,"BOP:pBetaS4");
  memory->create(pBetaS5,npairs,nr,"BOP:pBetaS5");
  memory->create(pBetaS6,npairs,nr,"BOP:pBetaS6");
  memory->create(pBetaP,npairs,nr,"BOP:pBetaP");
  memory->create(pBetaP1,npairs,nr,"BOP:pBetaP1");
  memory->create(pBetaP2,npairs,nr,"BOP:pBetaP2");
  memory->create(pBetaP3,npairs,nr,"BOP:pBetaP3");
  memory->create(pBetaP4,npairs,nr,"BOP:pBetaP4");
  memory->create(pBetaP5,npairs,nr,"BOP:pBetaP5");
  memory->create(pBetaP6,npairs,nr,"BOP:pBetaP6");
  memory->create(pRepul,npairs,nr,"BOP:pRepul");
  memory->create(pRepul1,npairs,nr,"BOP:pRepul1");
  memory->create(pRepul2,npairs,nr,"BOP:pRepul2");
  memory->create(pRepul3,npairs,nr,"BOP:pRepul3");
  memory->create(pRepul4,npairs,nr,"BOP:pRepul4");
  memory->create(pRepul5,npairs,nr,"BOP:pRepul5");
  memory->create(pRepul6,npairs,nr,"BOP:pRepul6");
  memory->create(FsigBO,npairs,nBOt,"BOP:FsigBO");
  memory->create(FsigBO1,npairs,nBOt,"BOP:FsigBO1");
  memory->create(FsigBO2,npairs,nBOt,"BOP:FsigBO2");
  memory->create(FsigBO3,npairs,nBOt,"BOP:FsigBO3");
  memory->create(FsigBO4,npairs,nBOt,"BOP:FsigBO4");
  memory->create(FsigBO5,npairs,nBOt,"BOP:FsigBO5");
  memory->create(FsigBO6,npairs,nBOt,"BOP:FsigBO6");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBOP::settings(int narg, char **arg)
{
  table = 0;
  otfly = 1;
  a_flag = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"table") == 0) {
      table = 1;
      iarg++;
    } else if (strcmp(arg[iarg],"save") == 0) {
      otfly = 0;
      iarg++;
    } else if (strcmp(arg[iarg],"sigmaoff") == 0) {
      a_flag = 1;
      iarg++;
    } else error->all(FLERR,"Illegal pair_style command");
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs(Updated: D.K. Ward 05/06/10)
------------------------------------------------------------------------- */

void PairBOP::coeff(int narg, char **arg)
{
  int i,j,n;
  MPI_Comm_rank(world,&me);
  map = new int[atom->ntypes+1];

  if (narg < 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // ensure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read the potential file

  nr=2000;
  nBOt=2000;
  bop_step=0;
  nb_pi=0;
  nb_sg=0;
  allocate_sigma=0;
  allocate_pi=0;
  allocate_neigh=0;
  update_list=0;

  if (table == 0) read_file(arg[2]);
  else read_table(arg[2]);

  if (table == 0) {
    setPbetaS();
    setPbetaP();
    setPrepul();
    setSign();
  }

  // match element names to BOP word types

  if (me == 0) {
    for (i = 3; i < narg; i++) {
      if (strcmp(arg[i],"NULL") == 0) {
        map[i-2] = -1;
        continue;
      }
      for (j = 0; j < bop_types; j++)
        if (strcmp(arg[i],words[j]) == 0) break;
      map[i-2] = j;
    }
  }

  MPI_Bcast(&map[1],atom->ntypes,MPI_INT,0,world);

  if (me == 0) {
    if (words) {
      for (i = 0; i < bop_types; i++) delete [] words[i];
      delete [] words;
    }
  }

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

void PairBOP::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style BOP requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style BOP requires newton pair on");

  // check that user sets comm->cutghostuser to 3x the max BOP cutoff

  if (comm->cutghostuser < 3.0*cutmax - EPSILON) {
    char str[128];
    sprintf(str,"Pair style bop requires comm ghost cutoff "
            "at least 3x larger than %g",cutmax);
    error->all(FLERR,str);
  }

  // need a full neighbor list and neighbors of ghosts

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->ghost = 1;
}

/* ---------------------------------------------------------------------- */

double PairBOP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  int ii = map[i]+1;
  int jj = map[j]+1;

  int ij;
  if (ii==jj) ij=ii-1;
  else if (ii<jj) ij=ii*bop_types-ii*(ii+1)/2+jj-1;
  else ij=jj*bop_types-jj*(jj+1)/2+ii-1;

  cutghost[i][j] = rcut[ij];
  cutghost[j][i] = cutghost[i][j];
  cutsq[i][j] = rcut[ij]*rcut[ij];
  cutsq[j][i] = cutsq[i][j];
  return rcut[ij];
}

/* ----------------------------------------------------------------------
   create BOP neighbor list from main neighbor list
   BOP neighbor list stores neighbors of ghost atoms
   BOP requires neighbor's of k if k is a neighbor of
   j and j is a neighbor of i
------------------------------------------------------------------------- */

void PairBOP::gneigh()
{
  int i,ii;
  int *ilist,*numneigh;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if(allocate_neigh==0) {
    memory->create (BOP_index,nall,"BOP_index");
    if (otfly==0) memory->create (cos_index,nall,"cos_index");
    allocate_neigh=1;
  }
  else {
    memory->grow (BOP_index,nall,"BOP_index");
    if (otfly==0) memory->grow (cos_index,nall,"cos_index");
    allocate_neigh=1;
  }
  ilist = list->ilist;
  numneigh = list->numneigh;
  if(bop_step==0) {
    maxneigh=0;
    maxnall=0;
  }
  neigh_total=0;
  cos_total=0;
  for (ii = 0; ii < nall; ii++) {
    if (ii < nlocal) {
      i=ilist[ii];
      if(numneigh[i]>maxneigh) maxneigh=numneigh[i];
    } else {
      i=ii;
      if(numneigh[i]>maxneigh) maxneigh=numneigh[i];
    }
    BOP_index[i]=neigh_total;
    neigh_total+=numneigh[i];
    if(otfly==0) {
      cos_index[i]=cos_total;
      cos_total+=numneigh[i]*(numneigh[i]-1)/2;
    }
  }
  maxnall=nall;
}

/* ---------------------------------------------------------------------- */

void PairBOP::theta()
{
  int i,j,ii,jj,kk;
  int itype,jtype,i12;
  int temp_ij,temp_ik,temp_ijk;
  int n,nlocal,nall,ks;
  int *ilist,*numneigh;
  int *iilist;
  int **firstneigh;
  double rj2,rk2,rsq,ps;
  double rj1k1,rj2k2;
  double **x = atom->x;
  int *type = atom->type;

  nlocal = atom->nlocal;
  nall = nlocal+atom->nghost;
  ilist = list->ilist;
  firstneigh = list->firstneigh;
  numneigh = list->numneigh;
  if(update_list!=0)
    memory_theta_grow();
  else
    memory_theta_create();
  for (ii = 0; ii < nall; ii++) {
    if(ii<nlocal)
      i= ilist[ii];
    else
      i=ii;
    itype = map[type[i]]+1;

    iilist=firstneigh[i];
    for(jj=0;jj<numneigh[i];jj++) {
      j=iilist[jj];
      temp_ij=BOP_index[i]+jj;
      jtype = map[type[j]]+1;

      if(itype==jtype)
        i12=itype-1;
      else if(itype<jtype)
        i12=itype*bop_types-itype*(itype+1)/2+jtype-1;
      else
        i12=jtype*bop_types-jtype*(jtype+1)/2+itype-1;
      if(i12>=npairs) {
        error->one(FLERR,"Too many atom pairs for pair bop");
      }
      disij[0][temp_ij]=x[j][0]-x[i][0];
      disij[1][temp_ij]=x[j][1]-x[i][1];
      disij[2][temp_ij]=x[j][2]-x[i][2];
      rsq=disij[0][temp_ij]*disij[0][temp_ij]
          +disij[1][temp_ij]*disij[1][temp_ij]
          +disij[2][temp_ij]*disij[2][temp_ij];
      rij[temp_ij]=sqrt(rsq);
      if(rij[temp_ij]<=rcut[i12])
        neigh_flag[temp_ij]=1;
      else
        neigh_flag[temp_ij]=0;
      ps=rij[temp_ij]*rdr[i12]+1.0;
      ks=(int)ps;

      if(nr-1<ks)
        ks=nr-1;
      ps=ps-ks;
      if(ps>1.0)
        ps=1.0;
      betaS[temp_ij]=((pBetaS3[i12][ks-1]*ps+pBetaS2[i12][ks-1])*ps
          +pBetaS1[i12][ks-1])*ps+pBetaS[i12][ks-1];
      dBetaS[temp_ij]=(pBetaS6[i12][ks-1]*ps+pBetaS5[i12][ks-1])*ps
          +pBetaS4[i12][ks-1];
      betaP[temp_ij]=((pBetaP3[i12][ks-1]*ps+pBetaP2[i12][ks-1])*ps
          +pBetaP1[i12][ks-1])*ps+pBetaP[i12][ks-1];
      dBetaP[temp_ij]=(pBetaP6[i12][ks-1]*ps+pBetaP5[i12][ks-1])*ps
          +pBetaP4[i12][ks-1];
      repul[temp_ij]=((pRepul3[i12][ks-1]*ps+pRepul2[i12][ks-1])*ps
          +pRepul1[i12][ks-1])*ps+pRepul[i12][ks-1];
      dRepul[temp_ij]=(pRepul6[i12][ks-1]*ps+pRepul5[i12][ks-1])*ps
          +pRepul4[i12][ks-1];
    }
  }
  for (ii = 0; ii < nall; ii++) {
    n=0;
    if(ii<nlocal)
      i= ilist[ii];
    else
      i=ii;
    iilist=firstneigh[i];
    for(jj=0;jj<numneigh[i];jj++) {
      j=iilist[jj];
      temp_ij=BOP_index[i]+jj;
      rj2=rij[temp_ij]*rij[temp_ij];
      for(kk=jj+1;kk<numneigh[i];kk++) {
        if(cos_index[i]+n>=cos_total) {
          error->one(FLERR,"Too many atom triplets for pair bop");
        }
        temp_ik=BOP_index[i]+kk;
        temp_ijk=cos_index[i]+n;
        if(temp_ijk>=cos_total) {
          error->one(FLERR,"Too many atom triplets for pair bop");
        }
        rk2=rij[temp_ik]*rij[temp_ik];
        rj1k1=rij[temp_ij]*rij[temp_ik];
        rj2k2=rj1k1*rj1k1;
        if(temp_ijk>=cos_total) {
          error->one(FLERR,"Too many atom triplets for pair bop");
        }
        cosAng[temp_ijk]=(disij[0][temp_ij]*disij[0][temp_ik]+disij[1][temp_ij]
            *disij[1][temp_ik]+disij[2][temp_ij]*disij[2][temp_ik])/rj1k1;
        dcAng[temp_ijk][0][0]=(disij[0][temp_ik]*rj1k1-cosAng[temp_ijk]
              *disij[0][temp_ij]*rk2)/(rj2k2);
        dcAng[temp_ijk][1][0]=(disij[1][temp_ik]*rj1k1-cosAng[temp_ijk]
            *disij[1][temp_ij]*rk2)/(rj2k2);
        dcAng[temp_ijk][2][0]=(disij[2][temp_ik]*rj1k1-cosAng[temp_ijk]
            *disij[2][temp_ij]*rk2)/(rj2k2);
        dcAng[temp_ijk][0][1]=(disij[0][temp_ij]*rj1k1-cosAng[temp_ijk]
            *disij[0][temp_ik]*rj2)/(rj2k2);
        dcAng[temp_ijk][1][1]=(disij[1][temp_ij]*rj1k1-cosAng[temp_ijk]
            *disij[1][temp_ik]*rj2)/(rj2k2);
        dcAng[temp_ijk][2][1]=(disij[2][temp_ij]*rj1k1-cosAng[temp_ijk]
            *disij[2][temp_ik]*rj2)/(rj2k2);
        n++;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairBOP::theta_mod()
{
  if(update_list!=0)
    memory_theta_grow();
  else
    memory_theta_create();
}

/* ---------------------------------------------------------------------- */

/*  The formulation differs slightly to avoid negative square roots
    in the calculation of Sigma^(1/2) of (a) Eq. 6 and (b) Eq. 11 */

void PairBOP::sigmaBo()
{
  int nb_t,new_n_tot;
  int n,i,j,k,kp,m,pp,kkp;
  int iij,ji,ki;
  int itmp,jtmp,ktmp,ltmp,mtmp;
  tagint i_tag,j_tag;
  int ngi,ngj,ngk,nglkp,ngli,nglj,ngl;
  int ngji,ngjk,nikj,ngki,ngkj,ngjkp;
  int ngkpk,ngkpj,ngkkp,nglk;
  int njik,nijk,nikkp,nkp,nijkp;
  int nkikp,njikp,nk0;
  int njkpk,nkjkp,njkkp;
  int jNeik,kNeii,kNeij,kNeikp;
  int kpNeij,kpNeik;
  int new1,new2,nlocal;
  int inum,*ilist,*iilist,*jlist,*klist,*kplist;
  int **firstneigh,*numneigh;
  int temp_ji,temp_ikp,temp_kkp;
  int temp_ij,temp_ik,temp_jkp,temp_kk,temp_jk;
  int ang_ijkp,ang_ikkp,ang_jkpk,ang_kjkp;
  int ang_ijk,ang_ikj,ang_jikp,ang_jkkp;
  int ang_jik,ang_kikp;
  int nb_ij,nb_ik,nb_ikp;
  int nb_jk,nb_jkp,nb_kkp;
  int nsearch;
  int sig_flag,setting,ncmp,ks;
  int itype,jtype,ktype,kptype;
  int bt_i,bt_j,bt_ij;
  int kp_index,same_ikp,same_jkp;
  int same_kkp;
  double AA,BB,CC,DD,EE,EE1,FF;
  double AAC,BBC,CCC,DDC,EEC,FFC,GGC;
  double AACFF,UT,bndtmp,UTcom;
  double amean,gmean0,gmean1,gmean2,ps;
  double gfactor1,gprime1,gsqprime;
  double gfactorsq,gfactor2,gprime2;
  double gfactorsq2,gsqprime2;
  double gfactor3,gprime3,gfactor,rfactor;
  double drfactor,gfactor4,gprime4,agpdpr3;
  double rfactor0,rfactorrt,rfactor1rt,rfactor1;
  double rcm1,rcm2,gcm1,gcm2,gcm3;
  double agpdpr1,agpdpr2,app1,app2,app3,app4;
  double dsigB1,dsigB2;
  double part0,part1,part2,part3,part4;
  double psign,bndtmp0,pp1;
  double bndtmp1,bndtmp2,bndtmp3,bndtmp4,bndtmp5;
  double ftmp[3];
  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int newton_pair = force->newton_pair;
  int *type = atom->type;

  nlocal = atom->nlocal;
  firstneigh = list->firstneigh;
  numneigh = list->numneigh;
  inum = list->inum;
  ilist = list->ilist;
  n=0;

//loop over all local atoms

  if(nb_sg>16) {
    nb_sg=16;
  }
  if(nb_sg==0) {
    nb_sg=(maxneigh)*(maxneigh/2);
  }
  if(allocate_sigma) {
    destroy_sigma();
  }
  create_sigma(nb_sg);
  for(itmp=0;itmp<inum;itmp++) {
    i = ilist[itmp];
    i_tag=tag[i];
    itype = map[type[i]]+1;

//j is loop over all neighbors of i

    for(jtmp=0;jtmp<numneigh[i];jtmp++) {
      temp_ij=BOP_index[i]+jtmp;
      if(neigh_flag[temp_ij]) {
        for(m=0;m<nb_sg;m++) {
          for(pp=0;pp<3;pp++) {
            bt_sg[m].dAA[pp]=0.0;
            bt_sg[m].dBB[pp]=0.0;
            bt_sg[m].dCC[pp]=0.0;
            bt_sg[m].dDD[pp]=0.0;
            bt_sg[m].dEE[pp]=0.0;
            bt_sg[m].dEE1[pp]=0.0;
            bt_sg[m].dFF[pp]=0.0;
            bt_sg[m].dAAC[pp]=0.0;
            bt_sg[m].dBBC[pp]=0.0;
            bt_sg[m].dCCC[pp]=0.0;
            bt_sg[m].dDDC[pp]=0.0;
            bt_sg[m].dEEC[pp]=0.0;
            bt_sg[m].dFFC[pp]=0.0;
            bt_sg[m].dGGC[pp]=0.0;
            bt_sg[m].dUT[pp]=0.0;
            bt_sg[m].dSigB1[pp]=0.0;
            bt_sg[m].dSigB[pp]=0.0;
          }
          bt_sg[m].i=-1;
          bt_sg[m].j=-1;
          bt_sg[m].temp=-1;
        }
        nb_t=0;
        iilist=firstneigh[i];
        j=iilist[jtmp];
        jlist=firstneigh[j];
        for(ki=0;ki<numneigh[j];ki++) {
          if(x[jlist[ki]][0]==x[i][0]) {
            if(x[jlist[ki]][1]==x[i][1]) {
              if(x[jlist[ki]][2]==x[i][2]) {
                break;
              }
            }
          }
        }
        j_tag=tag[j];
        jtype = map[type[j]]+1;
        nb_ij=nb_t;
        nb_t++;
        if(nb_t>nb_sg) {
          new_n_tot=nb_sg+maxneigh;
          grow_sigma(nb_sg,new_n_tot);
          nb_sg=new_n_tot;
        }
        bt_sg[nb_ij].temp=temp_ij;
        bt_sg[nb_ij].i=i;
        bt_sg[nb_ij].j=j;
        if(j_tag>=i_tag) {
          if(itype==jtype)
            iij=itype-1;
          else if(itype<jtype)
            iij=itype*bop_types-itype*(itype+1)/2+jtype-1;
          else
            iij=jtype*bop_types-jtype*(jtype+1)/2+itype-1;
          for(ji=0;ji<numneigh[j];ji++) {
            temp_ji=BOP_index[j]+ji;
            if(x[jlist[ji]][0]==x[i][0]) {
              if(x[jlist[ji]][1]==x[i][1]) {
                if(x[jlist[ji]][2]==x[i][2]) {
                  break;
                }
              }
            }
          }
          nSigBk[n]=0;

//AA-EE1 are the components making up Eq. 30 (a)

          AA=0.0;
          BB=0.0;
          CC=0.0;
          DD=0.0;
          EE=0.0;
          EE1=0.0;

//FF is the Beta_sigma^2 term

          FF=betaS[temp_ij]*betaS[temp_ij];

//agpdpr1 is derivative of FF w.r.t. r_ij

          agpdpr1=2.0*betaS[temp_ij]*dBetaS[temp_ij]/rij[temp_ij];

//dXX derivatives are taken with respect to all pairs contributing to the energy
//nb_ij is derivative w.r.t. ij pair

          bt_sg[nb_ij].dFF[0]=agpdpr1*disij[0][temp_ij];
          bt_sg[nb_ij].dFF[1]=agpdpr1*disij[1][temp_ij];
          bt_sg[nb_ij].dFF[2]=agpdpr1*disij[2][temp_ij];

//k is loop over all neighbors of i again with j neighbor of i

          for(ktmp=0;ktmp<numneigh[i];ktmp++) {
            temp_ik=BOP_index[i]+ktmp;
            if(neigh_flag[temp_ik]) {
              if(ktmp!=jtmp) {
                if(jtmp<ktmp) {
                  njik=jtmp*(2*numneigh[i]-jtmp-1)/2+(ktmp-jtmp)-1;
                  ngj=0;
                  ngk=1;
                }
                else {
                  njik=ktmp*(2*numneigh[i]-ktmp-1)/2+(jtmp-ktmp)-1;
                  ngj=1;
                  ngk=0;
                }
                k=iilist[ktmp];
                ktype = map[type[k]]+1;

//find neighbor of k that is equal to i

                klist=firstneigh[k];
                for(kNeii=0;kNeii<numneigh[k];kNeii++) {
                  if(x[klist[kNeii]][0]==x[i][0]) {
                    if(x[klist[kNeii]][1]==x[i][1]) {
                      if(x[klist[kNeii]][2]==x[i][2]) {
                        break;
                      }
                    }
                  }
                }

//find neighbor of i that is equal to k

                for(jNeik=0;jNeik<numneigh[j];jNeik++) {
                  temp_jk=BOP_index[j]+jNeik;
                  if(x[jlist[jNeik]][0]==x[k][0]) {
                    if(x[jlist[jNeik]][1]==x[k][1]) {
                      if(x[jlist[jNeik]][2]==x[k][2]) {
                        break;
                      }
                    }
                  }
                }

//find neighbor of k that is equal to j

                for(kNeij=0;kNeij<numneigh[k];kNeij++) {
                  if(x[klist[kNeij]][0]==x[j][0]) {
                    if(x[klist[kNeij]][1]==x[j][1]) {
                      if(x[klist[kNeij]][2]==x[j][2]) {
                        break;
                      }
                    }
                  }
                }
                sig_flag=0;
                for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                  ncmp=itypeSigBk[n][nsearch];
                  if(x[ncmp][0]==x[k][0]) {
                    if(x[ncmp][1]==x[k][1]) {
                      if(x[ncmp][2]==x[k][2]) {
                        nk0=nsearch;
                        sig_flag=1;
                        break;
                      }
                    }
                  }
                }
                if(sig_flag==0) {
                  nSigBk[n]=nSigBk[n]+1;
                  nk0=nSigBk[n]-1;
                  itypeSigBk[n][nk0]=k;
                }
                nb_ik=nb_t;
                nb_t++;
                if(nb_t>nb_sg) {
                  new_n_tot=nb_sg+maxneigh;
                  grow_sigma(nb_sg,new_n_tot);
                  nb_sg=new_n_tot;
                }
                bt_sg[nb_ik].temp=temp_ik;
                bt_sg[nb_ik].i=i;
                bt_sg[nb_ik].j=k;
                nb_jk=nb_t;
                nb_t++;
                if(nb_t>nb_sg) {
                  new_n_tot=nb_sg+maxneigh;
                  grow_sigma(nb_sg,new_n_tot);
                  nb_sg=new_n_tot;
                }
                bt_sg[nb_jk].temp=temp_jk;
                bt_sg[nb_jk].i=j;
                bt_sg[nb_jk].j=k;
                ang_jik=cos_index[i]+njik;
                gmean0=sigma_g0[jtype-1][itype-1][ktype-1];
                gmean1=sigma_g1[jtype-1][itype-1][ktype-1];
                gmean2=sigma_g2[jtype-1][itype-1][ktype-1];
                amean=cosAng[ang_jik];
                gfactor1=gmean0+gmean1*amean
                    +gmean2*amean*amean;
                gfactorsq=gfactor1*gfactor1;
                gprime1=gmean1+2.0*gmean2*amean;
                gsqprime=2.0*gfactor1*gprime1;

//AA is Eq. 34 (a) or Eq. 10 (c) for the i atom
//1st CC is Eq. 11 (c) for i atom where j & k=neighbor of i

                AA=AA+gfactorsq*betaS[temp_ik]*betaS[temp_ik];
                CC=CC+gfactorsq*betaS[temp_ik]*betaS[temp_ik]*betaS[temp_ik]*betaS[temp_ik];

//agpdpr1 is derivative of AA w.r.t. Beta(rik)
//agpdpr2 is derivative of CC 1st term w.r.t. Beta(rik)
//app1 is derivative of AA w.r.t. cos(theta_jik)
//app2 is derivative of CC 1st term w.r.t. cos(theta_jik)

                agpdpr1=2.0*gfactorsq*betaS[temp_ik]*dBetaS[temp_ik]/rij[temp_ik];
                agpdpr1=2.0*betaS[temp_ik]*betaS[temp_ik]*agpdpr1;
                app1=betaS[temp_ik]*betaS[temp_ik]*gsqprime;
                app1=betaS[temp_ik]*betaS[temp_ik]*app1;
                bt_sg[nb_ij].dAA[0]+=
                    app1*dcAng[ang_jik][0][ngj];
                bt_sg[nb_ij].dAA[1]+=
                    app1*dcAng[ang_jik][1][ngj];
                bt_sg[nb_ij].dAA[2]+=
                    app1*dcAng[ang_jik][2][ngj];
                bt_sg[nb_ij].dCC[0]+=
                    app2*dcAng[ang_jik][0][ngj];
                bt_sg[nb_ij].dCC[1]+=
                    app2*dcAng[ang_jik][1][ngj];
                bt_sg[nb_ij].dCC[2]+=
                    app2*dcAng[ang_jik][2][ngj];
                bt_sg[nb_ik].dAA[0]+=
                    app1*dcAng[ang_jik][0][ngk]
                    +agpdpr1*disij[0][temp_ik];
                bt_sg[nb_ik].dAA[1]+=
                    app1*dcAng[ang_jik][1][ngk]
                    +agpdpr1*disij[1][temp_ik];
                bt_sg[nb_ik].dAA[2]+=
                    app1*dcAng[ang_jik][2][ngk]
                    +agpdpr1*disij[2][temp_ik];
                bt_sg[nb_ik].dCC[0]+=
                    app2*dcAng[ang_jik][0][ngk]
                    +agpdpr2*disij[0][temp_ik];
                bt_sg[nb_ik].dCC[1]+=
                    app2*dcAng[ang_jik][1][ngk]
                    +agpdpr2*disij[1][temp_ik];
                bt_sg[nb_ik].dCC[2]+=
                    app2*dcAng[ang_jik][2][ngk]
                    +agpdpr2*disij[2][temp_ik];

//k' is loop over neighbors all neighbors of j with k a neighbor
//of i and j a neighbor of i and determine which k' is k

                kp_index=0;
                for(ltmp=0;ltmp<numneigh[j];ltmp++) {
                  temp_jkp=BOP_index[j]+ltmp;
                  kp=jlist[ltmp];
                  if(x[kp][0]==x[k][0]) {
                    if(x[kp][1]==x[k][1]) {
                      if(x[kp][2]==x[k][2]) {
                        kp_index=1;
                        break;
                      }
                    }
                  }
                }
                if(kp_index) {

//loop over neighbors of k

                  for(mtmp=0;mtmp<numneigh[k];mtmp++) {
                    kp=klist[mtmp];
                    if(x[kp][0]==x[j][0]) {
                      if(x[kp][1]==x[j][1]) {
                        if(x[kp][2]==x[j][2]) {
                          break;
                        }
                      }
                    }
                  }
                  if(ki<ltmp) {
                    nijk=ki*(2*numneigh[j]-ki-1)/2+(ltmp-ki)-1;
                    ngji=0;
                    ngjk=1;
                  }
                  else {
                    nijk=ltmp*(2*numneigh[j]-ltmp-1)/2+(ki-ltmp)-1;
                    ngji=1;
                    ngjk=0;
                  }
                  if(kNeii<mtmp) {
                    nikj=kNeii*(2*numneigh[k]-kNeii-1)/2+(mtmp-kNeii)-1;
                    ngki=0;
                    ngkj=1;
                  }
                  else {
                    nikj=mtmp*(2*numneigh[k]-mtmp-1)/2+(kNeii-mtmp)-1;
                    ngki=1;
                    ngkj=0;
                  }
                  ang_ijk=cos_index[j]+nijk;
                  gmean0=sigma_g0[itype-1][jtype-1][ktype-1];
                  gmean1=sigma_g1[itype-1][jtype-1][ktype-1];
                  gmean2=sigma_g2[itype-1][jtype-1][ktype-1];
                  amean=cosAng[ang_ijk];
                  gfactor2=gmean0+gmean1*amean
                      +gmean2*amean*amean;
                  gprime2=gmean1+2.0*gmean2*amean;
                  gmean0=sigma_g0[itype-1][ktype-1][jtype-1];
                  gmean1=sigma_g1[itype-1][ktype-1][jtype-1];
                  gmean2=sigma_g2[itype-1][ktype-1][jtype-1];
                  ang_ikj=cos_index[k]+nikj;
                  amean=cosAng[ang_ikj];
                  gfactor3=gmean0+gmean1*amean
                      +gmean2*amean*amean;
                  gprime3=gmean1+2.0*gmean2*amean;
                  gfactor=gfactor1*gfactor2*gfactor3;
                  rfactor=betaS[temp_ik]*betaS[temp_jkp];

//EE1 is (b) Eq. 12

                  EE1=EE1+gfactor*rfactor;

//rcm2 is derivative of EE1 w.r.t Beta(r_jk')
//gcm1 is derivative of EE1 w.r.t cos(theta_jik)
//gcm2 is derivative of EE1 w.r.t cos(theta_ijk)
//gcm3 is derivative of EE1 w.r.t cos(theta_ikj)

                  rcm1=gfactor*betaS[temp_jkp]*dBetaS[temp_ik]/rij[temp_ik];
                  rcm2=gfactor*betaS[temp_ik]*dBetaS[temp_jkp]/rij[temp_jkp];
                  gcm1=rfactor*gprime1*gfactor2*gfactor3;
                  gcm2=rfactor*gfactor1*gprime2*gfactor3;
                  gcm3=rfactor*gfactor1*gfactor2*gprime3;
                  bt_sg[nb_ij].dEE1[0]+=
                      gcm1*dcAng[ang_jik][0][ngj]
                      -gcm2*dcAng[ang_ijk][0][ngji];
                  bt_sg[nb_ij].dEE1[1]+=
                      gcm1*dcAng[ang_jik][1][ngj]
                      -gcm2*dcAng[ang_ijk][1][ngji];
                  bt_sg[nb_ij].dEE1[2]+=
                      gcm1*dcAng[ang_jik][2][ngj]
                      -gcm2*dcAng[ang_ijk][2][ngji];
                  bt_sg[nb_ik].dEE1[0]+=
                      gcm1*dcAng[ang_jik][0][ngk]
                      +rcm1*disij[0][temp_ik]
                      -gcm3*dcAng[ang_ikj][0][ngki];
                  bt_sg[nb_ik].dEE1[1]+=
                      gcm1*dcAng[ang_jik][1][ngk]
                      +rcm1*disij[1][temp_ik]
                      -gcm3*dcAng[ang_ikj][1][ngki];
                  bt_sg[nb_ik].dEE1[2]+=
                      gcm1*dcAng[ang_jik][2][ngk]
                      +rcm1*disij[2][temp_ik]
                      -gcm3*dcAng[ang_ikj][2][ngki];
                  bt_sg[nb_jk].dEE1[0]+=
                      gcm2*dcAng[ang_ijk][0][ngjk]
                      +rcm2*disij[0][temp_jkp]
                      -gcm3*dcAng[ang_ikj][0][ngkj];
                  bt_sg[nb_jk].dEE1[1]+=
                      gcm2*dcAng[ang_ijk][1][ngjk]
                      +rcm2*disij[1][temp_jkp]
                      -gcm3*dcAng[ang_ikj][1][ngkj];
                  bt_sg[nb_jk].dEE1[2]+=
                      gcm2*dcAng[ang_ijk][2][ngjk]
                      +rcm2*disij[2][temp_jkp]
                      -gcm3*dcAng[ang_ikj][2][ngkj];
                }

// k and k' and j are all different neighbors of i

                for(ltmp=0;ltmp<ktmp;ltmp++) {
                  if(ltmp!=jtmp) {
                    temp_ikp=BOP_index[i]+ltmp;
                    if(neigh_flag[temp_ikp]) {
                      kp=iilist[ltmp];
                      kptype = map[type[kp]]+1;
                      for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                        ncmp=itypeSigBk[n][nsearch];
                        if(x[ncmp][0]==x[kp][0]) {
                          if(x[ncmp][1]==x[kp][1]) {
                            if(x[ncmp][2]==x[kp][2]) {
                              break;
                            }
                          }
                        }
                      }
                      if(jtmp<ltmp) {
                        njikp=jtmp*(2*numneigh[i]-jtmp-1)/2+(ltmp-jtmp)-1;
                        nglj=0;
                        ngl=1;
                      }
                      else {
                        njikp=ltmp*(2*numneigh[i]-ltmp-1)/2+(jtmp-ltmp)-1;
                        nglj=1;
                        ngl=0;
                      }
                      if(ktmp<ltmp) {
                        nkikp=ktmp*(2*numneigh[i]-ktmp-1)/2+(ltmp-ktmp)-1;
                        nglk=0;
                        nglkp=1;
                      }
                      else {
                        nkikp=ltmp*(2*numneigh[i]-ltmp-1)/2+(ktmp-ltmp)-1;
                        nglk=1;
                        nglkp=0;
                      }
                      ang_jikp=cos_index[i]+njikp;
                      nb_ikp=nb_t;
                      nb_t++;
                      if(nb_t>nb_sg) {
                        new_n_tot=nb_sg+maxneigh;
                        grow_sigma(nb_sg,new_n_tot);
                        nb_sg=new_n_tot;
                      }
                      bt_sg[nb_ikp].temp=temp_ikp;
                      bt_sg[nb_ikp].i=i;
                      bt_sg[nb_ikp].j=kp;
                      gmean0=sigma_g0[jtype-1][itype-1][kptype-1];
                      gmean1=sigma_g1[jtype-1][itype-1][kptype-1];
                      gmean2=sigma_g2[jtype-1][itype-1][kptype-1];
                      amean=cosAng[ang_jikp];
                      gfactor2=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime2=gmean1+2.0*gmean2*amean;
                      gmean0=sigma_g0[ktype-1][itype-1][kptype-1];
                      gmean1=sigma_g1[ktype-1][itype-1][kptype-1];
                      gmean2=sigma_g2[ktype-1][itype-1][kptype-1];
                      ang_kikp=cos_index[i]+nkikp;
                      amean=cosAng[ang_kikp];
                      gfactor3=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime3=gmean1+2.0*gmean2*amean;
                      gfactor=gfactor1*gfactor2*gfactor3;
                      rfactorrt=betaS[temp_ik]*betaS[temp_ikp];
                      rfactor=rfactorrt*rfactorrt;

//2nd CC is second term of Eq. 11 (c) for i atom where j , k & k' =neighbor of i

                      CC=CC+2.0*gfactor*rfactor;

//agpdpr1 is derivative of CC 2nd term w.r.t. Beta(r_ik)
//agpdpr2 is derivative of CC 2nd term w.r.t. Beta(r_ik')
//app1 is derivative of CC 2nd term w.r.t. cos(theta_jik)
//app2 is derivative of CC 2nd term w.r.t. cos(theta_jik')
//app3 is derivative of CC 2nd term w.r.t. cos(theta_kik')

                      agpdpr1=4.0*gfactor*rfactorrt*betaS[temp_ikp]
                          *dBetaS[temp_ik]/rij[temp_ik];
                      agpdpr2=4.0*gfactor*rfactorrt*betaS[temp_ik]
                          *dBetaS[temp_ikp]/rij[temp_ikp];
                      app1=2.0*rfactor*gfactor2*gfactor3*gprime1;
                      app2=2.0*rfactor*gfactor1*gfactor3*gprime2;
                      app3=2.0*rfactor*gfactor1*gfactor2*gprime3;
                      bt_sg[nb_ij].dCC[0]+=
                          app1*dcAng[ang_jik][0][ngj]
                          +app2*dcAng[ang_jikp][0][nglj];
                      bt_sg[nb_ij].dCC[1]+=
                          app1*dcAng[ang_jik][1][ngj]
                          +app2*dcAng[ang_jikp][1][nglj];
                      bt_sg[nb_ij].dCC[2]+=
                          app1*dcAng[ang_jik][2][ngj]
                          +app2*dcAng[ang_jikp][2][nglj];
                      bt_sg[nb_ik].dCC[0]+=
                          app1*dcAng[ang_jik][0][ngk]
                          +app3*dcAng[ang_kikp][0][nglk]
                          +agpdpr1*disij[0][temp_ik];
                      bt_sg[nb_ik].dCC[1]+=
                          app1*dcAng[ang_jik][1][ngk]
                          +app3*dcAng[ang_kikp][1][nglk]
                          +agpdpr1*disij[1][temp_ik];
                      bt_sg[nb_ik].dCC[2]+=
                          app1*dcAng[ang_jik][2][ngk]
                          +app3*dcAng[ang_kikp][2][nglk]
                          +agpdpr1*disij[2][temp_ik];
                      bt_sg[nb_ikp].dCC[0]+=
                          app2*dcAng[ang_jikp][0][ngl]
                          +app3*dcAng[ang_kikp][0][nglkp]
                          +agpdpr2*disij[0][temp_ikp];
                      bt_sg[nb_ikp].dCC[1]+=
                          app2*dcAng[ang_jikp][1][ngl]
                          +app3*dcAng[ang_kikp][1][nglkp]
                          +agpdpr2*disij[1][temp_ikp];
                      bt_sg[nb_ikp].dCC[2]+=
                          app2*dcAng[ang_jikp][2][ngl]
                          +app3*dcAng[ang_kikp][2][nglkp]
                          +agpdpr2*disij[2][temp_ikp];
                    }
                  }
                }

// j and k are different neighbors of i and k' is a neighbor k not equal to i

                for(ltmp=0;ltmp<numneigh[k];ltmp++) {
                  temp_kkp=BOP_index[k]+ltmp;
                  if(neigh_flag[temp_kkp]) {
                    kp=klist[ltmp];;
                    kptype = map[type[kp]]+1;
                    same_ikp=0;
                    same_jkp=0;
                    if(x[i][0]==x[kp][0]) {
                      if(x[i][1]==x[kp][1]) {
                        if(x[i][2]==x[kp][2]) {
                          same_ikp=1;
                        }
                      }
                    }
                    if(x[j][0]==x[kp][0]) {
                      if(x[j][1]==x[kp][1]) {
                        if(x[j][2]==x[kp][2]) {
                          same_jkp=1;
                        }
                      }
                    }
                    if(!same_ikp&&!same_jkp) {
                      if(kNeii<ltmp) {
                        nikkp=kNeii*(2*numneigh[k]-kNeii-1)/2+(ltmp-kNeii)-1;
                        nglkp=1;
                        ngli=0;
                      }
                      else {
                        nikkp=ltmp*(2*numneigh[k]-ltmp-1)/2+(kNeii-ltmp)-1;
                        nglkp=0;
                        ngli=1;
                      }
                      sig_flag=0;
                      for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                        ncmp=itypeSigBk[n][nsearch];
                        if(x[ncmp][0]==x[kp][0]) {
                          if(x[ncmp][1]==x[kp][1]) {
                            if(x[ncmp][2]==x[kp][2]) {
                              sig_flag=1;
                              nkp=nsearch;
                              break;
                            }
                          }
                        }
                      }
                      if(sig_flag==0) {
                        nSigBk[n]=nSigBk[n]+1;
                        nkp=nSigBk[n]-1;
                        itypeSigBk[n][nkp]=kp;
                      }
                      ang_ikkp=cos_index[k]+nikkp;
                      nb_kkp=nb_t;
                      nb_t++;
                      if(nb_t>nb_sg) {
                        new_n_tot=nb_sg+maxneigh;
                        grow_sigma(nb_sg,new_n_tot);
                        nb_sg=new_n_tot;
                      }
                      bt_sg[nb_kkp].temp=temp_kkp;
                      bt_sg[nb_kkp].i=k;
                      bt_sg[nb_kkp].j=kp;
                      gmean0=sigma_g0[itype-1][ktype-1][kptype-1];
                      gmean1=sigma_g1[itype-1][ktype-1][kptype-1];
                      gmean2=sigma_g2[itype-1][ktype-1][kptype-1];
                      amean=cosAng[ang_ikkp];
                      gfactor2=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime2=gmean1+2.0*gmean2*amean;
                      gfactorsq2=gfactor2*gfactor2;
                      gsqprime2=2.0*gfactor2*gprime2;
                      gfactor=gfactorsq*gfactorsq2;
                      rfactorrt=betaS[temp_ik]*betaS[temp_kkp];
                      rfactor=rfactorrt*rfactorrt;

//3rd CC is third term of Eq. 11 (c) for i atom
//where j , k =neighbor of i & k' =neighbor of k

                      CC=CC+gfactor*rfactor;
                      agpdpr1=2.0*gfactor*rfactorrt*betaS[temp_kkp]
                          *dBetaS[temp_ik]/rij[temp_ik];
                      agpdpr2=2.0*gfactor*rfactorrt*betaS[temp_ik]
                          *dBetaS[temp_kkp]/rij[temp_kkp];
                      app1=rfactor*gfactorsq2*gsqprime;
                      app2=rfactor*gfactorsq*gsqprime2;
                      bt_sg[nb_ij].dCC[0]+=
                          app1*dcAng[ang_jik][0][ngj];
                      bt_sg[nb_ij].dCC[1]+=
                          app1*dcAng[ang_jik][1][ngj];
                      bt_sg[nb_ij].dCC[2]+=
                          app1*dcAng[ang_jik][2][ngj];
                      bt_sg[nb_ik].dCC[0]+=
                          app1*dcAng[ang_jik][0][ngk]
                          +agpdpr1*disij[0][temp_ik]
                          -app2*dcAng[ang_ikkp][0][ngli];
                      bt_sg[nb_ik].dCC[1]+=
                          app1*dcAng[ang_jik][1][ngk]
                          +agpdpr1*disij[1][temp_ik]
                          -app2*dcAng[ang_ikkp][1][ngli];
                      bt_sg[nb_ik].dCC[2]+=
                          app1*dcAng[ang_jik][2][ngk]
                          +agpdpr1*disij[2][temp_ik]
                          -app2*dcAng[ang_ikkp][2][ngli];
                      bt_sg[nb_kkp].dCC[0]+=
                          app2*dcAng[ang_ikkp][0][nglkp]
                          +agpdpr2*disij[0][temp_kkp];
                      bt_sg[nb_kkp].dCC[1]+=
                          app2*dcAng[ang_ikkp][1][nglkp]
                          +agpdpr2*disij[1][temp_kkp];
                      bt_sg[nb_kkp].dCC[2]+=
                          app2*dcAng[ang_ikkp][2][nglkp]
                          +agpdpr2*disij[2][temp_kkp];

                    }
                  }
                }

       //j and k are different neighbors of i and k' is a neighbor j not equal to k

                kplist=firstneigh[kp];
                for(ltmp=0;ltmp<numneigh[j];ltmp++) {
                  sig_flag=0;
                  temp_jkp=BOP_index[j]+ltmp;
                  if(neigh_flag[temp_jkp]) {
                    kp=jlist[ltmp];
                    kptype = map[type[kp]]+1;
                    same_jkp=0;
                    same_kkp=0;
                    for(kpNeij=0;kpNeij<numneigh[kp];kpNeij++) {
                      if(x[j][0]==x[kp][0]) {
                        if(x[j][1]==x[kp][1]) {
                          if(x[j][2]==x[kp][2]) {
                            same_jkp=1;
                            break;
                          }
                        }
                      }
                    }
                    for(kpNeik=0;kpNeik<numneigh[kp];kpNeik++) {
                      if(x[k][0]==x[kp][0]) {
                        if(x[k][1]==x[kp][1]) {
                          if(x[k][2]==x[kp][2]) {
                            same_kkp=1;
                            break;
                          }
                        }
                      }
                    }
                    if(!same_kkp&&!same_jkp) {
                      for(kNeikp=0;kNeikp<numneigh[k];kNeikp++) {
                        temp_kkp=BOP_index[k]+kNeikp;
                        kkp=klist[kNeikp];
                        if(x[kkp][0]==x[kp][0]) {
                          if(x[kkp][1]==x[kp][1]) {
                            if(x[kkp][2]==x[kp][2]) {
                              sig_flag=1;
                              break;
                            }
                          }
                        }
                      }
                      if(sig_flag==1) {
                        for(nsearch=0;nsearch<numneigh[kp];nsearch++) {
                          ncmp=kplist[nsearch];
                          if(x[ncmp][0]==x[j][0]) {
                            if(x[ncmp][1]==x[j][1]) {
                              if(x[ncmp][2]==x[j][2]) {
                                kpNeij=nsearch;
                              }
                            }
                          }
                          if(x[ncmp][0]==x[k][0]) {
                            if(x[ncmp][1]==x[k][1]) {
                              if(x[ncmp][2]==x[k][2]) {
                                kpNeik=nsearch;
                              }
                            }
                          }
                        }
                        if(ji<ltmp) {
                          nijkp=(ji)*numneigh[j]-(ji+1)*(ji+2)/2+ltmp;
                          ngji=0;
                          ngjkp=1;
                        }
                        else {
                          nijkp=(ltmp)*numneigh[j]-(ltmp+1)*(ltmp+2)/2+ji;
                          ngji=1;
                          ngjkp=0;
                        }
                        if(kNeii<kNeikp) {
                          nikkp=(kNeii)*numneigh[k]-(kNeii+1)*(kNeii+2)/2+kNeikp;
                          ngki=0;
                          ngkkp=1;
                        }
                        else {
                          nikkp=(kNeikp)*numneigh[k]-(kNeikp+1)*(kNeikp+2)/2+kNeii;
                          ngki=1;
                          ngkkp=0;
                        }
                        if(kpNeij<kpNeik) {
                          njkpk=(kpNeij)*numneigh[kp]-(kpNeij+1)*(kpNeij+2)/2+kpNeik;
                          ngkpj=0;
                          ngkpk=1;
                        }
                        else {
                          njkpk=(kpNeik)*numneigh[kp]-(kpNeik+1)*(kpNeik+2)/2+kpNeij;
                          ngkpj=1;
                          ngkpk=0;
                        }
                        sig_flag=0;
                        for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                          ncmp=itypeSigBk[n][nsearch];
                          if(x[ncmp][0]==x[kp][0]) {
                            if(x[ncmp][1]==x[kp][1]) {
                              if(x[ncmp][2]==x[kp][2]) {
                                nkp=nsearch;
                                sig_flag=1;
                                break;
                              }
                            }
                          }
                        }
                        if(sig_flag==0) {
                          nSigBk[n]=nSigBk[n]+1;
                          nkp=nSigBk[n]-1;
                          itypeSigBk[n][nkp]=kp;
                        }
                        ang_ijkp=cos_index[j]+nijkp;
                        ang_ikkp=cos_index[k]+nikkp;
                        ang_jkpk=cos_index[kp]+njkpk;
                        nb_jkp=nb_t;
                        nb_t++;
                        if(nb_t>nb_sg) {
                          new_n_tot=nb_sg+maxneigh;
                          grow_sigma(nb_sg,new_n_tot);
                          nb_sg=new_n_tot;
                        }
                        bt_sg[nb_jkp].temp=temp_jkp;
                        bt_sg[nb_jkp].i=j;
                        bt_sg[nb_jkp].j=kp;
                        nb_kkp=nb_t;
                        nb_t++;
                        if(nb_t>nb_sg) {
                          new_n_tot=nb_sg+maxneigh;
                          grow_sigma(nb_sg,new_n_tot);
                          nb_sg=new_n_tot;
                        }
                        bt_sg[nb_kkp].temp=temp_kkp;
                        bt_sg[nb_kkp].i=k;
                        bt_sg[nb_kkp].j=kp;
                        gmean0=sigma_g0[itype-1][jtype-1][kptype-1];
                        gmean1=sigma_g1[itype-1][jtype-1][kptype-1];
                        gmean2=sigma_g2[itype-1][jtype-1][kptype-1];
                        amean=cosAng[ang_ijkp];
                        gfactor2=gmean0+gmean1*amean
                            +gmean2*amean*amean;
                        gprime2=gmean1+2.0*gmean2*amean;
                        gmean0=sigma_g0[itype-1][ktype-1][kptype-1];
                        gmean1=sigma_g1[itype-1][ktype-1][kptype-1];
                        gmean2=sigma_g2[itype-1][ktype-1][kptype-1];
                        amean=cosAng[ang_ikkp];
                        gfactor3=gmean0+gmean1*amean
                            +gmean2*amean*amean;
                        gprime3=gmean1+2.0*gmean2*amean;
                        gmean0=sigma_g0[jtype-1][kptype-1][ktype-1];
                        gmean1=sigma_g1[jtype-1][kptype-1][ktype-1];
                        gmean2=sigma_g2[jtype-1][kptype-1][ktype-1];
                        amean=cosAng[ang_jkpk];
                        gfactor4=gmean0+gmean1*amean
                            +gmean2*amean*amean;
                        gprime4=gmean1+2.0*gmean2*amean;
                        gfactor=gfactor1*gfactor2*gfactor3*gfactor4;
                        rfactor0=(betaS[temp_ik]+small2)*(betaS[temp_jkp]+small2)
                            *(betaS[temp_kkp]+small2);
                        rfactor=pow(rfactor0,2.0/3.0);
                        drfactor=2.0/3.0*pow(rfactor0,-1.0/3.0);

//EE is Eq. 25(notes)

                        EE=EE+gfactor*rfactor;

//agpdpr1 is derivative of agpdpr1 w.r.t. Beta(r_ik)
//agpdpr2 is derivative of agpdpr1 w.r.t. Beta(r_jk')
//agpdpr3 is derivative of agpdpr1 w.r.t. Beta(r_kk')
//app1 is derivative of agpdpr1 w.r.t. cos(theta_jik)
//app2 is derivative of agpdpr1 w.r.t. cos(theta_ijk')
//app3 is derivative of agpdpr1 w.r.t. cos(theta_ikk')
//app4 is derivative of agpdpr1 w.r.t. cos(theta_jk'k)

                        agpdpr1=gfactor*drfactor*(betaS[temp_jkp]+small2)*(betaS[temp_kkp]
                            +small2)*dBetaS[temp_ik]/rij[temp_ik];
                        agpdpr2=gfactor*drfactor*(betaS[temp_ik]+small2)*(betaS[temp_kkp]
                            +small2)*dBetaS[temp_jkp]/rij[temp_jkp];
                        agpdpr3=gfactor*drfactor*(betaS[temp_ik]+small2)*(betaS[temp_jkp]
                            +small2)*dBetaS[temp_kkp]/rij[temp_kkp];
                        app1=rfactor*gfactor2*gfactor3*gfactor4*gprime1;
                        app2=rfactor*gfactor1*gfactor3*gfactor4*gprime2;
                        app3=rfactor*gfactor1*gfactor2*gfactor4*gprime3;
                        app4=rfactor*gfactor1*gfactor2*gfactor3*gprime4;
                        bt_sg[nb_ij].dEE[0]+=
                            app1*dcAng[ang_jik][0][ngj]
                            -app2*dcAng[ang_ijkp][0][ngji];
                        bt_sg[nb_ij].dEE[1]+=
                            app1*dcAng[ang_jik][1][ngj]
                            -app2*dcAng[ang_ijkp][1][ngji];
                        bt_sg[nb_ij].dEE[2]+=
                            app1*dcAng[ang_jik][2][ngj]
                            -app2*dcAng[ang_ijkp][2][ngji];
                        bt_sg[nb_ik].dEE[0]+=
                            app1*dcAng[ang_jik][0][ngk]
                            +agpdpr1*disij[0][temp_ik]
                            -app3*dcAng[ang_ikkp][0][ngki];
                        bt_sg[nb_ik].dEE[1]+=
                            app1*dcAng[ang_jik][1][ngk]
                            +agpdpr1*disij[1][temp_ik]
                            -app3*dcAng[ang_ikkp][1][ngki];
                        bt_sg[nb_ik].dEE[2]+=
                            app1*dcAng[ang_jik][2][ngk]
                            +agpdpr1*disij[2][temp_ik]
                            -app3*dcAng[ang_ikkp][2][ngki];
                        bt_sg[nb_jkp].dEE[0]+=
                            app2*dcAng[ang_ijkp][0][ngjkp]
                            +agpdpr2*disij[0][temp_jkp]
                            -app4*dcAng[ang_jkpk][0][ngkpj];
                        bt_sg[nb_jkp].dEE[1]+=
                            app2*dcAng[ang_ijkp][1][ngjkp]
                            +agpdpr2*disij[1][temp_jkp]
                            -app4*dcAng[ang_jkpk][1][ngkpj];
                        bt_sg[nb_jkp].dEE[2]+=
                            app2*dcAng[ang_ijkp][2][ngjkp]
                            +agpdpr2*disij[2][temp_jkp]
                            -app4*dcAng[ang_jkpk][2][ngkpj];
                        bt_sg[nb_kkp].dEE[0]+=
                            app3*dcAng[ang_ikkp][0][ngkkp]
                            +agpdpr3*disij[0][temp_kkp]
                            -app4*dcAng[ang_jkpk][0][ngkpk];
                        bt_sg[nb_kkp].dEE[1]+=
                            app3*dcAng[ang_ikkp][1][ngkkp]
                            +agpdpr3*disij[1][temp_kkp]
                            -app4*dcAng[ang_jkpk][1][ngkpk];
                        bt_sg[nb_kkp].dEE[2]+=
                            app3*dcAng[ang_ikkp][2][ngkkp]
                            +agpdpr3*disij[2][temp_kkp]
                            -app4*dcAng[ang_jkpk][2][ngkpk];
                      }
                    }
                  }
                }
              }
            }
          }

//j is a neighbor of i and k is a neighbor of j not equal to i

          for(ktmp=0;ktmp<numneigh[j];ktmp++) {
            if(ktmp!=ji) {
              if(ktmp<ji) {
                njik=ktmp*(2*numneigh[j]-ktmp-1)/2+(ji-ktmp)-1;
                ngi=1;
                ngk=0;
              }
              else {
                njik=ji*(2*numneigh[j]-ji-1)/2+(ktmp-ji)-1;
                ngi=0;
                ngk=1;
              }
              temp_jk=BOP_index[j]+ktmp;
              if(neigh_flag[temp_jk]) {
                k=jlist[ktmp];
                ktype=map[type[k]]+1;
                klist=firstneigh[k];

                for(kNeij=0;kNeij<numneigh[k];kNeij++) {
                  if(x[klist[kNeij]][0]==x[j][0]) {
                    if(x[klist[kNeij]][1]==x[j][1]) {
                      if(x[klist[kNeij]][2]==x[j][2]) {
                        break;
                      }
                    }
                  }
                }
                sig_flag=0;
                for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                  ncmp=itypeSigBk[n][nsearch];
                  if(x[ncmp][0]==x[k][0]) {
                    if(x[ncmp][1]==x[k][1]) {
                      if(x[ncmp][2]==x[k][2]) {
                        new1=nsearch;
                        sig_flag=1;
                        break;
                      }
                    }
                  }
                }
                if(sig_flag==0) {
                  nSigBk[n]=nSigBk[n]+1;
                  new1=nSigBk[n]-1;
                  itypeSigBk[n][new1]=k;
                }
                ang_ijk=cos_index[j]+njik;
                nb_jk=nb_t;
                nb_t++;
                if(nb_t>nb_sg) {
                  new_n_tot=nb_sg+maxneigh;
                  grow_sigma(nb_sg,new_n_tot);
                  nb_sg=new_n_tot;
                }
                bt_sg[nb_jk].temp=temp_jk;
                bt_sg[nb_jk].i=j;
                bt_sg[nb_jk].j=k;
                gmean0=sigma_g0[itype-1][jtype-1][ktype-1];
                gmean1=sigma_g1[itype-1][jtype-1][ktype-1];
                gmean2=sigma_g2[itype-1][jtype-1][ktype-1];
                amean=cosAng[ang_ijk];
                gfactor1=gmean0+gmean1*amean
                    +gmean2*amean*amean;
                gprime1=gmean1+2.0*gmean2*amean;
                gfactorsq=gfactor1*gfactor1;
                gsqprime=2.0*gfactor1*gprime1;
                rfactor1rt=betaS[temp_jk]*betaS[temp_jk];
                rfactor1=rfactor1rt*rfactor1rt;

//BB is Eq. 34 (a) or Eq. 10 (c) for the j atom
//1st DD is Eq. 11 (c) for j atom where i & k=neighbor of j

                BB=BB+gfactorsq*rfactor1rt;
                DD=DD+gfactorsq*rfactor1;

//agpdpr1 is derivative of BB  w.r.t. Beta(r_jk)
//app1 is derivative of BB w.r.t. cos(theta_ijk)

                agpdpr1=2.0*gfactorsq*betaS[temp_jk]*dBetaS[temp_jk]/rij[temp_jk];
                app1=rfactor1rt*gsqprime;
                bt_sg[nb_ij].dBB[0]-=
                    app1*dcAng[ang_ijk][0][ngi];
                bt_sg[nb_ij].dBB[1]-=
                    app1*dcAng[ang_ijk][1][ngi];
                bt_sg[nb_ij].dBB[2]-=
                    app1*dcAng[ang_ijk][2][ngi];
                bt_sg[nb_ij].dDD[0]-=
                    app2*dcAng[ang_ijk][0][ngi];
                bt_sg[nb_ij].dDD[1]-=
                    app2*dcAng[ang_ijk][1][ngi];
                bt_sg[nb_ij].dDD[2]-=
                    app2*dcAng[ang_ijk][2][ngi];
                bt_sg[nb_jk].dBB[0]+=
                    app1*dcAng[ang_ijk][0][ngk]
                    +agpdpr1*disij[0][temp_jk];
                bt_sg[nb_jk].dBB[1]+=
                    app1*dcAng[ang_ijk][1][ngk]
                    +agpdpr1*disij[1][temp_jk];
                bt_sg[nb_jk].dBB[2]+=
                    app1*dcAng[ang_ijk][2][ngk]
                    +agpdpr1*disij[2][temp_jk];
                bt_sg[nb_jk].dDD[0]+=
                    app2*dcAng[ang_ijk][0][ngk]
                    +agpdpr2*disij[0][temp_jk];
                bt_sg[nb_jk].dDD[1]+=
                    app2*dcAng[ang_ijk][1][ngk]
                    +agpdpr2*disij[1][temp_jk];
                bt_sg[nb_jk].dDD[2]+=
                    app2*dcAng[ang_ijk][2][ngk]
                    +agpdpr2*disij[2][temp_jk];

//j is a neighbor of i, k and k' prime different neighbors of j not equal to i

                for(ltmp=0;ltmp<ktmp;ltmp++) {
                  if(ltmp!=ji) {
                    temp_jkp=BOP_index[j]+ltmp;
                    if(neigh_flag[temp_jkp]) {
                      kp=jlist[ltmp];
                      kptype=map[type[kp]]+1;
                      for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                        ncmp=itypeSigBk[n][nsearch];
                        if(x[ncmp][0]==x[kp][0]) {
                          if(x[ncmp][1]==x[kp][1]) {
                            if(x[ncmp][2]==x[kp][2]) {
                              new2=nsearch;
                              break;
                            }
                          }
                        }
                      }
                      if(ji<ltmp) {
                        nijkp=ji*(2*numneigh[j]-ji-1)/2+(ltmp-ji)-1;
                        ngli=0;
                        ngl=1;
                      }
                      else {
                        nijkp=ltmp*(2*numneigh[j]-ltmp-1)/2+(ji-ltmp)-1;
                        ngli=1;
                        ngl=0;
                      }
                      if(ktmp<ltmp) {
                        nkjkp=ktmp*(2*numneigh[j]-ktmp-1)/2+(ltmp-ktmp)-1;
                        ngjk=0;
                        ngjkp=1;
                      }
                      else {
                        nkjkp=ltmp*(2*numneigh[j]-ltmp-1)/2+(ktmp-ltmp)-1;
                        ngjk=1;
                        ngjkp=0;
                      }
                      ang_ijkp=cos_index[j]+nijkp;
                      ang_kjkp=cos_index[j]+nkjkp;
                      gmean0=sigma_g0[itype-1][jtype-1][kptype-1];
                      gmean1=sigma_g1[itype-1][jtype-1][kptype-1];
                      gmean2=sigma_g2[itype-1][jtype-1][kptype-1];
                      amean=cosAng[ang_ijkp];
                      gfactor2=gmean0+gmean1*amean
                        +gmean2*amean*amean;
                      gprime2=gmean1+2.0*gmean2*amean;
                      gmean0=sigma_g0[ktype-1][jtype-1][kptype-1];
                      gmean1=sigma_g1[ktype-1][jtype-1][kptype-1];
                      gmean2=sigma_g2[ktype-1][jtype-1][kptype-1];
                      amean=cosAng[ang_kjkp];
                      gfactor3=gmean0+gmean1*amean
                        +gmean2*amean*amean;
                      gprime3=gmean1+2.0*gmean2*amean;
                      gfactor=gfactor1*gfactor2*gfactor3;
                      rfactorrt=betaS[temp_jk]*betaS[temp_jkp];
                      rfactor=rfactorrt*rfactorrt;

//2nd DD is Eq. 11 (c) for j atom where i , k & k'=neighbor of j

                      DD=DD+2.0*gfactor*rfactor;

//agpdpr1 is derivative of DD  w.r.t. Beta(r_jk)
//agpdpr2 is derivative of DD  w.r.t. Beta(r_jk')
//app1 is derivative of DD  w.r.t. cos(theta_ijk)
//app2 is derivative of DD  w.r.t. cos(theta_ijkp)
//app3 is derivative of DD  w.r.t. cos(theta_kjkp)

                      agpdpr1=4.0*gfactor*rfactorrt*betaS[temp_jkp]
                          *dBetaS[temp_jk]/rij[temp_jk];
                          agpdpr2=4.0*gfactor*rfactorrt*betaS[temp_jk]
                          *dBetaS[temp_jkp]/rij[temp_jkp];
                      app1=2.0*rfactor*gfactor2*gfactor3*gprime1;
                      app2=2.0*rfactor*gfactor1*gfactor3*gprime2;
                      app3=2.0*rfactor*gfactor1*gfactor2*gprime3;
                      bt_sg[nb_ij].dDD[0]-=
                          app1*dcAng[ang_ijk][0][ngi]
                          +app2*dcAng[ang_ijkp][0][ngli];
                      bt_sg[nb_ij].dDD[1]-=
                          app1*dcAng[ang_ijk][1][ngi]
                          +app2*dcAng[ang_ijkp][1][ngli];
                      bt_sg[nb_ij].dDD[2]-=
                          app1*dcAng[ang_ijk][2][ngi]
                          +app2*dcAng[ang_ijkp][2][ngli];
                      bt_sg[nb_jk].dDD[0]+=
                          app1*dcAng[ang_ijk][0][ngk]
                          +app3*dcAng[ang_kjkp][0][ngjk]
                          +agpdpr1*disij[0][temp_jk];
                      bt_sg[nb_jk].dDD[1]+=
                          app1*dcAng[ang_ijk][1][ngk]
                          +app3*dcAng[ang_kjkp][1][ngjk]
                          +agpdpr1*disij[1][temp_jk];
                      bt_sg[nb_jk].dDD[2]+=
                          app1*dcAng[ang_ijk][2][ngk]
                          +app3*dcAng[ang_kjkp][2][ngjk]
                          +agpdpr1*disij[2][temp_jk];
                      bt_sg[nb_jkp].dDD[0]+=
                          app2*dcAng[ang_ijkp][0][ngl]
                          +app3*dcAng[ang_kjkp][0][ngjkp]
                          +agpdpr2*disij[0][temp_jkp];
                      bt_sg[nb_jkp].dDD[1]+=
                          app2*dcAng[ang_ijkp][1][ngl]
                          +app3*dcAng[ang_kjkp][1][ngjkp]
                          +agpdpr2*disij[1][temp_jkp];
                      bt_sg[nb_jkp].dDD[2]+=
                          app2*dcAng[ang_ijkp][2][ngl]
                          +app3*dcAng[ang_kjkp][2][ngjkp]
                          +agpdpr2*disij[2][temp_jkp];
                    }
                  }
                }

//j is a neighbor of i, k is a neighbor of j not equal to i and k'
//is a neighbor of k not equal to j or i

                for(ltmp=0;ltmp<numneigh[k];ltmp++) {
                  temp_kkp=BOP_index[k]+ltmp;
                  if(neigh_flag[temp_kkp]) {
                    kp=klist[ltmp];
                    kptype=map[type[kp]]+1;
                    same_ikp=0;
                    same_jkp=0;
                    if(x[i][0]==x[kp][0]) {
                      if(x[i][1]==x[kp][1]) {
                        if(x[i][2]==x[kp][2]) {
                          same_ikp=1;
                        }
                      }
                    }
                    if(x[j][0]==x[kp][0]) {
                      if(x[j][1]==x[kp][1]) {
                        if(x[j][2]==x[kp][2]) {
                          same_jkp=1;
                        }
                      }
                    }
                    if(!same_ikp&&!same_jkp) {
                      for(kNeij=0;kNeij<numneigh[k];kNeij++) {
                        if(x[klist[kNeij]][0]==x[j][0]) {
                          if(x[klist[kNeij]][1]==x[j][1]) {
                            if(x[klist[kNeij]][2]==x[j][2]) {
                              break;
                            }
                          }
                        }
                      }
                      if(kNeij<ltmp) {
                        njkkp=kNeij*(2*numneigh[k]-kNeij-1)/2+(ltmp-kNeij)-1;
                        nglkp=1;
                        nglj=0;
                      }
                      else {
                        njkkp=ltmp*(2*numneigh[k]-ltmp-1)/2+(kNeij-ltmp)-1;
                        nglkp=0;
                        nglj=1;
                      }
                      sig_flag=0;
                      for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                        ncmp=itypeSigBk[n][nsearch];
                        if(x[ncmp][0]==x[kp][0]) {
                          if(x[ncmp][1]==x[kp][1]) {
                            if(x[ncmp][2]==x[kp][2]) {
                              new2=nsearch;
                              sig_flag=1;
                              break;
                            }
                          }
                        }
                      }
                      if(sig_flag==0) {
                        nSigBk[n]=nSigBk[n]+1;
                        new2=nSigBk[n]-1;
                        itypeSigBk[n][new2]=kp;
                      }
                      ang_jkkp=cos_index[k]+njkkp;
                      nb_kkp=nb_t;
                      nb_t++;
                      if(nb_t>nb_sg) {
                        new_n_tot=nb_sg+maxneigh;
                        grow_sigma(nb_sg,new_n_tot);
                        nb_sg=new_n_tot;
                      }
                      bt_sg[nb_kkp].temp=temp_kkp;
                      bt_sg[nb_kkp].i=k;
                      bt_sg[nb_kkp].j=kp;
                      gmean0=sigma_g0[jtype-1][ktype-1][kptype-1];
                      gmean1=sigma_g1[jtype-1][ktype-1][kptype-1];
                      gmean2=sigma_g2[jtype-1][ktype-1][kptype-1];
                      amean=cosAng[ang_jkkp];
                      gfactor2=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime2=gmean1+2.0*gmean2*amean;
                      gfactorsq2=gfactor2*gfactor2;
                      gsqprime2=2.0*gfactor2*gprime2;
                      gfactor=gfactorsq*gfactorsq2;
                      rfactorrt=betaS[temp_jk]*betaS[temp_kkp];
                      rfactor=rfactorrt*rfactorrt;

//3rd DD is Eq. 11 (c) for j atom where i & k=neighbor of j & k'=neighbor of k

                      DD=DD+gfactor*rfactor;

//agpdpr1 is derivative of DD  3rd term w.r.t. Beta(r_jk)
//agpdpr2 is derivative of DD  3rd term w.r.t. Beta(r_kk')
//app1 is derivative of DD  3rd term w.r.t. cos(theta_ijk)
//app2 is derivative of DD  3rd term w.r.t. cos(theta_jkkp)

                      agpdpr1=2.0*gfactor*rfactorrt*betaS[temp_kkp]
                          *dBetaS[temp_jk]/rij[temp_jk];
                      agpdpr2=2.0*gfactor*rfactorrt*betaS[temp_jk]
                          *dBetaS[temp_kkp]/rij[temp_kkp];
                      app1=rfactor*gfactorsq2*gsqprime;
                      app2=rfactor*gfactorsq*gsqprime2;
                      bt_sg[nb_ij].dDD[0]-=
                          app1*dcAng[ang_ijk][0][ngi];
                      bt_sg[nb_ij].dDD[1]-=
                          app1*dcAng[ang_ijk][1][ngi];
                      bt_sg[nb_ij].dDD[2]-=
                          app1*dcAng[ang_ijk][2][ngi];
                      bt_sg[nb_jk].dDD[0]+=
                          app1*dcAng[ang_ijk][0][ngk]
                          +agpdpr1*disij[0][temp_jk]
                          -app2*dcAng[ang_jkkp][0][nglj];
                      bt_sg[nb_jk].dDD[1]+=
                          app1*dcAng[ang_ijk][1][ngk]
                          +agpdpr1*disij[1][temp_jk]
                          -app2*dcAng[ang_jkkp][1][nglj];
                      bt_sg[nb_jk].dDD[2]+=
                          app1*dcAng[ang_ijk][2][ngk]
                          +agpdpr1*disij[2][temp_jk]
                          -app2*dcAng[ang_jkkp][2][nglj];
                      bt_sg[nb_kkp].dDD[0]+=
                          app2*dcAng[ang_jkkp][0][nglkp]
                          +agpdpr2*disij[0][temp_kkp];
                      bt_sg[nb_kkp].dDD[1]+=
                          app2*dcAng[ang_jkkp][1][nglkp]
                          +agpdpr2*disij[1][temp_kkp];
                      bt_sg[nb_kkp].dDD[2]+=
                          app2*dcAng[ang_jkkp][2][nglkp]
                          +agpdpr2*disij[2][temp_kkp];
                    }
                  }
                }
              }
            }
          }

          sig_flag=0;
          if(FF<=0.000001) {
            sigB[n]=0.0;
            sig_flag=1;
          }
          if(sig_flag==0) {
            if(AA<0.0)
              AA=0.0;
            if(BB<0.0)
              BB=0.0;
            if(CC<0.0)
              CC=0.0;
            if(DD<0.0)
              DD=0.0;

// AA and BB are the representations of (a) Eq. 34 and (b) Eq. 9
// for atoms i and j respectively

            AAC=AA+BB;
            BBC=AA*BB;
            CCC=AA*AA+BB*BB;
            DDC=CC+DD;

//EEC is a modified form of (a) Eq. 33

            EEC=(DDC-CCC)/(AAC+2.0*small1);
            AACFF=1.0/(AAC+2.0*small1);
            for(m=0;m<nb_t;m++) {
              if((bt_sg[m].i>-1)&&(bt_sg[m].j>-1)) {
                bt_i=bt_sg[m].i;
                bt_j=bt_sg[m].j;
                bt_sg[m].dAAC[0]=bt_sg[m].dAA[0]
                    +bt_sg[m].dBB[0];
                bt_sg[m].dAAC[1]=bt_sg[m].dAA[1]
                    +bt_sg[m].dBB[1];
                bt_sg[m].dAAC[2]=bt_sg[m].dAA[2]
                    +bt_sg[m].dBB[2];
                bt_sg[m].dBBC[0]=bt_sg[m].dAA[0]*BB
                    +AA*bt_sg[m].dBB[0];
                bt_sg[m].dBBC[1]=bt_sg[m].dAA[1]*BB
                    +AA*bt_sg[m].dBB[1];
                bt_sg[m].dBBC[2]=bt_sg[m].dAA[2]*BB
                    +AA*bt_sg[m].dBB[2];
                bt_sg[m].dCCC[0]=2.0*AA*bt_sg[m].dAA[0]
                    +2.0*BB*bt_sg[m].dBB[0];
                bt_sg[m].dCCC[1]=2.0*AA*bt_sg[m].dAA[1]
                    +2.0*BB*bt_sg[m].dBB[1];
                bt_sg[m].dCCC[2]=2.0*AA*bt_sg[m].dAA[2]
                    +2.0*BB*bt_sg[m].dBB[2];
                bt_sg[m].dDDC[0]=bt_sg[m].dCC[0]
                    +bt_sg[m].dDD[0];
                bt_sg[m].dDDC[1]=bt_sg[m].dCC[1]
                    +bt_sg[m].dDD[1];
                bt_sg[m].dDDC[2]=bt_sg[m].dCC[2]
                    +bt_sg[m].dDD[2];
                bt_sg[m].dEEC[0]=(bt_sg[m].dDDC[0]
                    -bt_sg[m].dCCC[0]
                    -EEC*bt_sg[m].dAAC[0])*AACFF;
                bt_sg[m].dEEC[1]=(bt_sg[m].dDDC[1]
                    -bt_sg[m].dCCC[1]
                    -EEC*bt_sg[m].dAAC[1])*AACFF;
                bt_sg[m].dEEC[2]=(bt_sg[m].dDDC[2]
                    -bt_sg[m].dCCC[2]
                    -EEC*bt_sg[m].dAAC[2])*AACFF;
              }
            }
            UT=EEC*FF+BBC+small3[iij];
            UT=1.0/sqrt(UT);

// FFC is slightly modified form of (a) Eq. 31
// GGC is slightly modified form of (a) Eq. 32
// bndtmp is a slightly modified form of (a) Eq. 30 and (b) Eq. 8

            FFC=BBC*UT;
            GGC=EEC*UT;
            bndtmp=(FF+sigma_delta[iij]*sigma_delta[iij])
                +sigma_c[iij]*AAC+small4;
            UTcom=-0.5*UT*UT*UT;
            for(m=0;m<nb_t;m++) {
              if((bt_sg[m].i>-1)&&(bt_sg[m].j>-1)) {
                bt_sg[m].dUT[0]=UTcom*(bt_sg[m].dEEC[0]*FF
                    +EEC*bt_sg[m].dFF[0]+bt_sg[m].dBBC[0]);
                bt_sg[m].dUT[1]=UTcom*(bt_sg[m].dEEC[1]*FF
                    +EEC*bt_sg[m].dFF[1]+bt_sg[m].dBBC[1]);
                bt_sg[m].dUT[2]=UTcom*(bt_sg[m].dEEC[2]*FF
                    +EEC*bt_sg[m].dFF[2]+bt_sg[m].dBBC[2]);
                bt_sg[m].dFFC[0]=bt_sg[m].dBBC[0]*UT
                    +BBC*bt_sg[m].dUT[0];
                bt_sg[m].dFFC[1]=bt_sg[m].dBBC[1]*UT
                    +BBC*bt_sg[m].dUT[1];
                bt_sg[m].dFFC[2]=bt_sg[m].dBBC[2]*UT
                    +BBC*bt_sg[m].dUT[2];
                bt_sg[m].dGGC[0]=bt_sg[m].dEEC[0]*UT
                    +EEC*bt_sg[m].dUT[0];
                bt_sg[m].dGGC[1]=bt_sg[m].dEEC[1]*UT
                    +EEC*bt_sg[m].dUT[1];
                bt_sg[m].dGGC[2]=bt_sg[m].dEEC[2]*UT
                    +EEC*bt_sg[m].dUT[2];
              }
            }
            psign=1.0;
            if(1.0+sigma_a[iij]*GGC<0.0)
              psign=-1.0;
            bndtmp0=1.0/sqrt(bndtmp);
            sigB1[n]=psign*betaS[temp_ij]*(1.0+sigma_a[iij]*GGC)*bndtmp0;
            bndtmp=-0.5*bndtmp0*bndtmp0*bndtmp0;
            bndtmp1=psign*(1.0+sigma_a[iij]*GGC)*bndtmp0+psign*betaS[temp_ij]
                *(1.0+sigma_a[iij]*GGC)*bndtmp*2.0*betaS[temp_ij]*(1.0
                +sigma_a[iij]*GGC)*(1.0+sigma_a[iij]*GGC);
            bndtmp1=bndtmp1*dBetaS[temp_ij]/rij[temp_ij];
            bndtmp2=psign*betaS[temp_ij]*(1.0+sigma_a[iij]*GGC)*bndtmp*sigma_c[iij];
            bndtmp3=psign*betaS[temp_ij]*(1.0+sigma_a[iij]*GGC)
                *bndtmp*sigma_c[iij]*sigma_a[iij];
            bndtmp4=psign*betaS[temp_ij]*(1.0+sigma_a[iij]*GGC)
                *bndtmp*sigma_c[iij]*sigma_a[iij]*(2.0+GGC);
            bndtmp5=sigma_a[iij]*psign*betaS[temp_ij]*bndtmp0
                +psign*betaS[temp_ij]*(1.0+sigma_a[iij]*GGC)*bndtmp
                *(2.0*(FF+sigma_delta[iij]*sigma_delta[iij])*(1.0
                +sigma_a[iij]*GGC)*sigma_a[iij]+sigma_c[iij]*sigma_a[iij]*FFC);
            setting=0;
            for(m=0;m<nb_t;m++) {
              if((bt_sg[m].i>-1)&&(bt_sg[m].j>-1)) {
                temp_kk=bt_sg[m].temp;
                if(temp_kk==temp_ij&&setting==0) {
                  bt_sg[m].dSigB1[0]=bndtmp1*disij[0][temp_ij]
                      +(bndtmp2*bt_sg[m].dAAC[0]
                      +bndtmp3*bt_sg[m].dEE[0]
                      +bndtmp4*bt_sg[m].dFFC[0]
                      +bndtmp5*bt_sg[m].dGGC[0]);
                  bt_sg[m].dSigB1[1]=bndtmp1*disij[1][temp_ij]
                      +(bndtmp2*bt_sg[m].dAAC[1]
                      +bndtmp3*bt_sg[m].dEE[1]
                      +bndtmp4*bt_sg[m].dFFC[1]
                      +bndtmp5*bt_sg[m].dGGC[1]);
                  bt_sg[m].dSigB1[2]=bndtmp1*disij[2][temp_ij]
                      +(bndtmp2*bt_sg[m].dAAC[2]
                      +bndtmp3*bt_sg[m].dEE[2]
                      +bndtmp4*bt_sg[m].dFFC[2]
                      +bndtmp5*bt_sg[m].dGGC[2]);
                  setting=1;
                }
                else if(temp_kk==temp_ji&&setting==0) {
                  bt_sg[m].dSigB1[0]=-bndtmp1*disij[0][temp_ij]
                      +(bndtmp2*bt_sg[m].dAAC[0]
                      +bndtmp3*bt_sg[m].dEE[0]
                      +bndtmp4*bt_sg[m].dFFC[0]
                      +bndtmp5*bt_sg[m].dGGC[0]);
                  bt_sg[m].dSigB1[1]=-bndtmp1*disij[1][temp_ij]
                      +(bndtmp2*bt_sg[m].dAAC[1]
                      +bndtmp3*bt_sg[m].dEE[1]
                      +bndtmp4*bt_sg[m].dFFC[1]
                      +bndtmp5*bt_sg[m].dGGC[1]);
                  bt_sg[m].dSigB1[2]=-bndtmp1*disij[2][temp_ij]
                      +(bndtmp2*bt_sg[m].dAAC[2]
                      +bndtmp3*bt_sg[m].dEE[2]
                      +bndtmp4*bt_sg[m].dFFC[2]
                      +bndtmp5*bt_sg[m].dGGC[2]);
                  setting=1;
                }
                else {
                  bt_sg[m].dSigB1[0]=(bndtmp2*bt_sg[m].dAAC[0]
                      +bndtmp3*bt_sg[m].dEE[0]
                      +bndtmp4*bt_sg[m].dFFC[0]
                      +bndtmp5*bt_sg[m].dGGC[0]);
                  bt_sg[m].dSigB1[1]=(bndtmp2*bt_sg[m].dAAC[1]
                      +bndtmp3*bt_sg[m].dEE[1]
                      +bndtmp4*bt_sg[m].dFFC[1]
                      +bndtmp5*bt_sg[m].dGGC[1]);
                  bt_sg[m].dSigB1[2]=(bndtmp2*bt_sg[m].dAAC[2]
                      +bndtmp3*bt_sg[m].dEE[2]
                      +bndtmp4*bt_sg[m].dFFC[2]
                      +bndtmp5*bt_sg[m].dGGC[2]);
                }
              }
            }

//This loop is to ensure there is not an error for atoms with no neighbors (deposition)

            if(nb_t==0) {
              if(j>i) {
                bt_sg[0].dSigB1[0]=bndtmp1*disij[0][temp_ij];
                bt_sg[0].dSigB1[1]=bndtmp1*disij[1][temp_ij];
                bt_sg[0].dSigB1[2]=bndtmp1*disij[2][temp_ij];
              }
              else {
                bt_sg[0].dSigB1[0]=-bndtmp1*disij[0][temp_ij];
                bt_sg[0].dSigB1[1]=-bndtmp1*disij[1][temp_ij];
                bt_sg[0].dSigB1[2]=-bndtmp1*disij[2][temp_ij];
              }
              for(pp=0;pp<3;pp++) {
                bt_sg[0].dAA[pp]=0.0;
                bt_sg[0].dBB[pp]=0.0;
                bt_sg[0].dCC[pp]=0.0;
                bt_sg[0].dDD[pp]=0.0;
                bt_sg[0].dEE[pp]=0.0;
                bt_sg[0].dEE1[pp]=0.0;
                bt_sg[0].dFF[pp]=0.0;
                bt_sg[0].dAAC[pp]=0.0;
                bt_sg[0].dBBC[pp]=0.0;
                bt_sg[0].dCCC[pp]=0.0;
                bt_sg[0].dDDC[pp]=0.0;
                bt_sg[0].dEEC[pp]=0.0;
                bt_sg[0].dFFC[pp]=0.0;
                bt_sg[0].dGGC[pp]=0.0;
                bt_sg[0].dUT[pp]=0.0;
                bt_sg[0].dSigB1[pp]=0.0;
                bt_sg[0].dSigB[pp]=0.0;
              }
              bt_sg[0].i=i;
              bt_sg[0].j=j;
              bt_sg[0].temp=temp_ij;
              nb_t++;
              if(nb_t>nb_sg) {
                new_n_tot=nb_sg+maxneigh;
                grow_sigma(nb_sg,new_n_tot);
                nb_sg=new_n_tot;
              }
            }
            ps=sigB1[n]*rdBO+1.0;
            ks=(int)ps;
            if(nBOt-1<ks)
              ks=nBOt-1;
            ps=ps-ks;
            if(ps>1.0)
              ps=1.0;
            dsigB1=((FsigBO3[iij][ks-1]*ps+FsigBO2[iij][ks-1])*ps
                +FsigBO1[iij][ks-1])*ps+FsigBO[iij][ks-1];
            dsigB2=(FsigBO6[iij][ks-1]*ps+FsigBO5[iij][ks-1])*ps+FsigBO4[iij][ks-1];
            part0=(FF+0.5*AAC+small5);
            part1=(sigma_f[iij]-0.5)*sigma_k[iij];
            part2=1.0-part1*EE1/part0;
            part3=dsigB1*part1/part0;
            part4=part3/part0*EE1;

// sigB is the final expression for (a) Eq. 6 and (b) Eq. 11

            sigB[n]=dsigB1*part2;
            pp1=2.0*betaS[temp_ij];
            for(m=0;m<nb_t;m++) {
              if((bt_sg[m].i>-1)&&(bt_sg[m].j>-1)) {
                temp_kk=bt_sg[m].temp;
                bt_ij=bt_sg[m].temp;
                bt_i=bt_sg[m].i;
                bt_j=bt_sg[m].j;
                for(pp=0;pp<3;pp++) {
                  bt_sg[m].dSigB[pp]=dsigB2*part2*bt_sg[m].dSigB1[pp]
                      -part3*bt_sg[m].dEE1[pp]
                      +part4*(bt_sg[m].dFF[pp]
                      +0.5*bt_sg[m].dAAC[pp]);
                }
                for(pp=0;pp<3;pp++) {
                  ftmp[pp]=pp1*bt_sg[m].dSigB[pp];
                  f[bt_i][pp]-=ftmp[pp];
                  f[bt_j][pp]+=ftmp[pp];
                }
                if(evflag) {
                  ev_tally_xyz(bt_i,bt_j,nlocal,newton_pair,0.0,0.0,ftmp[0],ftmp[1]
                      ,ftmp[2],disij[0][bt_ij],disij[1][bt_ij],disij[2][bt_ij]);
                }
              }
            }
          }
          n++;
        }
      }
    }
  }
  if(allocate_sigma)
    destroy_sigma();
}

/* ---------------------------------------------------------------------- */

void PairBOP::sigmaBo_noa()
{
  int nb_t,new_n_tot;
  int n,i,j,k,kp,m,pp;
  int iij,ji,ki;
  int itmp,jtmp,ktmp,ltmp,mtmp;
  tagint i_tag,j_tag;
  int ngi,ngj,ngk;
  int ngji,ngjk,nikj,ngki,ngkj;
  int njik,nijk,nikkp,nkp,nijkp;
  int nkikp,njikp,nk0,nkjkp,njkkp;
  int jNeik,kNeii,kNeij;
  int new1,new2,nlocal,nsearch;
  int inum,*ilist,*iilist,*jlist,*klist;
  int **firstneigh,*numneigh;
  int temp_ji,temp_ikp,temp_kkp;
  int temp_ij,temp_ik,temp_jkp,temp_kk,temp_jk;
  int ang_ijkp,ang_ikkp,ang_kjkp;
  int ang_ijk,ang_ikj,ang_jikp,ang_jkkp;
  int ang_jik,ang_kikp;
  int nb_ij,nb_ik,nb_jk;
  int sig_flag,setting,ncmp,ks;
  int itype,jtype,ktype,kptype;
  int bt_i,bt_j,bt_ij;
  int kp_index,same_ikp,same_jkp;
  double AA,BB,CC,DD,EE1,FF;
  double AAC,BBC,CCC,DDC,EEC;
  double UT,bndtmp;
  double amean,gmean0,gmean1,gmean2,ps;
  double gfactor1,gprime1,gsqprime;
  double gfactorsq,gfactor2,gprime2;
  double gfactorsq2;
  double gfactor3,gprime3,gfactor,rfactor;
  double rfactorrt,rfactor1rt,rfactor1;
  double rcm1,rcm2,gcm1,gcm2,gcm3;
  double agpdpr1,app1;
  double dsigB1,dsigB2;
  double part0,part1,part2,part3,part4;
  double psign,bndtmp0,pp1,bndtmp1,bndtmp2;
  double ftmp[3];
  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int newton_pair = force->newton_pair;
  int *type = atom->type;

  nlocal = atom->nlocal;
  firstneigh = list->firstneigh;
  numneigh = list->numneigh;
  inum = list->inum;
  ilist = list->ilist;
  n=0;

//loop over all local atoms

  if(nb_sg>16) {
    nb_sg=16;
  }
  if(nb_sg==0) {
    nb_sg=(maxneigh)*(maxneigh/2);
  }
  if(allocate_sigma) {
    destroy_sigma();
  }
  create_sigma(nb_sg);
  for(itmp=0;itmp<inum;itmp++) {
    i = ilist[itmp];
    i_tag=tag[i];
    itype = map[type[i]]+1;

//j is loop over all neighbors of i

    for(jtmp=0;jtmp<numneigh[i];jtmp++) {
      temp_ij=BOP_index[i]+jtmp;
      if(neigh_flag[temp_ij]) {
        for(m=0;m<nb_sg;m++) {
          for(pp=0;pp<3;pp++) {
            bt_sg[m].dAA[pp]=0.0;
            bt_sg[m].dBB[pp]=0.0;
            bt_sg[m].dEE1[pp]=0.0;
            bt_sg[m].dFF[pp]=0.0;
            bt_sg[m].dAAC[pp]=0.0;
            bt_sg[m].dSigB1[pp]=0.0;
            bt_sg[m].dSigB[pp]=0.0;
          }
          bt_sg[m].i=-1;
          bt_sg[m].j=-1;
        }
        nb_t=0;
        iilist=firstneigh[i];
        j=iilist[jtmp];
        jlist=firstneigh[j];
        for(ki=0;ki<numneigh[j];ki++) {
          if(x[jlist[ki]][0]==x[i][0]) {
            if(x[jlist[ki]][1]==x[i][1]) {
              if(x[jlist[ki]][2]==x[i][2]) {
                break;
              }
            }
          }
        }
        j_tag=tag[j];
        jtype = map[type[j]]+1;
        nb_ij=nb_t;
        nb_t++;
        if(nb_t>nb_sg) {
          new_n_tot=nb_sg+maxneigh;
          grow_sigma(nb_sg,new_n_tot);
          nb_sg=new_n_tot;
        }
        bt_sg[nb_ij].temp=temp_ij;
        bt_sg[nb_ij].i=i;
        bt_sg[nb_ij].j=j;
        if(j_tag>=i_tag) {
          if(itype==jtype)
            iij=itype-1;
          else if(itype<jtype)
            iij=itype*bop_types-itype*(itype+1)/2+jtype-1;
          else
            iij=jtype*bop_types-jtype*(jtype+1)/2+itype-1;
          for(ji=0;ji<numneigh[j];ji++) {
            temp_ji=BOP_index[j]+ji;
            if(x[jlist[ji]][0]==x[i][0]) {
              if(x[jlist[ji]][1]==x[i][1]) {
                if(x[jlist[ji]][2]==x[i][2]) {
                  break;
                }
              }
            }
          }
          nSigBk[n]=0;

//AA-EE1 are the components making up Eq. 30 (a)

          AA=0.0;
          BB=0.0;
          CC=0.0;
          DD=0.0;
          EE1=0.0;

//FF is the Beta_sigma^2 term

          FF=betaS[temp_ij]*betaS[temp_ij];

//agpdpr1 is derivative of FF w.r.t. r_ij

          agpdpr1=2.0*betaS[temp_ij]*dBetaS[temp_ij]/rij[temp_ij];

//dXX derivatives are taken with respect to all pairs contributing to the energy
//nb_ij is derivative w.r.t. ij pair

          bt_sg[nb_ij].dFF[0]=agpdpr1*disij[0][temp_ij];
          bt_sg[nb_ij].dFF[1]=agpdpr1*disij[1][temp_ij];
          bt_sg[nb_ij].dFF[2]=agpdpr1*disij[2][temp_ij];

//k is loop over all neighbors of i again with j neighbor of i
          for(ktmp=0;ktmp<numneigh[i];ktmp++) {
            temp_ik=BOP_index[i]+ktmp;
            if(neigh_flag[temp_ik]) {
              if(ktmp!=jtmp) {
                if(jtmp<ktmp) {
                  njik=jtmp*(2*numneigh[i]-jtmp-1)/2+(ktmp-jtmp)-1;
                  ngj=0;
                  ngk=1;
                }
                else {
                  njik=ktmp*(2*numneigh[i]-ktmp-1)/2+(jtmp-ktmp)-1;
                  ngj=1;
                  ngk=0;
                }
                k=iilist[ktmp];
                ktype = map[type[k]]+1;

//find neighbor of k that is equal to i

                klist=firstneigh[k];
                for(kNeii=0;kNeii<numneigh[k];kNeii++) {
                  if(x[klist[kNeii]][0]==x[i][0]) {
                    if(x[klist[kNeii]][1]==x[i][1]) {
                      if(x[klist[kNeii]][2]==x[i][2]) {
                        break;
                      }
                    }
                  }
                }

//find neighbor of i that is equal to k

                for(jNeik=0;jNeik<numneigh[j];jNeik++) {
                  temp_jk=BOP_index[j]+jNeik;
                  if(x[jlist[jNeik]][0]==x[k][0]) {
                    if(x[jlist[jNeik]][1]==x[k][1]) {
                      if(x[jlist[jNeik]][2]==x[k][2]) {
                        break;
                      }
                    }
                  }
                }

//find neighbor of k that is equal to j

                for(kNeij=0;kNeij<numneigh[k];kNeij++) {
                  if(x[klist[kNeij]][0]==x[j][0]) {
                    if(x[klist[kNeij]][1]==x[j][1]) {
                      if(x[klist[kNeij]][2]==x[j][2]) {
                        break;
                      }
                    }
                  }
                }
                sig_flag=0;
                for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                  ncmp=itypeSigBk[n][nsearch];
                  if(x[ncmp][0]==x[k][0]) {
                    if(x[ncmp][1]==x[k][1]) {
                      if(x[ncmp][2]==x[k][2]) {
                        nk0=nsearch;
                        sig_flag=1;
                        break;
                      }
                    }
                  }
                }
                if(sig_flag==0) {
                  nSigBk[n]=nSigBk[n]+1;
                  nk0=nSigBk[n]-1;
                  itypeSigBk[n][nk0]=k;
                }
                nb_ik=nb_t;
                nb_t++;
                if(nb_t>nb_sg) {
                  new_n_tot=nb_sg+maxneigh;
                  grow_sigma(nb_sg,new_n_tot);
                  nb_sg=new_n_tot;
                }
                bt_sg[nb_ik].temp=temp_ik;
                bt_sg[nb_ik].i=i;
                bt_sg[nb_ik].j=k;
                nb_jk=nb_t;
                nb_t++;
                if(nb_t>nb_sg) {
                  new_n_tot=nb_sg+maxneigh;
                  grow_sigma(nb_sg,new_n_tot);
                  nb_sg=new_n_tot;
                }
                bt_sg[nb_jk].temp=temp_jk;
                bt_sg[nb_jk].i=j;
                bt_sg[nb_jk].j=k;
                ang_jik=cos_index[i]+njik;
                if(ang_jik>=cos_total) {
                  error->one(FLERR,"Too many atom triplets for pair bop");
                }
                gmean0=sigma_g0[jtype-1][itype-1][ktype-1];
                gmean1=sigma_g1[jtype-1][itype-1][ktype-1];
                gmean2=sigma_g2[jtype-1][itype-1][ktype-1];
                amean=cosAng[ang_jik];
                gfactor1=gmean0+gmean1*amean
                    +gmean2*amean*amean;
                gfactorsq=gfactor1*gfactor1;
                gprime1=gmean1+2.0*gmean2*amean;
                gsqprime=2.0*gfactor1*gprime1;

//AA is Eq. 34 (a) or Eq. 10 (c) for the i atom
//1st CC is Eq. 11 (c) for i atom where j & k=neighbor of i

                AA=AA+gfactorsq*betaS[temp_ik]*betaS[temp_ik];
                CC=CC+gfactorsq*betaS[temp_ik]*betaS[temp_ik]*betaS[temp_ik]*betaS[temp_ik];
//agpdpr1 is derivative of AA w.r.t. Beta(rik)
//agpdpr2 is derivative of CC 1st term w.r.t. Beta(rik)
//app1 is derivative of AA w.r.t. cos(theta_jik)
//app2 is derivative of CC 1st term w.r.t. cos(theta_jik)

                agpdpr1=2.0*gfactorsq*betaS[temp_ik]*dBetaS[temp_ik]/rij[temp_ik];
                app1=betaS[temp_ik]*betaS[temp_ik]*gsqprime;
                bt_sg[nb_ij].dAA[0]+=
                    app1*dcAng[ang_jik][0][ngj];
                bt_sg[nb_ij].dAA[1]+=
                    app1*dcAng[ang_jik][1][ngj];
                bt_sg[nb_ij].dAA[2]+=
                    app1*dcAng[ang_jik][2][ngj];
                bt_sg[nb_ik].dAA[0]+=
                    app1*dcAng[ang_jik][0][ngk]
                    +agpdpr1*disij[0][temp_ik];
                bt_sg[nb_ik].dAA[1]+=
                    app1*dcAng[ang_jik][1][ngk]
                    +agpdpr1*disij[1][temp_ik];
                bt_sg[nb_ik].dAA[2]+=
                    app1*dcAng[ang_jik][2][ngk]
                    +agpdpr1*disij[2][temp_ik];

//k' is loop over neighbors all neighbors of j with k a neighbor
//of i and j a neighbor of i and determine which k' is k

                kp_index=0;
                for(ltmp=0;ltmp<numneigh[j];ltmp++) {
                  temp_jkp=BOP_index[j]+ltmp;
                  kp=jlist[ltmp];
                  if(x[kp][0]==x[k][0]) {
                    if(x[kp][1]==x[k][1]) {
                      if(x[kp][2]==x[k][2]) {
                        kp_index=1;
                        break;
                      }
                    }
                  }
                }
                if(kp_index) {

//loop over neighbors of k

                  for(mtmp=0;mtmp<numneigh[k];mtmp++) {
                    kp=klist[mtmp];
                    if(x[kp][0]==x[j][0]) {
                      if(x[kp][1]==x[j][1]) {
                        if(x[kp][2]==x[j][2]) {
                          break;
                        }
                      }
                    }
                  }
                  if(ki<ltmp) {
                    nijk=ki*(2*numneigh[j]-ki-1)/2+(ltmp-ki)-1;
                    ngji=0;
                    ngjk=1;
                  }
                  else {
                    nijk=ltmp*(2*numneigh[j]-ltmp-1)/2+(ki-ltmp)-1;
                    ngji=1;
                    ngjk=0;
                  }
                  if(kNeii<mtmp) {
                    nikj=kNeii*(2*numneigh[k]-kNeii-1)/2+(mtmp-kNeii)-1;
                    ngki=0;
                    ngkj=1;
                  }
                  else {
                    nikj=mtmp*(2*numneigh[k]-mtmp-1)/2+(kNeii-mtmp)-1;
                    ngki=1;
                    ngkj=0;
                  }
                  ang_ijk=cos_index[j]+nijk;
                  if(ang_ijk>=cos_total) {
                    error->one(FLERR,"Too many atom triplets for pair bop");
                  }
                  gmean0=sigma_g0[itype-1][jtype-1][ktype-1];
                  gmean1=sigma_g1[itype-1][jtype-1][ktype-1];
                  gmean2=sigma_g2[itype-1][jtype-1][ktype-1];
                  amean=cosAng[ang_ijk];
                  gfactor2=gmean0+gmean1*amean
                      +gmean2*amean*amean;
                  gprime2=gmean1+2.0*gmean2*amean;
                  gmean0=sigma_g0[itype-1][ktype-1][jtype-1];
                  gmean1=sigma_g1[itype-1][ktype-1][jtype-1];
                  gmean2=sigma_g2[itype-1][ktype-1][jtype-1];
                  ang_ikj=cos_index[k]+nikj;
                  if(ang_ikj>=cos_total) {
                    error->one(FLERR,"Too many atom triplets for pair bop");
                  }
                  amean=cosAng[ang_ikj];
                  gfactor3=gmean0+gmean1*amean
                      +gmean2*amean*amean;
                  gprime3=gmean1+2.0*gmean2*amean;
                  gfactor=gfactor1*gfactor2*gfactor3;
                  rfactor=betaS[temp_ik]*betaS[temp_jkp];

//EE1 is (b) Eq. 12

                  EE1=EE1+gfactor*rfactor;

//rcm2 is derivative of EE1 w.r.t Beta(r_jk')
//gcm1 is derivative of EE1 w.r.t cos(theta_jik)
//gcm2 is derivative of EE1 w.r.t cos(theta_ijk)
//gcm3 is derivative of EE1 w.r.t cos(theta_ikj)

                  rcm1=gfactor*betaS[temp_jkp]*dBetaS[temp_ik]/rij[temp_ik];
                  rcm2=gfactor*betaS[temp_ik]*dBetaS[temp_jkp]/rij[temp_jkp];
                  gcm1=rfactor*gprime1*gfactor2*gfactor3;
                  gcm2=rfactor*gfactor1*gprime2*gfactor3;
                  gcm3=rfactor*gfactor1*gfactor2*gprime3;
                  bt_sg[nb_ij].dEE1[0]+=
                      gcm1*dcAng[ang_jik][0][ngj]
                      -gcm2*dcAng[ang_ijk][0][ngji];
                  bt_sg[nb_ij].dEE1[1]+=
                      gcm1*dcAng[ang_jik][1][ngj]
                      -gcm2*dcAng[ang_ijk][1][ngji];
                  bt_sg[nb_ij].dEE1[2]+=
                      gcm1*dcAng[ang_jik][2][ngj]
                      -gcm2*dcAng[ang_ijk][2][ngji];
                  bt_sg[nb_ik].dEE1[0]+=
                      gcm1*dcAng[ang_jik][0][ngk]
                      +rcm1*disij[0][temp_ik]
                      -gcm3*dcAng[ang_ikj][0][ngki];
                  bt_sg[nb_ik].dEE1[1]+=
                      gcm1*dcAng[ang_jik][1][ngk]
                      +rcm1*disij[1][temp_ik]
                      -gcm3*dcAng[ang_ikj][1][ngki];
                  bt_sg[nb_ik].dEE1[2]+=
                      gcm1*dcAng[ang_jik][2][ngk]
                      +rcm1*disij[2][temp_ik]
                      -gcm3*dcAng[ang_ikj][2][ngki];
                  bt_sg[nb_jk].dEE1[0]+=
                      gcm2*dcAng[ang_ijk][0][ngjk]
                      +rcm2*disij[0][temp_jkp]
                      -gcm3*dcAng[ang_ikj][0][ngkj];
                  bt_sg[nb_jk].dEE1[1]+=
                      gcm2*dcAng[ang_ijk][1][ngjk]
                      +rcm2*disij[1][temp_jkp]
                      -gcm3*dcAng[ang_ikj][1][ngkj];
                  bt_sg[nb_jk].dEE1[2]+=
                      gcm2*dcAng[ang_ijk][2][ngjk]
                      +rcm2*disij[2][temp_jkp]
                      -gcm3*dcAng[ang_ikj][2][ngkj];
                }

// k and k' and j are all different neighbors of i

                for(ltmp=0;ltmp<ktmp;ltmp++) {
                  if(ltmp!=jtmp) {
                    temp_ikp=BOP_index[i]+ltmp;
                    if(neigh_flag[temp_ikp]) {
                      kp=iilist[ltmp];
                      kptype = map[type[kp]]+1;
                      for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                        ncmp=itypeSigBk[n][nsearch];
                        if(x[ncmp][0]==x[kp][0]) {
                          if(x[ncmp][1]==x[kp][1]) {
                            if(x[ncmp][2]==x[kp][2]) {
                              break;
                            }
                          }
                        }
                      }
                      if(jtmp<ltmp) {
                        njikp=jtmp*(2*numneigh[i]-jtmp-1)/2+(ltmp-jtmp)-1;
                      } else {
                        njikp=ltmp*(2*numneigh[i]-ltmp-1)/2+(jtmp-ltmp)-1;
                      }
                      if(ktmp<ltmp) {
                        nkikp=ktmp*(2*numneigh[i]-ktmp-1)/2+(ltmp-ktmp)-1;
                      } else {
                        nkikp=ltmp*(2*numneigh[i]-ltmp-1)/2+(ktmp-ltmp)-1;
                      }
                      ang_jikp=cos_index[i]+njikp;
                      if(ang_jikp>=cos_total) {
                        error->one(FLERR,"Too many atom triplets for pair bop");
                      }
                      gmean0=sigma_g0[jtype-1][itype-1][kptype-1];
                      gmean1=sigma_g1[jtype-1][itype-1][kptype-1];
                      gmean2=sigma_g2[jtype-1][itype-1][kptype-1];
                      amean=cosAng[ang_jikp];
                      gfactor2=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime2=gmean1+2.0*gmean2*amean;
                      gmean0=sigma_g0[ktype-1][itype-1][kptype-1];
                      gmean1=sigma_g1[ktype-1][itype-1][kptype-1];
                      gmean2=sigma_g2[ktype-1][itype-1][kptype-1];
                      ang_kikp=cos_index[i]+nkikp;
                      if(ang_kikp>=cos_total) {
                        error->one(FLERR,"Too many atom triplets for pair bop");
                      }
                      amean=cosAng[ang_kikp];
                      gfactor3=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime3=gmean1+2.0*gmean2*amean;
                      gfactor=gfactor1*gfactor2*gfactor3;
                      rfactorrt=betaS[temp_ik]*betaS[temp_ikp];
                      rfactor=rfactorrt*rfactorrt;

//2nd CC is second term of Eq. 11 (c) for i atom where j , k & k' =neighbor of i

                      CC=CC+2.0*gfactor*rfactor;
                    }
                  }
                }

// j and k are different neighbors of i and k' is a neighbor k not equal to i

                for(ltmp=0;ltmp<numneigh[k];ltmp++) {
                  temp_kkp=BOP_index[k]+ltmp;
                  if(neigh_flag[temp_kkp]) {
                    kp=klist[ltmp];;
                    kptype = map[type[kp]]+1;
                    same_ikp=0;
                    same_jkp=0;
                    if(x[i][0]==x[kp][0]) {
                      if(x[i][1]==x[kp][1]) {
                        if(x[i][2]==x[kp][2]) {
                          same_ikp=1;
                        }
                      }
                    }
                    if(x[j][0]==x[kp][0]) {
                      if(x[j][1]==x[kp][1]) {
                        if(x[j][2]==x[kp][2]) {
                          same_jkp=1;
                        }
                      }
                    }
                    if(!same_ikp&&!same_jkp) {
                      if(kNeii<ltmp) {
                        nikkp=kNeii*(2*numneigh[k]-kNeii-1)/2+(ltmp-kNeii)-1;
                      } else {
                        nikkp=ltmp*(2*numneigh[k]-ltmp-1)/2+(kNeii-ltmp)-1;
                      }
                      sig_flag=0;
                      for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                        ncmp=itypeSigBk[n][nsearch];
                        if(x[ncmp][0]==x[kp][0]) {
                          if(x[ncmp][1]==x[kp][1]) {
                            if(x[ncmp][2]==x[kp][2]) {
                              sig_flag=1;
                              nkp=nsearch;
                              break;
                            }
                          }
                        }
                      }
                      if(sig_flag==0) {
                        nSigBk[n]=nSigBk[n]+1;
                        nkp=nSigBk[n]-1;
                        itypeSigBk[n][nkp]=kp;
                      }
                      ang_ikkp=cos_index[k]+nikkp;
                      if(ang_ikkp>=cos_total) {
                        error->one(FLERR,"Too many atom triplets for pair bop");
                      }
                      gmean0=sigma_g0[itype-1][ktype-1][kptype-1];
                      gmean1=sigma_g1[itype-1][ktype-1][kptype-1];
                      gmean2=sigma_g2[itype-1][ktype-1][kptype-1];
                      amean=cosAng[ang_ikkp];
                      gfactor2=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime2=gmean1+2.0*gmean2*amean;
                      gfactorsq2=gfactor2*gfactor2;
                      gfactor=gfactorsq*gfactorsq2;
                      rfactorrt=betaS[temp_ik]*betaS[temp_kkp];
                      rfactor=rfactorrt*rfactorrt;

//3rd CC is third term of Eq. 11 (c) for i atom
//where j , k =neighbor of i & k' =neighbor of k

                      CC=CC+gfactor*rfactor;
                    }
                  }
                }
              }
            }
          }

//j is a neighbor of i and k is a neighbor of j not equal to i

          for(ktmp=0;ktmp<numneigh[j];ktmp++) {
            if(ktmp!=ji) {
              if(ktmp<ji) {
                njik=ktmp*(2*numneigh[j]-ktmp-1)/2+(ji-ktmp)-1;
                ngi=1;
                ngk=0;
              }
              else {
                njik=ji*(2*numneigh[j]-ji-1)/2+(ktmp-ji)-1;
                ngi=0;
                ngk=1;
              }
              temp_jk=BOP_index[j]+ktmp;
              if(neigh_flag[temp_jk]) {
                k=jlist[ktmp];
                ktype=map[type[k]]+1;
                klist=firstneigh[k];

                for(kNeij=0;kNeij<numneigh[k];kNeij++) {
                  if(x[klist[kNeij]][0]==x[j][0]) {
                    if(x[klist[kNeij]][1]==x[j][1]) {
                      if(x[klist[kNeij]][2]==x[j][2]) {
                        break;
                      }
                    }
                  }
                }
                sig_flag=0;
                for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                  ncmp=itypeSigBk[n][nsearch];
                  if(x[ncmp][0]==x[k][0]) {
                    if(x[ncmp][1]==x[k][1]) {
                      if(x[ncmp][2]==x[k][2]) {
                        new1=nsearch;
                        sig_flag=1;
                        break;
                      }
                    }
                  }
                }
                if(sig_flag==0) {
                  nSigBk[n]=nSigBk[n]+1;
                  new1=nSigBk[n]-1;
                  itypeSigBk[n][new1]=k;
                }
                ang_ijk=cos_index[j]+njik;
                if(ang_ijk>=cos_total) {
                  error->one(FLERR,"Too many atom triplets for pair bop");
                }
                nb_jk=nb_t;
                nb_t++;
                if(nb_t>nb_sg) {
                  new_n_tot=nb_sg+maxneigh;
                  grow_sigma(nb_sg,new_n_tot);
                  nb_sg=new_n_tot;
                }
                bt_sg[nb_jk].temp=temp_jk;
                bt_sg[nb_jk].i=j;
                bt_sg[nb_jk].j=k;
                gmean0=sigma_g0[itype-1][jtype-1][ktype-1];
                gmean1=sigma_g1[itype-1][jtype-1][ktype-1];
                gmean2=sigma_g2[itype-1][jtype-1][ktype-1];
                amean=cosAng[ang_ijk];
                gfactor1=gmean0+gmean1*amean
                    +gmean2*amean*amean;
                gprime1=gmean1+2.0*gmean2*amean;
                gfactorsq=gfactor1*gfactor1;
                gsqprime=2.0*gfactor1*gprime1;
                rfactor1rt=betaS[temp_jk]*betaS[temp_jk];
                rfactor1=rfactor1rt*rfactor1rt;

//BB is Eq. 34 (a) or Eq. 10 (c) for the j atom
//1st DD is Eq. 11 (c) for j atom where i & k=neighbor of j
                BB=BB+gfactorsq*rfactor1rt;
                DD=DD+gfactorsq*rfactor1;

//agpdpr1 is derivative of BB  w.r.t. Beta(r_jk)
//app1 is derivative of BB w.r.t. cos(theta_ijk)

                agpdpr1=2.0*gfactorsq*betaS[temp_jk]*dBetaS[temp_jk]/rij[temp_jk];
                app1=rfactor1rt*gsqprime;
                bt_sg[nb_ij].dBB[0]-=
                    app1*dcAng[ang_ijk][0][ngi];
                bt_sg[nb_ij].dBB[1]-=
                    app1*dcAng[ang_ijk][1][ngi];
                bt_sg[nb_ij].dBB[2]-=
                    app1*dcAng[ang_ijk][2][ngi];
                bt_sg[nb_jk].dBB[0]+=
                    app1*dcAng[ang_ijk][0][ngk]
                    +agpdpr1*disij[0][temp_jk];
                bt_sg[nb_jk].dBB[1]+=
                    app1*dcAng[ang_ijk][1][ngk]
                    +agpdpr1*disij[1][temp_jk];
                bt_sg[nb_jk].dBB[2]+=
                    app1*dcAng[ang_ijk][2][ngk]
                    +agpdpr1*disij[2][temp_jk];

//j is a neighbor of i, k and k' prime different neighbors of j not equal to i

                for(ltmp=0;ltmp<ktmp;ltmp++) {
                  if(ltmp!=ji) {
                    temp_jkp=BOP_index[j]+ltmp;
                    if(neigh_flag[temp_jkp]) {
                      kp=jlist[ltmp];
                      kptype=map[type[kp]]+1;
                      for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                        ncmp=itypeSigBk[n][nsearch];
                        if(x[ncmp][0]==x[kp][0]) {
                          if(x[ncmp][1]==x[kp][1]) {
                            if(x[ncmp][2]==x[kp][2]) {
                              new2=nsearch;
                              break;
                            }
                          }
                        }
                      }
                      if(ji<ltmp) {
                        nijkp=ji*(2*numneigh[j]-ji-1)/2+(ltmp-ji)-1;
                      } else {
                        nijkp=ltmp*(2*numneigh[j]-ltmp-1)/2+(ji-ltmp)-1;
                      }
                      if(ktmp<ltmp) {
                        nkjkp=ktmp*(2*numneigh[j]-ktmp-1)/2+(ltmp-ktmp)-1;
                        ngjk=0;
                      }
                      else {
                        nkjkp=ltmp*(2*numneigh[j]-ltmp-1)/2+(ktmp-ltmp)-1;
                        ngjk=1;
                      }
                      ang_ijkp=cos_index[j]+nijkp;
                      if(ang_ijkp>=cos_total) {
                        error->one(FLERR,"Too many atom triplets for pair bop");
                      }
                      ang_kjkp=cos_index[j]+nkjkp;
                      if(ang_kjkp>=cos_total) {
                        error->one(FLERR,"Too many atom triplets for pair bop");
                      }
                      gmean0=sigma_g0[itype-1][jtype-1][kptype-1];
                      gmean1=sigma_g1[itype-1][jtype-1][kptype-1];
                      gmean2=sigma_g2[itype-1][jtype-1][kptype-1];
                      amean=cosAng[ang_ijkp];
                      gfactor2=gmean0+gmean1*amean
                        +gmean2*amean*amean;
                      gprime2=gmean1+2.0*gmean2*amean;
                      gmean0=sigma_g0[ktype-1][jtype-1][kptype-1];
                      gmean1=sigma_g1[ktype-1][jtype-1][kptype-1];
                      gmean2=sigma_g2[ktype-1][jtype-1][kptype-1];
                      amean=cosAng[ang_kjkp];
                      gfactor3=gmean0+gmean1*amean
                        +gmean2*amean*amean;
                      gprime3=gmean1+2.0*gmean2*amean;
                      gfactor=gfactor1*gfactor2*gfactor3;
                      rfactorrt=betaS[temp_jk]*betaS[temp_jkp];
                      rfactor=rfactorrt*rfactorrt;

//2nd DD is Eq. 11 (c) for j atom where i , k & k'=neighbor of j

                      DD=DD+2.0*gfactor*rfactor;
                    }
                  }
                }

//j is a neighbor of i, k is a neighbor of j not equal to i and k'
//is a neighbor of k not equal to j or i

                for(ltmp=0;ltmp<numneigh[k];ltmp++) {
                  temp_kkp=BOP_index[k]+ltmp;
                  if(neigh_flag[temp_kkp]) {
                    kp=klist[ltmp];
                    kptype=map[type[kp]]+1;
                    same_ikp=0;
                    same_jkp=0;
                    if(x[i][0]==x[kp][0]) {
                      if(x[i][1]==x[kp][1]) {
                        if(x[i][2]==x[kp][2]) {
                          same_ikp=1;
                        }
                      }
                    }
                    if(x[j][0]==x[kp][0]) {
                      if(x[j][1]==x[kp][1]) {
                        if(x[j][2]==x[kp][2]) {
                          same_jkp=1;
                        }
                      }
                    }
                    if(!same_ikp&&!same_jkp) {
                      for(kNeij=0;kNeij<numneigh[k];kNeij++) {
                        if(x[klist[kNeij]][0]==x[j][0]) {
                          if(x[klist[kNeij]][1]==x[j][1]) {
                            if(x[klist[kNeij]][2]==x[j][2]) {
                              break;
                            }
                          }
                        }
                      }
                      if(kNeij<ltmp) {
                        njkkp=kNeij*(2*numneigh[k]-kNeij-1)/2+(ltmp-kNeij)-1;
                      } else {
                        njkkp=ltmp*(2*numneigh[k]-ltmp-1)/2+(kNeij-ltmp)-1;
                      }
                      sig_flag=0;
                      for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                        ncmp=itypeSigBk[n][nsearch];
                        if(x[ncmp][0]==x[kp][0]) {
                          if(x[ncmp][1]==x[kp][1]) {
                            if(x[ncmp][2]==x[kp][2]) {
                              new2=nsearch;
                              sig_flag=1;
                              break;
                            }
                          }
                        }
                      }
                      if(sig_flag==0) {
                        nSigBk[n]=nSigBk[n]+1;
                        new2=nSigBk[n]-1;
                        itypeSigBk[n][new2]=kp;
                      }
                      ang_jkkp=cos_index[k]+njkkp;
                      if(ang_jkkp>=cos_total) {
                        error->one(FLERR,"Too many atom triplets for pair bop");
                      }
                      gmean0=sigma_g0[jtype-1][ktype-1][kptype-1];
                      gmean1=sigma_g1[jtype-1][ktype-1][kptype-1];
                      gmean2=sigma_g2[jtype-1][ktype-1][kptype-1];
                      amean=cosAng[ang_jkkp];
                      gfactor2=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime2=gmean1+2.0*gmean2*amean;
                      gfactorsq2=gfactor2*gfactor2;
                      gfactor=gfactorsq*gfactorsq2;
                      rfactorrt=betaS[temp_jk]*betaS[temp_kkp];
                      rfactor=rfactorrt*rfactorrt;

//3rd DD is Eq. 11 (c) for j atom where i & k=neighbor of j & k'=neighbor of k

                      DD=DD+gfactor*rfactor;
                    }
                  }
                }
              }
            }
          }

          sig_flag=0;
          if(sig_flag==0) {

// AA and BB are the representations of (a) Eq. 34 and (b) Eq. 9
// for atoms i and j respectively

            AAC=AA+BB;
            BBC=AA*BB;
            CCC=AA*AA+BB*BB;
            DDC=CC+DD;

//EEC is a modified form of (a) Eq. 33

            EEC=(DDC-CCC)/(AAC+2.0*small1);
            for(m=0;m<nb_t;m++) {
              if((bt_sg[m].i>-1)&&(bt_sg[m].j>-1)) {
                bt_i=bt_sg[m].i;
                bt_j=bt_sg[m].j;
                bt_sg[m].dAAC[0]=bt_sg[m].dAA[0]
                    +bt_sg[m].dBB[0];
                bt_sg[m].dAAC[1]=bt_sg[m].dAA[1]
                    +bt_sg[m].dBB[1];
                bt_sg[m].dAAC[2]=bt_sg[m].dAA[2]
                    +bt_sg[m].dBB[2];
              }
            }
            UT=EEC*FF+BBC+small3[iij];
            UT=1.0/sqrt(UT);

// FFC is slightly modified form of (a) Eq. 31
// GGC is slightly modified form of (a) Eq. 32
// bndtmp is a slightly modified form of (a) Eq. 30 and (b) Eq. 8

            bndtmp=(FF+sigma_delta[iij]*sigma_delta[iij])
                +sigma_c[iij]*AAC+small4;
            psign=1.0;
            bndtmp0=1.0/sqrt(bndtmp);
            sigB1[n]=psign*betaS[temp_ij]*bndtmp0;
            bndtmp=-0.5*bndtmp0*bndtmp0*bndtmp0;
            bndtmp1=psign*bndtmp0+psign*betaS[temp_ij]
                *bndtmp*2.0*betaS[temp_ij];
            bndtmp1=bndtmp1*dBetaS[temp_ij]/rij[temp_ij];
            bndtmp2=psign*betaS[temp_ij]*bndtmp*sigma_c[iij];
            setting=0;
            for(m=0;m<nb_t;m++) {
              if((bt_sg[m].i>-1)&&(bt_sg[m].j>-1)) {
                temp_kk=bt_sg[m].temp;
                if(temp_kk==temp_ij&&setting==0) {
                  bt_sg[m].dSigB1[0]=bndtmp1*disij[0][temp_ij]
                      +(bndtmp2*bt_sg[m].dAAC[0]);
                  bt_sg[m].dSigB1[1]=bndtmp1*disij[1][temp_ij]
                      +(bndtmp2*bt_sg[m].dAAC[1]);
                  bt_sg[m].dSigB1[2]=bndtmp1*disij[2][temp_ij]
                      +(bndtmp2*bt_sg[m].dAAC[2]);
                  setting=1;
                }
                else if(temp_kk==temp_ji&&setting==0) {
                  bt_sg[m].dSigB1[0]=-bndtmp1*disij[0][temp_ij]
                      +(bndtmp2*bt_sg[m].dAAC[0]);
                  bt_sg[m].dSigB1[1]=-bndtmp1*disij[1][temp_ij]
                      +(bndtmp2*bt_sg[m].dAAC[1]);
                  bt_sg[m].dSigB1[2]=-bndtmp1*disij[2][temp_ij]
                      +(bndtmp2*bt_sg[m].dAAC[2]);
                  setting=1;
                }
                else {
                  bt_sg[m].dSigB1[0]=(bndtmp2*bt_sg[m].dAAC[0]);
                  bt_sg[m].dSigB1[1]=(bndtmp2*bt_sg[m].dAAC[1]);
                  bt_sg[m].dSigB1[2]=(bndtmp2*bt_sg[m].dAAC[2]);
                }
              }
            }

//This loop is to ensure there is not an error for atoms with no neighbors (deposition)

            if(nb_t==0) {
              if(j>i) {
                bt_sg[0].dSigB1[0]=bndtmp1*disij[0][temp_ij];
                bt_sg[0].dSigB1[1]=bndtmp1*disij[1][temp_ij];
                bt_sg[0].dSigB1[2]=bndtmp1*disij[2][temp_ij];
              }
              else {
                bt_sg[0].dSigB1[0]=-bndtmp1*disij[0][temp_ij];
                bt_sg[0].dSigB1[1]=-bndtmp1*disij[1][temp_ij];
                bt_sg[0].dSigB1[2]=-bndtmp1*disij[2][temp_ij];
              }
              for(pp=0;pp<3;pp++) {
                bt_sg[0].dAA[pp]=0.0;
                bt_sg[0].dBB[pp]=0.0;
                bt_sg[0].dEE1[pp]=0.0;
                bt_sg[0].dFF[pp]=0.0;
                bt_sg[0].dAAC[pp]=0.0;
                bt_sg[0].dSigB[pp]=0.0;
              }
              bt_sg[0].i=i;
              bt_sg[0].j=j;
              bt_sg[0].temp=temp_ij;
              nb_t++;
              if(nb_t>nb_sg) {
                new_n_tot=nb_sg+maxneigh;
                grow_sigma(nb_sg,new_n_tot);
                nb_sg=new_n_tot;
              }
            }
            ps=sigB1[n]*rdBO+1.0;
            ks=(int)ps;
            if(nBOt-1<ks)
              ks=nBOt-1;
            ps=ps-ks;
            if(ps>1.0)
              ps=1.0;
            dsigB1=((FsigBO3[iij][ks-1]*ps+FsigBO2[iij][ks-1])*ps
                +FsigBO1[iij][ks-1])*ps+FsigBO[iij][ks-1];
            dsigB2=(FsigBO6[iij][ks-1]*ps+FsigBO5[iij][ks-1])*ps+FsigBO4[iij][ks-1];
            part0=(FF+0.5*AAC+small5);
            part1=(sigma_f[iij]-0.5)*sigma_k[iij];
            part2=1.0-part1*EE1/part0;
            part3=dsigB1*part1/part0;
            part4=part3/part0*EE1;

// sigB is the final expression for (a) Eq. 6 and (b) Eq. 11

            sigB[n]=dsigB1*part2;
            pp1=2.0*betaS[temp_ij];
            for(m=0;m<nb_t;m++) {
              if((bt_sg[m].i>-1)&&(bt_sg[m].j>-1)) {
                temp_kk=bt_sg[m].temp;
                bt_ij=bt_sg[m].temp;
                bt_i=bt_sg[m].i;
                bt_j=bt_sg[m].j;
                for(pp=0;pp<3;pp++) {
                  bt_sg[m].dSigB[pp]=dsigB2*part2*bt_sg[m].dSigB1[pp]
                      -part3*bt_sg[m].dEE1[pp]
                      +part4*(bt_sg[m].dFF[pp]
                      +0.5*bt_sg[m].dAAC[pp]);
                }
                for(pp=0;pp<3;pp++) {
                  ftmp[pp]=pp1*bt_sg[m].dSigB[pp];
                  f[bt_i][pp]-=ftmp[pp];
                  f[bt_j][pp]+=ftmp[pp];
                }
                if(evflag) {
                  ev_tally_xyz(bt_i,bt_j,nlocal,newton_pair,0.0,0.0,ftmp[0],ftmp[1]
                      ,ftmp[2],disij[0][bt_ij],disij[1][bt_ij],disij[2][bt_ij]);
                }
              }
            }
          }
          n++;
        }
      }
    }
  }
  destroy_sigma();
}

/* ---------------------------------------------------------------------- */

/*  The formulation differs slightly to avoid negative square roots
    in the calculation of Theta_pi,ij of (a) Eq. 36 and (b) Eq. 18 */

void PairBOP::sigmaBo_otf()
{
  int nb_t,new_n_tot;
  int n,i,j,k,kp,m,pp,kpj,kpk,kkp;
  int itmp,jtmp,ktmp,ltmp,mtmp;
  tagint i_tag,j_tag;
  int kp1,kp2,kp1type;
  int iij,iik,ijk,ikkp,ji,iikp,ijkp;
  int nkp;
  int nk0;
  int jNeik,kNeii,kNeij,kNeikp;
  int kpNeij,kpNeik;
  int new1,new2,nlocal;
  int inum,*ilist,*iilist,*jlist,*klist,*kplist;
  int **firstneigh,*numneigh;
  int temp_ij,temp_ik,temp_jkp,temp_kk,temp_jk;
  int temp_ji,temp_kkp;
  int temp_ikp;
  int nb_ij,nb_ik,nb_ikp;
  int nb_jk,nb_jkp,nb_kkp;
  int nsearch;
  int sig_flag,setting,ncmp,ks;
  int itype,jtype,ktype,kptype;
  int bt_i,bt_j;
  int same_ikp,same_jkp,same_kpk;
  int same_jkpj,same_kkpk;
  double AA,BB,CC,DD,EE,EE1,FF;
  double AAC,BBC,CCC,DDC,EEC,FFC,GGC;
  double AACFF,UT,bndtmp,UTcom;
  double amean,gmean0,gmean1,gmean2,ps;
  double gfactor1,gprime1,gsqprime;
  double gfactorsq,gfactor2,gprime2;
  double gfactorsq2,gsqprime2;
  double gfactor3,gprime3,gfactor,rfactor;
  double drfactor,gfactor4,gprime4,agpdpr3;
  double rfactor0,rfactorrt,rfactor1rt,rfactor1;
  double rcm1,rcm2,gcm1,gcm2,gcm3;
  double agpdpr1,agpdpr2,app1,app2,app3,app4;
  double dsigB1,dsigB2;
  double part0,part1,part2,part3,part4;
  double psign,bndtmp0,pp1;
  double bndtmp1,bndtmp2,bndtmp3,bndtmp4,bndtmp5;
  double dis_ij[3],rsq_ij,r_ij;
  double betaS_ij,dBetaS_ij;
  double dis_ik[3],rsq_ik,r_ik;
  double betaS_ik,dBetaS_ik;
  double dis_ikp[3],rsq_ikp,r_ikp;
  double betaS_ikp,dBetaS_ikp;
  double dis_jk[3],rsq_jk,r_jk;
  double betaS_jk,dBetaS_jk;
  double dis_jkp[3],rsq_jkp,r_jkp;
  double betaS_jkp,dBetaS_jkp;
  double dis_kkp[3],rsq_kkp,r_kkp;
  double betaS_kkp,dBetaS_kkp;
  double cosAng_jik,dcA_jik[3][2];
  double cosAng_jikp,dcA_jikp[3][2];
  double cosAng_kikp,dcA_kikp[3][2];
  double cosAng_ijk,dcA_ijk[3][2];
  double cosAng_ijkp,dcA_ijkp[3][2];
  double cosAng_kjkp,dcA_kjkp[3][2];
  double cosAng_ikj,dcA_ikj[3][2];
  double cosAng_ikkp,dcA_ikkp[3][2];
  double cosAng_jkkp,dcA_jkkp[3][2];
  double cosAng_jkpk,dcA_jkpk[3][2];

  double ftmp[3],xtmp[3];
  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int newton_pair = force->newton_pair;
  int *type = atom->type;

  nlocal = atom->nlocal;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  n=0;
  if(nb_sg==0) {
    nb_sg=(maxneigh)*(maxneigh/2);
  }
  if(allocate_sigma) {
    destroy_sigma();
  }

  create_sigma(nb_sg);

  for(itmp=0;itmp<inum;itmp++) {
    i = ilist[itmp];
    i_tag=tag[i];
    itype = map[type[i]]+1;

//j is loop over all neighbors of i

    for(jtmp=0;jtmp<numneigh[i];jtmp++) {
      for(m=0;m<nb_sg;m++) {
        for(pp=0;pp<3;pp++) {
          bt_sg[m].dAA[pp]=0.0;
          bt_sg[m].dBB[pp]=0.0;
          bt_sg[m].dCC[pp]=0.0;
          bt_sg[m].dDD[pp]=0.0;
          bt_sg[m].dEE[pp]=0.0;
          bt_sg[m].dEE1[pp]=0.0;
          bt_sg[m].dFF[pp]=0.0;
          bt_sg[m].dAAC[pp]=0.0;
          bt_sg[m].dBBC[pp]=0.0;
          bt_sg[m].dCCC[pp]=0.0;
          bt_sg[m].dDDC[pp]=0.0;
          bt_sg[m].dEEC[pp]=0.0;
          bt_sg[m].dFFC[pp]=0.0;
          bt_sg[m].dGGC[pp]=0.0;
          bt_sg[m].dUT[pp]=0.0;
          bt_sg[m].dSigB1[pp]=0.0;
          bt_sg[m].dSigB[pp]=0.0;
        }
        bt_sg[m].i=-1;
        bt_sg[m].j=-1;
        bt_sg[m].temp=-1;
      }
      nb_t=0;
      iilist=firstneigh[i];
      temp_ij=BOP_index[i]+jtmp;
      j=iilist[jtmp];
      jlist=firstneigh[j];
      j_tag=tag[j];
      jtype = map[type[j]]+1;
      nb_ij=nb_t;
      nb_t++;
      if(nb_t>nb_sg) {
        new_n_tot=nb_sg+maxneigh;
        grow_sigma(nb_sg,new_n_tot);
        nb_sg=new_n_tot;
      }
      bt_sg[nb_ij].temp=temp_ij;
      bt_sg[nb_ij].i=i;
      bt_sg[nb_ij].j=j;
      if(j_tag>=i_tag) {
        if(itype==jtype)
          iij=itype-1;
        else if(itype<jtype)
          iij=itype*bop_types-itype*(itype+1)/2+jtype-1;
        else
          iij=jtype*bop_types-jtype*(jtype+1)/2+itype-1;
        for(ji=0;ji<numneigh[j];ji++) {
          temp_ji=BOP_index[j]+ji;
          if(x[jlist[ji]][0]==x[i][0]) {
            if(x[jlist[ji]][1]==x[i][1]) {
              if(x[jlist[ji]][2]==x[i][2]) {
                break;
              }
            }
          }
        }
        dis_ij[0]=x[j][0]-x[i][0];
        dis_ij[1]=x[j][1]-x[i][1];
        dis_ij[2]=x[j][2]-x[i][2];
        rsq_ij=dis_ij[0]*dis_ij[0]
            +dis_ij[1]*dis_ij[1]
            +dis_ij[2]*dis_ij[2];
        r_ij=sqrt(rsq_ij);

        if(r_ij<rcut[iij]) {
          ps=r_ij*rdr[iij]+1.0;
          ks=(int)ps;
          if(nr-1<ks)
            ks=nr-1;
          ps=ps-ks;
          if(ps>1.0)
            ps=1.0;
          betaS_ij=((pBetaS3[iij][ks-1]*ps+pBetaS2[iij][ks-1])*ps
              +pBetaS1[iij][ks-1])*ps+pBetaS[iij][ks-1];
          dBetaS_ij=(pBetaS6[iij][ks-1]*ps+pBetaS5[iij][ks-1])*ps
              +pBetaS4[iij][ks-1];
          nSigBk[n]=0;

//AA-EE1 are the components making up Eq. 30 (a)

          AA=0.0;
          BB=0.0;
          CC=0.0;
          DD=0.0;
          EE=0.0;
          EE1=0.0;

//FF is the Beta_sigma^2 term

          FF=betaS_ij*betaS_ij;

//agpdpr1 is derivative of FF w.r.t. r_ij

          agpdpr1=2.0*betaS_ij*dBetaS_ij/r_ij;

//dXX derivatives are taken with respect to all pairs contributing to the energy
//nb_ij is derivative w.r.t. ij pair

          bt_sg[nb_ij].dFF[0]=agpdpr1*dis_ij[0];
          bt_sg[nb_ij].dFF[1]=agpdpr1*dis_ij[1];
          bt_sg[nb_ij].dFF[2]=agpdpr1*dis_ij[2];

//k is loop over all neighbors of i again with j neighbor of i

          for(ktmp=0;ktmp<numneigh[i];ktmp++) {
            temp_ik=BOP_index[i]+ktmp;
            if(ktmp!=jtmp) {
              k=iilist[ktmp];
              klist=firstneigh[k];
              ktype = map[type[k]]+1;
              if(itype==ktype)
                iik=itype-1;
              else if(itype<ktype)
                iik=itype*bop_types-itype*(itype+1)/2+ktype-1;
              else
                iik=ktype*bop_types-ktype*(ktype+1)/2+itype-1;

//find neighbor of k that is equal to i

              for(kNeii=0;kNeii<numneigh[k];kNeii++) {
                if(x[klist[kNeii]][0]==x[i][0]) {
                  if(x[klist[kNeii]][1]==x[i][1]) {
                    if(x[klist[kNeii]][2]==x[i][2]) {
                      break;
                    }
                  }
                }
              }
              dis_ik[0]=x[k][0]-x[i][0];
              dis_ik[1]=x[k][1]-x[i][1];
              dis_ik[2]=x[k][2]-x[i][2];
              rsq_ik=dis_ik[0]*dis_ik[0]
                  +dis_ik[1]*dis_ik[1]
                  +dis_ik[2]*dis_ik[2];
              r_ik=sqrt(rsq_ik);
              if(r_ik<=rcut[iik]) {
                ps=r_ik*rdr[iik]+1.0;
                ks=(int)ps;
                if(nr-1<ks)
                  ks=nr-1;
                ps=ps-ks;
                if(ps>1.0)
                  ps=1.0;
                betaS_ik=((pBetaS3[iik][ks-1]*ps+pBetaS2[iik][ks-1])*ps
                    +pBetaS1[iik][ks-1])*ps+pBetaS[iik][ks-1];
                dBetaS_ik=(pBetaS6[iik][ks-1]*ps+pBetaS5[iik][ks-1])*ps
                    +pBetaS4[iik][ks-1];

//find neighbor of i that is equal to k

                for(jNeik=0;jNeik<numneigh[j];jNeik++) {
                  temp_jk=BOP_index[j]+jNeik;
                  if(x[jlist[jNeik]][0]==x[k][0]) {
                    if(x[jlist[jNeik]][1]==x[k][1]) {
                      if(x[jlist[jNeik]][2]==x[k][2]) {
                        break;
                      }
                    }
                  }
                }

//find neighbor of k that is equal to j

                for(kNeij=0;kNeij<numneigh[k];kNeij++) {
                  if(x[klist[kNeij]][0]==x[j][0]) {
                    if(x[klist[kNeij]][1]==x[j][1]) {
                      if(x[klist[kNeij]][2]==x[j][2]) {
                        break;
                      }
                    }
                  }
                }
                dis_jk[0]=x[k][0]-x[j][0];
                dis_jk[1]=x[k][1]-x[j][1];
                dis_jk[2]=x[k][2]-x[j][2];
                rsq_jk=dis_jk[0]*dis_jk[0]
                    +dis_jk[1]*dis_jk[1]
                    +dis_jk[2]*dis_jk[2];
                r_jk=sqrt(rsq_jk);

                sig_flag=0;
                for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                  ncmp=itypeSigBk[n][nsearch];
                  if(x[ncmp][0]==x[k][0]) {
                    if(x[ncmp][1]==x[k][1]) {
                      if(x[ncmp][2]==x[k][2]) {
                        nk0=nsearch;
                        sig_flag=1;
                        break;
                      }
                    }
                  }
                }
                if(sig_flag==0) {
                  nSigBk[n]=nSigBk[n]+1;
                  nk0=nSigBk[n]-1;
                  itypeSigBk[n][nk0]=k;
                }
                nb_ik=nb_t;
                nb_t++;
                if(nb_t>nb_sg) {
                  new_n_tot=nb_sg+maxneigh;
                  grow_sigma(nb_sg,new_n_tot);
                  nb_sg=new_n_tot;
                }
                bt_sg[nb_ik].temp=temp_ik;
                bt_sg[nb_ik].i=i;
                bt_sg[nb_ik].j=k;
                nb_jk=nb_t;
                nb_t++;
                if(nb_t>nb_sg) {
                  new_n_tot=nb_sg+maxneigh;
                  grow_sigma(nb_sg,new_n_tot);
                  nb_sg=new_n_tot;
                }
                bt_sg[nb_jk].temp=temp_jk;
                bt_sg[nb_jk].i=j;
                bt_sg[nb_jk].j=k;
                cosAng_jik=(dis_ij[0]*dis_ik[0]+dis_ij[1]*dis_ik[1]
                    +dis_ij[2]*dis_ik[2])/(r_ij*r_ik);
                dcA_jik[0][0]=(dis_ik[0]*r_ij*r_ik-cosAng_jik
                    *dis_ij[0]*r_ik*r_ik)/(r_ij*r_ij*r_ik*r_ik);
                dcA_jik[1][0]=(dis_ik[1]*r_ij*r_ik-cosAng_jik
                    *dis_ij[1]*r_ik*r_ik)/(r_ij*r_ij*r_ik*r_ik);
                dcA_jik[2][0]=(dis_ik[2]*r_ij*r_ik-cosAng_jik
                    *dis_ij[2]*r_ik*r_ik)/(r_ij*r_ij*r_ik*r_ik);
                dcA_jik[0][1]=(dis_ij[0]*r_ij*r_ik-cosAng_jik
                    *dis_ik[0]*r_ij*r_ij)/(r_ij*r_ij*r_ik*r_ik);
                dcA_jik[1][1]=(dis_ij[1]*r_ij*r_ik-cosAng_jik
                    *dis_ik[1]*r_ij*r_ij)/(r_ij*r_ij*r_ik*r_ik);
                dcA_jik[2][1]=(dis_ij[2]*r_ij*r_ik-cosAng_jik
                    *dis_ik[2]*r_ij*r_ij)/(r_ij*r_ij*r_ik*r_ik);
                gmean0=sigma_g0[jtype-1][itype-1][ktype-1];
                gmean1=sigma_g1[jtype-1][itype-1][ktype-1];
                gmean2=sigma_g2[jtype-1][itype-1][ktype-1];
                amean=cosAng_jik;
                gfactor1=gmean0+gmean1*amean
                    +gmean2*amean*amean;
                gfactorsq=gfactor1*gfactor1;
                gprime1=gmean1+2.0*gmean2*amean;
                gsqprime=2.0*gfactor1*gprime1;

//AA is Eq. 34 (a) or Eq. 10 (c) for the i atom
//1st CC is Eq. 11 (c) for i atom where j & k=neighbor of i

                AA=AA+gfactorsq*betaS_ik*betaS_ik;
                CC=CC+gfactorsq*betaS_ik*betaS_ik*betaS_ik*betaS_ik;

//agpdpr1 is derivative of AA w.r.t. Beta(rik)
//app1 is derivative of AA w.r.t. cos(theta_jik)

                agpdpr1=2.0*gfactorsq*betaS_ik*dBetaS_ik/r_ik;
                app1=betaS_ik*betaS_ik*gsqprime;
                bt_sg[nb_ij].dAA[0]+=
                    app1*dcA_jik[0][0];
                bt_sg[nb_ij].dAA[1]+=
                    app1*dcA_jik[1][0];
                bt_sg[nb_ij].dAA[2]+=
                    app1*dcA_jik[2][0];
                bt_sg[nb_ij].dCC[0]+=
                    app2*dcA_jik[0][0];
                bt_sg[nb_ij].dCC[1]+=
                    app2*dcA_jik[1][0];
                bt_sg[nb_ij].dCC[2]+=
                    app2*dcA_jik[2][0];
                bt_sg[nb_ik].dAA[0]+=
                    app1*dcA_jik[0][1]
                    +agpdpr1*dis_ik[0];
                bt_sg[nb_ik].dAA[1]+=
                    app1*dcA_jik[1][1]
                    +agpdpr1*dis_ik[1];
                bt_sg[nb_ik].dAA[2]+=
                    app1*dcA_jik[2][1]
                    +agpdpr1*dis_ik[2];
                bt_sg[nb_ik].dCC[0]+=
                    app2*dcA_jik[0][1]
                    +agpdpr2*dis_ik[0];
                bt_sg[nb_ik].dCC[1]+=
                    app2*dcA_jik[1][1]
                    +agpdpr2*dis_ik[1];
                bt_sg[nb_ik].dCC[2]+=
                    app2*dcA_jik[2][1]
                    +agpdpr2*dis_ik[2];

//k' is loop over neighbors all neighbors of j with k a neighbor
//of i and j a neighbor of i and determine which k' is k

                same_kpk=0;
                for(ltmp=0;ltmp<numneigh[j];ltmp++) {
                  temp_jkp=BOP_index[j]+ltmp;
                  kp1=jlist[ltmp];
                  kp1type=map[type[kp1]]+1;
                  if(x[kp1][0]==x[k][0]) {
                    if(x[kp1][1]==x[k][1]) {
                      if(x[kp1][2]==x[k][2]) {
                        same_kpk=1;
                        break;
                      }
                    }
                  }
                }
                if(same_kpk){

//loop over neighbors of k

                  for(mtmp=0;mtmp<numneigh[k];mtmp++) {
                    kp2=klist[mtmp];
                    if(x[kp2][0]==x[k][0]) {
                      if(x[kp2][1]==x[k][1]) {
                        if(x[kp2][2]==x[k][2]) {
                          break;
                        }
                      }
                    }
                  }
                  if(jtype==ktype)
                    ijk=jtype-1;
                  else if(jtype < ktype)
                    ijk=jtype*bop_types-jtype*(jtype+1)/2+ktype-1;
                  else
                    ijk=ktype*bop_types-ktype*(ktype+1)/2+jtype-1;
                  if(jtype==kp1type)
                    ijkp=jtype-1;
                  else if(jtype<kp1type)
                    ijkp=jtype*bop_types-jtype*(jtype+1)/2+kp1type-1;
                  else
                    ijkp=kp1type*bop_types-kp1type*(kp1type+1)/2+jtype-1;

                  dis_jkp[0]=x[kp1][0]-x[j][0];
                  dis_jkp[1]=x[kp1][1]-x[j][1];
                  dis_jkp[2]=x[kp1][2]-x[j][2];
                  rsq_jkp=dis_jkp[0]*dis_jkp[0]
                      +dis_jkp[1]*dis_jkp[1]
                      +dis_jkp[2]*dis_jkp[2];
                  r_jkp=sqrt(rsq_jkp);
                  if(r_jkp<=rcut[ijkp]) {
                    ps=r_jkp*rdr[ijkp]+1.0;
                    ks=(int)ps;
                    if(nr-1<ks)
                      ks=nr-1;
                    ps=ps-ks;
                    if(ps>1.0)
                      ps=1.0;
                    betaS_jkp=((pBetaS3[ijkp][ks-1]*ps+pBetaS2[ijkp][ks-1])*ps
                        +pBetaS1[ijkp][ks-1])*ps+pBetaS[ijkp][ks-1];
                    dBetaS_jkp=(pBetaS6[ijkp][ks-1]*ps+pBetaS5[ijkp][ks-1])*ps
                        +pBetaS4[ijkp][ks-1];
                    cosAng_ijk=(-dis_ij[0]*dis_jk[0]-dis_ij[1]*dis_jk[1]
                        -dis_ij[2]*dis_jk[2])/(r_ij*r_jk);
                    dcA_ijk[0][0]=(dis_jk[0]*r_ij*r_jk-cosAng_ijk
                        *-dis_ij[0]*r_jk*r_jk)/(r_ij*r_ij*r_jk*r_jk);
                    dcA_ijk[1][0]=(dis_jk[1]*r_ij*r_jk-cosAng_ijk
                        *-dis_ij[1]*r_jk*r_jk)/(r_ij*r_ij*r_jk*r_jk);
                    dcA_ijk[2][0]=(dis_jk[2]*r_ij*r_jk-cosAng_ijk
                        *-dis_ij[2]*r_jk*r_jk)/(r_ij*r_ij*r_jk*r_jk);
                    dcA_ijk[0][1]=(-dis_ij[0]*r_ij*r_jk-cosAng_ijk
                        *dis_jk[0]*r_ij*r_ij)/(r_ij*r_ij*r_jk*r_jk);
                    dcA_ijk[1][1]=(-dis_ij[1]*r_ij*r_jk-cosAng_ijk
                        *dis_jk[1]*r_ij*r_ij)/(r_ij*r_ij*r_jk*r_jk);
                    dcA_ijk[2][1]=(-dis_ij[2]*r_ij*r_jk-cosAng_ijk
                        *dis_jk[2]*r_ij*r_ij)/(r_ij*r_ij*r_jk*r_jk);
                    gmean0=sigma_g0[itype-1][jtype-1][ktype-1];
                    gmean1=sigma_g1[itype-1][jtype-1][ktype-1];
                    gmean2=sigma_g2[itype-1][jtype-1][ktype-1];
                    amean=cosAng_ijk;
                    gfactor2=gmean0+gmean1*amean
                        +gmean2*amean*amean;
                    gprime2=gmean1+2.0*gmean2*amean;
                    gmean0=sigma_g0[itype-1][ktype-1][jtype-1];
                    gmean1=sigma_g1[itype-1][ktype-1][jtype-1];
                    gmean2=sigma_g2[itype-1][ktype-1][jtype-1];
                    cosAng_ikj=(dis_ik[0]*dis_jk[0]+dis_ik[1]*dis_jk[1]
                        +dis_ik[2]*dis_jk[2])/(r_ik*r_jk);
                    dcA_ikj[0][0]=(-dis_jk[0]*r_ik*r_jk-cosAng_ikj
                        *-dis_ik[0]*r_jk*r_jk)/(r_ik*r_ik*r_jk*r_jk);
                    dcA_ikj[1][0]=(-dis_jk[1]*r_ik*r_jk-cosAng_ikj
                        *-dis_ik[1]*r_jk*r_jk)/(r_ik*r_ik*r_jk*r_jk);
                    dcA_ikj[2][0]=(-dis_jk[2]*r_ik*r_jk-cosAng_ikj
                        *-dis_ik[2]*r_jk*r_jk)/(r_ik*r_ik*r_jk*r_jk);
                    dcA_ikj[0][1]=(-dis_ik[0]*r_ik*r_jk-cosAng_ikj
                        *-dis_jk[0]*r_ik*r_ik)/(r_ik*r_ik*r_jk*r_jk);
                    dcA_ikj[1][1]=(-dis_ik[1]*r_ik*r_jk-cosAng_ikj
                        *-dis_jk[1]*r_ik*r_ik)/(r_ik*r_ik*r_jk*r_jk);
                    dcA_ikj[2][1]=(-dis_ik[2]*r_ik*r_jk-cosAng_ikj
                        *-dis_jk[2]*r_ik*r_ik)/(r_ik*r_ik*r_jk*r_jk);
                    amean=cosAng_ikj;
                    gfactor3=gmean0+gmean1*amean
                        +gmean2*amean*amean;
                    gprime3=gmean1+2.0*gmean2*amean;
                    gfactor=gfactor1*gfactor2*gfactor3;
                    rfactor=betaS_ik*betaS_jkp;

//EE1 is (b) Eq. 12

                    EE1=EE1+gfactor*rfactor;

//rcm1 is derivative of EE1 w.r.t Beta(r_ik)
//rcm2 is derivative of EE1 w.r.t Beta(r_jk')
//gcm1 is derivative of EE1 w.r.t cos(theta_jik)
//gcm2 is derivative of EE1 w.r.t cos(theta_ijk)
//gcm3 is derivative of EE1 w.r.t cos(theta_ikj)

                    rcm1=gfactor*betaS_jkp*dBetaS_ik/r_ik;
                    rcm2=gfactor*betaS_ik*dBetaS_jkp/r_jkp;
                    gcm1=rfactor*gprime1*gfactor2*gfactor3;
                    gcm2=rfactor*gfactor1*gprime2*gfactor3;
                    gcm3=rfactor*gfactor1*gfactor2*gprime3;
                    bt_sg[nb_ij].dEE1[0]+=
                        gcm1*dcA_jik[0][0]
                        -gcm2*dcA_ijk[0][0];
                    bt_sg[nb_ij].dEE1[1]+=
                        gcm1*dcA_jik[1][0]
                        -gcm2*dcA_ijk[1][0];
                    bt_sg[nb_ij].dEE1[2]+=
                        gcm1*dcA_jik[2][0]
                        -gcm2*dcA_ijk[2][0];
                    bt_sg[nb_ik].dEE1[0]+=
                        gcm1*dcA_jik[0][1]
                        +rcm1*dis_ik[0]
                        -gcm3*dcA_ikj[0][0];
                    bt_sg[nb_ik].dEE1[1]+=
                        gcm1*dcA_jik[1][1]
                        +rcm1*dis_ik[1]
                        -gcm3*dcA_ikj[1][0];
                    bt_sg[nb_ik].dEE1[2]+=
                        gcm1*dcA_jik[2][1]
                        +rcm1*dis_ik[2]
                        -gcm3*dcA_ikj[2][0];
                    bt_sg[nb_jk].dEE1[0]+=
                        gcm2*dcA_ijk[0][1]
                        +rcm2*dis_jkp[0]
                        -gcm3*dcA_ikj[0][1];
                    bt_sg[nb_jk].dEE1[1]+=
                        gcm2*dcA_ijk[1][1]
                        +rcm2*dis_jkp[1]
                        -gcm3*dcA_ikj[1][1];
                    bt_sg[nb_jk].dEE1[2]+=
                        gcm2*dcA_ijk[2][1]
                        +rcm2*dis_jkp[2]
                        -gcm3*dcA_ikj[2][1];
                  }
                }

// k and k' and j are all different neighbors of i

                for(ltmp=0;ltmp<ktmp;ltmp++) {
                  if(ltmp!=jtmp) {
                    temp_ikp=BOP_index[i]+ltmp;
                    kp=iilist[ltmp];;
                    kptype = map[type[kp]]+1;
                    if(itype==kptype)
                      iikp=itype-1;
                    else if(itype<kptype)
                      iikp=itype*bop_types-itype*(itype+1)/2+kptype-1;
                    else
                      iikp=kptype*bop_types-kptype*(kptype+1)/2+itype-1;
                    for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                      ncmp=itypeSigBk[n][nsearch];
                      if(x[ncmp][0]==x[kp][0]) {
                        if(x[ncmp][1]==x[kp][1]) {
                          if(x[ncmp][2]==x[kp][2]) {
                            break;
                          }
                        }
                      }
                    }
                    dis_ikp[0]=x[kp][0]-x[i][0];
                    dis_ikp[1]=x[kp][1]-x[i][1];
                    dis_ikp[2]=x[kp][2]-x[i][2];
                    rsq_ikp=dis_ikp[0]*dis_ikp[0]
                        +dis_ikp[1]*dis_ikp[1]
                        +dis_ikp[2]*dis_ikp[2];
                    r_ikp=sqrt(rsq_ikp);
                    if(r_ikp<=rcut[iikp]) {
                      ps=r_ikp*rdr[iikp]+1.0;
                      ks=(int)ps;
                      if(nr-1<ks)
                        ks=nr-1;
                      ps=ps-ks;
                      if(ps>1.0)
                        ps=1.0;
                      betaS_ikp=((pBetaS3[iikp][ks-1]*ps+pBetaS2[iikp][ks-1])*ps
                          +pBetaS1[iikp][ks-1])*ps+pBetaS[iikp][ks-1];
                      dBetaS_ikp=(pBetaS6[iikp][ks-1]*ps+pBetaS5[iikp][ks-1])*ps
                          +pBetaS4[iikp][ks-1];
                      nb_ikp=nb_t;
                      nb_t++;
                      if(nb_t>nb_sg) {
                        new_n_tot=nb_sg+maxneigh;
                        grow_sigma(nb_sg,new_n_tot);
                        nb_sg=new_n_tot;
                      }
                      bt_sg[nb_ikp].temp=temp_ikp;
                      bt_sg[nb_ikp].i=i;
                      bt_sg[nb_ikp].j=kp;
                      gmean0=sigma_g0[jtype-1][itype-1][kptype-1];
                      gmean1=sigma_g1[jtype-1][itype-1][kptype-1];
                      gmean2=sigma_g2[jtype-1][itype-1][kptype-1];
                      cosAng_jikp=(dis_ij[0]*dis_ikp[0]+dis_ij[1]*dis_ikp[1]
                          +dis_ij[2]*dis_ikp[2])/(r_ij*r_ikp);
                      dcA_jikp[0][0]=(dis_ikp[0]*r_ij*r_ikp-cosAng_jikp
                          *dis_ij[0]*r_ikp*r_ikp)/(r_ij*r_ij*r_ikp*r_ikp);
                      dcA_jikp[1][0]=(dis_ikp[1]*r_ij*r_ikp-cosAng_jikp
                          *dis_ij[1]*r_ikp*r_ikp)/(r_ij*r_ij*r_ikp*r_ikp);
                      dcA_jikp[2][0]=(dis_ikp[2]*r_ij*r_ikp-cosAng_jikp
                          *dis_ij[2]*r_ikp*r_ikp)/(r_ij*r_ij*r_ikp*r_ikp);
                      dcA_jikp[0][1]=(dis_ij[0]*r_ij*r_ikp-cosAng_jikp
                          *dis_ikp[0]*r_ij*r_ij)/(r_ij*r_ij*r_ikp*r_ikp);
                      dcA_jikp[1][1]=(dis_ij[1]*r_ij*r_ikp-cosAng_jikp
                          *dis_ikp[1]*r_ij*r_ij)/(r_ij*r_ij*r_ikp*r_ikp);
                      dcA_jikp[2][1]=(dis_ij[2]*r_ij*r_ikp-cosAng_jikp
                          *dis_ikp[2]*r_ij*r_ij)/(r_ij*r_ij*r_ikp*r_ikp);
                      cosAng_kikp=(dis_ik[0]*dis_ikp[0]+dis_ik[1]*dis_ikp[1]
                          +dis_ik[2]*dis_ikp[2])/(r_ik*r_ikp);
                      dcA_kikp[0][0]=(dis_ikp[0]*r_ik*r_ikp-cosAng_kikp
                          *dis_ik[0]*r_ikp*r_ikp)/(r_ik*r_ik*r_ikp*r_ikp);
                      dcA_kikp[1][0]=(dis_ikp[1]*r_ik*r_ikp-cosAng_kikp
                          *dis_ik[1]*r_ikp*r_ikp)/(r_ik*r_ik*r_ikp*r_ikp);
                      dcA_kikp[2][0]=(dis_ikp[2]*r_ik*r_ikp-cosAng_kikp
                          *dis_ik[2]*r_ikp*r_ikp)/(r_ik*r_ik*r_ikp*r_ikp);
                      dcA_kikp[0][1]=(dis_ik[0]*r_ik*r_ikp-cosAng_kikp
                          *dis_ikp[0]*r_ik*r_ik)/(r_ik*r_ik*r_ikp*r_ikp);
                      dcA_kikp[1][1]=(dis_ik[1]*r_ik*r_ikp-cosAng_kikp
                          *dis_ikp[1]*r_ik*r_ik)/(r_ik*r_ik*r_ikp*r_ikp);
                      dcA_kikp[2][1]=(dis_ik[2]*r_ik*r_ikp-cosAng_kikp
                          *dis_ikp[2]*r_ik*r_ik)/(r_ik*r_ik*r_ikp*r_ikp);
                      amean=cosAng_jikp;
                      gfactor2=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime2=gmean1+2.0*gmean2*amean;
                      gmean0=sigma_g0[ktype-1][itype-1][kptype-1];
                      gmean1=sigma_g1[ktype-1][itype-1][kptype-1];
                      gmean2=sigma_g2[ktype-1][itype-1][kptype-1];
                      amean=cosAng_kikp;
                      gfactor3=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime3=gmean1+2.0*gmean2*amean;
                      gfactor=gfactor1*gfactor2*gfactor3;
                      rfactorrt=betaS_ik*betaS_ikp;
                      rfactor=rfactorrt*rfactorrt;

//2nd CC is second term of Eq. 11 (c) for i atom where j , k & k' =neighbor of i

                      CC=CC+2.0*gfactor*rfactor;

//agpdpr1 is derivative of CC 2nd term w.r.t. Beta(r_ik)
//agpdpr2 is derivative of CC 2nd term w.r.t. Beta(r_ik')
//app1 is derivative of CC 2nd term w.r.t. cos(theta_jik)
//app2 is derivative of CC 2nd term w.r.t. cos(theta_jik')
//app3 is derivative of CC 2nd term w.r.t. cos(theta_kik')

                      agpdpr1=4.0*gfactor*rfactorrt*betaS_ikp
                          *dBetaS_ik/r_ik;
                      agpdpr2=4.0*gfactor*rfactorrt*betaS_ik
                          *dBetaS_ikp/r_ikp;
                      app1=2.0*rfactor*gfactor2*gfactor3*gprime1;
                      app2=2.0*rfactor*gfactor1*gfactor3*gprime2;
                      app3=2.0*rfactor*gfactor1*gfactor2*gprime3;
                      bt_sg[nb_ij].dCC[0]+=
                          app1*dcA_jik[0][0]
                          +app2*dcA_jikp[0][0];
                      bt_sg[nb_ij].dCC[1]+=
                          app1*dcA_jik[1][0]
                          +app2*dcA_jikp[1][0];
                      bt_sg[nb_ij].dCC[2]+=
                          app1*dcA_jik[2][0]
                          +app2*dcA_jikp[2][0];
                      bt_sg[nb_ik].dCC[0]+=
                          app1*dcA_jik[0][1]
                          +app3*dcA_kikp[0][0]
                          +agpdpr1*dis_ik[0];
                      bt_sg[nb_ik].dCC[1]+=
                          app1*dcA_jik[1][1]
                          +app3*dcA_kikp[1][0]
                          +agpdpr1*dis_ik[1];
                      bt_sg[nb_ik].dCC[2]+=
                          app1*dcA_jik[2][1]
                          +app3*dcA_kikp[2][0]
                          +agpdpr1*dis_ik[2];
                      bt_sg[nb_ikp].dCC[0]=
                          app2*dcA_jikp[0][1]
                          +app3*dcA_kikp[0][1]
                          +agpdpr2*dis_ikp[0];
                      bt_sg[nb_ikp].dCC[1]=
                          app2*dcA_jikp[1][1]
                          +app3*dcA_kikp[1][1]
                          +agpdpr2*dis_ikp[1];
                      bt_sg[nb_ikp].dCC[2]=
                          app2*dcA_jikp[2][1]
                          +app3*dcA_kikp[2][1]
                          +agpdpr2*dis_ikp[2];
                    }
                  }
                }

// j and k are different neighbors of i and k' is a neighbor k not equal to i

                for(ltmp=0;ltmp<numneigh[k];ltmp++) {
                  temp_kkp=BOP_index[k]+ltmp;
                  kp=klist[ltmp];;
                  kptype = map[type[kp]]+1;
                  same_ikp=0;
                  same_jkp=0;
                  if(x[i][0]==x[kp][0]) {
                    if(x[i][1]==x[kp][1]) {
                      if(x[i][2]==x[kp][2]) {
                        same_ikp=1;
                      }
                    }
                  }
                  if(x[j][0]==x[kp][0]) {
                    if(x[j][1]==x[kp][1]) {
                      if(x[j][2]==x[kp][2]) {
                        same_jkp=1;
                      }
                    }
                  }
                  if(!same_ikp&&!same_jkp) {
                    if(ktype==kptype)
                      ikkp=ktype-1;
                    else if(ktype<kptype)
                      ikkp=ktype*bop_types-ktype*(ktype+1)/2+kptype-1;
                    else
                      ikkp=kptype*bop_types-kptype*(kptype+1)/2+ktype-1;
                    dis_kkp[0]=x[kp][0]-x[k][0];
                    dis_kkp[1]=x[kp][1]-x[k][1];
                    dis_kkp[2]=x[kp][2]-x[k][2];
                    rsq_kkp=dis_kkp[0]*dis_kkp[0]
                        +dis_kkp[1]*dis_kkp[1]
                        +dis_kkp[2]*dis_kkp[2];
                    r_kkp=sqrt(rsq_kkp);
                    if(r_kkp<=rcut[ikkp]) {
                      ps=r_kkp*rdr[ikkp]+1.0;
                      ks=(int)ps;
                      if(nr-1<ks)
                        ks=nr-1;
                      ps=ps-ks;
                      if(ps>1.0)
                        ps=1.0;
                      betaS_kkp=((pBetaS3[ikkp][ks-1]*ps+pBetaS2[ikkp][ks-1])*ps
                          +pBetaS1[ikkp][ks-1])*ps+pBetaS[ikkp][ks-1];
                      dBetaS_kkp=(pBetaS6[ikkp][ks-1]*ps+pBetaS5[ikkp][ks-1])*ps
                          +pBetaS4[ikkp][ks-1];
                      sig_flag=0;
                      for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                        ncmp=itypeSigBk[n][nsearch];
                        if(x[ncmp][0]==x[kp][0]) {
                          if(x[ncmp][1]==x[kp][1]) {
                            if(x[ncmp][2]==x[kp][2]) {
                              sig_flag=1;
                              nkp=nsearch;
                              break;
                            }
                          }
                        }
                      }
                      if(sig_flag==0) {
                        nSigBk[n]=nSigBk[n]+1;
                        nkp=nSigBk[n]-1;
                        itypeSigBk[n][nkp]=kp;
                      }
                      cosAng_ikkp=(-dis_ik[0]*dis_kkp[0]-dis_ik[1]*dis_kkp[1]
                          -dis_ik[2]*dis_kkp[2])/(r_ik*r_kkp);
                      dcA_ikkp[0][0]=(dis_kkp[0]*r_ik*r_kkp-cosAng_ikkp
                          *-dis_ik[0]*r_kkp*r_kkp)/(r_ik*r_ik*r_kkp*r_kkp);
                      dcA_ikkp[1][0]=(dis_kkp[1]*r_ik*r_kkp-cosAng_ikkp
                          *-dis_ik[1]*r_kkp*r_kkp)/(r_ik*r_ik*r_kkp*r_kkp);
                      dcA_ikkp[2][0]=(dis_kkp[2]*r_ik*r_kkp-cosAng_ikkp
                          *-dis_ik[2]*r_kkp*r_kkp)/(r_ik*r_ik*r_kkp*r_kkp);
                      dcA_ikkp[0][1]=(-dis_ik[0]*r_ik*r_kkp-cosAng_ikkp
                          *dis_kkp[0]*r_ik*r_ik)/(r_ik*r_ik*r_kkp*r_kkp);
                      dcA_ikkp[1][1]=(-dis_ik[1]*r_ik*r_kkp-cosAng_ikkp
                          *dis_kkp[1]*r_ik*r_ik)/(r_ik*r_ik*r_kkp*r_kkp);
                      dcA_ikkp[2][1]=(-dis_ik[2]*r_ik*r_kkp-cosAng_ikkp
                          *dis_kkp[2]*r_ik*r_ik)/(r_ik*r_ik*r_kkp*r_kkp);
                      nb_kkp=nb_t;
                      nb_t++;
                      if(nb_t>nb_sg) {
                        new_n_tot=nb_sg+maxneigh;
                        grow_sigma(nb_sg,new_n_tot);
                        nb_sg=new_n_tot;
                      }
                      bt_sg[nb_kkp].temp=temp_kkp;
                      bt_sg[nb_kkp].i=k;
                      bt_sg[nb_kkp].j=kp;
                      gmean0=sigma_g0[itype-1][ktype-1][kptype-1];
                      gmean1=sigma_g1[itype-1][ktype-1][kptype-1];
                      gmean2=sigma_g2[itype-1][ktype-1][kptype-1];
                      amean=cosAng_ikkp;
                      gfactor2=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime2=gmean1+2.0*gmean2*amean;
                      gfactorsq2=gfactor2*gfactor2;
                      gsqprime2=2.0*gfactor2*gprime2;
                      gfactor=gfactorsq*gfactorsq2;
                      rfactorrt=betaS_ik*betaS_kkp;
                      rfactor=rfactorrt*rfactorrt;

//3rd CC is third term of Eq. 11 (c) for i atom
//where j , k =neighbor of i & k' =neighbor of k

                      CC=CC+gfactor*rfactor;

//agpdpr1 is derivative of CC 3rd term w.r.t. Beta(r_ik)
//agpdpr2 is derivative of CC 3rd term w.r.t. Beta(r_kk')
//app1 is derivative of CC 3rd term w.r.t. cos(theta_jik)
//app2 is derivative of CC 3rd term w.r.t. cos(theta_ikk')

                      agpdpr1=2.0*gfactor*rfactorrt*betaS_kkp
                          *dBetaS_ik/r_ik;
                      agpdpr2=2.0*gfactor*rfactorrt*betaS_ik
                          *dBetaS_kkp/r_kkp;
                      app1=rfactor*gfactorsq2*gsqprime;
                      app2=rfactor*gfactorsq*gsqprime2;
                      bt_sg[nb_ij].dCC[0]+=
                          app1*dcA_jik[0][0];
                      bt_sg[nb_ij].dCC[1]+=
                          app1*dcA_jik[1][0];
                      bt_sg[nb_ij].dCC[2]+=
                          app1*dcA_jik[2][0];
                      bt_sg[nb_ik].dCC[0]+=
                          app1*dcA_jik[0][1]
                          +agpdpr1*dis_ik[0]
                          -app2*dcA_ikkp[0][0];
                      bt_sg[nb_ik].dCC[1]+=
                          app1*dcA_jik[1][1]
                          +agpdpr1*dis_ik[1]
                          -app2*dcA_ikkp[1][0];
                      bt_sg[nb_ik].dCC[2]+=
                          app1*dcA_jik[2][1]
                          +agpdpr1*dis_ik[2]
                          -app2*dcA_ikkp[2][0];
                      bt_sg[nb_kkp].dCC[0]+=
                          app2*dcA_ikkp[0][1]
                          +agpdpr2*dis_kkp[0];
                      bt_sg[nb_kkp].dCC[1]+=
                          app2*dcA_ikkp[1][1]
                          +agpdpr2*dis_kkp[1];
                      bt_sg[nb_kkp].dCC[2]+=
                          app2*dcA_ikkp[2][1]
                          +agpdpr2*dis_kkp[2];
                    }
                  }
                }

//j and k are different neighbors of i and k' is a neighbor j not equal to k

                for(ltmp=0;ltmp<numneigh[j];ltmp++) {
                  sig_flag=0;
                  temp_jkp=BOP_index[j]+ltmp;
                  kp=jlist[ltmp];
                  kptype = map[type[kp]]+1;
                  kplist=firstneigh[kp];

                  same_kkpk=0;
                  same_jkpj=0;

                  for(kpNeij=0;kpNeij<numneigh[kp];kpNeij++) {
                    kpj=kplist[kpNeij];
                    if(x[j][0]==x[kpj][0]) {
                      if(x[j][1]==x[kpj][1]) {
                        if(x[j][2]==x[kpj][2]) {
                          same_jkpj=1;
                          break;
                        }
                      }
                    }
                  }
                  for(kpNeik=0;kpNeik<numneigh[kp];kpNeik++) {
                    kpk=kplist[kpNeik];
                    if(x[k][0]==x[kpk][0]) {
                      if(x[k][1]==x[kpk][1]) {
                        if(x[k][2]==x[kpk][2]) {
                          same_kkpk=1;
                          break;
                        }
                      }
                    }
                  }
                  if(!same_jkpj&&!same_kkpk) {
                    same_kkpk=0;
                    for(kNeikp=0;kNeikp<numneigh[k];kNeikp++) {
                      temp_kkp=BOP_index[k]+kNeikp;
                      kkp=kplist[kNeikp];
                      if(x[kp][0]==x[kkp][0]) {
                        if(x[kp][1]==x[kkp][1]) {
                          if(x[kp][2]==x[kkp][2]) {
                            sig_flag=1;
                            break;
                          }
                        }
                      }
                    }
                    if(sig_flag==1) {
                      for(nsearch=0;nsearch<numneigh[kp];nsearch++) {
                        ncmp=kplist[nsearch];
                        if(x[ncmp][0]==x[j][0]) {
                          if(x[ncmp][1]==x[j][1]) {
                            if(x[ncmp][2]==x[j][2]) {
                              kpNeij=nsearch;
                            }
                          }
                        }
                        if(x[ncmp][0]==x[k][0]) {
                          if(x[ncmp][1]==x[k][1]) {
                            if(x[ncmp][2]==x[k][2]) {
                              kpNeik=nsearch;
                            }
                          }
                        }
                      }
                      if(jtype==kptype)
                        ijkp=jtype-1;
                      else if(jtype<kptype)
                        ijkp=jtype*bop_types-jtype*(jtype+1)/2+kptype-1;
                      else
                        ijkp=kptype*bop_types-kptype*(kptype+1)/2+jtype-1;
                      if(ktype==kptype)
                        ikkp=ktype-1;
                      else if(ktype<kptype)
                        ikkp=ktype*bop_types-ktype*(ktype+1)/2+kptype-1;
                      else
                        ikkp=kptype*bop_types-kptype*(kptype+1)/2+ktype-1;

                      dis_jkp[0]=x[kp][0]-x[j][0];
                      dis_jkp[1]=x[kp][1]-x[j][1];
                      dis_jkp[2]=x[kp][2]-x[j][2];
                      rsq_jkp=dis_jkp[0]*dis_jkp[0]
                          +dis_jkp[1]*dis_jkp[1]
                          +dis_jkp[2]*dis_jkp[2];
                      r_jkp=sqrt(rsq_jkp);
                      ps=r_jkp*rdr[ijkp]+1.0;
                      ks=(int)ps;
                      if(nr-1<ks)
                        ks=nr-1;
                      ps=ps-ks;
                      if(ps>1.0)
                        ps=1.0;
                      betaS_jkp=((pBetaS3[ijkp][ks-1]*ps+pBetaS2[ijkp][ks-1])*ps
                          +pBetaS1[ijkp][ks-1])*ps+pBetaS[ijkp][ks-1];
                      dBetaS_jkp=(pBetaS6[ijkp][ks-1]*ps+pBetaS5[ijkp][ks-1])*ps
                          +pBetaS4[ijkp][ks-1];
                      dis_kkp[0]=x[kp][0]-x[k][0];
                      dis_kkp[1]=x[kp][1]-x[k][1];
                      dis_kkp[2]=x[kp][2]-x[k][2];
                      rsq_kkp=dis_kkp[0]*dis_kkp[0]
                          +dis_kkp[1]*dis_kkp[1]
                          +dis_kkp[2]*dis_kkp[2];
                      r_kkp=sqrt(rsq_kkp);
                      ps=r_kkp*rdr[ikkp]+1.0;
                      ks=(int)ps;
                      if(nr-1<ks)
                        ks=nr-1;
                      ps=ps-ks;
                      if(ps>1.0)
                        ps=1.0;
                      betaS_kkp=((pBetaS3[ikkp][ks-1]*ps+pBetaS2[ikkp][ks-1])*ps
                          +pBetaS1[ikkp][ks-1])*ps+pBetaS[ikkp][ks-1];
                      dBetaS_kkp=(pBetaS6[ikkp][ks-1]*ps+pBetaS5[ikkp][ks-1])*ps
                          +pBetaS4[ikkp][ks-1];
                      cosAng_ijkp=(-dis_ij[0]*dis_jkp[0]-dis_ij[1]*dis_jkp[1]
                          -dis_ij[2]*dis_jkp[2])/(r_ij*r_jkp);
                      dcA_ijkp[0][0]=(dis_jkp[0]*r_ij*r_jkp-cosAng_ijkp
                          *-dis_ij[0]*r_jkp*r_jkp)/(r_ij*r_ij*r_jkp*r_jkp);
                      dcA_ijkp[1][0]=(dis_jkp[1]*r_ij*r_jkp-cosAng_ijkp
                          *-dis_ij[1]*r_jkp*r_jkp)/(r_ij*r_ij*r_jkp*r_jkp);
                      dcA_ijkp[2][0]=(dis_jkp[2]*r_ij*r_jkp-cosAng_ijkp
                          *-dis_ij[2]*r_jkp*r_jkp)/(r_ij*r_ij*r_jkp*r_jkp);
                      dcA_ijkp[0][1]=(-dis_ij[0]*r_ij*r_jkp-cosAng_ijkp
                          *dis_jkp[0]*r_ij*r_ij)/(r_ij*r_ij*r_jkp*r_jkp);
                      dcA_ijkp[1][1]=(-dis_ij[1]*r_ij*r_jkp-cosAng_ijkp
                          *dis_jkp[1]*r_ij*r_ij)/(r_ij*r_ij*r_jkp*r_jkp);
                      dcA_ijkp[2][1]=(-dis_ij[2]*r_ij*r_jkp-cosAng_ijkp
                          *dis_jkp[2]*r_ij*r_ij)/(r_ij*r_ij*r_jkp*r_jkp);
                      cosAng_ikkp=(-dis_ik[0]*dis_kkp[0]-dis_ik[1]*dis_kkp[1]
                          -dis_ik[2]*dis_kkp[2])/(r_ik*r_kkp);
                      dcA_ikkp[0][0]=(dis_kkp[0]*r_ik*r_kkp-cosAng_ikkp
                          *-dis_ik[0]*r_kkp*r_kkp)/(r_ik*r_ik*r_kkp*r_kkp);
                      dcA_ikkp[1][0]=(dis_kkp[1]*r_ik*r_kkp-cosAng_ikkp
                          *-dis_ik[1]*r_kkp*r_kkp)/(r_ik*r_ik*r_kkp*r_kkp);
                      dcA_ikkp[2][0]=(dis_kkp[2]*r_ik*r_kkp-cosAng_ikkp
                          *-dis_ik[2]*r_kkp*r_kkp)/(r_ik*r_ik*r_kkp*r_kkp);
                      dcA_ikkp[0][1]=(-dis_ik[0]*r_ik*r_kkp-cosAng_ikkp
                          *dis_kkp[0]*r_ik*r_ik)/(r_ik*r_ik*r_kkp*r_kkp);
                      dcA_ikkp[1][1]=(-dis_ik[1]*r_ik*r_kkp-cosAng_ikkp
                          *dis_kkp[1]*r_ik*r_ik)/(r_ik*r_ik*r_kkp*r_kkp);
                      dcA_ikkp[2][1]=(-dis_ik[2]*r_ik*r_kkp-cosAng_ikkp
                          *dis_kkp[2]*r_ik*r_ik)/(r_ik*r_ik*r_kkp*r_kkp);
                      cosAng_jkpk=(dis_jkp[0]*dis_kkp[0]+dis_jkp[1]*dis_kkp[1]
                          +dis_jkp[2]*dis_kkp[2])/(r_jkp*r_kkp);
                      dcA_jkpk[0][0]=(-dis_kkp[0]*r_jkp*r_kkp-cosAng_jkpk
                          *-dis_jkp[0]*r_kkp*r_kkp)/(r_jkp*r_jkp*r_kkp*r_kkp);
                      dcA_jkpk[1][0]=(-dis_kkp[1]*r_jkp*r_kkp-cosAng_jkpk
                          *-dis_jkp[1]*r_kkp*r_kkp)/(r_jkp*r_jkp*r_kkp*r_kkp);
                      dcA_jkpk[2][0]=(-dis_kkp[2]*r_jkp*r_kkp-cosAng_jkpk
                          *-dis_jkp[2]*r_kkp*r_kkp)/(r_jkp*r_jkp*r_kkp*r_kkp);
                      dcA_jkpk[0][1]=(-dis_jkp[0]*r_jkp*r_kkp-cosAng_jkpk
                          *-dis_kkp[0]*r_jkp*r_jkp)/(r_jkp*r_jkp*r_kkp*r_kkp);
                      dcA_jkpk[1][1]=(-dis_jkp[1]*r_jkp*r_kkp-cosAng_jkpk
                          *-dis_kkp[1]*r_jkp*r_jkp)/(r_jkp*r_jkp*r_kkp*r_kkp);
                      dcA_jkpk[2][1]=(-dis_jkp[2]*r_jkp*r_kkp-cosAng_jkpk
                          *-dis_kkp[2]*r_jkp*r_jkp)/(r_jkp*r_jkp*r_kkp*r_kkp);
                      sig_flag=0;
                      for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                        ncmp=itypeSigBk[n][nsearch];
                        if(x[ncmp][0]==x[kp][0]) {
                          if(x[ncmp][1]==x[kp][1]) {
                            if(x[ncmp][2]==x[kp][2]) {
                              nkp=nsearch;
                              sig_flag=1;
                              break;
                            }
                          }
                        }
                      }
                      if(sig_flag==0) {
                        nSigBk[n]=nSigBk[n]+1;
                        nkp=nSigBk[n]-1;
                        itypeSigBk[n][nkp]=kp;
                      }
                      nb_jkp=nb_t;
                      nb_t++;
                      if(nb_t>nb_sg) {
                        new_n_tot=nb_sg+maxneigh;
                        grow_sigma(nb_sg,new_n_tot);
                        nb_sg=new_n_tot;
                      }
                      bt_sg[nb_jkp].temp=temp_jkp;
                      bt_sg[nb_jkp].i=j;
                      bt_sg[nb_jkp].j=kp;
                      nb_kkp=nb_t;
                      nb_t++;
                      if(nb_t>nb_sg) {
                        new_n_tot=nb_sg+maxneigh;
                        grow_sigma(nb_sg,new_n_tot);
                        nb_sg=new_n_tot;
                      }
                      bt_sg[nb_kkp].temp=temp_kkp;
                      bt_sg[nb_kkp].i=k;
                      bt_sg[nb_kkp].j=kp;
                      gmean0=sigma_g0[itype-1][jtype-1][kptype-1];
                      gmean1=sigma_g1[itype-1][jtype-1][kptype-1];
                      gmean2=sigma_g2[itype-1][jtype-1][kptype-1];
                      amean=cosAng_ijkp;
                      gfactor2=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime2=gmean1+2.0*gmean2*amean;
                      gmean0=sigma_g0[itype-1][ktype-1][kptype-1];
                      gmean1=sigma_g1[itype-1][ktype-1][kptype-1];
                      gmean2=sigma_g2[itype-1][ktype-1][kptype-1];
                      amean=cosAng_ikkp;
                      gfactor3=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime3=gmean1+2.0*gmean2*amean;
                      gmean0=sigma_g0[jtype-1][kptype-1][ktype-1];
                      gmean1=sigma_g1[jtype-1][kptype-1][ktype-1];
                      gmean2=sigma_g2[jtype-1][kptype-1][ktype-1];
                      amean=cosAng_jkpk;
                      gfactor4=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime4=gmean1+2.0*gmean2*amean;
                      gfactor=gfactor1*gfactor2*gfactor3*gfactor4;
                      rfactor0=(betaS_ik+small2)*(betaS_jkp+small2)
                          *(betaS_kkp+small2);
                      rfactor=pow(rfactor0,2.0/3.0);
                      drfactor=2.0/3.0*pow(rfactor0,-1.0/3.0);

//EE is Eq. 25(notes)

                      EE=EE+gfactor*rfactor;

//agpdpr1 is derivative of agpdpr1 w.r.t. Beta(r_ik)
//agpdpr2 is derivative of agpdpr1 w.r.t. Beta(r_jk')
//agpdpr3 is derivative of agpdpr1 w.r.t. Beta(r_kk')
//app1 is derivative of agpdpr1 w.r.t. cos(theta_jik)
//app2 is derivative of agpdpr1 w.r.t. cos(theta_ijk')
//app3 is derivative of agpdpr1 w.r.t. cos(theta_ikk')
//app4 is derivative of agpdpr1 w.r.t. cos(theta_jk'k)

                      agpdpr1=gfactor*drfactor*(betaS_jkp+small2)*(betaS_kkp
                          +small2)*dBetaS_ik/r_ik;
                      agpdpr2=gfactor*drfactor*(betaS_ik+small2)*(betaS_kkp
                          +small2)*dBetaS_jkp/r_jkp;
                      agpdpr3=gfactor*drfactor*(betaS_ik+small2)*(betaS_jkp
                          +small2)*dBetaS_kkp/r_kkp;
                      app1=rfactor*gfactor2*gfactor3*gfactor4*gprime1;
                      app2=rfactor*gfactor1*gfactor3*gfactor4*gprime2;
                      app3=rfactor*gfactor1*gfactor2*gfactor4*gprime3;
                      app4=rfactor*gfactor1*gfactor2*gfactor3*gprime4;
                      bt_sg[nb_ij].dEE[0]+=
                          app1*dcA_jik[0][0]
                          -app2*dcA_ijkp[0][0];
                      bt_sg[nb_ij].dEE[1]+=
                          app1*dcA_jik[1][0]
                          -app2*dcA_ijkp[1][0];
                      bt_sg[nb_ij].dEE[2]+=
                          app1*dcA_jik[2][0]
                          -app2*dcA_ijkp[2][0];
                      bt_sg[nb_ik].dEE[0]+=
                          app1*dcA_jik[0][1]
                          +agpdpr1*dis_ik[0]
                          -app3*dcA_ikkp[0][0];
                      bt_sg[nb_ik].dEE[1]+=
                          app1*dcA_jik[1][1]
                          +agpdpr1*dis_ik[1]
                          -app3*dcA_ikkp[1][0];
                      bt_sg[nb_ik].dEE[2]+=
                          app1*dcA_jik[2][1]
                          +agpdpr1*dis_ik[2]
                          -app3*dcA_ikkp[2][0];
                      bt_sg[nb_jkp].dEE[0]+=
                          app2*dcA_ijkp[0][1]
                          +agpdpr2*dis_jkp[0]
                          -app4*dcA_jkpk[0][0];
                      bt_sg[nb_jkp].dEE[1]+=
                          app2*dcA_ijkp[1][1]
                          +agpdpr2*dis_jkp[1]
                          -app4*dcA_jkpk[1][0];
                      bt_sg[nb_jkp].dEE[2]+=
                          app2*dcA_ijkp[2][1]
                          +agpdpr2*dis_jkp[2]
                          -app4*dcA_jkpk[2][0];
                      bt_sg[nb_kkp].dEE[0]+=
                          app3*dcA_ikkp[0][1]
                          +agpdpr3*dis_kkp[0]
                          -app4*dcA_jkpk[0][1];
                      bt_sg[nb_kkp].dEE[1]+=
                          app3*dcA_ikkp[1][1]
                          +agpdpr3*dis_kkp[1]
                          -app4*dcA_jkpk[1][1];
                      bt_sg[nb_kkp].dEE[2]+=
                          app3*dcA_ikkp[2][1]
                          +agpdpr3*dis_kkp[2]
                          -app4*dcA_jkpk[2][1];
                    }
                  }
                }
              }
            }
          }

//j is a neighbor of i and k is a neighbor of j not equal to i

          for(ktmp=0;ktmp<numneigh[j];ktmp++) {
            if(ktmp!=ji) {
              temp_jk=BOP_index[j]+ktmp;
              k=jlist[ktmp];
              klist=firstneigh[k];
              ktype=map[type[k]]+1;
              for(kNeij=0;kNeij<numneigh[k];kNeij++) {
                if(x[klist[kNeij]][0]==x[j][0]) {
                  if(x[klist[kNeij]][1]==x[j][1]) {
                    if(x[klist[kNeij]][2]==x[j][2]) {
                      break;
                    }
                  }
                }
              }
              if(jtype==ktype)
                ijk=jtype-1;
              else if(jtype<ktype)
                ijk=jtype*bop_types-jtype*(jtype+1)/2+ktype-1;
              else
                ijk=ktype*bop_types-ktype*(ktype+1)/2+jtype-1;
              sig_flag=0;
              for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                ncmp=itypeSigBk[n][nsearch];
                if(x[ncmp][0]==x[k][0]) {
                  if(x[ncmp][1]==x[k][1]) {
                    if(x[ncmp][2]==x[k][2]) {
                      new1=nsearch;
                      sig_flag=1;
                      break;
                    }
                  }
                }
              }
              if(sig_flag==0) {
                nSigBk[n]=nSigBk[n]+1;
                new1=nSigBk[n]-1;
                itypeSigBk[n][new1]=k;
              }
              dis_jk[0]=x[k][0]-x[j][0];
              dis_jk[1]=x[k][1]-x[j][1];
              dis_jk[2]=x[k][2]-x[j][2];
              rsq_jk=dis_jk[0]*dis_jk[0]
                  +dis_jk[1]*dis_jk[1]
                  +dis_jk[2]*dis_jk[2];
              r_jk=sqrt(rsq_jk);
              if(r_jk<=rcut[ijk]) {
                ps=r_jk*rdr[ijk]+1.0;
                ks=(int)ps;
                if(nr-1<ks)
                  ks=nr-1;
                ps=ps-ks;
                if(ps>1.0)
                  ps=1.0;
                betaS_jk=((pBetaS3[ijk][ks-1]*ps+pBetaS2[ijk][ks-1])*ps
                    +pBetaS1[ijk][ks-1])*ps+pBetaS[ijk][ks-1];
                dBetaS_jk=(pBetaS6[ijk][ks-1]*ps+pBetaS5[ijk][ks-1])*ps
                    +pBetaS4[ijk][ks-1];
                cosAng_ijk=(-dis_ij[0]*dis_jk[0]-dis_ij[1]*dis_jk[1]
                    -dis_ij[2]*dis_jk[2])/(r_ij*r_jk);
                dcA_ijk[0][0]=(dis_jk[0]*r_ij*r_jk-cosAng_ijk
                    *-dis_ij[0]*r_jk*r_jk)/(r_ij*r_ij*r_jk*r_jk);
                dcA_ijk[1][0]=(dis_jk[1]*r_ij*r_jk-cosAng_ijk
                    *-dis_ij[1]*r_jk*r_jk)/(r_ij*r_ij*r_jk*r_jk);
                dcA_ijk[2][0]=(dis_jk[2]*r_ij*r_jk-cosAng_ijk
                    *-dis_ij[2]*r_jk*r_jk)/(r_ij*r_ij*r_jk*r_jk);
                dcA_ijk[0][1]=(-dis_ij[0]*r_ij*r_jk-cosAng_ijk
                    *dis_jk[0]*r_ij*r_ij)/(r_ij*r_ij*r_jk*r_jk);
                dcA_ijk[1][1]=(-dis_ij[1]*r_ij*r_jk-cosAng_ijk
                    *dis_jk[1]*r_ij*r_ij)/(r_ij*r_ij*r_jk*r_jk);
                dcA_ijk[2][1]=(-dis_ij[2]*r_ij*r_jk-cosAng_ijk
                    *dis_jk[2]*r_ij*r_ij)/(r_ij*r_ij*r_jk*r_jk);
                nb_jk=nb_t;
                nb_t++;
                if(nb_t>nb_sg) {
                  new_n_tot=nb_sg+maxneigh;
                  grow_sigma(nb_sg,new_n_tot);
                  nb_sg=new_n_tot;
                }
                bt_sg[nb_jk].temp=temp_jk;
                bt_sg[nb_jk].i=j;
                bt_sg[nb_jk].j=k;
                gmean0=sigma_g0[itype-1][jtype-1][ktype-1];
                gmean1=sigma_g1[itype-1][jtype-1][ktype-1];
                gmean2=sigma_g2[itype-1][jtype-1][ktype-1];
                amean=cosAng_ijk;
                gfactor1=gmean0+gmean1*amean
                    +gmean2*amean*amean;
                gprime1=gmean1+2.0*gmean2*amean;
                gfactorsq=gfactor1*gfactor1;
                gsqprime=2.0*gfactor1*gprime1;
                rfactor1rt=betaS_jk*betaS_jk;
                rfactor1=rfactor1rt*rfactor1rt;

//BB is Eq. 34 (a) or Eq. 10 (c) for the j atom
//1st DD is Eq. 11 (c) for j atom where i & k=neighbor of j

                BB=BB+gfactorsq*rfactor1rt;
                DD=DD+gfactorsq*rfactor1;

//agpdpr1 is derivative of BB  w.r.t. Beta(r_jk)
//app1 is derivative of BB w.r.t. cos(theta_ijk)

                agpdpr1=2.0*gfactorsq*betaS_jk*dBetaS_jk/r_jk;
                agpdpr2=2.0*rfactor1rt*agpdpr1;
                app1=rfactor1rt*gsqprime;
                app2=rfactor1rt*app1;
                bt_sg[nb_ij].dBB[0]-=
                    app1*dcA_ijk[0][0];
                bt_sg[nb_ij].dBB[1]-=
                    app1*dcA_ijk[1][0];
                bt_sg[nb_ij].dBB[2]-=
                    app1*dcA_ijk[2][0];
                bt_sg[nb_ij].dDD[0]-=
                    app2*dcA_ijk[0][0];
                bt_sg[nb_ij].dDD[1]-=
                    app2*dcA_ijk[1][0];
                bt_sg[nb_ij].dDD[2]-=
                    app2*dcA_ijk[2][0];
                bt_sg[nb_jk].dBB[0]+=
                    app1*dcA_ijk[0][1]
                    +agpdpr1*dis_jk[0];
                bt_sg[nb_jk].dBB[1]+=
                    app1*dcA_ijk[1][1]
                    +agpdpr1*dis_jk[1];
                bt_sg[nb_jk].dBB[2]+=
                    app1*dcA_ijk[2][1]
                    +agpdpr1*dis_jk[2];
                bt_sg[nb_jk].dDD[0]+=
                    app2*dcA_ijk[0][1]
                    +agpdpr2*dis_jk[0];
                bt_sg[nb_jk].dDD[1]+=
                    app2*dcA_ijk[1][1]
                    +agpdpr2*dis_jk[1];
                bt_sg[nb_jk].dDD[2]+=
                    app2*dcA_ijk[2][1]
                    +agpdpr2*dis_jk[2];

//j is a neighbor of i, k and k' prime different neighbors of j not equal to i

                for(ltmp=0;ltmp<ktmp;ltmp++) {
                  if(ltmp!=ji) {
                    temp_jkp=BOP_index[j]+ltmp;
                    kp=jlist[ltmp];
                    kptype=map[type[kp]]+1;
                    if(jtype==kptype)
                      ijkp=jtype-1;
                    else if(jtype<kptype)
                      ijkp=jtype*bop_types-jtype*(jtype+1)/2+kptype-1;
                    else
                      ijkp=kptype*bop_types-kptype*(kptype+1)/2+jtype-1;
                    for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                      ncmp=itypeSigBk[n][nsearch];
                      if(x[ncmp][0]==x[kp][0]) {
                        if(x[ncmp][1]==x[kp][1]) {
                          if(x[ncmp][2]==x[kp][2]) {
                            new2=nsearch;
                            break;
                          }
                        }
                      }
                    }
                    dis_jkp[0]=x[kp][0]-x[j][0];
                    dis_jkp[1]=x[kp][1]-x[j][1];
                    dis_jkp[2]=x[kp][2]-x[j][2];
                    rsq_jkp=dis_jkp[0]*dis_jkp[0]
                        +dis_jkp[1]*dis_jkp[1]
                        +dis_jkp[2]*dis_jkp[2];
                    r_jkp=sqrt(rsq_jkp);
                    if(r_jkp<=rcut[ijkp]) {
                      ps=r_jkp*rdr[ijkp]+1.0;
                      ks=(int)ps;
                      if(nr-1<ks)
                        ks=nr-1;
                      ps=ps-ks;
                      if(ps>1.0)
                        ps=1.0;
                      betaS_jkp=((pBetaS3[ijkp][ks-1]*ps+pBetaS2[ijkp][ks-1])*ps
                        +pBetaS1[ijkp][ks-1])*ps+pBetaS[ijkp][ks-1];
                      dBetaS_jkp=(pBetaS6[ijkp][ks-1]*ps+pBetaS5[ijkp][ks-1])*ps
                        +pBetaS4[ijkp][ks-1];
                      cosAng_ijkp=(-dis_ij[0]*dis_jkp[0]-dis_ij[1]*dis_jkp[1]
                        -dis_ij[2]*dis_jkp[2])/(r_ij*r_jkp);
                      dcA_ijkp[0][0]=(dis_jkp[0]*r_ij*r_jkp-cosAng_ijkp
                        *-dis_ij[0]*r_jkp*r_jkp)/(r_ij*r_ij*r_jkp*r_jkp);
                      dcA_ijkp[1][0]=(dis_jkp[1]*r_ij*r_jkp-cosAng_ijkp
                        *-dis_ij[1]*r_jkp*r_jkp)/(r_ij*r_ij*r_jkp*r_jkp);
                      dcA_ijkp[2][0]=(dis_jkp[2]*r_ij*r_jkp-cosAng_ijkp
                        *-dis_ij[2]*r_jkp*r_jkp)/(r_ij*r_ij*r_jkp*r_jkp);
                      dcA_ijkp[0][1]=(-dis_ij[0]*r_ij*r_jkp-cosAng_ijkp
                        *dis_jkp[0]*r_ij*r_ij)/(r_ij*r_ij*r_jkp*r_jkp);
                      dcA_ijkp[1][1]=(-dis_ij[1]*r_ij*r_jkp-cosAng_ijkp
                        *dis_jkp[1]*r_ij*r_ij)/(r_ij*r_ij*r_jkp*r_jkp);
                      dcA_ijkp[2][1]=(-dis_ij[2]*r_ij*r_jkp-cosAng_ijkp
                        *dis_jkp[2]*r_ij*r_ij)/(r_ij*r_ij*r_jkp*r_jkp);
                      cosAng_kjkp=(dis_jk[0]*dis_jkp[0]+dis_jk[1]*dis_jkp[1]
                        +dis_jk[2]*dis_jkp[2])/(r_jk*r_jkp);
                      dcA_kjkp[0][0]=(dis_jkp[0]*r_jk*r_jkp-cosAng_kjkp
                        *dis_jk[0]*r_jkp*r_jkp)/(r_jk*r_jk*r_jkp*r_jkp);
                      dcA_kjkp[1][0]=(dis_jkp[1]*r_jk*r_jkp-cosAng_kjkp
                        *dis_jk[1]*r_jkp*r_jkp)/(r_jk*r_jk*r_jkp*r_jkp);
                      dcA_kjkp[2][0]=(dis_jkp[2]*r_jk*r_jkp-cosAng_kjkp
                        *dis_jk[2]*r_jkp*r_jkp)/(r_jk*r_jk*r_jkp*r_jkp);
                      dcA_kjkp[0][1]=(dis_jk[0]*r_jk*r_jkp-cosAng_kjkp
                        *dis_jkp[0]*r_jk*r_jk)/(r_jk*r_jk*r_jkp*r_jkp);
                      dcA_kjkp[1][1]=(dis_jk[1]*r_jk*r_jkp-cosAng_kjkp
                        *dis_jkp[1]*r_jk*r_jk)/(r_jk*r_jk*r_jkp*r_jkp);
                      dcA_kjkp[2][1]=(dis_jk[2]*r_jk*r_jkp-cosAng_kjkp
                        *dis_jkp[2]*r_jk*r_jk)/(r_jk*r_jk*r_jkp*r_jkp);
                      nb_jkp=nb_t;
                      nb_t++;
                      if(nb_t>nb_sg) {
                        new_n_tot=nb_sg+maxneigh;
                        grow_sigma(nb_sg,new_n_tot);
                        nb_sg=new_n_tot;
                      }
                      bt_sg[nb_jkp].temp=temp_jkp;
                      bt_sg[nb_jkp].i=j;
                      bt_sg[nb_jkp].j=kp;
                      gmean0=sigma_g0[itype-1][jtype-1][kptype-1];
                      gmean1=sigma_g1[itype-1][jtype-1][kptype-1];
                      gmean2=sigma_g2[itype-1][jtype-1][kptype-1];
                      amean=cosAng_ijkp;
                      gfactor2=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime2=gmean1+2.0*gmean2*amean;
                      gmean0=sigma_g0[ktype-1][jtype-1][kptype-1];
                      gmean1=sigma_g1[ktype-1][jtype-1][kptype-1];
                      gmean2=sigma_g2[ktype-1][jtype-1][kptype-1];
                      amean=cosAng_kjkp;
                      gfactor3=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime3=gmean1+2.0*gmean2*amean;
                      gfactor=gfactor1*gfactor2*gfactor3;
                      rfactorrt=betaS_jk*betaS_jkp;
                      rfactor=rfactorrt*rfactorrt;

//2nd DD is Eq. 11 (c) for j atom where i , k & k'=neighbor of j

                      DD=DD+2.0*gfactor*rfactor;

//agpdpr1 is derivative of DD  w.r.t. Beta(r_jk)
//agpdpr2 is derivative of DD  w.r.t. Beta(r_jk')
//app1 is derivative of DD  w.r.t. cos(theta_ijk)
//app2 is derivative of DD  w.r.t. cos(theta_ijkp)
//app3 is derivative of DD  w.r.t. cos(theta_kjkp)

                      agpdpr1=4.0*gfactor*rfactorrt*betaS_jkp
                          *dBetaS_jk/r_jk;
                      agpdpr2=4.0*gfactor*rfactorrt*betaS_jk
                          *dBetaS_jkp/r_jkp;
                      app1=2.0*rfactor*gfactor2*gfactor3*gprime1;
                      app2=2.0*rfactor*gfactor1*gfactor3*gprime2;
                      app3=2.0*rfactor*gfactor1*gfactor2*gprime3;
                      bt_sg[nb_ij].dDD[0]-=
                          app1*dcA_ijk[0][0]
                          +app2*dcA_ijkp[0][0];
                      bt_sg[nb_ij].dDD[1]-=
                          app1*dcA_ijk[1][0]
                          +app2*dcA_ijkp[1][0];
                      bt_sg[nb_ij].dDD[2]-=
                          app1*dcA_ijk[2][0]
                          +app2*dcA_ijkp[2][0];
                      bt_sg[nb_jk].dDD[0]+=
                          app1*dcA_ijk[0][1]
                          +app3*dcA_kjkp[0][0]
                          +agpdpr1*dis_jk[0];
                      bt_sg[nb_jk].dDD[1]+=
                          app1*dcA_ijk[1][1]
                          +app3*dcA_kjkp[1][0]
                          +agpdpr1*dis_jk[1];
                      bt_sg[nb_jk].dDD[2]+=
                          app1*dcA_ijk[2][1]
                          +app3*dcA_kjkp[2][0]
                          +agpdpr1*dis_jk[2];
                      bt_sg[nb_jkp].dDD[0]+=
                          app2*dcA_ijkp[0][1]
                          +app3*dcA_kjkp[0][1]
                          +agpdpr2*dis_jkp[0];
                      bt_sg[nb_jkp].dDD[1]+=
                          app2*dcA_ijkp[1][1]
                          +app3*dcA_kjkp[1][1]
                          +agpdpr2*dis_jkp[1];
                      bt_sg[nb_jkp].dDD[2]+=
                          app2*dcA_ijkp[2][1]
                          +app3*dcA_kjkp[2][1]
                          +agpdpr2*dis_jkp[2];

                    }
                  }
                }

//j is a neighbor of i, k is a neighbor of j not equal to i and k'
//is a neighbor of k not equal to j or i

                for(ltmp=0;ltmp<numneigh[k];ltmp++) {
                  temp_kkp=BOP_index[k]+ltmp;
                  kp=klist[ltmp];
                  kptype=map[type[kp]]+1;
                  same_ikp=0;
                  same_jkp=0;
                  if(x[i][0]==x[kp][0]) {
                    if(x[i][1]==x[kp][1]) {
                      if(x[i][2]==x[kp][2]) {
                        same_ikp=1;
                      }
                    }
                  }
                  if(x[j][0]==x[kp][0]) {
                    if(x[j][1]==x[kp][1]) {
                      if(x[j][2]==x[kp][2]) {
                        same_jkp=1;
                      }
                    }
                  }
                  if(!same_ikp&&!same_jkp) {
                    if(ktype==kptype)
                      ikkp=ktype-1;
                    else if(ktype<kptype)
                      ikkp=ktype*bop_types-ktype*(ktype+1)/2+kptype-1;
                    else
                      ikkp=kptype*bop_types-kptype*(kptype+1)/2+ktype-1;
                    for(kNeij=0;kNeij<numneigh[k];kNeij++) {
                      if(x[klist[kNeij]][0]==x[j][0]) {
                        if(x[klist[kNeij]][1]==x[j][1]) {
                          if(x[klist[kNeij]][2]==x[j][2]) {
                            break;
                          }
                        }
                      }
                    }
                    sig_flag=0;
                    for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                      ncmp=itypeSigBk[n][nsearch];
                      if(x[ncmp][0]==x[kp][0]) {
                        if(x[ncmp][1]==x[kp][1]) {
                          if(x[ncmp][2]==x[kp][2]) {
                            new2=nsearch;
                            sig_flag=1;
                            break;
                          }
                        }
                      }
                    }
                    if(sig_flag==0) {
                      nSigBk[n]=nSigBk[n]+1;
                      new2=nSigBk[n]-1;
                      itypeSigBk[n][new2]=kp;
                    }
                    dis_kkp[0]=x[kp][0]-x[k][0];
                    dis_kkp[1]=x[kp][1]-x[k][1];
                    dis_kkp[2]=x[kp][2]-x[k][2];
                    rsq_kkp=dis_kkp[0]*dis_kkp[0]
                        +dis_kkp[1]*dis_kkp[1]
                        +dis_kkp[2]*dis_kkp[2];
                    r_kkp=sqrt(rsq_kkp);
                    if(r_kkp<=rcut[ikkp]) {
                      ps=r_kkp*rdr[ikkp]+1.0;
                      ks=(int)ps;
                      if(nr-1<ks)
                        ks=nr-1;
                      ps=ps-ks;
                      if(ps>1.0)
                        ps=1.0;
                      betaS_kkp=((pBetaS3[ikkp][ks-1]*ps+pBetaS2[ikkp][ks-1])*ps
                          +pBetaS1[ikkp][ks-1])*ps+pBetaS[ikkp][ks-1];
                      dBetaS_kkp=(pBetaS6[ikkp][ks-1]*ps+pBetaS5[ikkp][ks-1])*ps
                          +pBetaS4[ikkp][ks-1];
                      cosAng_jkkp=(-dis_jk[0]*dis_kkp[0]-dis_jk[1]*dis_kkp[1]
                          -dis_jk[2]*dis_kkp[2])/(r_jk*r_kkp);
                      dcA_jkkp[0][0]=(dis_kkp[0]*r_jk*r_kkp-cosAng_jkkp
                          *-dis_jk[0]*r_kkp*r_kkp)/(r_jk*r_jk*r_kkp*r_kkp);
                      dcA_jkkp[1][0]=(dis_kkp[1]*r_jk*r_kkp-cosAng_jkkp
                          *-dis_jk[1]*r_kkp*r_kkp)/(r_jk*r_jk*r_kkp*r_kkp);
                      dcA_jkkp[2][0]=(dis_kkp[2]*r_jk*r_kkp-cosAng_jkkp
                          *-dis_jk[2]*r_kkp*r_kkp)/(r_jk*r_jk*r_kkp*r_kkp);
                      dcA_jkkp[0][1]=(-dis_jk[0]*r_jk*r_kkp-cosAng_jkkp
                          *dis_kkp[0]*r_jk*r_jk)/(r_jk*r_jk*r_kkp*r_kkp);
                      dcA_jkkp[1][1]=(-dis_jk[1]*r_jk*r_kkp-cosAng_jkkp
                          *dis_kkp[1]*r_jk*r_jk)/(r_jk*r_jk*r_kkp*r_kkp);
                      dcA_jkkp[2][1]=(-dis_jk[2]*r_jk*r_kkp-cosAng_jkkp
                          *dis_kkp[2]*r_jk*r_jk)/(r_jk*r_jk*r_kkp*r_kkp);
                      nb_kkp=nb_t;
                      nb_t++;
                      if(nb_t>nb_sg) {
                        new_n_tot=nb_sg+maxneigh;
                        grow_sigma(nb_sg,new_n_tot);
                        nb_sg=new_n_tot;
                      }
                      bt_sg[nb_kkp].temp=temp_kkp;
                      bt_sg[nb_kkp].i=k;
                      bt_sg[nb_kkp].j=kp;
                      gmean0=sigma_g0[jtype-1][ktype-1][kptype-1];
                      gmean1=sigma_g1[jtype-1][ktype-1][kptype-1];
                      gmean2=sigma_g2[jtype-1][ktype-1][kptype-1];
                      amean=cosAng_jkkp;
                      gfactor2=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime2=gmean1+2.0*gmean2*amean;
                      gfactorsq2=gfactor2*gfactor2;
                      gsqprime2=2.0*gfactor2*gprime2;
                      gfactor=gfactorsq*gfactorsq2;
                      rfactorrt=betaS_jk*betaS_kkp;
                      rfactor=rfactorrt*rfactorrt;

//3rd DD is Eq. 11 (c) for j atom where i & k=neighbor of j & k'=neighbor of k

                      DD=DD+gfactor*rfactor;

//agpdpr1 is derivative of DD  3rd term w.r.t. Beta(r_jk)
//agpdpr2 is derivative of DD  3rd term w.r.t. Beta(r_kk')
//app1 is derivative of DD  3rd term w.r.t. cos(theta_ijk)
//app2 is derivative of DD  3rd term w.r.t. cos(theta_jkkp)

                      agpdpr1=2.0*gfactor*rfactorrt*betaS_kkp
                          *dBetaS_jk/r_jk;
                      agpdpr2=2.0*gfactor*rfactorrt*betaS_jk
                          *dBetaS_kkp/r_kkp;
                      app1=rfactor*gfactorsq2*gsqprime;
                      app2=rfactor*gfactorsq*gsqprime2;
                      bt_sg[nb_ij].dDD[0]-=
                          app1*dcA_ijk[0][0];
                      bt_sg[nb_ij].dDD[1]-=
                          app1*dcA_ijk[1][0];
                      bt_sg[nb_ij].dDD[2]-=
                          app1*dcA_ijk[2][0];
                      bt_sg[nb_jk].dDD[0]+=
                          app1*dcA_ijk[0][1]
                          +agpdpr1*dis_jk[0]
                          -app2*dcA_jkkp[0][0];
                      bt_sg[nb_jk].dDD[1]+=
                          app1*dcA_ijk[1][1]
                          +agpdpr1*dis_jk[1]
                          -app2*dcA_jkkp[1][0];
                      bt_sg[nb_jk].dDD[2]+=
                          app1*dcA_ijk[2][1]
                          +agpdpr1*dis_jk[2]
                          -app2*dcA_jkkp[2][0];
                      bt_sg[nb_kkp].dDD[0]+=
                          app2*dcA_jkkp[0][1]
                          +agpdpr2*dis_kkp[0];
                      bt_sg[nb_kkp].dDD[1]+=
                          app2*dcA_jkkp[1][1]
                          +agpdpr2*dis_kkp[1];
                      bt_sg[nb_kkp].dDD[2]+=
                          app2*dcA_jkkp[2][1]
                          +agpdpr2*dis_kkp[2];

                    }
                  }
                }
              }
            }
          }

          sig_flag=0;
          if(FF<=0.000001) {
            sigB[n]=0.0;
            sig_flag=1;
          }
          if(sig_flag==0) {
            if(AA<0.0)
              AA=0.0;
            if(BB<0.0)
              BB=0.0;
            if(CC<0.0)
              CC=0.0;
            if(DD<0.0)
              DD=0.0;

// AA and BB are the representations of (a) Eq. 34 and (b) Eq. 9
// for atoms i and j respectively

            AAC=AA+BB;
            BBC=AA*BB;
            CCC=AA*AA+BB*BB;
            DDC=CC+DD;

//EEC is a modified form of (a) Eq. 33

            EEC=(DDC-CCC)/(AAC+2.0*small1);
            AACFF=1.0/(AAC+2.0*small1);
            for(m=0;m<nb_t;m++) {
              if((bt_sg[m].i>-1)&&(bt_sg[m].j>-1)) {
                bt_sg[m].dAAC[0]=bt_sg[m].dAA[0]
                    +bt_sg[m].dBB[0];
                bt_sg[m].dAAC[1]=bt_sg[m].dAA[1]
                    +bt_sg[m].dBB[1];
                bt_sg[m].dAAC[2]=bt_sg[m].dAA[2]
                    +bt_sg[m].dBB[2];
                bt_sg[m].dBBC[0]=bt_sg[m].dAA[0]*BB
                    +AA*bt_sg[m].dBB[0];
                bt_sg[m].dBBC[1]=bt_sg[m].dAA[1]*BB
                    +AA*bt_sg[m].dBB[1];
                bt_sg[m].dBBC[2]=bt_sg[m].dAA[2]*BB
                    +AA*bt_sg[m].dBB[2];
                bt_sg[m].dCCC[0]=2.0*AA*bt_sg[m].dAA[0]
                    +2.0*BB*bt_sg[m].dBB[0];
                bt_sg[m].dCCC[1]=2.0*AA*bt_sg[m].dAA[1]
                    +2.0*BB*bt_sg[m].dBB[1];
                bt_sg[m].dCCC[2]=2.0*AA*bt_sg[m].dAA[2]
                    +2.0*BB*bt_sg[m].dBB[2];
                bt_sg[m].dDDC[0]=bt_sg[m].dCC[0]
                    +bt_sg[m].dDD[0];
                bt_sg[m].dDDC[1]=bt_sg[m].dCC[1]
                    +bt_sg[m].dDD[1];
                bt_sg[m].dDDC[2]=bt_sg[m].dCC[2]
                    +bt_sg[m].dDD[2];
                bt_sg[m].dEEC[0]=(bt_sg[m].dDDC[0]
                    -bt_sg[m].dCCC[0]
                    -EEC*bt_sg[m].dAAC[0])*AACFF;
                bt_sg[m].dEEC[1]=(bt_sg[m].dDDC[1]
                    -bt_sg[m].dCCC[1]
                    -EEC*bt_sg[m].dAAC[1])*AACFF;
                bt_sg[m].dEEC[2]=(bt_sg[m].dDDC[2]
                    -bt_sg[m].dCCC[2]
                    -EEC*bt_sg[m].dAAC[2])*AACFF;
              }
            }
            UT=EEC*FF+BBC+small3[iij];
            UT=1.0/sqrt(UT);

// FFC is slightly modified form of (a) Eq. 31
// GGC is slightly modified form of (a) Eq. 32
// bndtmp is a slightly modified form of (a) Eq. 30 and (b) Eq. 8

            FFC=BBC*UT;
            GGC=EEC*UT;
            bndtmp=(FF+sigma_delta[iij]*sigma_delta[iij])*(1.0+sigma_a[iij]*GGC)
                *(1.0+sigma_a[iij]*GGC)+sigma_c[iij]*(AAC+sigma_a[iij]*EE
                +sigma_a[iij]*FFC*(2.0+GGC))+small4;
            UTcom=-0.5*UT*UT*UT;
            for(m=0;m<nb_t;m++) {
              if((bt_sg[m].i>-1)&&(bt_sg[m].j>-1)) {
                bt_sg[m].dUT[0]=UTcom*(bt_sg[m].dEEC[0]*FF
                    +EEC*bt_sg[m].dFF[0]+bt_sg[m].dBBC[0]);
                bt_sg[m].dUT[1]=UTcom*(bt_sg[m].dEEC[1]*FF
                    +EEC*bt_sg[m].dFF[1]+bt_sg[m].dBBC[1]);
                bt_sg[m].dUT[2]=UTcom*(bt_sg[m].dEEC[2]*FF
                    +EEC*bt_sg[m].dFF[2]+bt_sg[m].dBBC[2]);
                bt_sg[m].dFFC[0]=bt_sg[m].dBBC[0]*UT
                    +BBC*bt_sg[m].dUT[0];
                bt_sg[m].dFFC[1]=bt_sg[m].dBBC[1]*UT
                    +BBC*bt_sg[m].dUT[1];
                bt_sg[m].dFFC[2]=bt_sg[m].dBBC[2]*UT
                    +BBC*bt_sg[m].dUT[2];
                bt_sg[m].dGGC[0]=bt_sg[m].dEEC[0]*UT
                    +EEC*bt_sg[m].dUT[0];
                bt_sg[m].dGGC[1]=bt_sg[m].dEEC[1]*UT
                    +EEC*bt_sg[m].dUT[1];
                bt_sg[m].dGGC[2]=bt_sg[m].dEEC[2]*UT
                    +EEC*bt_sg[m].dUT[2];
              }
            }
            psign=1.0;
            if(1.0+sigma_a[iij]*GGC<0.0)
              psign=-1.0;
            bndtmp0=1.0/sqrt(bndtmp);
            sigB1[n]=psign*betaS_ij*(1.0+sigma_a[iij]*GGC)*bndtmp0;
            bndtmp=-0.5*bndtmp0*bndtmp0*bndtmp0;
            bndtmp1=psign*(1.0+sigma_a[iij]*GGC)*bndtmp0+psign*betaS_ij
                *(1.0+sigma_a[iij]*GGC)*bndtmp*2.0*betaS_ij*(1.0
                +sigma_a[iij]*GGC)*(1.0+sigma_a[iij]*GGC);
            bndtmp1=bndtmp1*dBetaS_ij/r_ij;
            bndtmp2=psign*betaS_ij*(1.0+sigma_a[iij]*GGC)*bndtmp*sigma_c[iij];
            bndtmp3=psign*betaS_ij*(1.0+sigma_a[iij]*GGC)
                *bndtmp*sigma_c[iij]*sigma_a[iij];
            bndtmp4=psign*betaS_ij*(1.0+sigma_a[iij]*GGC)
                *bndtmp*sigma_c[iij]*sigma_a[iij]*(2.0+GGC);
            bndtmp5=sigma_a[iij]*psign*betaS_ij*bndtmp0
                +psign*betaS_ij*(1.0+sigma_a[iij]*GGC)*bndtmp
                *(2.0*(FF+sigma_delta[iij]*sigma_delta[iij])*(1.0
                +sigma_a[iij]*GGC)*sigma_a[iij]+sigma_c[iij]*sigma_a[iij]*FFC);
            setting=0;
            for(m=0;m<nb_t;m++) {
              if((bt_sg[m].i>-1)&&(bt_sg[m].j>-1)) {
                temp_kk=bt_sg[m].temp;
                if(temp_kk==temp_ij&&setting==0) {
                  bt_sg[m].dSigB1[0]=bndtmp1*dis_ij[0]
                      +(bndtmp2*bt_sg[m].dAAC[0]
                      +bndtmp3*bt_sg[m].dEE[0]
                      +bndtmp4*bt_sg[m].dFFC[0]
                      +bndtmp5*bt_sg[m].dGGC[0]);
                  bt_sg[m].dSigB1[1]=bndtmp1*dis_ij[1]
                      +(bndtmp2*bt_sg[m].dAAC[1]
                      +bndtmp3*bt_sg[m].dEE[1]
                      +bndtmp4*bt_sg[m].dFFC[1]
                      +bndtmp5*bt_sg[m].dGGC[1]);
                  bt_sg[m].dSigB1[2]=bndtmp1*dis_ij[2]
                      +(bndtmp2*bt_sg[m].dAAC[2]
                      +bndtmp3*bt_sg[m].dEE[2]
                      +bndtmp4*bt_sg[m].dFFC[2]
                      +bndtmp5*bt_sg[m].dGGC[2]);
                  setting=1;
                }
                else if(temp_kk==temp_ji&&setting==0) {
                  bt_sg[m].dSigB1[0]=-bndtmp1*dis_ij[0]
                      +(bndtmp2*bt_sg[m].dAAC[0]
                      +bndtmp3*bt_sg[m].dEE[0]
                      +bndtmp4*bt_sg[m].dFFC[0]
                      +bndtmp5*bt_sg[m].dGGC[0]);
                  bt_sg[m].dSigB1[1]=-bndtmp1*dis_ij[1]
                      +(bndtmp2*bt_sg[m].dAAC[1]
                      +bndtmp3*bt_sg[m].dEE[1]
                      +bndtmp4*bt_sg[m].dFFC[1]
                      +bndtmp5*bt_sg[m].dGGC[1]);
                  bt_sg[m].dSigB1[2]=-bndtmp1*dis_ij[2]
                      +(bndtmp2*bt_sg[m].dAAC[2]
                      +bndtmp3*bt_sg[m].dEE[2]
                      +bndtmp4*bt_sg[m].dFFC[2]
                      +bndtmp5*bt_sg[m].dGGC[2]);
                  setting=1;
                }
                else {
                  bt_sg[m].dSigB1[0]=(bndtmp2*bt_sg[m].dAAC[0]
                      +bndtmp3*bt_sg[m].dEE[0]
                      +bndtmp4*bt_sg[m].dFFC[0]
                      +bndtmp5*bt_sg[m].dGGC[0]);
                  bt_sg[m].dSigB1[1]=(bndtmp2*bt_sg[m].dAAC[1]
                      +bndtmp3*bt_sg[m].dEE[1]
                      +bndtmp4*bt_sg[m].dFFC[1]
                      +bndtmp5*bt_sg[m].dGGC[1]);
                  bt_sg[m].dSigB1[2]=(bndtmp2*bt_sg[m].dAAC[2]
                      +bndtmp3*bt_sg[m].dEE[2]
                      +bndtmp4*bt_sg[m].dFFC[2]
                      +bndtmp5*bt_sg[m].dGGC[2]);
                }
              }
            }

//This loop is to ensure there is not an error for atoms with no neighbors (deposition)

            if(nb_t==0) {
              if(j>i) {
                bt_sg[0].dSigB1[0]=bndtmp1*dis_ij[0];
                bt_sg[0].dSigB1[1]=bndtmp1*dis_ij[1];
                bt_sg[0].dSigB1[2]=bndtmp1*dis_ij[2];
              }
              else {
                bt_sg[0].dSigB1[0]=-bndtmp1*dis_ij[0];
                bt_sg[0].dSigB1[1]=-bndtmp1*dis_ij[1];
                bt_sg[0].dSigB1[2]=-bndtmp1*dis_ij[2];
              }
              for(pp=0;pp<3;pp++) {
                bt_sg[0].dAA[pp]=0.0;
                bt_sg[0].dBB[pp]=0.0;
                bt_sg[0].dCC[pp]=0.0;
                bt_sg[0].dDD[pp]=0.0;
                bt_sg[0].dEE[pp]=0.0;
                bt_sg[0].dEE1[pp]=0.0;
                bt_sg[0].dFF[pp]=0.0;
                bt_sg[0].dAAC[pp]=0.0;
                bt_sg[0].dBBC[pp]=0.0;
                bt_sg[0].dCCC[pp]=0.0;
                bt_sg[0].dDDC[pp]=0.0;
                bt_sg[0].dEEC[pp]=0.0;
                bt_sg[0].dFFC[pp]=0.0;
                bt_sg[0].dGGC[pp]=0.0;
                bt_sg[0].dUT[pp]=0.0;
                bt_sg[0].dSigB1[pp]=0.0;
                bt_sg[0].dSigB[pp]=0.0;
              }
              bt_sg[0].i=i;
              bt_sg[0].j=j;
              bt_sg[0].temp=temp_ij;
              nb_t++;
              if(nb_t>nb_sg) {
                new_n_tot=nb_sg+maxneigh;
                grow_sigma(nb_sg,new_n_tot);
                nb_sg=new_n_tot;
              }
            }
            ps=sigB1[n]*rdBO+1.0;
            ks=(int)ps;
            if(nBOt-1<ks)
              ks=nBOt-1;
            ps=ps-ks;
            if(ps>1.0)
              ps=1.0;
            dsigB1=((FsigBO3[iij][ks-1]*ps+FsigBO2[iij][ks-1])*ps
                +FsigBO1[iij][ks-1])*ps+FsigBO[iij][ks-1];
            dsigB2=(FsigBO6[iij][ks-1]*ps+FsigBO5[iij][ks-1])*ps+FsigBO4[iij][ks-1];
            part0=(FF+0.5*AAC+small5);
            part1=(sigma_f[iij]-0.5)*sigma_k[iij];
            part2=1.0-part1*EE1/part0;
            part3=dsigB1*part1/part0;
            part4=part3/part0*EE1;

// sigB is the final expression for (a) Eq. 6 and (b) Eq. 11

            sigB[n]=dsigB1*part2;
            pp1=2.0*betaS_ij;
            for(m=0;m<nb_t;m++) {
              if((bt_sg[m].i>-1)&&(bt_sg[m].j>-1)) {
                temp_kk=bt_sg[m].temp;
                bt_i=bt_sg[m].i;
                bt_j=bt_sg[m].j;
                xtmp[0]=x[bt_j][0]-x[bt_i][0];
                xtmp[1]=x[bt_j][1]-x[bt_i][1];
                xtmp[2]=x[bt_j][2]-x[bt_i][2];
                for(pp=0;pp<3;pp++) {
                  bt_sg[m].dSigB[pp]=dsigB2*part2*bt_sg[m].dSigB1[pp]
                      -part3*bt_sg[m].dEE1[pp]
                      +part4*(bt_sg[m].dFF[pp]
                      +0.5*bt_sg[m].dAAC[pp]);
                }
                for(pp=0;pp<3;pp++) {
                  ftmp[pp]=pp1*bt_sg[m].dSigB[pp];
                  f[bt_i][pp]-=ftmp[pp];
                  f[bt_j][pp]+=ftmp[pp];
                }
                if(evflag) {
                  ev_tally_xyz(bt_i,bt_j,nlocal,newton_pair,0.0,0.0,ftmp[0],ftmp[1]
                      ,ftmp[2],xtmp[0],xtmp[1],xtmp[2]);
                }
              }
            }
          }
          n++;
        }
      }
    }
  }
  destroy_sigma();
}

/* ---------------------------------------------------------------------- */

/*  The formulation differs slightly to avoid negative square roots
    in the calculation of Theta_pi,ij of (a) Eq. 36 and (b) Eq. 18
    see (d) */

void PairBOP::sigmaBo_noa_otf()
{
  int nb_t,new_n_tot;
  int n,i,j,k,kp,m,pp;
  int itmp,jtmp,ktmp,ltmp,mtmp;
  tagint i_tag,j_tag;
  int kp1,kp2,kp1type;
  int iij,iik,ijk,ikkp,ji,iikp,ijkp;
  int nkp;
  int nk0;
  int jNeik,kNeii,kNeij;
  int new1,new2,nlocal;
  int inum,*ilist,*iilist,*jlist,*klist;
  int **firstneigh,*numneigh;
  int temp_ij,temp_ik,temp_jkp,temp_kk,temp_jk;
  int temp_ji,temp_kkp;
  int nb_ij,nb_ik;
  int nb_jk,nb_jkp,nb_kkp;
  int nsearch;
  int sig_flag,setting,ncmp,ks;
  int itype,jtype,ktype,kptype;
  int bt_i,bt_j;
  int same_ikp,same_jkp,same_kpk;
  double AA,BB,CC,DD,EE1,FF;
  double AAC,BBC,CCC,DDC,EEC;
  double UT,bndtmp;
  double amean,gmean0,gmean1,gmean2,ps;
  double gfactor1,gprime1,gsqprime;
  double gfactorsq,gfactor2,gprime2;
  double gfactorsq2;
  double gfactor3,gprime3,gfactor,rfactor;
  double rfactorrt,rfactor1rt,rfactor1;
  double rcm1,rcm2,gcm1,gcm2,gcm3;
  double agpdpr1,app1;
  double dsigB1,dsigB2;
  double part0,part1,part2,part3,part4;
  double psign,bndtmp0,pp1;
  double bndtmp1,bndtmp2;
  double dis_ij[3],rsq_ij,r_ij;
  double betaS_ij,dBetaS_ij;
  double dis_ik[3],rsq_ik,r_ik;
  double betaS_ik,dBetaS_ik;
  double dis_ikp[3],rsq_ikp,r_ikp;
  double betaS_ikp;
  double dis_jk[3],rsq_jk,r_jk;
  double betaS_jk,dBetaS_jk;
  double dis_jkp[3],rsq_jkp,r_jkp;
  double betaS_jkp,dBetaS_jkp;
  double dis_kkp[3],rsq_kkp,r_kkp;
  double betaS_kkp;
  double cosAng_jik,dcA_jik[3][2];
  double cosAng_jikp;
  double cosAng_kikp;
  double cosAng_ijk,dcA_ijk[3][2];
  double cosAng_ijkp;
  double cosAng_kjkp;
  double cosAng_ikj,dcA_ikj[3][2];
  double cosAng_ikkp;
  double cosAng_jkkp;


  double ftmp[3],xtmp[3];
  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int newton_pair = force->newton_pair;
  int *type = atom->type;

  nlocal = atom->nlocal;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  n=0;
  if(nb_sg==0) {
    nb_sg=4;
  }
  if(allocate_sigma) {
    destroy_sigma();
  }
    create_sigma(nb_sg);
  for(itmp=0;itmp<inum;itmp++) {

    i = ilist[itmp];
    i_tag=tag[i];
    itype = map[type[i]]+1;

//j is loop over all neighbors of i

    for(jtmp=0;jtmp<numneigh[i];jtmp++) {
      for(m=0;m<nb_sg;m++) {
        for(pp=0;pp<3;pp++) {
          bt_sg[m].dAA[pp]=0.0;
          bt_sg[m].dBB[pp]=0.0;
          bt_sg[m].dEE1[pp]=0.0;
          bt_sg[m].dFF[pp]=0.0;
          bt_sg[m].dAAC[pp]=0.0;
          bt_sg[m].dSigB1[pp]=0.0;
          bt_sg[m].dSigB[pp]=0.0;
        }
        bt_sg[m].i=-1;
        bt_sg[m].j=-1;
      }
      nb_t=0;
      iilist=firstneigh[i];
      temp_ij=BOP_index[i]+jtmp;
      j=iilist[jtmp];
      jlist=firstneigh[j];
      j_tag=tag[j];
      jtype = map[type[j]]+1;
      nb_ij=nb_t;
      nb_t++;
      if(nb_t>nb_sg) {
        new_n_tot=nb_sg+maxneigh;
        grow_sigma(nb_sg,new_n_tot);
        nb_sg=new_n_tot;
      }
      bt_sg[nb_ij].temp=temp_ij;
      bt_sg[nb_ij].i=i;
      bt_sg[nb_ij].j=j;
      if(j_tag>=i_tag) {
        if(itype==jtype)
          iij=itype-1;
        else if(itype<jtype)
          iij=itype*bop_types-itype*(itype+1)/2+jtype-1;
        else
          iij=jtype*bop_types-jtype*(jtype+1)/2+itype-1;
        for(ji=0;ji<numneigh[j];ji++) {
          temp_ji=BOP_index[j]+ji;
          if(x[jlist[ji]][0]==x[i][0]) {
            if(x[jlist[ji]][1]==x[i][1]) {
              if(x[jlist[ji]][2]==x[i][2]) {
                break;
              }
            }
          }
        }
        dis_ij[0]=x[j][0]-x[i][0];
        dis_ij[1]=x[j][1]-x[i][1];
        dis_ij[2]=x[j][2]-x[i][2];
        rsq_ij=dis_ij[0]*dis_ij[0]
            +dis_ij[1]*dis_ij[1]
            +dis_ij[2]*dis_ij[2];
        r_ij=sqrt(rsq_ij);

        if(r_ij<rcut[iij]) {

          ps=r_ij*rdr[iij]+1.0;
          ks=(int)ps;
          if(nr-1<ks)
            ks=nr-1;
          ps=ps-ks;
          if(ps>1.0)
            ps=1.0;
          betaS_ij=((pBetaS3[iij][ks-1]*ps+pBetaS2[iij][ks-1])*ps
              +pBetaS1[iij][ks-1])*ps+pBetaS[iij][ks-1];
          dBetaS_ij=(pBetaS6[iij][ks-1]*ps+pBetaS5[iij][ks-1])*ps
              +pBetaS4[iij][ks-1];
          nSigBk[n]=0;

//AA-EE1 are the components making up Eq. 30 (a)

          AA=0.0;
          BB=0.0;
          CC=0.0;
          DD=0.0;
          EE1=0.0;

//FF is the Beta_sigma^2 term

          FF=betaS_ij*betaS_ij;

//agpdpr1 is derivative of FF w.r.t. r_ij

          agpdpr1=2.0*betaS_ij*dBetaS_ij/r_ij;

//dXX derivatives are taken with respect to all pairs contributing to the energy
//nb_ij is derivative w.r.t. ij pair

          bt_sg[nb_ij].dFF[0]=agpdpr1*dis_ij[0];
          bt_sg[nb_ij].dFF[1]=agpdpr1*dis_ij[1];
          bt_sg[nb_ij].dFF[2]=agpdpr1*dis_ij[2];

//k is loop over all neighbors of i again with j neighbor of i

          for(ktmp=0;ktmp<numneigh[i];ktmp++) {
            temp_ik=BOP_index[i]+ktmp;
            if(ktmp!=jtmp) {
              k=iilist[ktmp];
              klist=firstneigh[k];
              ktype = map[type[k]]+1;
              if(itype==ktype)
                iik=itype-1;
              else if(itype<ktype)
                iik=itype*bop_types-itype*(itype+1)/2+ktype-1;
              else
                iik=ktype*bop_types-ktype*(ktype+1)/2+itype-1;

//find neighbor of k that is equal to i

              for(kNeii=0;kNeii<numneigh[k];kNeii++) {
                if(x[klist[kNeii]][0]==x[i][0]) {
                  if(x[klist[kNeii]][1]==x[i][1]) {
                    if(x[klist[kNeii]][2]==x[i][2]) {
                      break;
                    }
                  }
                }
              }
              dis_ik[0]=x[k][0]-x[i][0];
              dis_ik[1]=x[k][1]-x[i][1];
              dis_ik[2]=x[k][2]-x[i][2];
              rsq_ik=dis_ik[0]*dis_ik[0]
                  +dis_ik[1]*dis_ik[1]
                  +dis_ik[2]*dis_ik[2];
              r_ik=sqrt(rsq_ik);
              if(r_ik<=rcut[iik]) {
                ps=r_ik*rdr[iik]+1.0;
                ks=(int)ps;
                if(nr-1<ks)
                  ks=nr-1;
                ps=ps-ks;
                if(ps>1.0)
                  ps=1.0;
                betaS_ik=((pBetaS3[iik][ks-1]*ps+pBetaS2[iik][ks-1])*ps
                    +pBetaS1[iik][ks-1])*ps+pBetaS[iik][ks-1];
                dBetaS_ik=(pBetaS6[iik][ks-1]*ps+pBetaS5[iik][ks-1])*ps
                    +pBetaS4[iik][ks-1];

//find neighbor of i that is equal to k

                for(jNeik=0;jNeik<numneigh[j];jNeik++) {
                  temp_jk=BOP_index[j]+jNeik;
                  if(x[jlist[jNeik]][0]==x[k][0]) {
                    if(x[jlist[jNeik]][1]==x[k][1]) {
                      if(x[jlist[jNeik]][2]==x[k][2]) {
                        break;
                      }
                    }
                  }
                }

//find neighbor of k that is equal to j

                for(kNeij=0;kNeij<numneigh[k];kNeij++) {
                  if(x[klist[kNeij]][0]==x[j][0]) {
                    if(x[klist[kNeij]][1]==x[j][1]) {
                      if(x[klist[kNeij]][2]==x[j][2]) {
                        break;
                      }
                    }
                  }
                }
                dis_jk[0]=x[k][0]-x[j][0];
                dis_jk[1]=x[k][1]-x[j][1];
                dis_jk[2]=x[k][2]-x[j][2];
                rsq_jk=dis_jk[0]*dis_jk[0]
                    +dis_jk[1]*dis_jk[1]
                    +dis_jk[2]*dis_jk[2];
                r_jk=sqrt(rsq_jk);

                sig_flag=0;
                for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                  ncmp=itypeSigBk[n][nsearch];
                  if(x[ncmp][0]==x[k][0]) {
                    if(x[ncmp][1]==x[k][1]) {
                      if(x[ncmp][2]==x[k][2]) {
                        nk0=nsearch;
                        sig_flag=1;
                        break;
                      }
                    }
                  }
                }
                if(sig_flag==0) {
                  nSigBk[n]=nSigBk[n]+1;
                  nk0=nSigBk[n]-1;
                  itypeSigBk[n][nk0]=k;
                }
                nb_ik=nb_t;
                nb_t++;
                if(nb_t>nb_sg) {
                  new_n_tot=nb_sg+maxneigh;
                  grow_sigma(nb_sg,new_n_tot);
                  nb_sg=new_n_tot;
                }
                bt_sg[nb_ik].temp=temp_ik;
                bt_sg[nb_ik].i=i;
                bt_sg[nb_ik].j=k;
                nb_jk=nb_t;
                nb_t++;
                if(nb_t>nb_sg) {
                  new_n_tot=nb_sg+maxneigh;
                  grow_sigma(nb_sg,new_n_tot);
                  nb_sg=new_n_tot;
                }
                bt_sg[nb_jk].temp=temp_jk;
                bt_sg[nb_jk].i=j;
                bt_sg[nb_jk].j=k;
                cosAng_jik=(dis_ij[0]*dis_ik[0]+dis_ij[1]*dis_ik[1]
                    +dis_ij[2]*dis_ik[2])/(r_ij*r_ik);
                dcA_jik[0][0]=(dis_ik[0]*r_ij*r_ik-cosAng_jik
                    *dis_ij[0]*r_ik*r_ik)/(r_ij*r_ij*r_ik*r_ik);
                dcA_jik[1][0]=(dis_ik[1]*r_ij*r_ik-cosAng_jik
                    *dis_ij[1]*r_ik*r_ik)/(r_ij*r_ij*r_ik*r_ik);
                dcA_jik[2][0]=(dis_ik[2]*r_ij*r_ik-cosAng_jik
                    *dis_ij[2]*r_ik*r_ik)/(r_ij*r_ij*r_ik*r_ik);
                dcA_jik[0][1]=(dis_ij[0]*r_ij*r_ik-cosAng_jik
                    *dis_ik[0]*r_ij*r_ij)/(r_ij*r_ij*r_ik*r_ik);
                dcA_jik[1][1]=(dis_ij[1]*r_ij*r_ik-cosAng_jik
                    *dis_ik[1]*r_ij*r_ij)/(r_ij*r_ij*r_ik*r_ik);
                dcA_jik[2][1]=(dis_ij[2]*r_ij*r_ik-cosAng_jik
                    *dis_ik[2]*r_ij*r_ij)/(r_ij*r_ij*r_ik*r_ik);
                gmean0=sigma_g0[jtype-1][itype-1][ktype-1];
                gmean1=sigma_g1[jtype-1][itype-1][ktype-1];
                gmean2=sigma_g2[jtype-1][itype-1][ktype-1];
                amean=cosAng_jik;
                gfactor1=gmean0+gmean1*amean
                    +gmean2*amean*amean;
                gfactorsq=gfactor1*gfactor1;
                gprime1=gmean1+2.0*gmean2*amean;
                gsqprime=2.0*gfactor1*gprime1;

//AA is Eq. 34 (a) or Eq. 10 (c) for the i atom
//1st CC is Eq. 11 (c) for i atom where j & k=neighbor of i

                AA=AA+gfactorsq*betaS_ik*betaS_ik;
                CC=CC+gfactorsq*betaS_ik*betaS_ik*betaS_ik*betaS_ik;

//agpdpr1 is derivative of AA w.r.t. Beta(rik)
//app1 is derivative of AA w.r.t. cos(theta_jik)

                agpdpr1=2.0*gfactorsq*betaS_ik*dBetaS_ik/r_ik;
                app1=betaS_ik*betaS_ik*gsqprime;
                bt_sg[nb_ij].dAA[0]+=
                    app1*dcA_jik[0][0];
                bt_sg[nb_ij].dAA[1]+=
                    app1*dcA_jik[1][0];
                bt_sg[nb_ij].dAA[2]+=
                    app1*dcA_jik[2][0];
                bt_sg[nb_ik].dAA[0]+=
                    app1*dcA_jik[0][1]
                    +agpdpr1*dis_ik[0];
                bt_sg[nb_ik].dAA[1]+=
                    app1*dcA_jik[1][1]
                    +agpdpr1*dis_ik[1];
                bt_sg[nb_ik].dAA[2]+=
                    app1*dcA_jik[2][1]
                    +agpdpr1*dis_ik[2];

//k' is loop over neighbors all neighbors of j with k a neighbor
//of i and j a neighbor of i and determine which k' is k

                same_kpk=0;
                for(ltmp=0;ltmp<numneigh[j];ltmp++) {
                  temp_jkp=BOP_index[j]+ltmp;
                  kp1=jlist[ltmp];
                  kp1type=map[type[kp1]]+1;
                  if(x[kp1][0]==x[k][0]) {
                    if(x[kp1][1]==x[k][1]) {
                      if(x[kp1][2]==x[k][2]) {
                        same_kpk=1;
                        break;
                      }
                    }
                  }
                }
                if(same_kpk){

//loop over neighbors of k

                  for(mtmp=0;mtmp<numneigh[k];mtmp++) {
                    kp2=klist[mtmp];
                    if(x[kp2][0]==x[k][0]) {
                      if(x[kp2][1]==x[k][1]) {
                        if(x[kp2][2]==x[k][2]) {
                          break;
                        }
                      }
                    }
                  }
                  if(jtype==ktype)
                    ijk=jtype-1;
                  else if(jtype < ktype)
                    ijk=jtype*bop_types-jtype*(jtype+1)/2+ktype-1;
                  else
                    ijk=ktype*bop_types-ktype*(ktype+1)/2+jtype-1;
                  if(jtype==kp1type)
                    ijkp=jtype-1;
                  else if(jtype<kp1type)
                    ijkp=jtype*bop_types-jtype*(jtype+1)/2+kp1type-1;
                  else
                    ijkp=kp1type*bop_types-kp1type*(kp1type+1)/2+jtype-1;

                  dis_jkp[0]=x[kp1][0]-x[j][0];
                  dis_jkp[1]=x[kp1][1]-x[j][1];
                  dis_jkp[2]=x[kp1][2]-x[j][2];
                  rsq_jkp=dis_jkp[0]*dis_jkp[0]
                      +dis_jkp[1]*dis_jkp[1]
                      +dis_jkp[2]*dis_jkp[2];
                  r_jkp=sqrt(rsq_jkp);
                  if(r_jkp<=rcut[ijkp]) {
                    ps=r_jkp*rdr[ijkp]+1.0;
                    ks=(int)ps;
                    if(nr-1<ks)
                      ks=nr-1;
                    ps=ps-ks;
                    if(ps>1.0)
                      ps=1.0;
                    betaS_jkp=((pBetaS3[ijkp][ks-1]*ps+pBetaS2[ijkp][ks-1])*ps
                        +pBetaS1[ijkp][ks-1])*ps+pBetaS[ijkp][ks-1];
                    dBetaS_jkp=(pBetaS6[ijkp][ks-1]*ps+pBetaS5[ijkp][ks-1])*ps
                        +pBetaS4[ijkp][ks-1];
                    cosAng_ijk=(-dis_ij[0]*dis_jk[0]-dis_ij[1]*dis_jk[1]
                        -dis_ij[2]*dis_jk[2])/(r_ij*r_jk);
                    dcA_ijk[0][0]=(dis_jk[0]*r_ij*r_jk-cosAng_ijk
                        *-dis_ij[0]*r_jk*r_jk)/(r_ij*r_ij*r_jk*r_jk);
                    dcA_ijk[1][0]=(dis_jk[1]*r_ij*r_jk-cosAng_ijk
                        *-dis_ij[1]*r_jk*r_jk)/(r_ij*r_ij*r_jk*r_jk);
                    dcA_ijk[2][0]=(dis_jk[2]*r_ij*r_jk-cosAng_ijk
                        *-dis_ij[2]*r_jk*r_jk)/(r_ij*r_ij*r_jk*r_jk);
                    dcA_ijk[0][1]=(-dis_ij[0]*r_ij*r_jk-cosAng_ijk
                        *dis_jk[0]*r_ij*r_ij)/(r_ij*r_ij*r_jk*r_jk);
                    dcA_ijk[1][1]=(-dis_ij[1]*r_ij*r_jk-cosAng_ijk
                        *dis_jk[1]*r_ij*r_ij)/(r_ij*r_ij*r_jk*r_jk);
                    dcA_ijk[2][1]=(-dis_ij[2]*r_ij*r_jk-cosAng_ijk
                        *dis_jk[2]*r_ij*r_ij)/(r_ij*r_ij*r_jk*r_jk);
                    gmean0=sigma_g0[itype-1][jtype-1][ktype-1];
                    gmean1=sigma_g1[itype-1][jtype-1][ktype-1];
                    gmean2=sigma_g2[itype-1][jtype-1][ktype-1];
                    amean=cosAng_ijk;
                    gfactor2=gmean0+gmean1*amean
                        +gmean2*amean*amean;
                    gprime2=gmean1+2.0*gmean2*amean;
                    gmean0=sigma_g0[itype-1][ktype-1][jtype-1];
                    gmean1=sigma_g1[itype-1][ktype-1][jtype-1];
                    gmean2=sigma_g2[itype-1][ktype-1][jtype-1];
                    cosAng_ikj=(dis_ik[0]*dis_jk[0]+dis_ik[1]*dis_jk[1]
                        +dis_ik[2]*dis_jk[2])/(r_ik*r_jk);
                    dcA_ikj[0][0]=(-dis_jk[0]*r_ik*r_jk-cosAng_ikj
                        *-dis_ik[0]*r_jk*r_jk)/(r_ik*r_ik*r_jk*r_jk);
                    dcA_ikj[1][0]=(-dis_jk[1]*r_ik*r_jk-cosAng_ikj
                        *-dis_ik[1]*r_jk*r_jk)/(r_ik*r_ik*r_jk*r_jk);
                    dcA_ikj[2][0]=(-dis_jk[2]*r_ik*r_jk-cosAng_ikj
                        *-dis_ik[2]*r_jk*r_jk)/(r_ik*r_ik*r_jk*r_jk);
                    dcA_ikj[0][1]=(-dis_ik[0]*r_ik*r_jk-cosAng_ikj
                        *-dis_jk[0]*r_ik*r_ik)/(r_ik*r_ik*r_jk*r_jk);
                    dcA_ikj[1][1]=(-dis_ik[1]*r_ik*r_jk-cosAng_ikj
                        *-dis_jk[1]*r_ik*r_ik)/(r_ik*r_ik*r_jk*r_jk);
                    dcA_ikj[2][1]=(-dis_ik[2]*r_ik*r_jk-cosAng_ikj
                        *-dis_jk[2]*r_ik*r_ik)/(r_ik*r_ik*r_jk*r_jk);
                    amean=cosAng_ikj;
                    gfactor3=gmean0+gmean1*amean
                        +gmean2*amean*amean;
                    gprime3=gmean1+2.0*gmean2*amean;
                    gfactor=gfactor1*gfactor2*gfactor3;
                    rfactor=betaS_ik*betaS_jkp;

//EE1 is (b) Eq. 12

                    EE1=EE1+gfactor*rfactor;

//rcm1 is derivative of EE1 w.r.t Beta(r_ik)
//rcm2 is derivative of EE1 w.r.t Beta(r_jk')
//gcm1 is derivative of EE1 w.r.t cos(theta_jik)
//gcm2 is derivative of EE1 w.r.t cos(theta_ijk)
//gcm3 is derivative of EE1 w.r.t cos(theta_ikj)

                    rcm1=gfactor*betaS_jkp*dBetaS_ik/r_ik;
                    rcm2=gfactor*betaS_ik*dBetaS_jkp/r_jkp;
                    gcm1=rfactor*gprime1*gfactor2*gfactor3;
                    gcm2=rfactor*gfactor1*gprime2*gfactor3;
                    gcm3=rfactor*gfactor1*gfactor2*gprime3;
                    bt_sg[nb_ij].dEE1[0]+=
                        gcm1*dcA_jik[0][0]
                        -gcm2*dcA_ijk[0][0];
                    bt_sg[nb_ij].dEE1[1]+=
                        gcm1*dcA_jik[1][0]
                        -gcm2*dcA_ijk[1][0];
                    bt_sg[nb_ij].dEE1[2]+=
                        gcm1*dcA_jik[2][0]
                        -gcm2*dcA_ijk[2][0];
                    bt_sg[nb_ik].dEE1[0]+=
                        gcm1*dcA_jik[0][1]
                        +rcm1*dis_ik[0]
                        -gcm3*dcA_ikj[0][0];
                    bt_sg[nb_ik].dEE1[1]+=
                        gcm1*dcA_jik[1][1]
                        +rcm1*dis_ik[1]
                        -gcm3*dcA_ikj[1][0];
                    bt_sg[nb_ik].dEE1[2]+=
                        gcm1*dcA_jik[2][1]
                        +rcm1*dis_ik[2]
                        -gcm3*dcA_ikj[2][0];
                    bt_sg[nb_jk].dEE1[0]+=
                        gcm2*dcA_ijk[0][1]
                        +rcm2*dis_jkp[0]
                        -gcm3*dcA_ikj[0][1];
                    bt_sg[nb_jk].dEE1[1]+=
                        gcm2*dcA_ijk[1][1]
                        +rcm2*dis_jkp[1]
                        -gcm3*dcA_ikj[1][1];
                    bt_sg[nb_jk].dEE1[2]+=
                        gcm2*dcA_ijk[2][1]
                        +rcm2*dis_jkp[2]
                        -gcm3*dcA_ikj[2][1];
                  }
                }

// k and k' and j are all different neighbors of i

                for(ltmp=0;ltmp<ktmp;ltmp++) {
                  if(ltmp!=jtmp) {
                    kp=iilist[ltmp];;
                    kptype = map[type[kp]]+1;
                    if(itype==kptype)
                      iikp=itype-1;
                    else if(itype<kptype)
                      iikp=itype*bop_types-itype*(itype+1)/2+kptype-1;
                    else
                      iikp=kptype*bop_types-kptype*(kptype+1)/2+itype-1;
                    for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                      ncmp=itypeSigBk[n][nsearch];
                      if(x[ncmp][0]==x[kp][0]) {
                        if(x[ncmp][1]==x[kp][1]) {
                          if(x[ncmp][2]==x[kp][2]) {
                            break;
                          }
                        }
                      }
                    }
                    dis_ikp[0]=x[kp][0]-x[i][0];
                    dis_ikp[1]=x[kp][1]-x[i][1];
                    dis_ikp[2]=x[kp][2]-x[i][2];
                    rsq_ikp=dis_ikp[0]*dis_ikp[0]
                        +dis_ikp[1]*dis_ikp[1]
                        +dis_ikp[2]*dis_ikp[2];
                    r_ikp=sqrt(rsq_ikp);
                    if(r_ikp<=rcut[iikp]) {
                      ps=r_ikp*rdr[iikp]+1.0;
                      ks=(int)ps;
                      if(nr-1<ks)
                        ks=nr-1;
                      ps=ps-ks;
                      if(ps>1.0)
                        ps=1.0;
                      betaS_ikp=((pBetaS3[iikp][ks-1]*ps+pBetaS2[iikp][ks-1])*ps
                          +pBetaS1[iikp][ks-1])*ps+pBetaS[iikp][ks-1];
                      gmean0=sigma_g0[jtype-1][itype-1][kptype-1];
                      gmean1=sigma_g1[jtype-1][itype-1][kptype-1];
                      gmean2=sigma_g2[jtype-1][itype-1][kptype-1];
                      cosAng_jikp=(dis_ij[0]*dis_ikp[0]+dis_ij[1]*dis_ikp[1]
                          +dis_ij[2]*dis_ikp[2])/(r_ij*r_ikp);
                      cosAng_kikp=(dis_ik[0]*dis_ikp[0]+dis_ik[1]*dis_ikp[1]
                          +dis_ik[2]*dis_ikp[2])/(r_ik*r_ikp);
                      amean=cosAng_jikp;
                      gfactor2=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime2=gmean1+2.0*gmean2*amean;
                      gmean0=sigma_g0[ktype-1][itype-1][kptype-1];
                      gmean1=sigma_g1[ktype-1][itype-1][kptype-1];
                      gmean2=sigma_g2[ktype-1][itype-1][kptype-1];
                      amean=cosAng_kikp;
                      gfactor3=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime3=gmean1+2.0*gmean2*amean;
                      gfactor=gfactor1*gfactor2*gfactor3;
                      rfactorrt=betaS_ik*betaS_ikp;
                      rfactor=rfactorrt*rfactorrt;

//2nd CC is second term of Eq. 11 (c) for i atom where j , k & k' =neighbor of i

                      CC=CC+2.0*gfactor*rfactor;
                    }
                  }
                }

// j and k are different neighbors of i and k' is a neighbor k not equal to i

                for(ltmp=0;ltmp<numneigh[k];ltmp++) {
                  temp_kkp=BOP_index[k]+ltmp;
                  kp=klist[ltmp];;
                  kptype = map[type[kp]]+1;
                  same_ikp=0;
                  same_jkp=0;
                  if(x[i][0]==x[kp][0]) {
                    if(x[i][1]==x[kp][1]) {
                      if(x[i][2]==x[kp][2]) {
                        same_ikp=1;
                      }
                    }
                  }
                  if(x[j][0]==x[kp][0]) {
                    if(x[j][1]==x[kp][1]) {
                      if(x[j][2]==x[kp][2]) {
                        same_jkp=1;
                      }
                    }
                  }
                  if(!same_ikp&&!same_jkp) {
                    if(ktype==kptype)
                      ikkp=ktype-1;
                    else if(ktype<kptype)
                      ikkp=ktype*bop_types-ktype*(ktype+1)/2+kptype-1;
                    else
                      ikkp=kptype*bop_types-kptype*(kptype+1)/2+ktype-1;
                    dis_kkp[0]=x[kp][0]-x[k][0];
                    dis_kkp[1]=x[kp][1]-x[k][1];
                    dis_kkp[2]=x[kp][2]-x[k][2];
                    rsq_kkp=dis_kkp[0]*dis_kkp[0]
                        +dis_kkp[1]*dis_kkp[1]
                        +dis_kkp[2]*dis_kkp[2];
                    r_kkp=sqrt(rsq_kkp);
                    if(r_kkp<=rcut[ikkp]) {
                      ps=r_kkp*rdr[ikkp]+1.0;
                      ks=(int)ps;
                      if(nr-1<ks)
                        ks=nr-1;
                      ps=ps-ks;
                      if(ps>1.0)
                        ps=1.0;
                      betaS_kkp=((pBetaS3[ikkp][ks-1]*ps+pBetaS2[ikkp][ks-1])*ps
                          +pBetaS1[ikkp][ks-1])*ps+pBetaS[ikkp][ks-1];
                      sig_flag=0;
                      for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                        ncmp=itypeSigBk[n][nsearch];
                        if(x[ncmp][0]==x[kp][0]) {
                          if(x[ncmp][1]==x[kp][1]) {
                            if(x[ncmp][2]==x[kp][2]) {
                              sig_flag=1;
                              nkp=nsearch;
                              break;
                            }
                          }
                        }
                      }
                      if(sig_flag==0) {
                        nSigBk[n]=nSigBk[n]+1;
                        nkp=nSigBk[n]-1;
                        itypeSigBk[n][nkp]=kp;
                      }
                      cosAng_ikkp=(-dis_ik[0]*dis_kkp[0]-dis_ik[1]*dis_kkp[1]
                          -dis_ik[2]*dis_kkp[2])/(r_ik*r_kkp);
                      gmean0=sigma_g0[itype-1][ktype-1][kptype-1];
                      gmean1=sigma_g1[itype-1][ktype-1][kptype-1];
                      gmean2=sigma_g2[itype-1][ktype-1][kptype-1];
                      amean=cosAng_ikkp;
                      gfactor2=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime2=gmean1+2.0*gmean2*amean;
                      gfactorsq2=gfactor2*gfactor2;
                      gfactor=gfactorsq*gfactorsq2;
                      rfactorrt=betaS_ik*betaS_kkp;
                      rfactor=rfactorrt*rfactorrt;

//3rd CC is third term of Eq. 11 (c) for i atom
//where j , k =neighbor of i & k' =neighbor of k

                      CC=CC+gfactor*rfactor;
                    }
                  }
                }
              }
            }
          }

//j is a neighbor of i and k is a neighbor of j not equal to i

          for(ktmp=0;ktmp<numneigh[j];ktmp++) {
            if(ktmp!=ji) {
              temp_jk=BOP_index[j]+ktmp;
              k=jlist[ktmp];
              klist=firstneigh[k];
              ktype=map[type[k]]+1;
              for(kNeij=0;kNeij<numneigh[k];kNeij++) {
                if(x[klist[kNeij]][0]==x[j][0]) {
                  if(x[klist[kNeij]][1]==x[j][1]) {
                    if(x[klist[kNeij]][2]==x[j][2]) {
                      break;
                    }
                  }
                }
              }
              if(jtype==ktype)
                ijk=jtype-1;
              else if(jtype<ktype)
                ijk=jtype*bop_types-jtype*(jtype+1)/2+ktype-1;
              else
                ijk=ktype*bop_types-ktype*(ktype+1)/2+jtype-1;
              sig_flag=0;
              for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                ncmp=itypeSigBk[n][nsearch];
                if(x[ncmp][0]==x[k][0]) {
                  if(x[ncmp][1]==x[k][1]) {
                    if(x[ncmp][2]==x[k][2]) {
                      new1=nsearch;
                      sig_flag=1;
                      break;
                    }
                  }
                }
              }
              if(sig_flag==0) {
                nSigBk[n]=nSigBk[n]+1;
                new1=nSigBk[n]-1;
                itypeSigBk[n][new1]=k;
              }
              dis_jk[0]=x[k][0]-x[j][0];
              dis_jk[1]=x[k][1]-x[j][1];
              dis_jk[2]=x[k][2]-x[j][2];
              rsq_jk=dis_jk[0]*dis_jk[0]
                  +dis_jk[1]*dis_jk[1]
                  +dis_jk[2]*dis_jk[2];
              r_jk=sqrt(rsq_jk);
              if(r_jk<=rcut[ijk]) {
                ps=r_jk*rdr[ijk]+1.0;
                ks=(int)ps;
                if(nr-1<ks)
                  ks=nr-1;
                ps=ps-ks;
                if(ps>1.0)
                  ps=1.0;
                betaS_jk=((pBetaS3[ijk][ks-1]*ps+pBetaS2[ijk][ks-1])*ps
                    +pBetaS1[ijk][ks-1])*ps+pBetaS[ijk][ks-1];
                dBetaS_jk=(pBetaS6[ijk][ks-1]*ps+pBetaS5[ijk][ks-1])*ps
                    +pBetaS4[ijk][ks-1];
                cosAng_ijk=(-dis_ij[0]*dis_jk[0]-dis_ij[1]*dis_jk[1]
                    -dis_ij[2]*dis_jk[2])/(r_ij*r_jk);
                dcA_ijk[0][0]=(dis_jk[0]*r_ij*r_jk-cosAng_ijk
                    *-dis_ij[0]*r_jk*r_jk)/(r_ij*r_ij*r_jk*r_jk);
                dcA_ijk[1][0]=(dis_jk[1]*r_ij*r_jk-cosAng_ijk
                    *-dis_ij[1]*r_jk*r_jk)/(r_ij*r_ij*r_jk*r_jk);
                dcA_ijk[2][0]=(dis_jk[2]*r_ij*r_jk-cosAng_ijk
                    *-dis_ij[2]*r_jk*r_jk)/(r_ij*r_ij*r_jk*r_jk);
                dcA_ijk[0][1]=(-dis_ij[0]*r_ij*r_jk-cosAng_ijk
                    *dis_jk[0]*r_ij*r_ij)/(r_ij*r_ij*r_jk*r_jk);
                dcA_ijk[1][1]=(-dis_ij[1]*r_ij*r_jk-cosAng_ijk
                    *dis_jk[1]*r_ij*r_ij)/(r_ij*r_ij*r_jk*r_jk);
                dcA_ijk[2][1]=(-dis_ij[2]*r_ij*r_jk-cosAng_ijk
                    *dis_jk[2]*r_ij*r_ij)/(r_ij*r_ij*r_jk*r_jk);
                nb_jk=nb_t;
                nb_t++;
                if(nb_t>nb_sg) {
                  new_n_tot=nb_sg+maxneigh;
                  grow_sigma(nb_sg,new_n_tot);
                  nb_sg=new_n_tot;
                }
                bt_sg[nb_jk].temp=temp_jk;
                bt_sg[nb_jk].i=j;
                bt_sg[nb_jk].j=k;
                gmean0=sigma_g0[itype-1][jtype-1][ktype-1];
                gmean1=sigma_g1[itype-1][jtype-1][ktype-1];
                gmean2=sigma_g2[itype-1][jtype-1][ktype-1];
                amean=cosAng_ijk;
                gfactor1=gmean0+gmean1*amean
                    +gmean2*amean*amean;
                gprime1=gmean1+2.0*gmean2*amean;
                gfactorsq=gfactor1*gfactor1;
                gsqprime=2.0*gfactor1*gprime1;
                rfactor1rt=betaS_jk*betaS_jk;
                rfactor1=rfactor1rt*rfactor1rt;

//BB is Eq. 34 (a) or Eq. 10 (c) for the j atom
//1st DD is Eq. 11 (c) for j atom where i & k=neighbor of j

                BB=BB+gfactorsq*rfactor1rt;
                DD=DD+gfactorsq*rfactor1;

//agpdpr1 is derivative of BB  w.r.t. Beta(r_jk)
//app1 is derivative of BB w.r.t. cos(theta_ijk)

                agpdpr1=2.0*gfactorsq*betaS_jk*dBetaS_jk/r_jk;
                app1=rfactor1rt*gsqprime;
                bt_sg[nb_ij].dBB[0]-=
                    app1*dcA_ijk[0][0];
                bt_sg[nb_ij].dBB[1]-=
                    app1*dcA_ijk[1][0];
                bt_sg[nb_ij].dBB[2]-=
                    app1*dcA_ijk[2][0];
                bt_sg[nb_jk].dBB[0]+=
                    app1*dcA_ijk[0][1]
                    +agpdpr1*dis_jk[0];
                bt_sg[nb_jk].dBB[1]+=
                    app1*dcA_ijk[1][1]
                    +agpdpr1*dis_jk[1];
                bt_sg[nb_jk].dBB[2]+=
                    app1*dcA_ijk[2][1]
                    +agpdpr1*dis_jk[2];

//j is a neighbor of i, k and k' prime different neighbors of j not equal to i

                for(ltmp=0;ltmp<ktmp;ltmp++) {
                  if(ltmp!=ji) {
                    temp_jkp=BOP_index[j]+ltmp;
                    kp=jlist[ltmp];
                    kptype=map[type[kp]]+1;
                    if(jtype==kptype)
                      ijkp=jtype-1;
                    else if(jtype<kptype)
                      ijkp=jtype*bop_types-jtype*(jtype+1)/2+kptype-1;
                    else
                      ijkp=kptype*bop_types-kptype*(kptype+1)/2+jtype-1;
                    for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                      ncmp=itypeSigBk[n][nsearch];
                      if(x[ncmp][0]==x[kp][0]) {
                        if(x[ncmp][1]==x[kp][1]) {
                          if(x[ncmp][2]==x[kp][2]) {
                            new2=nsearch;
                            break;
                          }
                        }
                      }
                    }
                    dis_jkp[0]=x[kp][0]-x[j][0];
                    dis_jkp[1]=x[kp][1]-x[j][1];
                    dis_jkp[2]=x[kp][2]-x[j][2];
                    rsq_jkp=dis_jkp[0]*dis_jkp[0]
                        +dis_jkp[1]*dis_jkp[1]
                        +dis_jkp[2]*dis_jkp[2];
                    r_jkp=sqrt(rsq_jkp);
                    if(r_jkp<=rcut[ijkp]) {
                      ps=r_jkp*rdr[ijkp]+1.0;
                      ks=(int)ps;
                      if(nr-1<ks)
                        ks=nr-1;
                      ps=ps-ks;
                      if(ps>1.0)
                        ps=1.0;
                      betaS_jkp=((pBetaS3[ijkp][ks-1]*ps+pBetaS2[ijkp][ks-1])*ps
                        +pBetaS1[ijkp][ks-1])*ps+pBetaS[ijkp][ks-1];
                      dBetaS_jkp=(pBetaS6[ijkp][ks-1]*ps+pBetaS5[ijkp][ks-1])*ps
                        +pBetaS4[ijkp][ks-1];
                      cosAng_ijkp=(-dis_ij[0]*dis_jkp[0]-dis_ij[1]*dis_jkp[1]
                        -dis_ij[2]*dis_jkp[2])/(r_ij*r_jkp);
                      cosAng_kjkp=(dis_jk[0]*dis_jkp[0]+dis_jk[1]*dis_jkp[1]
                        +dis_jk[2]*dis_jkp[2])/(r_jk*r_jkp);
                      nb_jkp=nb_t;
                      nb_t++;
                      if(nb_t>nb_sg) {
                        new_n_tot=nb_sg+maxneigh;
                        grow_sigma(nb_sg,new_n_tot);
                        nb_sg=new_n_tot;
                      }
                      bt_sg[nb_jkp].temp=temp_jkp;
                      bt_sg[nb_jkp].i=j;
                      bt_sg[nb_jkp].j=kp;
                      gmean0=sigma_g0[itype-1][jtype-1][kptype-1];
                      gmean1=sigma_g1[itype-1][jtype-1][kptype-1];
                      gmean2=sigma_g2[itype-1][jtype-1][kptype-1];
                      amean=cosAng_ijkp;
                      gfactor2=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime2=gmean1+2.0*gmean2*amean;
                      gmean0=sigma_g0[ktype-1][jtype-1][kptype-1];
                      gmean1=sigma_g1[ktype-1][jtype-1][kptype-1];
                      gmean2=sigma_g2[ktype-1][jtype-1][kptype-1];
                      amean=cosAng_kjkp;
                      gfactor3=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime3=gmean1+2.0*gmean2*amean;
                      gfactor=gfactor1*gfactor2*gfactor3;
                      rfactorrt=betaS_jk*betaS_jkp;
                      rfactor=rfactorrt*rfactorrt;

//2nd DD is Eq. 11 (c) for j atom where i , k & k'=neighbor of j

                      DD=DD+2.0*gfactor*rfactor;
                    }
                  }
                }

//j is a neighbor of i, k is a neighbor of j not equal to i and k'
//is a neighbor of k not equal to j or i

                for(ltmp=0;ltmp<numneigh[k];ltmp++) {
                  temp_kkp=BOP_index[k]+ltmp;
                  kp=klist[ltmp];
                  kptype=map[type[kp]]+1;
                  same_ikp=0;
                  same_jkp=0;
                  if(x[i][0]==x[kp][0]) {
                    if(x[i][1]==x[kp][1]) {
                      if(x[i][2]==x[kp][2]) {
                        same_ikp=1;
                      }
                    }
                  }
                  if(x[j][0]==x[kp][0]) {
                    if(x[j][1]==x[kp][1]) {
                      if(x[j][2]==x[kp][2]) {
                        same_jkp=1;
                      }
                    }
                  }
                  if(!same_ikp&&!same_jkp) {
                    if(ktype==kptype)
                      ikkp=ktype-1;
                    else if(ktype<kptype)
                      ikkp=ktype*bop_types-ktype*(ktype+1)/2+kptype-1;
                    else
                      ikkp=kptype*bop_types-kptype*(kptype+1)/2+ktype-1;
                    for(kNeij=0;kNeij<numneigh[k];kNeij++) {
                      if(x[klist[kNeij]][0]==x[j][0]) {
                        if(x[klist[kNeij]][1]==x[j][1]) {
                          if(x[klist[kNeij]][2]==x[j][2]) {
                            break;
                          }
                        }
                      }
                    }
                    sig_flag=0;
                    for(nsearch=0;nsearch<nSigBk[n];nsearch++) {
                      ncmp=itypeSigBk[n][nsearch];
                      if(x[ncmp][0]==x[kp][0]) {
                        if(x[ncmp][1]==x[kp][1]) {
                          if(x[ncmp][2]==x[kp][2]) {
                            new2=nsearch;
                            sig_flag=1;
                            break;
                          }
                        }
                      }
                    }
                    if(sig_flag==0) {
                      nSigBk[n]=nSigBk[n]+1;
                      new2=nSigBk[n]-1;
                      itypeSigBk[n][new2]=kp;
                    }
                    dis_kkp[0]=x[kp][0]-x[k][0];
                    dis_kkp[1]=x[kp][1]-x[k][1];
                    dis_kkp[2]=x[kp][2]-x[k][2];
                    rsq_kkp=dis_kkp[0]*dis_kkp[0]
                        +dis_kkp[1]*dis_kkp[1]
                        +dis_kkp[2]*dis_kkp[2];
                    r_kkp=sqrt(rsq_kkp);
                    if(r_kkp<=rcut[ikkp]) {
                      ps=r_kkp*rdr[ikkp]+1.0;
                      ks=(int)ps;
                      if(nr-1<ks)
                        ks=nr-1;
                      ps=ps-ks;
                      if(ps>1.0)
                        ps=1.0;
                      betaS_kkp=((pBetaS3[ikkp][ks-1]*ps+pBetaS2[ikkp][ks-1])*ps
                          +pBetaS1[ikkp][ks-1])*ps+pBetaS[ikkp][ks-1];
                      cosAng_jkkp=(-dis_jk[0]*dis_kkp[0]-dis_jk[1]*dis_kkp[1]
                          -dis_jk[2]*dis_kkp[2])/(r_jk*r_kkp);
                      nb_kkp=nb_t;
                      nb_t++;
                      if(nb_t>nb_sg) {
                        new_n_tot=nb_sg+maxneigh;
                        grow_sigma(nb_sg,new_n_tot);
                        nb_sg=new_n_tot;
                      }
                      bt_sg[nb_kkp].temp=temp_kkp;
                      bt_sg[nb_kkp].i=k;
                      bt_sg[nb_kkp].j=kp;
                      gmean0=sigma_g0[jtype-1][ktype-1][kptype-1];
                      gmean1=sigma_g1[jtype-1][ktype-1][kptype-1];
                      gmean2=sigma_g2[jtype-1][ktype-1][kptype-1];
                      amean=cosAng_jkkp;
                      gfactor2=gmean0+gmean1*amean
                          +gmean2*amean*amean;
                      gprime2=gmean1+2.0*gmean2*amean;
                      gfactorsq2=gfactor2*gfactor2;
                      gfactor=gfactorsq*gfactorsq2;
                      rfactorrt=betaS_jk*betaS_kkp;
                      rfactor=rfactorrt*rfactorrt;

//3rd DD is Eq. 11 (c) for j atom where i & k=neighbor of j & k'=neighbor of k

                      DD=DD+gfactor*rfactor;
                    }
                  }
                }
              }
            }
          }

          sig_flag=0;
          if(FF<=0.000001) {
            sigB[n]=0.0;
            sig_flag=1;
          }
          if(sig_flag==0) {
            if(AA<0.0)
              AA=0.0;
            if(BB<0.0)
              BB=0.0;
            if(CC<0.0)
              CC=0.0;
            if(DD<0.0)
              DD=0.0;

// AA and BB are the representations of (a) Eq. 34 and (b) Eq. 9
// for atoms i and j respectively

            AAC=AA+BB;
            BBC=AA*BB;
            CCC=AA*AA+BB*BB;
            DDC=CC+DD;

//EEC is a modified form of (a) Eq. 33

            EEC=(DDC-CCC)/(AAC+2.0*small1);
            for(m=0;m<nb_t;m++) {
              if((bt_sg[m].i>-1)&&(bt_sg[m].j>-1)) {
                bt_sg[m].dAAC[0]=bt_sg[m].dAA[0]
                    +bt_sg[m].dBB[0];
                bt_sg[m].dAAC[1]=bt_sg[m].dAA[1]
                    +bt_sg[m].dBB[1];
                bt_sg[m].dAAC[2]=bt_sg[m].dAA[2]
                    +bt_sg[m].dBB[2];
              }
            }
            UT=EEC*FF+BBC+small3[iij];
            UT=1.0/sqrt(UT);

// bndtmp is a slightly modified form of (a) Eq. 30 and (b) Eq. 8

            bndtmp=(FF+sigma_delta[iij]*sigma_delta[iij])
                +sigma_c[iij]*AAC+small4;
            psign=1.0;
            bndtmp0=1.0/sqrt(bndtmp);
            sigB1[n]=psign*betaS_ij*bndtmp0;
            bndtmp=-0.5*bndtmp0*bndtmp0*bndtmp0;
            bndtmp1=psign*bndtmp0+psign*betaS_ij
                *bndtmp*2.0*betaS_ij;
            bndtmp1=bndtmp1*dBetaS_ij/r_ij;
            bndtmp2=psign*betaS_ij*bndtmp*sigma_c[iij];
            setting=0;
            for(m=0;m<nb_t;m++) {
              if((bt_sg[m].i>-1)&&(bt_sg[m].j>-1)) {
                temp_kk=bt_sg[m].temp;
                if(temp_kk==temp_ij&&setting==0) {
                  bt_sg[m].dSigB1[0]=bndtmp1*dis_ij[0]
                      +(bndtmp2*bt_sg[m].dAAC[0]);
                  bt_sg[m].dSigB1[1]=bndtmp1*dis_ij[1]
                      +(bndtmp2*bt_sg[m].dAAC[1]);
                  bt_sg[m].dSigB1[2]=bndtmp1*dis_ij[2]
                      +(bndtmp2*bt_sg[m].dAAC[2]);
                  setting=1;
                }
                else if(temp_kk==temp_ji&&setting==0) {
                  bt_sg[m].dSigB1[0]=-bndtmp1*dis_ij[0]
                      +(bndtmp2*bt_sg[m].dAAC[0]);
                  bt_sg[m].dSigB1[1]=-bndtmp1*dis_ij[1]
                      +(bndtmp2*bt_sg[m].dAAC[1]);
                  bt_sg[m].dSigB1[2]=-bndtmp1*dis_ij[2]
                      +(bndtmp2*bt_sg[m].dAAC[2]);
                  setting=1;
                }
                else {
                  bt_sg[m].dSigB1[0]=(bndtmp2*bt_sg[m].dAAC[0]);
                  bt_sg[m].dSigB1[1]=(bndtmp2*bt_sg[m].dAAC[1]);
                  bt_sg[m].dSigB1[2]=(bndtmp2*bt_sg[m].dAAC[2]);
                }
              }
            }

//This loop is to ensure there is not an error for atoms with no neighbors (deposition)

            if(nb_t==0) {
              if(j>i) {
                bt_sg[0].dSigB1[0]=bndtmp1*dis_ij[0];
                bt_sg[0].dSigB1[1]=bndtmp1*dis_ij[1];
                bt_sg[0].dSigB1[2]=bndtmp1*dis_ij[2];
              }
              else {
                bt_sg[0].dSigB1[0]=-bndtmp1*dis_ij[0];
                bt_sg[0].dSigB1[1]=-bndtmp1*dis_ij[1];
                bt_sg[0].dSigB1[2]=-bndtmp1*dis_ij[2];
              }
              for(pp=0;pp<3;pp++) {
                bt_sg[0].dAA[pp]=0.0;
                bt_sg[0].dBB[pp]=0.0;
                bt_sg[0].dEE1[pp]=0.0;
                bt_sg[0].dFF[pp]=0.0;
                bt_sg[0].dAAC[pp]=0.0;
                bt_sg[0].dSigB[pp]=0.0;
              }
              bt_sg[0].i=i;
              bt_sg[0].j=j;
              bt_sg[0].temp=temp_ij;
              nb_t++;
              if(nb_t>nb_sg) {
                new_n_tot=nb_sg+maxneigh;
                grow_sigma(nb_sg,new_n_tot);
                nb_sg=new_n_tot;
              }
            }
            ps=sigB1[n]*rdBO+1.0;
            ks=(int)ps;
            if(nBOt-1<ks)
              ks=nBOt-1;
            ps=ps-ks;
            if(ps>1.0)
              ps=1.0;
            dsigB1=((FsigBO3[iij][ks-1]*ps+FsigBO2[iij][ks-1])*ps
                +FsigBO1[iij][ks-1])*ps+FsigBO[iij][ks-1];
            dsigB2=(FsigBO6[iij][ks-1]*ps+FsigBO5[iij][ks-1])*ps+FsigBO4[iij][ks-1];
            part0=(FF+0.5*AAC+small5);
            part1=(sigma_f[iij]-0.5)*sigma_k[iij];
            part2=1.0-part1*EE1/part0;
            part3=dsigB1*part1/part0;
            part4=part3/part0*EE1;

// sigB is the final expression for (a) Eq. 6 and (b) Eq. 11

            sigB[n]=dsigB1*part2;
            pp1=2.0*betaS_ij;
            for(m=0;m<nb_t;m++) {
              if((bt_sg[m].i>-1)&&(bt_sg[m].j>-1)) {
                temp_kk=bt_sg[m].temp;
                bt_i=bt_sg[m].i;
                bt_j=bt_sg[m].j;
                xtmp[0]=x[bt_j][0]-x[bt_i][0];
                xtmp[1]=x[bt_j][1]-x[bt_i][1];
                xtmp[2]=x[bt_j][2]-x[bt_i][2];
                for(pp=0;pp<3;pp++) {
                  bt_sg[m].dSigB[pp]=dsigB2*part2*bt_sg[m].dSigB1[pp]
                      -part3*bt_sg[m].dEE1[pp]
                      +part4*(bt_sg[m].dFF[pp]
                      +0.5*bt_sg[m].dAAC[pp]);
                }
                for(pp=0;pp<3;pp++) {
                  ftmp[pp]=pp1*bt_sg[m].dSigB[pp];
                  f[bt_i][pp]-=ftmp[pp];
                  f[bt_j][pp]+=ftmp[pp];
                }
                if(evflag) {
                  ev_tally_xyz(bt_i,bt_j,nlocal,newton_pair,0.0,0.0,ftmp[0],ftmp[1]
                      ,ftmp[2],xtmp[0],xtmp[1],xtmp[2]);
                }
              }
            }
          }
          n++;
        }
      }
    }
  }
  destroy_sigma();
}

/* ---------------------------------------------------------------------- */

void PairBOP::PiBo()
{
  int new_n_tot;
  int i,j,k,kp,m,n,pp,nb_t;
  int iij,ji,ki;
  int nsearch,ncmp;
  tagint i_tag,j_tag;
  int njik,ngj,ngk,nglj,ngl,ngi;
  int nkjkp,nijkp,ngli,nkikp,njikp;
  int itmp,ltmp,jtmp,ktmp;
  int nlocal,pi_flag;
  int inum,*ilist,*iilist,*jlist;
  int **firstneigh,*numneigh;
  int itype,jtype;
  int temp_ij,temp_ik,temp_ikp;
  int temp_jk,temp_jkp;
  int ang_jikp,ang_kikp,ang_ijk;
  int ang_ijkp,ang_kjkp,ang_jik;
  int nb_ij,nb_ik,nb_jk,nb_ikp,nb_jkp;
  int bt_ij,bt_i,bt_j;
  double AA,BB,CC;
  double cosSq,sinFactor,cosFactor;
  double cosSq1,dotV,BBrt,AB1,AB2;
  double BBrtR,ABrtR1,ABrtR2;
  double angFactor,angFactor1,angFactor2;
  double angFactor3,angFactor4,angRfactor;
  double dAngR1,dAngR2,agpdpr3;
  double agpdpr1,agpdpr2,app1,app2,app3;
  double betaCapSq1,dbetaCapSq1;
  double betaCapSq2,dbetaCapSq2;
  double betaCapSum,ftmp[3];
  double dPiB1,dPiB2,dPiB3,pp2;
  double **f = atom->f;
  double **x = atom->x;
  int *type = atom->type;
  tagint *tag = atom->tag;
  int newton_pair = force->newton_pair;

  nlocal = atom->nlocal;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  n=0;

// Loop over all local atoms for i

  if(nb_pi>16) {
    nb_pi=16;
  }
  if(nb_pi==0) {
    nb_pi=(maxneigh)*(maxneigh/2);
  }
  if(allocate_pi) {
    destroy_pi();
  }
  create_pi(nb_pi);
  for(itmp=0;itmp<inum;itmp++) {
    nb_t=0;
    i = ilist[itmp];
    itype = map[type[i]]+1;
    i_tag=tag[i];

// j is a loop over all neighbors of i

    iilist=firstneigh[i];
    for(jtmp=0;jtmp<numneigh[i];jtmp++) {
      temp_ij=BOP_index[i]+jtmp;
      if(neigh_flag[temp_ij]) {
        for(m=0;m<nb_pi;m++) {
          for(pp=0;pp<3;pp++) {
            bt_pi[m].dAA[pp]=0.0;
            bt_pi[m].dBB[pp]=0.0;
            bt_pi[m].dPiB[pp]=0.0;
          }
          bt_pi[m].i=-1;
          bt_pi[m].j=-1;
        }
        j=iilist[jtmp];
        jlist=firstneigh[j];
        jtype=map[type[j]]+1;
        j_tag=tag[j];
        nb_t=0;
        ftmp[0]=0.0;
        ftmp[1]=0.0;
        ftmp[2]=0.0;
        nb_ij=nb_t;
        nb_t++;
        if(nb_t>nb_pi) {
          new_n_tot=nb_pi+maxneigh;
          grow_pi(nb_pi,new_n_tot);
          nb_pi=new_n_tot;
        }
        bt_pi[nb_ij].i=i;
        bt_pi[nb_ij].j=j;
        bt_pi[nb_ij].temp=temp_ij;
        if(j_tag>=i_tag) {
          if(itype==jtype)
            iij=itype-1;
          else if(itype<jtype)
            iij=itype*bop_types-itype*(itype+1)/2+jtype-1;
          else
            iij=jtype*bop_types-jtype*(jtype+1)/2+itype-1;
          AA=0.0;
          BB=0.0;
          nPiBk[n]=0;
          for(ji=0;ji<numneigh[j];ji++) {
            if(x[jlist[ji]][0]==x[i][0]) {
              if(x[jlist[ji]][1]==x[i][1]) {
                if(x[jlist[ji]][2]==x[i][2]) {
                  break;
                }
              }
            }
          }

// j and k are different neighbors of i

          for(ktmp=0;ktmp<numneigh[i];ktmp++) {
            if(ktmp!=jtmp) {
              temp_ik=BOP_index[i]+ktmp;
              if(neigh_flag[temp_ik]) {
                k=iilist[ktmp];
                if(jtmp<ktmp) {
                  njik=jtmp*(2*numneigh[i]-jtmp-1)/2+(ktmp-jtmp)-1;
                  ngj=0;
                  ngk=1;
                }
                else {
                  njik=ktmp*(2*numneigh[i]-ktmp-1)/2+(jtmp-ktmp)-1;
                  ngj=1;
                  ngk=0;
                }
                ang_jik=cos_index[i]+njik;
                if(ang_jik>=cos_total) {
                  error->one(FLERR,"Too many atom triplets for pair bop");
                }
                nb_ik=nb_t;
                nb_t++;
                if(nb_t>nb_pi) {
                  new_n_tot=nb_pi+maxneigh;
                  grow_pi(nb_pi,new_n_tot);
                  nb_pi=new_n_tot;
                }
                bt_pi[nb_ik].i=i;
                bt_pi[nb_ik].j=k;
                bt_pi[nb_ik].temp=temp_ik;
                cosSq=cosAng[ang_jik]*cosAng[ang_jik];
                sinFactor=.5*(1.0-cosSq)*pi_p[itype-1]*betaS[temp_ik];
                cosFactor=.5*(1.0+cosSq)*betaP[temp_ik];
                betaCapSq1=pi_p[itype-1]*betaS[temp_ik]*betaS[temp_ik]-betaP[temp_ik]
                    *betaP[temp_ik];
                dbetaCapSq1=2.0*pi_p[itype-1]*betaS[temp_ik]*dBetaS[temp_ik]
                    -2.0*betaP[temp_ik]*dBetaP[temp_ik];

//AA is Eq. 37 (a) and Eq. 19 (b) or i atoms
//1st BB is first term of Eq. 38 (a) where j and k =neighbors i

                AA=AA+sinFactor*betaS[temp_ik]+cosFactor*betaP[temp_ik];
                BB=BB+.25*(1.0-cosSq)*(1.0-cosSq)*betaCapSq1*betaCapSq1;

//agpdpr1 is derivative of AA w.r.t. for atom i w.r.t. Beta(r_ik)
//agpdpr2 is derivative of BB w.r.t. for atom i w.r.t. Beta(r_ik)
//app1 is derivative of AA w.r.t. for atom i w.r.t. cos(theta_jik)
//app2 is derivative of BB w.r.t. for atom i w.r.t. cos(theta_jik)

                agpdpr1=(2.0*sinFactor*dBetaS[temp_ik]+2.0*cosFactor
                    *dBetaP[temp_ik])/rij[temp_ik];
                app1=cosAng[ang_jik]*(-pi_p[itype-1]*betaS[temp_ik]*betaS[temp_ik]
                    +betaP[temp_ik]*betaP[temp_ik]);
                app2=-(1.0-cosSq)*cosAng[ang_jik]*betaCapSq1*betaCapSq1;
                agpdpr2=.5*(1.0-cosSq)*(1.0-cosSq)*betaCapSq1*dbetaCapSq1/rij[temp_ik];
                itypePiBk[n][nPiBk[n]]=k;
                bt_pi[nb_ij].dAA[0]+=
                    app1*dcAng[ang_jik][0][ngj];
                bt_pi[nb_ij].dAA[1]+=
                    app1*dcAng[ang_jik][1][ngj];
                bt_pi[nb_ij].dAA[2]+=
                    app1*dcAng[ang_jik][2][ngj];
                bt_pi[nb_ij].dBB[0]+=
                    app2*dcAng[ang_jik][0][ngj];
                bt_pi[nb_ij].dBB[1]+=
                    app2*dcAng[ang_jik][1][ngj];
                bt_pi[nb_ij].dBB[2]+=
                    app2*dcAng[ang_jik][2][ngj];
                bt_pi[nb_ik].dAA[0]+=
                    agpdpr1*disij[0][temp_ik]
                    +app1*dcAng[ang_jik][0][ngk];
                bt_pi[nb_ik].dAA[1]+=
                    agpdpr1*disij[1][temp_ik]
                    +app1*dcAng[ang_jik][1][ngk];
                bt_pi[nb_ik].dAA[2]+=
                    agpdpr1*disij[2][temp_ik]
                    +app1*dcAng[ang_jik][2][ngk];
                bt_pi[nb_ik].dBB[0]+=
                    app2*dcAng[ang_jik][0][ngk]
                    +agpdpr2*disij[0][temp_ik];
                bt_pi[nb_ik].dBB[1]+=
                    app2*dcAng[ang_jik][1][ngk]
                    +agpdpr2*disij[1][temp_ik];
                bt_pi[nb_ik].dBB[2]+=
                    app2*dcAng[ang_jik][2][ngk]
                    +agpdpr2*disij[2][temp_ik];

// j and k and k' are different neighbors of i

                for(ltmp=0;ltmp<ktmp;ltmp++) {
                  if(ltmp!=jtmp) {
                    temp_ikp=BOP_index[i]+ltmp;
                    if(neigh_flag[temp_ikp]) {
                      kp=iilist[ltmp];
                      for(nsearch=0;nsearch<nPiBk[n];nsearch++) {
                        ncmp=itypePiBk[n][nsearch];
                        if(x[ncmp][0]==x[kp][0]) {
                          if(x[ncmp][1]==x[kp][1]) {
                            if(x[ncmp][2]==x[kp][2]) {
                              break;
                            }
                          }
                        }
                      }
                      nkikp=ltmp*(2*numneigh[i]-ltmp-1)/2+(ktmp-ltmp)-1;
                      if(jtmp<ltmp) {
                        njikp=jtmp*(2*numneigh[i]-jtmp-1)/2+(ltmp-jtmp)-1;
                        nglj=0;
                        ngl=1;
                      }
                      else {
                        njikp=ltmp*(2*numneigh[i]-ltmp-1)/2+(jtmp-ltmp)-1;
                        nglj=1;
                        ngl=0;
                      }
                      ang_jikp=cos_index[i]+njikp;
                      if(ang_jikp>=cos_total) {
                        error->one(FLERR,"Too many atom triplets for pair bop");
                      }
                      nb_ikp=nb_t;
                      nb_t++;
                      if(nb_t>nb_pi) {
                        new_n_tot=nb_pi+maxneigh;
                        grow_pi(nb_pi,new_n_tot);
                        nb_pi=new_n_tot;
                      }
                      bt_pi[nb_ikp].i=i;
                      bt_pi[nb_ikp].j=kp;
                      bt_pi[nb_ikp].temp=temp_ikp;
                      ang_kikp=cos_index[i]+nkikp;
                      if(ang_kikp>=cos_total) {
                        error->one(FLERR,"Too many atom triplets for pair bop");
                      }
                      betaCapSq2=pi_p[itype-1]*betaS[temp_ikp]*betaS[temp_ikp]
                          -betaP[temp_ikp]*betaP[temp_ikp];
                      dbetaCapSq2=2.0*pi_p[itype-1]*betaS[temp_ikp]*dBetaS[temp_ikp]
                          -2.0*betaP[temp_ikp]*dBetaP[temp_ikp];
                      cosSq1=cosAng[ang_jikp]*cosAng[ang_jikp];
                      angFactor=cosAng[ang_kikp]-cosAng[ang_jikp]*cosAng[ang_jik];
                      angFactor1=4.0*angFactor;
                      angFactor2=-angFactor1*cosAng[ang_jikp]
                          +2.0*cosAng[ang_jik]*(1.0-cosSq1);
                      angFactor3=-angFactor1*cosAng[ang_jik]
                          +2.0*cosAng[ang_jikp]*(1.0-cosSq);
                      angFactor4=2.0*angFactor*angFactor-(1.0-cosSq)*(1.0-cosSq1);
                      betaCapSum=.5*betaCapSq1*betaCapSq2;

//2nd BB is third term of Eq. 38 (a) where j , k and k'=neighbors i

                      BB=BB+betaCapSum*angFactor4;

//agpdpr1 is derivative of BB w.r.t. for atom i w.r.t. Beta(r_ik)
//agpdpr2 is derivative of BB w.r.t. for atom i w.r.t. Beta(r_ik')
//app1 is derivative of BB 3rd term w.r.t. cos(theta_kik')
//app2 is derivative of BB 3rd term w.r.t. cos(theta_jik)
//app3 is derivative of BB 3rd term w.r.t. cos(theta_jik')

                      app1=betaCapSum*angFactor1;
                      app2=betaCapSum*angFactor2;
                      app3=betaCapSum*angFactor3;
                      agpdpr1=.5*angFactor4*dbetaCapSq1*betaCapSq2/rij[temp_ik];
                      agpdpr2=.5*angFactor4*betaCapSq1*dbetaCapSq2/rij[temp_ikp];

                      bt_pi[nb_ij].dBB[0]+=
                          app2*dcAng[ang_jik][0][ngj]
                          +app3*dcAng[ang_jikp][0][nglj];
                      bt_pi[nb_ij].dBB[1]+=
                          app2*dcAng[ang_jik][1][ngj]
                          +app3*dcAng[ang_jikp][1][nglj];
                      bt_pi[nb_ij].dBB[2]+=
                          app2*dcAng[ang_jik][2][ngj]
                          +app3*dcAng[ang_jikp][2][nglj];
                      bt_pi[nb_ik].dBB[0]+=
                          agpdpr1*disij[0][temp_ik]
                          +app1*dcAng[ang_kikp][0][1]
                          +app2*dcAng[ang_jik][0][ngk];
                      bt_pi[nb_ik].dBB[1]+=
                          agpdpr1*disij[1][temp_ik]
                          +app1*dcAng[ang_kikp][1][1]
                          +app2*dcAng[ang_jik][1][ngk];
                      bt_pi[nb_ik].dBB[2]+=
                          agpdpr1*disij[2][temp_ik]
                          +app1*dcAng[ang_kikp][2][1]
                          +app2*dcAng[ang_jik][2][ngk];
                      bt_pi[nb_ikp].dBB[0]+=
                          agpdpr2*disij[0][temp_ikp]
                          +app1*dcAng[ang_kikp][0][0]
                          +app3*dcAng[ang_jikp][0][ngl];
                      bt_pi[nb_ikp].dBB[1]+=
                          agpdpr2*disij[1][temp_ikp]
                          +app1*dcAng[ang_kikp][1][0]
                          +app3*dcAng[ang_jikp][1][ngl];
                      bt_pi[nb_ikp].dBB[2]+=
                          agpdpr2*disij[2][temp_ikp]
                          +app1*dcAng[ang_kikp][2][0]
                          +app3*dcAng[ang_jikp][2][ngl];
                    }
                  }
                }
                nPiBk[n]=nPiBk[n]+1;
              }
            }
          }

//j is a neighbor of i and k is a neighbor of j and equal to i

          for(ki=0;ki<numneigh[j];ki++) {
            k=jlist[ki];
            if(x[k][0]==x[i][0]) {
              if(x[k][1]==x[i][1]) {
                if(x[k][2]==x[i][2]) {
                  break;
                }
              }
            }
          }

//j is a neighbor of i and k is a neighbor of j not equal to i

          for(ktmp=0;ktmp<numneigh[j];ktmp++) {
            if(ktmp!=ki) {
              temp_jk=BOP_index[j]+ktmp;
              if(neigh_flag[temp_jk]) {
                k=jlist[ktmp];
                pi_flag=0;
                for(nsearch=0;nsearch<nPiBk[n];nsearch++) {
                  ncmp=itypePiBk[n][nsearch];
                  if(x[ncmp][0]==x[k][0]) {
                    if(x[ncmp][1]==x[k][1]) {
                      if(x[ncmp][2]==x[k][2]) {
                        pi_flag=1;
                        break;
                      }
                    }
                  }
                }
                if(pi_flag==0) {
                  itypePiBk[n][nPiBk[n]]=k;
                }
                if(ktmp<ki) {
                  njik=ktmp*(2*numneigh[j]-ktmp-1)/2+(ki-ktmp)-1;
                  ngi=1;
                  ngk=0;
                }
                else {
                  njik=ki*(2*numneigh[j]-ki-1)/2+(ktmp-ki)-1;
                  ngi=0;
                  ngk=1;
                }
                ang_ijk=cos_index[j]+njik;
                if(ang_ijk>=cos_total) {
                  error->one(FLERR,"Too many atom triplets for pair bop");
                }
                nb_jk=nb_t;
                nb_t++;
                if(nb_t>nb_pi) {
                  new_n_tot=nb_pi+maxneigh;
                  grow_pi(nb_pi,new_n_tot);
                  nb_pi=new_n_tot;
                }
                bt_pi[nb_jk].i=j;
                bt_pi[nb_jk].j=k;
                bt_pi[nb_jk].temp=temp_jk;
                cosSq=cosAng[ang_ijk]*cosAng[ang_ijk];
                sinFactor=.5*(1.0-cosSq)*pi_p[jtype-1]*betaS[temp_jk];
                cosFactor=.5*(1.0+cosSq)*betaP[temp_jk];
                betaCapSq1=pi_p[jtype-1]*betaS[temp_jk]*betaS[temp_jk]
                    -betaP[temp_jk]*betaP[temp_jk];
                dbetaCapSq1=2.0*pi_p[jtype-1]*betaS[temp_jk]*dBetaS[temp_jk]
                    -2.0*betaP[temp_jk]*dBetaP[temp_jk];

//AA is Eq. 37 (a) and Eq. 19 (b) for j atoms
//3rd BB is 2nd term of Eq. 38 (a) where i and k =neighbors j

                AA=AA+sinFactor*betaS[temp_jk]+cosFactor*betaP[temp_jk];
                BB=BB+.25*(1.0-cosSq)*(1.0-cosSq)*betaCapSq1*betaCapSq1;

//agpdpr1 is derivative of AA for atom j w.r.t. Beta(r_jk)
//agpdpr2 is derivative of BB for atom j w.r.t. Beta(r_jk)
//app1 is derivative of AA for j atom w.r.t. cos(theta_ijk)
//app2 is derivative of BB 2nd term w.r.t. cos(theta_ijk)

                agpdpr1=(2.0*sinFactor*dBetaS[temp_jk]+2.0*cosFactor
                    *dBetaP[temp_jk])/rij[temp_jk];
                agpdpr2=.5*(1.0-cosSq)*(1.0-cosSq)*betaCapSq1*dbetaCapSq1/rij[temp_jk];
                app1=cosAng[ang_ijk]*(-pi_p[jtype-1]*betaS[temp_jk]*betaS[temp_jk]
                    +betaP[temp_jk]*betaP[temp_jk]);
                app2=-(1.0-cosSq)*cosAng[ang_ijk]*betaCapSq1*betaCapSq1;
                bt_pi[nb_ij].dAA[0]-=
                    app1*dcAng[ang_ijk][0][ngi];
                bt_pi[nb_ij].dAA[1]-=
                    app1*dcAng[ang_ijk][1][ngi];
                bt_pi[nb_ij].dAA[2]-=
                    app1*dcAng[ang_ijk][2][ngi];
                bt_pi[nb_ij].dBB[0]-=
                    app2*dcAng[ang_ijk][0][ngi];
                bt_pi[nb_ij].dBB[1]-=
                    app2*dcAng[ang_ijk][1][ngi];
                bt_pi[nb_ij].dBB[2]-=
                    app2*dcAng[ang_ijk][2][ngi];
                bt_pi[nb_jk].dAA[0]+=
                    agpdpr1*disij[0][temp_jk]
                    +app1*dcAng[ang_ijk][0][ngk];
                bt_pi[nb_jk].dAA[1]+=
                    agpdpr1*disij[1][temp_jk]
                    +app1*dcAng[ang_ijk][1][ngk];
                bt_pi[nb_jk].dAA[2]+=
                    agpdpr1*disij[2][temp_jk]
                    +app1*dcAng[ang_ijk][2][ngk];
                bt_pi[nb_jk].dBB[0]+=
                    app2*dcAng[ang_ijk][0][ngk]
                    +agpdpr2*disij[0][temp_jk];
                bt_pi[nb_jk].dBB[1]+=
                    app2*dcAng[ang_ijk][1][ngk]
                    +agpdpr2*disij[1][temp_jk];
                bt_pi[nb_jk].dBB[2]+=
                    app2*dcAng[ang_ijk][2][ngk]
                    +agpdpr2*disij[2][temp_jk];

//j is a neighbor of i and k and k' are different neighbors of j not equal to i

                for(ltmp=0;ltmp<ktmp;ltmp++) {
                  if(ltmp!=ki) {
                    temp_jkp=BOP_index[j]+ltmp;
                    if(neigh_flag[temp_jkp]) {
                      kp=jlist[ltmp];
                      for(nsearch=0;nsearch<nPiBk[n];nsearch++) {
                        ncmp=itypePiBk[n][nsearch];
                        if(x[ncmp][0]==x[kp][0]) {
                          if(x[ncmp][1]==x[kp][1]) {
                            if(x[ncmp][2]==x[kp][2]) {
                              break;
                            }
                          }
                        }
                      }
                      nkjkp=ltmp*(2*numneigh[j]-ltmp-1)/2+(ktmp-ltmp)-1;
                      if(ki<ltmp) {
                        nijkp=ki*(2*numneigh[j]-ki-1)/2+(ltmp-ki)-1;
                        ngli=0;
                        ngl=1;
                      }
                      else {
                        nijkp=ltmp*(2*numneigh[j]-ltmp-1)/2+(ki-ltmp)-1;
                        ngli=1;
                        ngl=0;
                      }
                      ang_ijkp=cos_index[j]+nijkp;
                      if(ang_ijkp>=cos_total) {
                        error->one(FLERR,"Too many atom triplets for pair bop");
                      }
                      ang_kjkp=cos_index[j]+nkjkp;
                      if(ang_kjkp>=cos_total) {
                        error->one(FLERR,"Too many atom triplets for pair bop");
                      }
                      nb_jkp=nb_t;
                      nb_t++;
                      if(nb_t>nb_pi) {
                        new_n_tot=nb_pi+maxneigh;
                        grow_pi(nb_pi,new_n_tot);
                        nb_pi=new_n_tot;
                      }
                      bt_pi[nb_jkp].i=j;
                      bt_pi[nb_jkp].j=kp;
                      bt_pi[nb_jkp].temp=temp_jkp;
                      betaCapSq2=pi_p[jtype-1]*betaS[temp_jkp]*betaS[temp_jkp]
                          -betaP[temp_jkp]*betaP[temp_jkp];
                      dbetaCapSq2=2.0*pi_p[jtype-1]*betaS[temp_jkp]*dBetaS[temp_jkp]
                          -2.0*betaP[temp_jkp]*dBetaP[temp_jkp];
                      cosSq1=cosAng[ang_ijkp]*cosAng[ang_ijkp];
                      angFactor=cosAng[ang_kjkp]-cosAng[ang_ijkp]*cosAng[ang_ijk];
                      angFactor1=4.0*angFactor;
                      angFactor2=-angFactor1*cosAng[ang_ijkp]
                          +2.0*cosAng[ang_ijk]*(1.0-cosSq1);
                      angFactor3=-angFactor1*cosAng[ang_ijk]
                          +2.0*cosAng[ang_ijkp]*(1.0-cosSq);
                      angFactor4=2.0*angFactor*angFactor-(1.0-cosSq)*(1.0-cosSq1);
                      betaCapSum=.5*betaCapSq1*betaCapSq2;

//4th BB is 4th term of Eq. 38 (a) where i , k and k' =neighbors j

                      BB=BB+betaCapSum*angFactor4;

//app1 is derivative of BB 4th term w.r.t. cos(theta_kjk')
//app2 is derivative of BB 4th term w.r.t. cos(theta_ijk)
//app3 is derivative of BB 4th term w.r.t. cos(theta_ijk')
//agpdpr1 is derivative of BB 4th term for atom j w.r.t. Beta(r_jk)
//agpdpr2 is derivative of BB 4th term for atom j w.r.t. Beta(r_jk')

                      app1=betaCapSum*angFactor1;
                      app2=betaCapSum*angFactor2;
                      app3=betaCapSum*angFactor3;
                      agpdpr1=.5*angFactor4*dbetaCapSq1*betaCapSq2/rij[temp_jk];
                      agpdpr2=.5*angFactor4*betaCapSq1*dbetaCapSq2/rij[temp_jkp];

                      bt_pi[nb_ij].dBB[0]-=
                          app3*dcAng[ang_ijkp][0][ngli]
                          +app2*dcAng[ang_ijk][0][ngi];
                      bt_pi[nb_ij].dBB[1]-=
                          app3*dcAng[ang_ijkp][1][ngli]
                          +app2*dcAng[ang_ijk][1][ngi];
                      bt_pi[nb_ij].dBB[2]-=
                          app3*dcAng[ang_ijkp][2][ngli]
                          +app2*dcAng[ang_ijk][2][ngi];
                      bt_pi[nb_jk].dBB[0]+=
                          agpdpr1*disij[0][temp_jk]
                          +app1*dcAng[ang_kjkp][0][1]
                          +app2*dcAng[ang_ijk][0][ngk];
                      bt_pi[nb_jk].dBB[1]+=
                          agpdpr1*disij[1][temp_jk]
                          +app1*dcAng[ang_kjkp][1][1]
                          +app2*dcAng[ang_ijk][1][ngk];
                      bt_pi[nb_jk].dBB[2]+=
                          agpdpr1*disij[2][temp_jk]
                          +app1*dcAng[ang_kjkp][2][1]
                          +app2*dcAng[ang_ijk][2][ngk];
                      bt_pi[nb_jkp].dBB[0]+=
                          agpdpr2*disij[0][temp_jkp]
                          +app1*dcAng[ang_kjkp][0][0]
                          +app3*dcAng[ang_ijkp][0][ngl];
                      bt_pi[nb_jkp].dBB[1]+=
                          agpdpr2*disij[1][temp_jkp]
                          +app1*dcAng[ang_kjkp][1][0]
                          +app3*dcAng[ang_ijkp][1][ngl];
                      bt_pi[nb_jkp].dBB[2]+=
                          agpdpr2*disij[2][temp_jkp]
                          +app1*dcAng[ang_kjkp][2][0]
                          +app3*dcAng[ang_ijkp][2][ngl];
                    }
                  }
                }

//j and k' are different neighbors of i and k is a neighbor of j not equal to i

                for(ltmp=0;ltmp<numneigh[i];ltmp++) {
                  if(ltmp!=jtmp) {
                    temp_ikp=BOP_index[i]+ltmp;
                    if(neigh_flag[temp_ikp]) {
                      kp=iilist[ltmp];
                      for(nsearch=0;nsearch<nPiBk[n];nsearch++) {
                        ncmp=itypePiBk[n][nsearch];
                        if(x[ncmp][0]==x[kp][0]) {
                          if(x[ncmp][1]==x[kp][1]) {
                            if(x[ncmp][2]==x[kp][2]) {
                              break;
                            }
                          }
                        }
                      }
                      if(ltmp<jtmp) {
                        njikp=ltmp*(2*numneigh[i]-ltmp-1)/2+(jtmp-ltmp)-1;
                        ngl=1;
                        nglj=0;
                      }
                      else {
                        njikp=jtmp*(2*numneigh[i]-jtmp-1)/2+(ltmp-jtmp)-1;
                        ngl=0;
                        nglj=1;
                      }
                      ang_jikp=cos_index[i]+njikp;
                      if(ang_jikp>=cos_total) {
                        error->one(FLERR,"Too many atom triplets for pair bop");
                      }
                      nb_ikp=nb_t;
                      nb_t++;
                      if(nb_t>nb_pi) {
                        new_n_tot=nb_pi+maxneigh;
                        grow_pi(nb_pi,new_n_tot);
                        nb_pi=new_n_tot;
                      }
                      bt_pi[nb_ikp].i=i;
                      bt_pi[nb_ikp].j=kp;
                      bt_pi[nb_ikp].temp=temp_ikp;
                      betaCapSq2=pi_p[itype-1]*betaS[temp_ikp]*betaS[temp_ikp]
                          -betaP[temp_ikp]*betaP[temp_ikp];
                      dbetaCapSq2=2.0*pi_p[itype-1]*betaS[temp_ikp]*dBetaS[temp_ikp]
                          -2.0*betaP[temp_ikp]*dBetaP[temp_ikp];
                      dotV=(disij[0][temp_jk]*disij[0][temp_ikp]+disij[1][temp_jk]
                          *disij[1][temp_ikp]+disij[2][temp_jk]*disij[2][temp_ikp])
                          /(rij[temp_jk]*rij[temp_ikp]);
                      cosSq1=cosAng[ang_jikp]*cosAng[ang_jikp];
                      angFactor=dotV+cosAng[ang_jikp]*cosAng[ang_ijk];
                      angRfactor=4.0*angFactor*dotV;
                      dAngR1=-angRfactor/rij[temp_jk];
                      dAngR2=-angRfactor/rij[temp_ikp];
                      angFactor1=4.0*angFactor*cosAng[ang_jikp]
                          +2.0*cosAng[ang_ijk]*(1.0-cosSq1);
                      angFactor2=4.0*angFactor*cosAng[ang_ijk]
                          +2.0*cosAng[ang_jikp]*(1.0-cosSq);
                      angFactor3=2.0*angFactor*angFactor-(1.0-cosSq)*(1.0-cosSq1);
                      betaCapSum=.5*betaCapSq1*betaCapSq2;

//5th BB is 5th term of Eq. 38 (a) Eq. 21 (b) where i , k and k' =neighbors j

                      BB=BB+betaCapSum*angFactor3;

//app1 is derivative of BB 5th term w.r.t. cos(theta_ijk)
//app2 is derivative of BB 5th term w.r.t. cos(theta_jik')
//agpdpr1 is derivative of BB 5th term for atom j w.r.t. Beta(r_jk)
//agpdpr2 is derivative of BB 5th term for atom j w.r.t. Beta(r_ik')
//agpdpr3 is derivative of BB 5th term for atom j w.r.t. dot(r_ik',r_ij)

                      app1=betaCapSum*angFactor1;
                      app2=betaCapSum*angFactor2;
                      agpdpr1=(.5*angFactor3*dbetaCapSq1*betaCapSq2
                          +betaCapSum*dAngR1)/rij[temp_jk];
                      agpdpr2=(.5*angFactor3*betaCapSq1*dbetaCapSq2
                          +betaCapSum*dAngR2)/rij[temp_ikp];
                      agpdpr3=4.0*betaCapSum*angFactor/(rij[temp_ikp]*rij[temp_jk]);

                      bt_pi[nb_ij].dBB[0]+=
                          +app2*dcAng[ang_jikp][0][ngl]
                          -app1*dcAng[ang_ijk][0][ngi];
                      bt_pi[nb_ij].dBB[1]+=
                          +app2*dcAng[ang_jikp][1][ngl]
                          -app1*dcAng[ang_ijk][1][ngi];
                      bt_pi[nb_ij].dBB[2]+=
                          +app2*dcAng[ang_jikp][2][ngl]
                          -app1*dcAng[ang_ijk][2][ngi];
                      bt_pi[nb_ikp].dBB[0]+=
                          agpdpr2*disij[0][temp_ikp]
                          +agpdpr3*disij[0][temp_jk]
                          +app2*dcAng[ang_jikp][0][nglj];
                      bt_pi[nb_ikp].dBB[1]+=
                          agpdpr2*disij[1][temp_ikp]
                          +agpdpr3*disij[1][temp_jk]
                          +app2*dcAng[ang_jikp][1][nglj];
                      bt_pi[nb_ikp].dBB[2]+=
                          agpdpr2*disij[2][temp_ikp]
                          +agpdpr3*disij[2][temp_jk]
                          +app2*dcAng[ang_jikp][2][nglj];
                      bt_pi[nb_jk].dBB[0]+=
                          agpdpr1*disij[0][temp_jk]
                          +agpdpr3*disij[0][temp_ikp]
                          +app1*dcAng[ang_ijk][0][ngk];
                      bt_pi[nb_jk].dBB[1]+=
                          agpdpr1*disij[1][temp_jk]
                          +agpdpr3*disij[1][temp_ikp]
                          +app1*dcAng[ang_ijk][1][ngk];
                      bt_pi[nb_jk].dBB[2]+=
                          agpdpr1*disij[2][temp_jk]
                          +agpdpr3*disij[2][temp_ikp]
                          +app1*dcAng[ang_ijk][2][ngk];
                    }
                  }
                }
                if(pi_flag==0)
                  nPiBk[n]=nPiBk[n]+1;
              }
            }
          }
          CC=betaP[temp_ij]*betaP[temp_ij]+pi_delta[iij]*pi_delta[iij];
          BBrt=sqrt(BB+small6);
          AB1=CC+pi_c[iij]*(AA+BBrt)+small7;
          AB2=CC+pi_c[iij]*(AA-BBrt+sqrt(small6))+small7;
          BBrtR=1.0/BBrt;
          ABrtR1=1.0/sqrt(AB1);
          ABrtR2=1.0/sqrt(AB2);

// piB is similary formulation to (a) Eq. 36 and (b) Eq. 18

          piB[n]=(ABrtR1+ABrtR2)*pi_a[iij]*betaP[temp_ij];
          dPiB1=-.5*(cube(ABrtR1)+cube(ABrtR2))*pi_c[iij]*pi_a[iij]*betaP[temp_ij];
          dPiB2=.25*BBrtR*(cube(ABrtR2)-cube(ABrtR1))*pi_c[iij]*pi_a[iij]*betaP[temp_ij];
          dPiB3=((ABrtR1+ABrtR2)*pi_a[iij]-(cube(ABrtR1)+cube(ABrtR2))*pi_a[iij]
              *betaP[temp_ij]*betaP[temp_ij])*dBetaP[temp_ij]/rij[temp_ij];
          n++;
          pp2=2.0*betaP[temp_ij];
          for(m=0;m<nb_t;m++) {
            bt_ij=bt_pi[m].temp;
            bt_i=bt_pi[m].i;
            bt_j=bt_pi[m].j;
            for(pp=0;pp<3;pp++) {
              bt_pi[m].dPiB[pp]=
                  +dPiB1*bt_pi[m].dAA[pp]
                  +dPiB2*bt_pi[m].dBB[pp];
              ftmp[pp]=pp2*bt_pi[m].dPiB[pp];
              f[bt_i][pp]-=ftmp[pp];
              f[bt_j][pp]+=ftmp[pp];

            }
            if(evflag) {
              ev_tally_xyz(bt_i,bt_j,nlocal,newton_pair,0.0,0.0,ftmp[0],ftmp[1]
                  ,ftmp[2],disij[0][bt_ij],disij[1][bt_ij],disij[2][bt_ij]);
            }
          }
          for(pp=0;pp<3;pp++) {
            ftmp[pp]=pp2*dPiB3*disij[pp][temp_ij];
            f[i][pp]-=ftmp[pp];
            f[j][pp]+=ftmp[pp];
          }
          if(evflag) {
            ev_tally_xyz(i,j,nlocal,newton_pair,0.0,0.0,ftmp[0],ftmp[1]
                ,ftmp[2],disij[0][temp_ij],disij[1][temp_ij],disij[2][temp_ij]);
          }
        }
      }
    }
  }
  destroy_pi();
}

/* ---------------------------------------------------------------------- */

void PairBOP::PiBo_otf()
{
  int new_n_tot;
  int i,j,k,kp,m,n,pp,nb_t;
  int iij,iik,iikp,ji,ki,ijkp,ijk;
  int nsearch,ncmp;
  tagint i_tag,j_tag;
  int itmp,ltmp,jtmp,ktmp;
  int pi_flag,ks;
  int nlocal;
  int inum,*ilist,*iilist,*jlist;
  int **firstneigh,*numneigh;
  int itype,jtype,ktype,kptype;
  int temp_ij,temp_ik,temp_ikp;
  int temp_jk,temp_jkp;
  int nb_ij,nb_ik,nb_jk,nb_ikp,nb_jkp;
  int bt_i,bt_j;
  double AA,BB,CC;
  double cosSq,sinFactor,cosFactor;
  double cosSq1,dotV,BBrt,AB1,AB2;
  double BBrtR,ABrtR1,ABrtR2;
  double angFactor,angFactor1,angFactor2;
  double angFactor3,angFactor4,angRfactor;
  double dAngR1,dAngR2,agpdpr3;
  double agpdpr1,agpdpr2,app1,app2,app3;
  double betaCapSq1,dbetaCapSq1;
  double betaCapSq2,dbetaCapSq2;
  double betaCapSum,ps;
  double ftmp[3],xtmp[3];
  double dPiB1,dPiB2,dPiB3,pp2;

  double dis_ij[3],rsq_ij,r_ij;
  double betaP_ij,dBetaP_ij;
  double dis_ik[3],rsq_ik,r_ik;
  double betaS_ik,dBetaS_ik;
  double betaP_ik,dBetaP_ik;
  double dis_ikp[3],rsq_ikp,r_ikp;
  double betaS_ikp,dBetaS_ikp;
  double betaP_ikp,dBetaP_ikp;
  double dis_jk[3],rsq_jk,r_jk;
  double betaS_jk,dBetaS_jk;
  double betaP_jk,dBetaP_jk;
  double dis_jkp[3],rsq_jkp,r_jkp;
  double betaS_jkp,dBetaS_jkp;
  double betaP_jkp,dBetaP_jkp;

  double cosAng_jik,dcA_jik[3][2];
  double cosAng_jikp,dcA_jikp[3][2];
  double cosAng_kikp,dcA_kikp[3][2];
  double cosAng_ijk,dcA_ijk[3][2];
  double cosAng_ijkp,dcA_ijkp[3][2];
  double cosAng_kjkp,dcA_kjkp[3][2];

  int newton_pair = force->newton_pair;

  double **f = atom->f;
  double **x = atom->x;
  int *type = atom->type;
  tagint *tag = atom->tag;

  nlocal = atom->nlocal;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  inum = list->inum;
  ilist = list->ilist;
  n=0;
  if(nb_pi>16) {
    nb_pi=16;
  }
  if(nb_pi==0) {
    nb_pi=(maxneigh)*(maxneigh/2);
  }

// Loop over all local atoms for i

  if(allocate_pi) {
    destroy_pi();
  }
  create_pi(nb_pi);

  for(itmp=0;itmp<inum;itmp++) {
    nb_t=0;
    i = ilist[itmp];
    itype = map[type[i]]+1;
    i_tag=tag[i];

// j is a loop over all neighbors of i

    iilist=firstneigh[i];
    for(jtmp=0;jtmp<numneigh[i];jtmp++) {
      for(m=0;m<nb_pi;m++) {
        for(pp=0;pp<3;pp++) {
          bt_pi[m].dAA[pp]=0.0;
          bt_pi[m].dBB[pp]=0.0;
          bt_pi[m].dPiB[pp]=0.0;
        }
        bt_pi[m].i=-1;
        bt_pi[m].j=-1;
      }
      temp_ij=BOP_index[i]+jtmp;
      j=iilist[jtmp];
      jlist=firstneigh[j];
      jtype=map[type[j]]+1;
      j_tag=tag[j];
      nb_t=0;
      ftmp[0]=0.0;
      ftmp[1]=0.0;
      ftmp[2]=0.0;
      if(j_tag>=i_tag) {
        if(itype==jtype)
          iij=itype-1;
        else if(itype<jtype)
          iij=itype*bop_types-itype*(itype+1)/2+jtype-1;
        else
          iij=jtype*bop_types-jtype*(jtype+1)/2+itype-1;
        AA=0.0;
        BB=0.0;
        nPiBk[n]=0;
        for(ji=0;ji<numneigh[j];ji++) {
          if(x[jlist[ji]][0]==x[i][0]) {
            if(x[jlist[ji]][1]==x[i][1]) {
              if(x[jlist[ji]][2]==x[i][2]) {
                  break;
              }
            }
          }
        }
        nb_ij=nb_t;
        nb_t++;
        if(nb_t>nb_pi) {
          new_n_tot=nb_pi+maxneigh;
          grow_pi(nb_pi,new_n_tot);
          nb_pi=new_n_tot;
        }
        bt_pi[nb_ij].i=i;
        bt_pi[nb_ij].j=j;
        bt_pi[nb_ij].temp=temp_ij;
        dis_ij[0]=x[j][0]-x[i][0];
        dis_ij[1]=x[j][1]-x[i][1];
        dis_ij[2]=x[j][2]-x[i][2];
        rsq_ij=dis_ij[0]*dis_ij[0]
            +dis_ij[1]*dis_ij[1]
            +dis_ij[2]*dis_ij[2];
        r_ij=sqrt(rsq_ij);
        if(r_ij<=rcut[iij]) {
          ps=r_ij*rdr[iij]+1.0;
          ks=(int)ps;
          if(nr-1<ks)
            ks=nr-1;
          ps=ps-ks;
          if(ps>1.0)
            ps=1.0;
          betaP_ij=((pBetaP3[iij][ks-1]*ps+pBetaP2[iij][ks-1])*ps
              +pBetaP1[iij][ks-1])*ps+pBetaP[iij][ks-1];
          dBetaP_ij=(pBetaP6[iij][ks-1]*ps+pBetaP5[iij][ks-1])*ps
              +pBetaP4[iij][ks-1];

// j and k are different neighbors of i

          for(ktmp=0;ktmp<numneigh[i];ktmp++) {
            if(ktmp!=jtmp) {
              temp_ik=BOP_index[i]+ktmp;
              k=iilist[ktmp];
              ktype=map[type[k]]+1;
              if(itype==ktype)
                iik=itype-1;
              else if(itype<ktype)
                iik=itype*bop_types-itype*(itype+1)/2+ktype-1;
              else
                iik=ktype*bop_types-ktype*(ktype+1)/2+itype-1;
              dis_ik[0]=x[k][0]-x[i][0];
              dis_ik[1]=x[k][1]-x[i][1];
              dis_ik[2]=x[k][2]-x[i][2];
              rsq_ik=dis_ik[0]*dis_ik[0]
                  +dis_ik[1]*dis_ik[1]
                  +dis_ik[2]*dis_ik[2];
              r_ik=sqrt(rsq_ik);
              if(r_ik<=rcut[iik]) {
                ps=r_ik*rdr[iik]+1.0;
                ks=(int)ps;
                if(nr-1<ks)
                  ks=nr-1;
                ps=ps-ks;
                if(ps>1.0)
                  ps=1.0;
                betaS_ik=((pBetaS3[iik][ks-1]*ps+pBetaS2[iik][ks-1])*ps
                    +pBetaS1[iik][ks-1])*ps+pBetaS[iik][ks-1];
                dBetaS_ik=(pBetaS6[iik][ks-1]*ps+pBetaS5[iik][ks-1])*ps
                    +pBetaS4[iik][ks-1];
                betaP_ik=((pBetaP3[iik][ks-1]*ps+pBetaP2[iik][ks-1])*ps
                    +pBetaP1[iik][ks-1])*ps+pBetaP[iik][ks-1];
                dBetaP_ik=(pBetaP6[iik][ks-1]*ps+pBetaP5[iik][ks-1])*ps
                    +pBetaP4[iik][ks-1];
                cosAng_jik=(dis_ij[0]*dis_ik[0]+dis_ij[1]*dis_ik[1]
                    +dis_ij[2]*dis_ik[2])/(r_ij*r_ik);
                dcA_jik[0][0]=(dis_ik[0]*r_ij*r_ik-cosAng_jik
                    *dis_ij[0]*r_ik*r_ik)/(r_ij*r_ij*r_ik*r_ik);
                dcA_jik[1][0]=(dis_ik[1]*r_ij*r_ik-cosAng_jik
                    *dis_ij[1]*r_ik*r_ik)/(r_ij*r_ij*r_ik*r_ik);
                dcA_jik[2][0]=(dis_ik[2]*r_ij*r_ik-cosAng_jik
                    *dis_ij[2]*r_ik*r_ik)/(r_ij*r_ij*r_ik*r_ik);
                dcA_jik[0][1]=(dis_ij[0]*r_ij*r_ik-cosAng_jik
                    *dis_ik[0]*r_ij*r_ij)/(r_ij*r_ij*r_ik*r_ik);
                dcA_jik[1][1]=(dis_ij[1]*r_ij*r_ik-cosAng_jik
                    *dis_ik[1]*r_ij*r_ij)/(r_ij*r_ij*r_ik*r_ik);
                dcA_jik[2][1]=(dis_ij[2]*r_ij*r_ik-cosAng_jik
                    *dis_ik[2]*r_ij*r_ij)/(r_ij*r_ij*r_ik*r_ik);
                nb_ik=nb_t;
                nb_t++;
                if(nb_t>nb_pi) {
                  new_n_tot=nb_pi+maxneigh;
                  grow_pi(nb_pi,new_n_tot);
                  nb_pi=new_n_tot;
                }
                bt_pi[nb_ik].i=i;
                bt_pi[nb_ik].j=k;
                bt_pi[nb_ik].temp=temp_ik;
                cosSq=cosAng_jik*cosAng_jik;
                sinFactor=.5*(1.0-cosSq)*pi_p[itype-1]*betaS_ik;
                cosFactor=.5*(1.0+cosSq)*betaP_ik;
                betaCapSq1=pi_p[itype-1]*betaS_ik*betaS_ik-betaP_ik
                    *betaP_ik;
                dbetaCapSq1=2.0*pi_p[itype-1]*betaS_ik*dBetaS_ik
                    -2.0*betaP_ik*dBetaP_ik;

//AA is Eq. 37 (a) and Eq. 19 (b) or i atoms
//1st BB is first term of Eq. 38 (a) where j and k =neighbors i

                AA=AA+sinFactor*betaS_ik+cosFactor*betaP_ik;
                BB=BB+.25*(1.0-cosSq)*(1.0-cosSq)*betaCapSq1*betaCapSq1;

//agpdpr1 is derivative of AA w.r.t. for atom i w.r.t. Beta(r_ik)
//agpdpr2 is derivative of BB w.r.t. for atom i w.r.t. Beta(r_ik)
//app1 is derivative of AA w.r.t. for atom i w.r.t. cos(theta_jik)
//app2 is derivative of BB w.r.t. for atom i w.r.t. cos(theta_jik)

                agpdpr1=(2.0*sinFactor*dBetaS_ik+2.0*cosFactor
                    *dBetaP_ik)/r_ik;
                app1=cosAng_jik*(-pi_p[itype-1]*betaS_ik*betaS_ik
                    +betaP_ik*betaP_ik);
                app2=-(1.0-cosSq)*cosAng_jik*betaCapSq1*betaCapSq1;
                agpdpr2=.5*(1.0-cosSq)*(1.0-cosSq)*betaCapSq1*dbetaCapSq1/r_ik;
                itypePiBk[n][nPiBk[n]]=k;
                bt_pi[nb_ij].dAA[0]+=
                    app1*dcA_jik[0][0];
                bt_pi[nb_ij].dAA[1]+=
                    app1*dcA_jik[1][0];
                bt_pi[nb_ij].dAA[2]+=
                    app1*dcA_jik[2][0];
                bt_pi[nb_ij].dBB[0]+=
                    app2*dcA_jik[0][0];
                bt_pi[nb_ij].dBB[1]+=
                    app2*dcA_jik[1][0];
                bt_pi[nb_ij].dBB[2]+=
                    app2*dcA_jik[2][0];
                bt_pi[nb_ik].dAA[0]+=
                    agpdpr1*dis_ik[0]
                    +app1*dcA_jik[0][1];
                bt_pi[nb_ik].dAA[1]+=
                    agpdpr1*dis_ik[1]
                    +app1*dcA_jik[1][1];
                bt_pi[nb_ik].dAA[2]+=
                    agpdpr1*dis_ik[2]
                    +app1*dcA_jik[2][1];
                bt_pi[nb_ik].dBB[0]+=
                    app2*dcA_jik[0][1]
                    +agpdpr2*dis_ik[0];
                bt_pi[nb_ik].dBB[1]+=
                    app2*dcA_jik[1][1]
                    +agpdpr2*dis_ik[1];
                bt_pi[nb_ik].dBB[2]+=
                    app2*dcA_jik[2][1]
                    +agpdpr2*dis_ik[2];

// j and k and k' are different neighbors of i

                for(ltmp=0;ltmp<ktmp;ltmp++) {
                  if(ltmp!=jtmp) {
                    temp_ikp=BOP_index[i]+ltmp;
                    kp=iilist[ltmp];
                    kptype=map[type[kp]]+1;
                    for(nsearch=0;nsearch<nPiBk[n];nsearch++) {
                      ncmp=itypePiBk[n][nsearch];
                      if(x[ncmp][0]==x[kp][0]) {
                        if(x[ncmp][1]==x[kp][1]) {
                          if(x[ncmp][2]==x[kp][2]) {
                            break;
                          }
                        }
                      }
                    }
                    if(itype==kptype)
                      iikp=itype-1;
                    else if(itype<kptype)
                      iikp=itype*bop_types-itype*(itype+1)/2+kptype-1;
                    else
                      iikp=kptype*bop_types-kptype*(kptype+1)/2+itype-1;
                    dis_ikp[0]=x[kp][0]-x[i][0];
                    dis_ikp[1]=x[kp][1]-x[i][1];
                    dis_ikp[2]=x[kp][2]-x[i][2];
                    rsq_ikp=dis_ikp[0]*dis_ikp[0]
                        +dis_ikp[1]*dis_ikp[1]
                        +dis_ikp[2]*dis_ikp[2];
                    r_ikp=sqrt(rsq_ikp);
                    if(r_ikp<=rcut[iikp]) {
                      ps=r_ikp*rdr[iikp]+1.0;
                      ks=(int)ps;
                      if(nr-1<ks)
                        ks=nr-1;
                      ps=ps-ks;
                      if(ps>1.0)
                        ps=1.0;
                      betaS_ikp=((pBetaS3[iikp][ks-1]*ps+pBetaS2[iikp][ks-1])*ps
                          +pBetaS1[iikp][ks-1])*ps+pBetaS[iikp][ks-1];
                      dBetaS_ikp=(pBetaS6[iikp][ks-1]*ps+pBetaS5[iikp][ks-1])*ps
                          +pBetaS4[iikp][ks-1];
                      betaP_ikp=((pBetaP3[iikp][ks-1]*ps+pBetaP2[iikp][ks-1])*ps
                          +pBetaP1[iikp][ks-1])*ps+pBetaP[iikp][ks-1];
                      dBetaP_ikp=(pBetaP6[iikp][ks-1]*ps+pBetaP5[iikp][ks-1])*ps
                          +pBetaP4[iikp][ks-1];
                      cosAng_jikp=(dis_ij[0]*dis_ikp[0]+dis_ij[1]*dis_ikp[1]
                          +dis_ij[2]*dis_ikp[2])/(r_ij*r_ikp);
                      dcA_jikp[0][0]=(dis_ikp[0]*r_ij*r_ikp-cosAng_jikp
                          *dis_ij[0]*r_ikp*r_ikp)/(r_ij*r_ij*r_ikp*r_ikp);
                      dcA_jikp[1][0]=(dis_ikp[1]*r_ij*r_ikp-cosAng_jikp
                          *dis_ij[1]*r_ikp*r_ikp)/(r_ij*r_ij*r_ikp*r_ikp);
                      dcA_jikp[2][0]=(dis_ikp[2]*r_ij*r_ikp-cosAng_jikp
                          *dis_ij[2]*r_ikp*r_ikp)/(r_ij*r_ij*r_ikp*r_ikp);
                      dcA_jikp[0][1]=(dis_ij[0]*r_ij*r_ikp-cosAng_jikp
                          *dis_ikp[0]*r_ij*r_ij)/(r_ij*r_ij*r_ikp*r_ikp);
                      dcA_jikp[1][1]=(dis_ij[1]*r_ij*r_ikp-cosAng_jikp
                          *dis_ikp[1]*r_ij*r_ij)/(r_ij*r_ij*r_ikp*r_ikp);
                      dcA_jikp[2][1]=(dis_ij[2]*r_ij*r_ikp-cosAng_jikp
                          *dis_ikp[2]*r_ij*r_ij)/(r_ij*r_ij*r_ikp*r_ikp);
                      cosAng_kikp=(dis_ik[0]*dis_ikp[0]+dis_ik[1]*dis_ikp[1]
                          +dis_ik[2]*dis_ikp[2])/(r_ik*r_ikp);
                      dcA_kikp[0][0]=(dis_ikp[0]*r_ik*r_ikp-cosAng_kikp
                          *dis_ik[0]*r_ikp*r_ikp)/(r_ik*r_ik*r_ikp*r_ikp);
                      dcA_kikp[1][0]=(dis_ikp[1]*r_ik*r_ikp-cosAng_kikp
                          *dis_ik[1]*r_ikp*r_ikp)/(r_ik*r_ik*r_ikp*r_ikp);
                      dcA_kikp[2][0]=(dis_ikp[2]*r_ik*r_ikp-cosAng_kikp
                          *dis_ik[2]*r_ikp*r_ikp)/(r_ik*r_ik*r_ikp*r_ikp);
                      dcA_kikp[0][1]=(dis_ik[0]*r_ik*r_ikp-cosAng_kikp
                          *dis_ikp[0]*r_ik*r_ik)/(r_ik*r_ik*r_ikp*r_ikp);
                      dcA_kikp[1][1]=(dis_ik[1]*r_ik*r_ikp-cosAng_kikp
                          *dis_ikp[1]*r_ik*r_ik)/(r_ik*r_ik*r_ikp*r_ikp);
                      dcA_kikp[2][1]=(dis_ik[2]*r_ik*r_ikp-cosAng_kikp
                          *dis_ikp[2]*r_ik*r_ik)/(r_ik*r_ik*r_ikp*r_ikp);
                      nb_ikp=nb_t;
                      nb_t++;
                      if(nb_t>nb_pi) {
                        new_n_tot=nb_pi+maxneigh;
                        grow_pi(nb_pi,new_n_tot);
                        nb_pi=new_n_tot;
                      }
                      bt_pi[nb_ikp].i=i;
                      bt_pi[nb_ikp].j=kp;
                      bt_pi[nb_ikp].temp=temp_ikp;
                      betaCapSq2=pi_p[itype-1]*betaS_ikp*betaS_ikp
                          -betaP_ikp*betaP_ikp;
                      dbetaCapSq2=2.0*pi_p[itype-1]*betaS_ikp*dBetaS_ikp
                          -2.0*betaP_ikp*dBetaP_ikp;
                      cosSq1=cosAng_jikp*cosAng_jikp;
                      angFactor=cosAng_kikp-cosAng_jikp*cosAng_jik;
                      angFactor1=4.0*angFactor;
                      angFactor2=-angFactor1*cosAng_jikp
                          +2.0*cosAng_jik*(1.0-cosSq1);
                      angFactor3=-angFactor1*cosAng_jik
                          +2.0*cosAng_jikp*(1.0-cosSq);
                      angFactor4=2.0*angFactor*angFactor-(1.0-cosSq)*(1.0-cosSq1);
                      betaCapSum=.5*betaCapSq1*betaCapSq2;

//2nd BB is third term of Eq. 38 (a) where j , k and k'=neighbors i

                      BB=BB+betaCapSum*angFactor4;

//agpdpr1 is derivative of BB w.r.t. for atom i w.r.t. Beta(r_ik)
//agpdpr2 is derivative of BB w.r.t. for atom i w.r.t. Beta(r_ik')
//app1 is derivative of BB 3rd term w.r.t. cos(theta_kik')
//app2 is derivative of BB 3rd term w.r.t. cos(theta_jik)
//app3 is derivative of BB 3rd term w.r.t. cos(theta_jik')

                      app1=betaCapSum*angFactor1;
                      app2=betaCapSum*angFactor2;
                      app3=betaCapSum*angFactor3;
                      agpdpr1=.5*angFactor4*dbetaCapSq1*betaCapSq2/r_ik;
                      agpdpr2=.5*angFactor4*betaCapSq1*dbetaCapSq2/r_ikp;
                      bt_pi[nb_ij].dBB[0]+=
                          app2*dcA_jik[0][0]
                          +app3*dcA_jikp[0][0];
                      bt_pi[nb_ij].dBB[1]+=
                          app2*dcA_jik[1][0]
                          +app3*dcA_jikp[1][0];
                      bt_pi[nb_ij].dBB[2]+=
                          app2*dcA_jik[2][0]
                          +app3*dcA_jikp[2][0];
                      bt_pi[nb_ik].dBB[0]+=
                          agpdpr1*dis_ik[0]
                          +app1*dcA_kikp[0][0]
                          +app2*dcA_jik[0][1];
                      bt_pi[nb_ik].dBB[1]+=
                          agpdpr1*dis_ik[1]
                          +app1*dcA_kikp[1][0]
                          +app2*dcA_jik[1][1];
                      bt_pi[nb_ik].dBB[2]+=
                          agpdpr1*dis_ik[2]
                          +app1*dcA_kikp[2][0]
                          +app2*dcA_jik[2][1];
                      bt_pi[nb_ikp].dBB[0]+=
                          agpdpr2*dis_ikp[0]
                          +app1*dcA_kikp[0][1]
                          +app3*dcA_jikp[0][1];
                      bt_pi[nb_ikp].dBB[1]+=
                          agpdpr2*dis_ikp[1]
                          +app1*dcA_kikp[1][1]
                          +app3*dcA_jikp[1][1];
                      bt_pi[nb_ikp].dBB[2]+=
                          agpdpr2*dis_ikp[2]
                          +app1*dcA_kikp[2][1]
                          +app3*dcA_jikp[2][1];
                      }
                    }
                  }
                nPiBk[n]=nPiBk[n]+1;
                }
              }
            }

//j is a neighbor of i and k is a neighbor of j and equal to i

          for(ki=0;ki<numneigh[j];ki++) {
            k=jlist[ki];
            if(x[k][0]==x[i][0]) {
              if(x[k][1]==x[i][1]) {
                if(x[k][2]==x[i][2]) {
                  break;
                }
              }
            }
          }

//j is a neighbor of i and k is a neighbor of j not equal to i

          for(ktmp=0;ktmp<numneigh[j];ktmp++) {
            if(ktmp!=ki) {
              temp_jk=BOP_index[j]+ktmp;
              k=jlist[ktmp];
              ktype=map[type[k]]+1;
              pi_flag=0;
              for(nsearch=0;nsearch<nPiBk[n];nsearch++) {
                ncmp=itypePiBk[n][nsearch];
                if(x[ncmp][0]==x[k][0]) {
                  if(x[ncmp][1]==x[k][1]) {
                    if(x[ncmp][2]==x[k][2]) {
                      pi_flag=1;
                      break;
                    }
                  }
                }
              }
              if(pi_flag==0) {
                itypePiBk[n][nPiBk[n]]=k;
              }
              if(jtype==ktype)
                ijk=jtype-1;
              else if(jtype<ktype)
                ijk=jtype*bop_types-jtype*(jtype+1)/2+ktype-1;
              else
                ijk=ktype*bop_types-ktype*(ktype+1)/2+jtype-1;
              dis_jk[0]=x[k][0]-x[j][0];
              dis_jk[1]=x[k][1]-x[j][1];
              dis_jk[2]=x[k][2]-x[j][2];
              rsq_jk=dis_jk[0]*dis_jk[0]
                  +dis_jk[1]*dis_jk[1]
                  +dis_jk[2]*dis_jk[2];
              r_jk=sqrt(rsq_jk);
              if(r_jk<=rcut[ijk]) {
                ps=r_jk*rdr[ijk]+1.0;
                ks=(int)ps;
                if(nr-1<ks)
                  ks=nr-1;
                ps=ps-ks;
                if(ps>1.0)
                  ps=1.0;
                betaS_jk=((pBetaS3[ijk][ks-1]*ps+pBetaS2[ijk][ks-1])*ps
                    +pBetaS1[ijk][ks-1])*ps+pBetaS[ijk][ks-1];
                dBetaS_jk=(pBetaS6[ijk][ks-1]*ps+pBetaS5[ijk][ks-1])*ps
                    +pBetaS4[ijk][ks-1];
                betaP_jk=((pBetaP3[ijk][ks-1]*ps+pBetaP2[ijk][ks-1])*ps
                    +pBetaP1[ijk][ks-1])*ps+pBetaP[ijk][ks-1];
                dBetaP_jk=(pBetaP6[ijk][ks-1]*ps+pBetaP5[ijk][ks-1])*ps
                    +pBetaP4[ijk][ks-1];
                cosAng_ijk=(-dis_ij[0]*dis_jk[0]-dis_ij[1]*dis_jk[1]
                    -dis_ij[2]*dis_jk[2])/(r_ij*r_jk);
                dcA_ijk[0][0]=(dis_jk[0]*r_ij*r_jk-cosAng_ijk
                    *-dis_ij[0]*r_jk*r_jk)/(r_ij*r_ij*r_jk*r_jk);
                dcA_ijk[1][0]=(dis_jk[1]*r_ij*r_jk-cosAng_ijk
                    *-dis_ij[1]*r_jk*r_jk)/(r_ij*r_ij*r_jk*r_jk);
                dcA_ijk[2][0]=(dis_jk[2]*r_ij*r_jk-cosAng_ijk
                    *-dis_ij[2]*r_jk*r_jk)/(r_ij*r_ij*r_jk*r_jk);
                dcA_ijk[0][1]=(-dis_ij[0]*r_ij*r_jk-cosAng_ijk
                    *dis_jk[0]*r_ij*r_ij)/(r_ij*r_ij*r_jk*r_jk);
                dcA_ijk[1][1]=(-dis_ij[1]*r_ij*r_jk-cosAng_ijk
                    *dis_jk[1]*r_ij*r_ij)/(r_ij*r_ij*r_jk*r_jk);
                dcA_ijk[2][1]=(-dis_ij[2]*r_ij*r_jk-cosAng_ijk
                    *dis_jk[2]*r_ij*r_ij)/(r_ij*r_ij*r_jk*r_jk);
                nb_jk=nb_t;
                nb_t++;
                if(nb_t>nb_pi) {
                  new_n_tot=nb_pi+maxneigh;
                  grow_pi(nb_pi,new_n_tot);
                  nb_pi=new_n_tot;
                }
                bt_pi[nb_jk].i=j;
                bt_pi[nb_jk].j=k;
                bt_pi[nb_jk].temp=temp_jk;
                cosSq=cosAng_ijk*cosAng_ijk;
                sinFactor=.5*(1.0-cosSq)*pi_p[jtype-1]*betaS_jk;
                cosFactor=.5*(1.0+cosSq)*betaP_jk;
                betaCapSq1=pi_p[jtype-1]*betaS_jk*betaS_jk
                    -betaP_jk*betaP_jk;
                dbetaCapSq1=2.0*pi_p[jtype-1]*betaS_jk*dBetaS_jk
                    -2.0*betaP_jk*dBetaP_jk;

//AA is Eq. 37 (a) and Eq. 19 (b) for j atoms
//3rd BB is 2nd term of Eq. 38 (a) where i and k =neighbors j

                AA=AA+sinFactor*betaS_jk+cosFactor*betaP_jk;
                BB=BB+.25*(1.0-cosSq)*(1.0-cosSq)*betaCapSq1*betaCapSq1;

                agpdpr1=(2.0*sinFactor*dBetaS_jk+2.0*cosFactor
                    *dBetaP_jk)/r_jk;

//agpdpr1 is derivative of AA for atom j w.r.t. Beta(r_jk)
//agpdpr2 is derivative of BB for atom j w.r.t. Beta(r_jk)
//app1 is derivative of AA for j atom w.r.t. cos(theta_ijk)
//app2 is derivative of BB 2nd term w.r.t. cos(theta_ijk)

                agpdpr2=.5*(1.0-cosSq)*(1.0-cosSq)*betaCapSq1*dbetaCapSq1/r_jk;
                app1=cosAng_ijk*(-pi_p[jtype-1]*betaS_jk*betaS_jk
                    +betaP_jk*betaP_jk);
                app2=-(1.0-cosSq)*cosAng_ijk*betaCapSq1*betaCapSq1;
                bt_pi[nb_ij].dAA[0]-=
                    app1*dcA_ijk[0][0];
                bt_pi[nb_ij].dAA[1]-=
                    app1*dcA_ijk[1][0];
                bt_pi[nb_ij].dAA[2]-=
                    app1*dcA_ijk[2][0];
                bt_pi[nb_ij].dBB[0]-=
                    app2*dcA_ijk[0][0];
                bt_pi[nb_ij].dBB[1]-=
                    app2*dcA_ijk[1][0];
                bt_pi[nb_ij].dBB[2]-=
                    app2*dcA_ijk[2][0];
                bt_pi[nb_jk].dAA[0]+=
                    agpdpr1*dis_jk[0]
                    +app1*dcA_ijk[0][1];
                bt_pi[nb_jk].dAA[1]+=
                    agpdpr1*dis_jk[1]
                    +app1*dcA_ijk[1][1];
                bt_pi[nb_jk].dAA[2]+=
                    agpdpr1*dis_jk[2]
                    +app1*dcA_ijk[2][1];
                bt_pi[nb_jk].dBB[0]+=
                    app2*dcA_ijk[0][1]
                    +agpdpr2*dis_jk[0];
                bt_pi[nb_jk].dBB[1]+=
                    app2*dcA_ijk[1][1]
                    +agpdpr2*dis_jk[1];
                bt_pi[nb_jk].dBB[2]+=
                    app2*dcA_ijk[2][1]
                    +agpdpr2*dis_jk[2];

//j is a neighbor of i and k and k' are different neighbors of j not equal to i

                for(ltmp=0;ltmp<ktmp;ltmp++) {
                  if(ltmp!=ki) {
                    temp_jkp=BOP_index[j]+ltmp;
                    kp=jlist[ltmp];
                    kptype=map[type[kp]]+1;
                    for(nsearch=0;nsearch<nPiBk[n];nsearch++) {
                      ncmp=itypePiBk[n][nsearch];
                      if(x[ncmp][0]==x[kp][0]) {
                        if(x[ncmp][1]==x[kp][1]) {
                          if(x[ncmp][2]==x[kp][2]) {
                            break;
                          }
                        }
                      }
                    }
                    if(jtype==kptype)
                      ijkp=jtype-1;
                    else if(jtype<kptype)
                      ijkp=jtype*bop_types-jtype*(jtype+1)/2+kptype-1;
                    else
                      ijkp=kptype*bop_types-kptype*(kptype+1)/2+jtype-1;
                    dis_jkp[0]=x[kp][0]-x[j][0];
                    dis_jkp[1]=x[kp][1]-x[j][1];
                    dis_jkp[2]=x[kp][2]-x[j][2];
                    rsq_jkp=dis_jkp[0]*dis_jkp[0]
                        +dis_jkp[1]*dis_jkp[1]
                        +dis_jkp[2]*dis_jkp[2];
                    r_jkp=sqrt(rsq_jkp);
                    if(r_jkp<=rcut[ijkp]) {
                      ps=r_jkp*rdr[ijkp]+1.0;
                      ks=(int)ps;
                      if(nr-1<ks)
                        ks=nr-1;
                      ps=ps-ks;
                      if(ps>1.0)
                        ps=1.0;
                      betaS_jkp=((pBetaS3[ijkp][ks-1]*ps+pBetaS2[ijkp][ks-1])*ps
                          +pBetaS1[ijkp][ks-1])*ps+pBetaS[ijkp][ks-1];
                      dBetaS_jkp=(pBetaS6[ijkp][ks-1]*ps+pBetaS5[ijkp][ks-1])*ps
                          +pBetaS4[ijkp][ks-1];
                      betaP_jkp=((pBetaP3[ijkp][ks-1]*ps+pBetaP2[ijkp][ks-1])*ps
                          +pBetaP1[ijkp][ks-1])*ps+pBetaP[ijkp][ks-1];
                      dBetaP_jkp=(pBetaP6[ijkp][ks-1]*ps+pBetaP5[ijkp][ks-1])*ps
                          +pBetaP4[ijkp][ks-1];
                      cosAng_ijkp=(-dis_ij[0]*dis_jkp[0]-dis_ij[1]*dis_jkp[1]
                          -dis_ij[2]*dis_jkp[2])/(r_ij*r_jkp);
                      dcA_ijkp[0][0]=(dis_jkp[0]*r_ij*r_jkp-cosAng_ijkp
                          *-dis_ij[0]*r_jkp*r_jkp)/(r_ij*r_ij*r_jkp*r_jkp);
                      dcA_ijkp[1][0]=(dis_jkp[1]*r_ij*r_jkp-cosAng_ijkp
                          *-dis_ij[1]*r_jkp*r_jkp)/(r_ij*r_ij*r_jkp*r_jkp);
                      dcA_ijkp[2][0]=(dis_jkp[2]*r_ij*r_jkp-cosAng_ijkp
                          *-dis_ij[2]*r_jkp*r_jkp)/(r_ij*r_ij*r_jkp*r_jkp);
                      dcA_ijkp[0][1]=(-dis_ij[0]*r_ij*r_jkp-cosAng_ijkp
                          *dis_jkp[0]*r_ij*r_ij)/(r_ij*r_ij*r_jkp*r_jkp);
                      dcA_ijkp[1][1]=(-dis_ij[1]*r_ij*r_jkp-cosAng_ijkp
                          *dis_jkp[1]*r_ij*r_ij)/(r_ij*r_ij*r_jkp*r_jkp);
                      dcA_ijkp[2][1]=(-dis_ij[2]*r_ij*r_jkp-cosAng_ijkp
                          *dis_jkp[2]*r_ij*r_ij)/(r_ij*r_ij*r_jkp*r_jkp);
                      cosAng_kjkp=(dis_jk[0]*dis_jkp[0]+dis_jk[1]*dis_jkp[1]
                          +dis_jk[2]*dis_jkp[2])/(r_jk*r_jkp);
                      dcA_kjkp[0][0]=(dis_jkp[0]*r_jk*r_jkp-cosAng_kjkp
                          *dis_jk[0]*r_jkp*r_jkp)/(r_jk*r_jk*r_jkp*r_jkp);
                      dcA_kjkp[1][0]=(dis_jkp[1]*r_jk*r_jkp-cosAng_kjkp
                          *dis_jk[1]*r_jkp*r_jkp)/(r_jk*r_jk*r_jkp*r_jkp);
                      dcA_kjkp[2][0]=(dis_jkp[2]*r_jk*r_jkp-cosAng_kjkp
                          *dis_jk[2]*r_jkp*r_jkp)/(r_jk*r_jk*r_jkp*r_jkp);
                      dcA_kjkp[0][1]=(dis_jk[0]*r_jk*r_jkp-cosAng_kjkp
                          *dis_jkp[0]*r_jk*r_jk)/(r_jk*r_jk*r_jkp*r_jkp);
                      dcA_kjkp[1][1]=(dis_jk[1]*r_jk*r_jkp-cosAng_kjkp
                          *dis_jkp[1]*r_jk*r_jk)/(r_jk*r_jk*r_jkp*r_jkp);
                      dcA_kjkp[2][1]=(dis_jk[2]*r_jk*r_jkp-cosAng_kjkp
                          *dis_jkp[2]*r_jk*r_jk)/(r_jk*r_jk*r_jkp*r_jkp);
                      nb_jkp=nb_t;
                      nb_t++;
                      if(nb_t>nb_pi) {
                        new_n_tot=nb_pi+maxneigh;
                        grow_pi(nb_pi,new_n_tot);
                        nb_pi=new_n_tot;
                      }
                    bt_pi[nb_jkp].i=j;
                    bt_pi[nb_jkp].j=kp;
                    bt_pi[nb_jkp].temp=temp_jkp;
                    betaCapSq2=pi_p[jtype-1]*betaS_jkp*betaS_jkp
                      -betaP_jkp*betaP_jkp;
                    dbetaCapSq2=2.0*pi_p[jtype-1]*betaS_jkp*dBetaS_jkp
                      -2.0*betaP_jkp*dBetaP_jkp;
                    cosSq1=cosAng_ijkp*cosAng_ijkp;
                    angFactor=cosAng_kjkp-cosAng_ijkp*cosAng_ijk;
                    angFactor1=4.0*angFactor;
                    angFactor2=-angFactor1*cosAng_ijkp
                      +2.0*cosAng_ijk*(1.0-cosSq1);
                    angFactor3=-angFactor1*cosAng_ijk
                      +2.0*cosAng_ijkp*(1.0-cosSq);
                    angFactor4=2.0*angFactor*angFactor-(1.0-cosSq)*(1.0-cosSq1);
                      betaCapSum=.5*betaCapSq1*betaCapSq2;

//4th BB is 4th term of Eq. 38 (a) where i , k and k' =neighbors j

                    BB=BB+betaCapSum*angFactor4;

//app1 is derivative of BB 4th term w.r.t. cos(theta_kjk')
//app2 is derivative of BB 4th term w.r.t. cos(theta_ijk)
//app3 is derivative of BB 4th term w.r.t. cos(theta_ijk')
//agpdpr1 is derivative of BB 4th term for atom j w.r.t. Beta(r_jk)
//agpdpr2 is derivative of BB 4th term for atom j w.r.t. Beta(r_jk')

                    app1=betaCapSum*angFactor1;
                    app2=betaCapSum*angFactor2;
                    app3=betaCapSum*angFactor3;
                    agpdpr1=.5*angFactor4*dbetaCapSq1*betaCapSq2/r_jk;
                    agpdpr2=.5*angFactor4*betaCapSq1*dbetaCapSq2/r_jkp;
                    bt_pi[nb_ij].dBB[0]-=
                      app3*dcA_ijkp[0][0]
                      +app2*dcA_ijk[0][0];
                    bt_pi[nb_ij].dBB[1]-=
                      app3*dcA_ijkp[1][0]
                      +app2*dcA_ijk[1][0];
                    bt_pi[nb_ij].dBB[2]-=
                      app3*dcA_ijkp[2][0]
                      +app2*dcA_ijk[2][0];
                    bt_pi[nb_jk].dBB[0]+=
                      agpdpr1*dis_jk[0]
                      +app1*dcA_kjkp[0][0]
                      +app2*dcA_ijk[0][1];
                    bt_pi[nb_jk].dBB[1]+=
                      agpdpr1*dis_jk[1]
                      +app1*dcA_kjkp[1][0]
                      +app2*dcA_ijk[1][1];
                    bt_pi[nb_jk].dBB[2]+=
                      agpdpr1*dis_jk[2]
                      +app1*dcA_kjkp[2][0]
                      +app2*dcA_ijk[2][1];
                    bt_pi[nb_jkp].dBB[0]+=
                      agpdpr2*dis_jkp[0]
                      +app1*dcA_kjkp[0][1]
                      +app3*dcA_ijkp[0][1];
                    bt_pi[nb_jkp].dBB[1]+=
                      agpdpr2*dis_jkp[1]
                      +app1*dcA_kjkp[1][1]
                      +app3*dcA_ijkp[1][1];
                    bt_pi[nb_jkp].dBB[2]+=
                      agpdpr2*dis_jkp[2]
                      +app1*dcA_kjkp[2][1]
                      +app3*dcA_ijkp[2][1];
                    }
                  }
                }

//j and k' are different neighbors of i and k is a neighbor of j not equal to i

                for(ltmp=0;ltmp<numneigh[i];ltmp++) {
                  if(ltmp!=jtmp) {
                    temp_ikp=BOP_index[i]+ltmp;
                    kp=iilist[ltmp];
                    kptype=map[type[kp]]+1;
                    for(nsearch=0;nsearch<nPiBk[n];nsearch++) {
                      ncmp=itypePiBk[n][nsearch];
                      if(x[ncmp][0]==x[kp][0]) {
                        if(x[ncmp][1]==x[kp][1]) {
                          if(x[ncmp][2]==x[kp][2]) {
                            break;
                          }
                        }
                      }
                    }
                    if(itype==kptype)
                      iikp=itype-1;
                    else if(itype<kptype)
                      iikp=itype*bop_types-itype*(itype+1)/2+kptype-1;
                    else
                      iikp=kptype*bop_types-kptype*(kptype+1)/2+itype-1;
                    dis_ikp[0]=x[kp][0]-x[i][0];
                    dis_ikp[1]=x[kp][1]-x[i][1];
                    dis_ikp[2]=x[kp][2]-x[i][2];
                    rsq_ikp=dis_ikp[0]*dis_ikp[0]
                        +dis_ikp[1]*dis_ikp[1]
                        +dis_ikp[2]*dis_ikp[2];
                    r_ikp=sqrt(rsq_ikp);
                    if(r_ikp<=rcut[iikp]) {
                      ps=r_ikp*rdr[iikp]+1.0;
                      ks=(int)ps;
                      if(nr-1<ks)
                        ks=nr-1;
                      ps=ps-ks;
                      if(ps>1.0)
                        ps=1.0;
                      betaS_ikp=((pBetaS3[iikp][ks-1]*ps+pBetaS2[iikp][ks-1])*ps
                          +pBetaS1[iikp][ks-1])*ps+pBetaS[iikp][ks-1];
                      dBetaS_ikp=(pBetaS6[iikp][ks-1]*ps+pBetaS5[iikp][ks-1])*ps
                          +pBetaS4[iikp][ks-1];
                      betaP_ikp=((pBetaP3[iikp][ks-1]*ps+pBetaP2[iikp][ks-1])*ps
                          +pBetaP1[iikp][ks-1])*ps+pBetaP[iikp][ks-1];
                      dBetaP_ikp=(pBetaP6[iikp][ks-1]*ps+pBetaP5[iikp][ks-1])*ps
                          +pBetaP4[iikp][ks-1];
                      cosAng_jikp=(dis_ij[0]*dis_ikp[0]+dis_ij[1]*dis_ikp[1]
                          +dis_ij[2]*dis_ikp[2])/(r_ij*r_ikp);
                      dcA_jikp[0][0]=(dis_ikp[0]*r_ij*r_ikp-cosAng_jikp
                          *dis_ij[0]*r_ikp*r_ikp)/(r_ij*r_ij*r_ikp*r_ikp);
                      dcA_jikp[1][0]=(dis_ikp[1]*r_ij*r_ikp-cosAng_jikp
                          *dis_ij[1]*r_ikp*r_ikp)/(r_ij*r_ij*r_ikp*r_ikp);
                      dcA_jikp[2][0]=(dis_ikp[2]*r_ij*r_ikp-cosAng_jikp
                          *dis_ij[2]*r_ikp*r_ikp)/(r_ij*r_ij*r_ikp*r_ikp);
                      dcA_jikp[0][1]=(dis_ij[0]*r_ij*r_ikp-cosAng_jikp
                          *dis_ikp[0]*r_ij*r_ij)/(r_ij*r_ij*r_ikp*r_ikp);
                      dcA_jikp[1][1]=(dis_ij[1]*r_ij*r_ikp-cosAng_jikp
                          *dis_ikp[1]*r_ij*r_ij)/(r_ij*r_ij*r_ikp*r_ikp);
                      dcA_jikp[2][1]=(dis_ij[2]*r_ij*r_ikp-cosAng_jikp
                          *dis_ikp[2]*r_ij*r_ij)/(r_ij*r_ij*r_ikp*r_ikp);
                      nb_ikp=nb_t;
                      nb_t++;
                      if(nb_t>nb_pi) {
                        new_n_tot=nb_pi+maxneigh;
                        grow_pi(nb_pi,new_n_tot);
                        nb_pi=new_n_tot;
                      }
                      bt_pi[nb_ikp].i=i;
                      bt_pi[nb_ikp].j=kp;
                      bt_pi[nb_ikp].temp=temp_ikp;

                      betaCapSq2=pi_p[itype-1]*betaS_ikp*betaS_ikp
                          -betaP_ikp*betaP_ikp;
                      dbetaCapSq2=2.0*pi_p[itype-1]*betaS_ikp*dBetaS_ikp
                          -2.0*betaP_ikp*dBetaP_ikp;
                      dotV=(dis_jk[0]*dis_ikp[0]+dis_jk[1]
                          *dis_ikp[1]+dis_jk[2]*dis_ikp[2])
                          /(r_jk*r_ikp);
                      cosSq1=cosAng_jikp*cosAng_jikp;
                      angFactor=dotV+cosAng_jikp*cosAng_ijk;
                      angRfactor=4.0*angFactor*dotV;
                      dAngR1=-angRfactor/r_jk;
                      dAngR2=-angRfactor/r_ikp;
                      angFactor1=4.0*angFactor*cosAng_jikp
                          +2.0*cosAng_ijk*(1.0-cosSq1);
                      angFactor2=4.0*angFactor*cosAng_ijk
                          +2.0*cosAng_jikp*(1.0-cosSq);
                      angFactor3=2.0*angFactor*angFactor-(1.0-cosSq)*(1.0-cosSq1);
                      betaCapSum=.5*betaCapSq1*betaCapSq2;

//5th BB is 5th term of Eq. 38 (a) Eq. 21 (b) where i , k and k' =neighbors j

                      BB=BB+betaCapSum*angFactor3;

//app1 is derivative of BB 5th term w.r.t. cos(theta_ijk)
//app2 is derivative of BB 5th term w.r.t. cos(theta_jik')
//agpdpr1 is derivative of BB 5th term for atom j w.r.t. Beta(r_jk)
//agpdpr2 is derivative of BB 5th term for atom j w.r.t. Beta(r_ik')
//agpdpr3 is derivative of BB 5th term for atom j w.r.t. dot(r_ik',r_ij)

                      app1=betaCapSum*angFactor1;
                      app2=betaCapSum*angFactor2;
                      agpdpr1=(.5*angFactor3*dbetaCapSq1*betaCapSq2
                          +betaCapSum*dAngR1)/r_jk;
                      agpdpr2=(.5*angFactor3*betaCapSq1*dbetaCapSq2
                          +betaCapSum*dAngR2)/r_ikp;
                      agpdpr3=4.0*betaCapSum*angFactor/(r_ikp*r_jk);
                      bt_pi[nb_ij].dBB[0]+=
                          +app2*dcA_jikp[0][0]
                          -app1*dcA_ijk[0][0];
                      bt_pi[nb_ij].dBB[1]+=
                          +app2*dcA_jikp[1][0]
                          -app1*dcA_ijk[1][0];
                      bt_pi[nb_ij].dBB[2]+=
                          +app2*dcA_jikp[2][0]
                          -app1*dcA_ijk[2][0];
                      bt_pi[nb_ikp].dBB[0]+=
                          agpdpr2*dis_ikp[0]
                          +agpdpr3*dis_jk[0]
                          +app2*dcA_jikp[0][1];
                      bt_pi[nb_ikp].dBB[1]+=
                          agpdpr2*dis_ikp[1]
                          +agpdpr3*dis_jk[1]
                          +app2*dcA_jikp[1][1];
                      bt_pi[nb_ikp].dBB[2]+=
                          agpdpr2*dis_ikp[2]
                          +agpdpr3*dis_jk[2]
                          +app2*dcA_jikp[2][1];
                      bt_pi[nb_jk].dBB[0]+=
                          agpdpr1*dis_jk[0]
                          +agpdpr3*dis_ikp[0]
                          +app1*dcA_ijk[0][1];
                      bt_pi[nb_jk].dBB[1]+=
                          agpdpr1*dis_jk[1]
                          +agpdpr3*dis_ikp[1]
                          +app1*dcA_ijk[1][1];
                      bt_pi[nb_jk].dBB[2]+=
                          agpdpr1*dis_jk[2]
                          +agpdpr3*dis_ikp[2]
                          +app1*dcA_ijk[2][1];
                    }
                  }
                }
                if(pi_flag==0)
                  nPiBk[n]=nPiBk[n]+1;
              }
            }
          }
          CC=betaP_ij*betaP_ij+pi_delta[iij]*pi_delta[iij];
          BBrt=sqrt(BB+small6);
          AB1=CC+pi_c[iij]*(AA+BBrt)+small7;
          AB2=CC+pi_c[iij]*(AA-BBrt+sqrt(small6))+small7;
          BBrtR=1.0/BBrt;
          ABrtR1=1.0/sqrt(AB1);
          ABrtR2=1.0/sqrt(AB2);

// piB is similary formulation to (a) Eq. 36 and (b) Eq. 18

          piB[n]=(ABrtR1+ABrtR2)*pi_a[iij]*betaP_ij;
          dPiB1=-.5*(cube(ABrtR1)+cube(ABrtR2))*pi_c[iij]*pi_a[iij]*betaP_ij;
          dPiB2=.25*BBrtR*(cube(ABrtR2)-cube(ABrtR1))*pi_c[iij]*pi_a[iij]*betaP_ij;
          dPiB3=((ABrtR1+ABrtR2)*pi_a[iij]-(cube(ABrtR1)+cube(ABrtR2))*pi_a[iij]
              *betaP_ij*betaP_ij)*dBetaP_ij/r_ij;
          n++;

          pp2=2.0*betaP_ij;
          for(m=0;m<nb_t;m++) {
            bt_i=bt_pi[m].i;
            bt_j=bt_pi[m].j;
            xtmp[0]=x[bt_j][0]-x[bt_i][0];
            xtmp[1]=x[bt_j][1]-x[bt_i][1];
            xtmp[2]=x[bt_j][2]-x[bt_i][2];
            for(pp=0;pp<3;pp++) {
              bt_pi[m].dPiB[pp]=
                  +dPiB1*bt_pi[m].dAA[pp]
                  +dPiB2*bt_pi[m].dBB[pp];
              ftmp[pp]=pp2*bt_pi[m].dPiB[pp];
              f[bt_i][pp]-=ftmp[pp];
              f[bt_j][pp]+=ftmp[pp];
            }
            if(evflag) {
              ev_tally_xyz(bt_i,bt_j,nlocal,newton_pair,0.0,0.0,ftmp[0],ftmp[1]
                  ,ftmp[2],xtmp[0],xtmp[1],xtmp[2]);
            }
          }
          for(pp=0;pp<3;pp++) {
            ftmp[pp]=pp2*dPiB3*dis_ij[pp];
            f[i][pp]-=ftmp[pp];
            f[j][pp]+=ftmp[pp];
          }
          if(evflag) {
            ev_tally_xyz(i,j,nlocal,newton_pair,0.0,0.0,ftmp[0],ftmp[1]
                ,ftmp[2],dis_ij[0],dis_ij[1],dis_ij[2]);
          }
        }
      }
    }
  }
  destroy_pi();
}

/* ----------------------------------------------------------------------
   read BOP potential file
------------------------------------------------------------------------- */

void PairBOP::read_file(char *filename)
{
  int i,j,k;
  int ij,ii,jj;
  int buf1;
  int n;
  double buf2;
  char s[MAXLINE];
  char buf[2];

  MPI_Comm_rank(world,&me);

// read file on proc 0

  rcore=0.1;

  if (me == 0) {
    FILE *fp = force->open_potential(filename);
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open BOP potential file %s",filename);
      error->one(FLERR,str);
    }

// read parameters

    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&bop_types);
    fclose(fp);
    npairs=bop_types*(bop_types+1)/2;
  }
  MPI_Bcast(&bop_types,1,MPI_INT,0,world);
  MPI_Bcast(&npairs,1,MPI_INT,0,world);
  allocate();
  memory->create(pi_a,npairs,"BOP:pi_a");
  memory->create(pro_delta,bop_types,"BOP:pro_delta");
  memory->create(pi_delta,npairs,"BOP:pi_delta");
  memory->create(pi_p,bop_types,"BOP:pi_p");
  memory->create(pi_c,npairs,"BOP:pi_c");
  memory->create(sigma_r0,npairs,"BOP:sigma_r0");
  memory->create(pi_r0,npairs,"BOP:pi_r0");
  memory->create(phi_r0,npairs,"BOP:phi_r0");
  memory->create(sigma_rc,npairs,"BOP:sigma_rc");
  memory->create(pi_rc,npairs,"BOP:pi_rc");
  memory->create(phi_rc,npairs,"BOP:phi_rc");
  memory->create(r1,npairs,"BOP:r1");
  memory->create(sigma_beta0,npairs,"BOP:sigma_beta0");
  memory->create(pi_beta0,npairs,"BOP:pi_beta0");
  memory->create(phi0,npairs,"BOP:phi0");
  memory->create(sigma_n,npairs,"BOP:sigma_n");
  memory->create(pi_n,npairs,"BOP:pi_n");
  memory->create(phi_m,npairs,"BOP:phi_m");
  memory->create(sigma_nc,npairs,"BOP:sigma_nc");
  memory->create(pi_nc,npairs,"BOP:pi_nc");
  memory->create(phi_nc,npairs,"BOP:phi_nc");
  memory->create(pro,bop_types,"BOP:pro");
  memory->create(sigma_delta,npairs,"BOP:sigma_delta");
  memory->create(sigma_c,npairs,"BOP:sigma_c");
  memory->create(sigma_a,npairs,"BOP:sigma_a");
  memory->create(sigma_g0,bop_types
      ,bop_types,bop_types,"BOP:sigma_g0");
  memory->create(sigma_g1,bop_types
      ,bop_types,bop_types,"BOP:sigma_g1");
  memory->create(sigma_g2,bop_types
      ,bop_types,bop_types,"BOP:sigma_g2");
  memory->create(sigma_g3,bop_types
      ,bop_types,bop_types,"BOP:sigma_g3");
  memory->create(sigma_g4,bop_types
      ,bop_types,bop_types,"BOP:sigma_g4");
  memory->create(sigma_f,npairs,"BOP:sigma_f");
  memory->create(sigma_k,npairs,"BOP:sigma_k");
  memory->create(small3,npairs,"BOP:small3");

  if (me == 0) {
    words = new char*[bop_types];
    for(i=0;i<bop_types;i++) words[i]=NULL;
    FILE *fp = force->open_potential(filename);
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open BOP potential file %s",filename);
      error->one(FLERR,str);
    }
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    for(i=0;i<bop_types;i++) {
      fgets(s,MAXLINE,fp);
      sscanf(s,"%d %lf %s",&buf1,&buf2,buf);
      n= strlen(buf)+1;
      words[i] = new char[n];
      strcpy(words[i],buf);
    }
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lf%lf%lf%lf%lf%lf%lf",&small1,&small2,&small3g,&small4
        ,&small5,&small6,&small7);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d%lf%lf",&ncutoff,&rbig,&rsmall);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lf%lf%d",&which,&alpha,&nfunc);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lf%lf%lf",&alpha1,&beta1,&gamma1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lf%lf",&alpha2,&beta2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lf%lf",&alpha3,&beta3);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    for(i=0;i<bop_types;i++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lf%lf%lf",&pro[i],&pro_delta[i],&pi_p[i]);
    }
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    cutmax=0;

    for(i=0;i<bop_types;i++) {
      ii=i+1;
      for(j=i;j<bop_types;j++) {
        jj=j+1;
        if(ii==jj)
          ij=ii-1;
        else if(ii<jj)
          ij=ii*bop_types-ii*(ii+1)/2+jj-1;
        else
          ij=jj*bop_types-jj*(jj+1)/2+ii-1;
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lf%lf%lf%lf",&sigma_r0[ij],&sigma_rc[ij],&r1[ij],&rcut[ij]);
        if(rcut[ij]>cutmax)
          cutmax=rcut[ij];
        pi_r0[ij]=sigma_r0[ij];
        phi_r0[ij]=sigma_r0[ij];
        pi_rc[ij]=sigma_rc[ij];
        phi_rc[ij]=sigma_rc[ij];
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lf%lf%lf",&phi_m[ij],&sigma_n[ij],&sigma_nc[ij]);
        pi_n[ij]=sigma_n[ij];
        pi_nc[ij]=sigma_nc[ij];
        phi_nc[ij]=sigma_nc[ij];
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lf%lf%lf",&phi0[ij],&sigma_beta0[ij],&pi_beta0[ij]);
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lf%lf%lf",&sigma_a[ij],&sigma_c[ij],&sigma_delta[ij]);
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lf%lf%lf",&pi_a[ij],&pi_c[ij],&pi_delta[ij]);
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lf%lf%lf",&sigma_f[ij],&sigma_k[ij],&small3[ij]);
      }
    }
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    for(i=0;i<bop_types;i++) {
      for(j=0;j<bop_types;j++) {
        for(k=j;k<bop_types;k++) {
          fgets(s,MAXLINE,fp);
          sscanf(s,"%lf%lf%lf",&sigma_g0[j][i][k],&sigma_g1[j][i][k]
              ,&sigma_g2[j][i][k]);
          sigma_g0[k][i][j]=sigma_g0[j][i][k];
          sigma_g1[k][i][j]=sigma_g1[j][i][k];
          sigma_g2[k][i][j]=sigma_g2[j][i][k];
        }
      }
    }
    for(i=0;i<npairs;i++) {
      dr[i]=rcut[i]/(nr-1.0);
      rdr[i]=1.0/dr[i];
    }
    fclose(fp);
  }
  MPI_Bcast(&small1,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small2,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small3g,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small4,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small5,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small6,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small7,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&ncutoff,1,MPI_INT,0,world);
  MPI_Bcast(&rbig,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&rsmall,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&which,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&alpha,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&nfunc,1,MPI_INT,0,world);
  MPI_Bcast(&alpha1,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&beta1,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&gamma1,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&alpha2,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&beta2,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&alpha3,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&beta3,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&pro[0],bop_types,MPI_DOUBLE,0,world);
  MPI_Bcast(&pro_delta[0],bop_types,MPI_DOUBLE,0,world);
  MPI_Bcast(&pi_p[0],bop_types,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_r0[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_rc[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&r1[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&rcut[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&cutmax,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&pi_r0[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&phi_r0[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&pi_rc[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&phi_rc[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&phi_m[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_n[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_nc[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&pi_n[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&pi_nc[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&phi_nc[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&phi0[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_beta0[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&pi_beta0[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_a[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_c[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_delta[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&pi_a[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&pi_c[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&pi_delta[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_f[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_k[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&small3[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_g0[0][0][0],bop_types*bop_types*bop_types,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_g1[0][0][0],bop_types*bop_types*bop_types,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_g2[0][0][0],bop_types*bop_types*bop_types,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_g3[0][0][0],bop_types*bop_types*bop_types,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_g4[0][0][0],bop_types*bop_types*bop_types,MPI_DOUBLE,0,world);
  MPI_Bcast(&dr[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&rdr[0],npairs,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

void PairBOP::read_table(char *filename)
{
  int i,j,k,n;
  int buf1;
  double buf2;
  char s[MAXLINE],buf[2];

  MPI_Comm_rank(world,&me);

  if (me == 0) {
    FILE *fp = force->open_potential(filename);
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open BOP potential file %s",filename);
      error->one(FLERR,str);
    }
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&bop_types);
    words = new char*[bop_types];
    for(i=0;i<bop_types;i++) words[i]=NULL;
    for(i=0;i<bop_types;i++) {
      fgets(s,MAXLINE,fp);
      sscanf(s,"%d %lf %s",&buf1,&buf2,buf);
      n= strlen(buf)+1;
      words[i] = new char[n];
      strcpy(words[i],buf);
    }
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d %d",&nr,&nBOt);
    fclose(fp);
    npairs=bop_types*(bop_types+1)/2;
  }

  MPI_Bcast(&nr,1,MPI_INT,0,world);
  MPI_Bcast(&nBOt,1,MPI_INT,0,world);
  MPI_Bcast(&bop_types,1,MPI_INT,0,world);
  MPI_Bcast(&npairs,1,MPI_INT,0,world);
  memory->create(pi_a,npairs,"BOP:pi_a");
  memory->create(pro_delta,bop_types,"BOP:pro_delta");
  memory->create(pi_delta,npairs,"BOP:pi_delta");
  memory->create(pi_p,bop_types,"BOP:pi_p");
  memory->create(pi_c,npairs,"BOP:pi_c");
  memory->create(r1,npairs,"BOP:r1");
  memory->create(pro,bop_types,"BOP:pro");
  memory->create(sigma_delta,npairs,"BOP:sigma_delta");
  memory->create(sigma_c,npairs,"BOP:sigma_c");
  memory->create(sigma_a,npairs,"BOP:sigma_a");
  memory->create(sigma_g0,bop_types
      ,bop_types,bop_types,"BOP:sigma_g0");
  memory->create(sigma_g1,bop_types
      ,bop_types,bop_types,"BOP:sigma_g1");
  memory->create(sigma_g2,bop_types
      ,bop_types,bop_types,"BOP:sigma_g2");
  memory->create(sigma_f,npairs,"BOP:sigma_f");
  memory->create(sigma_k,npairs,"BOP:sigma_k");
  memory->create(small3,npairs,"BOP:small3");
  allocate();

  if (me == 0) {
    FILE *fp = force->open_potential(filename);
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open BOP potential file %s",filename);
      error->one(FLERR,str);
    }
    for(i=0;i<bop_types+2;i++) {
      fgets(s,MAXLINE,fp);
    }
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lf%lf%lf%lf%lf%lf%lf",&small1,&small2,&small3g
        ,&small4,&small5,&small6,&small7);
    for(i=0;i<bop_types;i++) {
      fgets(s,MAXLINE,fp);
      sscanf(s,"%lf",&pi_p[i]);
    }
    cutmax=0;
    for(i=0;i<npairs;i++) {
      fgets(s,MAXLINE,fp);
      sscanf(s,"%lf",&rcut[i]);
      if(rcut[i]>cutmax)
        cutmax=rcut[i];
      fgets(s,MAXLINE,fp);
      sscanf(s,"%lf%lf%lf%lf",&sigma_c[i],&sigma_a[i],&pi_c[i],&pi_a[i]);
      fgets(s,MAXLINE,fp);
      sscanf(s,"%lf%lf",&sigma_delta[i],&pi_delta[i]);
      fgets(s,MAXLINE,fp);
      sscanf(s,"%lf%lf%lf",&sigma_f[i],&sigma_k[i],&small3[i]);
    }
    for(i=0;i<bop_types;i++)
      for(j=0;j<bop_types;j++)
        for(k=0;k<bop_types;k++) {
          fgets(s,MAXLINE,fp);
          sscanf(s,"%lf%lf%lf",&sigma_g0[i][j][k],&sigma_g1[i][j][k],&sigma_g2[i][j][k]);
        }
    for(i=0;i<npairs;i++) {
      for(j=0;j<nr;j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lf%lf%lf%lf%lf",&pRepul[i][j],&pRepul[i][j+1]
            ,&pRepul[i][j+2],&pRepul[i][j+3],&pRepul[i][j+4]);
        j+=4;
      }
    }
    for(i=0;i<npairs;i++) {
      for(j=0;j<nr;j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lf%lf%lf%lf%lf",&pBetaS[i][j],&pBetaS[i][j+1]
            ,&pBetaS[i][j+2],&pBetaS[i][j+3],&pBetaS[i][j+4]);
        j+=4;
      }
    }
    for(i=0;i<npairs;i++) {
      for(j=0;j<nr;j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lf%lf%lf%lf%lf",&pBetaP[i][j],&pBetaP[i][j+1]
            ,&pBetaP[i][j+2],&pBetaP[i][j+3],&pBetaP[i][j+4]);
        j+=4;
      }
    }
    for(i=0;i<npairs;i++) {
      for(j=0;j<nBOt;j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lf%lf%lf%lf%lf",&FsigBO[i][j],&FsigBO[i][j+1]
            ,&FsigBO[i][j+2],&FsigBO[i][j+3],&FsigBO[i][j+4]);
        j+=4;
      }
    }
    for(i=0;i<bop_types;i++) {
      fgets(s,MAXLINE,fp);
      sscanf(s,"%lf",&pro_delta[i]);
    }
    for(i=0;i<bop_types;i++) {
      fgets(s,MAXLINE,fp);
      sscanf(s,"%lf",&pro[i]);
    }
    for(i=0;i<npairs;i++) {
      dr[i]=rcut[i]/((double)nr-1.0);
      rdr[i]=1.0/dr[i];
    }
    dBO=1.0/((double)nBOt-1.0);
    rdBO=1.0/(double)dBO;

    for(i=0;i<npairs;i++) {
      pBetaS1[i][0]=pBetaS[i][1]-pBetaS[i][0];
      pBetaS1[i][1]=0.5*(pBetaS[i][2]-pBetaS[i][0]);
      pBetaS1[i][nr-2]=0.5*(pBetaS[i][nr-1]-pBetaS[i][nr-3]);
      pBetaS1[i][nr-1]=pBetaS[i][nr-1]-pBetaS[i][nr-2];
      pBetaP1[i][0]=pBetaP[i][1]-pBetaP[i][0];
      pBetaP1[i][1]=0.5*(pBetaP[i][2]-pBetaP[i][0]);
      pBetaP1[i][nr-2]=0.5*(pBetaP[i][nr-1]-pBetaP[i][nr-3]);
      pBetaP1[i][nr-1]=pBetaP[i][nr-1]-pBetaP[i][nr-2];
      pRepul1[i][0]=pRepul[i][1]-pRepul[i][0];
      pRepul1[i][1]=0.5*(pRepul[i][2]-pRepul[i][0]);
      pRepul1[i][nr-2]=0.5*(pRepul[i][nr-1]-pRepul[i][nr-3]);
      pRepul1[i][nr-1]=pRepul[i][nr-1]-pRepul[i][nr-2];
      FsigBO1[i][0]=FsigBO[i][1]-FsigBO[i][0];
      FsigBO1[i][1]=0.5*(FsigBO[i][2]-FsigBO[i][0]);
      FsigBO1[i][nBOt-2]=0.5*(FsigBO[i][nBOt-1]-FsigBO[i][nBOt-3]);
      FsigBO1[i][nBOt-1]=FsigBO[i][nBOt-1]-FsigBO[i][nBOt-2];
      for(k=2;k<nr-2;k++) {
        pBetaS1[i][k]=((pBetaS[i][k-2]-pBetaS[i][k+2])
            +8.0*(pBetaS[i][k+1]-pBetaS[i][k-1]))/12.0;
        pBetaP1[i][k]=((pBetaP[i][k-2]-pBetaP[i][k+2])
            +8.0*(pBetaP[i][k+1]-pBetaP[i][k-1]))/12.0;
        pRepul1[i][k]=((pRepul[i][k-2]-pRepul[i][k+2])
            +8.0*(pRepul[i][k+1]-pRepul[i][k-1]))/12.0;
      }
      for(k=2;k<nr-2;k++) {
        FsigBO1[i][k]=((FsigBO[i][k-2]-FsigBO[i][k+2])
            +8.0*(FsigBO[i][k+1]-FsigBO[i][k-1]))/12.0;
      }
      for(k=0;k<nr-1;k++) {
        pBetaS2[i][k]=3.0*(pBetaS[i][k+1]-pBetaS[i][k])
            -2.0*pBetaS1[i][k]-pBetaS1[i][k+1];
        pBetaS3[i][k]=pBetaS1[i][k]+pBetaS1[i][k+1]
            -2.0*(pBetaS[i][k+1]-pBetaS[i][k]);
        pBetaP2[i][k]=3.0*(pBetaP[i][k+1]-pBetaP[i][k])
            -2.0*pBetaP1[i][k]-pBetaP1[i][k+1];
        pBetaP3[i][k]=pBetaP1[i][k]+pBetaP1[i][k+1]
            -2.0*(pBetaP[i][k+1]-pBetaP[i][k]);
        pRepul2[i][k]=3.0*(pRepul[i][k+1]-pRepul[i][k])
            -2.0*pRepul1[i][k]-pRepul1[i][k+1];
        pRepul3[i][k]=pRepul1[i][k]+pRepul1[i][k+1]
            -2.0*(pRepul[i][k+1]-pRepul[i][k]);
      }
      for(k=0;k<nBOt-1;k++) {
        FsigBO2[i][k]=3.0*(FsigBO[i][k+1]-FsigBO[i][k])
            -2.0*FsigBO1[i][k]-FsigBO1[i][k+1];
        FsigBO3[i][k]=FsigBO1[i][k]+FsigBO1[i][k+1]
            -2.0*(FsigBO[i][k+1]-FsigBO[i][k]);
      }
      pBetaS2[i][nr-1]=0.0;
      pBetaS3[i][nr-1]=0.0;
      pBetaP2[i][nr-1]=0.0;
      pBetaP3[i][nr-1]=0.0;
      pRepul2[i][nr-1]=0.0;
      pRepul3[i][nr-1]=0.0;
      FsigBO2[i][nBOt-1]=0.0;
      FsigBO3[i][nBOt-1]=0.0;
      for(k=0;k<nr;k++) {
        pBetaS4[i][k]=pBetaS1[i][k]/dr[i];
        pBetaS5[i][k]=2.0*pBetaS2[i][k]/dr[i];
        pBetaS6[i][k]=3.0*pBetaS3[i][k]/dr[i];
        pBetaP4[i][k]=pBetaP1[i][k]/dr[i];
        pBetaP5[i][k]=2.0*pBetaP2[i][k]/dr[i];
        pBetaP6[i][k]=3.0*pBetaP3[i][k]/dr[i];
        pRepul4[i][k]=pRepul1[i][k]/dr[i];
        pRepul5[i][k]=2.0*pRepul2[i][k]/dr[i];
        pRepul6[i][k]=3.0*pRepul3[i][k]/dr[i];
      }
      for(k=0;k<nBOt;k++) {
        FsigBO4[i][k]=FsigBO1[i][k]/dBO;
        FsigBO5[i][k]=2.0*FsigBO2[i][k]/dBO;
        FsigBO6[i][k]=3.0*FsigBO3[i][k]/dBO;
      }
    }
    fclose(fp);
  }
  MPI_Bcast(&rdBO,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&dBO,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&bop_types,1,MPI_INT,0,world);
  MPI_Bcast(&small1,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small2,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small3g,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small4,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small5,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small6,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&small7,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&pro[0],bop_types,MPI_DOUBLE,0,world);
  MPI_Bcast(&pro_delta[0],bop_types,MPI_DOUBLE,0,world);
  MPI_Bcast(&pi_p[0],bop_types,MPI_DOUBLE,0,world);
  MPI_Bcast(&r1[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&rcut[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&cutmax,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_a[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_c[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_delta[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&pi_a[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&pi_c[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&pi_delta[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_f[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_k[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&small3[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_g0[0][0][0],bop_types*bop_types*bop_types,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_g1[0][0][0],bop_types*bop_types*bop_types,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma_g2[0][0][0],bop_types*bop_types*bop_types,MPI_DOUBLE,0,world);
  MPI_Bcast(&dr[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&rdr[0],npairs,MPI_DOUBLE,0,world);
  MPI_Bcast(&pBetaS[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pBetaS1[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pBetaS2[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pBetaS3[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pBetaS4[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pBetaS5[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pBetaS6[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pBetaP[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pBetaP1[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pBetaP2[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pBetaP3[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pBetaP4[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pBetaP5[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pBetaP6[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pRepul[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pRepul1[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pRepul2[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pRepul3[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pRepul4[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pRepul5[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&pRepul6[0][0],npairs*nr,MPI_DOUBLE,0,world);
  MPI_Bcast(&FsigBO[0][0],npairs*nBOt,MPI_DOUBLE,0,world);
  MPI_Bcast(&FsigBO1[0][0],npairs*nBOt,MPI_DOUBLE,0,world);
  MPI_Bcast(&FsigBO2[0][0],npairs*nBOt,MPI_DOUBLE,0,world);
  MPI_Bcast(&FsigBO3[0][0],npairs*nBOt,MPI_DOUBLE,0,world);
  MPI_Bcast(&FsigBO4[0][0],npairs*nBOt,MPI_DOUBLE,0,world);
  MPI_Bcast(&FsigBO5[0][0],npairs*nBOt,MPI_DOUBLE,0,world);
  MPI_Bcast(&FsigBO6[0][0],npairs*nBOt,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

void PairBOP::setPbetaS()
{
  int i,j,k;
  double r,value,dvalue;

  for(i=0;i<npairs;i++) {
    for(j=0;j<nr;j++) {
      r=(double)j*dr[i];
      if(r<rcore)
        r=rcore;
      if(ncutoff==3) {
        if(r>=rcut[i])
          pBetaS[i][j]=0.0;
        else if(r<=r1[i]) {
          value=betaSfunc(i,r);
          dvalue=dBetaSfunc(i,r,value,1.0);
          pBetaS[i][j]=value;
        }
        else {
          value=betaSfunc(i,r1[i]);
          dvalue=dBetaSfunc(i,r1[i],value,1.0);
          pBetaS[i][j]=-(r-rcut[i])*(r-rcut[i])*(value*(2.0*r-3.0*r1[i]+rcut[i])
              -dvalue*(r-r1[i])*(r1[i]-rcut[i]))/((r1[i]-rcut[i])
              *(r1[i]-rcut[i])*(r1[i]-rcut[i]));
        }
      }
      else {
        if(r>=rcut[i])
          pBetaS[i][j]=0.0;
        else {
          value=betaSfunc(i,r);
          dvalue=dBetaSfunc(i,r,value,0.0);
          pBetaS[i][j]=value*cutoff(r1[i],rcut[i],ncutoff,r);
        }
      }
    }
    pBetaS[i][nr-1]=0.0;
    pBetaS1[i][0]=pBetaS[i][1]-pBetaS[i][0];
    pBetaS1[i][1]=0.5*(pBetaS[i][2]-pBetaS[i][0]);
    pBetaS1[i][nr-2]=0.5*(pBetaS[i][nr-1]-pBetaS[i][nr-3]);
    pBetaS1[i][nr-1]=pBetaS[i][nr-1]-pBetaS[i][nr-2];

    for(k=2;k<nr-2;k++) {
      pBetaS1[i][k]=((pBetaS[i][k-2]-pBetaS[i][k+2])+8.0*(pBetaS[i][k+1]
          -pBetaS[i][k-1]))/12.0;
    }
    for(k=0;k<nr-1;k++) {
      pBetaS2[i][k]=3.0*(pBetaS[i][k+1]-pBetaS[i][k])-2.0*pBetaS1[i][k]-pBetaS1[i][k+1];
      pBetaS3[i][k]=pBetaS1[i][k]+pBetaS1[i][k+1]-2.0*(pBetaS[i][k+1]-pBetaS[i][k]);
    }
    pBetaS2[i][nr-1]=0.0;
    pBetaS3[i][nr-1]=0.0;
    for(k=0;k<nr;k++) {
      pBetaS4[i][k]=pBetaS1[i][k]/dr[i];
      pBetaS5[i][k]=2.0*pBetaS2[i][k]/dr[i];
      pBetaS6[i][k]=3.0*pBetaS3[i][k]/dr[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairBOP::setPbetaP()
{
  int i,j,k;
  double r,value,dvalue;

  for(i=0;i<npairs;i++) {
    for(j=0;j<nr;j++) {
      r=(double)j*dr[i];
      if(r<rcore)
        r=rcore;
      if(ncutoff==3) {
        if(r>=rcut[i])
          pBetaP[i][j]=0.0;
        else if(r<=r1[i]) {
          value=betaPfunc(i,r);
          dvalue=dBetaPfunc(i,r,value,0.0);
          pBetaP[i][j]=value;
        }
        else {
          value=betaPfunc(i,r1[i]);
          dvalue=dBetaPfunc(i,r1[i],value,1.0);
          pBetaP[i][j]=-(r-rcut[i])*(r-rcut[i])*(value*(2.0*r-3.0*r1[i]
              +rcut[i])-dvalue*(r-r1[1])*(r1[i]-rcut[i]))/((r1[i]-rcut[i])
              *(r1[i]-rcut[i])*(r1[i]-rcut[i]));
        }
      }
      else {
        if(r>=rcut[i])
          pBetaP[i][j]=0.0;
        else {
          value=betaPfunc(i,r);
          dvalue=dBetaPfunc(i,r,value,0.0);
          pBetaP[i][j]=value*cutoff(r1[i],rcut[i],ncutoff,r);
        }
      }
    }
    pBetaP[i][nr-1]=0.0;
    pBetaP1[i][0]=pBetaP[i][1]-pBetaP[i][0];
    pBetaP1[i][1]=0.5*(pBetaP[i][2]-pBetaP[i][0]);
    pBetaP1[i][nr-2]=0.5*(pBetaP[i][nr-1]-pBetaP[i][nr-3]);
    pBetaP1[i][nr-1]=pBetaP[i][nr-1]-pBetaP[i][nr-2];
    for(k=2;k<nr-2;k++)
      pBetaP1[i][k]=((pBetaP[i][k-2]-pBetaP[i][k+2])+8.0*(pBetaP[i][k+1]
          -pBetaP[i][k-1]))/12.0;
    for(k=0;k<nr-1;k++) {
      pBetaP2[i][k]=3.0*(pBetaP[i][k+1]-pBetaP[i][k])-2.0*pBetaP1[i][k]-pBetaP1[i][k+1];
      pBetaP3[i][k]=pBetaP1[i][k]+pBetaP1[i][k+1]-2.0*(pBetaP[i][k+1]-pBetaP[i][k]);
    }
    pBetaP2[i][nr-1]=0.0;
    pBetaP3[i][nr-1]=0.0;
    for(k=0;k<nr;k++) {
      pBetaP4[i][k]=pBetaP1[i][k]/dr[i];
      pBetaP5[i][k]=2.0*pBetaP2[i][k]/dr[i];
      pBetaP6[i][k]=3.0*pBetaP3[i][k]/dr[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairBOP::setPrepul()
{
  int i,j,k;
  double r,value,dvalue;

  for(i=0;i<npairs;i++) {
    for(j=0;j<nr;j++) {
      r=(double)j*dr[i];
      if(r<rcore)
        r=rcore;
      if(ncutoff==3) {
        if(r>=rcut[i])
          pRepul[i][j]=0.0;
        else if(r<=r1[i]) {
          value=repulfunc(i,r);
          dvalue=dRepulfunc(i,r,value,0.0);
          pRepul[i][j]=value;
        }
        else {
          value=repulfunc(i,r1[i]);
          dvalue=dRepulfunc(i,r1[i],value,1.0);
          pRepul[i][j]=-(r-rcut[i])*(r-rcut[i])*(value*(2.0*r-3.0*r1[i]+rcut[i])
              -dvalue*(r-r1[i])*(r1[i]-rcut[i]))/((r1[i]-rcut[i])
              *(r1[i]-rcut[i])*(r1[i]-rcut[i]));
        }
      }
      else {
        if(r>=rcut[i])
          pRepul[i][j]=0.0;
        else {
          value=repulfunc(i,r);
          dvalue=dRepulfunc(i,r,value,0.0);
          pRepul[i][j]=value*cutoff(r1[i],rcut[i],ncutoff,r);
        }
      }
    }
    pRepul[i][nr-1]=0.0;
    pRepul1[i][0]=pRepul[i][1]-pRepul[i][0];
    pRepul1[i][1]=0.5*(pRepul[i][2]-pRepul[i][0]);
    pRepul1[i][nr-2]=0.5*(pRepul[i][nr-1]-pRepul[i][nr-3]);
    pRepul1[i][nr-1]=pRepul[i][nr-1]-pRepul[i][nr-2];
    for(k=2;k<nr-2;k++)
      pRepul1[i][k]=((pRepul[i][k-2]-pRepul[i][k+2])+8.0*(pRepul[i][k+1]
          -pRepul[i][k-1]))/12.0;
    for(k=0;k<nr-1;k++) {
      pRepul2[i][k]=3.0*(pRepul[i][k+1]-pRepul[i][k])-2.0*pRepul1[i][k]-pRepul1[i][k+1];
      pRepul3[i][k]=pRepul1[i][k]+pRepul1[i][k+1]-2.0*(pRepul[i][k+1]-pRepul[i][k]);
    }
    pRepul2[i][nr-1]=0.0;
    pRepul3[i][nr-1]=0.0;
    for(k=0;k<nr;k++) {
      pRepul4[i][k]=pRepul1[i][k]/dr[i];
      pRepul5[i][k]=2.0*pRepul2[i][k]/dr[i];
      pRepul6[i][k]=3.0*pRepul3[i][k]/dr[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairBOP::betaSfunc(int i,double r)
{
  double temp_value;

  if(nfunc==1) {
    temp_value=pow(sigma_r0[i]/r,sigma_n[i])*exp(sigma_n[i]*pow(sigma_r0[i]
        /sigma_rc[i],sigma_nc[i])-sigma_n[i]*pow(r/sigma_rc[i],sigma_nc[i]));
    temp_value=sigma_beta0[i]*temp_value;
  }
  if(nfunc==2)
    temp_value=sigma_beta0[i]*exp(-sigma_n[i]*r);
  if(nfunc==3)
    temp_value=sigma_beta0[i]/pow(r,sigma_n[i]);
  return(temp_value);
}

/* ---------------------------------------------------------------------- */

double PairBOP::dBetaSfunc(int i,double r,double value,double dmore)
{
  double temp_dvalue;

  if(nfunc==1)
    if(dmore==1.0)
      temp_dvalue=-sigma_n[i]*value/r*(1.0+sigma_nc[i]
          *pow(r/sigma_rc[i],sigma_nc[i]));
  if(nfunc==2)
    if(dmore==1.0)
      temp_dvalue=-sigma_n[i]*value;
  if(nfunc==3)
    if(dmore==1.0)
      temp_dvalue=-sigma_n[i]*value/r;
  return(temp_dvalue);
}

/* ---------------------------------------------------------------------- */

double PairBOP::betaPfunc(int i,double r)
{
  double temp_value;

  if(nfunc==1) {
    temp_value=pow(pi_r0[i]/r,pi_n[i])*exp(pi_n[i]*pow(pi_r0[i]
        /pi_rc[i],pi_nc[i])-pi_n[i]*pow(r/pi_rc[i],pi_nc[i]));
        temp_value=pi_beta0[i]*temp_value;
  }
  if(nfunc==2)
    temp_value=pi_beta0[i]*exp(-pi_n[i]*r);
  if(nfunc==3)
    temp_value=pi_beta0[i]/pow(r,pi_n[i]);
  return(temp_value);
}

/* ---------------------------------------------------------------------- */

double PairBOP::dBetaPfunc(int i,double r,double value,double dmore)
{
  double temp_dvalue;

  if(nfunc==1)
    if(dmore==1.0)
      temp_dvalue=-pi_n[i]*value/r*(1.0+pi_nc[i]*pow(r/pi_rc[i],pi_nc[i]));
  if(nfunc==2)
    if(dmore==1.0)
      temp_dvalue=-pi_n[i]*value;
  if(nfunc==3)
    if(dmore==1.0)
      temp_dvalue=-pi_n[i]*value/r;
  return(temp_dvalue);
}

/* ---------------------------------------------------------------------- */

double PairBOP::repulfunc(int i,double r)
{
  double temp_value;

  if(nfunc==1) {
    temp_value=pow(phi_r0[i]/r,phi_m[i])*exp(phi_m[i]*pow(phi_r0[i]/phi_rc[i]
        ,phi_nc[i])-phi_m[i]*pow(r/phi_rc[i],phi_nc[i]));
    temp_value=phi0[i]*temp_value;
  }
  if(nfunc==2)
    temp_value=phi0[i]*exp(-phi_m[i]*r);
  if(nfunc==3)
    temp_value=phi0[i]/pow(r,phi_m[i]);
  return(temp_value);
}

/* ---------------------------------------------------------------------- */

double PairBOP::dRepulfunc(int i,double r,double value,double dmore)
{
  double temp_dvalue;

  if(nfunc==1)
    if(dmore==1.0)
      temp_dvalue=-phi_m[i]*value/r*(1.0+phi_nc[i]*pow(r/phi_rc[i],phi_nc[i]));
  if(nfunc==2)
    if(dmore==1.0)
      temp_dvalue=-phi_m[i]*value;
  if(nfunc==3)
    if(dmore==1.0)
      temp_dvalue=-phi_m[i]*value/r;
  return(temp_dvalue);
}

/* ---------------------------------------------------------------------- */

void PairBOP::setSign()
{
  int i,j,k;
  double y0,tmp,xBO,fth,cs,bigF;
  double epsilon,fsigma1,slope,sat;

  dBO=1.0/(nBOt-1.0);
  rdBO=1.0/dBO;
  for(i=0;i<npairs;i++) {
    for(j=0;j<nBOt;j++) {
      xBO=(double)j*dBO;
      if(which==1.0) {
        fth=0.0;
        if(xBO>alpha)
          fth=4.0/3.0*(xBO-alpha);
        if(sigma_f[i]<=fth)
          FsigBO[i][j]=2.0*sigma_f[i];
        else if(sigma_f[i]>=1.0-fth)
          FsigBO[i][j]=2.0*(1.0-sigma_f[i]);
        else {
          cs=0.0;
          if(xBO<alpha)
            cs=32.0*(alpha-xBO);
          bigF=(sigma_f[i]*(1.0-sigma_f[i])-fth*(1.0-fth))/square(1.0-2.0*fth);
          FsigBO[i][j]=2.0*fth+2.0*bigF*(1.0-2.0*fth)*(1.0+bigF*(1.0-cs*bigF));
        }
      }
      else if(which==2.0) {
        epsilon=0.0000000001;
        fsigma1=sigma_f[i];
        if(fsigma1>0.5)
          fsigma1=1.0-fsigma1;
        y0=alpha1*pow(fsigma1,beta1)*pow(0.5-fsigma1,gamma1);
        slope=(1.0-exp(-alpha2*pow(fsigma1,beta2)))/(1.0-exp(-alpha2*pow(0.5,beta2)));
        sat=alpha3*fsigma1+beta3;
        tmp=y0+slope*xBO+sat;
        FsigBO[i][j]=(tmp-sqrt(tmp*tmp-4.0*(-epsilon*sqrt(1.0+slope*slope)
            +y0*sat+slope*sat*xBO)))/2.0;
      }
    }
    FsigBO1[i][0]=FsigBO[i][1]-FsigBO[i][0];
    FsigBO1[i][1]=0.5*(FsigBO[i][2]-FsigBO[i][0]);
    FsigBO1[i][nBOt-2]=0.5*(FsigBO[i][nBOt-1]-FsigBO[i][nBOt-3]);
    FsigBO1[i][nBOt-1]=FsigBO[i][nBOt-1]-FsigBO[i][nBOt-2];
    for(k=2;k<nBOt-2;k++)
      FsigBO1[i][k]=((FsigBO[i][k-2]-FsigBO[i][k+2])+8.0*(FsigBO[i][k+1]
          -FsigBO[i][k-1]))/12.0;
    for(k=0;k<nBOt-1;k++) {
      FsigBO2[i][k]=3.0*(FsigBO[i][k+1]-FsigBO[i][k])-2.0*FsigBO1[i][k]-FsigBO1[i][k+1];
      FsigBO3[i][k]=FsigBO1[i][k]+FsigBO1[i][k+1]-2.0*(FsigBO[i][k+1]-FsigBO[i][k]);
    }
    FsigBO2[i][nBOt-1]=0.0;
    FsigBO3[i][nBOt-1]=0.0;
    for(k=0;k<nBOt;k++) {
      FsigBO4[i][k]=FsigBO1[i][k]/dBO;
      FsigBO5[i][k]=2.0*FsigBO2[i][k]/dBO;
      FsigBO6[i][k]=3.0*FsigBO3[i][k]/dBO;
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairBOP::cutoff(double rp,double vrcut,int mode,double r)
{
  double tmp,tmp_beta,tmp_alpha,cut_store;

  if(mode==1) {
    tmp=(rsmall-rbig)*(r-rp)/(vrcut-rp)+rbig;
    cut_store=(erfc(tmp)-erfc(rsmall))/(erfc(rbig)-erfc(rsmall));
  }
  else {
    tmp_beta=log(log(rbig)/log(rsmall))/log(rp/vrcut);
    tmp_alpha=-log(rbig)/pow(rp,tmp_beta);
    cut_store=(exp(-tmp_alpha*pow(r,tmp_beta))-exp(-tmp_alpha*pow(vrcut
        ,tmp_beta)))/(exp(-tmp_alpha*pow(rp,tmp_beta))-exp(-tmp_alpha
        *pow(vrcut,tmp_beta)));
  }
  return(cut_store);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairBOP::memory_usage()
{
  int nlocal,nghost,nall;
  int n = atom->ntypes;
  nlocal = atom->nlocal;
  nghost = atom->nghost;
  nall = nlocal + nghost;
  double bytes = 0.0;

// rcut
  bytes += npairs * sizeof (double);
// dr
  bytes += npairs * sizeof (double);
// rdr
  bytes += npairs * sizeof (double);
// setflag
  bytes += (n+1) * (n+1) * sizeof (int);
// cutsq
  bytes += (n+1) * (n+1) * sizeof (double);
// cutghost
  bytes += (n+1) * (n+1) * sizeof (double);
// cutghost
  bytes += (n+1) * (n+1) * sizeof (double);
// pBetaS
  bytes += npairs * nr * sizeof (double);
// pBetaS1
  bytes += npairs * nr * sizeof (double);
// pBetaS2
  bytes += npairs * nr * sizeof (double);
// pBetaS3
  bytes += npairs * nr * sizeof (double);
// pBetaS4
  bytes += npairs * nr * sizeof (double);
// pBetaS5
  bytes += npairs * nr * sizeof (double);
// pBetaS6
  bytes += npairs * nr * sizeof (double);
// pBetaP
  bytes += npairs * nr * sizeof (double);
// pBetaP1
  bytes += npairs * nr * sizeof (double);
// pBetaP2
  bytes += npairs * nr * sizeof (double);
// pBetaP3
  bytes += npairs * nr * sizeof (double);
// pBetaP4
  bytes += npairs * nr * sizeof (double);
// pBetaP5
  bytes += npairs * nr * sizeof (double);
// pBetaP6
  bytes += npairs * nr * sizeof (double);
// pRepul
  bytes += npairs * nr * sizeof (double);
// pRepul1
  bytes += npairs * nr * sizeof (double);
// pRepul2
  bytes += npairs * nr * sizeof (double);
// pRepul3
  bytes += npairs * nr * sizeof (double);
// pRepul4
  bytes += npairs * nr * sizeof (double);
// pRepul5
  bytes += npairs * nr * sizeof (double);
// pRepul6
  bytes += npairs * nr * sizeof (double);
// FsigBO
  bytes += npairs * nr * sizeof (double);
// FsigBO1
  bytes += npairs * nr * sizeof (double);
// FsigBO2
  bytes += npairs * nr * sizeof (double);
// FsigBO3
  bytes += npairs * nr * sizeof (double);
// FsigBO4
  bytes += npairs * nr * sizeof (double);
// FsigBO5
  bytes += npairs * nr * sizeof (double);
// FsigBO6
  bytes += npairs * nr * sizeof (double);
// itypeSigBk
  bytes += neigh_total *neigh_ct* sizeof(int);
// nSigBk
  bytes += neigh_total * sizeof(int);
// sigB
  bytes += neigh_total * sizeof(int);
// sigB1
  bytes += neigh_total * sizeof(int);
// nPiBk
  bytes += neigh_total * sizeof(int);
// piB
  bytes += neigh_total * sizeof(int);
// itypePiBk
  bytes += neigh_total *neigh_ct* sizeof(int);
// BOP_index
    bytes += nall * sizeof(double);
  if(otfly==0) {
// cosAng
    bytes += cos_total* sizeof(double);
// dcAng
    bytes += cos_total * 3 * 2 * sizeof(double);
// disij
    bytes += neigh_total * 3 * sizeof(double);
// rij
    bytes += neigh_total * sizeof(double);
// betaS
    bytes += neigh_total * sizeof(double);
// dBetaS
    bytes += neigh_total * sizeof(double);
// betaP
    bytes += neigh_total * sizeof(double);
// dBetaP
    bytes += neigh_total * sizeof(double);
// repul
    bytes += neigh_total * sizeof(double);
// dRepul
    bytes += neigh_total * sizeof(double);
// cos_index
    bytes += nall * sizeof(double);
  }
// pi_a
  bytes += npairs * sizeof(double);
// pro_delta
  bytes += npairs * sizeof(double);
// pi_delta
  bytes += npairs * sizeof(double);
// pi_p
  bytes += npairs * sizeof(double);
// pi_c
  bytes += npairs * sizeof(double);
// sigma_r0
  bytes += npairs * sizeof(double);
// pi_r0
  bytes += npairs * sizeof(double);
// phi_r0
  bytes += npairs * sizeof(double);
// sigma_rc
  bytes += npairs * sizeof(double);
// pi_rc
  bytes += npairs * sizeof(double);
// pi_a
  bytes += npairs * sizeof(double);
// pro_delta
  bytes += npairs * sizeof(double);
// pi_delta
  bytes += npairs * sizeof(double);
// pi_p
  bytes += npairs * sizeof(double);
// pi_c
  bytes += npairs * sizeof(double);
// sigma_r0
  bytes += npairs * sizeof(double);
// pi_r0
  bytes += npairs * sizeof(double);
// phi_r0
  bytes += npairs * sizeof(double);
// sigma_rc
  bytes += npairs * sizeof(double);
// pi_rc
  bytes += npairs * sizeof(double);
// phi_rc
  bytes += npairs * sizeof(double);
// r1
  bytes += npairs * sizeof(double);
// sigma_beta0
  bytes += npairs * sizeof(double);
// pi_beta0
  bytes += npairs * sizeof(double);
// phi0
  bytes += npairs * sizeof(double);
// sigma_n
  bytes += npairs * sizeof(double);
// pi_n
  bytes += npairs * sizeof(double);
// phi_m
  bytes += npairs * sizeof(double);
// sigma_nc
  bytes += npairs * sizeof(double);
// pi_nc
  bytes += npairs * sizeof(double);
// phi_nc
  bytes += npairs * sizeof(double);
// pro
  bytes += npairs * sizeof(double);
// sigma_delta
  bytes += npairs * sizeof(double);
// sigma_c
  bytes += npairs * sizeof(double);
// sigma_a
  bytes += npairs * sizeof(double);
// sigma_g0
  bytes += bop_types * bop_types *bop_types * sizeof(double);
// sigma_g1
  bytes += bop_types * bop_types *bop_types * sizeof(double);
// sigma_g2
  bytes += bop_types * bop_types *bop_types * sizeof(double);
// sigma_g3
  bytes += bop_types * bop_types *bop_types * sizeof(double);
// sigma_g4
  bytes += bop_types * bop_types *bop_types * sizeof(double);
// sigma_f
  bytes += npairs * sizeof(double);
// sigma_k
  bytes += npairs * sizeof(double);
// small3
  bytes += npairs * sizeof(double);
// bt_pi
  bytes += maxneigh*(maxneigh/2) *sizeof(B_PI);
// bt_sigma
  bytes += maxneigh*(maxneigh/2) *sizeof(B_SG);

  return bytes;
}

/* ---------------------------------------------------------------------- */

void PairBOP::memory_theta_create()
{
  if(maxneigh<8)
    neigh_ct=(maxneigh-1)*(maxneigh-1)*(maxneigh-1);
  else
    neigh_ct=(maxneigh-1)*(maxneigh-1);
  memory->create(itypeSigBk,neigh_total
      ,neigh_ct,"itypeSigBk");
  memory->create(nSigBk,neigh_total,"nSigBk");
  memory->create(sigB,neigh_total,"sigB");
  memory->create(sigB1,neigh_total,"sigB1");
  memory->create(itypePiBk,neigh_total
      ,neigh_ct,"itypePiBk");
  memory->create(nPiBk,neigh_total,"nPiBk");
  memory->create(piB,neigh_total,"piB");
  memory->create(neigh_flag,neigh_total,"neigh_flag");
  if(otfly==0) {
    memory->create(cosAng,cos_total,"BOP:cosAng");
    memory->create(dcAng,cos_total*2,3,2,"BOP:dcAng");
    memory->create(disij,3,neigh_total,"disij");
    memory->create(rij,neigh_total,"rij");
    memory->create(betaS,neigh_total,"betaS");
    memory->create(dBetaS,neigh_total,"dBetaS");
    memory->create(betaP,neigh_total,"betaP");
    memory->create(dBetaP,neigh_total,"dBetaP");
    memory->create(repul,neigh_total,"repul");
    memory->create(dRepul,neigh_total,"dRepul");
  }
  update_list=1;
}

/* ---------------------------------------------------------------------- */

void PairBOP::memory_theta_grow()
{
  if(maxneigh<8)
    neigh_ct=(maxneigh-1)*(maxneigh-1)*(maxneigh-1);
  else
    neigh_ct=(maxneigh-1)*(maxneigh-1);
  memory->grow(itypeSigBk,neigh_total
      ,neigh_ct,"itypeSigBk");
  memory->grow(nSigBk,neigh_total,"nSigBk");
  memory->grow(sigB,neigh_total,"sigB");
  memory->grow(sigB1,neigh_total,"sigB1");
  memory->grow(itypePiBk,neigh_total
      ,neigh_ct,"itypePiBk");
  memory->grow(nPiBk,neigh_total,"nPiBk");
  memory->grow(piB,neigh_total,"piB");
  memory->grow(neigh_flag,neigh_total,"neigh_flag");
  if(otfly==0) {
    memory->grow(cosAng,cos_total,"BOP:cosAng");
    memory->grow(dcAng,cos_total*2,3,2,"BOP:dcAng");
    memory->grow(disij,3,neigh_total,"disij");
    memory->grow(rij,neigh_total,"rij");
    memory->grow(betaS,neigh_total,"betaS");
    memory->grow(dBetaS,neigh_total,"dBetaS");
    memory->grow(betaP,neigh_total,"betaP");
    memory->grow(dBetaP,neigh_total,"dBetaP");
    memory->grow(repul,neigh_total,"repul");
    memory->grow(dRepul,neigh_total,"dRepul");
  }
  update_list=1;
}

/* ---------------------------------------------------------------------- */

void PairBOP::memory_theta_destroy()
{

  memory->destroy(itypeSigBk);
  memory->destroy(nSigBk);
  memory->destroy(sigB);
  memory->destroy(sigB1);
  memory->destroy(itypePiBk);
  memory->destroy(nPiBk);
  memory->destroy(piB);
  memory->destroy(neigh_flag);
  if(otfly==0) {
    memory->destroy(cosAng);
    memory->destroy(dcAng);
    memory->destroy(disij);
    memory->destroy(rij);
    memory->destroy(betaS);
    memory->destroy(dBetaS);
    memory->destroy(betaP);
    memory->destroy(dBetaP);
    memory->destroy(repul);
    memory->destroy(dRepul);
  }
 update_list=0;
}

/* ---------------------------------------------------------------------- */

void PairBOP::create_pi(int n_tot)
{
  bt_pi = (B_PI *) memory->smalloc(n_tot*sizeof(B_PI),"BOP:bt_pi");
  allocate_pi=1;
}

void PairBOP::create_sigma(int n_tot)
{
  bt_sg = (B_SG *) memory->smalloc(n_tot*sizeof(B_SG),"BOP:bt_sg");
  allocate_sigma=1;
}

void PairBOP::destroy_pi()
{
  memory->destroy(bt_pi);
  allocate_pi=0;
}

void PairBOP::destroy_sigma()
{
  memory->destroy(bt_sg);
  allocate_sigma=0;
}

/* ---------------------------------------------------------------------- */

void PairBOP::grow_pi(int n1, int n2)
{
  int i,j;
  B_PI *bt_temp;
  bt_temp = (B_PI *) memory->smalloc(n1*sizeof(B_PI),"BOP:b_temp");
  for(i=0;i<n1;i++) {
    bt_temp[i].temp = bt_pi[i].temp;
    bt_temp[i].i = bt_pi[i].i;
    bt_temp[i].j = bt_pi[i].j;
    for(j=0;j<3;j++) {
      bt_temp[i].dAA[j] = bt_pi[i].dAA[j];
      bt_temp[i].dBB[j] = bt_pi[i].dBB[j];
      bt_temp[i].dPiB[j] = bt_pi[i].dPiB[j];
    }
  }
  memory->destroy(bt_pi);
  bt_pi=NULL;
  bt_pi = (B_PI *) memory->smalloc(n2*sizeof(B_PI),"BOP:bt_pi");
  for(i=0;i<n1;i++) {
    bt_pi[i].temp = bt_temp[i].temp;
    bt_pi[i].i = bt_temp[i].i;
    bt_pi[i].j = bt_temp[i].j;
    for(j=0;j<3;j++) {
      bt_pi[i].dAA[j] = bt_temp[i].dAA[j];
      bt_pi[i].dBB[j] = bt_temp[i].dBB[j];
      bt_pi[i].dPiB[j] = bt_temp[i].dPiB[j];
    }
  }
  for(i=n1;i<n2;i++) {
    bt_pi[i].i = -1;
    bt_pi[i].j = -1;
    for(j=0;j<3;j++) {
      bt_pi[i].dAA[j] = 0.0;
      bt_pi[i].dBB[j] = 0.0;
      bt_pi[i].dPiB[j] = 0.0;
    }
  }
  memory->destroy(bt_temp);
}

/* ---------------------------------------------------------------------- */

void PairBOP::grow_sigma(int n1,int n2)
{
  int i,j;
  B_SG *bt_temp;
  bt_temp = (B_SG *) memory->smalloc(n1*sizeof(B_SG),"BOP:bt_temp");
  for(i=0;i<n1;i++) {
    bt_temp[i].temp = bt_sg[i].temp;
    bt_temp[i].i = bt_sg[i].i;
    bt_temp[i].j = bt_sg[i].j;
    for(j=0;j<3;j++) {
      bt_temp[i].dAA[j] = bt_sg[i].dAA[j];
      bt_temp[i].dBB[j] = bt_sg[i].dBB[j];
      bt_temp[i].dCC[j] = bt_sg[i].dCC[j];
      bt_temp[i].dDD[j] = bt_sg[i].dDD[j];
      bt_temp[i].dEE[j] = bt_sg[i].dEE[j];
      bt_temp[i].dEE1[j] = bt_sg[i].dEE1[j];
      bt_temp[i].dFF[j] = bt_sg[i].dFF[j];
      bt_temp[i].dAAC[j] = bt_sg[i].dAAC[j];
      bt_temp[i].dBBC[j] = bt_sg[i].dBBC[j];
      bt_temp[i].dCCC[j] = bt_sg[i].dCCC[j];
      bt_temp[i].dDDC[j] = bt_sg[i].dDDC[j];
      bt_temp[i].dEEC[j] = bt_sg[i].dEEC[j];
      bt_temp[i].dFFC[j] = bt_sg[i].dFFC[j];
      bt_temp[i].dGGC[j] = bt_sg[i].dGGC[j];
      bt_temp[i].dUT[j] = bt_sg[i].dUT[j];
      bt_temp[i].dSigB1[j] = bt_sg[i].dSigB1[j];
      bt_temp[i].dSigB[j] = bt_sg[i].dSigB[j];
    }
  }
  memory->destroy(bt_sg);
  bt_sg=NULL;
  bt_sg = (B_SG *) memory->smalloc(n2*sizeof(B_SG),"BOP:bt_sg");
  for(i=0;i<n1;i++) {
    bt_sg[i].temp = bt_temp[i].temp;
    bt_sg[i].i = bt_temp[i].i;
    bt_sg[i].j = bt_temp[i].j;
    for(j=0;j<3;j++) {
      bt_sg[i].dAA[j] = bt_temp[i].dAA[j];
      bt_sg[i].dBB[j] = bt_temp[i].dBB[j];
      bt_sg[i].dCC[j] = bt_temp[i].dCC[j];
      bt_sg[i].dDD[j] = bt_temp[i].dDD[j];
      bt_sg[i].dEE[j] = bt_temp[i].dEE[j];
      bt_sg[i].dEE1[j] = bt_temp[i].dEE1[j];
      bt_sg[i].dFF[j] = bt_temp[i].dFF[j];
      bt_sg[i].dAAC[j] = bt_temp[i].dAAC[j];
      bt_sg[i].dBBC[j] = bt_temp[i].dBBC[j];
      bt_sg[i].dCCC[j] = bt_temp[i].dCCC[j];
      bt_sg[i].dDDC[j] = bt_temp[i].dDDC[j];
      bt_sg[i].dEEC[j] = bt_temp[i].dEEC[j];
      bt_sg[i].dFFC[j] = bt_temp[i].dFFC[j];
      bt_sg[i].dGGC[j] = bt_temp[i].dGGC[j];
      bt_sg[i].dUT[j] = bt_temp[i].dUT[j];
      bt_sg[i].dSigB1[j] = bt_temp[i].dSigB1[j];
      bt_sg[i].dSigB[j] = bt_temp[i].dSigB[j];
    }
  }
  for(i=n1;i<n2;i++) {
    bt_sg[i].i = -1;
    bt_sg[i].j = -1;
    for(j=0;j<3;j++) {
      bt_sg[i].dAA[j] = 0.0;
      bt_sg[i].dBB[j] = 0.0;
      bt_sg[i].dCC[j] = 0.0;
      bt_sg[i].dDD[j] = 0.0;
      bt_sg[i].dEE[j] = 0.0;
      bt_sg[i].dEE1[j] = 0.0;
      bt_sg[i].dFF[j] = 0.0;
      bt_sg[i].dAAC[j] = 0.0;
      bt_sg[i].dBBC[j] = 0.0;
      bt_sg[i].dCCC[j] = 0.0;
      bt_sg[i].dDDC[j] = 0.0;
      bt_sg[i].dEEC[j] = 0.0;
      bt_sg[i].dFFC[j] = 0.0;
      bt_sg[i].dGGC[j] = 0.0;
      bt_sg[i].dUT[j] = 0.0;
      bt_sg[i].dSigB1[j] = 0.0;
      bt_sg[i].dSigB[j] = 0.0;
    }
  }
  memory->destroy(bt_temp);
}
