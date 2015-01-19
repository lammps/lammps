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
   Contributing author: Dominik Wójt (Wroclaw University of Technology)
     based on pair_airebo by Ase Henry (MIT)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include "pair_lcbop.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "my_page.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 1024
#define TOL 1.0e-9
#define PGDELTA 1

/* ---------------------------------------------------------------------- */

PairLCBOP::PairLCBOP(LAMMPS *lmp) : Pair(lmp) {
  single_enable = 0;
  one_coeff = 1;
  manybody_flag = 1;
  ghostneigh = 1;

  maxlocal = 0;
  SR_numneigh = NULL;
  SR_firstneigh = NULL;
  ipage = NULL;
  pgsize = oneatom = 0;

  N = NULL;
  M = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairLCBOP::~PairLCBOP()
{
  memory->destroy(SR_numneigh);
  memory->sfree(SR_firstneigh);
  delete [] ipage;
  memory->destroy(N);
  memory->destroy(M);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cutghost);

    delete [] map;
  }
}

/* ---------------------------------------------------------------------- */

void PairLCBOP::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = vflag_atom = 0;

  SR_neigh();
  FSR(eflag,vflag);
  FLR(eflag,vflag);

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLCBOP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cutghost,n+1,n+1,"pair:cutghost");

  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLCBOP::settings(int narg, char **arg) {
  if( narg != 0 ) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLCBOP::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom types to C and NULL
  // map[i] = which element (0 for C) the Ith atom type is, -1 if NULL

  for (int i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
    } else if (strcmp(arg[i],"C") == 0) {
      map[i-2] = 0;
    } else error->all(FLERR,"Incorrect args for pair coefficients");
  }

  // read potential file and initialize fitting splines

  read_file(arg[2]);
  spline_init();

  // clear setflag since coeff() called once with I,J = * *

  int n = atom->ntypes;
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

void PairLCBOP::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style LCBOP requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style LCBOP requires newton pair on");

  // need a full neighbor list, including neighbors of ghosts

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->ghost = 1;

  // local SR neighbor list
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
      ipage[i].init(oneatom,pgsize,PGDELTA);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLCBOP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  // cut3rebo = 3 SR distances

  cut3rebo = 3.0 * r_2;

  // cutmax = furthest distance from an owned atom
  //          at which another atom will feel force, i.e. the ghost cutoff
  // for SR term in potential:
  //   interaction = M-K-I-J-L-N with I = owned and J = ghost
  //   I to N is max distance = 3 SR distances
  // for V_LR term in potential:
  //   r_2_LR
  // cutghost = SR cutoff used in SR_neigh() for neighbors of ghosts

  double cutmax = MAX( cut3rebo,r_2_LR );

  cutghost[i][j] = r_2;
  cutLRsq = r_2_LR*r_2_LR;

  cutghost[j][i] = cutghost[i][j];

  r_2_sq = r_2*r_2;

  return cutmax;
}

/* ----------------------------------------------------------------------
   create SR neighbor list from main neighbor list
   SR neighbor list stores neighbors of ghost atoms
------------------------------------------------------------------------- */

void PairLCBOP::SR_neigh()
{
  int i,j,ii,jj,n,allnum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,dS;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *neighptr;

  double **x = atom->x;

  if (atom->nmax > maxlocal) {  // ensure ther is enough space
    maxlocal = atom->nmax;      // for atoms and ghosts allocated
    memory->destroy(SR_numneigh);
    memory->sfree(SR_firstneigh);
    memory->destroy(N);
    memory->destroy(M);
    memory->create(SR_numneigh,maxlocal,"LCBOP:numneigh");
    SR_firstneigh = (int **) memory->smalloc(maxlocal*sizeof(int *),
                           "LCBOP:firstneigh");
    memory->create(N,maxlocal,"LCBOP:N");
    memory->create(M,maxlocal,"LCBOP:M");
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

    n = 0;
    neighptr = ipage->vget();

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    N[i] = 0.0;
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < r_2_sq) {
        neighptr[n++] = j;
        N[i] += f_c(sqrt(rsq),r_1,r_2,&dS);
      }
    }

    SR_firstneigh[i] = neighptr;
    SR_numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  // calculate M_i

  for (ii = 0; ii < allnum; ii++) {
    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    M[i] = 0.0;

    jlist = SR_firstneigh[i];
    jnum = SR_numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < r_2_sq) {
        double f_c_ij = f_c(sqrt(rsq),r_1,r_2,&dS);
        double Nji = N[j]-f_c_ij;
        // F(xij) = 1-f_c_LR(Nji, 2,3,&dummy)
        M[i] += f_c_ij * ( 1-f_c_LR(Nji, 2,3,&dS) );
      }
    }
  }
}

/* ----------------------------------------------------------------------
  Short range forces and energy
------------------------------------------------------------------------- */

void PairLCBOP::FSR(int eflag, int vflag)
{
  int i,j,jj,ii,inum;
  tagint itag,jtag;
  double delx,dely,delz,fpair,xtmp,ytmp,ztmp;
  double r_sq,rijmag,f_c_ij,df_c_ij;
  double VR,dVRdi,VA,Bij,dVAdi,dVA;
  double del[3];
  int *ilist,*SR_neighs;

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;

  // two-body interactions from SR neighbor list, skip half of them

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    SR_neighs = SR_firstneigh[i];

    for (jj = 0; jj < SR_numneigh[i]; jj++) {
      j = SR_neighs[jj];
      jtag = tag[j];

      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];
      r_sq = delx*delx + dely*dely + delz*delz;
      rijmag = sqrt(r_sq);
      f_c_ij = f_c( rijmag,r_1,r_2,&df_c_ij );
      if( f_c_ij <= TOL ) continue;

      VR = A*exp(-alpha*rijmag);
      dVRdi = -alpha*VR;
      dVRdi = dVRdi*f_c_ij + df_c_ij*VR; // VR -> VR * f_c_ij
      VR *= f_c_ij;

      VA = dVA = 0.0;
      {
        double term = B_1 * exp(-beta_1*rijmag);
        VA += term;
        dVA += -beta_1 * term;
        term = B_2 * exp(-beta_2*rijmag);
        VA += term;
        dVA += -beta_2 * term;
      }
      dVA = dVA*f_c_ij + df_c_ij*VA; // VA -> VA * f_c_ij
      VA *= f_c_ij;
      del[0] = delx;
      del[1] = dely;
      del[2] = delz;
      Bij = bondorder(i,j,del,rijmag,VA,f,vflag_atom);
      dVAdi = Bij*dVA;

      // F = (dVRdi+dVAdi)*(-grad rijmag)
      // grad_i rijmag =  \vec{rij} /rijmag
      // grad_j rijmag = -\vec{rij} /rijmag
      fpair = -(dVRdi-dVAdi) / rijmag;
      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      double evdwl=0.0;
      if (eflag) evdwl = VR - Bij*VA;
      if (evflag) ev_tally(i,j,nlocal,newton_pair,
        evdwl,0.0,fpair,delx,dely,delz);
    }
  }
}

/* ----------------------------------------------------------------------
   compute long range forces and energy
------------------------------------------------------------------------- */

void PairLCBOP::FLR(int eflag, int vflag)
{
  int i,j,jj,ii;
  tagint itag,jtag;
  double delx,dely,delz,fpair,xtmp,ytmp,ztmp;
  double r_sq,rijmag,f_c_ij,df_c_ij;
  double V,dVdi;

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // two-body interactions from full neighbor list, skip half of them

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    int *neighs = firstneigh[i];

    for (jj = 0; jj < numneigh[i]; jj++) {
      j = neighs[jj];
      j &= NEIGHMASK;
      jtag = tag[j];

      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];
      r_sq = delx*delx + dely*dely + delz*delz;
      rijmag = sqrt(r_sq);
      f_c_ij = 1-f_c( rijmag,r_1,r_2,&df_c_ij );
      df_c_ij = -df_c_ij;
      // derivative may be inherited from previous call, see f_c_LR definition
      f_c_ij *= f_c_LR( rijmag, r_1_LR, r_2_LR, &df_c_ij );
      if( f_c_ij <= TOL ) continue;

      V = dVdi = 0;
      if( rijmag<r_0 ) {
        double exp_part = exp( -lambda_1*(rijmag-r_0) );
        V = eps_1*( exp_part*exp_part - 2*exp_part) + v_1;
        dVdi = 2*eps_1*lambda_1*exp_part*( 1-exp_part );
      } else {
        double exp_part = exp( -lambda_2*(rijmag-r_0) );
        V = eps_2*( exp_part*exp_part - 2*exp_part) + v_2;
        dVdi = 2*eps_2*lambda_2*exp_part*( 1-exp_part );
      }
      dVdi = dVdi*f_c_ij + df_c_ij*V; // V -> V * f_c_ij
      V *= f_c_ij;

      // F = (dVdi)*(-grad rijmag)
      // grad_i rijmag =  \vec{rij} /rijmag
      // grad_j rijmag = -\vec{rij} /rijmag
      fpair = -dVdi / rijmag;
      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      double evdwl=0.0;
      if (eflag) evdwl = V;
      if (evflag) ev_tally(i,j,nlocal,newton_pair,
        evdwl,0.0,fpair,delx,dely,delz);
    }
  }
}

/* ----------------------------------------------------------------------
   forces for Nij and Mij
------------------------------------------------------------------------- */

void PairLCBOP::FNij( int i, int j, double factor, double **f, int vflag_atom ) {
  int atomi = i;
  int atomj = j;
  int *SR_neighs = SR_firstneigh[i];
  double **x = atom->x;
  for( int k=0; k<SR_numneigh[i]; k++ ) {
    int atomk = SR_neighs[k];
    if (atomk != atomj) {
      double rik[3];
      rik[0] = x[atomi][0]-x[atomk][0];
      rik[1] = x[atomi][1]-x[atomk][1];
      rik[2] = x[atomi][2]-x[atomk][2];
      double riksq = (rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]);
      if( riksq > r_1*r_1 ) { // &&  riksq < r_2*r_2, if second condition not fulfilled neighbor would not be in the list
        double rikmag = sqrt(riksq);
        double df_c_ik;
        f_c( rikmag, r_1, r_2, &df_c_ik );

        // F = factor*df_c_ik*(-grad rikmag)
        // grad_i rikmag =  \vec{rik} /rikmag
        // grad_k rikmag = -\vec{rik} /rikmag
        double fpair = -factor*df_c_ik / rikmag;
        f[atomi][0] += rik[0]*fpair;
        f[atomi][1] += rik[1]*fpair;
        f[atomi][2] += rik[2]*fpair;
        f[atomk][0] -= rik[0]*fpair;
        f[atomk][1] -= rik[1]*fpair;
        f[atomk][2] -= rik[2]*fpair;

        if (vflag_atom) v_tally2(atomi,atomk,fpair,rik);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLCBOP::FMij( int i, int j, double factor, double **f, int vflag_atom ) {
  int atomi = i;
  int atomj = j;
  int *SR_neighs = SR_firstneigh[i];
  double **x = atom->x;
  for( int k=0; k<SR_numneigh[i]; k++ ) {
    int atomk = SR_neighs[k];
    if (atomk != atomj) {
      double rik[3];
      rik[0] = x[atomi][0]-x[atomk][0];
      rik[1] = x[atomi][1]-x[atomk][1];
      rik[2] = x[atomi][2]-x[atomk][2];
      double rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      double df_c_ik;
      double f_c_ik = f_c( rikmag, r_1, r_2, &df_c_ik );
      double Nki = N[k]-(f_c_ik);
//      double Mij = M[i] - f_c_ij*( 1-f_c(Nji, 2,3,&dummy) );
      double dF=0;
      double Fx = 1-f_c_LR(Nki, 2,3,&dF);
      dF = -dF;

      if( df_c_ik > TOL ) {
        double factor2 = factor*df_c_ik*Fx;
        // F = factor2*(-grad rikmag)
        // grad_i rikmag =  \vec{rik} /rikmag
        // grad_k rikmag = -\vec{rik} /rikmag
        double fpair = -factor2 / rikmag;
        f[atomi][0] += rik[0]*fpair;
        f[atomi][1] += rik[1]*fpair;
        f[atomi][2] += rik[2]*fpair;
        f[atomk][0] -= rik[0]*fpair;
        f[atomk][1] -= rik[1]*fpair;
        f[atomk][2] -= rik[2]*fpair;
        if (vflag_atom) v_tally2(atomi,atomk,fpair,rik);
      }

      if( dF > TOL ) {
        double factor2 = factor*f_c_ik*dF;
        FNij( atomk, atomi, factor2, f, vflag_atom );
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Bij function
------------------------------------------------------------------------- */

double PairLCBOP::bondorder(int i, int j, double rij[3],
    double rijmag, double VA,
    double **f, int vflag_atom)
{

  double bij, bji;
  /* bij & bji */{
    double rji[3];
    rji[0] = -rij[0]; rji[1] = -rij[1]; rji[2] = -rij[2];
    bij = b(i,j,rij,rijmag,VA,f,vflag_atom);
    bji = b(j,i,rji,rijmag,VA,f,vflag_atom);
  }

  double Fij_conj;
  /* F_conj */{
    double dummy;

    double df_c_ij;
    double f_c_ij = f_c( rijmag, r_1, r_2, &df_c_ij );
    double Nij = MIN( 3, N[i]-(f_c_ij) );
    double Nji = MIN( 3, N[j]-(f_c_ij) );

    // F(xij) = 1-f_c(Nji, 2,3,&dummy)
    double Mij = M[i] - f_c_ij*( 1-f_c(Nji, 2,3,&dummy) );
    double Mji = M[j] - f_c_ij*( 1-f_c(Nij, 2,3,&dummy) );
    Mij = MIN( Mij, 3 );
    Mji = MIN( Mji, 3 );

    double Nij_el, dNij_el_dNij, dNij_el_dMij;
    double Nji_el, dNji_el_dNji, dNji_el_dMji;
    {
      double num_Nij_el = 4 - Mij;
      double num_Nji_el = 4 - Mji;
      double den_Nij_el = Nij + 1 - Mij;
      double den_Nji_el = Nji + 1 - Mji;
      Nij_el = num_Nij_el / den_Nij_el;
      Nji_el = num_Nji_el / den_Nji_el;
      dNij_el_dNij = -Nij_el/den_Nij_el;
      dNji_el_dNji = -Nji_el/den_Nji_el;
      dNij_el_dMij = ( -1 + Nij_el ) /den_Nij_el;
      dNji_el_dMji = ( -1 + Nji_el ) /den_Nji_el;
    }

    double Nconj;
    double dNconj_dNij;
    double dNconj_dNji;
    double dNconj_dNel;
    {
      double num_Nconj = ( Nij+1 )*( Nji+1 )*( Nij_el+Nji_el ) - 4*( Nij+Nji+2);
      double den_Nconj = Nij*( 3-Nij )*( Nji+1 ) + Nji*( 3-Nji )*( Nij+1 ) + eps;
      Nconj = num_Nconj / den_Nconj;
      if( Nconj <= 0 ) {
        Nconj = 0;
        dNconj_dNij = 0;
        dNconj_dNji = 0;
        dNconj_dNel = 0;
      } else if( Nconj >= 1 ) {
        Nconj = 1;
        dNconj_dNij = 0;
        dNconj_dNji = 0;
        dNconj_dNel = 0;
      } else {
        dNconj_dNij = (
            ( (Nji+1)*(Nij_el + Nji_el)-4)
            - Nconj*( (Nji+1)*(3-2*Nij) + Nji*(3-Nji) )
          ) /den_Nconj;
        dNconj_dNji = (
            ( (Nij+1)*(Nji_el + Nij_el)-4)
            - Nconj*( (Nij+1)*(3-2*Nji) + Nij*(3-Nij) )
          ) /den_Nconj;
        dNconj_dNel = (Nij+1)*(Nji+1) / den_Nconj;
      }
    }

    double dF_dNij, dF_dNji, dF_dNconj;
    Fij_conj = F_conj( Nij, Nji, Nconj, &dF_dNij, &dF_dNji, &dF_dNconj );

    /*forces for Nij*/
    if( 3-Nij > TOL ) {
      double factor = -VA*0.5*( dF_dNij + dF_dNconj*( dNconj_dNij + dNconj_dNel*dNij_el_dNij ) );
      FNij( i, j, factor, f, vflag_atom );
    }
    /*forces for Nji*/
    if( 3-Nji > TOL ) {
      double factor = -VA*0.5*( dF_dNji + dF_dNconj*( dNconj_dNji + dNconj_dNel*dNji_el_dNji ) );
      FNij( j, i, factor, f, vflag_atom );
    }
    /*forces for Mij*/
    if( 3-Mij > TOL ) {
      double factor = -VA*0.5*( dF_dNconj*dNconj_dNel*dNij_el_dMij );
      FMij( i, j, factor, f, vflag_atom );
    }
    if( 3-Mji > TOL ) {
      double factor = -VA*0.5*( dF_dNconj*dNconj_dNel*dNji_el_dMji );
      FMij( j, i, factor, f, vflag_atom );
    }
  }


  double Bij = 0.5*( bij + bji + Fij_conj );
  return Bij;
}

/* ----------------------------------------------------------------------
  bij function
------------------------------------------------------------------------- */

double PairLCBOP::b(int i, int j, double rij[3],
                 double rijmag, double VA,
                 double **f, int vflag_atom) {
  int *SR_neighs = SR_firstneigh[i];
  double **x = atom->x;
  int atomi = i;
  int atomj = j;

  //calculate bij magnitude
  double bij = 1.0;
  for (int k = 0; k < SR_numneigh[i]; k++) {
    int atomk = SR_neighs[k];
    if (atomk != atomj) {
      double rik[3];
      rik[0] = x[atomi][0]-x[atomk][0];
      rik[1] = x[atomi][1]-x[atomk][1];
      rik[2] = x[atomi][2]-x[atomk][2];
      double rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      double delta_ijk = rijmag-rikmag;
      double dummy;
      double f_c_ik = f_c( rikmag, r_1, r_2, &dummy );
      double cos_ijk = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2]))
                / (rijmag*rikmag);
      cos_ijk = MIN(cos_ijk,1.0);
      cos_ijk = MAX(cos_ijk,-1.0);

      double G = gSpline(cos_ijk,   &dummy);
      double H = hSpline(delta_ijk, &dummy);
      bij += (f_c_ik*G*H);
    }
  }
  bij = pow( bij, -delta );

  // bij forces

  for (int k = 0; k < SR_numneigh[i]; k++) {
    int atomk = SR_neighs[k];
    if (atomk != atomj) {
      double rik[3];
      rik[0] = x[atomi][0]-x[atomk][0];
      rik[1] = x[atomi][1]-x[atomk][1];
      rik[2] = x[atomi][2]-x[atomk][2];
      double rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      double delta_ijk = rijmag-rikmag;
      double df_c_ik;
      double f_c_ik = f_c( rikmag, r_1, r_2, &df_c_ik );
      double cos_ijk = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2]))
                / (rijmag*rikmag);
      cos_ijk = MIN(cos_ijk,1.0);
      cos_ijk = MAX(cos_ijk,-1.0);

      double dcos_ijk_dri[3],dcos_ijk_drj[3],dcos_ijk_drk[3];
      dcos_ijk_drj[0] = -rik[0] / (rijmag*rikmag)
             + cos_ijk * rij[0] / (rijmag*rijmag);
      dcos_ijk_drj[1] = -rik[1] / (rijmag*rikmag)
             + cos_ijk * rij[1] / (rijmag*rijmag);
      dcos_ijk_drj[2] = -rik[2] / (rijmag*rikmag)
             + cos_ijk * rij[2] / (rijmag*rijmag);

      dcos_ijk_drk[0] = -rij[0] / (rijmag*rikmag)
             + cos_ijk * rik[0] / (rikmag*rikmag);
      dcos_ijk_drk[1] = -rij[1] / (rijmag*rikmag)
             + cos_ijk * rik[1] / (rikmag*rikmag);
      dcos_ijk_drk[2] = -rij[2] / (rijmag*rikmag)
             + cos_ijk * rik[2] / (rikmag*rikmag);

      dcos_ijk_dri[0] = -dcos_ijk_drk[0] - dcos_ijk_drj[0];
      dcos_ijk_dri[1] = -dcos_ijk_drk[1] - dcos_ijk_drj[1];
      dcos_ijk_dri[2] = -dcos_ijk_drk[2] - dcos_ijk_drj[2];

      double dG, dH;
      double G = gSpline( cos_ijk,   &dG );
      double H = hSpline( delta_ijk, &dH );
      double tmp = -VA*0.5*(-0.5*bij*bij*bij);

      double fi[3], fj[3], fk[3];

      double tmp2 = -tmp*df_c_ik*G*H/rikmag;
      // F = tmp*df_c_ik*G*H*(-grad rikmag)
      // grad_i rikmag =  \vec{rik} /rikmag
      // grad_k rikmag = -\vec{rik} /rikmag
      fi[0] =  tmp2*rik[0];
      fi[1] =  tmp2*rik[1];
      fi[2] =  tmp2*rik[2];
      fk[0] = -tmp2*rik[0];
      fk[1] = -tmp2*rik[1];
      fk[2] = -tmp2*rik[2];


      tmp2 = -tmp*f_c_ik*dG*H;
      // F = tmp*f_c_ik*dG*H*(-grad cos_ijk)
      // grad_i cos_ijk = dcos_ijk_dri
      // grad_j cos_ijk = dcos_ijk_drj
      // grad_k cos_ijk = dcos_ijk_drk
      fi[0] += tmp2*dcos_ijk_dri[0];
      fi[1] += tmp2*dcos_ijk_dri[1];
      fi[2] += tmp2*dcos_ijk_dri[2];
      fj[0] =  tmp2*dcos_ijk_drj[0];
      fj[1] =  tmp2*dcos_ijk_drj[1];
      fj[2] =  tmp2*dcos_ijk_drj[2];
      fk[0] += tmp2*dcos_ijk_drk[0];
      fk[1] += tmp2*dcos_ijk_drk[1];
      fk[2] += tmp2*dcos_ijk_drk[2];

      tmp2 = -tmp*f_c_ik*G*dH;
      // F = tmp*f_c_ik*G*dH*(-grad delta_ijk)
      // grad_i delta_ijk =  \vec{rij} /rijmag - \vec{rik} /rijmag
      // grad_j delta_ijk = -\vec{rij} /rijmag
      // grad_k delta_ijk =  \vec{rik} /rikmag
      fi[0] += tmp2*( rij[0]/rijmag - rik[0]/rikmag );
      fi[1] += tmp2*( rij[1]/rijmag - rik[1]/rikmag );
      fi[2] += tmp2*( rij[2]/rijmag - rik[2]/rikmag );
      fj[0] += tmp2*( -rij[0]/rijmag );
      fj[1] += tmp2*( -rij[1]/rijmag );
      fj[2] += tmp2*( -rij[2]/rijmag );
      fk[0] += tmp2*( rik[0]/rikmag );
      fk[1] += tmp2*( rik[1]/rikmag );
      fk[2] += tmp2*( rik[2]/rikmag );

      f[atomi][0] += fi[0]; f[atomi][1] += fi[1]; f[atomi][2] += fi[2];
      f[atomj][0] += fj[0]; f[atomj][1] += fj[1]; f[atomj][2] += fj[2];
      f[atomk][0] += fk[0]; f[atomk][1] += fk[1]; f[atomk][2] += fk[2];

      if (vflag_atom) {
        double rji[3], rki[3];
        rji[0] = -rij[0]; rji[1] = -rij[1]; rji[2] = -rij[2];
        rki[0] = -rik[0]; rki[1] = -rik[1]; rki[2] = -rik[2];
        v_tally3(atomi,atomj,atomk,fj,fk,rji,rki);
      }
    }
  }

  return bij;
}

/* ----------------------------------------------------------------------
   spline interpolation for G
------------------------------------------------------------------------- */

void PairLCBOP::g_decompose_x( double x, size_t *field_idx, double *offset ) {
  size_t i=0;
  while( i<(6-1) && !( x<gX[i+1] ) )
    i++;
  *field_idx = i;
  *offset = ( x - gX[i] );
}

/* ---------------------------------------------------------------------- */

double PairLCBOP::gSpline( double x, double *dgdc ) {
  size_t i;
  double x_n;
  g_decompose_x( x, &i, &x_n );
  double sum = 0;
  *dgdc = 0;
  double pow_x_n = 1.0;
  for( size_t j=0; j<5; j++ ) {
      sum += gC[j][i]*pow_x_n;
      *dgdc += gC[j+1][i]*(j+1)*pow_x_n;
      pow_x_n *= x_n;
  }
  sum += gC[5][i]*pow_x_n;
  return sum;
}

/* ---------------------------------------------------------------------- */

double PairLCBOP::hSpline( double x, double *dhdx ) {
  if( x < -d ) {
      double z = kappa*( x+d );
      double y = pow(z, 10.0);
      double w = pow( 1+y, -0.1 );
      *dhdx = kappa*L*w/(1+y);
      return L*( 1 + z*w );
    }
    if( x > d ) {
      *dhdx = R_1;
      return R_0 + R_1*( x-d );
    }

      double result = 1 + C_1*x;
      *dhdx    = C_1*result;
    double pow_x = x*x;
      result  += 0.5*C_1*C_1*pow_x;
    pow_x *= x;// == x^3
      *dhdx   += 4*C_4*pow_x;
    pow_x *= x;// == x^4
      result  += C_4*pow_x;
    pow_x *= x;// == x^5
      *dhdx   += 6*C_6*pow_x;
    pow_x *= x;// == x^5
      result += C_6*pow_x;
    return result;
}

/* ---------------------------------------------------------------------- */

double PairLCBOP::F_conj( double N_ij, double N_ji, double N_conj_ij, double *dFN_ij, double *dFN_ji, double *dFN_ij_conj ) {
  size_t N_ij_int         = MIN( static_cast<size_t>( floor( N_ij ) ), 2 ); // 2 is the highest number of field
  size_t N_ji_int         = MIN( static_cast<size_t>( floor( N_ji ) ), 2 ); // cast to suppress warning
  double x                = N_ij - N_ij_int;
  double y                = N_ji - N_ji_int;
  const TF_conj_field &f0 = F_conj_field[N_ij_int][N_ji_int][0];
  const TF_conj_field &f1 = F_conj_field[N_ij_int][N_ji_int][1];
  double F_0 = 0;
  double F_1 = 0;
  double dF_0_dx = 0, dF_0_dy = 0;
  double dF_1_dx = 0, dF_1_dy = 0;
  double l, r;
  if( N_conj_ij < 1 ) {
    l = (1-y)* (1-x);   r = ( f0.f_00 + x*     x*   f0.f_x_10   + y*     y*   f0.f_y_01 );    F_0 += l*r;   dF_0_dx += -(1-y)*r +l*2*x*    f0.f_x_10;    dF_0_dy += -(1-x)*r +l*2*y*    f0.f_y_01;
    l = (1-y)*  x;      r = ( f0.f_10 + (1-x)*(1-x)*f0.f_x_00   + y*     y*   f0.f_y_11 );    F_0 += l*r;   dF_0_dx +=  (1-y)*r -l*2*(1-x)*f0.f_x_00;    dF_0_dy += -x*    r +l*2*y*    f0.f_y_11;
    l = y*     (1-x);   r = ( f0.f_01 + x*     x*   f0.f_x_11   + (1-y)*(1-y)*f0.f_y_00 );    F_0 += l*r;   dF_0_dx += -y*    r +l*2*x*    f0.f_x_11;    dF_0_dy +=  (1-x)*r -l*2*(1-y)*f0.f_y_00;
    l = y*      x;      r = ( f0.f_11 + (1-x)*(1-x)*f0.f_x_01   + (1-y)*(1-y)*f0.f_y_10 );    F_0 += l*r;   dF_0_dx +=  y*    r -l*2*(1-x)*f0.f_x_01;    dF_0_dy +=  x*    r -l*2*(1-y)*f0.f_y_10;
  }
  if( N_conj_ij > 0 ) {
    l = (1-y)* (1-x);   r = ( f0.f_00 + x*     x*   f1.f_x_10   + y*     y*   f1.f_y_01 );    F_1 += l*r;   dF_1_dx += -(1-y)*r +l*2*x*    f1.f_x_10;    dF_1_dy += -(1-x)*r +l*2*y*    f1.f_y_01;
    l = (1-y)*  x;      r = ( f1.f_10 + (1-x)*(1-x)*f1.f_x_00   + y*     y*   f1.f_y_11 );    F_1 += l*r;   dF_1_dx +=  (1-y)*r -l*2*(1-x)*f1.f_x_00;    dF_1_dy += -x*    r +l*2*y*    f1.f_y_11;
    l = y*     (1-x);   r = ( f1.f_01 + x*     x*   f1.f_x_11   + (1-y)*(1-y)*f1.f_y_00 );    F_1 += l*r;   dF_1_dx += -y*    r +l*2*x*    f1.f_x_11;    dF_1_dy +=  (1-x)*r -l*2*(1-y)*f1.f_y_00;
    l = y*      x;      r = ( f1.f_11 + (1-x)*(1-x)*f1.f_x_01   + (1-y)*(1-y)*f1.f_y_10 );    F_1 += l*r;   dF_1_dx +=  y*    r -l*2*(1-x)*f1.f_x_01;    dF_1_dy +=  x*    r -l*2*(1-y)*f1.f_y_10;
  }
  double result = (1-N_conj_ij)*F_0 + N_conj_ij*F_1;
  *dFN_ij = (1-N_conj_ij)*dF_0_dx + N_conj_ij*dF_1_dx;
  *dFN_ji = (1-N_conj_ij)*dF_0_dy + N_conj_ij*dF_1_dy;
  *dFN_ij_conj = -F_0 + F_1;

  return result;
}

/* ----------------------------------------------------------------------
   read LCBOP potential file
------------------------------------------------------------------------- */

void PairLCBOP::read_file(char *filename)
{
  int i,k,l;
  char s[MAXLINE];

  MPI_Comm_rank(world,&me);

  // read file on proc 0

  if (me == 0) {
    FILE *fp = force->open_potential(filename);
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open LCBOP potential file %s",filename);
      error->one(FLERR,str);
    }

    // skip initial comment lines

    while (1) {
      fgets(s,MAXLINE,fp);
      if (s[0] != '#') break;
    }

    // read parameters

    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&r_1);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&r_2);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&gamma_1);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&A);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&B_1);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&B_2);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&alpha);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&beta_1);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&beta_2);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&d);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&C_1);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&C_4);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&C_6);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&L);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&kappa);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&R_0);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&R_1);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&r_0);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&r_1_LR);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&r_2_LR);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&v_1);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&v_2);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&eps_1);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&eps_2);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&lambda_1);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&lambda_2);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&eps);
    fgets(s,MAXLINE,fp);    sscanf(s,"%lg",&delta);

    while (1) {
      fgets(s,MAXLINE,fp);
      if (s[0] != '#') break;
    }

    // F_conj spline

    for (k = 0; k < 2; k++) { // 2 values of N_ij_conj
      for (l = 0; l < 3; l++) { // 3 types of data: f, dfdx, dfdy
        for (i = 0; i < 4; i++) { // 4x4 matrix
          fgets(s,MAXLINE,fp);
          sscanf(s,"%lg %lg %lg %lg",
            &F_conj_data[i][0][k][l],
            &F_conj_data[i][1][k][l],
            &F_conj_data[i][2][k][l],
            &F_conj_data[i][3][k][l]);
        }
        while (1) { fgets(s,MAXLINE,fp); if (s[0] != '#') break; }
      }
    }

    // G spline

    // x coordinates of mesh points
    fgets(s,MAXLINE,fp);
    sscanf( s,"%lg %lg %lg %lg %lg %lg",
      &gX[0], &gX[1], &gX[2],
      &gX[3], &gX[4], &gX[5] );

    for (i = 0; i < 6; i++) { // for each power in polynomial
      fgets(s,MAXLINE,fp);
      sscanf( s,"%lg %lg %lg %lg %lg",
        &gC[i][0], &gC[i][1], &gC[i][2],
        &gC[i][3], &gC[i][4] );
    }

    fclose(fp);
  }

  // broadcast read-in and setup values

  MPI_Bcast(&r_1      ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&r_2      ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&gamma_1  ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&A        ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&B_1      ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&B_2      ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&alpha    ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&beta_1   ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&beta_2   ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&d        ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&C_1      ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&C_4      ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&C_6      ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&L        ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&kappa    ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&R_0      ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&R_1      ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&r_0      ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&r_1_LR   ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&r_2_LR   ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&v_1      ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&v_2      ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&eps_1    ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&eps_2    ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&lambda_1 ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&lambda_2 ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&eps      ,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&delta    ,1,MPI_DOUBLE,0,world);

  MPI_Bcast(&gX[0]    ,6,MPI_DOUBLE,0,world);
  MPI_Bcast(&gC[0][0] ,(6-1)*(5+1),MPI_DOUBLE,0,world);

  MPI_Bcast(&F_conj_data[0],6*4*4,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   init coefficients for TF_conj
------------------------------------------------------------------------- */

#include <iostream>
#include <fstream>
#include <functional>
template< class function > void print_function( double x_0, double x_1, size_t n, function f, std::ostream &stream ) {
  double dx = (x_1-x_0)/n;
  for( double x=x_0; x<=x_1+0.0001; x+=dx ) {
    double f_val, df;
    f_val = f(x, &df);
    stream << x << " " << f_val << "   " << df << std::endl;
  }
  stream << std::endl;
}

void PairLCBOP::spline_init() {
  for( size_t N_conj_ij=0; N_conj_ij<2; N_conj_ij++ ) // N_conj_ij
  for( size_t N_ij=0; N_ij<4-1; N_ij++ )
  for( size_t N_ji=0; N_ji<4-1; N_ji++ ) {
    TF_conj_field &field = F_conj_field[N_ij][N_ji][N_conj_ij];
    field.f_00 = F_conj_data[N_ij  ][N_ji  ][N_conj_ij][0];
    field.f_01 = F_conj_data[N_ij  ][N_ji+1][N_conj_ij][0];
    field.f_10 = F_conj_data[N_ij+1][N_ji  ][N_conj_ij][0];
    field.f_11 = F_conj_data[N_ij+1][N_ji+1][N_conj_ij][0];

    field.f_x_00 =   F_conj_data[N_ij  ][N_ji  ][N_conj_ij][1] - field.f_10 + field.f_00;
    field.f_x_01 =   F_conj_data[N_ij  ][N_ji+1][N_conj_ij][1] - field.f_11 + field.f_01;
    field.f_x_10 = -(F_conj_data[N_ij+1][N_ji  ][N_conj_ij][1] - field.f_10 + field.f_00);
    field.f_x_11 = -(F_conj_data[N_ij+1][N_ji+1][N_conj_ij][1] - field.f_11 + field.f_01);

    field.f_y_00 =   F_conj_data[N_ij  ][N_ji  ][N_conj_ij][2] - field.f_01 + field.f_00;
    field.f_y_01 = -(F_conj_data[N_ij  ][N_ji+1][N_conj_ij][2] - field.f_01 + field.f_00);
    field.f_y_10 =   F_conj_data[N_ij+1][N_ji  ][N_conj_ij][2] - field.f_11 + field.f_10;
    field.f_y_11 = -(F_conj_data[N_ij+1][N_ji+1][N_conj_ij][2] - field.f_11 + field.f_10);
  }

  //some testing:
//  std::ofstream file( "test.txt" );
//    file << "gX:\n";
//    file  << gX[0] << " "
//          << gX[1] << " "
//          << gX[2] << " "
//          << gX[3] << " "
//          << gX[4] << " "
//          << gX[5] << std::endl;
//    file << "gC:\n";
//    for( int i=0; i<6; i++ )
//      file  << gC[i][0] << " "
//            << gC[i][1] << " "
//            << gC[i][2] << " "
//            << gC[i][3] << " "
//            << gC[i][4] << std::endl;
//    file << std::endl;
//
//    file << "gamma_1 = " << gamma_1 << std::endl;
//    file << "r_1 = " << r_1 << std::endl;
//    file << "r_2 = " << r_2 << std::endl;
//    file << "A = " << A << std::endl;
//    file << "B_1 = " << B_1 << std::endl;
//    file << "B_2 = " << B_2 << std::endl;
//    file << "alpha = " << alpha << std::endl;
//    file << "beta_1 = " << beta_1 << std::endl;
//    file << "beta_2 = " << beta_2 << std::endl;
//    file << "d = " << d << std::endl;
//    file << "C_1 = " << C_1 << std::endl;
//    file << "C_4 = " << C_4 << std::endl;
//    file << "C_6 = " << C_6 << std::endl;
//    file << "L = " << L << std::endl;
//    file << "kappa = " << kappa << std::endl;
//    file << "R_0 = " << R_0 << std::endl;
//    file << "R_1 = " << R_1 << std::endl;
//    file << "r_0 = " << r_0 << std::endl;
//    file << "r_1_LR = " << r_1_LR << std::endl;
//    file << "r_2_LR = " << r_2_LR << std::endl;
//    file << "v_1 = " << v_1 << std::endl;
//    file << "v_2 = " << v_2 << std::endl;
//    file << "eps_1 = " << eps_1 << std::endl;
//    file << "eps_2 = " << eps_2 << std::endl;
//    file << "lambda_1 = " << lambda_1 << std::endl;
//    file << "lambda_2 = " << lambda_2 << std::endl;
//    file << "eps = " << eps << std::endl;
//    file << "delta = " << delta << std::endl;
//    file << "r_2_sq = " << r_2_sq << std::endl;
//    file << std::endl;
//
//
//    file << "gSpline:" << std::endl;
//    double x_1 = 1, x_0 = -1;
//    int n=1000;
//    double dx = (x_1-x_0)/n;
//    for( double x=x_0; x<=x_1+0.0001; x+=dx ) {
//      double g, dg;
//      g = gSpline(x, &dg);
//      file << x << " " << g << " " << dg << std::endl;
//    }
//    file << std::endl;
//
//  file << "hSpline:" << std::endl;
//  double x_1 = 1, x_0 = -1;
//  int n=1000;
//  double dx = (x_1-x_0)/n;
//  for( double x=x_0; x<=x_1+0.0001; x+=dx ) {
//    double h, dh;
//    h = hSpline(x, &dh);
//    file << x << " " << h << " " << dh << std::endl;
//  }
//  file << std::endl;
//
//
//  file << "f_c:" << std::endl;
//  double x_1 = 4, x_0 = 0;
//  int n=1000;
//  double dx = (x_1-x_0)/n;
//  for( double x=x_0; x<=x_1+0.0001; x+=dx ) {
//    double f, df;
//    f = f_c(x, r_1, r_2, &df);
//    file << x << " " << f << " " << df << std::endl;
//  }
//  file << std::endl;

//  file << "F_conj_data\n";
//  for (int k = 0; k < 2; k++) { // 2 values of N_ij_conj
//    for (int l = 0; l < 3; l++) { // 3 types of data: f, dfdx, dfdy
//      for (int i = 0; i < 4; i++) { // 4x4 matrix
//        file
//          << F_conj_data[i][0][k][l] << " "
//          << F_conj_data[i][1][k][l] << " "
//          << F_conj_data[i][2][k][l] << " "
//          << F_conj_data[i][3][k][l] << std::endl;
//      }
//    file << std::endl;
//    }
//  }
//
//
//  file << "F_conj_0 ";
//  double dummy;
//  for( double y=0; y<=3.0+0.0001; y+=0.1 )
//    file << y << " ";
//  file << std::endl;
//  for( double x=0; x<=3.0+0.0001; x+=0.1 ){
//    file << x << " ";
//    for( double y=0; y<=3.0+0.0001; y+=0.1 )
//      file << F_conj( x, y, 0, &dummy, &dummy, &dummy ) << " ";
//    file << std::endl;
//  }
//
//  file << "dF0_dx ";
//  for( double y=0; y<=3.0+0.0001; y+=0.1 )
//    file << y << " ";
//  file << std::endl;
//  for( double x=0; x<=3.0+0.0001; x+=0.1 ){
//    file << x << " ";
//    for( double y=0; y<=3.0+0.0001; y+=0.1 ) {
//      double dF_dx;
//      F_conj( x, y, 0, &dF_dx, &dummy, &dummy );
//      file << dF_dx << " ";
//    }
//    file << std::endl;
//  }
//
//
//
//  file << "F_conj_1 ";
//  for( double y=0; y<=3.0+0.0001; y+=0.1 )
//    file << y << " ";
//  file << std::endl;
//  for( double x=0; x<=3.0+0.0001; x+=0.1 ){
//    file << x << " ";
//    for( double y=0; y<=3.0+0.0001; y+=0.1 )
//      file << F_conj( x, y, 0, &dummy, &dummy, &dummy ) << " ";
//    file << std::endl;
//  }

}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairLCBOP::memory_usage()
{
  double bytes = 0.0;
  bytes += maxlocal * sizeof(int);
  bytes += maxlocal * sizeof(int *);

  for (int i = 0; i < comm->nthreads; i++)
    bytes += ipage[i].size();

  bytes += 3*maxlocal * sizeof(double);
  return bytes;
}
