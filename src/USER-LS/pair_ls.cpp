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
   Contributing authors: Alexey Lipnitskiy (BSU), Daniil Poletaev (BSU) 
------------------------------------------------------------------------- */

#include "pair_ls.h"

#include <cmath>

#include <cstring>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "update.h"

#include "tokenizer.h"
#include "potential_file_reader.h"

// #include <Accelerate/Accelerate.h>
// #include <lapacke.h>

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

PairLS::PairLS(LAMMPS *lmp) : Pair(lmp)
{
  std::cout << "!!!!! PairLS debug mode !!!!! " << "PairLS constructor started working"  << std::endl;

  shag_sp_fi  = nullptr;
  shag_sp_ro  = nullptr;
  shag_sp_emb = nullptr;
  shag_sp_f   = nullptr;

  R_sp_fi  = nullptr; 
  R_sp_ro  = nullptr;
  R_sp_emb = nullptr;
  R_sp_f   = nullptr;
  R_sp_g   = nullptr;
  
  a_sp_fi = nullptr;
  b_sp_fi = nullptr;
  c_sp_fi = nullptr;
  d_sp_fi = nullptr;

  a_sp_ro = nullptr;
  b_sp_ro = nullptr;
  c_sp_ro = nullptr;
  d_sp_ro = nullptr;

  a_sp_emb = nullptr;
  b_sp_emb = nullptr;
  c_sp_emb = nullptr;
  d_sp_emb = nullptr;

  a_sp_f3 = nullptr;
  b_sp_f3 = nullptr;
  c_sp_f3 = nullptr;
  d_sp_f3 = nullptr;

  a_sp_g3 = nullptr;
  b_sp_g3 = nullptr;
  c_sp_g3 = nullptr;
  d_sp_g3 = nullptr; 

  a_sp_f4 = nullptr;
  b_sp_f4 = nullptr;
  c_sp_f4 = nullptr;
  d_sp_f4 = nullptr;

  a_sp_g4 = nullptr;
  b_sp_g4 = nullptr;
  c_sp_g4 = nullptr;
  d_sp_g4 = nullptr;

  fip_rmin = nullptr;

  z_ion  = nullptr;
  c_ZBL  = nullptr;
  d_ZBL  = nullptr;
  zz_ZBL = nullptr;
  a_ZBL  = nullptr;
  e0_ZBL = nullptr;

  n_sp_fi  = nullptr;
  n_sp_ro  = nullptr;
  n_sp_emb = nullptr;
  n_sp_f   = nullptr;
  n_sp_g   = nullptr;

  n_sort = atom->ntypes;

  comm_forward = 1;
  comm_reverse = 1;

  periodic[0] = domain->xperiodic;
  periodic[1] = domain->yperiodic;
  periodic[2] = domain->zperiodic;
  
  std::cout << "!!!!! PairLS debug mode !!!!! " << "PairLS constructor end working"  << std::endl;


}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairLS::~PairLS()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(shag_sp_fi);
    memory->destroy(shag_sp_ro);
    memory->destroy(shag_sp_emb);
    memory->destroy(shag_sp_f);

    memory->destroy(R_sp_fi);
    memory->destroy(R_sp_ro);
    memory->destroy(R_sp_emb);
    memory->destroy(R_sp_f);
    memory->destroy(R_sp_g);

    memory->destroy(a_sp_fi);
    memory->destroy(b_sp_fi);
    memory->destroy(c_sp_fi);
    memory->destroy(d_sp_fi);

    memory->destroy(a_sp_ro);
    memory->destroy(b_sp_ro);
    memory->destroy(c_sp_ro);
    memory->destroy(d_sp_ro);

    memory->destroy(a_sp_emb);
    memory->destroy(b_sp_emb);
    memory->destroy(c_sp_emb);
    memory->destroy(d_sp_emb);

    memory->destroy(a_sp_f3);
    memory->destroy(b_sp_f3);
    memory->destroy(c_sp_f3);
    memory->destroy(d_sp_f3);

    memory->destroy(a_sp_g3);
    memory->destroy(b_sp_g3);
    memory->destroy(c_sp_g3);
    memory->destroy(d_sp_g3);

    memory->destroy(a_sp_f4);
    memory->destroy(b_sp_f4);
    memory->destroy(c_sp_f4);
    memory->destroy(d_sp_f4);

    memory->destroy(a_sp_g4);
    memory->destroy(b_sp_g4);
    memory->destroy(c_sp_g4);
    memory->destroy(d_sp_g4);

    memory->destroy(fip_rmin);

    memory->destroy(z_ion);
    memory->destroy(c_ZBL);
    memory->destroy(d_ZBL);
    memory->destroy(zz_ZBL);
    memory->destroy(a_ZBL);
    memory->destroy(e0_ZBL);

    memory->destroy(Rmin_fi_ZBL);
    memory->destroy(c_fi_ZBL);

    memory->destroy(n_sp_fi);
    memory->destroy(n_sp_ro);
    memory->destroy(n_sp_emb);
    memory->destroy(n_sp_f);
    memory->destroy(n_sp_g);
  }


}

/* ---------------------------------------------------------------------- */

void PairLS::compute(int eflag, int vflag)
{
  int i,j,ii,jj,m,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,p,rhoip,rhojp,z2,z2p,recip,phip,psip,phi;
  double *coeff;
  int *ilist,*jlist,*numneigh,**firstneigh;

  // ev_init(eflag,vflag);


  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double sizex, sizey, sizez;
  int n_at = atom->natoms;

  memory->create(if_true_i,,"PairLS:if_true_i");

  sizex = domain->xprd; // global x box dimension
  sizey = domain->yprd; // global y box dimension
  sizez = domain->zprd; // global z box dimension

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLS::allocate()
{
  int n = atom->ntypes;

  std::cout << "!!!!! PairLS debug mode !!!!! " << " Start allocating memory"  << std::endl;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(shag_sp_fi,mi,mi,"PairLS:shag_sp_fi");
  memory->create(shag_sp_ro,mi,mi,"PairLS:shag_sp_ro");
  memory->create(shag_sp_emb,mi,"PairLS:shag_sp_emb");
  memory->create(shag_sp_f,mi,mi,"PairLS:shag_sp_f");

  memory->create(R_sp_fi,mfi,mi,mi,"PairLS:R_sp_fi");
  memory->create(R_sp_ro,mfi,mi,mi,"PairLS:R_sp_ro");
  memory->create(R_sp_emb,memb,mi,"PairLS:R_sp_emb");
  memory->create(R_sp_f,mf,mi,mi,"PairLS:R_sp_f");
  memory->create(R_sp_g,mg,"PairLS:R_sp_g");

  memory->create(a_sp_fi,mfi,mi,mi,"PairLS:a_sp_fi");
  memory->create(b_sp_fi,mfi,mi,mi,"PairLS:b_sp_fi");
  memory->create(c_sp_fi,mfi,mi,mi,"PairLS:c_sp_fi");
  memory->create(d_sp_fi,mfi,mi,mi,"PairLS:d_sp_fi");

  memory->create(a_sp_ro,mro,mi,mi,"PairLS:a_sp_ro");
  memory->create(b_sp_ro,mro,mi,mi,"PairLS:b_sp_ro");
  memory->create(c_sp_ro,mro,mi,mi,"PairLS:c_sp_ro");
  memory->create(d_sp_ro,mro,mi,mi,"PairLS:d_sp_ro");

  memory->create(a_sp_emb,memb,mi,"PairLS:a_sp_emb");
  memory->create(b_sp_emb,memb,mi,"PairLS:b_sp_emb");
  memory->create(c_sp_emb,memb,mi,"PairLS:c_sp_emb");
  memory->create(d_sp_emb,memb,mi,"PairLS:d_sp_emb");

  memory->create(a_sp_f3,mf,mf3,mi,mi,"PairLS:a_sp_f3");
  memory->create(b_sp_f3,mf,mf3,mi,mi,"PairLS:b_sp_f3");
  memory->create(c_sp_f3,mf,mf3,mi,mi,"PairLS:c_sp_f3");
  memory->create(d_sp_f3,mf,mf3,mi,mi,"PairLS:d_sp_f3");

  memory->create(a_sp_g3,mg,mf3,mf3,mi,"PairLS:a_sp_g3");
  memory->create(b_sp_g3,mg,mf3,mf3,mi,"PairLS:b_sp_g3");
  memory->create(c_sp_g3,mg,mf3,mf3,mi,"PairLS:c_sp_g3");
  memory->create(d_sp_g3,mg,mf3,mf3,mi,"PairLS:d_sp_g3");

  memory->create(a_sp_f4,mf,mi,mi,"PairLS:a_sp_f4");
  memory->create(b_sp_f4,mf,mi,mi,"PairLS:b_sp_f4");
  memory->create(c_sp_f4,mf,mi,mi,"PairLS:c_sp_f4");
  memory->create(d_sp_f4,mf,mi,mi,"PairLS:d_sp_f4");

  memory->create(a_sp_g4,mi,mi,"PairLS:a_sp_g4");
  memory->create(b_sp_g4,mi,mi,"PairLS:b_sp_g4");
  memory->create(c_sp_g4,mi,mi,"PairLS:c_sp_g4");
  memory->create(d_sp_g4,mi,mi,"PairLS:d_sp_g4");

  memory->create(fip_rmin,mi,mi,"PairLS:fip_rmin");

  memory->create(z_ion,mi,"PairLS:z_ion");
  memory->create(c_ZBL,4,"PairLS:c_ZBL");
  // memory->create(c_ZBL,4,mi,mi,"PairLS:c_ZBL");
  memory->create(d_ZBL,4,"PairLS:d_ZBL");
  // memory->create(d_ZBL,4,mi,mi,"PairLS:d_ZBL");
  memory->create(zz_ZBL,mi,mi,"PairLS:zz_ZBL");
  memory->create(a_ZBL,mi,mi,"PairLS:a_ZBL");
  memory->create(e0_ZBL,mi,mi,"PairLS:e0_ZBL");

  memory->create(Rmin_fi_ZBL,mi,mi,"PairLS:Rmin_fi_ZBL");
  memory->create(c_fi_ZBL,6,mi,mi,"PairLS:c_fi_ZBL");

  memory->create(n_sp_fi,mi,mi,"PairLS:n_sp_fi");
  memory->create(n_sp_ro,mi,mi,"PairLS:n_sp_ro");
  memory->create(n_sp_emb,mi,mi,"PairLS:n_sp_emb");
  memory->create(n_sp_f,mi,mi,"PairLS:n_sp_f");
  memory->create(n_sp_g,mi,mi,"PairLS:n_sp_g");
  memory->create(n_f3,mi,"PairLS:n_f3");

  allocated = 1;
  std::cout << "!!!!! PairLS debug mode !!!!! " << " End allocating memory"  << std::endl;

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLS::settings(int narg, char **/*arg*/)
{
  if (narg > 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
  Pair coeffs for each pair of atoms are the names of files with potential functions.
  They should be written in a row in the following order (example for 3-component system):
  pair_coeff pot_1 pot_2 pot_3 pot_1_2 pot_1_3 pot_2_3

  Here pot_i are the names of files containing potential functions for one sort of atom i,
  pot_i_j are the names of files containing cross potential functions for pair of atoms i and j.
   
  The total number of potential files should be equal to N(N+1)/2 where N is a number of the atom types in the system
------------------------------------------------------------------------- */

void PairLS::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  // int n_sort = atom->ntypes;
  int n_pot = n_sort*(n_sort+1)/2;
  double r_pot, w;
  int i, j;
  char name_32;

  std::cout << "!!!!! PairLS debug mode !!!!! " << " Number of atomtypes is " << n_sort << std::endl;

  // spelling that the number of arguments (files with potentials) is correspond to the number of atom types in the system
  if (narg != n_pot) error->all(FLERR,"Incorrect number of args for pair coefficients");

  // Start reading monoatomic potentials

  // read_pot_ls.f (14-18):
  // do i=1,n_sort
  //   name_32=name_pot_is(i)
  //   call wr_pot_ls_is(-1,name_32,i,info_pot)
  //   call par2pot_is(i)
  // enddo
  std::cout << "!!!!! PairLS debug mode !!!!! " << " Number of potentials is " << n_pot << std::endl;

  for (i = 1; i <= (n_sort); i++)
    {
      std::cout << "!!!!! PairLS debug mode !!!!! " << " Start reading potential for atom " << i << std::endl;
      // name_32 = arg[i];
      // r_pot_ls_is(name_32, i);
      std::cout << "!!!!! PairLS debug mode !!!!! " << " The arg with potential name is " << arg[i-1] << std::endl;
      r_pot_ls_is(arg[i-1], i);
      std::cout << "!!!!! PairLS debug mode !!!!! " << " End reading potential for atom " << i << std::endl;
      par2pot_is(i);
      std::cout << "!!!!! PairLS debug mode !!!!! " << " End parametrizing potential for atom " << i << std::endl;
      setflag[i][i] = 1;
    }

  // read_pot_ls.f (20-25):
	// w=R_sp_fi(n_sp_fi,1,1)
	// do i=1,n_sort
	// 	if(R_sp_fi(n_sp_fi,i,i) > w) w=R_sp_fi(n_sp_fi,i,i)
	// enddo
	// Rc_fi=w
	// r_pot=Rc_fi

	w = R_sp_fi[n_sp_fi[1][1]-1][1][1];
	for (i = 1; i <= n_sort; i++)
    {
      if (R_sp_fi[n_sp_fi[i][i]-1][i][i] > w) w = R_sp_fi[n_sp_f[i][i]-1][i][i];
    }
	Rc_fi = w;
	r_pot = Rc_fi;

  // read_pot_ls.f (27-32):
	// w=R_sp_f(n_sp_f,1,1)
	// do i=1,n_sort
	// 	if(R_sp_f(n_sp_f,i,i) > w) w=R_sp_f(n_sp_f,i,i)
	// enddo
	// Rc_f=w

	w = R_sp_f[n_sp_f[1][1]-1][1][1];
	for (i = 1; i <= n_sort; i++)
    {
      if (R_sp_f[n_sp_f[i][i]-1][i][i] > w) w = R_sp_f[n_sp_f[i][i]-1][i][i];
    }
	Rc_f = w;

  // End reading monoatomic potentials

  // Start reading cross-potentials

  // read_pot_ls.f (35-44):
	// if(n_sort>1) then
	// 	do i=1,n_sort-1
	// 		do j=i+1,n_sort
	// 			name_32=name_pot_is1_is2(i,j)
	// 			call wr_pot_ls_is1_is2(-1,name_32,i,j,info_pot)
	// 			call par2pot_is1_is2(i,j)
	// 			call par2pot_is1_is2(j,i)
	// 		enddo
	// 	enddo
	// endif


  if (n_sort > 1)
    {
      int ij = n_sort+1;
      for (i = 1; i <= n_sort-1; i++)
      {
        for (j = i + 1; j <= n_sort; j++)
        {
          // name_32 = arg[ij];
          // r_pot_ls_is1_is2(name_32, i, j);
          std::cout << "!!!!! PairLS debug mode !!!!! " << " Start reading cross potential for atoms " << i << " and " << j << std::endl;
          std::cout << "!!!!! PairLS debug mode !!!!! " << " The arg with potential name is " << arg[ij-1] << std::endl;
          r_pot_ls_is1_is2(arg[ij-1], i, j);
          std::cout << "!!!!! PairLS debug mode !!!!! " << " End reading cross potential for atoms " << i << " and " << j << std::endl;
          par2pot_is1_is2(i,j);
          std::cout << "!!!!! PairLS debug mode !!!!! " << " End parametrizing cross potential for atoms " << i << " and " << j << std::endl;
          par2pot_is1_is2(j,i);          
          std::cout << "!!!!! PairLS debug mode !!!!! " << " End parametrizing cross potential for atoms " << j << " and " << i << std::endl;
          setflag[i][j] = 1;
          ij++;
        }
      }
    }








}

 

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLS::init_style()
{
  // convert read-in file(s) to arrays and spline them

  // file2array();
  // array2spline();

  neighbor->request(this,instance_me);
  // embedstep = -1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLS::init_one(int i, int j)
{
  // single global cutoff = max of cut from all files read in
  // for funcfl could be multiple files
  // for setfl or fs, just one file

  // if (setflag[i][j] == 0) scale[i][j] = 1.0;
  // scale[j][i] = scale[i][j];

  // if (funcfl) {
  //   cutmax = 0.0;
  //   for (int m = 0; m < nfuncfl; m++)
  //     cutmax = MAX(cutmax,funcfl[m].cut);
  // } else if (setfl) cutmax = setfl->cut;
  // else if (fs) cutmax = fs->cut;

  // cutforcesq = cutmax*cutmax;

  // return cutmax;
  return 3.14;
}


// Specific functions for this pair style

void PairLS::r_pot_ls_is(char *name_32, int is)
{
  // FILE *potfile;

  // potfile = fopen(name_32, "r")

  // char info_pot;
  // getline(potfile, info_pot);
  // getline(potfile, if_g3_pot);
  // getline(potfile, if_g4_pot);
  std::cout << "!!!!! PairLS debug mode !!!!! " << " I am in a void PairLS::r_pot_ls_is" << std::endl;
  std::cout << "!!!!! PairLS debug mode !!!!! " << " name_32 is " << name_32 << std::endl;

  if(comm->me == 0) {
  //  {
    std::cout << "!!!!! PairLS debug mode !!!!! " << " if(comm->me == 0) condition is valid " << std::endl;
    PotentialFileReader reader(lmp, name_32, "ls");
    std::cout << "!!!!! PairLS debug mode !!!!! " << " PotentialFileReader instance is created" << std::endl;

    try {

      // wr_pot_ls_is.f (75):
      // read(23,'(A128)') info_pot

      reader.skip_line();  // the unnessessary information in the 1st line
      std::cout << " The 1st line is skipped" << std::endl;

      // wr_pot_ls_is.f (76-78):
      // read(23,*) if_g3_pot
      // read(23,*) if_g4_pot
      // read(23,*) if_gp0_pot(is)


      if_g3_pot = strcmp(reader.next_string().c_str(),".true.") == 0;
      std::cout << " if_g3_pot = " << if_g3_pot << std::endl;
      if_g4_pot = strcmp(reader.next_string().c_str(),".true.") == 0;
      std::cout << " if_g4_pot = " << if_g4_pot << std::endl;
      reader.skip_line(); // the unnessessary information in the 4th line
      std::cout << " The 4th line is skipped" << std::endl;

      // wr_pot_ls_is.f (80-83):
      // read(23,*) n_sp_fi
      // do i=1,n_sp_fi
      // read(23,*) R_sp_fi(i,is,is),a_sp_fi(i,is,is)
      // enddo

      n_sp_fi[is][is] = reader.next_int();
      std::cout << " n_sp_fi[" << is << "][" << is << "] = " << n_sp_fi[is][is] << std::endl;
      for (int n = 0; n < n_sp_fi[is][is]; n++)
      {
        ValueTokenizer sp_fi_values = reader.next_values(2);
        R_sp_fi[n][is][is] = sp_fi_values.next_double();
        a_sp_fi[n][is][is] = sp_fi_values.next_double();
        std::cout << " R_sp_fi[" << n << "][" << is << "][" << is << "] = " << R_sp_fi[n][is][is] << "   ";
        std::cout << " a_sp_fi[" << n << "][" << is << "][" << is << "] = " << a_sp_fi[n][is][is] << std::endl;
      }

      // wr_pot_ls_is.f (84-86):
      // read(23,*) fip_Rmin(is,is)
      // read(23,*) Rmin_fi_ZBL(is,is)
      // read(23,*) e0_ZBL(is,is)

      fip_rmin[is][is] = reader.next_double();
      std::cout << " fip_rmin[" << is << "][" << is << "] = " << fip_rmin[is][is] << std::endl;
      Rmin_fi_ZBL[is][is] = reader.next_double();
      std::cout << " Rmin_fi_ZBL[" << is << "][" << is << "] = " << Rmin_fi_ZBL[is][is] << std::endl;
      e0_ZBL[is][is] = reader.next_double();
      std::cout << " e0_ZBL[" << is << "][" << is << "] = " << e0_ZBL[is][is] << std::endl;


      // wr_pot_ls_is.f (88-91):
      // read(23,*) n_sp_ro
      // do i=1,n_sp_ro
      // read(23,*) R_sp_ro(i,is,is),a_sp_ro(i,is,is)
      // enddo

      n_sp_ro[is][is] = reader.next_int();
      std::cout << " n_sp_ro[" << is << "][" << is << "] = " << n_sp_ro[is][is] << std::endl;
      for (int n = 0; n < n_sp_ro[is][is]; n++)
      {
        ValueTokenizer sp_ro_values = reader.next_values(2);
        R_sp_ro[n][is][is] = sp_ro_values.next_double();
        a_sp_ro[n][is][is] = sp_ro_values.next_double();
        std::cout << " R_sp_ro[" << n << "][" << is << "][" << is << "] = " << R_sp_ro[n][is][is] << "   ";
        std::cout << " a_sp_ro[" << n << "][" << is << "][" << is << "] = " << a_sp_ro[n][is][is] << std::endl;
      }    

      // wr_pot_ls_is.f (92-95):
      // read(23,*) n_sp_emb
      // do i=1,n_sp_emb
      // read(23,*) R_sp_emb(i,is),a_sp_emb(i,is)
      // enddo

      n_sp_emb[is][is] = reader.next_int();
      std::cout << " n_sp_emb[" << is << "][" << is << "] = " << n_sp_emb[is][is] << std::endl;
      for (int n = 0; n < n_sp_emb[is][is]; n++)
      {
        ValueTokenizer sp_emb_values = reader.next_values(2);
        R_sp_emb[n][is] = sp_emb_values.next_double();
        a_sp_emb[n][is] = sp_emb_values.next_double();
        std::cout << " R_sp_emb[" << n << "][" << is << "] = " << R_sp_emb[n][is] << "   ";
        std::cout << " a_sp_emb[" << n << "][" << is << "] = " << a_sp_emb[n][is] << std::endl;
      }         

      // wr_pot_ls_is.f (97-101):
      // read(23,*) n_sp_f,n_f3(is)
      // do i=1,n_sp_f
      // read(23,*) R_sp_f(i,is,is),(a_sp_f3(i,i1,is,is),i1=1,n_f3(is))
      // enddo

      ValueTokenizer values = reader.next_values(2);
      n_sp_f[is][is] = values.next_int();
      std::cout << " n_sp_f[" << is << "][" << is << "] = " << n_sp_f[is][is] << std::endl;
      n_f3[is] = values.next_int();
      std::cout << " n_f3[" << is << "] = " << n_f3[is] << std::endl;
      for (int n = 0; n < n_sp_f[is][is]; n++)
      {
        ValueTokenizer sp_f_values = reader.next_values(n_f3[is]+1);
        R_sp_f[n][is][is] = sp_f_values.next_double();
        std::cout << " R_sp_f[" << n << "][" << is << "][" << is << "] = " << R_sp_f[n][is][is] << "   ";
        for (int n1 = 0; n1 < n_f3[is]; n1++)
        {
          a_sp_f3[n][n1][is][is] = sp_f_values.next_double();
          std::cout << " a_sp_f[" << n << "][" << n1 << "]["<< is << "][" << is << "] = " << a_sp_f3[n][n1][is][is] << "  ";
        }
        std::cout << std::endl;

      }

      // wr_pot_ls_is.f (102-108):
      // read(23,*) n_sp_g
      // read(23,*) (R_sp_g(i),i=1,n_sp_g)
      // do i1=1,n_f3(is)
      // do i2=1,i1
      // read(23,*) (a_sp_g3(i,i1,i2,is),i=1,n_sp_g)
      // enddo
      // enddo

      n_sp_g[is][is] = reader.next_int();
      std::cout << " n_sp_g[" << is << "][" << is << "] = " << n_sp_g[is][is] << "   ";
      values = reader.next_values(n_sp_g[is][is]);
      for (int n = 0; n < n_sp_g[is][is]; n++)
      {
        R_sp_g[n] = values.next_double();  // R_sp_g actually are the values of cos(ijk) from -1 (180 grad) to 1 (0 grad)
        std::cout << " R_sp_g[" << n << "] = " << R_sp_g[n]<< "   ";
      }
      std::cout << std::endl;
      
      for (int n = 0; n < n_f3[is]; n++)
      {
        for (int n1 = 0; n1 <= n; n1++)
        {
          ValueTokenizer sp_g_values = reader.next_values(n_sp_g[is][is]);
          for (int n2 = 0; n2 < n_sp_g[is][is]; n2++)
          {
            a_sp_g3[n2][n][n1][is] = sp_g_values.next_double();
            std::cout << " a_sp_f[" << n2 << "][" << n << "]["<< n1 << "][" << is << "] = " << a_sp_g3[n2][n][n1][is] << "  ";
          }
          std::cout << std::endl;
        }
      }

      // wr_pot_ls_is.f (120-126):
      // read(23,*) z_ion(is)
      // do i=1,4
      // read(23,*) c_ZBL(i)
      // enddo
      // do i=1,4
      // read(23,*) d_ZBL(i)
      // enddo

      z_ion[is] = reader.next_double();
      std::cout << " z_ion[" << is << "] = " << z_ion[is]<< std::endl;
      for (int n = 0; n < 4; n++)
      {
        c_ZBL[n] = reader.next_double();
        std::cout << " c_ZBL[" << n << "] = " << c_ZBL[n]<< std::endl;
      }

      for (int n = 0; n < 4; n++)
      {
        // d_ZBL[n][is][is] = reader.next_double();
        d_ZBL[n] = reader.next_double();
        // std::cout << " d_ZBL[" << n << "]["<< is << "][" << is << "] = " << d_ZBL[n][is][is]<< std::endl;
        std::cout << " d_ZBL[" << n <<  "] = " << d_ZBL[n]<< std::endl;
      }
    }
    catch (TokenizerException &e) {
      error->one(FLERR, e.what());
    }
  }
}



void PairLS::r_pot_ls_is1_is2(char *name_32, int is1, int is2)
{
  // FILE *potfile;

  // potfile = fopen(name_32, "r")

  // char info_pot;
  // getline(potfile, info_pot);
  // getline(potfile, if_g3_pot);
  // getline(potfile, if_g4_pot);
  std::cout << "!!!!! PairLS debug mode !!!!! " << " I am in a void PairLS::r_pot_ls_is1_is2" << std::endl;
  std::cout << "!!!!! PairLS debug mode !!!!! " << " name_32 is " << name_32 << std::endl;

  if(comm->me == 0) {
  //  {
    std::cout << "!!!!! PairLS debug mode !!!!! " << " if(comm->me == 0) condition is valid " << std::endl;
    PotentialFileReader reader(lmp, name_32, "ls");
    std::cout << "!!!!! PairLS debug mode !!!!! " << " PotentialFileReader instance is created" << std::endl;

    try {

      // wr_pot_ls_is1_is2.f (54):
      // read(23,'(A128)') info_pot

      reader.skip_line();  // the unnessessary information in the 1st line
      std::cout << " The 1st line is skipped" << std::endl;

      // wr_pot_ls_is1_is2.f (57-60):
      // read(23,*) n_sp_fi
      // do i=1,n_sp_fi
      // read(23,*) R_sp_fi(i,is1,is2),a_sp_fi(i,is1,is2)
      // enddo

      n_sp_fi[is1][is2] = reader.next_int();
      std::cout << " n_sp_fi[" << is1 << "][" << is2 << "] = " << n_sp_fi[is1][is2] << std::endl;
      for (int n = 0; n < n_sp_fi[is1][is2]; n++)
      {
        ValueTokenizer sp_fi_values = reader.next_values(2);
        R_sp_fi[n][is1][is2] = sp_fi_values.next_double();
        a_sp_fi[n][is1][is2] = sp_fi_values.next_double();
        std::cout << " R_sp_fi[" << n << "][" << is1 << "][" << is2 << "] = " << R_sp_fi[n][is1][is2] << "   ";
        std::cout << " a_sp_fi[" << n << "][" << is1 << "][" << is2 << "] = " << a_sp_fi[n][is1][is2] << std::endl;
      }

      // wr_pot_ls_is1_is2.f (61-63):
      // read(23,*) fip_Rmin(is1,is2)
      // read(23,*) Rmin_fi_ZBL(is1,is2)
      // read(23,*) e0_ZBL(is1,is2)

      fip_rmin[is1][is2] = reader.next_double();
      std::cout << " fip_rmin[" << is1 << "][" << is2 << "] = " << fip_rmin[is1][is2] << std::endl;
      Rmin_fi_ZBL[is1][is2] = reader.next_double();
      std::cout << " Rmin_fi_ZBL[" << is1 << "][" << is2 << "] = " << Rmin_fi_ZBL[is1][is2] << std::endl;
      e0_ZBL[is1][is2] = reader.next_double();
      std::cout << " e0_ZBL[" << is1 << "][" << is2 << "] = " << e0_ZBL[is1][is2] << std::endl;

      // wr_pot_ls_is1_is2.f (65-73):
      // do i=1,n_sp_fi
      // R_sp_fi(i,is2,is1)=R_sp_fi(i,is1,is2)
      // enddo
      // do i=1,n_sp_fi
      // a_sp_fi(i,is2,is1)=a_sp_fi(i,is1,is2)
      // enddo
      // fip_Rmin(is2,is1)=fip_Rmin(is1,is2)
      // Rmin_fi_ZBL(is2,is1)=Rmin_fi_ZBL(is1,is2)
      // e0_ZBL(is2,is1)=e0_ZBL(is1,is2)

      n_sp_fi[is2][is1] = n_sp_fi[is1][is2];
      for (int n = 0; n < n_sp_fi[is1][is2]; n++)
      {
        R_sp_fi[n][is2][is1] = R_sp_fi[n][is1][is2];
        a_sp_fi[n][is2][is1] = a_sp_fi[n][is1][is2];
        std::cout << " R_sp_fi[" << n << "][" << is2 << "][" << is1 << "] = " << R_sp_fi[n][is2][is1] << "   ";
        std::cout << " a_sp_fi[" << n << "][" << is2 << "][" << is1 << "] = " << a_sp_fi[n][is2][is1] << std::endl;
      }

      fip_rmin[is2][is1] = fip_rmin[is1][is2];
      std::cout << " fip_rmin[" << is2 << "][" << is1 << "] = " << fip_rmin[is2][is1] << std::endl;
      Rmin_fi_ZBL[is2][is1] = Rmin_fi_ZBL[is1][is2];
      std::cout << " Rmin_fi_ZBL[" << is2 << "][" << is1 << "] = " << Rmin_fi_ZBL[is2][is1] << std::endl;
      e0_ZBL[is2][is1] = e0_ZBL[is1][is2];
      std::cout << " e0_ZBL[" << is2 << "][" << is1 << "] = " << e0_ZBL[is2][is1] << std::endl;

      // wr_pot_ls_is1_is2.f (76-82):
      // read(23,*) n_sp_ro
      // do i=1,n_sp_ro
      // read(23,*) R_sp_ro(i,is1,is2),a_sp_ro(i,is1,is2)
      // enddo
      // do i=1,n_sp_ro
      // read(23,*) R_sp_ro(i,is2,is1),a_sp_ro(i,is2,is1)
      // enddo

      n_sp_ro[is1][is2] = reader.next_int();
      n_sp_ro[is2][is1] = n_sp_ro[is1][is2];
      std::cout << " n_sp_ro[" << is1 << "][" << is2 << "] = " << n_sp_ro[is1][is2] << std::endl;
      for (int n = 0; n < n_sp_ro[is1][is2]; n++)
      {
        ValueTokenizer sp_ro_is1_values = reader.next_values(2);
        R_sp_ro[n][is1][is2] = sp_ro_is1_values.next_double();
        a_sp_ro[n][is1][is2] = sp_ro_is1_values.next_double();
        std::cout << " R_sp_ro[" << n << "][" << is1 << "][" << is2 << "] = " << R_sp_ro[n][is1][is2] << "   ";
        std::cout << " a_sp_ro[" << n << "][" << is1 << "][" << is2 << "] = " << a_sp_ro[n][is1][is2] << std::endl;
      }    

      for (int n = 0; n < n_sp_ro[is2][is1]; n++)
      {
        ValueTokenizer sp_ro_is2_values = reader.next_values(2);
        R_sp_ro[n][is2][is1] = sp_ro_is2_values.next_double();
        a_sp_ro[n][is2][is1] = sp_ro_is2_values.next_double();
        std::cout << " R_sp_ro[" << n << "][" << is2 << "][" << is1 << "] = " << R_sp_ro[n][is2][is1] << "   ";
        std::cout << " a_sp_ro[" << n << "][" << is2 << "][" << is1 << "] = " << a_sp_ro[n][is2][is1] << std::endl;
      }   

      // wr_pot_ls_is1_is2.f (85-88):
      // read(23,*) n_sp_f,n_f3(is1)
      // do i=1,n_sp_f
      // read(23,*) R_sp_f(i,is2,is1),(a_sp_f3(i,i1,is2,is1),i1=1,n_f3(is1))
      // enddo

      ValueTokenizer n_sp_f_is1_values = reader.next_values(2);
      n_sp_f[is2][is1] = n_sp_f_is1_values.next_int();
      std::cout << " n_sp_f[" << is2 << "][" << is1 << "] = " << n_sp_f[is2][is1] << std::endl;
      n_f3[is1] = n_sp_f_is1_values.next_int();
      std::cout << " n_f3[" << is1 << "] = " << n_f3[is1] << std::endl;
      for (int n = 0; n < n_sp_f[is2][is1]; n++)
      {
        ValueTokenizer sp_f_is1_values = reader.next_values(n_f3[is1]+1);
        R_sp_f[n][is2][is1] = sp_f_is1_values.next_double();
        std::cout << " R_sp_f[" << n << "][" << is2 << "][" << is1 << "] = " << R_sp_f[n][is2][is1] << "   ";
        for (int n1 = 0; n1 < n_f3[is1]; n1++)
        {
          a_sp_f3[n][n1][is2][is1] = sp_f_is1_values.next_double();
          std::cout << " a_sp_f[" << n << "][" << n1 << "]["<< is2 << "][" << is1 << "] = " << a_sp_f3[n][n1][is2][is1] << "  ";
        }
        std::cout << std::endl;

      }

      // wr_pot_ls_is1_is2.f (89-92):
      // read(23,*) n_sp_f,n_f3(is2)
      // do i=1,n_sp_f
      // read(23,*) R_sp_f(i,is1,is2),(a_sp_f3(i,i1,is1,is2),i1=1,n_f3(is2))
      // enddo

      ValueTokenizer n_sp_f_is2_values = reader.next_values(2);
      n_sp_f[is1][is2] = n_sp_f_is2_values.next_int();
      std::cout << " n_sp_f[" << is1 << "][" << is2 << "] = " << n_sp_f[is1][is2] << std::endl;
      n_f3[is2] = n_sp_f_is2_values.next_int();
      std::cout << " n_f3[" << is2 << "] = " << n_f3[is2] << std::endl;
      for (int n = 0; n < n_sp_f[is1][is2]; n++)
      {
        ValueTokenizer sp_f_is2_values = reader.next_values(n_f3[is2]+1);
        R_sp_f[n][is1][is2] = sp_f_is2_values.next_double();
        std::cout << " R_sp_f[" << n << "][" << is1 << "][" << is2 << "] = " << R_sp_f[n][is1][is2] << "   ";
        for (int n1 = 0; n1 < n_f3[is2]; n1++)
        {
          a_sp_f3[n][n1][is1][is2] = sp_f_is2_values.next_double();
          std::cout << " a_sp_f[" << n << "][" << n1 << "]["<< is1 << "][" << is2 << "] = " << a_sp_f3[n][n1][is1][is2] << "  ";
        }
        std::cout << std::endl;

      }

    }
    catch (TokenizerException &e) {
      error->one(FLERR, e.what());
    }
  }
}



void PairLS::par2pot_is(int is)
{
  // Start pot_ls_a_sp.h
  int n_sp;
  // double R_sp[mfi];
  // double a_sp[mfi], b_sp[mfi], c_sp[mfi], d_sp[mfi], e_sp[mfi];
  double *R_sp;
  double *a_sp, *b_sp, *c_sp, *d_sp, *e_sp; 
  // End pot_ls_a_sp.h

  int i, j, i1, i2, n;
  double p1, p2, pn;
  // double B6[6][1];
  double *B6;
  double r1, r2, f1, fp1, fpp1, f2, fp2, fpp2;
  std::cout << "!!!!! PairLS debug mode !!!!! " << " entering par2pot_is" << std::endl;
  
  memory->destroy(R_sp);
  memory->destroy(a_sp);
  memory->destroy(b_sp);
  memory->destroy(c_sp);
  memory->destroy(d_sp);
  memory->destroy(e_sp);
  memory->destroy(B6);
  // std::cout << "!!!!! PairLS debug mode !!!!! " << " destroyed memory for R_sp, a_sp, etc" << std::endl;

  memory->create(R_sp, mfi, "PairLS:a_sp_w");
  memory->create(a_sp, mfi, "PairLS:a_sp_w");
  memory->create(b_sp, mfi, "PairLS:b_sp_w");
  memory->create(c_sp, mfi, "PairLS:c_sp_w");
  memory->create(d_sp, mfi, "PairLS:d_sp_w");
  memory->create(e_sp, mfi, "PairLS:e_sp_w");
  memory->create(B6, 6, "PairLS:B6_w");

  // par2pot_is.f(15-18):
  //   zz_ZBL(is,is)=z_ion(is)*z_ion(is)*(3.795D0**2)
  //   a_ZBL(is,is)=0.8853D0
  //  :	*0.5291772083D0/(z_ion(is)**0.23D0 + z_ion(is)**0.23D0)
  // std::cout << "!!!!! PairLS debug mode !!!!! " << " created memory for R_sp, a_sp, etc" << std::endl;

  zz_ZBL[is][is] = z_ion[is]*z_ion[is]*pow(3.795,2);
  a_ZBL[is][is] = 0.8853*0.5291772083/(pow(z_ion[is],0.23) + pow(z_ion[is],0.23));

  // std::cout << "!!!!! PairLS debug mode !!!!! " << " zz_ZBL[is][is] and a_ZBL[is][is] assigned" << std::endl;

  // par2pot_is.f(21-28):
	// n_sp=n_sp_fi
	// do i=1,n_sp
	// R_sp(i)=R_sp_fi(i,is,is)
	// enddo
	// do i=1,n_sp-1
	// a_sp(i)=a_sp_fi(i,is,is)
	// enddo
	// a_sp(n_sp)=0.0D0

	n_sp = n_sp_fi[is][is];
	for (i = 0; i < n_sp; i++)
  {
    R_sp[i]=R_sp_fi[i][is][is];
  }
	for (i = 0; i < n_sp-1; i++)
  {
    a_sp[i]=a_sp_fi[i][is][is];
  }  
	a_sp[n_sp-1] = 0.0;

  // par2pot_is.f(29-36):
	// call SPL(n_sp, R_sp, a_sp, 1, fip_Rmin(is,is),0.0D0, b_sp,c_sp,d_sp)
		// do i=1,n_sp_fi
		// a_sp_fi(i,is,is)=a_sp(i)
		// b_sp_fi(i,is,is)=b_sp(i)
		// c_sp_fi(i,is,is)=c_sp(i)
		// d_sp_fi(i,is,is)=d_sp(i)
		// enddo
	// shag_sp_fi(is,is)=1.0D0/((R_sp_fi(n_sp,is,is)-R_sp_fi(1,is,is))/dfloat(n_sp-1))

	// SPL(n_sp, &R_sp, &a_sp, 1, fip_rmin[is][is], 0.0, &b_sp, &c_sp, &d_sp);
  // std::cout << "!!!!! PairLS debug mode !!!!! " << " now SPL function will be called" << std::endl;

	SPL(n_sp, R_sp, a_sp, 1, fip_rmin[is][is], 0.0, b_sp, c_sp, d_sp);

  // std::cout << "!!!!! PairLS debug mode !!!!! " << " SPL function worked well" << std::endl;

  std::cout << "!!!!!!!!!!!! Fi" << std::endl;
	for (i = 0; i < n_sp; i++)
  {
    a_sp_fi[i][is][is] = a_sp[i];
    b_sp_fi[i][is][is] = b_sp[i];
    c_sp_fi[i][is][is] = c_sp[i];
    d_sp_fi[i][is][is] = d_sp[i];
    std::cout << " R_sp_fi[" << i << "][" << is << "][" << is << "] = " << R_sp_fi[i][is][is] << "   ";
    std::cout << " a_sp_fi[" << i << "][" << is << "][" << is << "] = " << a_sp_fi[i][is][is] << "   ";    
    std::cout << " b_sp_fi[" << i << "][" << is << "][" << is << "] = " << b_sp_fi[i][is][is] << "   ";    
    std::cout << " c_sp_fi[" << i << "][" << is << "][" << is << "] = " << c_sp_fi[i][is][is] << "   ";    
    std::cout << " d_sp_fi[" << i << "][" << is << "][" << is << "] = " << d_sp_fi[i][is][is] << std::endl;    
  }
	
	shag_sp_fi[is][is] = 1.0/((R_sp_fi[n_sp-1][is][is]-R_sp_fi[0][is][is])/(n_sp-1));

  // par2pot_is.f(39-49):
  // c fi_ZBL
	// r1=Rmin_fi_ZBL(is,is)
	// f1=v_ZBL(r1,is,is)+e0_ZBL(is,is)
	// fp1=vp_ZBL(r1,is,is)
	// fpp1=vpp_ZBL(r1,is,is)
	//     r2=R_sp_fi(1,is,is)
	//     f2=a_sp_fi(1,is,is)
	//     fp2=b_sp_fi(1,is,is)
	//     fpp2=2.0D0*c_sp_fi(1,is,is)	
	// call smooth_zero_22 (B6,r1,r2,f1,fp1,fpp1,
  //    :                         f2,fp2,fpp2)
	// c_fi_ZBL(1:6,is,is)=B6(1:6)

	r1 = Rmin_fi_ZBL[is][is];
	f1 = v_ZBL(r1, is, is) + e0_ZBL[is][is];
	fp1 = vp_ZBL(r1, is, is);
	fpp1 = vpp_ZBL(r1, is, is);
  r2 = R_sp_fi[0][is][is];
  f2 = a_sp_fi[0][is][is];
  fp2 = b_sp_fi[0][is][is];
  fpp2 = 2.0*c_sp_fi[0][is][is];	
	
  // std::cout << "!!!!! PairLS debug mode !!!!! " << " now smooth_zero_22 function will be called" << std::endl;

  smooth_zero_22(B6, r1, r2, f1, fp1, fpp1, f2, fp2, fpp2);
  // smooth_zero_22(B6, &r1, &r2, &f1, &fp1, &fpp1, &f2, &fp2, &fpp2);
  // std::cout << "!!!!! PairLS debug mode !!!!! " << " smooth_zero_22 function worked well" << std::endl;

  std::cout << "!!!!!!!!!!!! c_fi_ZBL" << std::endl;
  std::cout << " r1 = " << r1 <<std::endl;
  std::cout << " r2 = " << r2 <<std::endl;
  std::cout << " f1 = " << f1 <<std::endl;
  std::cout << " fp1 = " << fp1 <<std::endl;
  std::cout << " fpp1 = " << fpp1 <<std::endl;
  std::cout << " f2 = " << f2 <<std::endl;
  std::cout << " fp2 = " << fp2 <<std::endl;
  std::cout << " fpp2 = " << fpp2 <<std::endl;  

  for (i = 0; i < 6; i++)
  {
    c_fi_ZBL[i][is][is] = B6[i];
    std::cout << c_fi_ZBL[i][is][is] << "  ";
  }
  std::cout << std::endl;
  // std::cout << "!!!!! PairLS debug mode !!!!! " << " for (i = 0; i < 6; i++) worked well" << std::endl;

  // par2pot_is.f(51-68):
  // ! ro
	// n_sp=n_sp_ro
	// do i=1,n_sp
	// R_sp(i)=R_sp_ro(i,is,is)
	// enddo
	// do i=1,n_sp
	// a_sp(i)=a_sp_ro(i,is,is)
	// enddo
	// p1=0.0D0
  // c	p1=(a_sp(2)-a_sp(1))/(R_sp(2)-R_sp(1))
	// call SPL(n_sp, R_sp, a_sp, 1, p1,0.0D0, b_sp,c_sp,d_sp)
	// do i=1,n_sp
	// a_sp_ro(i,is,is)=a_sp(i)
	// b_sp_ro(i,is,is)=b_sp(i)
	// c_sp_ro(i,is,is)=c_sp(i)
	// d_sp_ro(i,is,is)=d_sp(i)
	// enddo
	// shag_sp_ro(is,is)=1.0D0/((R_sp_ro(n_sp,is,is)-R_sp_ro(1,is,is))/dfloat(n_sp-1))

  n_sp = n_sp_ro[is][is];
  std::cout << " n_sp = " << n_sp_ro[is][is] << std::endl;
  for (i = 0; i < n_sp; i++)
  {
    std::cout << " i = " << i << " R_sp_ro[i][is][is] = " << R_sp_ro[i][is][is] << std::endl;
    R_sp[i] = R_sp_ro[i][is][is];
    a_sp[i] = a_sp_ro[i][is][is];
  }
  p1=0.0;

  // std::cout << "!!!!! PairLS debug mode !!!!! " << " for (i = 0; i < n_sp; i++) worked well" << std::endl;


	// SPL(n_sp, &R_sp, &a_sp, 1, p1, 0.0, &b_sp, &c_sp, &d_sp);
	SPL(n_sp, R_sp, a_sp, 1, p1, 0.0, b_sp, c_sp, d_sp);
  // std::cout << "!!!!! PairLS debug mode !!!!! " << " SPL(n_sp, R_sp, a_sp, 1, p1, 0.0, b_sp, c_sp, d_sp) worked well" << std::endl;

  std::cout << "!!!!!!!!!!!! Rho" << std::endl;
	for (i = 0; i < n_sp; i++)
  {  
    a_sp_ro[i][is][is] = a_sp[i];
    b_sp_ro[i][is][is] = b_sp[i];
    c_sp_ro[i][is][is] = c_sp[i];
    d_sp_ro[i][is][is] = d_sp[i];
    std::cout << " R_sp_ro[" << i << "][" << is << "][" << is << "] = " << R_sp_ro[i][is][is] << "   ";
    std::cout << " a_sp_ro[" << i << "][" << is << "][" << is << "] = " << a_sp_ro[i][is][is] << "   ";    
    std::cout << " b_sp_ro[" << i << "][" << is << "][" << is << "] = " << b_sp_ro[i][is][is] << "   ";    
    std::cout << " c_sp_ro[" << i << "][" << is << "][" << is << "] = " << c_sp_ro[i][is][is] << "   ";    
    std::cout << " d_sp_ro[" << i << "][" << is << "][" << is << "] = " << d_sp_ro[i][is][is] << std::endl;     
  }
	shag_sp_ro[is][is] = 1.0/((R_sp_ro[n_sp-1][is][is]-R_sp_ro[0][is][is])/(n_sp-1));

  // std::cout << "!!!!! PairLS debug mode !!!!! " << " for (i = 0; i < n_sp; i++) worked well" << std::endl;

  // par2pot_is.f(70-91):
  // ! emb
  // 	n_sp=n_sp_emb
  // 	n=n_sp_emb
  // 	do i=1,n_sp
  // 	R_sp(i)=R_sp_emb(i,is)
  // 	enddo
  // 	do i=1,n_sp
  // 	a_sp(i)=a_sp_emb(i,is)
  // 	enddo
  // 	a_sp(1)=0.0D0
  // 	p1=(a_sp(2)-a_sp(1))/(R_sp(2)-R_sp(1))
  // c	pn=(a_sp(n)-a_sp(n-1))/(R_sp(n)-R_sp(n-1))
  // 	pn=0.0D0
  // 	call SPL(n_sp, R_sp, a_sp, 1, p1,pn, b_sp,c_sp,d_sp)
  // 		do i=1,n_sp
  // 		a_sp_emb(i,is)=a_sp(i)
  // 		b_sp_emb(i,is)=b_sp(i)
  // 		c_sp_emb(i,is)=c_sp(i)
  // 		d_sp_emb(i,is)=d_sp(i)
  // 		enddo
  // 	shag_sp_emb(is)=1.0D0/((R_sp_emb(n_sp,is)-R_sp_emb(1,is))/dfloat(n_sp-1))

	n_sp = n_sp_emb[is][is];
	n = n_sp_emb[is][is];
  for (i = 0; i < n_sp; i++)
  {
    R_sp[i] = R_sp_emb[i][is];
    a_sp[i] = a_sp_emb[i][is];
  }  
	a_sp[0] = 0.0;

	p1 = (a_sp[1]-a_sp[0])/(R_sp[1]-R_sp[0]);
	pn = 0.0;
	// SPL(n_sp, &R_sp, &a_sp, 1, p1, pn, &b_sp, &c_sp, &d_sp);
	SPL(n_sp, R_sp, a_sp, 1, p1, pn, b_sp, c_sp, d_sp);
  // std::cout << "!!!!! PairLS debug mode !!!!! " << " SPL(n_sp, R_sp, a_sp, 1, p1, pn, b_sp, c_sp, d_sp) worked well" << std::endl;
  std::cout << "!!!!!!!!!!!! Emb" << std::endl;
	for (i = 0; i < n_sp; i++)
  {  
    a_sp_emb[i][is] = a_sp[i];
    b_sp_emb[i][is] = b_sp[i];
    c_sp_emb[i][is] = c_sp[i];
    d_sp_emb[i][is] = d_sp[i];
    std::cout << " R_sp_emb[" << i << "][" << is << "] = " << R_sp_emb[i][is] << "   ";
    std::cout << " a_sp_emb[" << i << "][" << is << "] = " << a_sp_emb[i][is] << "   ";    
    std::cout << " b_sp_emb[" << i << "][" << is << "] = " << b_sp_emb[i][is] << "   ";    
    std::cout << " c_sp_emb[" << i << "][" << is << "] = " << c_sp_emb[i][is] << "   ";    
    std::cout << " d_sp_emb[" << i << "][" << is << "] = " << d_sp_emb[i][is] << std::endl;     
  }

	shag_sp_emb[is] = 1.0/((R_sp_emb[n_sp-1][is] - R_sp_emb[0][is])/(n_sp - 1)); 

  // par2pot_is.f(97-115):
  // ! f3
  // 	n_sp=n_sp_f
  // 	do i=1,n_sp
  // 	R_sp(i)=R_sp_f(i,is,is)
  // 	enddo
  // 	do i1=1,n_f3(is)
  // 	    do i=1,n_sp
  // 	    a_sp(i)=a_sp_f3(i,i1,is,is)
  // 	    enddo
  // 	    p1=0.0D0
  // c	    p1=(a_sp(2)-a_sp(1))/(R_sp(2)-R_sp(1))
  // 	call SPL(n_sp, R_sp, a_sp, 1, p1,0.0D0, b_sp,c_sp,d_sp)
  // 		do i=1,n_sp
  // 		a_sp_f3(i,i1,is,is)=a_sp(i)
  // 		b_sp_f3(i,i1,is,is)=b_sp(i)
  // 		c_sp_f3(i,i1,is,is)=c_sp(i)
  // 		d_sp_f3(i,i1,is,is)=d_sp(i)
  // 		enddo
  // 	enddo
  // 	shag_sp_f(is,is)=1.0D0/((R_sp_f(n_sp,is,is)-R_sp_f(1,is,is))/dfloat(n_sp-1))

	n_sp = n_sp_f[is][is];
  for (i = 0; i < n_sp; i++)
  {
    R_sp[i] = R_sp_f[i][is][is];
  }

  std::cout << "!!!!!!!!!!!! f3" << std::endl;
  for (i1 = 0; i1 < n_f3[is]; i1++)
  {
    for (i = 0; i < n_sp; i++)
    {
      a_sp[i] = a_sp_f3[i][i1][is][is];
    }
    p1 = 0.0;
    // SPL(n_sp, &R_sp, &a_sp, 1, p1, 0.0, &b_sp, &c_sp, &d_sp);
    SPL(n_sp, R_sp, a_sp, 1, p1, 0.0, b_sp, c_sp, d_sp);
    // std::cout << "!!!!! PairLS debug mode !!!!! " << " SPL(n_sp, R_sp, a_sp, 1, p1, 0.0, b_sp, c_sp, d_sp) worked well" << std::endl;
    for (i = 0; i < n_sp; i++)
    { 
      a_sp_f3[i][i1][is][is] = a_sp[i];
      b_sp_f3[i][i1][is][is] = b_sp[i];
      c_sp_f3[i][i1][is][is] = c_sp[i];
      d_sp_f3[i][i1][is][is] = d_sp[i];
      std::cout << " R_sp_f3[" << i << "][" << is << "][" << is << "] = " << R_sp_f[i][is][is] << "   ";
      std::cout << " a_sp_f3[" << i << "][" << i1 << "][" << is << "][" << is << "] = " << a_sp_f3[i][i1][is][is] << "   ";    
      std::cout << " b_sp_f3[" << i << "][" << i1 << "][" << is << "][" << is << "] = " << b_sp_f3[i][i1][is][is] << "   ";    
      std::cout << " c_sp_f3[" << i << "][" << i1 << "][" << is << "][" << is << "] = " << c_sp_f3[i][i1][is][is] << "   ";    
      std::cout << " d_sp_f3[" << i << "][" << i1 << "][" << is << "][" << is << "] = " << d_sp_f3[i][i1][is][is] << std::endl;      
    }
  }
	shag_sp_f[is][is] = 1.0/((R_sp_f[n_sp-1][is][is]-R_sp_f[0][is][is])/(n_sp-1));

  // par2pot_is.f(117-150):
  // ! g3
  // 	n_sp=n_sp_g
  // 	do i=1,n_sp
  // 	R_sp(i)=R_sp_g(i)
  // 	enddo
  // c						    if_diag(is)=.false.
  // 	do i1=1,n_f3(is)
  // 	do i2=1,i1
  // 	    do i=1,n_sp
  // 	    a_sp(i)=a_sp_g3(i,i1,i2,is)
  // 	    enddo
  // c	    if(i2.NE.i1) then
  // c	    if(abs(a_sp(2))<0.00000001.AND.abs(a_sp(7))<0.00000001) if_diag(is)=.true.
  // c	    endif
  // 	p1=0.0D0
  // 	p2=0.0D0
  // 	    if(.NOT.if_gp0_pot(is)) then
  // 	    p1=(a_sp(2)-a_sp(1))/(R_sp(2)-R_sp(1))
  // 	    p2=(a_sp(n_sp)-a_sp(n_sp-1))/(R_sp(n_sp)-R_sp(n_sp-1))
  // 	    endif
  // 	call SPL(n_sp, R_sp, a_sp, 1, p1,p2, b_sp,c_sp,d_sp)
  // 		do i=1,n_sp
  // 		a_sp_g3(i,i1,i2,is)=a_sp(i)
  // 		a_sp_g3(i,i2,i1,is)=a_sp(i)
  // 		b_sp_g3(i,i1,i2,is)=b_sp(i)
  // 		b_sp_g3(i,i2,i1,is)=b_sp(i)
  // 		c_sp_g3(i,i1,i2,is)=c_sp(i)
  // 		c_sp_g3(i,i2,i1,is)=c_sp(i)
  // 		d_sp_g3(i,i1,i2,is)=d_sp(i)
  // 		d_sp_g3(i,i2,i1,is)=d_sp(i)
  // 		enddo
  // 	enddo
  // 	enddo
  // 	shag_sp_g=1.0D0/((R_sp_g(n_sp)-R_sp_g(1))/dfloat(n_sp-1))

  std::cout << "!!!!!!!!!!!! g3" << std::endl;
	n_sp = n_sp_g[is][is];
  for (i = 0; i < n_sp; i++)
  {
    R_sp[i] = R_sp_g[i];
  }
  for (i1 = 0; i1 < n_f3[is]; i1++)
  {
    for (i2 = 0; i2 <= i1; i2++)
    {
      for (i = 0; i < n_sp; i++)
      {
        a_sp[i] = a_sp_g3[i][i1][i2][is];
      }
      p1 = 0.0;
      p2 = 0.0;
      // SPL(n_sp, &R_sp, &a_sp, 1, p1, p2, &b_sp, &c_sp, &d_sp);
      SPL(n_sp, R_sp, a_sp, 1, p1, p2, b_sp, c_sp, d_sp);
      // std::cout << "!!!!! PairLS debug mode !!!!! " << " SPL(n_sp, R_sp, a_sp, 1, p1, p2, b_sp, c_sp, d_sp) worked well" << std::endl;
      for (i = 0; i < n_sp; i++)
      { 
        a_sp_g3[i][i1][i2][is] = a_sp[i];
        a_sp_g3[i][i2][i1][is] = a_sp[i];
        b_sp_g3[i][i1][i2][is] = b_sp[i];
        b_sp_g3[i][i2][i1][is] = b_sp[i];
        c_sp_g3[i][i1][i2][is] = c_sp[i];
        c_sp_g3[i][i2][i1][is] = c_sp[i];
        d_sp_g3[i][i1][i2][is] = d_sp[i];
        d_sp_g3[i][i2][i1][is] = d_sp[i];
        std::cout << " R_sp_g3[" << i <<  "] = " << R_sp_g[i] << "   ";
        std::cout << " a_sp_g3[" << i << "][" << i1 << "][" << i2 << "][" << is << "] = " << a_sp_g3[i][i1][i2][is] << "   ";    
        std::cout << " b_sp_g3[" << i << "][" << i1 << "][" << i2 << "][" << is << "] = " << b_sp_g3[i][i1][i2][is] << "   ";    
        std::cout << " c_sp_g3[" << i << "][" << i1 << "][" << i2 << "][" << is << "] = " << c_sp_g3[i][i1][i2][is] << "   ";    
        std::cout << " d_sp_g3[" << i << "][" << i1 << "][" << i2 << "][" << is << "] = " << d_sp_g3[i][i1][i2][is] << std::endl;   
        // 
        std::cout << " R_sp_g3[" << i <<  "] = " << R_sp_g[i] << "   ";
        std::cout << " a_sp_g3[" << i << "][" << i2 << "][" << i1 << "][" << is << "] = " << a_sp_g3[i][i2][i1][is] << "   ";    
        std::cout << " b_sp_g3[" << i << "][" << i2 << "][" << i1 << "][" << is << "] = " << b_sp_g3[i][i2][i1][is] << "   ";    
        std::cout << " c_sp_g3[" << i << "][" << i2 << "][" << i1 << "][" << is << "] = " << c_sp_g3[i][i2][i1][is] << "   ";    
        std::cout << " d_sp_g3[" << i << "][" << i2 << "][" << i1 << "][" << is << "] = " << d_sp_g3[i][i2][i1][is] << std::endl;     
      }
    }  
  }
	shag_sp_g = 1.0/((R_sp_g[n_sp-1]-R_sp_g[0])/(n_sp-1));

  memory->destroy(R_sp);
  memory->destroy(a_sp);
  memory->destroy(b_sp);
  memory->destroy(c_sp);
  memory->destroy(d_sp);
  memory->destroy(e_sp);
  memory->destroy(B6);
}



void PairLS::par2pot_is1_is2(int is1, int is2)
{
  // Start pot_ls_a_sp.h
  int n_sp;
  // double R_sp[mfi];
  // double a_sp[mfi], b_sp[mfi], c_sp[mfi], d_sp[mfi], e_sp[mfi];
  double *R_sp;
  double *a_sp, *b_sp, *c_sp, *d_sp, *e_sp; 
  // End pot_ls_a_sp.h

  int i, j, i1, i2, n;
  double p1, p2, pn;
  // double B6[6][1];
  double *B6;
  double r1, r2, f1, fp1, fpp1, f2, fp2, fpp2;
  std::cout << "!!!!! PairLS debug mode !!!!! " << " entering par2pot_is1_is2" << std::endl;
  
  memory->destroy(R_sp);
  memory->destroy(a_sp);
  memory->destroy(b_sp);
  memory->destroy(c_sp);
  memory->destroy(d_sp);
  memory->destroy(e_sp);
  memory->destroy(B6);
  // std::cout << "!!!!! PairLS debug mode !!!!! " << " destroyed memory for R_sp, a_sp, etc" << std::endl;

  memory->create(R_sp, mfi, "PairLS:a_sp_w");
  memory->create(a_sp, mfi, "PairLS:a_sp_w");
  memory->create(b_sp, mfi, "PairLS:b_sp_w");
  memory->create(c_sp, mfi, "PairLS:c_sp_w");
  memory->create(d_sp, mfi, "PairLS:d_sp_w");
  memory->create(e_sp, mfi, "PairLS:e_sp_w");
  memory->create(B6, 6, "PairLS:B6_w");  

  // std::cout << "!!!!! PairLS debug mode !!!!! " << " created memory for R_sp, a_sp, etc" << std::endl;

  zz_ZBL[is1][is2] = z_ion[is1]*z_ion[is2]*pow(3.795,2);
  a_ZBL[is1][is2] = 0.8853*0.5291772083/(pow(z_ion[is1],0.23) + pow(z_ion[is2],0.23));

  // std::cout << "!!!!! PairLS debug mode !!!!! " << " zz_ZBL[is][is] and a_ZBL[is][is] assigned" << std::endl;

	n_sp = n_sp_fi[is1][is2];
	for (i = 0; i < n_sp; i++)
  {
    R_sp[i]=R_sp_fi[i][is1][is2];
  }
	for (i = 0; i < n_sp-1; i++)
  {
    a_sp[i]=a_sp_fi[i][is1][is2];
  }  
	a_sp[n_sp-1] = 0.0;

	SPL(n_sp, R_sp, a_sp, 1, fip_rmin[is1][is2], 0.0, b_sp, c_sp, d_sp);

  // std::cout << "!!!!! PairLS debug mode !!!!! " << " SPL function worked well" << std::endl;

  std::cout << "!!!!!!!!!!!! Fi" << std::endl;
	for (i = 0; i < n_sp; i++)
  {
    a_sp_fi[i][is1][is2] = a_sp[i];
    b_sp_fi[i][is1][is2] = b_sp[i];
    c_sp_fi[i][is1][is2] = c_sp[i];
    d_sp_fi[i][is1][is2] = d_sp[i];
    std::cout << " R_sp_fi[" << i << "][" << is1 << "][" << is2 << "] = " << R_sp_fi[i][is1][is2] << "   ";
    std::cout << " a_sp_fi[" << i << "][" << is1 << "][" << is2 << "] = " << a_sp_fi[i][is1][is2] << "   ";    
    std::cout << " b_sp_fi[" << i << "][" << is1 << "][" << is2 << "] = " << b_sp_fi[i][is1][is2] << "   ";    
    std::cout << " c_sp_fi[" << i << "][" << is1 << "][" << is2 << "] = " << c_sp_fi[i][is1][is2] << "   ";    
    std::cout << " d_sp_fi[" << i << "][" << is1 << "][" << is2 << "] = " << d_sp_fi[i][is1][is2] << std::endl;    
  }
	
	shag_sp_fi[is1][is2] = 1.0/((R_sp_fi[n_sp-1][is1][is2]-R_sp_fi[0][is1][is2])/(n_sp-1));


	r1 = Rmin_fi_ZBL[is1][is2];
	f1 = v_ZBL(r1, is1, is2) + e0_ZBL[is1][is2];
	fp1 = vp_ZBL(r1, is1, is2);
	fpp1 = vpp_ZBL(r1, is1, is2);
  r2 = R_sp_fi[0][is1][is2];
  f2 = a_sp_fi[0][is1][is2];
  fp2 = b_sp_fi[0][is1][is2];
  fpp2 = 2.0*c_sp_fi[0][is1][is2];	
	
  // std::cout << "!!!!! PairLS debug mode !!!!! " << " now smooth_zero_22 function will be called" << std::endl;

  smooth_zero_22(B6, r1, r2, f1, fp1, fpp1, f2, fp2, fpp2);
  // smooth_zero_22(B6, &r1, &r2, &f1, &fp1, &fpp1, &f2, &fp2, &fpp2);
  // std::cout << "!!!!! PairLS debug mode !!!!! " << " smooth_zero_22 function worked well" << std::endl;

  std::cout << "!!!!!!!!!!!! c_fi_ZBL" << std::endl;
  std::cout << " r1 = " << r1 <<std::endl;
  std::cout << " r2 = " << r2 <<std::endl;
  std::cout << " f1 = " << f1 <<std::endl;
  std::cout << " fp1 = " << fp1 <<std::endl;
  std::cout << " fpp1 = " << fpp1 <<std::endl;
  std::cout << " f2 = " << f2 <<std::endl;
  std::cout << " fp2 = " << fp2 <<std::endl;
  std::cout << " fpp2 = " << fpp2 <<std::endl;

  for (i = 0; i < 6; i++)
  {
    c_fi_ZBL[i][is1][is2] = B6[i];
    std::cout << c_fi_ZBL[i][is1][is2] << "  ";
  }
  std::cout << std::endl;
  // std::cout << "!!!!! PairLS debug mode !!!!! " << " for (i = 0; i < 6; i++) worked well" << std::endl;


  // ro
  n_sp = n_sp_ro[is1][is2];
  std::cout << " n_sp = " << n_sp_ro[is1][is2] << std::endl;
  for (i = 0; i < n_sp; i++)
  {
    // std::cout << " i = " << i << " R_sp_ro[i][is][is] = " << R_sp_ro[i][is1][is2] << std::endl;
    R_sp[i] = R_sp_ro[i][is1][is2];
    a_sp[i] = a_sp_ro[i][is1][is2];
  }
  p1=0.0;

  // std::cout << "!!!!! PairLS debug mode !!!!! " << " for (i = 0; i < n_sp; i++) worked well" << std::endl;


	// SPL(n_sp, &R_sp, &a_sp, 1, p1, 0.0, &b_sp, &c_sp, &d_sp);
	SPL(n_sp, R_sp, a_sp, 1, p1, 0.0, b_sp, c_sp, d_sp);
  // std::cout << "!!!!! PairLS debug mode !!!!! " << " SPL(n_sp, R_sp, a_sp, 1, p1, 0.0, b_sp, c_sp, d_sp) worked well" << std::endl;

  std::cout << "!!!!!!!!!!!! Rho" << std::endl;
	for (i = 0; i < n_sp; i++)
  {  
    a_sp_ro[i][is1][is2] = a_sp[i];
    b_sp_ro[i][is1][is2] = b_sp[i];
    c_sp_ro[i][is1][is2] = c_sp[i];
    d_sp_ro[i][is1][is2] = d_sp[i];
    std::cout << " R_sp_ro[" << i << "][" << is1 << "][" << is2 << "] = " << R_sp_ro[i][is1][is2] << "   ";
    std::cout << " a_sp_ro[" << i << "][" << is1 << "][" << is2 << "] = " << a_sp_ro[i][is1][is2] << "   ";    
    std::cout << " b_sp_ro[" << i << "][" << is1 << "][" << is2 << "] = " << b_sp_ro[i][is1][is2] << "   ";    
    std::cout << " c_sp_ro[" << i << "][" << is1 << "][" << is2 << "] = " << c_sp_ro[i][is1][is2] << "   ";    
    std::cout << " d_sp_ro[" << i << "][" << is1 << "][" << is2 << "] = " << d_sp_ro[i][is1][is2] << std::endl;     
  }
	shag_sp_ro[is1][is2] = 1.0/((R_sp_ro[n_sp-1][is1][is2]-R_sp_ro[0][is1][is2])/(n_sp-1));


  // f3

	n_sp = n_sp_f[is1][is2];
  for (i = 0; i < n_sp; i++)
  {
    R_sp[i] = R_sp_f[i][is1][is2];
  }

  std::cout << "!!!!!!!!!!!! f3" << std::endl;
  for (i1 = 0; i1 < n_f3[is2]; i1++)
  {
    for (i = 0; i < n_sp; i++)
    {
      a_sp[i] = a_sp_f3[i][i1][is1][is2];
    }
    p1 = 0.0;
    // SPL(n_sp, &R_sp, &a_sp, 1, p1, 0.0, &b_sp, &c_sp, &d_sp);
    SPL(n_sp, R_sp, a_sp, 1, p1, 0.0, b_sp, c_sp, d_sp);
    // std::cout << "!!!!! PairLS debug mode !!!!! " << " SPL(n_sp, R_sp, a_sp, 1, p1, 0.0, b_sp, c_sp, d_sp) worked well" << std::endl;
    for (i = 0; i < n_sp; i++)
    { 
      a_sp_f3[i][i1][is1][is2] = a_sp[i];
      b_sp_f3[i][i1][is1][is2] = b_sp[i];
      c_sp_f3[i][i1][is1][is2] = c_sp[i];
      d_sp_f3[i][i1][is1][is2] = d_sp[i];
      std::cout << " R_sp_f3[" << i << "][" << is1 << "][" << is2 << "] = " << R_sp_f[i][is1][is2] << "   ";
      std::cout << " a_sp_f3[" << i << "][" << i1 << "][" << is1 << "][" << is2 << "] = " << a_sp_f3[i][i1][is1][is2] << "   ";    
      std::cout << " b_sp_f3[" << i << "][" << i1 << "][" << is1 << "][" << is2 << "] = " << b_sp_f3[i][i1][is1][is2] << "   ";    
      std::cout << " c_sp_f3[" << i << "][" << i1 << "][" << is1 << "][" << is2 << "] = " << c_sp_f3[i][i1][is1][is2] << "   ";    
      std::cout << " d_sp_f3[" << i << "][" << i1 << "][" << is1 << "][" << is2 << "] = " << d_sp_f3[i][i1][is1][is2] << std::endl;      
    }
  }
	shag_sp_f[is1][is2] = 1.0/((R_sp_f[n_sp-1][is1][is2]-R_sp_f[0][is1][is2])/(n_sp-1));

}


// Subroutines for spline creation written by A.G. Lipnitskii and translated from Fortran to C++ 

// smooth_zero_22.f
void PairLS::smooth_zero_22(double *B, double R1, double R2, double f1, double fp1, double fpp1, double f2, double fp2, double fpp2)
{
  //c == calc sqear delta  ==>
  int N = 6, NRHS = 1, LDA = 6, LDB = 7, INFO = 1;
  // int IPIV[N];
  // double A[LDA][N];
  int *IPIV;
  double *A;
  memory->create(IPIV, N, "PairLS:smooth_zero_22_IPIV_w");
  // memory->create(A, LDA, N, "PairLS:smooth_zero_22_A_w");
  memory->create(A, LDA*N, "PairLS:smooth_zero_22_A_w");

  // double B[6][1];
  // double R1, R2, f1, fp1, fpp1, f2, fp2, fpp2;

  A[0] = 1;            // A[0][0] = 1;
  A[1] = 0;            // A[1][0] = 0;
  A[2] = 0;            // A[2][0] = 0;
  A[3] = 1;            // A[3][0] = 1;
  A[4] = 0;            // A[4][0] = 0;
  A[5] = 0;            // A[5][0] = 0;

  A[6] = R1;           // A[0][1] = R1;
  A[7] = 1;            // A[1][1] = 1;
  A[8] = 0;            // A[2][1] = 0;
  A[9] = R2;           // A[3][1] = R2;
  A[10] = 1;            // A[4][1] = 1;
  A[11] = 0;            // A[5][1] = 0;

  A[12] = pow(R1,2);    // A[0][2] = pow(R1,2);
  A[13] = 2*R1;         // A[1][2] = 2*R1;
  A[14] = 2;            // A[2][2] = 2;
  A[15] = pow(R2,2);    // A[3][2] = pow(R2,2);
  A[16] = 2*R2;         // A[4][2] = 2*R2;
  A[17] = 2;            // A[5][2] = 2;

  A[18] = pow(R1,3);    // A[0][3] = pow(R1,3);
  A[19] = 3*pow(R1,2);  // A[1][3] = 3*pow(R1,2);
  A[20] = 6*R1;         // A[2][3] = 6*R1;
  A[21] = pow(R2,3);    // A[3][3] = pow(R2,3);
  A[22] = 3*pow(R2,2);  // A[4][3] = 3*pow(R2,2);
  A[23] = 6*R2;         // A[5][3] = 6*R2;

  A[24] = pow(R1,4);    // A[0][4] = pow(R1,4);
  A[25] = 4*pow(R1,3);  // A[1][4] = 4*pow(R1,3);
  A[26] = 12*pow(R1,2); // A[2][4] = 12*pow(R1,2);
  A[27] = pow(R2,4);    // A[3][4] = pow(R2,4);
  A[28] = 4*pow(R2,3);  // A[4][4] = 4*pow(R2,3);
  A[29] = 12*pow(R2,2); // A[5][4] = 12*pow(R2,2);

  A[30] = pow(R1,5);    // A[0][5] = pow(R1,5);
  A[31] = 5*pow(R1,4);  // A[1][5] = 5*pow(R1,4);
  A[32] = 20*pow(R1,3); // A[2][5] = 20*pow(R1,3);
  A[33] = pow(R2,5);    // A[3][5] = pow(R2,5);
  A[34] = 5*pow(R2,4);  // A[4][5] = 5*pow(R2,4);
  A[35] = 20*pow(R2,3); // A[5][5] = 20*pow(R2,3);

  B[0] = f1;
  B[1] = fp1;
  B[2] = fpp1;
  B[3] = f2;
  B[4] = fp2;
  B[5] = fpp2;

  // std::cout << "!!!!! PairLS debug mode !!!!! " << " entering dgesv" << std::endl;

  dgesv_(&N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO);
  // dgesv_(LAPACK_COL_MAJOR, N, NRHS, A, LDA, IPIV, B, LDB);
  // std::cout << "!!!!! PairLS debug mode !!!!! " << " _dgesv worked well" << std::endl;

  memory->destroy(IPIV);
  memory->destroy(A);
};


// SPL.f90
void PairLS::SPL(int n, double *X, double *Y, int ib, double D1, double DN, double *B, double *C, double *D)
{
  // int ib, n;         // intent(in)
  // double X[n], Y[n]; // intent(in)
  // double D1, DN;     // intent(in)
  // double B[n], C[n], D[n];
  double *A, *S;
  double t1, t2, t3;
  int i, n1, nn, Err;
  //begin
  // n = size(X)
  // std::cout << "!!!!! PairLS debug mode !!!!! " << " entering SPL" << std::endl;
  if (n == 1) 
  {
      B[0] = 0.0; C[0] = 0.0; D[0] = 0.0;
      return;
  } else
  if (n == 2) 
  {
      B[0] = (Y[1] - Y[0])/(X[1] - X[0]);
      C[0] = 0.0; D[0] = 0.0;
      B[1] = 0.0; C[1] = 0.0; D[1] = 0.0;
      return;
  }
  n1 = n - 1;
  B[0] = X[1] - X[0]; B[n-1] = 0.0;
  C[0] = 0.0; C[1] = B[0];
  D[0] = (Y[1] - Y[0])/B[0]; D[1] = D[0];
  memory->create(A, n, "PairLS:SPL_A_w");
  memory->create(S, n, "PairLS:SPL_S_w");  
  // std::cout << "!!!!! PairLS debug mode !!!!! " << " A and S arrays created" << std::endl;

  // SPL.f90(27-32):
  // do i=2, n1
  //   B(i)=X(i+1)-X(i); C(i+1)=B(i)
  //   A(i)=2.0*(X(i+1)-X(i-1))
  //   D(i+1)=(Y(i+1)-Y(i))/B(i)
  //   D(i)=D(i+1)-D(i)
  // end do

  for (i = 1; i < n1; i++)
  {
    // std::cout << i << std::endl;
    B[i] = X[i + 1] - X[i]; 
    C[i + 1] = B[i];
    A[i] = 2.0*(X[i + 1] - X[i - 1]);
    D[i + 1] = (Y[i + 1] - Y[i])/B[i];
    D[i] = D[i + 1] - D[i];    
  }

  // std::cout << "!!!!! PairLS debug mode !!!!! " << " Cycle: for (i = 1; i < n1; i++) done" << std::endl;

  switch (ib)
  {
    case 1:
      A[0] = 2.0*B[0]; A[n-1] = 2.0*B[n1-1];
      D[0] = D[0] - D1; D[n-1] = DN - D[n-1];
      nn = n;
      break;
    case 2:
      A[0] = 6.0; A[n-1] = 6.0; B[0] = 0.0; C[n-1] = 0.0;
      D[0] = D1; D[n-1] = DN;
      nn = n;
      break;
    case 3:
      D[0] = D[0] - D[n-1];   
      if (n == 3) 
      {
        A[0] = X[2] - X[0]; 
        A[1] = A[0]; 
        A[2] = A[0]; 
        D[2] = D[0];
        B[0] = 0.0; 
        B[1] = 0.0; 
        C[1] = 0.0; 
        C[2] = 0.0;
        nn = n;   // maybe it should be nn = n - 1
      } 
      else
      {
        A[0] = 2.0*(B[0] + B[n1-1]); C[0] = B[n1-1];
        nn = n1;  // maybe it should be nn = n1 - 1
      }     
      break;
    default:
      A[0] = -B[0]; A[n] = -B[n1-1];
      if (n == 3) 
      {
        D[0] = 0.0; D[2] = 0.0;
      } 
      else
      {
        D[0] = D[2]/(X[3] - X[1]) - D[1]/(X[2] - X[0]);
        D[n-1] = D[n1-1]/(X[n-1] - X[n-3]) - D[n-3]/(X[n1-1] - X[n-4]);
        D[0] = -D[0]*B[0]*B[0]/(X[3] - X[0]);
        D[n-1] = D[n-1]*B[n1-1]*B[n1-1]/(X[n-1] - X[n-4]);
      }
      nn = n;
      break;
  }    

  // std::cout << "!!!!! PairLS debug mode !!!!! " << " now LA30 function will be called" << std::endl;
  LA30(nn, A, B, C, D, S, &Err);
  // LA30(nn, A[0:nn-1], B[0:nn-1], C[0:nn-1], D[0:nn-1], S[0:nn-1], Err);
  // std::cout << "!!!!! PairLS debug mode !!!!! " << " LA30 worked well" << std::endl;

  B[0] = X[1] - X[0];
  if (ib == 3) 
  {
    S[n-1] = S[0]; B[1] = X[2] - X[1];
  }
  for (i = 0; i < n1; i++)
  {
    D[i] = (S[i + 1] - S[i])/B[i];
    C[i] = 3.0*S[i];
    B[i] = (Y[i + 1] - Y[i])/B[i] - B[i]*(S[i + 1] + 2.0*S[i]);
  }
  D[n-1] = D[n1-1]; C[n-1] = 3.0*S[n-1]; B[n-1] = B[n1-1];

  memory->destroy(A);
  memory->destroy(S);
}

// LA30.f
void PairLS::LA30(int n, double *A, double *B, double *C, double *D, double *X, int *Error)
{
  //
  // int(4), intent(in):: n;
  // float(8), intent(in):: A[n], B[n], C[n], D[n];
  // float(8):: X[n];
  // int(4):: Error;
  // float(8), allocatable:: P[:], Q[:], R[:], S[:], T[:];
  double *P, *Q, *R, *S, *T;
  double W;
  int i, ii;
  //begin
  // n = size(A)
  if (n < 3) 
  {
    *Error = 1; 
    return;
  }

  memory->create(P, n, "PairLS:LA30_P_w");
  memory->create(Q, n, "PairLS:LA30_Q_w");
  memory->create(R, n, "PairLS:LA30_R_w");
  memory->create(S, n, "PairLS:LA30_S_w");
  memory->create(T, n, "PairLS:LA30_T_w");
  P[0] = 0.0; 
  Q[0] = 0.0; 
  R[0] = 1.0;

  for (i = 0; i < n-1; i++)
  {
    ii = i + 1;
    W = A[i] + Q[i]*C[i];
    if (1.0 + W == 1.0) 
    {
      memory->destroy(P);
      memory->destroy(Q);
      memory->destroy(R);
      memory->destroy(S);
      memory->destroy(T);
      *Error = 65; 
      return;
    }
    P[ii] = (D[i] - P[i]*C[i])/W;
    Q[ii] = -B[i]/W;
    R[ii] = -R[i]*C[i]/W;
  }

  S[n-1] = 1.0; 
  T[n-1] = 0.0;
  for (i = n-2; i >= 0; i--) // check for consistency with LA30.f in testing
  {
    ii = i + 1;
    S[i] = Q[ii]*S[ii] + R[ii];
    T[i] = Q[ii]*T[ii] + P[ii];
  }


  W = A[n-1] + B[n-1]*S[0] + C[n-1]*S[n-2];

  if (1.0 + W == 1.0) 
  {
    memory->destroy(P);
    memory->destroy(Q);
    memory->destroy(R);
    memory->destroy(S);
    memory->destroy(T);
    *Error = 65; 
    return;
  }

  X[n-1] = (D[n-1] - B[n-1]*T[0] - C[n-1]*T[n-2])/W;
  for (i = 0; i < n-1; i++)
  {
    X[i] = S[i]*X[n-1] + T[i];
  }

  memory->destroy(P);
  memory->destroy(Q);
  memory->destroy(R);
  memory->destroy(S);
  memory->destroy(T);
  *Error = 0; 
  return;
}

// functions for calculating energies and forces

// e_force_fi_emb.f
// void PairLS::e_force_fi_emb(double *e_at, double **f_at, double *px_at, double *py_at, double *pz_at,  double **r_at, int *i_sort_at, int n_at, double *sizex, double *sizey, double *sizez)
void PairLS::e_force_fi_emb(double *e_at, double **f_at, double *px_at, double *py_at, double *pz_at,  double **r_at, int *i_sort_at, int n_at, double sizex, double sizey, double sizez)
{

  // Array Arguments
  // int *i_sort_at = atom->type // aary with atom types that must range from 1 to specified # of types.
  // Scalar Arguments
  double e_sum, pressure, sxx, syy, szz;
  
  // Local Scalars
  int i, j, ii, jj, jnum, is, js;
  int n_i, in_list;
  double rr_pot[10][10], size[3];
  double sizex05, sizey05, sizez05, x, y, z, xx, yy, zz, rr, r, r1;
  double roi, roj, w, ropi, ropj, w1, w2, w3;
  int i1, i2, i3, iw;

  // Local Arrays
  double *rosum;
  memory->create(rosum, max_at, "PairLS:rosum");

  // pointers to LAMMPS arrays 
  int *ilist,*jlist,*numneigh,**firstneigh;
  // double **x = atom->x;
  // double **f = atom->f;
  // int *type = atom->type; // Atom types must range from 1 to specified # of types.
  // int nlocal = atom->nlocal;
  // int newton_pair = force->newton_pair;
  int inum = list->inum; // # of I atoms neighbors are stored for

  ilist = list->ilist;           // local indices of I atoms
  numneigh = list->numneigh;     // # of J neighbors for each I atom
  firstneigh = list->firstneigh; // ptr to 1st J int value of each I atom

  // time1 = timesec(idum);

  sizex05 = sizex*0.5D0;
  sizey05 = sizey*0.5D0;
  sizez05 = sizez*0.5D0;

  for (is = 1; is <= n_sort; is++)
  {
    for (js = 1; js <= n_sort; js++)
    {
      rr_pot[js][is] = pow(R_sp_fi[n_sp_fi-1][js][is],2);
    }
  }
  
  e_at = {0.0};
  rosum = {0.0};

  // == calc rosum(n_at) =========================>

//////////////////////////////////////////////////
  // This is how an original loop over the md_ls Verlet's list is look 
  // i = 0; // start atom i;
  // in_list = 0; // start list of atom i;
  // while (i < n_at)
  // {
  //   // if (i > n_at) break;
  //   n_i = neighbour_list(in_list) // number of neighbours of atom i;
  //   if (n_i == 0) continue;
  //   x = r_at(1, i); 
  //   y = r_at(2, i); 
  //   z = r_at(3, i); 
  //   is = i_sort_at(i);

  //   // Loop over neighbours of atom i
  //   for (jj = in_list + 1; jj = in_list + n_i; jj++)
  //   {
  //     j = neighbour_list(jj);
  //     js = i_sort_at(j);
  //     xx = r_at(1, j) - x; 
  //     yy = r_at(2, j) - y; 
  //     zz = r_at(3, j) - z;;
  //     if (periodic(1)) && (xx > sizex05) xx = xx - sizex;
  //     if (periodic(1)) && (xx < - sizex05) xx = xx + sizex;
  //     if (periodic(2)) && (yy > sizey05) yy = yy - sizey;
  //     if (periodic(2)) && (yy < - sizey05) yy = yy + sizey;
  //     if (periodic(3)) && (zz > sizez05) zz = zz - sizez;
  //     if (periodic(3)) && (zz < - sizez05) zz = zz + sizez;
  //     rr = xx**2 + yy**2 + zz**2;
  //     if (rr < rr_pot[js, is]) 
  //     {
  //     r = dsqrt(rr);
  //     roj = fun_ro(r, js, is);
  //     rosum[i] = rosum[i] + roj;
  //     if (if_true_i(i)) 
  //     {
  //     w = fun_fi(r, is, js);
  //     e_at(i) = e_at(i) + w;
  //   }// end loop over jj;
  //   // g111:
  //   // Set starting position of next atom i
  //   i = i + 1;
  //   in_list = in_list + n_i + 1;
  // };

//////////////////////////////////////////////////

  // Loop over the LAMMPS neighbour list 
  std::cout << " i indexes of atoms in the LAMMPS neighbour list" << "  "; 
  for (ii = 0; ii < inum; ii++) // inum is the # of I atoms neighbors are stored for
  {
    i = ilist[ii];  // index of atom i
    std::cout << i << "  ";
    x = r_at[i][0];
    y = r_at[i][1];
    z = r_at[i][2];
    is = i_sort_at[i];
    jlist = firstneigh[i]; // ptr to 1st J int value of the atom i
    jnum = numneigh[i];    // # of J neighbors for the atom i

    for (jj = 0; jj < jnum; jj++) 
    {
      j = jlist[jj];  // index of atom j
      js = type[j];
      xx = r_at[j][1] - x;
      yy = r_at[j][2] - y;
      zz = r_at[j][3] - z;
      // if (domain->xperiodic && xx >  domain->xprd_half) xx = xx - sizex;
      // if (domain->xperiodic && xx < -domain->xprd_half) xx = xx + sizex;
      // if (domain->yperiodic && yy >  domain->yprd_half) yy = yy - sizey;
      // if (domain->yperiodic && yy < -domain->yprd_half) yy = yy + sizey;
      // if (domain->zperiodic && zz >  domain->zprd_half) zz = zz - sizez;
      // if (domain->zperiodic && zz < -domain->zprd_half) zz = zz + sizez;
      /////////////////////////////////////
      // if (periodic[0] && xx >  sizex05) xx = xx - sizex;
      // if (periodic[0] && xx < -sizex05) xx = xx + sizex;
      // if (periodic[1] && yy >  sizey05) yy = yy - sizey;
      // if (periodic[1] && yy < -sizey05) yy = yy + sizey;
      // if (periodic[2] && zz >  sizez05) zz = zz - sizez;
      // if (periodic[2] && zz < -sizez05) zz = zz + sizez;
      // Start debugging version
      // if (domain->xperiodic && xx >  domain->xprd_half) 
      if (periodic[0] && xx > sizex05) 
      {
        std::cout << " domain->xperiodic && xx > domain->xprd_half is true" << std::endl;
        xx = xx - sizex;
      }
      // if (domain->xperiodic && xx < -domain->xprd_half) 
      if (periodic[0] && xx < -sizex05) 
      {
        std::cout << " domain->xperiodic && xx < -domain->xprd_half is true" << std::endl;
        xx = xx + sizex;
      }
      // if (domain->yperiodic && yy >  domain->yprd_half) 
      if (periodic[1] && yy > sizey05) 
      {
        std::cout << " domain->yperiodic && yy >  domain->yprd_half is true" << std::endl;
        yy = yy - sizey;
      }
      // if (domain->yperiodic && yy < -domain->yprd_half) 
      if (periodic[1] && yy < -sizey05) 
      {
        std::cout << " domain->yperiodic && yy < -domain->yprd_half is true" << std::endl;
        yy = yy + sizey;
      }
      // if (domain->zperiodic && zz >  domain->zprd_half) 
      if (periodic[2] && zz > sizez05) 
      {
        std::cout << " domain->zperiodic && zz >  domain->zprd_half is true" << std::endl;
        zz = zz - sizez;
      }
      // if (domain->zperiodic && zz < -domain->zprd_half) 
      if (periodic[2] && zz < -sizez05) 
      {
        std::cout << " domain->zperiodic && zz < -domain->zprd_half is true" << std::endl;
        zz = zz + sizez;  
      }
      // End debugging version
      rr = xx**2 + yy**2 + zz**2;
      if (rr < rr_pot[js][is]) 
      {
        r = sqrt(rr);
        roj = fun_ro(r, js, is);
        rosum[i] = rosum[i] + roj;
        if (if_true_i[i])     // check what it mean
        {
          w = fun_fi(r, is, js);
          e_at[i] = e_at[i] + w; 
        }
      }
    }
  }
//////////////////////////////////////////////////

  // == calc energies: e_at(n_at) ==>
  for (i = 0; i < n_at; i++)
  {
    if (.not.if_true_i(i)) cycle;
    w = fun_emb(rosum[i], i_sort_at(i)) + 0.5D0*e_at(i);
    e_at(i) = w;
  }


  //c
  //c
  //c
  //c == calc funp_emb(n_at) and put one into rosum(n_at) ==>
  do i = 1, n_at
  //c	    if(.not.if_true_i(i)) cycle
  rosum[i] = funp_emb(rosum[i], i_sort_at(i));
  }do;


  //c
  //c == calc forces: f_at(3,n_at) ==============================>
  //c
  f_at = 0.0D0;
  px_at = 0.0D0;
  py_at = 0.0D0;
  pz_at = 0.0D0;
  pressure = 0.0D0;

  // Loop over the Verlet"s listundefinedundefined
  i = 1 // start atom i;
  in_list = 1 // start list of atom i;
  do;
  if (i > n_at) exit;
  n_i = neighbour_list(in_list) // number of neighbours of atom i;
  if (n_i =  = 0) goto g222;
  if (.not.if_true_i(i)) goto g222;
  x = r_at(1, i); y = r_at(2, i); z = r_at(3, i); is = i_sort_at(i);;
  // Loop over neighbours of atom i
  do jj = in_list + 1, in_list + n_i
  j = neighbour_list(jj);
  js = i_sort_at(j);
  xx = r_at(1, j) - x; yy = r_at(2, j) - y; zz = r_at(3, j) - z;;
  if (periodic(1)) && (xx > sizex05) xx = xx - sizex;
  if (periodic(1)) && (xx < - sizex05) xx = xx + sizex;
  if (periodic(2)) && (yy > sizey05) yy = yy - sizey;
  if (periodic(2)) && (yy < - sizey05) yy = yy + sizey;
  if (periodic(3)) && (zz > sizez05) zz = zz - sizez;
  if (periodic(3)) && (zz < - sizez05) zz = zz + sizez;
  rr = xx**2 + yy**2 + zz**2;
  if (rr < rr_pot[js, is]) 
  {
  r = dsqrt(rr);
  r1 = 1.0D0/r;
  ropi = funp_ro(r, is, js);
  if (js =  = is) 
  {
  ropj = ropi;
  } else
  {
  ropj = funp_ro(r, js, is);
  }
  w = ((rosum[i]*ropj + rosum[j]*ropi) + funp_fi(r, is, js))*r1;
  w1 = w*xx; w2 = w*yy; w3 = w*zz;;
  f_at(1, i) = f_at(1, i) + w1;
  f_at(2, i) = f_at(2, i) + w2 // add the fi force j = >i;
  f_at(3, i) = f_at(3, i) + w3;
  w1 = w1*xx;
  w2 = w2*yy;
  w3 = w3*zz;
  px_at(i) = px_at(i) + w1;
  py_at(i) = py_at(i) + w2;
  pz_at(i) = pz_at(i) + w3;
  //c			px_at(j)=px_at(j)+w1
  //c		        py_at(j)=py_at(j)+w2
  //c			pz_at(j)=pz_at(j)+w3
  }
  }do // end loop over jj;



  g222:

  // Set starting position of next atom i
  i = i + 1;
  in_list = in_list + n_i + 1;
  }do;


  //c
  //c	sxx=0.0D0
  //c	syy=0.0D0
  //c	szz=0.0D0
  //c	    e_sum=0.0D0

  w = 0.5D0*(-1.0D0)/(sizex*sizey*sizez);
  do i = 1, n_at
  //c	if(.not.if_true_i(i)) cycle
  //c	    e_sum=e_sum+e_at(i)
  //c	sxx=sxx + px_at(i)
  //c	syy=syy + py_at(i)
  //c	szz=szz + pz_at(i)
  px_at(i) = w*px_at(i);
  py_at(i) = w*py_at(i);
  pz_at(i) = w*pz_at(i);
  }do;


  ;
  //c	sxx=0.5D0*(-1.0D0)*sxx/(sizex*sizey*sizez)
  //c	syy=0.5D0*(-1.0D0)*syy/(sizex*sizey*sizez)
  //c	szz=0.5D0*(-1.0D0)*szz/(sizex*sizey*sizez)
  ;
  //c	pressure=(1.0D0/3.0D0)*(sxx+syy+szz)
  //c
  // g777:

  // time2 = timesec(idum);
  // t_e_force_fi_emb = t_e_force_fi_emb + (time2 - time1);
  //c
  memory->destroy(rosum);
  return;
}

void PairLS::e_force_g3(double *, double *, double *, double *, double *,  double **, int *, int, double *, double *, double *)
{

}







// Potential functions from fun_pot_ls.f

// fun_pot_ls.f(3-36):
double PairLS::fun_fi(double r , int is, int js)
{
  int i;
  double fun_fi, r, dr, r0_min;

  if (r >= R_sp_fi[n_sp_fi-1][is][js]) 
  {
    fun_fi = 0.0;
    return fun_fi;
  }

  if (r < Rmin_fi_ZBL[is][js]) 
  {
    fun_fi = v_ZBL(r, is, js) + e0_ZBL[is][js];
    return fun_fi;
  }

  r0_min = R_sp_fi[0][is][js];

  if (r < r0_min) 
  {
    fun_fi = fun_fi_ZBL(r, is, js);
    return fun_fi;
  }

  i = int((r - R_sp_fi[0][is][js])*shag_sp_fi[is][js]);
  i = i + 1;
  if (i < 1) i = 1;
  dr = r - R_sp_fi[i-1][is][js];
  fun_fi = a_sp_fi[i-1][is][js] + dr*(b_sp_fi[i-1][is][js] + dr*(c_sp_fi[i-1][is][js] + dr*(d_sp_fi[i-1][is][js])));

  return fun_fi;

}

// fun_pot_ls.f(39-70):
double PairLS::funp_fi(double r, int is, int js)
{
  int i;
  double funp_fi, r, dr, r0_min;

  if (r >= R_sp_fi[n_sp_fi-1][is][js]) 
  {
    funp_fi = 0.0;
    return funp_fi;
  }

  if (r < Rmin_fi_ZBL[is][js]) 
  {
    funp_fi = vp_ZBL(r, is, js);
    return funp_fi;
  }

  r0_min = R_sp_fi[0][is][js];

  if (r < r0_min) 
  {
    funp_fi = funp_fi_ZBL(r, is, js);
    return funp_fi;
  }

  i = int((r - R_sp_fi[0][is][js])*shag_sp_fi[is][js]);
  i = i + 1;
  if (i < 1) i = 1;
  dr = r - R_sp_fi[i-1][is][js];
  funp_fi = b_sp_fi[i-1][is][js] + dr*(2.0*c_sp_fi[i-1][is][js] + dr*(3.0*d_sp_fi[i-1][is][js]));

  return funp_fi;
}

// fun_pot_ls.f(74-106):
double PairLS::funpp_fi(double r, int is, int js)
{
  int i;
  double funpp_fi, r, dr, r0_min;

  if (r >= R_sp_fi[n_sp_fi-1][is][js]) 
  {
    funpp_fi = 0.0;
    return funpp_fi;
  }

  if (r < Rmin_fi_ZBL[is][js]) 
  {
    funpp_fi = vpp_ZBL(r, is, js);
    return funpp_fi;
  }

  r0_min = R_sp_fi[0][is][js];

  if (r < r0_min) 
  {
    funpp_fi = funpp_fi_ZBL(r, is, js);
    return funpp_fi;
  }

  i = int((r - r0_min)*shag_sp_fi[is][js]);
  i = i + 1;
  if (i < 1) i = 1;
  dr = r - R_sp_fi[i-1][is][js];
  funpp_fi = 2.0*c_sp_fi[i-1][is][js] + dr*(6.0*d_sp_fi[i-1][is][js]);

  return funpp_fi;
}



// fun_pot_ls.f(112-139):
double PairLS::fun_ro(double r, int is, int js)
{
  int i;
  double fun_ro, r, dr, r0_min;

  if (r >= R_sp_ro[n_sp_ro-1][is][js]) 
  {
    fun_ro = 0.0;
    return fun_ro;
  }

  r0_min = R_sp_ro[0][is][js];

  if (r < r0_min) 
  {
    fun_ro = a_sp_ro[0][is][js];
    return fun_ro;
  }

  i = int((r - r0_min)*shag_sp_ro[is][js]);
  i = i + 1;
  dr = r - R_sp_ro[i-1][is][js];
  fun_ro = a_sp_ro[i-1][is][js] + dr*(b_sp_ro[i-1][is][js] + dr*(c_sp_ro[i-1][is][js] + dr*(d_sp_ro[i-1][is][js])));

  return fun_ro;
}

// fun_pot_ls.f(142-169):
double PairLS::funp_ro(double r, int is, int js)
{
  int i;
  double funp_ro, r, dr, r0_min;

  if (r >= R_sp_ro[n_sp_ro-1][is][js]) 
  {
    funp_ro = 0.0;
    return funp_ro;
  }

  r0_min = R_sp_ro[0][is][js];

  if (r < r0_min) 
  {
    funp_ro = 0.0;
    return funp_ro;
  }

  i = int((r - r0_min)*shag_sp_ro[is][js]);
  i = i + 1;
  dr = r - R_sp_ro[i-1][is][js];
  funp_ro = b_sp_ro[i-1][is][js] + dr*(2.0*c_sp_ro[i-1][is][js] + dr*(3.0*d_sp_ro[i-1][is][js]));

  return funp_ro;  
}

// fun_pot_ls.f(173-200):
double PairLS::funpp_ro(double r, int is, int js)
{
  int i;
  double funpp_ro, r, dr, r0_min;

  if (r >= R_sp_ro[n_sp_ro-1][is][js]) 
  {
    funpp_ro = 0.0;
    return funpp_ro;
  }

  r0_min = R_sp_ro[0][is][js];

  if (r < r0_min) 
  {
    funpp_ro = 0.0;
    return funpp_ro;
  }

  i = int((r - r0_min)*shag_sp_ro[is][js]);
  i = i + 1;
  if (i <= 0) i = 1;
  dr = r - R_sp_ro[i-1][is][js];
  funpp_ro = 2.0*c_sp_ro[i-1][is][js] + dr*(6.0*d_sp_ro[i-1][is][js]);

  return funpp_ro;  
}  


// fun_pot_ls.f(209-245):
double PairLS::fun_emb(double r, int is)
{
  int i;
  double fun_emb, r, dr, r0_min;

  if (r >= R_sp_emb(n_sp_emb, is)) 
  {
      fun_emb = a_sp_emb[n_sp_emb-1][is];
      return fun_emb;
  }

  r0_min = R_sp_emb[0][is];

  if (r <= r0_min) 
  {
      fun_emb = b_sp_emb[0][is]*(r - r0_min);
      return fun_emb;
  }

  //c    	
  i = int((r - R_sp_emb[0][is])*shag_sp_emb[is]);
  i = i + 1;
  dr = r - R_sp_emb[i-1][is];
  fun_emb = a_sp_emb[i-1][is] + dr*(b_sp_emb[i-1][is] + dr*(c_sp_emb[i-1][is] + dr*(d_sp_emb[i-1][is])));

  return fun_emb;  
}

// fun_pot_ls.f(248-273):
double PairLS::funp_emb(double r, int is)
{
  int i;
  double funp_emb, r, dr, r0_min;

  if (r >= R_sp_emb(n_sp_emb, is)) 
  {
      funp_emb = 0.0;
      return funp_emb;
  }

  r0_min = R_sp_emb[0][is];

  if (r <= r0_min) 
  {
      funp_emb = b_sp_emb[0][is];
      return funp_emb;
  }

  //c    	
  i = int((r - r0_min)*shag_sp_emb[is]);
  i = i + 1;
  dr = r - R_sp_emb[i-1][is];
  funp_emb = b_sp_emb[i-1][is] + dr*(2.0*c_sp_emb[i-1][is] + dr*(3.0*d_sp_emb[i-1][is]));

  return funp_emb;  
}

// fun_pot_ls.f(285-312):
double PairLS::funpp_emb(double r, int is)
{
  int i;
  double funpp_emb, r, dr, r0_min;

  if (r >= R_sp_emb(n_sp_emb, is)) 
  {
      funpp_emb = 0.0;
      return funpp_emb;
  }

  r0_min = R_sp_emb[0][is];

  if (r <= r0_min) 
  {
      funpp_emb = 0.0;
      return funpp_emb;
  }

  //c    	
  i = int((r - r0_min)*shag_sp_emb[is]);
  i = i + 1;
  if(i <= 0) i = 1; // maybe this condition should be added also for fun_emb and funp_emb?
  dr = r - R_sp_emb[i-1][is];
  funpp_emb = 2.0*c_sp_emb[i-1][is] + dr*(6.0*d_sp_emb[i-1][is]));

  return funpp_emb;    
} 


// fun_pot_ls.f(319-347):
double PairLS::fun_f3(double r, int i_f3, int js, int is)
{
  int i;
  double fun_f3, r, dr, r0_min;

  if (r >= R_sp_f[n_sp_f-1][js][is]) 
  {
      fun_f3 = 0.0;
      return fun_f3;
  }

  r0_min = R_sp_f[0][js][is];

  if (r <= r0_min) 
  {
      fun_f3 = a_sp_f3[0][i_f3-1][js][is];
      return fun_f3;
  }

  i = int((r - r0_min)*shag_sp_f[js][is]);
  i = i + 1;
  dr = r - R_sp_f[i-1][js][is];
  // fun_f3 = a_sp_f3[i-1][i_f3-1][js][is] + dr*(b_sp_f3[i-1][i_f3-1][js][is] + dr*(c_sp_f3[i-1][i_f3-1][js][is] + dr*(d_sp_f3[i-1][i_f3-1][js][is])));
  fun_f3 = a_sp_f3[i-1][i_f3][js][is] + dr*(b_sp_f3[i-1][i_f3][js][is] + dr*(c_sp_f3[i-1][i_f3][js][is] + dr*(d_sp_f3[i-1][i_f3][js][is])));

  return fun_f3;
}

// fun_pot_ls.f(350-377):
double PairLS::funp_f3(double r, int i_f3, int js, int is)
{
  int i;
  double funp_f3, r, dr, r0_min;

  if (r >= R_sp_f[n_sp_f-1][js][is]) 
  {
      funp_f3 = 0.0;
      return funp_f3;
  }

  r0_min = R_sp_f[0][js][is];

  if (r <= r0_min) 
  {
      funp_f3 = 0.0;
      return funp_f3;
  }

  i = int((r - r0_min)*shag_sp_f[js][is]);
  i = i + 1;
  dr = r - R_sp_f[i-1][js][is];
  // funp_f3 = b_sp_f3[i-1][i_f3-1][js][is] + dr*(2.0*c_sp_f3[i-1][i_f3-1][js][is] + dr*(3.0*d_sp_f3[i-1][i_f3-1][js][is]));
  funp_f3 = b_sp_f3[i-1][i_f3][js][is] + dr*(2.0*c_sp_f3[i-1][i_f3][js][is] + dr*(3.0*d_sp_f3[i-1][i_f3][js][is]));

  return funp_f3;  
}

// fun_pot_ls.f(381-406):
double PairLS::funpp_f3(double r, int i_f3, int js, int is)
{
  int i;
  double funpp_f3, r, dr, r0_min;

  if (r >= R_sp_f[n_sp_f-1][js][is]) 
  {
      funpp_f3 = 0.0;
      return funpp_f3;
  }

  r0_min = R_sp_f[0][js][is];

  if (r <= r0_min) 
  {
      funpp_f3 = 0.0;
      return funpp_f3;
  }

  i = int((r - r0_min)*shag_sp_f[js][is]);
  i = i + 1;
  if (i <= 0) i = 1;
  dr = r - R_sp_f[i-1][js][is];
  // funpp_f3 = 2.0*c_sp_f3[i-1][i_f3-1][js][is] + dr*(6.0*d_sp_f3[i-1][i_f3-1][js][is]);
  funpp_f3 = 2.0*c_sp_f3[i-1][i_f3][js][is] + dr*(6.0*d_sp_f3[i-1][i_f3][js][is]);

  return funpp_f3;    
}


// fun_pot_ls.f(412-425):
double PairLS::fun_g3(double r, int i1, int i2, int is)
{
  int i;
  double fun_g3, r, dr;

  i = int((r - R_sp_g[0])*shag_sp_g);
  i = i + 1;
  if (i >= n_sp_g) i = n_sp_g - 1;
  dr = r - R_sp_g[i-1];
  // fun_g3 = a_sp_g3[i-1][i1-1][i2-1][is] + dr*(b_sp_g3[i-1][i1-1][i2-1][is] + dr*(c_sp_g3[i-1][i1-1][i2-1][is] + dr*(d_sp_g3[i-1][i1-1][i2-1][is])));
  fun_g3 = a_sp_g3[i-1][i1][i2][is] + dr*(b_sp_g3[i-1][i1][i2][is] + dr*(c_sp_g3[i-1][i1][i2][is] + dr*(d_sp_g3[i-1][i1][i2][is])));
  return fun_g3;
}

// fun_pot_ls.f(428-442):
double PairLS::funp_g3(double r, int i1, int i2, int is)
{
  int i;
  double funp_g3, r, dr;

  i = int((r - R_sp_g[0])*shag_sp_g);
  i = i + 1;
  if (i >= n_sp_g) i = n_sp_g - 1;
  dr = r - R_sp_g[i-1];
  // funp_g3 = b_sp_g3[i-1][i1-1][i2-1][is] + dr*(2.0*c_sp_g3[i-1][i1-1][i2-1][is] + dr*(3.0*d_sp_g3[i-1][i1-1][i2-1][is]));
  funp_g3 = b_sp_g3[i-1][i1][i2][is] + dr*(2.0*c_sp_g3[i-1][i1][i2][is] + dr*(3.0*d_sp_g3[i-1][i1][i2][is]));
  return funp_g3;
}

// fun_pot_ls.f(446-459):
double PairLS::funpp_g3(double r, int i1, int i2, int is)
{
  int i;
  double funpp_g3, r, dr;

  i = int((r - R_sp_g[0])*shag_sp_g);
  i = i + 1;
  if (i >= n_sp_g) i = n_sp_g - 1;
  dr = r - R_sp_g[i-1];
  // funpp_g3 = 2.0*c_sp_g3[i-1][i1-1][i2-1][is] + dr*(6.0*d_sp_g3[i-1][i1-1][i2-1][is]);
  funpp_g3 = 2.0*c_sp_g3[i-1][i1][i2][is] + dr*(6.0*d_sp_g3[i-1][i1][i2][is]);
  return funpp_g3;
}


// fun_pot_ls.f(603-623):
double PairLS::v_ZBL(double r, int is, int js)
{
  int i;
  double v_ZBL;
  double w, sum, zz_r;

  zz_r = zz_ZBL[is][js]/r;

  w = r/a_ZBL[is][js];

  sum = 0.0;
  for (i = 0; i < 4; i++)
  {
    // sum = sum + c_ZBL[i][is][js]*exp(-d_ZBL[i][is][js]*w);
    sum = sum + c_ZBL[i]*exp(-d_ZBL[i]*w);
  }

  v_ZBL = zz_r*sum;
        
  return v_ZBL;
}

// fun_pot_ls.f(627-655):
double PairLS::vp_ZBL(double r, int is, int js)
{
  int i;
  double vp_ZBL;
  double w, sum, sump, zz_r, zzp_r;

  zz_r = zz_ZBL[is][js]/r;
  zzp_r = -zz_r/r;

  w = r/a_ZBL[is][js];

  sum = 0.0;
  for (i = 0; i < 4; i++)
  {
    // sum = sum + c_ZBL[i][is][js]*exp(-d_ZBL[i][is][js]*w);
    sum = sum + c_ZBL[i]*exp(-d_ZBL[i]*w);
  }

  sump = 0.0;
  for (i = 0; i < 4; i++)
  {
    // sump = sump + c_ZBL[i][is][js]*exp(-d_ZBL[i][is][js]*w)*(-d_ZBL[i][is][js]/a_ZBL[is][js]);
    sump = sump + c_ZBL[i]*exp(-d_ZBL[i]*w)*(-d_ZBL[i]/a_ZBL[is][js]);
  }

  vp_ZBL = zzp_r*sum + zz_r*sump;
        
  return vp_ZBL;
}

// fun_pot_ls.f(659-694):
double PairLS::vpp_ZBL(double r, int is, int js)
{
  int i;
  double vpp_ZBL;
  double w, sum, sump, sumpp, zz_r, zzp_r, zzpp_r;

  zz_r = zz_ZBL[is][js]/r;
  zzp_r = -zz_r/r;
  zzpp_r = -2.0*zzp_r/r;

  w = r/a_ZBL[is][js];

  sum = 0.0;
  for (i = 0; i < 4; i++)
  {
    // sum = sum + c_ZBL[i][is][js]*exp(-d_ZBL[i][is][js]*w);
    sum = sum + c_ZBL[i]*exp(-d_ZBL[i]*w);
  }

  sump = 0.0;
  for (i = 0; i < 4; i++)
  {
    // sump = sump + c_ZBL[i][is][js]*exp(-d_ZBL[i][is][js]*w)*(-d_ZBL[i][is][js]/a_ZBL[is][js]);
    sump = sump + c_ZBL[i]*exp(-d_ZBL[i]*w)*(-d_ZBL[i]/a_ZBL[is][js]);
  }

  sumpp = 0.0;
  for (i = 0; i < 4; i++)
  {
    // sumpp = sumpp + c_ZBL[i][is][js]*exp(-d_ZBL[i][is][js]*w)*pow((d_ZBL[i][is][js]/a_ZBL[is][js]),2);
    sumpp = sumpp + c_ZBL[i]*exp(-d_ZBL[i]*w)*pow((d_ZBL[i]/a_ZBL[is][js]),2);
  }

  vpp_ZBL = zzpp_r*sum + 2.0*zzp_r*sump + zz_r*sumpp;
        
  return vpp_ZBL;
}


// fun_pot_ls.f(698-711):
double PairLS::fun_fi_ZBL(double r, int is, int js)
{
  double fun_fi_ZBL;

  fun_fi_ZBL = c_fi_ZBL[0][is][js] + r*(c_fi_ZBL[1][is][js] + r*(c_fi_ZBL[2][is][js] + r*(c_fi_ZBL[3][is][js] + r*(c_fi_ZBL[4][is][js] + r*(c_fi_ZBL[5][is][js])))));

  return fun_fi_ZBL;
}

// fun_pot_ls.f(715-727):
double PairLS::funp_fi_ZBL(double r, int is, int js)
{
  double funp_fi_ZBL;

  funp_fi_ZBL = c_fi_ZBL[1][is][js] + r*(2.0*c_fi_ZBL[2][is][js] + r*(3.0*c_fi_ZBL[3][is][js] + r*(4.0*c_fi_ZBL[4][is][js] + r*(5.0*c_fi_ZBL[5][is][js]))));

  return funp_fi_ZBL;
}

// fun_pot_ls.f(731-742):
double PairLS::funpp_fi_ZBL(double r, int is, int js)
{
  double funpp_fi_ZBL;

  funpp_fi_ZBL = 2.0*c_fi_ZBL[2][is][js] + r*(6.0*c_fi_ZBL[3][is][js] + r*(12.0*c_fi_ZBL[4][is][js] + r*(20.0*c_fi_ZBL[5][is][js])));

  return funpp_fi_ZBL;
}
















































// // LAPACK subroutines and functions translated from Fortran to C++
// // dgesv.f
// void PairLS::DGESV(int N, int NRHS, double **A, int LDA, int *IPIV, double **B, int LDB, int *INFO)
// {

//   // *  DGESV computes the solution to a real system of linear equations
//   // *     A * X = B,
//   // *  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
//   // *
//   // *  The LU decomposition with partial pivoting and row interchanges is
//   // *  used to factor A as
//   // *     A = P * L * U,
//   // *  where P is a permutation matrix, L is unit lower triangular, and U is
//   // *  upper triangular.  The factored form of A is then used to solve the
//   // *  system of equations A * X = B.
//   // 
//   // *  Arguments
//   // *  =========
//   // *
//   // *  N       (input) INTEGER
//   // *          The number of linear equations, i.e., the order of the
//   // *          matrix A.  N >= 0.
//   // *
//   // *  NRHS    (input) INTEGER
//   // *          The number of right hand sides, i.e., the number of columns
//   // *          of the matrix B.  NRHS >= 0.
//   // *
//   // *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
//   // *          On entry, the N-by-N coefficient matrix A.
//   // *          On exit, the factors L and U from the factorization
//   // *          A = P*L*U; the unit diagonal elements of L are not stored.
//   // *
//   // *  LDA     (input) INTEGER
//   // *          The leading dimension of the array A.  LDA >= max(1,N).
//   // *
//   // *  IPIV    (output) INTEGER array, dimension (N)
//   // *          The pivot indices that define the permutation matrix P;
//   // *          row i of the matrix was interchanged with row IPIV(i).
//   // *
//   // *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
//   // *          On entry, the N-by-NRHS matrix of right hand side matrix B.
//   // *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
//   // *
//   // *  LDB     (input) INTEGER
//   // *          The leading dimension of the array B.  LDB >= max(1,N).
//   // *
//   // *  INFO    (output) INTEGER
//   // *          = 0:  successful exit
//   // *          < 0:  if INFO = -i, the i-th argument had an illegal value
//   // *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
//   // *                has been completed, but the factor U is exactly
//   // *                singular, so the solution could not be computed.
//   // *
//   // *  =====================================================================
//   // 
//   // *     .. Executable Statements ..
//   // *
//   // *     Test the input parameters.
//   INFO = 0;
//   if (N < 0) 
//   {
//     INFO = -1;
//     std::cout << "PairLS::DGESV Illegal value of parameter N = " << N << std::endl;  
//   } 
//   else if (NRHS < 0) 
//   {
//     INFO = -2;
//     std::cout << "PairLS::DGESV Illegal value of parameter NRHS = " << NRHS << std::endl;  
//   } 
//   else if (LDA < MAX(1, N)) 
//   {
//     INFO = -4;
//     std::cout << "PairLS::DGESV Illegal values of parameters LDA = " << LDA << "and/or N = " << N << std::endl;  
//   } 
//   else if (LDB < MAX(1, N)) 
//   {
//     INFO = -7;
//     std::cout << "PairLS::DGESV Illegal values of parameters LDB = " << LDB << "and/or N = " << N << std::endl;  
//   }

//   if (INFO != 0) 
//   {
//       // XERBLA("DGESV ", -INFO);
//       return;
//   }
//   //*
//   //*     Compute the LU factorization of A.
//   //*
//   DGETRF(N, N, A, LDA, IPIV, INFO);
//   if (INFO == 0) 
//   {
//       //*
//       //*        Solve the system A*X = B, overwriting B with X.
//       //*
//       DGETRS("No transpose", N, NRHS, A, LDA, IPIV, B, LDB, INFO);
//   }
//   return;
// }

// // dgetrf.f
// void PairLS::DGETRF(int M, int N, double **A, int LDA, int *IPIV, int *INFO)
// {
//   // *  DGETRF computes an LU factorization of a general M-by-N matrix A
//   // *  using partial pivoting with row interchanges.
//   // *
//   // *  The factorization has the form
//   // *     A = P * L * U
//   // *  where P is a permutation matrix, L is lower triangular with unit
//   // *  diagonal elements (lower trapezoidal if m > n), and U is upper
//   // *  triangular (upper trapezoidal if m < n).
//   // *
//   // *  This is the right-looking Level 3 BLAS version of the algorithm.
//   // *
//   // *  Arguments
//   // *  =========
//   // *
//   // *  M       (input) INTEGER
//   // *          The number of rows of the matrix A.  M >= 0.
//   // *
//   // *  N       (input) INTEGER
//   // *          The number of columns of the matrix A.  N >= 0.
//   // *
//   // *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
//   // *          On entry, the M-by-N matrix to be factored.
//   // *          On exit, the factors L and U from the factorization
//   // *          A = P*L*U; the unit diagonal elements of L are not stored.
//   // *
//   // *  LDA     (input) INTEGER
//   // *          The leading dimension of the array A.  LDA >= max(1,M).
//   // *
//   // *  IPIV    (output) INTEGER array, dimension (min(M,N))
//   // *          The pivot indices; for 1 <= i <= min(M,N), row i of the
//   // *          matrix was interchanged with row IPIV(i).
//   // *
//   // *  INFO    (output) INTEGER
//   // *          = 0:  successful exit
//   // *          < 0:  if INFO = -i, the i-th argument had an illegal value
//   // *          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
//   // *                has been completed, but the factor U is exactly
//   // *                singular, and division by zero will occur if it is used
//   // *                to solve a system of equations.
//   // *
//   // *  =====================================================================

//   //*     .. Parameters ..
//   double const ONE = 1.0;

//   //*     ..
//   //*     .. Local Scalars ..
//   int I, IINFO, J, JB, NB;
//   //*     .. Additional Local Scalars for C-like indexing of arrays
//   int II, JJ; 
//   //*     ..
//   //*     .. External Subroutines ..
//   // EXTERNAL DGEMM, DGETF2, DLASWP, DTRSM, XERBLA;
//   //*     ..
//   //*     .. External Functions ..
//   // int ILAENV;
//   // EXTERNAL ILAENV;
//   //*     ..
//   //*     .. Intrinsic Functions ..
//   // INTRINSIC MAX, MIN;
//   //*     ..
//   //*     .. Executable Statements ..
//   //*
//   //*     Test the input parameters.
//   //*
//   INFO = 0;
//   if (M < 0) 
//   {
//     INFO = -1;
//     std::cout << "PairLS::DGETRF Illegal value of parameter M = " << M << std::endl;  
//   } 
//   else if (N < 0) 
//   {
//     INFO = -2;
//     std::cout << "PairLS::DGETRF Illegal value of parameter N = " << N << std::endl;  
//   } 
//   else if (LDA < MAX(1, M)) 
//   {
//     INFO = -4;
//     std::cout << "PairLS::DGETRF Illegal values of parameters LDA = " << LDA << "and/or M = " << M << std::endl;  
//   }

//   if (INFO != 0) 
//   {
//     // XERBLA("DGETRF", -INFO);
//     return;
//   }
//   //*
//   //*     Quick return if possible
//   //*
//   if (M == 0 || N == 0) return;
//   //*
//   //*     Determine the block size for this environment.
//   //*
//   NB = ILAENV(1, "DGETRF", " ", M, N, -1, -1); // check what this means
//   if (NB <= 1 || NB >= MIN(M, N)) 
//   {
//     //*
//     //*        Use unblocked code.
//     //*
//     DGETF2(M, N, A, LDA, IPIV, INFO);
//   } 
//   else
//   {
//     //*
//     //*        Use blocked code.
//     //*
//     for (J = 1, J <= MIN(M, N), J += NB)    // decide what to use
//     // for (JJ = 1, JJ <= MIN(M, N), JJ += NB)    // decide what to use
//     // for (J = 0, J < MIN(M, N), J += NB)  // decide what to use
//     {
//       // J = JJ - 1;
//       JB = MIN(MIN(M, N) - J + 1, NB);  // decide what to use
//       // JB = MIN(MIN(M, N) - (J+1) + 1, NB); // decide what to use
//       //*
//       //*           Factor diagonal and subdiagonal blocks and test for exact
//       //*           singularity.
//       //*
//       DGETF2(M - J + 1, JB, A[J-1][J-1], LDA, IPIV[0:J-1], IINFO); // decide what to use
//       // DGETF2(M - J + 1, JB, A(J, J), LDA, IPIV(J), IINFO);        // decide what to use
//       // DGETF2(M - (J+1) + 1, JB, A[0:J][0:J], LDA, IPIV[0:J], IINFO); // decide what to use
//       //*
//       //*           Adjust INFO and the pivot indices.
//       //*
//       if (INFO == 0 && IINFO > 0) INFO = IINFO + J - 1;
//       for(I = J, I <= MIN(M, J + JB - 1), I++)
//       {
//         IPIV[I-1] = J - 1 + IPIV[I-1];     // decide what to use
//         // IPIV[I] = (J+1) - 1 + IPIV[I];   // decide what to use
//       }
//       //*
//       //*           Apply interchanges to columns 1:J-1.
//       //*
//       DLASWP(J - 1, A, LDA, J, J + JB - 1, IPIV, 1);                 // decide what to use
//       // DLASWP((J+1) - 1, A, LDA, (J+1), (J+1) + JB - 1, IPIV, 1);  // decide what to use
//       //*
//       if (J + JB <= N)         // decide what to use
//       // if (J+1 + JB <= N)    // decide what to use
//       {
//         //*
//         //*              Apply interchanges to columns J+JB:N.
//         //*
//         DLASWP(N - J - JB + 1, A[1][J + JB], LDA, J, J + JB - 1,  IPIV, 1);      // decide what to use 
//         // DLASWP(N - (J+1) - JB + 1, A[1][(J+1) + JB], LDA, J, J + JB - 1,  IPIV, 1); // decide what to use
//         //*
//         //*              Compute block row of U.
//         //*
//         // DTRSM("Left", "Lower", "No transpose", "Unit", JB, N - J - JB + 1, ONE, A(J, J), LDA, A(J, J + JB),                     LDA);
//         DTRSM("Left", "Lower", "No transpose", "Unit", JB, N - (J+1) - JB + 1, ONE, A[0:J][0:J], LDA, A(J, J + JB),                     LDA);
//         if (J+1 + JB <= M) 
//         {
//           //*
//           //*                 Update trailing submatrix.
//           //*
//           DGEMM("No transpose", "No transpose", M - J - JB + 1,    N - J - JB + 1, JB, -ONE, A(J + JB, J), LDA,                        A(J, J + JB), LDA, ONE, A(J + JB, J + JB),                        LDA);
//         }
//       }
          
//     }
//   }
//   return;

// }

// // ilaenv.f
// int PairLS::ILAENV()
// {

// }

// // ieeeck.f
// int PairLS::IEEECK()
// {

// }

// // lsame.f
// bool PairLS::LSAME()
// {

// }



// // dgetf2.f
// void PairLS::DGETF2(int M, int N, double **A, int LDA, int *IPIV, int *INFO)
// {

// }

// // dlaswp.f
// void PairLS::DLASWP(int N, double **A, int LDA, int K1, int K2, int *IPIV, int INCX)
// {
//   //*
//   //*  -- LAPACK auxiliary routine (version 2.0) --
//   //*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//   //*     Courant Institute, Argonne National Lab, and Rice University
//   //*     October 31, 1992
//   //*
//   //*     .. Scalar Arguments ..
//   // int INCX, K1, K2, LDA, N;
//   //*     ..
//   //*     .. Array Arguments ..
//   // int IPIV[*];
//   // double A(LDA, *);
//   //*     ..
//   //*
//   //*  Purpose
//   //*  =======
//   //*
//   //*  DLASWP performs a series of row interchanges on the matrix A.
//   //*  One row interchange is initiated for each of rows K1 through K2 of A.
//   //*
//   //*  Arguments
//   //*  =========
//   //*
//   //*  N       (input) INTEGER
//   //*          The number of columns of the matrix A.
//   //*
//   //*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
//   //*          On entry, the matrix of column dimension N to which the row
//   //*          interchanges will be applied.
//   //*          On exit, the permuted matrix.
//   //*
//   //*  LDA     (input) INTEGER
//   //*          The leading dimension of the array A.
//   //*
//   //*  K1      (input) INTEGER
//   //*          The first element of IPIV for which a row interchange will
//   //*          be done.
//   //*
//   //*  K2      (input) INTEGER
//   //*          The last element of IPIV for which a row interchange will
//   //*          be done.
//   //*
//   //*  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
//   //*          The vector of pivot indices.  Only the elements in positions
//   //*          K1 through K2 of IPIV are accessed.
//   //*          IPIV(K) = L implies rows K and L are to be interchanged.
//   //*
//   //*  INCX    (input) INTEGER
//   //*          The increment between successive values of IPIV.  If IPIV
//   //*          is negative, the pivots are applied in reverse order.
//   //*
//   //* =====================================================================
//   //*
//   //*     .. Local Scalars ..
//   int I, IP, IX;
//   //*     ..
//   //*     .. External Subroutines ..
//   // EXTERNAL DSWAP;
//   //*     ..
//   //*     .. Executable Statements ..
//   //*
//   //*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
//   //*
//   if (INCX == 0) return;
//   if (INCX > 0) 
//   {
//     IX = K1;
//   } 
//   else
//   {
//     IX = 1 + (1 - K2)*INCX;
//   }

//   if (INCX == 1) 
//   {
//     for(I = K1; I <= K2; I++)
//     {
//       IP = IPIV[I-1];
//       if (IP != I) 
//       {
//         // DSWAP(N, A(I, 1), LDA, A(IP, 1), LDA);
//         DSWAP(N, A[0:I][0], LDA, A(IP, 1), LDA);
//       }
//     }
//   } 
//   else if (INCX > 1) 
//   {
//     for(I = K1; I <= K2; I++)
//     {
//       IP = IPIV[IX];
//       if (IP != I) 
//       {
//         DSWAP(N, A(I, 1), LDA, A(IP, 1), LDA);
//       }
//       IX = IX + INCX;
//     }
//   } 
//   else if (INCX < 0) 
//   {
//       for(I = K2, K1, I -= 1)
//       {
//         IP = IPIV[IX];
//         if (IP != I) 
//         {
//           DSWAP(N, A(I, 1), LDA, A(IP, 1), LDA);
//         }
//         IX = IX + INCX;
//       }
//   }
//   return;
// }

// void PairLS::DSWAP(int N, double **A, int LDA, int K1, int K2, int *IPIV, int INCX)
// {

// }

// // dtrsm.f
// void PairLS::DTRSM()
// {

// }


// // dgemm.f
// void PairLS::DGEMM()
// {

// }

// // dgetrs.f
// void PairLS::DGETRS()
// {

// }

// void PairLS::DGER()
// {
  
// }

// // idamax.f
// int PairLS::idamax()
// {

// }

// // dscal.f
// int PairLS::dscal()
// {

// }
