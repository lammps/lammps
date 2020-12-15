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

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

PairLS::PairLS(LAMMPS *lmp) : Pair(lmp)
{

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

  // memory allocation in heap for arrays with variables and spline coefficients
  // maybe it is not nessessary to allocate all this memory
  // shag_sp_fi  = new double[mi][mi];
  // shag_sp_ro  = new double[mi][mi];
  // shag_sp_emb = new double[mi]; 
  // shag_sp_f   = new double[mi][mi];
  // shag_sp_g   = new double; 

  // R_sp_fi  = new double[mfi][mi][mi]; 
  // R_sp_ro  = new double[mfi][mi][mi]; 
  // R_sp_emb = new double[memb][mi]; 
  // R_sp_f   = new double[mf][mi][mi]; 
  // R_sp_g   = new double[mg];

  // a_sp_fi = new double[mfi][mi][mi]; 
  // b_sp_fi = new double[mfi][mi][mi];
  // c_sp_fi = new double[mfi][mi][mi]; 
  // d_sp_fi = new double[mfi][mi][mi];

  // a_sp_ro = new double[mro][mi][mi];
  // b_sp_ro = new double[mro][mi][mi];
  // c_sp_ro = new double[mro][mi][mi];
  // d_sp_ro = new double[mro][mi][mi];
  
  // a_sp_emb = new double[memb][mi];
  // b_sp_emb = new double[memb][mi];
  // c_sp_emb = new double[memb][mi];
  // d_sp_emb = new double[memb][mi];
  
  // a_sp_f3 = new double[mf][mf3][mi][mi];
  // b_sp_f3 = new double[mf][mf3][mi][mi];
  // c_sp_f3 = new double[mf][mf3][mi][mi];
  // d_sp_f3 = new double[mf][mf3][mi][mi];
  
  // a_sp_g3 = new double[mg][mf3][mf3][mi];
  // b_sp_g3 = new double[mg][mf3][mf3][mi];
  // c_sp_g3 = new double[mg][mf3][mf3][mi];
  // d_sp_g3 = new double[mg][mf3][mf3][mi]; 
  
  // a_sp_f4 = new double[mf][mi][mi];
  // b_sp_f4 = new double[mf][mi][mi];
  // c_sp_f4 = new double[mf][mi][mi];
  // d_sp_f4 = new double[mf][mi][mi];
  
  // a_sp_g4 = new double[mi][mi];
  // b_sp_g4 = new double[mi][mi];
  // c_sp_g4 = new double[mi][mi];
  // d_sp_g4 = new double[mi][mi];
  
  
  // fip_rmin = new double[mi][mi];
  
  
  // z_ion  = new double[mi]; 
  // c_ZBL  = new double[4]; 
  // d_ZBL  = new double[4]; 
  // zz_ZBL = new double[mi][mi]; 
  // a_ZBL  = new double[mi][mi]; 
  // e0_ZBL = new double[mi][mi];
  
  
  // Rmin_fi_ZBL = new double[mi][mi]; 
  // c_fi_ZBL    = new double[6][mi][mi];
  
  // Rc_fi = new double; 
  // Rc_f  = new double;

  n_sort = atom->ntypes;


  // EAM inheritance
  // restartinfo = 0;
  // manybody_flag = 1;
  // embedstep = -1;
  // unit_convert_flag = utils::get_supported_conversions(utils::ENERGY);

  // nmax = 0;
  // rho = nullptr;
  // fp = nullptr;
  // numforce = nullptr;
  // map = nullptr;
  // type2frho = nullptr;

  // nfuncfl = 0;
  // funcfl = nullptr;

  // setfl = nullptr;
  // fs = nullptr;

  // frho = nullptr;
  // rhor = nullptr;
  // z2r = nullptr;
  // scale = nullptr;

  // frho_spline = nullptr;
  // rhor_spline = nullptr;
  // z2r_spline = nullptr;

  // set comm size needed by this Pair

  comm_forward = 1;
  comm_reverse = 1;
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

  ev_init(eflag,vflag);


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
  memory->create(d_ZBL,4,"PairLS:d_ZBL");
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
   (For future implementations. Coeffs for each pair of atoms are the names of files with potential functions should be written separately (example for 3-component system):
   pair_coeff 1 1 pot_1
   pair_coeff 2 2 pot_2
   pair_coeff 3 3 pot_3
   pair_coeff 1 2 pot_1_2
   pair_coeff 1 3 pot_1_3
   pair_coeff 2 3 pot_2_3)

  coeffs for each pair of atoms are the names of files with potential functions should be in a row in the following order (example for 3-component system):
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
      setflag[i][i] = 1;
      // print_pot_arrays_is(i); // check if potential was correctly read
      // par2pot_is(i);
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
          // par2pot_is1_is2(i,j);
          // par2pot_is1_is2(j,i);
          std::cout << "!!!!! PairLS debug mode !!!!! " << " Start reading cross potential for atoms " << i << " and " << j << std::endl;
          std::cout << "!!!!! PairLS debug mode !!!!! " << " The arg with potential name is " << arg[ij-1] << std::endl;
          r_pot_ls_is1_is2(arg[ij-1], i, j);
          std::cout << "!!!!! PairLS debug mode !!!!! " << " End reading cross potential for atoms " << i << " and " << j << std::endl;
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
        d_ZBL[n] = reader.next_double();
        std::cout << " d_ZBL[" << n << "] = " << d_ZBL[n]<< std::endl;
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



// void PairLS::par2pot_is(int is)
// {
  
// }



// void PairLS::par2pot_is1_is2(int is1, int is2)
// {
  
// }


// utility functions

// bool PairLS::string2bool(const std::string & v)
// {
//     return !v.empty () &&
//         (strcasecmp (v.c_str (), "true") == 0 ||
//          atoi (v.c_str ()) != 0);
// }

