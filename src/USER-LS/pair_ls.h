/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

 
#ifdef PAIR_CLASS

PairStyle(ls,PairLS)

#else

#ifndef LMP_PAIR_LS_H
#define LMP_PAIR_LS_H

#include "pair.h"

/******************************************************************************
* BLAS and LAPACK definitions
******************************************************************************/
#ifdef MKL
#include "mkl.h"
#define dgesv_ dgesv
#else
// #include <lapacke.h>
// #define dgesv_ LAPACKE_dgesv_work
// extern "C"
// {
extern "C" void dgesv_(int *, int *, double *, int *, int *, double *, int *, int *);
// extern void dgesv_(int *, int *, double *, int *, int *, double *, int *, int *);
// }
#endif

namespace LAMMPS_NS {

class PairLS : public Pair {
 public:

  PairLS(class LAMMPS *);
  virtual ~PairLS();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  void init_style();
  double init_one(int, int);

  // maybe not needed public methods
  // double single(int, int, int, int, double, double, double, double &);
  // virtual void *extract(const char *, int &);

  // virtual int pack_forward_comm(int, int *, double *, int, int *);
  // virtual void unpack_forward_comm(int, int, double *);
  // int pack_reverse_comm(int, int, double *);
  // void unpack_reverse_comm(int, int *, double *);
  // double memory_usage();



 protected:
 
  int *map;                   // which element each atom type maps to

  // Begin max_at.h (may be it is not needed in the LAMMPS implementation)
  const int  n_mark_at=10;
  const long  max_at=100000;
  const int  max_at_def=100;
  const int  max_neighb=200;
  const long  max_pair_at=max_at*max_neighb;
  const long  max_pair_at_def=max_at_def*max_neighb;
  const long  max_list=max_pair_at+max_at;
  const long  max_cell=100000;
  const long  max_at_group=max_at/5;
  const int  max_group=10000;
  const int  max_p=1;
  const int  max_hole=max_p/100;
  const long  max_seat=max_at;
  // End max_at.h

  // Begin pot_ls.h
  const int mf3=4;            // max_f3
  const int mfi=30;           // max_sp_fi
  const int mro=25;           // max_sp_ro
  const int memb=10;          // max_sp_emb
  const int mf=15;            // max_sp_f
  const int mg=15;            // max_sp_g
  const int mi=6;             // max_sort_at (actually is equal to mi-1 since C-arrays are started form 0 while the atom indexes are started from 1)



 
 // common/abcd_sp/
 // pointers to arrays with spline coefficients and data about steps 
 
  double **shag_sp_fi, **shag_sp_ro, *shag_sp_emb, **shag_sp_f, shag_sp_g; 
  double ***R_sp_fi, ***R_sp_ro, **R_sp_emb, ***R_sp_f, *R_sp_g; 
  // spline coefficients determining functions for pair interactions: the first term in the equation for total energy fi(R_ij) 
  double ***a_sp_fi, ***b_sp_fi, ***c_sp_fi, ***d_sp_fi;     
  // spline coefficients determining basis functions rho(R_ij)
  double ***a_sp_ro, ***b_sp_ro, ***c_sp_ro, ***d_sp_ro;     
  // spline coefficients determining functions for manybody interactions within the centro-symmetric approximation (CSA): the last term in the equation for total energy F(sum_i!=j(rho_ij)) 
  double **a_sp_emb, **b_sp_emb, **c_sp_emb, **d_sp_emb; 
  // spline coefficients determining basis functions f3(R_ij) in the second term describing three-body interactions in the equation for total energy
  double ****a_sp_f3, ****b_sp_f3, ****c_sp_f3, ****d_sp_f3;  
  // spline coefficients determining expansion coefficients g3(cos(theta_ijk)) for basis functions f3(R_ij) in the third term describing three-body interactions in the equation for total energy    
  double ****a_sp_g3, ****b_sp_g3, ****c_sp_g3, ****d_sp_g3;   
  // spline coefficients determining basis functions f4(R_ij) in the third term describing four-body interactions in the equation for total energy   
  double ***a_sp_f4, ***b_sp_f4, ***c_sp_f4, ***d_sp_f4;
  // spline coefficients determining expansion coefficients g4(cos(theta_ijk),cos(theta_jkl),cos(theta_ilj)) for basis functions f4(R_ij) in the third term describing four-body interactions in the equation for total energy    
  double **a_sp_g4, **b_sp_g4, **c_sp_g4, **d_sp_g4;    
  double **fip_rmin;
  double *z_ion, *c_ZBL, *d_ZBL, **zz_ZBL, **a_ZBL, **e0_ZBL; // maybe it is useful to make c_ZBL and d_ZBL three-dimensional arrays 
  double **Rmin_fi_ZBL, ***c_fi_ZBL;
  double Rc_fi, Rc_f;


 // common/abcd_i/

  bool if_g3_pot;  // if two- and three-body interactions are accurately described by the potential
  bool if_g4_pot;  // if two-, three- and four-body interactions are accurately described by the potential
  bool if_F2_pot;  // if two-body interactions are accurately described by the potential (like within the EAM potential)
  bool *if_gp0_pot;
  // bool *if_diag;   // probably deprecated variable
  int n_sort;     // number of sorts of atoms
  int **n_sp_fi, **n_sp_ro, **n_sp_emb, **n_sp_f, **n_sp_g; // numbers of spline nodes  for differen 
  int *n_f3;    // array with numbers of basis functions for each sort of atom

  // common/md_relax_manager/ 
  int which_pot;
  int i_shift_fixed_boundary;
  bool fixed_boundary_y; 
  bool periodic[3]; 
  bool if_Px_Py_Pz;
  bool if_relax_z_only;
  bool if_calc_pressure;
  bool if_print;
  bool *if_true_i;

  /* 
   Fortran subroutines translated into the void or double functions (class methods) in C++. 
   In the original program md1_temp, all these subrotines use the variables stored in the common blocks which are declared in the file pot_ls.h. 
   In the LAMMPS implementation, the variables from the Fortran common blocks are declared as protected members of PairLS class.
   The protected members of PairLS class can be inherited by any of its child classes.   
  */

  // subroutines for reading potential files
  void r_pot_ls_is(char *, int, double, double);
  void r_pot_ls_is1_is2(char *, int, int, double, double);
  void allocate();

  void par2pot_is(int);
  void par2pot_is1_is2(int, int);
  void smooth_zero_22(double *, double, double, double, double, double, double, double, double);
  void SPL(int, double *, double *, int, double, double, double *, double *, double *);
  void LA30(int, double *, double *, double *, double *, double *, int *);

  // subroutines for calculating energies and forces
  void e_force_fi_emb(int, double *, double **, double *, double *, double *,  double **, int *, int, double *, double *, double *);
  void e_force_g3(int, double *, double **, double *, double *, double *,  double **, int *, int, double *, double *, double *);

  // potential functions and their derivatives
  double fun_fi(double, int, int);
  double funp_fi(double, int, int);
  double funpp_fi(double, int, int);

  double fun_ro(double, int, int);
  double funp_ro(double, int, int);
  double funpp_ro(double, int, int);  

  double fun_emb(double, int);
  double funp_emb(double, int);
  double funpp_emb(double, int); 

  double fun_f3(double, int, int, int);
  double funp_f3(double, int, int, int);
  double funpp_f3(double, int, int, int);

  double fun_g3(double, int, int, int);
  double funp_g3(double, int, int, int);
  double funpp_g3(double, int, int, int);

  double v_ZBL(double, int, int);
  double vp_ZBL(double, int, int);
  double vpp_ZBL(double, int, int);

  double fun_fi_ZBL(double, int, int);
  double funp_fi_ZBL(double, int, int);
  double funpp_fi_ZBL(double, int, int);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Cannot open LS potential file %s

The specified LS potential file cannot be opened.  Check that the
path and name are correct.

E: Invalid LS potential file

UNDOCUMENTED

*/
