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

namespace LAMMPS_NS {


class PairLS : public Pair {
 public:
  friend class FixSemiGrandCanonicalMC;   // Alex Stukowski option

  // public variables so USER-ATC package can access them

  double cutmax;

  // potentials as array data

  int nrho,nr;
  int nfrho,nrhor,nz2r;
  double **frho,**rhor,**z2r;
  int *type2frho,**type2rhor,**type2z2r;

  // potentials in spline form used for force computation

  double dr,rdr,drho,rdrho,rhomax;
  double ***rhor_spline,***frho_spline,***z2r_spline;

  PairLS(class LAMMPS *);
  virtual ~PairLS();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double single(int, int, int, int, double, double, double, double &);
  virtual void *extract(const char *, int &);

  // virtual int pack_forward_comm(int, int *, double *, int, int *);
  // virtual void unpack_forward_comm(int, int, double *);
  // int pack_reverse_comm(int, int, double *);
  // void unpack_reverse_comm(int, int *, double *);
  double memory_usage();
  // void swap_ls(double *, double **);

 protected:
 
  const int  n_mark_at=10;
  const int  max_at=100000;
  const int  max_at_def=100;
  const int  max_neighb=200;
  const int  max_pair_at=max_at*max_neighb;
  const int  max_pair_at_def=max_at_def*max_neighb;
  const int  max_list=max_pair_at+max_at;
  const int  max_cell=100000;
  const int  max_at_group=max_at/5;
  const int  max_group=10000;
  const int  max_p=1;
  const int  max_hole= max_p/100;
  const int  max_seat= max_at;

  const int mf3=3;            // max_f3
  const int mfi=30;           // max_sp_fi
  const int mro=25;           // max_sp_ro
  const int memb=10;          // max_sp_emb
  const int mf=15;            // max_sp_f
  const int mg=15;            // max_sp_g
  const int mi=6;             // max_sort_at
 
  bool*1 if_g3_pot, if_g4_pot, if_F2_pot, if_gp0_pot;
  bool*1 if_diag;
  int*4 n_f3, n_sort;
  int*4 n_sp_fi, n_sp_ro, n_sp_emb, n_sp_f, n_sp_g;

  // arrays are inverted in comparison as they declared in common blocks in pot_ls.h
  double shag_sp_fi[mi][mi], shag_sp_ro[mi][mi], shag_sp_emb[mi], shag_sp_f[mi][mi], shag_sp_g;
  
  double R_sp_fi[mi][mi][mfi], R_sp_ro[mi][mi][mro], R_sp_emb[mi][memb], R_sp_f[mi][mi][mf], R_sp_g[mg]; 
  double a_sp_fi[mi][mi][mfi], b_sp_fi[mi][mi][mfi], c_sp_fi[mi][mi][mfi], d_sp_fi[mi][mi][mfi];
  double a_sp_ro, b_sp_ro, c_sp_ro, d_sp_ro, a_sp_emb, b_sp_emb, c_sp_emb, d_sp_emb, a_sp_f3, b_sp_f3, c_sp_f3, d_sp_f3,     a_sp_g3, b_sp_g3, c_sp_g3, d_sp_g3,       a_sp_f4, b_sp_f4, c_sp_f4, d_sp_f4,       a_sp_g4, b_sp_g4, c_sp_g4, d_sp_g4;
  double fip_rmin;
  double z_ion, c_ZBL, d_ZBL, zz_ZBL, a_ZBL, e0_ZBL;
  double Rmin_fi_ZBL, c_fi_ZBL;
  double Rc_fi, Rc_f;

  virtual void read_file(char *);
  virtual void file2array();

  // int nmax;                   // allocated size of per-atom arrays
  // double cutforcesq;
  // double **scale;
  // bigint embedstep;           // timestep, the embedding term was computed

  // // per-atom arrays

  // double *rho,*fp;
  // int *numforce;

  // // potentials as file data

  // int *map;                   // which element each atom type maps to

  // struct Funcfl {
  //   char *file;
  //   int nrho,nr;
  //   double drho,dr,cut,mass;
  //   double *frho,*rhor,*zr;
  // };
  // Funcfl *funcfl;
  // int nfuncfl;

  // struct Setfl {
  //   char **elements;
  //   int nelements,nrho,nr;
  //   double drho,dr,cut;
  //   double *mass;
  //   double **frho,**rhor,***z2r;
  // };
  // Setfl *setfl;

  // struct Fs {
  //   char **elements;
  //   int nelements,nrho,nr;
  //   double drho,dr,cut;
  //   double *mass;
  //   double **frho,***rhor,***z2r;
  // };
  // Fs *fs;

  // virtual void allocate();
  // virtual void array2spline();
  // void interpolate(int, double, double *, double **);

  // virtual void read_file(char *);
  // virtual void file2array();
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
