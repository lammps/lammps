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

#ifdef DIHEDRAL_CLASS

DihedralStyle(table/cut,DihedralTableCut)

#else

#ifndef LMP_DIHEDRAL_TABLE_CUT_H
#define LMP_DIHEDRAL_TABLE_CUT_H

#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralTableCut : public Dihedral {
 public:
  DihedralTableCut(class LAMMPS *);
  virtual ~DihedralTableCut();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int type, int i1, int i2, int i3, int i4);

 protected:
  double *aat_k,*aat_theta0_1,*aat_theta0_2;

  void allocate();

  int tabstyle,tablength;
  char *checkU_fname;
  char *checkF_fname;

  struct Table {
    int ninput;
    int f_unspecified; // boolean (but MPI does not like type "bool")
    int use_degrees;   // boolean (but MPI does not like type "bool")
    double *phifile,*efile,*ffile;
    double *e2file,*f2file;
    double delta,invdelta,deltasq6;
    double *phi,*e,*de,*f,*df,*e2,*f2;
  };

  int ntables;
  Table *tables;
  int *tabindex;

  void null_table(Table *);
  void free_table(Table *);
  void read_table(Table *, char *, char *);
  void bcast_table(Table *);
  void spline_table(Table *);
  void compute_table(Table *);

  void param_extract(Table *, char *);

  // --------------------------------------------
  // ------------ inline functions --------------
  // --------------------------------------------

  // -----------------------------------------------------------
  //   uf_lookup()
  //   quickly calculate the potential u and force f at angle x,
  //   using the internal tables tb->e and tb->f (evenly spaced)
  // -----------------------------------------------------------
  enum{LINEAR,SPLINE};

  inline void uf_lookup(int type, double x, double &u, double &f)
  {
    Table *tb = &tables[tabindex[type]];
    double x_over_delta = x*tb->invdelta;
    int i = static_cast<int> (x_over_delta);
    double a;
    double b = x_over_delta - i;
    // Apply periodic boundary conditions to indices i and i+1
    if (i >= tablength) i -= tablength;
    int ip1 = i+1; if (ip1 >= tablength) ip1 -= tablength;

    switch(tabstyle) {
      case LINEAR:
        u = tb->e[i] + b * tb->de[i];
        f = -(tb->f[i] + b * tb->df[i]); //<--works even if tb->f_unspecified==true
        break;
      case SPLINE:
        a = 1.0 - b;
        u = a * tb->e[i] + b * tb->e[ip1] +
          ((a*a*a-a)*tb->e2[i] + (b*b*b-b)*tb->e2[ip1]) *
          tb->deltasq6;
        if (tb->f_unspecified)
          //Formula below taken from equation3.3.5 of "numerical recipes in c"
          //"f"=-derivative of e with respect to x (or "phi" in this case)
          f = -((tb->e[i]-tb->e[ip1])*tb->invdelta +
            ((3.0*a*a-1.0)*tb->e2[i]+(1.0-3.0*b*b)*tb->e2[ip1])*tb->delta/6.0);
        else
          f = -(a * tb->f[i] + b * tb->f[ip1] +
            ((a*a*a-a)*tb->f2[i] + (b*b*b-b)*tb->f2[ip1]) *
            tb->deltasq6);
        break;
    } // switch(tabstyle)
  } // uf_lookup()


  // ----------------------------------------------------------
  //    u_lookup()
  //  quickly calculate the potential u at angle x using tb->e
  //-----------------------------------------------------------

  inline void u_lookup(int type, double x, double &u)
  {
    Table *tb = &tables[tabindex[type]];
    int N = tablength;

    //  i = static_cast<int> ((x - tb->lo) * tb->invdelta); <-general version
    double x_over_delta = x*tb->invdelta;
    int    i = static_cast<int> (x_over_delta);
    double b = x_over_delta - i;

    // Apply periodic boundary conditions to indices i and i+1
    if (i >= N) i -= N;
    int ip1 = i+1; if (ip1 >= N) ip1 -= N;

    if (tabstyle == LINEAR) {
      u = tb->e[i] + b * tb->de[i];
    }
    else if (tabstyle == SPLINE) {
      double a = 1.0 - b;
      u = a * tb->e[i] + b * tb->e[ip1] +
        ((a*a*a-a)*tb->e2[i] + (b*b*b-b)*tb->e2[ip1]) *
        tb->deltasq6;
    }
  } // u_lookup()


};

}

#endif
#endif

/* ERROR/WARNING messages:

W: Dihedral problem: %d %ld %d %d %d %d

Conformation of the 4 listed dihedral atoms is extreme; you may want
to check your simulation geometry.

E: Incorrect args for dihedral coefficients

Self-explanatory.  Check the input script or data file.

*/
