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
   Contributing author: Andrew Jewett (jewett.aij at g mail)
------------------------------------------------------------------------- */
#ifdef DIHEDRAL_CLASS

DihedralStyle(table,DihedralTable)

#else

#ifndef LMP_DIHEDRAL_TABLE_H
#define LMP_DIHEDRAL_TABLE_H
#include <cassert>
#include <cmath>
#include "domain.h"
#include "dihedral.h"
using namespace std;


namespace LAMMPS_NS {



class DihedralTable : public Dihedral {
 public:
  DihedralTable(class LAMMPS *);
  virtual ~DihedralTable();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int type, int i1, int i2, int i3, int i4);

 protected:
  int tabstyle,tablength;
  // double *phi0;       <- equilibrium angles not supported

  struct Table {
    int ninput;
    //double phi0;      <-equilibrium angles not supported
    bool f_unspecified;
    bool use_degrees;
    double *phifile,*efile,*ffile;
    double *e2file,*f2file;
    double delta,invdelta,deltasq6;
    double *phi,*e,*de,*f,*df,*e2,*f2;
  };

  int ntables;
  Table *tables;
  int *tabindex;
  
  void allocate();
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
	f = tb->f[i] + b * tb->df[i]; //<--works even if tb->f_unspecified==true
	break;
      case SPLINE:
	a = 1.0 - b;
	u = a * tb->e[i] + b * tb->e[ip1] + 
	  ((a*a*a-a)*tb->e2[i] + (b*b*b-b)*tb->e2[ip1]) * 
	  tb->deltasq6;
	if (tb->f_unspecified)
	  //Formula below taken from equation3.3.5 of "numerical recipes in c"
	  //"f"=-derivative of e with respect to x (or "phi" in this case)
	  f = (tb->e[i]-tb->e[ip1])*tb->invdelta +
	    ((3.0*a*a-1.0)*tb->e2[i]+(1.0-3.0*b*b)*tb->e2[ip1])*tb->delta/6.0;
	else
	  f = a * tb->f[i] + b * tb->f[ip1] + 
	    ((a*a*a-a)*tb->f2[i] + (b*b*b-b)*tb->f2[ip1]) * 
	    tb->deltasq6;
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


  // Pre-allocated strings to store file names for debugging splines.  (One day
  // I would really like to rewrite everything and use C++ strings instead.)
  static const int MAXLINE=2048;
  char checkU_fname[MAXLINE];
  char checkF_fname[MAXLINE];

}; //class DihedralTable













// ------------------------------------------------------------------------
// The following auxiliary functions were left out of the 
// DihedralTable class either because they require template parameters, 
// or because they have nothing to do with dihedral angles.
// ------------------------------------------------------------------------

namespace DIHEDRAL_TABLE_NS {

static const double PI    = 3.1415926535897931;
static const double TWOPI = 6.2831853071795862;

// Determine the array of "y2" parameters of a cyclic spline from its control
// points at positions x[] and y[]. (spline() must be invoked before splint())
// The x[] positions should be sorted in order and not exceed period.
void    cyc_spline(double const *xa, double const *ya, int n, 
                   double period, double *y2a);

// Evaluate a cyclic spline at position x with n control points at xa[], ya[],
// (The y2a array must be calculated using cyc_spline() above in advance.)
// x (and all the xa[] positions) should lie in the range from 0 to period.
// (Typically period = 2*PI, but this is optional.)
double  cyc_splint(double const *xa, double const *ya, double const *y2a, 
                   int n, double period, double x);

// Evaluate the deriviative of a cyclic spline at position x:
double cyc_splintD(double const *xa, double const *ya, double const *y2a, 
                   int n, double period, double x);

// -----------------------------------------------------------
// ----  some simple vector operations are defined below. ----
// -----------------------------------------------------------

//  --- g_dim ---   As elsewhere in the LAMMPS code, coordinates here are 
// represented as entries in an array, not as named variables "x" "y" "z".
// (I like this style.)  In this spirit, the vector operations here are 
// defined for vectors of arbitrary size.  For this to work, the number 
// of dimensions, "g_dim", must be known at compile time:
const int g_dim = 3;
// In LAMMPS at least, this constant is always 3, and is only used inside 
// the dihedral code here.  (It should not conflict with 2-D simulations.)
// Note: Compiler optimizations should eliminate any performance overhead
//       associated with loops like "for (int i=0; i<g_dim; i++)"
// If having a constant named "g_dim" is confusing people, I suppose
// we can replace it with "3".  Note: the supplemental file 
// "nd/dihedral_table_nd_mod.h" shows how to generalize the dihedral 
// code in higher dimensions.

template<class _Real>
inline _Real
DotProduct(_Real const *A, _Real const *B)
{
  _Real AdotB = 0.0;
  for (int d=0; d < g_dim; ++d)
    AdotB += A[d]*B[d];
  return AdotB;
}

// Normalize() divides the components of the vector "v" by it's length.
// Normalize() silently ignores divide-by-zero errors but does not
// crash.  (If "v" has length 0, then we replace v with the unit vector in 
// an arbitrary direction,(1,0,...).)
// It returns the length of v (useful for checking if the operation succeeded).
template<class _Real>
inline _Real
Normalize(_Real *v)
{
  _Real length = sqrt(DotProduct(v,v));
  if (length != 0.0)
  {
    _Real one_over_length = 1.0 / length;
    for (int d=0; d < g_dim; ++d) 
      v[d] *= one_over_length;
  }
  else {
    v[0] = 1.0;
    for (int d=1; d < g_dim; ++d) 
      v[d] = 0.0;
  }
  return length;
}


// CrossProduct(A,B,dest) computes the cross-product (A x B) 
// and stores the result in "dest".
template<class _Real>
inline void
CrossProduct(_Real const *A, _Real const *B, _Real *dest)
{
  dest[0] = A[1]*B[2] - A[2]*B[1];
  dest[1] = A[2]*B[0] - A[0]*B[2];
  dest[2] = A[0]*B[1] - A[1]*B[0];
}


// --------------------------------------------
// ------- Calculate the dihedral angle -------
// --------------------------------------------

inline double Phi(double const *x1, //array holding x,y,z coords atom 1
                  double const *x2, // :       :      :      :        2
                  double const *x3, // :       :      :      :        3
                  double const *x4, // :       :      :      :        4
		  Domain *domain, //<-periodic boundary information
                  // The following arrays are of doubles with g_dim elements.
                  // (g_dim is a constant known at compile time, usually 3).
                  // Their contents is calculated by this function.
                  // Space for these vectors must be allocated in advance.
                  // (This is not hidden internally because these vectors
                  //  may be needed outside the function, later on.)
                  double *vb12, // will store x2-x1
                  double *vb23, // will store x3-x2
                  double *vb34, // will store x4-x3
                  double *n123, // will store normal to plane x1,x2,x3
                  double *n234) // will store normal to plane x2,x3,x4
{

  for (int d=0; d < g_dim; ++d) {
    vb12[d] = x2[d] - x1[d]; // 1st bond
    vb23[d] = x3[d] - x2[d]; // 2nd bond
    vb34[d] = x4[d] - x3[d]; // 3rd bond
  }

  //Consider periodic boundary conditions:
  domain->minimum_image(vb12[0],vb12[1],vb12[2]);
  domain->minimum_image(vb23[0],vb23[1],vb23[2]);
  domain->minimum_image(vb34[0],vb34[1],vb34[2]);

  //--- Compute the normal to the planes formed by atoms 1,2,3 and 2,3,4 ---

  CrossProduct(vb12, vb23, n123);        // <- n123=vb12 x vb23
  CrossProduct(vb34, vb23, n234);        // <- n234=vb34 x vb23

  Normalize(n123);
  Normalize(n234);

  double cos_phi = -DotProduct(n123, n234);

  if (cos_phi > 1.0)
    cos_phi = 1.0;
  else if (cos_phi < -1.0)
    cos_phi = -1.0;

  double phi = acos(cos_phi);

  if (DotProduct(n123, vb34) > 0.0) {
    phi = -phi;   //(Note: Negative dihedral angles are possible only in 3-D.)
    phi += TWOPI; //<- This insure phi is always in the range 0 to 2*PI
  }
  return phi;
} // DihedralTable::Phi()


} // namespace DIHEDRAL_TABLE_NS

} // namespace LAMMPS_NS


#endif //#ifndef LMP_DIHEDRAL_TABLE_H
#endif //#ifdef DIHEDRAL_CLASS ... #else
