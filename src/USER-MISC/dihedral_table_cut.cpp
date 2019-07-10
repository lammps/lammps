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
   Contributing author: K. Michael Salerno (NRL)
   Based on tabulated dihedral (dihedral_table.cpp) by Andrew Jewett
------------------------------------------------------------------------- */

#include "dihedral_table_cut.h"
#include <mpi.h>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>  // IWYU pragma: keep
#include <sstream>  // IWYU pragma: keep

#include "atom.h"
#include "neighbor.h"
#include "update.h"
#include "comm.h"
#include "force.h"
#include "citeme.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;

static const char cite_dihedral_tablecut[] =
  "dihedral_style  table/cut  command:\n\n"
  "@Article{Salerno17,\n"
  " author =  {K. M. Salerno and N. Bernstein},\n"
  " title =   {Persistence Length, End-to-End Distance, and Structure of Coarse-Grained Polymers},\n"
  " journal = {J.~Chem.~Theory Comput.},\n"
  " year =    2018,\n"
  " DOI = 10.1021/acs.jctc.7b01229"
  "}\n\n";

/* ---------------------------------------------------------------------- */

#define TOLERANCE 0.05
#define SMALL     0.0000001

// ------------------------------------------------------------------------
// The following auxiliary functions were left out of the
// DihedralTable class either because they require template parameters,
// or because they have nothing to do with dihedral angles.
// ------------------------------------------------------------------------

// -------------------------------------------------------------------
// ---------    The function was taken verbatim from the    ---------
// ---------    GNU Scientific Library (GSL, version 1.15)   ---------
// -------------------------------------------------------------------

/* Author: Gerard Jungman */
/* for description of method see [Engeln-Mullges + Uhlig, p. 96]
 *
 *      diag[0]  offdiag[0]             0   .....  offdiag[N-1]
 *   offdiag[0]     diag[1]    offdiag[1]   .....
 *            0  offdiag[1]       diag[2]
 *            0           0    offdiag[2]   .....
 *          ...         ...
 * offdiag[N-1]         ...
 *
 */
// -- (A non-symmetric version of this function is also available.) --

enum { //GSL status return codes.
  GSL_FAILURE  = -1,
  GSL_SUCCESS  = 0,
  GSL_ENOMEM   = 8,
  GSL_EZERODIV = 12,
  GSL_EBADLEN  = 19
};


static int solve_cyc_tridiag( const double diag[], size_t d_stride,
                              const double offdiag[], size_t o_stride,
                              const double b[], size_t b_stride,
                              double x[], size_t x_stride,
                              size_t N, bool warn)
{
  int status = GSL_SUCCESS;
  double * delta = (double *) malloc (N * sizeof (double));
  double * gamma = (double *) malloc (N * sizeof (double));
  double * alpha = (double *) malloc (N * sizeof (double));
  double * c = (double *) malloc (N * sizeof (double));
  double * z = (double *) malloc (N * sizeof (double));

  if (delta == 0 || gamma == 0 || alpha == 0 || c == 0 || z == 0) {
    if (warn)
      fprintf(stderr,"Internal Cyclic Spline Error: failed to allocate working space\n");

    if (delta) free(delta);
    if (gamma) free(gamma);
    if (alpha) free(alpha);
    if (c) free(c);
    if (z) free(z);
    return GSL_ENOMEM;
  }
  else
    {
      size_t i, j;
      double sum = 0.0;

      /* factor */

      if (N == 1)
        {
          x[0] = b[0] / diag[0];
          free(delta);
          free(gamma);
          free(alpha);
          free(c);
          free(z);
          return GSL_SUCCESS;
        }

      alpha[0] = diag[0];
      gamma[0] = offdiag[0] / alpha[0];
      delta[0] = offdiag[o_stride * (N-1)] / alpha[0];

      if (alpha[0] == 0) {
        status = GSL_EZERODIV;
      }

      for (i = 1; i < N - 2; i++)
        {
          alpha[i] = diag[d_stride * i] - offdiag[o_stride * (i-1)] * gamma[i - 1];
          gamma[i] = offdiag[o_stride * i] / alpha[i];
          delta[i] = -delta[i - 1] * offdiag[o_stride * (i-1)] / alpha[i];
          if (alpha[i] == 0) {
            status = GSL_EZERODIV;
          }
        }

      for (i = 0; i < N - 2; i++)
        {
          sum += alpha[i] * delta[i] * delta[i];
        }

      alpha[N - 2] = diag[d_stride * (N - 2)] - offdiag[o_stride * (N - 3)] * gamma[N - 3];

      gamma[N - 2] = (offdiag[o_stride * (N - 2)] - offdiag[o_stride * (N - 3)] * delta[N - 3]) / alpha[N - 2];

      alpha[N - 1] = diag[d_stride * (N - 1)] - sum - alpha[(N - 2)] * gamma[N - 2] * gamma[N - 2];

      /* update */
      z[0] = b[0];
      for (i = 1; i < N - 1; i++)
        {
          z[i] = b[b_stride * i] - z[i - 1] * gamma[i - 1];
        }
      sum = 0.0;
      for (i = 0; i < N - 2; i++)
        {
          sum += delta[i] * z[i];
        }
      z[N - 1] = b[b_stride * (N - 1)] - sum - gamma[N - 2] * z[N - 2];
      for (i = 0; i < N; i++)
        {
          c[i] = z[i] / alpha[i];
        }

      /* backsubstitution */
      x[x_stride * (N - 1)] = c[N - 1];
      x[x_stride * (N - 2)] = c[N - 2] - gamma[N - 2] * x[x_stride * (N - 1)];
      if (N >= 3)
        {
          for (i = N - 3, j = 0; j <= N - 3; j++, i--)
            {
              x[x_stride * i] = c[i] - gamma[i] * x[x_stride * (i + 1)] - delta[i] * x[x_stride * (N - 1)];
            }
        }
    }

  free (z);
  free (c);
  free (alpha);
  free (gamma);
  free (delta);

  if ((status == GSL_EZERODIV) && warn)
      fprintf(stderr, "Internal Cyclic Spline Error: Matrix must be positive definite.\n");

  return status;
} //solve_cyc_tridiag()

/* ----------------------------------------------------------------------
   spline and splint routines modified from Numerical Recipes
------------------------------------------------------------------------- */

static int cyc_spline(double const *xa,
                      double const *ya,
                      int n,
                      double period,
                      double *y2a, bool warn)
{

  double *diag    = new double[n];
  double *offdiag = new double[n];
  double *rhs     = new double[n];
  double xa_im1, xa_ip1;

  // In the cyclic case, there are n equations with n unknows.
  // The for loop sets up the equations we need to solve.
  // Later we invoke the GSL tridiagonal matrix solver to solve them.

  for(int i=0; i < n; i++) {

    // I have to lookup xa[i+1] and xa[i-1].  This gets tricky because of
    // periodic boundary conditions.  We handle that now.
    int im1 = i-1;
    if (im1<0) {
      im1 += n;
      xa_im1 = xa[im1] - period;
    }
    else
      xa_im1 = xa[im1];

    int ip1 = i+1;
    if (ip1>=n) {
      ip1 -= n;
      xa_ip1 = xa[ip1] + period;
    }
    else
      xa_ip1 = xa[ip1];

    // Recall that we want to find the y2a[] parameters (there are n of them).
    // To solve for them, we have a linear equation with n unknowns
    // (in the cyclic case that is).  For details, the non-cyclic case is
    // explained in equation 3.3.7 in Numerical Recipes in C, p. 115.
    diag[i]    = (xa_ip1 - xa_im1) / 3.0;
    offdiag[i] = (xa_ip1 - xa[i]) / 6.0;
    rhs[i]     = ((ya[ip1] - ya[i]) / (xa_ip1 - xa[i])) -
                 ((ya[i] - ya[im1]) / (xa[i] - xa_im1));
  }

  // Because this matrix is tridiagonal (and cyclic), we can use the following
  // cheap method to invert it.
  if (solve_cyc_tridiag(diag, 1,
                    offdiag, 1,
                    rhs, 1,
                    y2a, 1,
                    n, warn) != GSL_SUCCESS) {
    if (warn)
      fprintf(stderr,"Error in inverting matrix for splines.\n");

    delete [] diag;
    delete [] offdiag;
    delete [] rhs;
    return 1;
  }
  delete [] diag;
  delete [] offdiag;
  delete [] rhs;
  return 0;
} // cyc_spline()

/* ---------------------------------------------------------------------- */

// cyc_splint(): Evaluates a spline at position x, with n control
//           points located at xa[], ya[], and with parameters y2a[]
//           The xa[] must be monotonically increasing and their
//           range should not exceed period (ie xa[n-1] < xa[0] + period).
//           x must lie in the range:  [(xa[n-1]-period), (xa[0]+period)]
//           "period" is typically 2*PI.
static double cyc_splint(double const *xa,
                         double const *ya,
                         double const *y2a,
                         int n,
                         double period,
                         double x)
{
  int klo = -1;
  int khi = n;
  int k;
  double xlo = xa[n-1] - period;
  double xhi = xa[0] + period;
  while (khi-klo > 1) {
    k = (khi+klo) >> 1; //(k=(khi+klo)/2)
    if (xa[k] > x) {
      khi = k;
      xhi = xa[k];
    }
    else {
      klo = k;
      xlo = xa[k];
    }
  }

  if (khi == n) khi = 0;
  if (klo ==-1) klo = n-1;

  double h = xhi-xlo;
  double a = (xhi-x) / h;
  double b = (x-xlo) / h;
  double y = a*ya[klo] + b*ya[khi] +
    ((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi]) * (h*h)/6.0;

  return y;

} // cyc_splint()


static double cyc_lin(double const *xa,
                      double const *ya,
                      int n,
                      double period,
                      double x)
{
  int klo = -1;
  int khi = n;
  int k;
  double xlo = xa[n-1] - period;
  double xhi = xa[0] + period;
  while (khi-klo > 1) {
    k = (khi+klo) >> 1; //(k=(khi+klo)/2)
    if (xa[k] > x) {
      khi = k;
      xhi = xa[k];
    }
    else {
      klo = k;
      xlo = xa[k];
    }
  }

  if (khi == n) khi = 0;
  if (klo ==-1) klo = n-1;

  double h = xhi-xlo;
  double a = (xhi-x) / h;
  double b = (x-xlo) / h;
  double y = a*ya[klo] + b*ya[khi];

  return y;

} // cyc_lin()




// cyc_splintD(): Evaluate the deriviative of a cyclic spline at position x,
//           with n control points at xa[], ya[], with parameters y2a[].
//           The xa[] must be monotonically increasing and their
//           range should not exceed period (ie xa[n-1] < xa[0] + period).
//           x must lie in the range:  [(xa[n-1]-period), (xa[0]+period)]
//           "period" is typically 2*PI.
static double cyc_splintD(double const *xa,
                          double const *ya,
                          double const *y2a,
                          int n,
                          double period,
                          double x)
{
  int klo = -1;
  int khi = n; // (not n-1)
  int k;
  double xlo = xa[n-1] - period;
  double xhi = xa[0] + period;
  while (khi-klo > 1) {
    k = (khi+klo) >> 1; //(k=(khi+klo)/2)
    if (xa[k] > x) {
      khi = k;
      xhi = xa[k];
    }
    else {
      klo = k;
      xlo = xa[k];
    }
  }

  if (khi == n) khi = 0;
  if (klo ==-1) klo = n-1;

  double yhi = ya[khi];
  double ylo = ya[klo];
  double h = xhi-xlo;
  double g = yhi-ylo;
  double a = (xhi-x) / h;
  double b = (x-xlo) / h;
  // Formula below taken from equation 3.3.5 of "numerical recipes in c"
  // "yD" = the derivative of y
  double yD = g/h - ( (3.0*a*a-1.0)*y2a[klo] - (3.0*b*b-1.0)*y2a[khi] ) * h/6.0;
  // For rerefence: y = a*ylo + b*yhi +
  //                  ((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi]) * (h*h)/6.0;

  return yD;

} // cyc_splintD()


/* ---------------------------------------------------------------------- */

DihedralTableCut::DihedralTableCut(LAMMPS *lmp) : Dihedral(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_dihedral_tablecut);
  ntables = 0;
  tables = NULL;
  checkU_fname = checkF_fname = NULL;
}

/* ---------------------------------------------------------------------- */

DihedralTableCut::~DihedralTableCut()
{
  if (allocated) {
    memory->destroy(aat_k);
    memory->destroy(aat_theta0_1);
    memory->destroy(aat_theta0_2);

    for (int m = 0; m < ntables; m++) free_table(&tables[m]);
    memory->sfree(tables);
    memory->sfree(checkU_fname);
    memory->sfree(checkF_fname);

    memory->destroy(setflag);
    memory->destroy(tabindex);

  }
}

/* ---------------------------------------------------------------------- */

void DihedralTableCut::compute(int eflag, int vflag)
{

  int i1,i2,i3,i4,i,j,k,n,type;
  double edihedral;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
  double fphi,fpphi;
  double r1mag2,r1,r2mag2,r2,r3mag2,r3;
  double sb1,rb1,sb2,rb2,sb3,rb3,c0,r12c1;
  double r12c2,costh12,costh13,costh23,sc1,sc2,s1,s2,c;
  double phi,sinphi,a11,a22,a33,a12,a13,a23,sx1,sx2;
  double sx12,sy1,sy2,sy12,sz1,sz2,sz12;
  double t1,t2,t3,t4;
  double da1,da2;
  double s12,sin2;
  double dcosphidr[4][3],dphidr[4][3],dthetadr[2][4][3];
  double fabcd[4][3];

  edihedral = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist = neighbor->ndihedrallist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < ndihedrallist; n++) {
    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = dihedrallist[n][4];

    // 1st bond

    vb1x = x[i1][0] - x[i2][0];
    vb1y = x[i1][1] - x[i2][1];
    vb1z = x[i1][2] - x[i2][2];

    // 2nd bond

    vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];

    vb2xm = -vb2x;
    vb2ym = -vb2y;
    vb2zm = -vb2z;

    // 3rd bond

    vb3x = x[i4][0] - x[i3][0];
    vb3y = x[i4][1] - x[i3][1];
    vb3z = x[i4][2] - x[i3][2];

    // distances

    r1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
    r1 = sqrt(r1mag2);
    r2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
    r2 = sqrt(r2mag2);
    r3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
    r3 = sqrt(r3mag2);

    sb1 = 1.0/r1mag2;
    rb1 = 1.0/r1;
    sb2 = 1.0/r2mag2;
    rb2 = 1.0/r2;
    sb3 = 1.0/r3mag2;
    rb3 = 1.0/r3;

    c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3;

    // angles

    r12c1 = rb1*rb2;
    r12c2 = rb2*rb3;
    costh12 = (vb1x*vb2x + vb1y*vb2y + vb1z*vb2z) * r12c1;
    costh13 = c0;
    costh23 = (vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z) * r12c2;

    // cos and sin of 2 angles and final c

    sin2 = MAX(1.0 - costh12*costh12,0.0);
    sc1 = sqrt(sin2);
    if (sc1 < SMALL) sc1 = SMALL;
    sc1 = 1.0/sc1;

    sin2 = MAX(1.0 - costh23*costh23,0.0);
    sc2 = sqrt(sin2);
    if (sc2 < SMALL) sc2 = SMALL;
    sc2 = 1.0/sc2;

    s1 = sc1 * sc1;
    s2 = sc2 * sc2;
    s12 = sc1 * sc2;
    c = (c0 + costh12*costh23) * s12;

    // error check

    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
      int me;
      MPI_Comm_rank(world,&me);
      if (screen) {
        char str[128];
        sprintf(str,"Dihedral problem: %d " BIGINT_FORMAT " "
                TAGINT_FORMAT " " TAGINT_FORMAT " "
                TAGINT_FORMAT " " TAGINT_FORMAT,
                me,update->ntimestep,
                atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4]);
        error->warning(FLERR,str,0);
        fprintf(screen,"  1st atom: %d %g %g %g\n",
                me,x[i1][0],x[i1][1],x[i1][2]);
        fprintf(screen,"  2nd atom: %d %g %g %g\n",
                me,x[i2][0],x[i2][1],x[i2][2]);
        fprintf(screen,"  3rd atom: %d %g %g %g\n",
                me,x[i3][0],x[i3][1],x[i3][2]);
        fprintf(screen,"  4th atom: %d %g %g %g\n",
                me,x[i4][0],x[i4][1],x[i4][2]);
      }
    }

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    double phil = acos(c);
    phi = acos(c);

    sinphi = sqrt(1.0 - c*c);
    sinphi = MAX(sinphi,SMALL);

    // n123 = vb1 x vb2

    double n123x = vb1y*vb2z - vb1z*vb2y;
    double n123y = vb1z*vb2x - vb1x*vb2z;
    double n123z = vb1x*vb2y - vb1y*vb2x;
    double n123_dot_vb3 = n123x*vb3x + n123y*vb3y + n123z*vb3z;
    if (n123_dot_vb3 > 0.0) {
      phil = -phil;
      phi = -phi;
      sinphi = -sinphi;
    }

    a11 = -c*sb1*s1;
    a22 = sb2 * (2.0*costh13*s12 - c*(s1+s2));
    a33 = -c*sb3*s2;
    a12 = r12c1 * (costh12*c*s1 + costh23*s12);
    a13 = rb1*rb3*s12;
    a23 = r12c2 * (-costh23*c*s2 - costh12*s12);

    sx1  = a11*vb1x + a12*vb2x + a13*vb3x;
    sx2  = a12*vb1x + a22*vb2x + a23*vb3x;
    sx12 = a13*vb1x + a23*vb2x + a33*vb3x;
    sy1  = a11*vb1y + a12*vb2y + a13*vb3y;
    sy2  = a12*vb1y + a22*vb2y + a23*vb3y;
    sy12 = a13*vb1y + a23*vb2y + a33*vb3y;
    sz1  = a11*vb1z + a12*vb2z + a13*vb3z;
    sz2  = a12*vb1z + a22*vb2z + a23*vb3z;
    sz12 = a13*vb1z + a23*vb2z + a33*vb3z;

    // set up d(cos(phi))/d(r) and dphi/dr arrays

    dcosphidr[0][0] = -sx1;
    dcosphidr[0][1] = -sy1;
    dcosphidr[0][2] = -sz1;
    dcosphidr[1][0] = sx2 + sx1;
    dcosphidr[1][1] = sy2 + sy1;
    dcosphidr[1][2] = sz2 + sz1;
    dcosphidr[2][0] = sx12 - sx2;
    dcosphidr[2][1] = sy12 - sy2;
    dcosphidr[2][2] = sz12 - sz2;
    dcosphidr[3][0] = -sx12;
    dcosphidr[3][1] = -sy12;
    dcosphidr[3][2] = -sz12;

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
        dphidr[i][j] = -dcosphidr[i][j] / sinphi;


    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
        fabcd[i][j] = 0;
    edihedral = 0;


    // set up d(theta)/d(r) array
    // dthetadr(i,j,k) = angle i, atom j, coordinate k

    for (i = 0; i < 2; i++)
      for (j = 0; j < 4; j++)
        for (k = 0; k < 3; k++)
          dthetadr[i][j][k] = 0.0;

    t1 = costh12 / r1mag2;
    t2 = costh23 / r2mag2;
    t3 = costh12 / r2mag2;
    t4 = costh23 / r3mag2;

    // angle12

    dthetadr[0][0][0] = sc1 * ((t1 * vb1x) - (vb2x * r12c1));
    dthetadr[0][0][1] = sc1 * ((t1 * vb1y) - (vb2y * r12c1));
    dthetadr[0][0][2] = sc1 * ((t1 * vb1z) - (vb2z * r12c1));

    dthetadr[0][1][0] = sc1 * ((-t1 * vb1x) + (vb2x * r12c1) +
                               (-t3 * vb2x) + (vb1x * r12c1));
    dthetadr[0][1][1] = sc1 * ((-t1 * vb1y) + (vb2y * r12c1) +
                               (-t3 * vb2y) + (vb1y * r12c1));
    dthetadr[0][1][2] = sc1 * ((-t1 * vb1z) + (vb2z * r12c1) +
                               (-t3 * vb2z) + (vb1z * r12c1));

    dthetadr[0][2][0] = sc1 * ((t3 * vb2x) - (vb1x * r12c1));
    dthetadr[0][2][1] = sc1 * ((t3 * vb2y) - (vb1y * r12c1));
    dthetadr[0][2][2] = sc1 * ((t3 * vb2z) - (vb1z * r12c1));

    // angle23

    dthetadr[1][1][0] = sc2 * ((t2 * vb2x) + (vb3x * r12c2));
    dthetadr[1][1][1] = sc2 * ((t2 * vb2y) + (vb3y * r12c2));
    dthetadr[1][1][2] = sc2 * ((t2 * vb2z) + (vb3z * r12c2));

    dthetadr[1][2][0] = sc2 * ((-t2 * vb2x) - (vb3x * r12c2) +
                               (t4 * vb3x) + (vb2x * r12c2));
    dthetadr[1][2][1] = sc2 * ((-t2 * vb2y) - (vb3y * r12c2) +
                               (t4 * vb3y) + (vb2y * r12c2));
    dthetadr[1][2][2] = sc2 * ((-t2 * vb2z) - (vb3z * r12c2) +
                               (t4 * vb3z) + (vb2z * r12c2));

    dthetadr[1][3][0] = -sc2 * ((t4 * vb3x) + (vb2x * r12c2));
    dthetadr[1][3][1] = -sc2 * ((t4 * vb3y) + (vb2y * r12c2));
    dthetadr[1][3][2] = -sc2 * ((t4 * vb3z) + (vb2z * r12c2));

    // angle/angle/torsion cutoff

    da1 = acos(costh12) - aat_theta0_1[type] ;
    da2 = acos(costh23) - aat_theta0_1[type] ;
    double dtheta = aat_theta0_2[type]-aat_theta0_1[type];

    fphi = 0.0;
    fpphi = 0.0;
    if (phil < 0) phil +=MY_2PI;
    uf_lookup(type, phil, fphi, fpphi);

    double gt = aat_k[type];
    double gtt = aat_k[type];
    double gpt = 0;
    double gptt = 0;

    if ( acos(costh12) > aat_theta0_1[type]) {
      gt *= 1-da1*da1/dtheta/dtheta;
      gpt = -aat_k[type]*2*da1/dtheta/dtheta;
    }

    if ( acos(costh23) > aat_theta0_1[type]) {
      gtt *= 1-da2*da2/dtheta/dtheta;
      gptt = -aat_k[type]*2*da2/dtheta/dtheta;
    }

    if (eflag) edihedral = gt*gtt*fphi;

      for (i = 0; i < 4; i++)
        for (j = 0; j < 3; j++)
          fabcd[i][j] -=  - gt*gtt*fpphi*dphidr[i][j]
            - gt*gptt*fphi*dthetadr[1][i][j] + gpt*gtt*fphi*dthetadr[0][i][j];

    // apply force to each of 4 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += fabcd[0][0];
      f[i1][1] += fabcd[0][1];
      f[i1][2] += fabcd[0][2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += fabcd[1][0];
      f[i2][1] += fabcd[1][1];
      f[i2][2] += fabcd[1][2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += fabcd[2][0];
      f[i3][1] += fabcd[2][1];
      f[i3][2] += fabcd[2][2];
    }

    if (newton_bond || i4 < nlocal) {
      f[i4][0] += fabcd[3][0];
      f[i4][1] += fabcd[3][1];
      f[i4][2] += fabcd[3][2];
    }

    if (evflag)
      ev_tally(i1,i2,i3,i4,nlocal,newton_bond,edihedral,
               fabcd[0],fabcd[2],fabcd[3],
               vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralTableCut::allocate()
{
  allocated = 1;
  int n = atom->ndihedraltypes;

  memory->create(aat_k,n+1,"dihedral:aat_k");
  memory->create(aat_theta0_1,n+1,"dihedral:aat_theta0_1");
  memory->create(aat_theta0_2,n+1,"dihedral:aat_theta0_2");

  memory->create(tabindex,n+1,"dihedral:tabindex");
  memory->create(setflag,n+1,"dihedral:setflag");

  for (int i = 1; i <= n; i++)
    setflag[i] = 0;
}

void DihedralTableCut::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal dihedral_style command");

  if (strcmp(arg[0],"linear") == 0) tabstyle = LINEAR;
  else if (strcmp(arg[0],"spline") == 0) tabstyle = SPLINE;
  else error->all(FLERR,"Unknown table style in dihedral style table_cut");

  tablength = force->inumeric(FLERR,arg[1]);
  if (tablength < 3)
    error->all(FLERR,"Illegal number of dihedral table entries");
  // delete old tables, since cannot just change settings

  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(tabindex);
  }
  allocated = 0;

  ntables = 0;
  tables = NULL;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
   arg1 = "aat" -> AngleAngleTorsion coeffs
   arg1 -> Dihedral coeffs
------------------------------------------------------------------------- */

void DihedralTableCut::coeff(int narg, char **arg)
{

  if (narg != 7) error->all(FLERR,"Incorrect args for dihedral coefficients");
  if (!allocated) allocate();
  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->ndihedraltypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[2]);
  double theta0_1_one = force->numeric(FLERR,arg[3]);
  double theta0_2_one = force->numeric(FLERR,arg[4]);

  // convert theta0's from degrees to radians

  for (int i = ilo; i <= ihi; i++) {
    aat_k[i] = k_one;
    aat_theta0_1[i] = theta0_1_one/180.0 * MY_PI;
    aat_theta0_2[i] = theta0_2_one/180.0 * MY_PI;
  }

  int me;
  MPI_Comm_rank(world,&me);
  tables = (Table *)
    memory->srealloc(tables,(ntables+1)*sizeof(Table), "dihedral:tables");
  Table *tb = &tables[ntables];
  null_table(tb);
  if (me == 0) read_table(tb,arg[5],arg[6]);
  bcast_table(tb);

  // --- check the angle data for range errors ---
  // ---  and resolve issues with periodicity  ---

  if (tb->ninput < 2) {
    string err_msg;
    err_msg = string("Invalid dihedral table length (")
      + string(arg[5]) + string(").");
    error->one(FLERR,err_msg.c_str());
  }
  else if ((tb->ninput == 2) && (tabstyle == SPLINE)) {
    string err_msg;
    err_msg = string("Invalid dihedral spline table length. (Try linear)\n (")
      + string(arg[5]) + string(").");
    error->one(FLERR,err_msg.c_str());
  }

  // check for monotonicity
  for (int i=0; i < tb->ninput-1; i++) {
    if (tb->phifile[i] >= tb->phifile[i+1]) {
      stringstream i_str;
      i_str << i+1;
      string err_msg =
        string("Dihedral table values are not increasing (") +
        string(arg[5]) + string(", ")+i_str.str()+string("th entry)");
      if (i==0)
        err_msg += string("\n(This is probably a mistake with your table format.)\n");
      error->all(FLERR,err_msg.c_str());
    }
  }

  // check the range of angles
  double philo = tb->phifile[0];
  double phihi = tb->phifile[tb->ninput-1];
  if (tb->use_degrees) {
    if ((phihi - philo) >= 360) {
      string err_msg;
      err_msg = string("Dihedral table angle range must be < 360 degrees (")
        +string(arg[5]) + string(").");
      error->all(FLERR,err_msg.c_str());
    }
  }
  else {
    if ((phihi - philo) >= MY_2PI) {
      string err_msg;
      err_msg = string("Dihedral table angle range must be < 2*PI radians (")
        + string(arg[5]) + string(").");
      error->all(FLERR,err_msg.c_str());
    }
  }

  // convert phi from degrees to radians
  if (tb->use_degrees) {
    for (int i=0; i < tb->ninput; i++) {
      tb->phifile[i] *= MY_PI/180.0;
      // I assume that if angles are in degrees, then the forces (f=dU/dphi)
      // are specified with "phi" in degrees as well.
      tb->ffile[i] *= 180.0/MY_PI;
    }
  }

  // We want all the phi dihedral angles to lie in the range from 0 to 2*PI.
  // But I don't want to restrict users to input their data in this range.
  // We also want the angles to be sorted in increasing order.
  // This messy code fixes these problems with the user's data:
  {
    double *phifile_tmp = new double [tb->ninput];  //temporary arrays
    double *ffile_tmp = new double [tb->ninput];  //used for sorting
    double *efile_tmp = new double [tb->ninput];

    // After re-imaging, does the range of angles cross the 0 or 2*PI boundary?
    // If so, find the discontinuity:
    int i_discontinuity = tb->ninput;
    for (int i=0; i < tb->ninput; i++) {
      double phi   = tb->phifile[i];
      // Add a multiple of 2*PI to phi until it lies in the range [0, 2*PI).
      phi -= MY_2PI * floor(phi/MY_2PI);
      phifile_tmp[i] = phi;
      efile_tmp[i] = tb->efile[i];
      ffile_tmp[i] = tb->ffile[i];
      if ((i>0) && (phifile_tmp[i] < phifile_tmp[i-1])) {
        //There should only be at most one discontinuity, because we have
        //insured that the data was sorted before imaging, and because the
        //range of angle values does not exceed 2*PI.
        i_discontinuity = i;
      }
    }

    int I = 0;
    for (int i = i_discontinuity; i < tb->ninput; i++) {
      tb->phifile[I] = phifile_tmp[i];
      tb->efile[I] = efile_tmp[i];
      tb->ffile[I] = ffile_tmp[i];
      I++;
    }
    for (int i = 0; i < i_discontinuity; i++) {
      tb->phifile[I] = phifile_tmp[i];
      tb->efile[I] = efile_tmp[i];
      tb->ffile[I] = ffile_tmp[i];
      I++;
    }

    // clean up temporary storage
    delete[] phifile_tmp;
    delete[] ffile_tmp;
    delete[] efile_tmp;
  }

  // spline read-in and compute r,e,f vectors within table

  spline_table(tb);
  compute_table(tb);

  // Optional: allow the user to print out the interpolated spline tables

  if (me == 0) {
    if (checkU_fname && (strlen(checkU_fname) != 0))
    {
      ofstream checkU_file;
      checkU_file.open(checkU_fname, ios::out);
      for (int i=0; i < tablength; i++) {
        double phi = i*MY_2PI/tablength;
        double u = tb->e[i];
        if (tb->use_degrees)
          phi *= 180.0/MY_PI;
        checkU_file << phi << " " << u << "\n";
      }
      checkU_file.close();
    }
    if (checkF_fname && (strlen(checkF_fname) != 0))
    {
      ofstream checkF_file;
      checkF_file.open(checkF_fname, ios::out);
      for (int i=0; i < tablength; i++)
      {
        double phi = i*MY_2PI/tablength;
        double f;
        if ((tabstyle == SPLINE) && (tb->f_unspecified)) {
          double dU_dphi =
            // (If the user did not specify the forces now, AND the user
            //  selected the "spline" option, (as opposed to "linear")
            //  THEN the tb->f array is uninitialized, so there's
            //  no point to print out the contents of the tb->f[] array.
            //  Instead, later on, we will calculate the force using the
            //  -cyc_splintD() routine to calculate the derivative of the
            //  energy spline, using the energy data (tb->e[]).
            //  To be nice and report something, I do the same thing here.)
            cyc_splintD(tb->phi, tb->e, tb->e2, tablength, MY_2PI,phi);
          f = -dU_dphi;
        } else
          // Otherwise we calculated the tb->f[] array.  Report its contents.
          f = tb->f[i];
        if (tb->use_degrees) {
          phi *= 180.0/MY_PI;
          // If the user wants degree angle units, we should convert our
          // internal force tables (in energy/radians) to (energy/degrees)
          f *= MY_PI/180.0;
        }
        checkF_file << phi << " " << f << "\n";
      }
      checkF_file.close();
    } // if (checkF_fname && (strlen(checkF_fname) != 0))
  } // if (me == 0)

  // store ptr to table in tabindex
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    tabindex[i] = ntables;
    //phi0[i] = tb->phi0; <- equilibrium dihedral angles not supported
    setflag[i] = 1;
    count++;
  }
  ntables++;

  if (count == 0) error->all(FLERR,"Incorrect args for dihedral coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void DihedralTableCut::write_restart(FILE *fp)
{
  fwrite(&tabstyle,sizeof(int),1,fp);
  fwrite(&tablength,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void DihedralTableCut::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&tabstyle,sizeof(int),1,fp);
    fread(&tablength,sizeof(int),1,fp);
  }

  MPI_Bcast(&tabstyle,1,MPI_INT,0,world);
  MPI_Bcast(&tablength,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

void DihedralTableCut::null_table(Table *tb)
{
  tb->phifile = tb->efile = tb->ffile = NULL;
  tb->e2file = tb->f2file = NULL;
  tb->phi = tb->e = tb->de = NULL;
  tb->f = tb->df = tb->e2 = tb->f2 = NULL;
}

/* ---------------------------------------------------------------------- */

void DihedralTableCut::free_table(Table *tb)
{
  memory->destroy(tb->phifile);
  memory->destroy(tb->efile);
  memory->destroy(tb->ffile);
  memory->destroy(tb->e2file);
  memory->destroy(tb->f2file);

  memory->destroy(tb->phi);
  memory->destroy(tb->e);
  memory->destroy(tb->de);
  memory->destroy(tb->f);
  memory->destroy(tb->df);
  memory->destroy(tb->e2);
  memory->destroy(tb->f2);
}

/* ----------------------------------------------------------------------
   read table file, only called by proc 0
------------------------------------------------------------------------- */

static const int MAXLINE=2048;

void DihedralTableCut::read_table(Table *tb, char *file, char *keyword)
{
  char line[MAXLINE];

  // open file

  FILE *fp = force->open_potential(file);
  if (fp == NULL) {
    string err_msg = string("Cannot open file ") + string(file);
    error->one(FLERR,err_msg.c_str());
  }

  // loop until section found with matching keyword

  while (1) {
    if (fgets(line,MAXLINE,fp) == NULL) {
      string err_msg=string("Did not find keyword \"")
        +string(keyword)+string("\" in dihedral table file.");
      error->one(FLERR, err_msg.c_str());
    }
    if (strspn(line," \t\n\r") == strlen(line)) continue;  // blank line
    if (line[0] == '#') continue;                          // comment
    char *word = strtok(line," \t\n\r");
    if (strcmp(word,keyword) == 0) break;           // matching keyword
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);                         // no match, skip section
    param_extract(tb,line);
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
    for (int i = 0; i < tb->ninput; i++)
      utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
  }

  // read args on 2nd line of section
  // allocate table arrays for file values

  utils::sfgets(FLERR,line,MAXLINE,fp,file,error);
  param_extract(tb,line);
  memory->create(tb->phifile,tb->ninput,"dihedral:phifile");
  memory->create(tb->efile,tb->ninput,"dihedral:efile");
  memory->create(tb->ffile,tb->ninput,"dihedral:ffile");

  // read a,e,f table values from file

  int itmp;
  for (int i = 0; i < tb->ninput; i++) {
    // Read the next line.  Make sure the file is long enough.
    utils::sfgets(FLERR,line,MAXLINE,fp,file,error);

    // Skip blank lines and delete text following a '#' character
    char *pe = strchr(line, '#');
    if (pe != NULL) *pe = '\0'; //terminate string at '#' character
    char *pc = line;
    while ((*pc != '\0') && isspace(*pc))
      pc++;
    if (*pc != '\0') { //If line is not a blank line
      stringstream line_ss(line);
      if (tb->f_unspecified) {
        line_ss >> itmp;
        line_ss >> tb->phifile[i];
        line_ss >> tb->efile[i];
      } else {
        line_ss >> itmp;
        line_ss >> tb->phifile[i];
        line_ss >> tb->efile[i];
        line_ss >> tb->ffile[i];
      }
      if (! line_ss) {
        stringstream err_msg;
        err_msg << "Read error in table "<< keyword<<", near line "<<i+1<<"\n"
                << "   (Check to make sure the number of colums is correct.)";
        if ((! tb->f_unspecified) && (i==0))
          err_msg << "\n   (This sometimes occurs if users forget to specify the \"NOF\" option.)\n";
        error->one(FLERR, err_msg.str().c_str());
      }
    } else //if it is a blank line, then skip it.
      i--;
  } //for (int i = 0; (i < tb->ninput) && fp; i++) {

  fclose(fp);
}

/* ----------------------------------------------------------------------
   build spline representation of e,f over entire range of read-in table
   this function sets these values in e2file,f2file.
   I also perform a crude check for force & energy consistency.
------------------------------------------------------------------------- */

void DihedralTableCut::spline_table(Table *tb)
{
  memory->create(tb->e2file,tb->ninput,"dihedral:e2file");
  memory->create(tb->f2file,tb->ninput,"dihedral:f2file");

  if (cyc_spline(tb->phifile, tb->efile, tb->ninput,
                 MY_2PI,tb->e2file,comm->me == 0))
    error->one(FLERR,"Error computing dihedral spline tables");

  if (! tb->f_unspecified) {
    if (cyc_spline(tb->phifile, tb->ffile, tb->ninput,
                   MY_2PI, tb->f2file, comm->me == 0))
      error->one(FLERR,"Error computing dihedral spline tables");
  }

  // CHECK to help make sure the user calculated forces in a way
  // which is grossly numerically consistent with the energy table.
  if (! tb->f_unspecified) {
    int num_disagreements = 0;
    for (int i=0; i<tb->ninput; i++) {

      // Calculate what the force should be at the control points
      // by using linear interpolation of the derivatives of the energy:

      double phi_i = tb->phifile[i];

      // First deal with periodicity
      double phi_im1, phi_ip1;
      int im1 = i-1;
      if (im1 < 0) {
        im1 += tb->ninput;
        phi_im1 = tb->phifile[im1] - MY_2PI;
      }
      else
        phi_im1 = tb->phifile[im1];
      int ip1 = i+1;
      if (ip1 >= tb->ninput) {
        ip1 -= tb->ninput;
        phi_ip1 = tb->phifile[ip1] + MY_2PI;
      }
      else
        phi_ip1 = tb->phifile[ip1];

      // Now calculate the midpoints above and below phi_i = tb->phifile[i]
      double phi_lo= 0.5*(phi_im1 + phi_i); //midpoint between phi_im1 and phi_i
      double phi_hi= 0.5*(phi_i + phi_ip1); //midpoint between phi_i and phi_ip1

      // Use a linear approximation to the derivative at these two midpoints
      double dU_dphi_lo =
        (tb->efile[i] - tb->efile[im1])
        /
        (phi_i - phi_im1);
      double dU_dphi_hi =
        (tb->efile[ip1] - tb->efile[i])
        /
        (phi_ip1 - phi_i);

      // Now calculate the derivative at position
      // phi_i (=tb->phifile[i]) using linear interpolation

      double a = (phi_i - phi_lo) / (phi_hi - phi_lo);
      double b = (phi_hi - phi_i) / (phi_hi - phi_lo);
      double dU_dphi = a*dU_dphi_lo  +  b*dU_dphi_hi;
      double f = -dU_dphi;
      // Alternately, we could use spline interpolation instead:
      // double f = - splintD(tb->phifile, tb->efile, tb->e2file,
      //                      tb->ninput, MY_2PI, tb->phifile[i]);
      // This is the way I originally did it, but I trust
      // the ugly simple linear way above better.
      // Recall this entire block of code doess not calculate
      // anything important.  It does not have to be perfect.
      // We are only checking for stupid user errors here.

      if ((f != 0.0) &&
          (tb->ffile[i] != 0.0) &&
          ((f/tb->ffile[i] < 0.5) || (f/tb->ffile[i] > 2.0))) {
        num_disagreements++;
      }
    } // for (int i=0; i<tb->ninput; i++)

    if ((num_disagreements > tb->ninput/2) && (num_disagreements > 2)) {
      string msg("Dihedral table has inconsistent forces and energies. (Try \"NOF\".)\n");
      error->all(FLERR,msg.c_str());
    }

  } // check for consistency if (! tb->f_unspecified)

} // DihedralTable::spline_table()


/* ----------------------------------------------------------------------
   compute a,e,f vectors from splined values
------------------------------------------------------------------------- */

void DihedralTableCut::compute_table(Table *tb)
{
  //delta = table spacing in dihedral angle for tablength (cyclic) bins
  tb->delta = MY_2PI / tablength;
  tb->invdelta = 1.0/tb->delta;
  tb->deltasq6 = tb->delta*tb->delta / 6.0;

  // N evenly spaced bins in dihedral angle from 0 to 2*PI
  // phi,e,f = value at lower edge of bin
  // de,df values = delta values of e,f (cyclic, in this case)
  // phi,e,f,de,df are arrays containing "tablength" number of entries

  memory->create(tb->phi,tablength,"dihedral:phi");
  memory->create(tb->e,tablength,"dihedral:e");
  memory->create(tb->de,tablength,"dihedral:de");
  memory->create(tb->f,tablength,"dihedral:f");
  memory->create(tb->df,tablength,"dihedral:df");
  memory->create(tb->e2,tablength,"dihedral:e2");
  memory->create(tb->f2,tablength,"dihedral:f2");

  if (tabstyle == SPLINE) {
    // Use cubic spline interpolation to calculate the entries in the
    // internal table. (This is true regardless...even if tabstyle!=SPLINE.)
    for (int i = 0; i < tablength; i++) {
      double phi = i*tb->delta;
      tb->phi[i] = phi;
      tb->e[i]= cyc_splint(tb->phifile,tb->efile,tb->e2file,tb->ninput,MY_2PI,phi);
      if (! tb->f_unspecified)
        tb->f[i] = cyc_splint(tb->phifile,tb->ffile,tb->f2file,tb->ninput,MY_2PI,phi);
    }
  } // if (tabstyle == SPLINE)
  else if (tabstyle == LINEAR) {
    if (! tb->f_unspecified) {
      for (int i = 0; i < tablength; i++) {
        double phi = i*tb->delta;
        tb->phi[i] = phi;
        tb->e[i]= cyc_lin(tb->phifile,tb->efile,tb->ninput,MY_2PI,phi);
        tb->f[i]= cyc_lin(tb->phifile,tb->ffile,tb->ninput,MY_2PI,phi);
      }
    }
    else {
      for (int i = 0; i < tablength; i++) {
        double phi = i*tb->delta;
        tb->phi[i] = phi;
        tb->e[i]= cyc_lin(tb->phifile,tb->efile,tb->ninput,MY_2PI,phi);
      }
      // In the linear case, if the user did not specify the forces, then we
      // must generate the "f" array. Do this using linear interpolation
      // of the e array (which itself was generated above)
      for (int i = 0; i < tablength; i++) {
        int im1 = i-1; if (im1 < 0) im1 += tablength;
        int ip1 = i+1; if (ip1 >= tablength) ip1 -= tablength;
        double dedx = (tb->e[ip1] - tb->e[im1]) / (2.0 * tb->delta);
        // (This is the average of the linear slopes on either side of the node.
        //  Note that the nodes in the internal table are evenly spaced.)
        tb->f[i] = -dedx;
      }
    }


    // Fill the linear interpolation tables (de, df)
    for (int i = 0; i < tablength; i++) {
      int ip1 = i+1; if (ip1 >= tablength) ip1 -= tablength;
      tb->de[i] = tb->e[ip1] - tb->e[i];
      tb->df[i] = tb->f[ip1] - tb->f[i];
    }
  } // else if (tabstyle == LINEAR)



  cyc_spline(tb->phi, tb->e, tablength, MY_2PI, tb->e2, comm->me == 0);
  if (! tb->f_unspecified)
    cyc_spline(tb->phi, tb->f, tablength, MY_2PI, tb->f2, comm->me == 0);
}


/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: N value NOF DEGREES RADIANS
   N is required, other params are optional
------------------------------------------------------------------------- */

void DihedralTableCut::param_extract(Table *tb, char *line)
{
  //tb->theta0 = 180.0; <- equilibrium angles not supported
  tb->ninput = 0;
  tb->f_unspecified = false; //default
  tb->use_degrees   = true;  //default

  char *word = strtok(line," \t\n\r\f");
  while (word) {
    if (strcmp(word,"N") == 0) {
      word = strtok(NULL," \t\n\r\f");
      tb->ninput = atoi(word);
    }
    else if (strcmp(word,"NOF") == 0) {
      tb->f_unspecified = true;
    }
    else if ((strcmp(word,"DEGREES") == 0) || (strcmp(word,"degrees") == 0)) {
      tb->use_degrees = true;
    }
    else if ((strcmp(word,"RADIANS") == 0) || (strcmp(word,"radians") == 0)) {
      tb->use_degrees = false;
    }
    else if (strcmp(word,"CHECKU") == 0) {
      word = strtok(NULL," \t\n\r\f");
      memory->sfree(checkU_fname);
      memory->create(checkU_fname,strlen(word)+1,"dihedral_table:checkU");
      strcpy(checkU_fname, word);
    }
    else if (strcmp(word,"CHECKF") == 0) {
      word = strtok(NULL," \t\n\r\f");
      memory->sfree(checkF_fname);
      memory->create(checkF_fname,strlen(word)+1,"dihedral_table:checkF");
      strcpy(checkF_fname, word);
    }
    // COMMENTING OUT:  equilibrium angles are not supported
    //else if (strcmp(word,"EQ") == 0) {
    //  word = strtok(NULL," \t\n\r\f");
    //  tb->theta0 = atof(word);
    //}
    else {
      string err_msg("Invalid keyword in dihedral angle table parameters");
      err_msg += string(" (") + string(word) + string(")");
      error->one(FLERR,err_msg.c_str());
    }
    word = strtok(NULL," \t\n\r\f");
  }

  if (tb->ninput == 0)
    error->one(FLERR,"Dihedral table parameters did not set N");

} // DihedralTable::param_extract()


/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,phifile,efile,ffile,
     f_unspecified,use_degrees
------------------------------------------------------------------------- */

void DihedralTableCut::bcast_table(Table *tb)
{
  MPI_Bcast(&tb->ninput,1,MPI_INT,0,world);

  int me;
  MPI_Comm_rank(world,&me);
  if (me > 0) {
    memory->create(tb->phifile,tb->ninput,"dihedral:phifile");
    memory->create(tb->efile,tb->ninput,"dihedral:efile");
    memory->create(tb->ffile,tb->ninput,"dihedral:ffile");
  }

  MPI_Bcast(tb->phifile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->efile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->ffile,tb->ninput,MPI_DOUBLE,0,world);

  MPI_Bcast(&tb->f_unspecified,1,MPI_INT,0,world);
  MPI_Bcast(&tb->use_degrees,1,MPI_INT,0,world);

  // COMMENTING OUT:  equilibrium angles are not supported
  //MPI_Bcast(&tb->theta0,1,MPI_DOUBLE,0,world);
}


