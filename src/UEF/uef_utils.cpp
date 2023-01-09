// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: David Nicholson (MIT)
-------------------------------------------------------------------------

   This class contains functions to calculate the evolution of the periodic
   simulation box under elongational flow as described by Matthew Dobson
   in the arXiv preprint at https://arxiv.org/abs/1408.7078

   Additionally, there are methods to do a lattice reduction to further
   reduce the simulation box using the method of Igor Semaev at
   https://link.springer.com/chapter/10.1007%2F3-540-44670-2_13
*/

#include <cmath>
#include "uef_utils.h"

namespace LAMMPS_NS {
  namespace UEF_utils{

UEFBox::UEFBox()
{

  // initial box (also an inverse eigenvector matrix of automorphisms)

  double x = 0.327985277605681;
  double y = 0.591009048506103;
  double z = 0.736976229099578;
  l0[0][0]= z; l0[0][1]= y; l0[0][2]= x;
  l0[1][0]=-x; l0[1][1]= z; l0[1][2]=-y;
  l0[2][0]=-y; l0[2][1]= x; l0[2][2]= z;

  // spectra of the two automorpisms (log of eigenvalues)

  w1[0]=-1.177725211523360;
  w1[1]=-0.441448620566067;
  w1[2]= 1.619173832089425;
  w2[0]= w1[1];
  w2[1]= w1[2];
  w2[2]= w1[0];

  // initialize theta
  // strain = w1 * theta1 + w2 * theta2

  theta[0]=theta[1]=0;

  //set up the initial box l and change of basis matrix r

  for (int k=0;k<3;k++)
    for (int j=0;j<3;j++) {
      l[k][j] = l0[k][j];
      r[j][k]=(j==k);
      ri[j][k]=(j==k);
    }

  // get the initial rotation and upper triangular matrix

  rotation_matrix(rot, lrot ,l);

  // this is just a way to calculate the automorphisms
  // themselves, which play a minor role in the calculations
  // it's overkill, but only called once

  double t1[3][3];
  double t1i[3][3];
  double t2[3][3];
  double t2i[3][3];
  double l0t[3][3];
  for (int k=0; k<3; ++k)
    for (int j=0; j<3; ++j) {
      t1[k][j] = exp(w1[k])*l0[k][j];
      t1i[k][j] = exp(-w1[k])*l0[k][j];
      t2[k][j] = exp(w2[k])*l0[k][j];
      t2i[k][j] = exp(-w2[k])*l0[k][j];
      l0t[k][j] = l0[j][k];
    }
  mul_m2(l0t,t1);
  mul_m2(l0t,t1i);
  mul_m2(l0t,t2);
  mul_m2(l0t,t2i);
  for (int k=0; k<3; ++k)
    for (int j=0; j<3; ++j) {
      a1[k][j] = round(t1[k][j]);
      a1i[k][j] = round(t1i[k][j]);
      a2[k][j] = round(t2[k][j]);
      a2i[k][j] = round(t2i[k][j]);
    }

  // winv used to transform between
  // strain increments and theta increments

  winv[0][0] = w2[1];
  winv[0][1] = -w2[0];
  winv[1][0] = -w1[1];
  winv[1][1] = w1[0];
  double d = w1[0]*w2[1] - w2[0]*w1[1];
  for (int k=0;k<2;k++)
    for (int j=0;j<2;j++)
      winv[k][j] /= d;
}

/* ----------------------------------------------------------------------
   get volume-correct r basis in: basis*cbrt(vol) = q*r
------------------------------------------------------------------------- */
void UEFBox::get_box(double x[3][3], double v)
{
  v = cbrtf(v);
  for (int k=0;k<3;k++)
    for (int j=0;j<3;j++)
      x[k][j] = lrot[k][j]*v;
}

/* ----------------------------------------------------------------------
   get rotation matrix q in: basis = q*r
------------------------------------------------------------------------- */
void UEFBox::get_rot(double x[3][3])
{
  for (int k=0;k<3;k++)
    for (int j=0;j<3;j++)
      x[k][j]=rot[k][j];
}

/* ----------------------------------------------------------------------
   get inverse change of basis matrix
------------------------------------------------------------------------- */
void UEFBox::get_inverse_cob(int x[3][3])
{
  for (int k=0;k<3;k++)
    for (int j=0;j<3;j++)
      x[k][j]=ri[k][j];
}

/* ----------------------------------------------------------------------
   apply diagonal, incompressible deformation
------------------------------------------------------------------------- */
void UEFBox::step_deform(const double ex, const double ey)
{
  // increment theta values used in the reduction

  theta[0] +=winv[0][0]*ex + winv[0][1]*ey;
  theta[1] +=winv[1][0]*ex + winv[1][1]*ey;

  // deformation of the box. reduce() needs to be called regularly or
  // calculation will become unstable

  double eps[3];
  eps[0]=ex; eps[1] = ey; eps[2] = -ex-ey;
  for (int k=0;k<3;k++) {
    eps[k] = exp(eps[k]);
    l[k][0] = eps[k]*l[k][0];
    l[k][1] = eps[k]*l[k][1];
    l[k][2] = eps[k]*l[k][2];
  }
  rotation_matrix(rot,lrot, l);
}

/* ----------------------------------------------------------------------
   reduce the current basis
------------------------------------------------------------------------- */
bool UEFBox::reduce()
{
  // determine how many times to apply the automorphisms and find new theta
  // values

  int f1 = round(theta[0]);
  int f2 = round(theta[1]);
  theta[0] -= f1;
  theta[1] -= f2;

  // store old change or basis matrix to determine if it changes

  int r0[3][3];
  for (int k=0;k<3;k++)
    for (int j=0;j<3;j++)
      r0[k][j]=r[k][j];

  // this modifies the old change basis matrix to handle the case where the
  // automorphism transforms the box but the reduced basis doesn't change
  // (r0 should still equal r at the end)

  if (f1 > 0) for (int k=0;k<f1;k++) mul_m2 (a1,r0);
  if (f1 < 0) for (int k=0;k<-f1;k++) mul_m2 (a1i,r0);
  if (f2 > 0) for (int k=0;k<f2;k++) mul_m2 (a2,r0);
  if (f2 < 0) for (int k=0;k<-f2;k++) mul_m2 (a2i,r0);

  // robust reduction to the box defined by Dobson

  for (int k=0;k<3;k++) {
    double eps = exp(theta[0]*w1[k]+theta[1]*w2[k]);
    l[k][0] = eps*l0[k][0];
    l[k][1] = eps*l0[k][1];
    l[k][2] = eps*l0[k][2];
  }

  // further reduce the box using greedy reduction and check
  // if it changed from the last step using the change of basis
  // matrices r and r0

  greedy(l,r,ri);

  // multiplying the inverse by the old change of basis matrix gives
  // the inverse of the transformation itself (should be identity if
  // no reduction takes place). This is used for image flags only.

  mul_m1(ri,r0);
  rotation_matrix(rot,lrot, l);
  return !mat_same(r,r0);
}

/* ----------------------------------------------------------------------
   set the strain to a specific value
------------------------------------------------------------------------- */
void UEFBox::set_strain(const double ex, const double ey)
{
  theta[0]  = winv[0][0]*ex + winv[0][1]*ey;
  theta[1]  = winv[1][0]*ex + winv[1][1]*ey;
  theta[0] -= round(theta[0]);
  theta[1] -= round(theta[1]);

  for (int k=0;k<3;k++) {
    double eps = exp(theta[0]*w1[k]+theta[1]*w2[k]);
    l[k][0] = eps*l0[k][0];
    l[k][1] = eps*l0[k][1];
    l[k][2] = eps*l0[k][2];
  }
  greedy(l,r,ri);
  rotation_matrix(rot,lrot, l);
}

/* ----------------------------------------------------------------------
   qr reduction using householder reflections
   q*m = r. q is orthogonal. m is input matrix. r is upper triangular
------------------------------------------------------------------------- */
void rotation_matrix(double q[3][3], double r[3][3], const double m[3][3])
{
  for (int k=0;k<3;k++)
    for (int j=0;j<3;j++)
      r[k][j] = m[k][j];

  double a = -sqrt(col_prod(r,0,0))*r[0][0]/fabs(r[0][0]);
  double v[3];
  v[0] = r[0][0]-a;
  v[1] = r[1][0];
  v[2] = r[2][0];
  a = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] /= a; v[1] /= a; v[2] /= a;
  double qt[3][3];
  for (int k=0;k<3;k++)
    for (int j=0;j<3;j++) {
      qt[k][j] = (k==j) - 2*v[k]*v[j];
      q[k][j]= qt[k][j];
    }
  mul_m2(qt,r);
  a = -sqrt(r[1][1]*r[1][1] + r[2][1]*r[2][1])*r[1][1]/fabs(r[1][1]);
  v[0] = 0;
  v[1] = r[1][1] - a;
  v[2] = r[2][1];
  a = sqrt(v[1]*v[1]+v[2]*v[2]);
  v[1] /= a;
  v[2] /= a;
  for (int k=0;k<3;k++)
    for (int j=0;j<3;j++)
      qt[k][j] = (k==j) - 2*v[k]*v[j];
  mul_m2(qt,r);
  mul_m2(qt,q);

  // this makes r have positive diagonals
  // q*m = r <==> (-q)*m = (-r) will hold row-wise

  if (r[0][0] < 0) { neg_row(q,0); neg_row(r,0); }
  if (r[1][1] < 0) { neg_row(q,1); neg_row(r,1); }
  if (r[2][2] < 0) { neg_row(q,2); neg_row(r,2); }
}

/* ----------------------------------------------------------------------
   sort columns of b in order of increasing length
   mimic column operations on ri and r
------------------------------------------------------------------------- */
void col_sort(double b[3][3],int r[3][3],int ri[3][3])
{
  if (col_prod(b,0,0)>col_prod(b,1,1)) {
    col_swap(b,0,1);
    col_swap(r,0,1);
    col_swap(ri,0,1);
  }
  if (col_prod(b,0,0)>col_prod(b,2,2)) {
    col_swap(b,0,2);
    col_swap(r,0,2);
    col_swap(ri,0,2);
  }
  if (col_prod(b,1,1)>col_prod(b,2,2)) {
    col_swap(b,1,2);
    col_swap(r,1,2);
    col_swap(ri,1,2);
  }
}

/* ----------------------------------------------------------------------
   1-2 reduction (Graham-Schmidt)
------------------------------------------------------------------------- */
void red12(double b[3][3],int r[3][3],int ri[3][3])
{
  int y = round(col_prod(b,0,1)/col_prod(b,0,0));
  b[0][1] -= y*b[0][0];
  b[1][1] -= y*b[1][0];
  b[2][1] -= y*b[2][0];

  r[0][1] -= y*r[0][0];
  r[1][1] -= y*r[1][0];
  r[2][1] -= y*r[2][0];

  ri[0][0] += y*ri[0][1];
  ri[1][0] += y*ri[1][1];
  ri[2][0] += y*ri[2][1];

  if (col_prod(b,1,1) < col_prod(b,0,0)) {
    col_swap(b,0,1);
    col_swap(r,0,1);
    col_swap(ri,0,1);
    red12(b,r,ri);
  }
}

/* ----------------------------------------------------------------------
   Apply the Semaev condition for a 3-reduced basis
------------------------------------------------------------------------- */
void red3(double b[3][3], int r[3][3], int ri[3][3])
{
  double b11 = col_prod(b,0,0);
  double b22 = col_prod(b,1,1);
  double b12 = col_prod(b,0,1);
  double b13 = col_prod(b,0,2);
  double b23 = col_prod(b,1,2);

  double y2 =-(b23/b22-b12/b22*b13/b11)/(1-b12/b11*b12/b22);
  double y1 =-(b13/b11-b12/b11*b23/b22)/(1-b12/b11*b12/b22);

  int x1=0;
  int x2=0;
  double min = col_prod(b,2,2);
  int x1v[2];
  int x2v[2];
  x1v[0] = floor(y1); x1v[1] = x1v[0]+1;
  x2v[0] = floor(y2); x2v[1] = x2v[0]+1;
  for (int k=0;k<2;k++)
    for (int j=0;j<2;j++) {
      double a[3];
      a[0] = b[0][2] + x1v[k]*b[0][0] + x2v[j]*b[0][1];
      a[1] = b[1][2] + x1v[k]*b[1][0] + x2v[j]*b[1][1];
      a[2] = b[2][2] + x1v[k]*b[2][0] + x2v[j]*b[2][1];
      double val=a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
      if (val<min) {
        min = val;
        x1 = x1v[k];
        x2 = x2v[j];
      }
    }
  if (x1 || x2) {
    b[0][2] += x1*b[0][0] + x2*b[0][1];
    b[1][2] += x1*b[1][0] + x2*b[1][1];
    b[2][2] += x1*b[2][0] + x2*b[2][1];
    r[0][2] += x1*r[0][0] + x2*r[0][1];
    r[1][2] += x1*r[1][0] + x2*r[1][1];
    r[2][2] += x1*r[2][0] + x2*r[2][1];
    ri[0][0] += -x1*ri[0][2];
    ri[1][0] += -x1*ri[1][2];
    ri[2][0] += -x1*ri[2][2];
    ri[0][1] += -x2*ri[0][2];
    ri[1][1] += -x2*ri[1][2];
    ri[2][1] += -x2*ri[2][2];
    greedy_recurse(b,r,ri); // note the recursion step is here
  }
}

/* ----------------------------------------------------------------------
   the meat of the greedy reduction algorithm
------------------------------------------------------------------------- */
void greedy_recurse(double b[3][3], int r[3][3], int ri[3][3])
{
  col_sort(b,r,ri);
  red12(b,r,ri);
  red3(b,r,ri); // recursive caller
}

/* ----------------------------------------------------------------------
   reduce the basis b. also output the change of basis matrix r and its
   inverse ri
------------------------------------------------------------------------- */
void greedy(double b[3][3],int r[3][3],int ri[3][3])
{
  r[0][1]=r[0][2]=r[1][0]=r[1][2]=r[2][0]=r[2][1]=0;
  r[0][0]=r[1][1]=r[2][2]=1;
  ri[0][1]=ri[0][2]=ri[1][0]=ri[1][2]=ri[2][0]=ri[2][1]=0;
  ri[0][0]=ri[1][1]=ri[2][2]=1;
  greedy_recurse(b,r,ri);
  make_unique(b,r,ri);
  transpose(ri);
}

/* ----------------------------------------------------------------------
   A reduced basis isn't unique. This procedure will make it
   "more" unique. Degenerate cases are possible, but unlikely
   with floating point math.
------------------------------------------------------------------------- */
void make_unique(double b[3][3], int r[3][3], int ri[3][3])
{
  if (fabs(b[0][0]) < fabs(b[0][1])) {
    col_swap(b,0,1);
    col_swap(r,0,1);
    col_swap(ri,0,1);
  }
  if (fabs(b[0][0]) < fabs(b[0][2])) {
    col_swap(b,0,2);
    col_swap(r,0,2);
    col_swap(ri,0,2);
  }
  if (fabs(b[1][1]) < fabs(b[1][2])) {
    col_swap(b,1,2);
    col_swap(r,1,2);
    col_swap(ri,1,2);
  }

  if (b[0][0] < 0) {
    neg_col(b,0);
    neg_col(r,0);
    neg_col(ri,0);
  }
  if (b[1][1] < 0) {
    neg_col(b,1);
    neg_col(r,1);
    neg_col(ri,1);
  }
  if (det(b) < 0) {
    neg_col(b,2);
    neg_col(r,2);
    neg_col(ri,2);
  }
}
}}
