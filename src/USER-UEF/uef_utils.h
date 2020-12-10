/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: David Nicholson (MIT)
------------------------------------------------------------------------- */

#ifndef LMP_UEF_UTILS_H
#define LMP_UEF_UTILS_H

namespace LAMMPS_NS{ namespace UEF_utils {

class UEFBox
{
  public:
    UEFBox();
    void set_strain(const double, const double);
    void step_deform(const double,const double);
    bool reduce();
    void get_box(double[3][3], double);
    void get_rot(double[3][3]);
    void get_inverse_cob(int[3][3]);
  private:
    double l0[3][3]; // initial basis
    double w1[3],w2[3],winv[3][3];//omega1 and omega2 (spectra of automorphisms)
    double theta[2];
    double l[3][3], rot[3][3], lrot[3][3];
    int r[3][3],ri[3][3],a1[3][3],a2[3][3],a1i[3][3],a2i[3][3];
};

// lattice reduction routines

void greedy(double[3][3],int[3][3],int[3][3]);
void col_sort(double[3][3],int[3][3],int[3][3]);
void red12(double[3][3],int[3][3],int[3][3]);
void greedy_recurse(double[3][3],int[3][3],int[3][3]);
void red3(double [3][3],int r[3][3],int[3][3]);
void make_unique(double[3][3],int[3][3],int[3][3]);
void rotation_matrix(double[3][3],double[3][3],const double [3][3]);

// A few utility functions for 3x3 arrays

template<typename T>
T col_prod(T x[3][3], int c1, int c2)
{
  return x[0][c1]*x[0][c2]+x[1][c1]*x[1][c2]+x[2][c1]*x[2][c2];
}

template<typename T>
void col_swap(T x[3][3], int c1, int c2)
{
  for (int k=0;k<3;k++) {
    T t = x[k][c2];
    x[k][c2]=x[k][c1];
    x[k][c1]=t;
  }
}

template<typename T>
void neg_col(T x[3][3], int c1)
{
  x[0][c1] = -x[0][c1];
  x[1][c1] = -x[1][c1];
  x[2][c1] = -x[2][c1];
}

template<typename T>
void neg_row(T x[3][3], int c1)
{
  x[c1][0] = -x[c1][0];
  x[c1][1] = -x[c1][1];
  x[c1][2] = -x[c1][2];
}

template<typename T>
T det(T x[3][3])
{
  double val;
  val  = x[0][0]*(x[1][1]*x[2][2] - x[1][2]*x[2][1]);
  val -= x[0][1]*(x[1][0]*x[2][2] - x[1][2]*x[2][0]);
  val += x[0][2]*(x[1][0]*x[2][1] - x[1][1]*x[2][0]);
  return val;
}

template<typename T>
bool mat_same(T x1[3][3], T x2[3][3])
{
  bool v = true;
  for (int k=0;k<3;k++)
    for (int j=0;j<3;j++)
      v &= (x1[k][j]==x2[k][j]);
  return v;
}

template<typename T>
void transpose(T m[3][3])
{
  for (int k=0;k<3;k++)
    for (int j=k+1;j<3;j++) {
      T x = m[k][j];
      m[k][j] = m[j][k];
      m[j][k] = x;
    }
}

template<typename T1,typename T2>
void mul_m1(T1 m1[3][3], const T2 m2[3][3])
{
  T1 t[3][3];
  for (int k=0;k<3;k++)
    for (int j=0;j<3;j++)
      t[k][j]=m1[k][j];

  for (int k=0;k<3;k++)
    for (int j=0;j<3;j++)
      m1[k][j] = t[k][0]*m2[0][j] + t[k][1]*m2[1][j] + t[k][2]*m2[2][j];
}

template<typename T1, typename T2>
void mul_m2(const T1 m1[3][3], T2 m2[3][3])
{
  T2 t[3][3];
  for (int k=0;k<3;k++)
    for (int j=0;j<3;j++)
      t[k][j]=m2[k][j];

  for (int k=0;k<3;k++)
    for (int j=0;j<3;j++)
      m2[k][j] = m1[k][0]*t[0][j] + m1[k][1]*t[1][j] + m1[k][2]*t[2][j];
}

}
}
#endif
