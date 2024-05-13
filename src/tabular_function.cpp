/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "tabular_function.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

TabularFunction::TabularFunction() :
    size(0), xmin(0.0), xmax(0.0), xmaxsq(0.0), rdx(0.0), vmax(0.0), xs(nullptr), ys(nullptr),
    ys1(nullptr), ys2(nullptr), ys3(nullptr), ys4(nullptr), ys5(nullptr), ys6(nullptr)
{
}

TabularFunction::~TabularFunction()
{
  delete[] xs;
  delete[] ys;
  delete[] ys1;
  delete[] ys2;
  delete[] ys3;
  delete[] ys4;
  delete[] ys5;
  delete[] ys6;
}

void TabularFunction::set_values(int n, double x1, double x2, double *values)
{
  reset_size(n);
  set_xrange(x1, x2);
  memcpy(ys, values, n * sizeof(double));
  initialize();
}

void TabularFunction::set_xrange(double x1, double x2)
{
  xmin = x1;
  xmax = x2;
  xmaxsq = x2 * x2;
}

void TabularFunction::reset_size(int n)
{
  if (n != size) {
    size = n;
    delete[] xs;
    delete[] ys;
    delete[] ys1;
    delete[] ys2;
    delete[] ys3;
    delete[] ys4;
    delete[] ys5;
    delete[] ys6;
    xs = new double[n];
    ys = new double[n];
    ys1 = new double[n];
    ys2 = new double[n];
    ys3 = new double[n];
    ys4 = new double[n];
    ys5 = new double[n];
    ys6 = new double[n];
  }
}

void TabularFunction::initialize()
{
  int i;
  rdx = (xmax - xmin) / (size - 1.0);
  vmax = 0.0;
  for (i = 0; i < size; i++)
    if (fabs(ys[i]) > vmax) vmax = fabs(ys[i]);
  for (i = 0; i < size; i++) xs[i] = xmin + i * rdx;
  rdx = 1.0 / rdx;
  ys1[0] = ys[1] - ys[0];
  ys1[1] = 0.5 * (ys[2] - ys[0]);
  ys1[size - 2] = 0.5 * (ys[size - 1] - ys[size - 3]);
  ys1[size - 1] = ys[size - 1] - ys[size - 2];
  for (i = 2; i < size - 2; i++)
    ys1[i] = ((ys[i - 2] - ys[i + 2]) + 8.0 * (ys[i + 1] - ys[i - 1])) / 12.0;
  for (i = 0; i < size - 1; i++) {
    ys2[i] = 3.0 * (ys[i + 1] - ys[i]) - 2.0 * ys1[i] - ys1[i + 1];
    ys3[i] = ys1[i] + ys1[i + 1] - 2.0 * (ys[i + 1] - ys[i]);
  }
  ys2[size - 1] = 0.0;
  ys3[size - 1] = 0.0;
  for (i = 0; i < size; i++) {
    ys4[i] = ys1[i] * rdx;
    ys5[i] = 2.0 * ys2[i] * rdx;
    ys6[i] = 3.0 * ys3[i] * rdx;
  }
}
