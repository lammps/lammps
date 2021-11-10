/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <cmath>
#include <direct.h>

// some symbols have different names in Windows

#undef ATOBIGINT
#define ATOBIGINT _atoi64

#define pclose _pclose
#define strdup _strdup

// the following functions are defined to get rid of
// 'ambiguous call to overloaded function' error in VSS for mismatched type arguments
#if !defined(__MINGW32__)
inline double pow(int i, int j)
{
  return pow((double) i, j);
}
inline double fabs(int i)
{
  return fabs((double) i);
}
inline double sqrt(int i)
{
  return sqrt((double) i);
}
#endif

inline double trunc(double x)
{
  return x > 0 ? floor(x) : ceil(x);
}

// Windows version of mkdir function does not have permission flags
#ifndef S_IRWXU
#define S_IRWXU 0
#endif
#ifndef S_IRGRP
#define S_IRGRP 0
#endif
#ifndef S_IXGRP
#define S_IXGRP 0
#endif
inline int mkdir(const char *path, int)
{
  return _mkdir(path);
}
