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

#ifndef LMP_SUFFIX_H
#define LMP_SUFFIX_H

namespace LAMMPS_NS {

namespace Suffix {
  static const int NONE = 0;
  static const int OPT  = 1<<0;
  static const int GPU  = 1<<1;
  static const int CUDA = 1<<2;
  static const int OMP  = 1<<3;
  static const int INTEL  = 1<<4;
}

}

#endif
