/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef CUDA_MODIFY_FLAGS_H
#define CUDA_MODIFY_FLAGS_H

#include "fix.h"

namespace LAMMPS_NS {
namespace FixConstCuda {
  static const int INITIAL_INTEGRATE_CUDA = FixConst::FIX_CONST_LAST << 0;
  static const int POST_INTEGRATE_CUDA    = FixConst::FIX_CONST_LAST << 1;
  static const int PRE_EXCHANGE_CUDA      = FixConst::FIX_CONST_LAST << 2;
  static const int PRE_NEIGHBOR_CUDA      = FixConst::FIX_CONST_LAST << 3;
  static const int PRE_FORCE_CUDA         = FixConst::FIX_CONST_LAST << 4;
  static const int POST_FORCE_CUDA        = FixConst::FIX_CONST_LAST << 5;
  static const int FINAL_INTEGRATE_CUDA   = FixConst::FIX_CONST_LAST << 6;
  static const int END_OF_STEP_CUDA       = FixConst::FIX_CONST_LAST << 7;
  static const int THERMO_ENERGY_CUDA     = FixConst::FIX_CONST_LAST << 8;
  static const int MIN_POST_FORCE_CUDA    = FixConst::FIX_CONST_LAST << 9;
}
}
// remember not to shift over 31 bits

#endif // CUDA_MODIFY_FLAGS_H
