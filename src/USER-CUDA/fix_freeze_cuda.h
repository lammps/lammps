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

#ifdef FIX_CLASS

FixStyle(freeze/cuda,FixFreezeCuda)

#else

#ifndef LMP_FIX_FREEZE_CUDA_H
#define LMP_FIX_FREEZE_CUDA_H

#include "fix.h"
#include "cuda_data.h"

namespace LAMMPS_NS {

class FixFreezeCuda : public Fix {
 public:
  FixFreezeCuda(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  double compute_vector(int);

 private:
  class Cuda *cuda;
  double foriginal[3],foriginal_all[3];
  cCudaData<double     , F_CFLOAT                   , x>* cu_foriginal;
  int force_flag;
};

}

#endif
#endif
