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

FixStyle(addforce/cuda,FixAddForceCuda)

#else

#ifndef LMP_FIX_ADD_FORCE_CUDA_H
#define LMP_FIX_ADD_FORCE_CUDA_H

#include "fix.h"
#include "cuda_data.h"

namespace LAMMPS_NS {

class FixAddForceCuda : public Fix {
 public:
  FixAddForceCuda(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

 private:
  class Cuda *cuda;
  int iregion;
  double xvalue,yvalue,zvalue;
  double foriginal[4],foriginal_all[4];
  cCudaData<double     , F_CFLOAT                   , x>* cu_foriginal;
  int force_flag;
  int nlevels_respa;
};

}

#endif
#endif
