/* ----------------------------------------------------------------------
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

FixStyle(viscous/cuda,FixViscousCuda)

#else

#ifndef LMP_FIX_VISCOUS_CUDA_H
#define LMP_FIX_VISCOUS_CUDA_H

#include "fix_viscous.h"
#include "cuda_data.h"

namespace LAMMPS_NS {

class FixViscousCuda : public FixViscous {
 public:
  FixViscousCuda(class LAMMPS *, int, char **);
  ~FixViscousCuda();
  int setmask();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  cCudaData<double, F_FLOAT, x>* cu_gamma;

  private:
  class Cuda *cuda;
};

}

#endif
#endif
