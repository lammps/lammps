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

FixStyle(temp/rescale/limit/cuda,FixTempRescaleLimitCuda)

#else

#ifndef FIX_TEMP_RESCALE_LIMIT_CUDA_H
#define FIX_TEMP_RESCALE_LIMIT_CUDA_H

#include "fix.h"

namespace LAMMPS_NS {
class FixTempRescaleLimitCuda : public Fix {
 public:
  FixTempRescaleLimitCuda(class LAMMPS *, int, char **);
  ~FixTempRescaleLimitCuda();
  int setmask();
  void init();
  void end_of_step();
  int modify_param(int, char **);
  void reset_target(double);
  double compute_scalar();

 private:
  class Cuda *cuda;
  int which;
  double t_start,t_stop,t_window;
  double fraction,energy,efactor;
  double limit;
  char *id_temp;
  class Compute *temperature;
  int tflag;
};

}

#endif
#endif
