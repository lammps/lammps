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

FixStyle(gravity/cuda,FixGravityCuda)

#else

#ifndef LMP_FIX_GRAVITY_CUDA_H
#define LMP_FIX_GRAVITY_CUDA_H

#include "fix.h"
#include "cuda_data.h"

namespace LAMMPS_NS {

class FixGravityCuda : public Fix {
 public:
  FixGravityCuda(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void setup(int);
  void post_force(int);

 private:
  class Cuda *cuda;
  int style;
  double magnitude,dt;
  double phi,theta,phigrad,thetagrad;
  double xdir,ydir,zdir;
  double xgrav,ygrav,zgrav,xacc,yacc,zacc;
  double degree2rad;
  int time_origin;
};

}

#endif
#endif
