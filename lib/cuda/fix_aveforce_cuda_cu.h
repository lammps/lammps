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

#include "cuda_shared.h"

extern "C" void Cuda_FixAveForceCuda_Init(cuda_shared_data* sdata);
extern "C" void Cuda_FixAveForceCuda_PostForce_FOrg(cuda_shared_data* sdata, int groupbit, F_CFLOAT* aforiginal);
extern "C" void Cuda_FixAveForceCuda_PostForce_Set(cuda_shared_data* sdata, int groupbit, int xflag, int yflag, int zflag, F_CFLOAT axvalue, F_CFLOAT ayvalue, F_CFLOAT azvalue);
