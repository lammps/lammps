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

extern "C" void Cuda_FixNVECuda_Init(cuda_shared_data* sdata, X_FLOAT dtv, V_FLOAT dtf);
extern "C" void Cuda_FixNVECuda_InitialIntegrate(cuda_shared_data* sdata, int groupbit, int mynlocal);
extern "C" void Cuda_FixNVECuda_FinalIntegrate(cuda_shared_data* sdata, int groupbit, int mynlocal);
