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

extern "C" void Cuda_FixNHCuda_Init(cuda_shared_data* sdata, X_CFLOAT dtv, V_CFLOAT dtf);
extern "C" void Cuda_FixNHCuda_nh_v_press(cuda_shared_data* sdata, int groupbit, double* factor_h, int mynlocal, int p_triclinic); //mynlocal can be nfirst if firstgroup==igroup  see cpp
extern "C" void Cuda_FixNHCuda_nh_v_temp(cuda_shared_data* sdata, int groupbit, F_CFLOAT factor_eta, int mynlocal); //mynlocal can be nfirst if firstgroup==igroup  see cpp
extern "C" void Cuda_FixNHCuda_nh_v_press_and_nve_v_NoBias(cuda_shared_data* sdata, int groupbit, double* factor_h, int mynlocal, int p_triclinic); //mynlocal can be nfirst if firstgroup==igroup  see cpp
extern "C" void Cuda_FixNHCuda_nve_v(cuda_shared_data* sdata, int groupbit, int mynlocal); //mynlocal can be nfirst if firstgroup==igroup  see cpp
extern "C" void Cuda_FixNHCuda_nve_x(cuda_shared_data* sdata, int groupbit, int mynlocal); //mynlocal can be nfirst if firstgroup==igroup  see cpp
extern "C" void Cuda_FixNHCuda_nve_v_and_nh_v_press_NoBias(cuda_shared_data* sdata, int groupbit, double* factor_h, int mynlocal, int p_triclinic); //mynlocal can be nfirst if firstgroup==igroup  see cpp
