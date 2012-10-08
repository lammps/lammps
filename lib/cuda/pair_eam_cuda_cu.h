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
extern "C" void Cuda_PairEAMCuda_Init(cuda_shared_data* sdata, double rdr, double rdrho, int nfrho, int nrhor, int nr, int nrho, int nz2r,
                                      void* frho_spline, void* rhor_spline, void* z2r_spline, void* rho, void* fp,
                                      int* type2frho, int** type2z2r, int** type2rhor);
extern "C" void Cuda_PairEAM1Cuda(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist, int eflag, int vflag, int eflag_atom, int vflag_atom);
extern "C" void Cuda_PairEAM2Cuda(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist, int eflag, int vflag, int eflag_atom, int vflag_atom);
extern "C" void Cuda_PairEAMCuda_PackComm(cuda_shared_data* sdata, int n, int iswap, void* buf_send);
extern "C" void Cuda_PairEAMCuda_UnpackComm(cuda_shared_data* sdata, int n, int first, void* buf_recv, void* fp);

#define EAM_COEFF_LENGTH 8
