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

extern "C" int Cuda_CommCuda_PackComm(cuda_shared_data* sdata, int n, int iswap, void* buf_send, int* pbc, int pbcflag);
extern "C" int Cuda_CommCuda_PackCommVel(cuda_shared_data* sdata, int n, int iswap, void* buf_send, int* pbc, int pbcflag);
extern "C" int Cuda_CommCuda_PackComm_Self(cuda_shared_data* sdata, int n, int iswap, int first, int* pbc, int pbcflag);
extern "C" int Cuda_CommCuda_PackCommVel_Self(cuda_shared_data* sdata, int n, int iswap, int first, int* pbc, int pbcflag);
extern "C" void Cuda_CommCuda_UnpackComm(cuda_shared_data* sdata, int n, int first, void* buf_recv, int iswap = -1);
extern "C" void Cuda_CommCuda_UnpackCommVel(cuda_shared_data* sdata, int n, int first, void* buf_recv, int iswap = -1);
extern "C" int Cuda_CommCuda_PackReverse(cuda_shared_data* sdata, int n, int first, void* buf_send);
extern "C" void Cuda_CommCuda_UnpackReverse(cuda_shared_data* sdata, int n, int iswap, void* buf_recv);
extern "C" void Cuda_CommCuda_UnpackReverse_Self(cuda_shared_data* sdata, int n, int iswap, int first);
extern "C" int Cuda_CommCuda_BuildSendlist(cuda_shared_data* sdata, int bordergroup, int ineed, int style, int atom_nfirst, int nfirst, int nlast, int dim, int iswap);
