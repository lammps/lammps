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

struct Param_Float {
  F_FLOAT lam1, lam2, lam3;
  F_FLOAT c, d, h;
  F_FLOAT gamma, powerm;
  F_FLOAT powern, beta;
  F_FLOAT biga, bigb, bigd, bigr;
  F_FLOAT cut, cutsq;
  F_FLOAT c1, c2, c3, c4;
  int ielement, jelement, kelement;
  int powermint;
  //F_FLOAT Z_i,Z_j;
  F_FLOAT ZBLcut, ZBLexpscale;
  F_FLOAT a_ij, premult;
};

extern "C" void Cuda_PairTersoffCuda_Init(cuda_shared_data* sdata, Param_Float* params_host, void* map_host, void* elem2param_host, int nelements_h, bool zbl);
extern "C" void Cuda_PairTersoffCuda(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist, int eflag, int vflag, int eflag_atom, int vflag_atom);
