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

struct ParamSW_Float {
  F_FLOAT epsilon, sigma;
  F_FLOAT littlea, lambda, gamma, costheta;
  F_FLOAT biga, bigb;
  F_FLOAT powerp, powerq;
  F_FLOAT tol;
  F_FLOAT cut, cutsq;
  F_FLOAT sigma_gamma, lambda_epsilon, lambda_epsilon2;
  F_FLOAT c1, c2, c3, c4, c5, c6;
  int ielement, jelement, kelement;
};

extern "C" void Cuda_PairSWCuda_Init(cuda_shared_data* sdata, ParamSW_Float* params_host, void* map_host, void* elem2param_host, int nelements_h);
extern "C" void Cuda_PairSWCuda(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist, int eflag, int vflag, int eflag_atom, int vflag_atom);
