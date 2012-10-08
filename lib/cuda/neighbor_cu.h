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

#ifndef NEIGHBOR_CU_H_
#define NEIGHBOR_CU_H_
#include "cuda_shared.h"

extern "C" int Cuda_BinAtoms(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist);
extern "C" int Cuda_NeighborBuildFullBin(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist);
extern "C" int Cuda_NeighborBuildFullNsq(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist);

#endif /*NEIGHBOR_CU_H_*/
