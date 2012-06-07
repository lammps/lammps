/* -*- c++ -*- ----------------------------------------------------------
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

#ifndef LMP_NEIGH_LIST_CUDA_H
#define LMP_NEIGH_LIST_CUDA_H

#include "pointers.h"
#include "cuda_data.h"
#include "neigh_list.h"

namespace LAMMPS_NS
{

class CudaNeighList : protected Pointers
{
        public:
                cCudaData<int , int , x>*  cu_ilist;
                cCudaData<int , int , x>*  cu_numneigh;
                cCudaData<int , int , x>*  cu_inum_border;
                cCudaData<int , int , x>*  cu_ilist_border;
                cCudaData<int , int , x>*  cu_numneigh_border;
                cCudaData<int , int , x>*  cu_numneigh_inner;
                cCudaData<int , int , x>*  cu_neighbors;
                cCudaData<int , int , x>*  cu_neighbors_border;
                cCudaData<int , int , x>*  cu_neighbors_inner;
                cCudaData<int , int , x>*  cu_ex_type;
                cCudaData<int , int , x>*  cu_ex1_bit;
                cCudaData<int , int , x>*  cu_ex2_bit;
                cCudaData<int , int , x>*  cu_ex_mol_bit;


                cuda_shared_neighlist sneighlist;

                int* neighbors;
                int* neighbors_inner;
                int* neighbors_border;
                int inum_border;
                int* ilist_border;
                int* numneigh_border;
                int* numneigh_inner;
                int nex_type;
                int nex_group;
                int nex_mol;

                bool build_cuda;

                CudaNeighList(class LAMMPS *, class NeighList* neigh_list);
                ~CudaNeighList();
                void grow_device(); // will grow pages memory on device, keeping old pages. will grow lists memory on device, deleting old lists
                void nl_upload(bool will_be_changed=true);
                void nl_download(bool will_be_changed=true);
                NeighList* neigh_list;

                void dev_alloc();
                void dev_free();

 private:
  class Cuda *cuda;
};

}

#endif
