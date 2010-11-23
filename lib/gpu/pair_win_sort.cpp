/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
 
/* ----------------------------------------------------------------------
   Contributing authors: Mike Brown (ORNL), brownw@ornl.gov
------------------------------------------------------------------------- */

#ifndef USE_OPENCL
#include "cudpp.h"
#endif

class PairWinSort {
 public:
  inline PairWinSort() : _allocated(false) {
    #ifndef USE_OPENCL
    sort_config.op = CUDPP_ADD;
    sort_config.datatype = CUDPP_UINT;
    sort_config.algorithm = CUDPP_SORT_RADIX;
    sort_config.options = CUDPP_OPTION_KEY_VALUE_PAIRS;
    #endif
  }
  inline ~PairWinSort() { clear(); }

  /// Free all memory on host and device
  inline void clear() {
    #ifndef USE_OPENCL
    if (_allocated) { cudppDestroyPlan(sort_plan); _allocated=false; }
    #endif
  }
 
  inline bool alloc(const int max_atoms) {
    #ifndef USE_OPENCL
    clear();
    CUDPPResult result = cudppPlan(&sort_plan, sort_config, max_atoms, 1, 0);  
    if (CUDPP_SUCCESS != result)
      return false;
    #endif
    return true;
  }

  /// Sort arrays for neighbor list calculation
  void sort_neighbor(const int num_atoms, unsigned *cell_begin, int *particle_begin) {
    #ifndef USE_OPENCL
    CUDPPResult result = cudppSort(sort_plan, cell_begin, particle_begin, 
                                   8*sizeof(unsigned), num_atoms);
    if (CUDPP_SUCCESS != result) {
      printf("Error in cudppSort\n");
      assert(1==0);
    }
    #endif
  }
  
 private:
  
  bool allocated;

  #ifndef USE_OPENCL
  CUDPPConfiguration sort_config;
  CUDPPHandle sort_plan;
  #endif
};

static PairWinSort win_sort;

extern "C" __declspec(dllexport) bool _win_sort_alloc(const int max_atoms) {
  win_sort.alloc(max_atoms);
}

extern "C" __declspec(dllexport) bool _win_sort(const int max_atoms, unsigned *cell_begin,
                                                int *particle_begin) {
  win_sort.sort(num_atoms,cell_begin,particle_begin);
}
