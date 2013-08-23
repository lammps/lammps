/***************************************************************************
                             neighbor_shared.cpp
                             -------------------
                            W. Michael Brown (ORNL)

  Class for management of data shared by all neighbor lists

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : brownw@ornl.gov
 ***************************************************************************/

#include "lal_precision.h"
#include "lal_neighbor_shared.h"

#if defined(USE_OPENCL)
#include "neighbor_cpu_cl.h"
#include "neighbor_gpu_cl.h"
#elif defined(USE_CUDART)
const char *neighbor_cpu=0;
const char *neighbor_gpu=0;
#else
#include "neighbor_cpu_cubin.h"
#include "neighbor_gpu_cubin.h"
#endif

using namespace LAMMPS_AL;

void NeighborShared::clear() {
  if (_compiled) {
    if (_gpu_nbor>0) {
      if (_gpu_nbor==1) {
        k_cell_id.clear();
        k_cell_counts.clear();
      }
      k_build_nbor.clear();
      k_transpose.clear();
      k_special.clear();
      delete build_program;
    } else {
      k_nbor.clear();
      delete nbor_program;
    }
    _compiled=false;
  }
}

void NeighborShared::compile_kernels(UCL_Device &dev, const int gpu_nbor,
                                     const std::string flags) {
  if (_compiled)
  	return;
  	
  _gpu_nbor=gpu_nbor;
  if (_gpu_nbor==0) {
    nbor_program=new UCL_Program(dev);
    nbor_program->load_string(neighbor_cpu,flags.c_str());
    k_nbor.set_function(*nbor_program,"kernel_unpack");
  } else {
    build_program=new UCL_Program(dev);
    build_program->load_string(neighbor_gpu,flags.c_str());

    if (_gpu_nbor==1) {
      k_cell_id.set_function(*build_program,"calc_cell_id");
      k_cell_counts.set_function(*build_program,"kernel_calc_cell_counts");
    }
    k_build_nbor.set_function(*build_program,"calc_neigh_list_cell");
    k_transpose.set_function(*build_program,"transpose");
    k_special.set_function(*build_program,"kernel_special");
    neigh_tex.get_texture(*build_program,"pos_tex");
  }
  _compiled=true;
}
