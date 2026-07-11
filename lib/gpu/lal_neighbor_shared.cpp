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

#include <cmath>
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

double NeighborShared::best_cell_size(const double subx, const double suby,
                                      const double subz, const int nlocal,
                                      const double cut) {
  if (_cached_cell_size && _cut_sort==cut) {
    _cached_cell_size=false;
    return _cell_size;
  }

  const double box_density = static_cast<double>(nlocal) / (subx*suby*subz);
  const double density=box_density*cut*cut*cut;
  if (density >= 4.0 * _simd_size) return cut*0.5;
  else if (density >= 0.5 * _simd_size) return cut;

  const double iters = 60;
  const double inc = cut/(iters-1);
  const double iss = 1.0 / _simd_size;
  double test_size = cut;
  double best_iters = 1e200;
  double best_size;
  for (int i = 0; i < iters; i++) {
    const double i_test_size = 1.0/test_size;
    const int ncellx = static_cast<int>(ceil(subx*i_test_size));
    const int ncelly = static_cast<int>(ceil(suby*i_test_size));
    const int ncellz = static_cast<int>(ceil(subz*i_test_size));
    const double density = box_density*test_size*test_size*test_size;
    const double iters_per_cell = ceil(iss*density);
    const double iters = ncellx*ncelly*ncellz*iters_per_cell*
      ceil(density*27.0*iss);
    if (iters < best_iters) {
      best_iters = iters;
      best_size = test_size;
    }
    test_size += inc;
  }
  const int cells_in_cutoff=static_cast<int>(ceil(cut/best_size));
  if (cells_in_cutoff > 2) best_size=cut*0.5;
  return best_size;
}

void NeighborShared::compile_kernels(UCL_Device &dev, const int gpu_nbor,
                                     const std::string &flags) {
  if (_compiled)
          return;

  _gpu_nbor=gpu_nbor;
  if (_gpu_nbor==0) {
    nbor_program=new UCL_Program(dev);
    nbor_program->load_string(neighbor_cpu,flags.c_str(),nullptr,stderr);
    k_nbor.set_function(*nbor_program,"kernel_unpack");
  } else {
    build_program=new UCL_Program(dev);
    build_program->load_string(neighbor_gpu,flags.c_str(),nullptr,stderr);

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
