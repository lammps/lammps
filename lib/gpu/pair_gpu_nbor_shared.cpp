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

#include "pair_gpu_nbor_shared.h"

#ifdef USE_OPENCL
#include "pair_gpu_nbor_cl.h"
#else
#include "pair_gpu_nbor_ptx.h"
#include "pair_gpu_build_ptx.h"
#endif
  
void PairGPUNborShared::clear() {
  if (_compiled) {
    if (_gpu_nbor) {
      k_cell_id.clear();
      k_cell_counts.clear();
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

void PairGPUNborShared::compile_kernels(UCL_Device &dev, const bool gpu_nbor) {
  if (_compiled)
  	return;
  	
  _gpu_nbor=gpu_nbor;
  std::string flags="-cl-fast-relaxed-math -cl-mad-enable";

  if (gpu_nbor==false) {
    nbor_program=new UCL_Program(dev);
    nbor_program->load_string(pair_gpu_nbor_kernel,flags.c_str());
    k_nbor.set_function(*nbor_program,"kernel_unpack");
  } else {
    build_program=new UCL_Program(dev);
    #ifdef USE_OPENCL
    std::cerr << "CANNOT CURRENTLY USE GPU NEIGHBORING WITH OPENCL\n";
    exit(1);
    #else
    build_program->load_string(pair_gpu_build_kernel,flags.c_str());
    #endif
    k_cell_id.set_function(*build_program,"calc_cell_id");
    k_cell_counts.set_function(*build_program,"kernel_calc_cell_counts");
    k_build_nbor.set_function(*build_program,"calc_neigh_list_cell");
    k_transpose.set_function(*build_program,"transpose");
    k_special.set_function(*build_program,"kernel_special");
    neigh_tex.get_texture(*build_program,"neigh_tex");
  }
  _compiled=true;
}
