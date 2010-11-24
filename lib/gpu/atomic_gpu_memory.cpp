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
 
#include "atomic_gpu_memory.h"
#define AtomicGPUMemoryT AtomicGPUMemory<numtyp, acctyp>

extern PairGPUDevice<PRECISION,ACC_PRECISION> pair_gpu_device;

template <class numtyp, class acctyp>
AtomicGPUMemoryT::AtomicGPUMemory() : _compiled(false), _max_bytes(0)  {
  device=&pair_gpu_device;
}

template <class numtyp, class acctyp>
AtomicGPUMemoryT::~AtomicGPUMemory() {
}

template <class numtyp, class acctyp>
int AtomicGPUMemoryT::bytes_per_atom_atomic(const int max_nbors) const {
  return device->atom.bytes_per_atom()+device->nbor.bytes_per_atom(max_nbors);
}

template <class numtyp, class acctyp>
bool AtomicGPUMemoryT::init_atomic(const int nlocal, const int nall,
                                   const int max_nbors, const int maxspecial,
                                   const double cell_size,
                                   const double gpu_split, FILE *_screen,
                                   const char *pair_program) {
  nbor_time_avail=false;
  screen=_screen;

  bool gpu_nbor=false;
  if (device->gpu_mode()==PairGPUDevice<numtyp,acctyp>::GPU_NEIGH)
    gpu_nbor=true;

  int _gpu_host=0;
  int host_nlocal=hd_balancer.first_host_count(nlocal,gpu_nbor,gpu_split);
  if (host_nlocal>0)
    _gpu_host=1;

  if (!device->init(false,false,nlocal,host_nlocal,nall,maxspecial,gpu_nbor,
                    _gpu_host,max_nbors,cell_size,false))
    return false;
  ucl_device=device->gpu;
  atom=&device->atom;
  nbor=&device->nbor;

  _block_size=BLOCK_1D;
  if (static_cast<size_t>(_block_size)>ucl_device->group_size())
    _block_size=ucl_device->group_size();
  compile_kernels(*ucl_device,pair_program);

  // Initialize host-device load balancer
  hd_balancer.init(device,gpu_split);

  // Initialize timers for the selected GPU
  time_pair.init(*ucl_device);
  time_pair.zero();

  pos_tex.bind_float(atom->dev_x,4);

  _max_an_bytes=atom->gpu_bytes()+nbor->gpu_bytes();

  return true;
}

template <class numtyp, class acctyp>
void AtomicGPUMemoryT::clear_atomic() {
  // Output any timing information
  acc_timers();
  double avg_split=hd_balancer.all_avg_split();
  device->output_times(time_pair,avg_split,_max_bytes+_max_an_bytes,screen);

  if (_compiled) {
    k_pair_fast.clear();
    k_pair.clear();
    delete pair_program;
    _compiled=false;
  }

  time_pair.clear();
  hd_balancer.clear();

  device->clear();
}

// ---------------------------------------------------------------------------
// Copy neighbor list from host
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int * AtomicGPUMemoryT::reset_nbors(const int nall, const int inum, int *ilist,
                                   int *numj, int **firstneigh, bool &success) {
  success=true;

  nbor_time_avail=true;

  int mn=nbor->max_nbor_loop(inum,numj);
  resize_atom(inum,nall,success);
  resize_local(inum,mn,success);
  if (!success)
    return false;

  nbor->get_host(inum,ilist,numj,firstneigh,block_size());

  double bytes=atom->gpu_bytes()+nbor->gpu_bytes();
  if (bytes>_max_an_bytes)
    _max_an_bytes=bytes;
  
  return ilist;
}

// ---------------------------------------------------------------------------
// Build neighbor list on device
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
inline void AtomicGPUMemoryT::build_nbor_list(const int inum,
                                              const int host_inum,
                                              const int nall, double **host_x,
                                              int *host_type, double *boxlo,
                                              double *boxhi, int *tag,
                                              int **nspecial, int **special,
                                              bool &success) {
  nbor_time_avail=true;

  success=true;
  resize_atom(inum,nall,success);
  resize_local(inum,host_inum,nbor->max_nbors(),success);
  if (!success)
    return;
  atom->cast_copy_x(host_x,host_type);

  int mn;
  nbor->build_nbor_list(inum, host_inum, nall, *atom, boxlo, boxhi, tag,
                        nspecial, special, success, mn);

  double bytes=atom->gpu_bytes()+nbor->gpu_bytes();
  if (bytes>_max_an_bytes)
    _max_an_bytes=bytes;
}

// ---------------------------------------------------------------------------
// Copy nbor list from host if necessary and then calculate forces, virials,..
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void AtomicGPUMemoryT::compute(const int timestep, const int f_ago,
			      const int inum_full, const int nall,
                              double **host_x, int *host_type,
                              int *ilist, int *numj, int **firstneigh,
                              const bool eflag, const bool vflag,
                              const bool eatom, const bool vatom,
                              int &host_start, const double cpu_time,
                              bool &success) {
  acc_timers();
  if (inum_full==0) {
    zero_timers();
    return;
  }
  
  int ago=hd_balancer.ago_first(f_ago);
  int inum=hd_balancer.balance(timestep,ago,inum_full,cpu_time,
		               nbor->gpu_nbor());
  atom->inum(inum);
  host_start=inum;

  if (ago==0) {
    reset_nbors(nall, inum, ilist, numj, firstneigh, success);
    if (!success)
      return;
  }

  atom->cast_x_data(host_x,host_type);
  hd_balancer.start_timer();
  atom->add_x_data(host_x,host_type);

  loop(eflag,vflag);
  atom->copy_answers(eflag,vflag,eatom,vatom,ilist);
  hd_balancer.stop_timer();
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU if necessary and then compute forces, virials, energies
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int * AtomicGPUMemoryT::compute(const int timestep, const int ago,
                                const int inum_full, const int nall,
                                double **host_x, int *host_type, double *boxlo,
                                double *boxhi, int *tag, int **nspecial,
                                int **special, const bool eflag, 
                                const bool vflag, const bool eatom,
                                const bool vatom, int &host_start,
                                const double cpu_time, bool &success) {
  acc_timers();
  if (inum_full==0) {
    zero_timers();
    return NULL;
  }
  
  hd_balancer.balance(cpu_time,nbor->gpu_nbor());
  int inum=hd_balancer.get_gpu_count(timestep,ago,inum_full);
  atom->inum(inum);
  host_start=inum;
 
  // Build neighbor list on GPU if necessary
  if (ago==0) {
    build_nbor_list(inum, inum_full-inum, nall, host_x, host_type,
                    boxlo, boxhi, tag, nspecial, special, success);
    if (!success)
      return NULL;
    hd_balancer.start_timer();
  } else {
    atom->cast_x_data(host_x,host_type);
    hd_balancer.start_timer();
    atom->add_x_data(host_x,host_type);
  }

  loop(eflag,vflag);
  atom->copy_answers(eflag,vflag,eatom,vatom);
  hd_balancer.stop_timer();
  
  return device->nbor.host_nbor.begin();
}

template <class numtyp, class acctyp>
double AtomicGPUMemoryT::host_memory_usage_atomic() const {
  return device->atom.host_memory_usage()+
         device->nbor.host_memory_usage()+4*sizeof(numtyp)+
         sizeof(AtomicGPUMemory<numtyp,acctyp>);
}

template <class numtyp, class acctyp>
void AtomicGPUMemoryT::compile_kernels(UCL_Device &dev, const char *pair_str) {
  if (_compiled)
    return;

  std::string flags="-cl-fast-relaxed-math -cl-mad-enable "+
                    std::string(OCL_PRECISION_COMPILE);

  pair_program=new UCL_Program(dev);
  pair_program->load_string(pair_str,flags.c_str());
  k_pair_fast.set_function(*pair_program,"kernel_pair_fast");
  k_pair.set_function(*pair_program,"kernel_pair");
  pos_tex.get_texture(*pair_program,"pos_tex");

  _compiled=true;
}

template class AtomicGPUMemory<PRECISION,ACC_PRECISION>;

