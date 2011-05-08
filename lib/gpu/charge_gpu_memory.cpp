/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Charge/Molecular Massively Parallel Simulator
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

#include "charge_gpu_memory.h"
#define ChargeGPUMemoryT ChargeGPUMemory<numtyp, acctyp>

extern PairGPUDevice<PRECISION,ACC_PRECISION> pair_gpu_device;

template <class numtyp, class acctyp>
ChargeGPUMemoryT::ChargeGPUMemory() : _compiled(false), _max_bytes(0) {
  device=&pair_gpu_device;
  ans=new PairGPUAns<numtyp,acctyp>();
  nbor=new PairGPUNbor();
}

template <class numtyp, class acctyp>
ChargeGPUMemoryT::~ChargeGPUMemory() {
  delete ans;
  delete nbor;
}

template <class numtyp, class acctyp>
int ChargeGPUMemoryT::bytes_per_atom_atomic(const int max_nbors) const {
  return device->atom.bytes_per_atom()+ans->bytes_per_atom()+
         nbor->bytes_per_atom(max_nbors);
}

template <class numtyp, class acctyp>
int ChargeGPUMemoryT::init_atomic(const int nlocal, const int nall,
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
  int host_nlocal=hd_balancer.first_host_count(nlocal,gpu_split,gpu_nbor);
  if (host_nlocal>0)
    _gpu_host=1;

  _threads_per_atom=device->threads_per_charge();
  if (_threads_per_atom>1 && gpu_nbor==false) {
    nbor->packing(true);
    _nbor_data=&(nbor->dev_packed);
  } else
    _nbor_data=&(nbor->dev_nbor);
    
  int success=device->init(*ans,true,false,nlocal,host_nlocal,nall,nbor,
                           maxspecial,_gpu_host,max_nbors,cell_size,false);
  if (success!=0)
    return success;

  ucl_device=device->gpu;
  atom=&device->atom;

  _block_size=device->pair_block_size();
  _block_bio_size=device->block_bio_pair();
  compile_kernels(*ucl_device,pair_program);

  // Initialize host-device load balancer
  hd_balancer.init(device,gpu_nbor,gpu_split);

  // Initialize timers for the selected GPU
  time_pair.init(*ucl_device);
  time_pair.zero();

  pos_tex.bind_float(atom->dev_x,4);
  q_tex.bind_float(atom->dev_q,1);

  _max_an_bytes=ans->gpu_bytes()+nbor->gpu_bytes();

  return success;
}

template <class numtyp, class acctyp>
void ChargeGPUMemoryT::estimate_gpu_overhead() {
  device->estimate_gpu_overhead(1,_gpu_overhead,_driver_overhead);
}

template <class numtyp, class acctyp>
void ChargeGPUMemoryT::clear_atomic() {
  // Output any timing information
  acc_timers();
  double avg_split=hd_balancer.all_avg_split();
  _gpu_overhead*=hd_balancer.timestep();
  _driver_overhead*=hd_balancer.timestep();
  device->output_times(time_pair,*ans,*nbor,avg_split,_max_bytes+_max_an_bytes,
                       _gpu_overhead,_driver_overhead,_threads_per_atom,screen);

  if (_compiled) {
    k_pair_fast.clear();
    k_pair.clear();
    delete pair_program;
    _compiled=false;
  }

  time_pair.clear();
  hd_balancer.clear();

  nbor->clear();
  ans->clear();
  device->clear();
}

// ---------------------------------------------------------------------------
// Copy neighbor list from host
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int * ChargeGPUMemoryT::reset_nbors(const int nall, const int inum, int *ilist,
                                   int *numj, int **firstneigh, bool &success) {
  success=true;

  nbor_time_avail=true;

  int mn=nbor->max_nbor_loop(inum,numj,ilist);
  resize_atom(inum,nall,success);
  resize_local(inum,mn,success);
  if (!success)
    return false;

  nbor->get_host(inum,ilist,numj,firstneigh,block_size());

  double bytes=ans->gpu_bytes()+nbor->gpu_bytes();
  if (bytes>_max_an_bytes)
    _max_an_bytes=bytes;

  return ilist;
}

// ---------------------------------------------------------------------------
// Build neighbor list on device
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
inline void ChargeGPUMemoryT::build_nbor_list(const int inum,
                                              const int host_inum,
                                              const int nall, double **host_x,
                                              int *host_type, double *sublo,
                                              double *subhi, int *tag, 
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
  nbor->build_nbor_list(inum, host_inum, nall, *atom, sublo, subhi, tag,
                        nspecial, special, success, mn);

  double bytes=ans->gpu_bytes()+nbor->gpu_bytes();
  if (bytes>_max_an_bytes)
    _max_an_bytes=bytes;
}

// ---------------------------------------------------------------------------
// Copy nbor list from host if necessary and then calculate forces, virials,..
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void ChargeGPUMemoryT::compute(const int f_ago, const int inum_full,
                               const int nall, double **host_x, int *host_type,
                               int *ilist, int *numj, int **firstneigh,
                               const bool eflag, const bool vflag,
                               const bool eatom, const bool vatom,
                               int &host_start, const double cpu_time,
                               bool &success, double *host_q,
                               const int nlocal, double *boxlo, double *prd) {
  acc_timers();
  if (inum_full==0) {
    host_start=0;
    // Make sure textures are correct if realloc by a different hybrid style
    resize_atom(0,nall,success);
    zero_timers();
    return;
  }
  
  int ago=hd_balancer.ago_first(f_ago);
  int inum=hd_balancer.balance(ago,inum_full,cpu_time);
  ans->inum(inum);
  host_start=inum;

  if (ago==0) {
    reset_nbors(nall, inum, ilist, numj, firstneigh, success);
    if (!success)
      return;
  }

  atom->cast_x_data(host_x,host_type);
  atom->cast_q_data(host_q);
  hd_balancer.start_timer();
  atom->add_x_data(host_x,host_type);
  atom->add_q_data();

  device->precompute(f_ago,nlocal,nall,host_x,host_type,success,host_q,
                     boxlo, prd);

  loop(eflag,vflag);
  ans->copy_answers(eflag,vflag,eatom,vatom,ilist);
  device->add_ans_object(ans);
  hd_balancer.stop_timer();
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU if necessary and then compute forces, virials, energies
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int** ChargeGPUMemoryT::compute(const int ago, const int inum_full,
                                const int nall, double **host_x, int *host_type,
                                double *sublo, double *subhi, int *tag,
                                int **nspecial, int **special, const bool eflag, 
                                const bool vflag, const bool eatom,
                                const bool vatom, int &host_start,
                                int **ilist, int **jnum,
                                const double cpu_time, bool &success,
                                double *host_q, double *boxlo, double *prd) {
  acc_timers();
  if (inum_full==0) {
    host_start=0;
    // Make sure textures are correct if realloc by a different hybrid style
    resize_atom(0,nall,success);
    zero_timers();
    return NULL;
  }
  
  hd_balancer.balance(cpu_time);
  int inum=hd_balancer.get_gpu_count(ago,inum_full);
  ans->inum(inum);
  host_start=inum;
 
  // Build neighbor list on GPU if necessary
  if (ago==0) {
    build_nbor_list(inum, inum_full-inum, nall, host_x, host_type,
                    sublo, subhi, tag, nspecial, special, success);
    if (!success)
      return NULL;
    atom->cast_q_data(host_q);
    hd_balancer.start_timer();
  } else {
    atom->cast_x_data(host_x,host_type);
    atom->cast_q_data(host_q);
    hd_balancer.start_timer();
    atom->add_x_data(host_x,host_type);
  }
  atom->add_q_data();
  *ilist=nbor->host_ilist.begin();
  *jnum=nbor->host_acc.begin();

  device->precompute(ago,inum_full,nall,host_x,host_type,success,host_q,
                     boxlo, prd);

  loop(eflag,vflag);
  ans->copy_answers(eflag,vflag,eatom,vatom);
  device->add_ans_object(ans);
  hd_balancer.stop_timer();
  
  return nbor->host_jlist.begin()-host_start;
}

template <class numtyp, class acctyp>
double ChargeGPUMemoryT::host_memory_usage_atomic() const {
  return device->atom.host_memory_usage()+nbor->host_memory_usage()+
         4*sizeof(numtyp)+sizeof(ChargeGPUMemory<numtyp,acctyp>);
}

template <class numtyp, class acctyp>
void ChargeGPUMemoryT::compile_kernels(UCL_Device &dev, const char *pair_str) {
  if (_compiled)
    return;

  std::string flags="-cl-fast-relaxed-math -cl-mad-enable "+
                    std::string(OCL_PRECISION_COMPILE);

  pair_program=new UCL_Program(dev);
  pair_program->load_string(pair_str,flags.c_str());
  k_pair_fast.set_function(*pair_program,"kernel_pair_fast");
  k_pair.set_function(*pair_program,"kernel_pair");
  pos_tex.get_texture(*pair_program,"pos_tex");
  q_tex.get_texture(*pair_program,"q_tex");

  _compiled=true;
}

template class ChargeGPUMemory<PRECISION,ACC_PRECISION>;

