/***************************************************************************
                               base_dpd.cpp
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Base class for pair styles needing per-particle data for position,
  dipole, and type.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : Jan 15, 2014
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#include "lal_base_dpd.h"
using namespace LAMMPS_AL;
#define BaseDPDT BaseDPD<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> global_device;

template <class numtyp, class acctyp>
BaseDPDT::BaseDPD() : _compiled(false), _max_bytes(0) {
  device=&global_device;
  ans=new Answer<numtyp,acctyp>();
  nbor=new Neighbor();
}

template <class numtyp, class acctyp>
BaseDPDT::~BaseDPD() {
  delete ans;
  delete nbor;
}

template <class numtyp, class acctyp>
int BaseDPDT::bytes_per_atom_atomic(const int max_nbors) const {
  return device->atom.bytes_per_atom()+ans->bytes_per_atom()+
         nbor->bytes_per_atom(max_nbors);
}

template <class numtyp, class acctyp>
int BaseDPDT::init_atomic(const int nlocal, const int nall,
                          const int max_nbors, const int maxspecial,
                          const double cell_size,
                          const double gpu_split, FILE *_screen,
                          const void *pair_program, const char *k_name) {
  screen=_screen;

  int gpu_nbor=0;
  if (device->gpu_mode()==Device<numtyp,acctyp>::GPU_NEIGH)
    gpu_nbor=1;
  else if (device->gpu_mode()==Device<numtyp,acctyp>::GPU_HYB_NEIGH)
    gpu_nbor=2;

  int _gpu_host=0;
  int host_nlocal=hd_balancer.first_host_count(nlocal,gpu_split,gpu_nbor);
  if (host_nlocal>0)
    _gpu_host=1;

  _threads_per_atom=device->threads_per_charge();
  if (_threads_per_atom>1 && gpu_nbor==0) {
    nbor->packing(true);
    _nbor_data=&(nbor->dev_packed);
  } else
    _nbor_data=&(nbor->dev_nbor);

  int success=device->init(*ans,false,false,nlocal,host_nlocal,nall,nbor,
                           maxspecial,_gpu_host,max_nbors,cell_size,false,
                           _threads_per_atom,true);
  if (success!=0)
    return success;

  success = device->init_nbor(nbor,nlocal,host_nlocal,nall,maxspecial,_gpu_host,
                  max_nbors,cell_size,false,_threads_per_atom);

  if (success!=0)
    return success;

  ucl_device=device->gpu;
  atom=&device->atom;

  _block_size=device->pair_block_size();
  _block_bio_size=device->block_bio_pair();
  compile_kernels(*ucl_device,pair_program,k_name);

  // Initialize host-device load balancer
  hd_balancer.init(device,gpu_nbor,gpu_split);

  // Initialize timers for the selected GPU
  time_pair.init(*ucl_device);
  time_pair.zero();

  pos_tex.bind_float(atom->x,4);
  vel_tex.bind_float(atom->v,4);

  _max_an_bytes=ans->gpu_bytes()+nbor->gpu_bytes();

  return success;
}

template <class numtyp, class acctyp>
void BaseDPDT::estimate_gpu_overhead() {
  device->estimate_gpu_overhead(1,_gpu_overhead,_driver_overhead);
}

template <class numtyp, class acctyp>
void BaseDPDT::clear_atomic() {
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
int * BaseDPDT::reset_nbors(const int nall, const int inum, int *ilist,
                            int *numj, int **firstneigh, bool &success) {
  success=true;

  int mn=nbor->max_nbor_loop(inum,numj,ilist);
  resize_atom(inum,nall,success);
  resize_local(inum,mn,success);
  if (!success)
    return NULL;

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
inline void BaseDPDT::build_nbor_list(const int inum, const int host_inum,
                                      const int nall, double **host_x,
                                      int *host_type, double *sublo,
                                      double *subhi, tagint *tag,
                                      int **nspecial, tagint **special,
                                      bool &success) {
  success=true;
  resize_atom(inum,nall,success);
  resize_local(inum,host_inum,nbor->max_nbors(),success);
  if (!success)
    return;
  atom->cast_copy_x(host_x,host_type);

  int mn;
  nbor->build_nbor_list(host_x, inum, host_inum, nall, *atom, sublo, subhi, tag,
                        nspecial, special, success, mn);

  double bytes=ans->gpu_bytes()+nbor->gpu_bytes();
  if (bytes>_max_an_bytes)
    _max_an_bytes=bytes;
}

// ---------------------------------------------------------------------------
// Copy nbor list from host if necessary and then calculate forces, virials,..
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void BaseDPDT::compute(const int f_ago, const int inum_full,
                       const int nall, double **host_x, int *host_type,
                       int *ilist, int *numj, int **firstneigh,
                       const bool eflag, const bool vflag,
                       const bool eatom, const bool vatom,
                       int &host_start, const double cpu_time,
                       bool &success, tagint *tag, double **host_v,
                       const double dtinvsqrt, const int seed, const int timestep,
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
  atom->cast_v_data(host_v,tag);
  hd_balancer.start_timer();
  atom->add_x_data(host_x,host_type);
  atom->add_v_data(host_v,tag);

  _dtinvsqrt = dtinvsqrt;
  _seed = seed;
  _timestep = timestep;

  loop(eflag,vflag);
  ans->copy_answers(eflag,vflag,eatom,vatom,ilist);
  device->add_ans_object(ans);
  hd_balancer.stop_timer();
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU if necessary and then compute forces, virials, energies
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int** BaseDPDT::compute(const int ago, const int inum_full,
                        const int nall, double **host_x, int *host_type,
                        double *sublo, double *subhi, tagint *tag,
                        int **nspecial, tagint **special, const bool eflag,
                        const bool vflag, const bool eatom,
                        const bool vatom, int &host_start,
                        int **ilist, int **jnum,
                        const double cpu_time, bool &success,
                        double **host_v, const double dtinvsqrt,
                        const int seed, const int timestep,
                        double *boxlo, double *prd) {
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
    atom->cast_v_data(host_v,tag);
    hd_balancer.start_timer();
  } else {
    atom->cast_x_data(host_x,host_type);
    atom->cast_v_data(host_v,tag);
    hd_balancer.start_timer();
    atom->add_x_data(host_x,host_type);
  }
  atom->add_v_data(host_v,tag);
  *ilist=nbor->host_ilist.begin();
  *jnum=nbor->host_acc.begin();

  _dtinvsqrt = dtinvsqrt;
  _seed = seed;
  _timestep = timestep;

  loop(eflag,vflag);
  ans->copy_answers(eflag,vflag,eatom,vatom);
  device->add_ans_object(ans);
  hd_balancer.stop_timer();

  return nbor->host_jlist.begin()-host_start;
}

template <class numtyp, class acctyp>
double BaseDPDT::host_memory_usage_atomic() const {
  return device->atom.host_memory_usage()+nbor->host_memory_usage()+
         4*sizeof(numtyp)+sizeof(BaseDPD<numtyp,acctyp>);
}

template <class numtyp, class acctyp>
void BaseDPDT::compile_kernels(UCL_Device &dev, const void *pair_str,
                                  const char *kname) {
  if (_compiled)
    return;

  std::string s_fast=std::string(kname)+"_fast";
  pair_program=new UCL_Program(dev);
  pair_program->load_string(pair_str,device->compile_string().c_str());
  k_pair_fast.set_function(*pair_program,s_fast.c_str());
  k_pair.set_function(*pair_program,kname);
  pos_tex.get_texture(*pair_program,"pos_tex");
  vel_tex.get_texture(*pair_program,"vel_tex");

  _compiled=true;
}

template class BaseDPD<PRECISION,ACC_PRECISION>;

