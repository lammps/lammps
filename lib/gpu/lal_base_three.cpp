/***************************************************************************
                               base_three.cpp
                             -------------------
                            W. Michael Brown (ORNL)

  Base class for pair styles with per-particle data for position and type

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : Tue April 2, 2013
    email                : brownw@ornl.gov
 ***************************************************************************/

#include "lal_base_three.h"
namespace LAMMPS_AL {
#define BaseThreeT BaseThree<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> global_device;

template <class numtyp, class acctyp>
BaseThreeT::BaseThree() : _compiled(false), _max_bytes(0), _onetype(-1) {
  device=&global_device;
  ans=new Answer<numtyp,acctyp>();
  nbor=new Neighbor();
  #ifdef THREE_CONCURRENT
  ans2=new Answer<numtyp,acctyp>();
  #endif
  pair_program=nullptr;
  ucl_device=nullptr;
  #if defined(LAL_OCL_EV_JIT)
  pair_program_noev=nullptr;
  #endif
}

template <class numtyp, class acctyp>
BaseThreeT::~BaseThree() {
  delete ans;
  delete nbor;
  #ifdef THREE_CONCURRENT
  delete ans2;
  #endif
  k_three_center.clear();
  k_three_end.clear();
  k_three_end_vatom.clear();
  k_pair.clear();
  k_short_nbor.clear();
  if (pair_program) delete pair_program;
  #if defined(LAL_OCL_EV_JIT)
  k_three_center_noev.clear();
  k_three_end_noev.clear();
  k_pair_noev.clear();
  if (pair_program_noev) delete pair_program_noev;
  #endif
}

template <class numtyp, class acctyp>
int BaseThreeT::bytes_per_atom_atomic(const int max_nbors) const {
  int b=device->atom.bytes_per_atom()+ans->bytes_per_atom()+
    nbor->bytes_per_atom(max_nbors);
  #ifdef THREE_CONCURRENT
  b+=ans2->bytes_per_atom();
  #endif
  return b;
}

template <class numtyp, class acctyp>
int BaseThreeT::init_three(const int nlocal, const int nall,
                           const int max_nbors, const int maxspecial,
                           const double cell_size, const double gpu_split,
                           FILE *_screen, const void *pair_program,
                           const char *two, const char *three_center,
                           const char *three_end, const char *short_nbor,
                           const int onetype, const int onetype3,
                           const int spq, const int tpa_override) {
  screen=_screen;

  int gpu_nbor=0;
  if (device->gpu_mode()==Device<numtyp,acctyp>::GPU_NEIGH)
    gpu_nbor=1;
  else if (device->gpu_mode()==Device<numtyp,acctyp>::GPU_HYB_NEIGH)
    gpu_nbor=2;
  _gpu_nbor=gpu_nbor;

  int _gpu_host=0;
  int host_nlocal=hd_balancer.first_host_count(nlocal,gpu_split,gpu_nbor);
  if (host_nlocal>0)
    _gpu_host=1;

  // Allow forcing threads per atom to 1 for tersoff due to subg sync issue
  if (tpa_override)
    _threads_per_atom=tpa_override;
  else
    _threads_per_atom=device->threads_per_three();

  int success=device->init(*ans,false,false,nlocal,nall,maxspecial);
  if (success!=0)
    return success;

  if (ucl_device!=device->gpu) _compiled=false;

  ucl_device=device->gpu;
  atom=&device->atom;

  #ifdef THREE_CONCURRENT
  _end_command_queue=ucl_device->num_queues();
  ucl_device->push_command_queue();
  if (!ans2->init(ans->max_inum(),false,false,*(device->gpu)))
    return -3;
  ans2->cq(_end_command_queue);
  #endif

  _block_pair=device->pair_block_size();
  _block_size=device->block_ellipse();
  compile_kernels(*ucl_device,pair_program,two,three_center,three_end,
                  short_nbor,onetype,onetype3,spq);

  while (_threads_per_atom*_threads_per_atom>device->simd_size())
    _threads_per_atom = _threads_per_atom / 2;

  if (_threads_per_atom*_threads_per_atom>device->simd_size())
    return -10;

  success = device->init_nbor(nbor,nall,host_nlocal,nall,maxspecial,
                              _gpu_host,max_nbors,cell_size,true,1,true);
  if (success!=0)
    return success;

  // Initialize host-device load balancer
  hd_balancer.init(device,gpu_nbor,gpu_split);

  // Initialize timers for the selected GPU
  time_pair.init(*ucl_device);
  time_pair.zero();

  pos_tex.bind_float(atom->x,4);

  int ef_nall=nall;
  if (ef_nall==0)
    ef_nall=2000;

  _max_an_bytes=ans->gpu_bytes()+nbor->gpu_bytes();
  #ifdef THREE_CONCURRENT
  _max_an_bytes+=ans2->gpu_bytes();
  #endif

  return 0;
}

template <class numtyp, class acctyp>
void BaseThreeT::estimate_gpu_overhead(const int add_kernels) {
  device->estimate_gpu_overhead(4+add_kernels,_gpu_overhead,_driver_overhead);
}

template <class numtyp, class acctyp>
void BaseThreeT::clear_atomic() {
  // Output any timing information
  acc_timers();
  double avg_split=hd_balancer.all_avg_split();
  _gpu_overhead*=hd_balancer.timestep();
  _driver_overhead*=hd_balancer.timestep();
  device->output_times(time_pair,*ans,*nbor,avg_split,_max_bytes+_max_an_bytes,
                       _gpu_overhead,_driver_overhead,_threads_per_atom,screen);

  time_pair.clear();
  hd_balancer.clear();

  nbor->clear();
  ans->clear();
  #ifdef THREE_CONCURRENT
  ans2->clear();
  assert(ucl_device->num_queues()==_end_command_queue+1);
  // ucl_device will clean up the command queue in its destructor
//  ucl_device->pop_command_queue();
  #endif
}

// ---------------------------------------------------------------------------
// Copy neighbor list from host
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int * BaseThreeT::reset_nbors(const int nall, const int inum, const int nlist,
                              int *ilist, int *numj, int **firstneigh,
                              bool &success) {
  success=true;

  int mn=nbor->max_nbor_loop(nlist,numj,ilist);
  resize_atom(inum,nall,success);
  resize_local(nall,mn,success);
  if (!success)
    return nullptr;

  _nall = nall;

  // originally the requirement that nall == nlist was enforced
  // to allow direct indexing neighbors of neighbors after re-arrangement
//  nbor->get_host3(nall,nlist,ilist,numj,firstneigh,block_size());

  // now the requirement is removed, allowing to work within pair hybrid
  nbor->get_host(nlist,ilist,numj,firstneigh,block_size());
  nbor->copy_unpacked(nlist,mn);

  double bytes=ans->gpu_bytes()+nbor->gpu_bytes();
  #ifdef THREE_CONCURRENT
  bytes+=ans2->gpu_bytes();
  #endif
  if (bytes>_max_an_bytes)
    _max_an_bytes=bytes;

  return ilist;
}

// ---------------------------------------------------------------------------
// Build neighbor list on device
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
inline void BaseThreeT::build_nbor_list(const int inum, const int host_inum,
                                        const int nall, double **host_x,
                                        int *host_type, double *sublo,
                                        double *subhi, tagint *tag,
                                        int **nspecial, tagint **special,
                                        bool &success) {
  success=true;
  resize_atom(inum,nall,success);
  resize_local(nall,host_inum,nbor->max_nbors(),success);
  if (!success)
    return;
  atom->cast_copy_x(host_x,host_type);

  _nall = nall;

  // Increase the effective sub-domain size for neighbors of ghosts
  // This is still inefficient because we are calculating neighbors for more
  // ghosts than necessary due to increased ghost cutoff
  const double ncut=nbor->cutoff()*2.0;
  for (int i=0; i<3; i++) sublo[i]-=ncut;
  for (int i=0; i<3; i++) subhi[i]+=ncut;

  int mn;
  nbor->build_nbor_list(host_x, nall, host_inum, nall, *atom, sublo, subhi,
                        tag, nspecial, special, success, mn, ans->error_flag);
  nbor->copy_unpacked(nall,mn);

  double bytes=ans->gpu_bytes()+nbor->gpu_bytes();
  #ifdef THREE_CONCURRENT
  bytes+=ans2->gpu_bytes();
  #endif
  if (bytes>_max_an_bytes)
    _max_an_bytes=bytes;
}

// ---------------------------------------------------------------------------
// Copy nbor list from host if necessary and then calculate forces, virials,..
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void BaseThreeT::compute(const int f_ago, const int inum_full, const int nall,
                         const int nlist, double **host_x, int *host_type,
                         int *ilist, int *numj, int **firstneigh,
                         const bool eflag_in, const bool vflag_in,
                         const bool eatom, const bool vatom, int &host_start,
                         const double cpu_time, bool &success) {
  acc_timers();
  int eflag, vflag;
  if (eatom) eflag=2;
  else if (eflag_in) eflag=1;
  else eflag=0;
  if (vatom) vflag=2;
  else if (vflag_in) vflag=1;
  else vflag=0;

  #ifdef LAL_NO_BLOCK_REDUCE
  if (eflag) eflag=2;
  if (vflag) vflag=2;
  #endif

  set_kernel(eflag,vflag);
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
  #ifdef THREE_CONCURRENT
  ans2->inum(inum);
  #endif
  host_start=inum;

  if (ago==0) {
    reset_nbors(nall, inum, nlist, ilist, numj, firstneigh, success);
    if (!success)
      return;
  }

  atom->cast_x_data(host_x,host_type);
  hd_balancer.start_timer();
  atom->add_x_data(host_x,host_type);

  // _ainum to be used in loop() for short neighbor list build
  _ainum = nlist;

  int evatom=0;
  if (eatom || vatom)
    evatom=1;
  #ifdef THREE_CONCURRENT
  ucl_device->sync();
  #endif
  const int red_blocks=loop(eflag,vflag,evatom,success);
  ans->copy_answers(eflag_in,vflag_in,eatom,vatom,ilist,red_blocks);
  device->add_ans_object(ans);
  #ifdef THREE_CONCURRENT
  ans2->copy_answers(eflag_in,vflag_in,eatom,vatom,ilist,red_blocks);
  device->add_ans_object(ans2);
  #endif
  hd_balancer.stop_timer();
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU if necessary and then compute forces, virials, energies
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int ** BaseThreeT::compute(const int ago, const int inum_full, const int nall,
                           double **host_x, int *host_type, double *sublo,
                           double *subhi, tagint *tag, int **nspecial,
                           tagint **special, const bool eflag_in,
                           const bool vflag_in, const bool eatom,
                           const bool vatom, int &host_start,
                           int **ilist, int **jnum,
                           const double cpu_time, bool &success) {
  acc_timers();
  int eflag, vflag;
  if (eatom) eflag=2;
  else if (eflag_in) eflag=1;
  else eflag=0;
  if (vatom) vflag=2;
  else if (vflag_in) vflag=1;
  else vflag=0;

  #ifdef LAL_NO_BLOCK_REDUCE
  if (eflag) eflag=2;
  if (vflag) vflag=2;
  #endif

  set_kernel(eflag,vflag);
  if (inum_full==0) {
    host_start=0;
    // Make sure textures are correct if realloc by a different hybrid style
    resize_atom(0,nall,success);
    zero_timers();
    return nullptr;
  }

  hd_balancer.balance(cpu_time);
  int inum=hd_balancer.get_gpu_count(ago,inum_full);
  ans->inum(inum);
  #ifdef THREE_CONCURRENT
  ans2->inum(inum);
  #endif
  host_start=inum;

  // Build neighbor list on GPU if necessary
  if (ago==0) {
    build_nbor_list(inum, inum_full-inum, nall, host_x, host_type,
                    sublo, subhi, tag, nspecial, special, success);
    if (!success)
      return nullptr;
    hd_balancer.start_timer();
  } else {
    atom->cast_x_data(host_x,host_type);
    hd_balancer.start_timer();
    atom->add_x_data(host_x,host_type);
  }
  *ilist=nbor->host_ilist.begin();
  *jnum=nbor->host_acc.begin();

  // _ainum to be used in loop() for short neighbor list build
  _ainum = nall;

  int evatom=0;
  if (eatom || vatom)
    evatom=1;
  #ifdef THREE_CONCURRENT
  ucl_device->sync();
  #endif
  const int red_blocks=loop(eflag,vflag,evatom,success);
  ans->copy_answers(eflag_in,vflag_in,eatom,vatom,red_blocks);
  device->add_ans_object(ans);
  #ifdef THREE_CONCURRENT
  ans2->copy_answers(eflag_in,vflag_in,eatom,vatom,red_blocks);
  device->add_ans_object(ans2);
  #endif
  hd_balancer.stop_timer();

  return nbor->host_jlist.begin()-host_start;
}

template <class numtyp, class acctyp>
double BaseThreeT::host_memory_usage_atomic() const {
  return device->atom.host_memory_usage()+nbor->host_memory_usage()+
         4*sizeof(numtyp)+sizeof(BaseThree<numtyp,acctyp>);
}

template <class numtyp, class acctyp>
void BaseThreeT::compile_kernels(UCL_Device &dev, const void *pair_str,
                                 const char *two, const char *three_center,
                                 const char *three_end, const char* short_nbor,
                                 const int onetype, const int onetype3,
                                 const int spq) {
  if (_compiled && _onetype==onetype && _onetype3==onetype3 && _spq==spq)
    return;

  _onetype=onetype;
  _onetype3=onetype3;
  _spq=spq;

  std::string vatom_name=std::string(three_end)+"_vatom";
  if (pair_program) delete pair_program;
  pair_program=new UCL_Program(dev);
  std::string oclstring = device->compile_string()+" -DEVFLAG=1";
  if (_onetype>=0) oclstring+=" -DONETYPE="+device->toa(_onetype)+
                     " -DONETYPE3="+device->toa(_onetype3);
  if (_spq) oclstring+=" -DSPQ="+device->toa(_spq);
  pair_program->load_string(pair_str,oclstring.c_str(),nullptr,screen);
  k_three_center.set_function(*pair_program,three_center);
  k_three_end.set_function(*pair_program,three_end);
  k_three_end_vatom.set_function(*pair_program,vatom_name.c_str());
  k_pair.set_function(*pair_program,two);
  k_short_nbor.set_function(*pair_program,short_nbor);
  pos_tex.get_texture(*pair_program,"pos_tex");

  #if defined(LAL_OCL_EV_JIT)
  oclstring = device->compile_string()+" -DEVFLAG=0";
  if (_onetype>=0) oclstring+=" -DONETYPE="+device->toa(_onetype)+
                     " -DONETYPE3="+device->toa(_onetype3);
  if (_spq) oclstring+=" -DSPQ="+device->toa(_spq);
  if (pair_program_noev) delete pair_program_noev;
  pair_program_noev=new UCL_Program(dev);
  pair_program_noev->load_string(pair_str,oclstring.c_str(),nullptr,screen);
  k_three_center_noev.set_function(*pair_program_noev,three_center);
  k_three_end_noev.set_function(*pair_program_noev,three_end);
  k_pair_noev.set_function(*pair_program_noev,two);
  #else
  k_sel = &k_pair;
  k_3center_sel = &k_three_center;
  k_3end_sel = &k_three_end;
  #endif

  #ifdef THREE_CONCURRENT
  k_three_end.cq(ucl_device->cq(_end_command_queue));
  k_three_end_vatom.cq(ucl_device->cq(_end_command_queue));
  #if defined(LAL_OCL_EV_JIT)
  k_three_end_noev.cq(ucl_device->cq(_end_command_queue));
  #endif
  #endif

  _compiled=true;

  #if defined(USE_OPENCL) && (defined(CL_VERSION_2_1) || defined(CL_VERSION_3_0))
  if (dev.has_subgroup_support()) {
    size_t mx_subgroup_sz = k_pair.max_subgroup_size(_block_size);
    mx_subgroup_sz = std::min(mx_subgroup_sz, k_three_center.max_subgroup_size(_block_size));
    mx_subgroup_sz = std::min(mx_subgroup_sz, k_three_end.max_subgroup_size(_block_size));
    mx_subgroup_sz = std::min(mx_subgroup_sz, k_three_end_vatom.max_subgroup_size(_block_size));
    #if defined(LAL_OCL_EV_JIT)
    mx_subgroup_sz = std::min(mx_subgroup_sz, k_pair_noev.max_subgroup_size(_block_size));
    mx_subgroup_sz = std::min(mx_subgroup_sz, k_three_center_noev.max_subgroup_size(_block_size));
    mx_subgroup_sz = std::min(mx_subgroup_sz, k_three_end_noev.max_subgroup_size(_block_size));
    #endif
    if (_threads_per_atom > mx_subgroup_sz)
      _threads_per_atom = mx_subgroup_sz;
    device->set_simd_size(mx_subgroup_sz);
  }
  #endif

}

template class BaseThree<PRECISION,ACC_PRECISION>;
}
