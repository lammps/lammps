/***************************************************************************
                                lal_eam.cpp
                             -------------------
                      W. Michael Brown, Trung Dac Nguyen (ORNL)

  Class for acceleration of the eam pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : brownw@ornl.gov nguyentd@ornl.gov
 ***************************************************************************/
 
#ifdef USE_OPENCL
#include "eam_cl.h"
#else
#include "eam_ptx.h"
#endif

#include "lal_eam.h"
#include <cassert>
using namespace LAMMPS_AL;
#define EAMT EAM<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> global_device;

template <class numtyp, class acctyp>
EAMT::EAM() : _compiled(false), _max_bytes(0), _allocated(false) {
  device=&global_device;
  ans=new Answer<numtyp,acctyp>();
  nbor=new Neighbor();
}

template <class numtyp, class acctyp>
EAMT::~EAM() {
  delete ans;
  delete nbor;
  clear();
}
 
template <class numtyp, class acctyp>
int EAMT::bytes_per_atom(const int max_nbors) const {
  return device->atom.bytes_per_atom()+ans->bytes_per_atom()+
         nbor->bytes_per_atom(max_nbors);
}

template <class numtyp, class acctyp>
int EAMT::init_atomic(const int nlocal, const int nall,
                      const int max_nbors, const int maxspecial,
                      const double cell_size,
                      const double gpu_split, FILE *_screen,
                      const char *pair_program) {
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
    
  bool charge = true;
  bool rot = false;
  bool pre_cut = false;
  int success=device->init(*ans,charge,rot,nlocal,host_nlocal,nall,nbor,
                           maxspecial,_gpu_host,max_nbors,cell_size,pre_cut,
                           _threads_per_atom);
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

  _max_an_bytes=ans->gpu_bytes()
    +nbor->gpu_bytes();

  return success;
}

template <class numtyp, class acctyp>
int EAMT::init(const int ntypes, double host_cutforcesq,
              int **host_type2rhor, int **host_type2z2r, int *host_type2frho,
              double ***host_rhor_spline, double ***host_z2r_spline,
              double ***host_frho_spline,
              double rdr, double rdrho, int nrhor, int nrho, 
              int nz2r, int nfrho, int nr,
              const int nlocal, const int nall, const int max_nbors,
              const int maxspecial, const double cell_size,
              const double gpu_split, FILE *_screen) 
{
  int success;
  success=init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,gpu_split,
                            _screen,eam);
  
  if (success!=0)
    return success;
    
  // If atom type constants fit in shared memory use fast kernel
  int lj_types=ntypes;
  shared_types=false;

  int max_shared_types=this->device->max_shared_types();
  if (lj_types<=max_shared_types && this->_block_size>=max_shared_types) {
    lj_types=max_shared_types;
    shared_types=true;
  }
  
  _ntypes=lj_types;
  _cutforcesq=host_cutforcesq;
  _rdr=rdr;
  _rdrho = rdrho;
  _nrhor=nrhor;
  _nrho=nrho;
  _nz2r=nz2r;
  _nfrho=nfrho;
  _nr=nr;
  
  UCL_H_Vec<numtyp> dview_type(lj_types*lj_types*2,*(this->ucl_device),
                               UCL_WRITE_OPTIMIZED);
  
  for (int i=0; i<lj_types*lj_types*2; i++)
    dview_type[i]=(numtyp)0.0; 
                                
  // pack type2rhor and type2z2r
  type2rhor_z2r.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  
  this->atom->type_pack2(ntypes,lj_types,type2rhor_z2r,dview_type,
                        host_type2rhor,
                        host_type2z2r);
  
  // pack type2frho
  UCL_H_Vec<numtyp> dview_type2frho(ntypes,*(this->ucl_device),
                               UCL_WRITE_OPTIMIZED);

  type2frho.alloc(ntypes,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<ntypes; i++)
    dview_type2frho[i]=(numtyp)host_type2frho[i];
  ucl_copy(type2frho,dview_type2frho,false);
                        
  // pack frho_spline
  UCL_H_Vec<numtyp> dview_frho_spline(nfrho*(nr+1)*7,*(this->ucl_device),
                               UCL_WRITE_OPTIMIZED);
                               
  for (int ix=0; ix<nfrho; ix++)
    for (int iy=0; iy<nr+1; iy++)
      for (int iz=0; iz<7; iz++) 
    dview_frho_spline[ix*(nr+1)*7+iy*7+iz]=host_frho_spline[ix][iy][iz];
  
  frho_spline.alloc(nfrho*(nr+1)*7,*(this->ucl_device),UCL_READ_ONLY);
  ucl_copy(frho_spline,dview_frho_spline,false);
  
  // pack rhor_spline
  UCL_H_Vec<numtyp> dview_rhor_spline(nrhor*(nr+1)*7,*(this->ucl_device),
                               UCL_WRITE_OPTIMIZED);
                               
  for (int ix=0; ix<nrhor; ix++)
    for (int iy=0; iy<nr+1; iy++)
      for (int iz=0; iz<7; iz++) 
    dview_rhor_spline[ix*(nr+1)*7+iy*7+iz]=host_rhor_spline[ix][iy][iz];
  
  rhor_spline.alloc(nrhor*(nr+1)*7,*(this->ucl_device),UCL_READ_ONLY);
  ucl_copy(rhor_spline,dview_rhor_spline,false);
  
  // pack z2r_spline
  UCL_H_Vec<numtyp> dview_z2r_spline(nz2r*(nr+1)*7,*(this->ucl_device),
                               UCL_WRITE_OPTIMIZED);
                               
  for (int ix=0; ix<nz2r; ix++)
    for (int iy=0; iy<nr+1; iy++)
      for (int iz=0; iz<7; iz++) 
    dview_z2r_spline[ix*(nr+1)*7+iy*7+iz]=host_z2r_spline[ix][iy][iz];
  
  z2r_spline.alloc(nz2r*(nr+1)*7,*(this->ucl_device),UCL_READ_ONLY);
  ucl_copy(z2r_spline,dview_z2r_spline,false);

  _allocated=true;
  this->_max_bytes=type2rhor_z2r.row_bytes()
        + type2frho.row_bytes()
        + rhor_spline.row_bytes()+z2r_spline.row_bytes()
        + frho_spline.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void EAMT::estimate_gpu_overhead() {
  device->estimate_gpu_overhead(1,_gpu_overhead,_driver_overhead);
}

template <class numtyp, class acctyp>
void EAMT::clear() {
  if (!_allocated)
    return;
  _allocated=false;
  
  type2rhor_z2r.clear();
  type2frho.clear();
  rhor_spline.clear();
  z2r_spline.clear();
  frho_spline.clear();
  
  host_fp.clear();
  dev_fp.clear();

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
    k_energy.clear();
    delete pair_program;
    _compiled=false;
  }
  
  time_pair.clear();
  hd_balancer.clear();

  nbor->clear();
  ans->clear();
  device->clear();
}

template <class numtyp, class acctyp>
double EAMT::host_memory_usage() const {
  return device->atom.host_memory_usage()+nbor->host_memory_usage()+
         4*sizeof(numtyp)+sizeof(EAM<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Copy neighbor list from host
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int * EAMT::reset_nbors(const int nall, const int inum, int *ilist,
                                   int *numj, int **firstneigh, bool &success) {
  success=true;

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
inline void EAMT::build_nbor_list(const int inum, const int host_inum,
                                         const int nall, double **host_x,
                                         int *host_type, double *sublo,
                                         double *subhi, int *tag, 
                                         int **nspecial, int **special,
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
// Copy nbor list from host if necessary and then compute atom energies/forces
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void EAMT::compute_energy(const int f_ago, const int inum_full,
                   const int nall, double **host_x, int *host_type,
                   int *ilist, int *numj, int **firstneigh,
                   const bool eflag, const bool vflag,
                   const bool eatom, const bool vatom,
                   int &host_start, const double cpu_time,
                   bool &success, double *fp,
                   const int nlocal, double *boxlo, double *prd, 
                   double *evdwl) {
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
  hd_balancer.start_timer();
  atom->add_x_data(host_x,host_type);

  energy(eflag,vflag);
  
  // copy fp from device to host for comm
  
  ucl_copy(host_fp,dev_fp,false);
  
  acctyp *ap=host_fp.begin();
  for (int i=0; i<inum; i++) {
    fp[i]=*ap;
    ap++;
  }

  hd_balancer.stop_timer();
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU and then compute per-atom densities
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int** EAMT::compute_energy(const int ago, const int inum_full,
                    const int nall, double **host_x, int *host_type,
                    double *sublo, double *subhi, int *tag,
                    int **nspecial, int **special, const bool eflag, 
                    const bool vflag, const bool eatom,
                    const bool vatom, int &host_start,
                    int **ilist, int **jnum,
                    const double cpu_time, bool &success,
                    double *fp, double *boxlo, double *prd,
                    double *evdwl, int &inum) {
  acc_timers();
  if (inum_full==0) {
    host_start=0;
    // Make sure textures are correct if realloc by a different hybrid style
    resize_atom(0,nall,success);
    zero_timers();
    return NULL;
  }
  
  // load balance, returning the atom count on the device (inum)
  hd_balancer.balance(cpu_time);
  inum=hd_balancer.get_gpu_count(ago,inum_full);
  ans->inum(inum);
  host_start=inum;
 
  // Build neighbor list on GPU if necessary 
  if (ago==0) {
    build_nbor_list(inum, inum_full-inum, nall, host_x, host_type,
                    sublo, subhi, tag, nspecial, special, success);
    if (!success)
      return NULL;
    hd_balancer.start_timer();
  } else {
    atom->cast_x_data(host_x,host_type);
    hd_balancer.start_timer();
    atom->add_x_data(host_x,host_type);
  }
  *ilist=nbor->host_ilist.begin();
  *jnum=nbor->host_acc.begin();

  energy(eflag,vflag);
  
  // copy fp from device to host for comm
  
  ucl_copy(host_fp,dev_fp,false);
  
  acctyp *ap=host_fp.begin();
  for (int i=0; i<inum; i++) {
    fp[i]=*ap;
    ap++;
  }

  hd_balancer.stop_timer();
  
  return nbor->host_jlist.begin()-host_start;
}

// ---------------------------------------------------------------------------
// Copy nbor list from host if necessary and then calculate forces, virials,..
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void EAMT::compute(const int f_ago, const int inum_full,
                   const int nall, double **host_x, int *host_type,
                   int *ilist, int *numj, int **firstneigh,
                   const bool eflag, const bool vflag,
                   const bool eatom, const bool vatom,
                   int &host_start, const double cpu_time,
                   bool &success, double *host_q,
                   const int nlocal, double *boxlo, double *prd) {
  acc_timers();
  
  // compute density already took care of the neighbor list

  atom->cast_q_data(host_q);
  atom->add_q_data();
  
  loop(eflag,vflag);
  ans->copy_answers(eflag,vflag,eatom,vatom,ilist);
  device->add_ans_object(ans);
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU if necessary and then compute forces, virials, energies
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int** EAMT::compute(const int ago, const int inum_full,
                    const int nall, double **host_x, int *host_type,
                    double *sublo, double *subhi, int *tag,
                    int **nspecial, int **special, const bool eflag, 
                    const bool vflag, const bool eatom,
                    const bool vatom, int &host_start,
                    int **ilist, int **jnum,
                    const double cpu_time, bool &success,
                    double *host_q, double *boxlo, double *prd, int inum) {
  acc_timers();
  
  // use the atom count returned from load balance invoked in compute energy
  ans->inum(inum);
  host_start=inum;

  atom->cast_q_data(host_q);
  hd_balancer.start_timer();
  atom->add_q_data();
  *ilist=nbor->host_ilist.begin();
  *jnum=nbor->host_acc.begin();

  loop(eflag,vflag);
  ans->copy_answers(eflag,vflag,eatom,vatom);
  device->add_ans_object(ans);
  hd_balancer.stop_timer();
  
  return nbor->host_jlist.begin()-host_start;
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void EAMT::loop(const bool _eflag, const bool _vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int eflag, vflag;
  if (_eflag)
    eflag=1;
  else
    eflag=0;

  if (_vflag)
    vflag=1;
  else
    vflag=0;
  
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  
  if (shared_types) {
    this->k_pair_fast.set_size(GX,BX);
    this->k_pair_fast.run(&this->atom->dev_x.begin(), &this->atom->dev_q.begin(), 
                   &type2rhor_z2r.begin(),
                   &rhor_spline.begin(), &z2r_spline.begin(),
                   &this->nbor->dev_nbor.begin(),
                   &this->_nbor_data->begin(), &this->ans->dev_ans.begin(),
                   &this->ans->dev_engv.begin(), &eflag, &vflag, &ainum,
                   &nbor_pitch, &_cutforcesq, &_rdr, &_nr,
                   &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->dev_x.begin(), &this->atom->dev_q.begin(), 
                   &type2rhor_z2r.begin(),
                   &rhor_spline.begin(), &z2r_spline.begin(),
                   &this->nbor->dev_nbor.begin(),
                   &this->_nbor_data->begin(), &this->ans->dev_ans.begin(),
                   &this->ans->dev_engv.begin(), &eflag, &vflag, &ainum,
                   &nbor_pitch, &_ntypes, &_cutforcesq, &_rdr, &_nr,
                   &this->_threads_per_atom);
  }

  this->time_pair.stop();
}

// ---------------------------------------------------------------------------
// Calculate per-atom energies and forces
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void EAMT::energy(const bool _eflag, const bool _vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int eflag, vflag;
  if (_eflag)
    eflag=1;
  else
    eflag=0;

  if (_vflag)
    vflag=1;
  else
    vflag=0;
  
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  
  
  this->k_energy.set_size(GX,BX);
  this->k_energy.run(&this->atom->dev_x.begin(), 
                 &type2rhor_z2r.begin(), &type2frho.begin(),
                 &rhor_spline.begin(), &frho_spline.begin(),
                 &this->nbor->dev_nbor.begin(), &this->_nbor_data->begin(),
                 &dev_fp.begin(), 
                 &ans->dev_engv.begin(),
                 &eflag, &vflag, &ainum,
                 &nbor_pitch, 
                 &_ntypes, &_cutforcesq, 
                 &_rdr, &_rdrho,
                 &_nrho, &_nr,
                 &this->_threads_per_atom);

  this->time_pair.stop();
}

template <class numtyp, class acctyp>
void EAMT::compile_kernels(UCL_Device &dev, const char *pair_str) {
  if (_compiled)
    return;

  std::string flags="-cl-fast-relaxed-math -cl-mad-enable "+
                    std::string(OCL_PRECISION_COMPILE)+" -D"+
                    std::string(OCL_VENDOR);

  pair_program=new UCL_Program(dev);
  pair_program->load_string(pair_str,flags.c_str());
  k_pair_fast.set_function(*pair_program,"kernel_pair_fast");
  k_pair.set_function(*pair_program,"kernel_pair");
  k_energy.set_function(*pair_program,"kernel_energy");
  pos_tex.get_texture(*pair_program,"pos_tex");
  q_tex.get_texture(*pair_program,"q_tex");

  _compiled=true;
}

template class EAM<PRECISION,ACC_PRECISION>;
