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

#ifdef USE_OPENCL
#include "pppm_gpu_cl.h"
#else
#include "pppm_gpu_ptx.h"
#endif
#include "pppm_gpu_memory.h"

#define PPPMGPUMemoryT PPPMGPUMemory<numtyp, acctyp>

extern PairGPUDevice<PRECISION,ACC_PRECISION> pair_gpu_device;

template <class numtyp, class acctyp>
PPPMGPUMemoryT::PPPMGPUMemory() : _allocated(false), _compiled(false),
                                  _max_bytes(0) {
  device=&pair_gpu_device;
  ans=new PairGPUAns<numtyp,acctyp>();
}

template <class numtyp, class acctyp>
PPPMGPUMemoryT::~PPPMGPUMemory() {
  clear();
  delete ans;
}

template <class numtyp, class acctyp>
int PPPMGPUMemoryT::bytes_per_atom() const {
  return device->atom.bytes_per_atom()+ans->bytes_per_atom();
}

template <class numtyp, class acctyp>
bool PPPMGPUMemoryT::init(const int nlocal, const int nall, FILE *_screen) {
  screen=_screen;

  if (!device->init(*ans,true,false,nlocal,nall))
    return false;
  ucl_device=device->gpu;
  atom=&device->atom;

  _block_size=BLOCK_1D;
  if (static_cast<size_t>(_block_size)>ucl_device->group_size())
    _block_size=ucl_device->group_size();
  compile_kernels(*ucl_device);

  // Initialize timers for the selected GPU
  time_in.init(*ucl_device);
  time_in.zero();

  pos_tex.bind_float(atom->dev_x,4);
  q_tex.bind_float(atom->dev_q,1);

  _allocated=true;
  _max_bytes=0;
  _max_an_bytes=ans->gpu_bytes();

  return true;
}

template <class numtyp, class acctyp>
void PPPMGPUMemoryT::clear() {
  if (!_allocated)
    return;
  _allocated=false;
  
  acc_timers();
  device->output_kspace_times(time_in,*ans,_max_bytes+_max_an_bytes,
                              screen);

  if (_compiled) {
    k_compute.clear();
    delete pppm_program;
    _compiled=false;
  }

  time_in.clear();

  device->clear();
}

// ---------------------------------------------------------------------------
// Copy nbor list from host if necessary and then calculate forces, virials,..
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void PPPMGPUMemoryT::compute(const int ago, const int nlocal, const int nall,
                             double **host_x, int *host_type,
                             bool &success, double *host_q) {
  acc_timers();
  if (nlocal==0) {
    zero_timers();
    return;
  }
  
  ans->inum(nlocal);

  if (ago==0) {
    resize_atom(nlocal,nall,success);
    resize_local(nlocal,success);
    if (!success)
      return;

    double bytes=ans->gpu_bytes();
    if (bytes>_max_an_bytes)
      _max_an_bytes=bytes;
  }

  atom->cast_x_data(host_x,host_type);
  atom->cast_q_data(host_q);
  atom->add_x_data(host_x,host_type);
  atom->add_q_data();


  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/BX));

  int ainum=this->ans->inum();
  int anall=this->atom->nall();

//  this->time_pair.start();
//  this->k_pair.set_size(GX,BX);
//  this->k_pair.run(&this->atom->dev_x.begin(), &lj1.begin(), &lj3.begin(),
//                   &_lj_types, &sp_lj.begin(), &this->nbor->dev_nbor.begin(),
//                   &this->ans->dev_ans.begin(),
//                   &this->ans->dev_engv.begin(), &eflag, &vflag, &ainum,
//                   &anall, &nbor_pitch, &this->atom->dev_q.begin(),
//                   &cutsq.begin(), &_qqrd2e);
//  this->time_pair.stop();


//  ans->copy_answers(eflag,vflag,eatom,vatom,ilist);
//  device->add_ans_object(ans);
}

template <class numtyp, class acctyp>
double PPPMGPUMemoryT::host_memory_usage() const {
  return device->atom.host_memory_usage()+sizeof(PPPMGPUMemory<numtyp,acctyp>);
}

template <class numtyp, class acctyp>
void PPPMGPUMemoryT::compile_kernels(UCL_Device &dev) {
  if (_compiled)
    return;

  std::string flags="-cl-fast-relaxed-math -cl-mad-enable "+
                    std::string(OCL_PRECISION_COMPILE);

  pppm_program=new UCL_Program(dev);
  pppm_program->load_string(pppm_gpu_kernel,flags.c_str());
//  k_pair_fast.set_function(*pair_program,"kernel_pair_fast");
//  k_pair.set_function(*pair_program,"kernel_pair");
  pos_tex.get_texture(*pppm_program,"pos_tex");
  q_tex.get_texture(*pppm_program,"q_tex");

  _compiled=true;
}

template class PPPMGPUMemory<PRECISION,ACC_PRECISION>;

