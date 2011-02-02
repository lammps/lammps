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
bool PPPMGPUMemoryT::init(const int nlocal, const int nall, FILE *_screen,
                          const int order, const int nxlo_out,
                          const int nylo_out, const int nzlo_out,
                          const int nxhi_out, const int nyhi_out,
                          const int nzhi_out, double **rho_coeff) {
  clear();
  
  _max_bytes=0;
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
  time_kernel.init(*ucl_device);
  time_kernel.zero();

  pos_tex.bind_float(atom->dev_x,4);
  q_tex.bind_float(atom->dev_q,1);

  _allocated=true;
  _max_bytes=0;
  _max_an_bytes=ans->gpu_bytes();
  
  _order=order;
  _nlower=-(_order-1)/2;
  _nxlo_out=nxlo_out;
  _nylo_out=nylo_out;
  _nzlo_out=nzlo_out;
  _nxhi_out=nxhi_out;
  _nyhi_out=nyhi_out;
  _nzhi_out=nzhi_out;
  
  // Get rho_coeff on device
  int n2lo=(1-order)/2;
  int numel=order*( order/2 - n2lo + 1 );
  d_rho_coeff.alloc(numel,*ucl_device,UCL_READ_ONLY);
  UCL_H_Vec<double> view;
  view.view(rho_coeff[0]+n2lo,numel,*ucl_device);
  ucl_copy(d_rho_coeff,view,true);
  _max_bytes+=d_rho_coeff.row_bytes();
  
  // Allocate vector with count of atoms assigned to each grid point
  _npts_x=nxhi_out-nxlo_out+1;
  _npts_y=nyhi_out-nylo_out+1;
  _npts_z=nzhi_out-nzlo_out+1;
  numel=_npts_x*_npts_y*_npts_z;
  d_brick_counts.alloc(numel,*ucl_device);
  _max_bytes+=d_brick_counts.row_bytes();

  // Allocate error flags for checking out of bounds atoms
  h_error_flag.alloc(1,*ucl_device);
  d_error_flag.alloc(1,*ucl_device,UCL_WRITE_ONLY);
  d_error_flag.zero();
  
  return true;
}

template <class numtyp, class acctyp>
void PPPMGPUMemoryT::clear() {
  if (!_allocated)
    return;
  _allocated=false;
  
  d_brick_counts.clear();
  h_error_flag.clear();
  d_error_flag.clear();
  acc_timers();
  device->output_kspace_times(time_in,time_kernel,*ans,_max_bytes+_max_an_bytes,
                              screen);

  if (_compiled) {
    k_particle_map.clear();
    delete pppm_program;
    _compiled=false;
  }

  time_in.clear();
  time_kernel.clear();

  device->clear();
}

// ---------------------------------------------------------------------------
// Copy nbor list from host if necessary and then calculate forces, virials,..
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int PPPMGPUMemoryT::compute(const int ago, const int nlocal, const int nall,
                            double **host_x, int *host_type, bool &success,
                            double *host_q, double *boxlo, 
                            const double delxinv, const double delyinv,
                            const double delzinv) {
  acc_timers();
  if (nlocal==0) {
    zero_timers();
    return 0;
  }
  
  ans->inum(nlocal);

  if (ago==0) {
    resize_atom(nlocal,nall,success);
    resize_local(nlocal,success);
    if (!success)
      return 0;

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

  int _max_atoms=10;
  int ainum=this->ans->inum();
  
  // Boxlo adjusted to include ghost cells and shift for even stencil order
  double lo_shift_x=static_cast<double>(_nlower);
  double lo_shift_y=static_cast<double>(_nlower);
  double lo_shift_z=static_cast<double>(_nlower);
  if (_order % 2) {
    lo_shift_x-=0.5;
    lo_shift_y-=0.5;
    lo_shift_z-=0.5;
  }
  lo_shift_x/=delxinv;
  lo_shift_y/=delyinv;
  lo_shift_z/=delzinv;
  
  numtyp f_boxlo_x=boxlo[0]+lo_shift_x;
  numtyp f_boxlo_y=boxlo[1]+lo_shift_y;
  numtyp f_boxlo_z=boxlo[2]+lo_shift_z;
  numtyp f_delxinv=delxinv;
  numtyp f_delyinv=delyinv;
  numtyp f_delzinv=delzinv;

  time_kernel.start();
  d_brick_counts.zero();
  k_particle_map.set_size(GX,BX);
  k_particle_map.run(&atom->dev_x.begin(), &ainum, &d_brick_counts.begin(),
                     &d_brick_counts.begin(), &f_boxlo_x, &f_boxlo_y, 
                     &f_boxlo_z, &f_delxinv, &f_delyinv, &f_delzinv, &_npts_x,
                     &_npts_y, &_npts_z, &_max_atoms, &d_error_flag.begin());
  time_kernel.stop();
  ucl_copy(h_error_flag,d_error_flag,false);
  
  if (h_error_flag[0]==2)
    std::cerr << "NEED TO RESIZE!\n";
  return h_error_flag[0];
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
  k_particle_map.set_function(*pppm_program,"particle_map");
  pos_tex.get_texture(*pppm_program,"pos_tex");
  q_tex.get_texture(*pppm_program,"q_tex");

  _compiled=true;
}

template class PPPMGPUMemory<PRECISION,ACC_PRECISION>;

