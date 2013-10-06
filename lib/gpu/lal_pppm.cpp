/***************************************************************************
                                  pppm.cpp
                             -------------------
                            W. Michael Brown (ORNL)

  Class for PPPM acceleration

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : brownw@ornl.gov
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "pppm_cl.h"
#elif defined(USE_CUDART)
const char *pppm_f=0;
const char *pppm_d=0;
#else
#include "pppm_f_cubin.h"
#include "pppm_d_cubin.h"
#endif
#include "lal_pppm.h"
#include <cassert>

using namespace LAMMPS_AL;
#define PPPMT PPPM<numtyp, acctyp, grdtyp, grdtyp4>

extern Device<PRECISION,ACC_PRECISION> global_device;

template <class numtyp, class acctyp, class grdtyp, class grdtyp4>
PPPMT::PPPM() : _allocated(false), _compiled(false),
                                  _max_bytes(0) {
  device=&global_device;
  ans=new Answer<numtyp,acctyp>();
}

template <class numtyp, class acctyp, class grdtyp, class grdtyp4>
PPPMT::~PPPM() {
  clear(0.0);
  delete ans;
}

template <class numtyp, class acctyp, class grdtyp, class grdtyp4>
int PPPMT::bytes_per_atom() const {
  return device->atom.bytes_per_atom()+ans->bytes_per_atom()+1;
}

template <class numtyp, class acctyp, class grdtyp, class grdtyp4>
grdtyp * PPPMT::init(const int nlocal, const int nall, FILE *_screen,
                              const int order, const int nxlo_out,
                              const int nylo_out, const int nzlo_out,
                              const int nxhi_out, const int nyhi_out,
                              const int nzhi_out, grdtyp **rho_coeff,
                              grdtyp **vd_brick_p, const double slab_volfactor, 
                              const int nx_pppm, const int ny_pppm,
                              const int nz_pppm, const bool split, int &flag) {
  _max_bytes=10;
  screen=_screen;
  _kspace_split=split;
  bool success=true;

  flag=device->init(*ans,nlocal,nall);
  if (flag!=0)
    return 0;
  if (sizeof(grdtyp)==sizeof(double) && device->double_precision()==false) {
    flag=-5;
    return 0;
  }
  if (device->ptx_arch()>0.0 && device->ptx_arch()<1.1) {
    flag=-4;
    return 0;
  }

  ucl_device=device->gpu;
  atom=&device->atom;

  _block_size=device->pppm_block();
  _pencil_size=device->num_mem_threads();
  _block_pencils=_block_size/_pencil_size;

  compile_kernels(*ucl_device);

  // Initialize timers for the selected GPU
  time_in.init(*ucl_device);
  time_in.zero();
  time_out.init(*ucl_device);
  time_out.zero();
  time_map.init(*ucl_device);
  time_map.zero();
  time_rho.init(*ucl_device);
  time_rho.zero();
  time_interp.init(*ucl_device);
  time_interp.zero();

  pos_tex.bind_float(atom->x,4);
  q_tex.bind_float(atom->q,1);

  _allocated=true;
  _max_bytes=0;
  _max_an_bytes=ans->gpu_bytes();
  
  _order=order;
  _order_m_1=order-1;
  _order2=_order_m_1*_order;
  _nlower=-(_order-1)/2;
  _nupper=order/2;
  _nxlo_out=nxlo_out;
  _nylo_out=nylo_out;
  _nzlo_out=nzlo_out;
  _nxhi_out=nxhi_out;
  _nyhi_out=nyhi_out;
  _nzhi_out=nzhi_out;

  _slab_volfactor=slab_volfactor;
  _nx_pppm=nx_pppm;
  _ny_pppm=ny_pppm;
  _nz_pppm=nz_pppm;

  _max_brick_atoms=10;

  // Get rho_coeff on device
  int n2lo=(1-order)/2;
  int numel=order*( order/2 - n2lo + 1 );
  success=success && (d_rho_coeff.alloc(numel,*ucl_device,UCL_READ_ONLY)==
                      UCL_SUCCESS);
  UCL_H_Vec<grdtyp> view;
  view.view(rho_coeff[0]+n2lo,numel,*ucl_device);
  ucl_copy(d_rho_coeff,view,true);
  _max_bytes+=d_rho_coeff.row_bytes();
  
  // Allocate storage for grid
  _npts_x=nxhi_out-nxlo_out+1;
  _npts_y=nyhi_out-nylo_out+1;
  _npts_z=nzhi_out-nzlo_out+1;
  _npts_yx=_npts_x*_npts_y;
  success=success && (brick.alloc(_npts_x*_npts_y*_npts_z,*ucl_device,
                                  UCL_READ_ONLY,UCL_WRITE_ONLY)==UCL_SUCCESS);
  success=success && (vd_brick.alloc(_npts_x*_npts_y*_npts_z*4,*ucl_device,
                                    UCL_READ_WRITE,UCL_READ_ONLY)==UCL_SUCCESS);
  *vd_brick_p=vd_brick.host.begin();
  _max_bytes+=brick.device.row_bytes()+vd_brick.device.row_bytes();

  // Allocate vector with count of atoms assigned to each grid point
  _nlocal_x=_npts_x+_nlower-_nupper;
  _nlocal_y=_npts_y+_nlower-_nupper;
  _nlocal_z=_npts_z+_nlower-_nupper;
  _nlocal_yx=_nlocal_x*_nlocal_y;
  _atom_stride=_nlocal_x*_nlocal_y*_nlocal_z;
  success=success && (d_brick_counts.alloc(_atom_stride,*ucl_device)==
                      UCL_SUCCESS);
  _max_bytes+=d_brick_counts.row_bytes();

  // Allocate storage for atoms assigned to each grid point
  success=success && (d_brick_atoms.alloc(_atom_stride*_max_brick_atoms,
                                          *ucl_device)==UCL_SUCCESS);
  _max_bytes+=d_brick_atoms.row_bytes();

  // Allocate error flags for checking out of bounds atoms
  success=success && (error_flag.alloc(1,*ucl_device,UCL_READ_ONLY,
                                       UCL_READ_WRITE)==UCL_SUCCESS);
  if (!success) {
    flag=-3;
    return 0;
  }
  
  error_flag.device.zero();
  _max_bytes+=1;
  
  _cpu_idle_time=0.0;

  return brick.host.begin();
}

template <class numtyp, class acctyp, class grdtyp, class grdtyp4>
void PPPMT::clear(const double cpu_time) {
  if (!_allocated)
    return;
  _allocated=false;
  _precompute_done=false;
  
  brick.clear();
  vd_brick.clear();
  d_brick_counts.clear();
  error_flag.clear();
  d_brick_atoms.clear();
  
  acc_timers();
  device->output_kspace_times(time_in,time_out,time_map,time_rho,time_interp,
                              *ans,_max_bytes+_max_an_bytes,cpu_time,
                              _cpu_idle_time,screen);

  if (_compiled) {
    k_particle_map.clear();
    k_make_rho.clear();
    k_interp.clear();
    delete pppm_program;
    _compiled=false;
  }

  time_in.clear();
  time_out.clear();
  time_map.clear();
  time_rho.clear();
  time_interp.clear();

  ans->clear();
  device->clear();
}

// ---------------------------------------------------------------------------
// Charge assignment that can be performed asynchronously
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp, class grdtyp, class grdtyp4>
void PPPMT::_precompute(const int ago, const int nlocal, const int nall,
                                 double **host_x, int *host_type, bool &success,
                                 double *host_q, double *boxlo, 
                                 const double delxinv, const double delyinv,
                                 const double delzinv) {
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

  time_map.start();

  // Compute the block size and grid size to keep all cores busy
  int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/BX));

  int ainum=this->ans->inum();
  
  // Boxlo adjusted to be upper left brick and shift for even spline order
  double shift=0.0;
  if (_order % 2)
    shift=0.5;
  _brick_x=boxlo[0]+(_nxlo_out-_nlower-shift)/delxinv;
  _brick_y=boxlo[1]+(_nylo_out-_nlower-shift)/delyinv;
  _brick_z=boxlo[2]+(_nzlo_out-_nlower-shift)/delzinv;
  
  _delxinv=delxinv;
  _delyinv=delyinv;
  _delzinv=delzinv;
  double delvolinv = delxinv*delyinv*delzinv;
  grdtyp f_delvolinv = delvolinv;

  device->zero(d_brick_counts,d_brick_counts.numel());
  k_particle_map.set_size(GX,BX);
  k_particle_map.run(&atom->x, &atom->q, &f_delvolinv, &ainum,
                     &d_brick_counts, &d_brick_atoms, &_brick_x, &_brick_y, 
                     &_brick_z, &_delxinv, &_delyinv, &_delzinv, &_nlocal_x,
                     &_nlocal_y, &_nlocal_z, &_atom_stride, &_max_brick_atoms,
                     &error_flag);
  time_map.stop();

  time_rho.start();
  BX=block_size();

  GX=static_cast<int>(ceil(static_cast<double>(_npts_y*_npts_z)/
                      _block_pencils));
  k_make_rho.set_size(GX,BX);
  k_make_rho.run(&d_brick_counts, &d_brick_atoms, &brick, &d_rho_coeff,
                 &_atom_stride, &_npts_x, &_npts_y, &_npts_z, &_nlocal_x,
                 &_nlocal_y, &_nlocal_z, &_order_m_1, &_order, &_order2);
  time_rho.stop();

  time_out.start();
  brick.update_host(_npts_yx*_npts_z,true);
  error_flag.update_host(true);
  time_out.stop();

  _precompute_done=true;
}

// ---------------------------------------------------------------------------
// Charge spreading stuff
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp, class grdtyp, class grdtyp4>
int PPPMT::spread(const int ago, const int nlocal, const int nall,
                           double **host_x, int *host_type, bool &success,
                           double *host_q, double *boxlo, 
                           const double delxinv, const double delyinv,
                           const double delzinv) {
  if (_precompute_done==false) {
    atom->acc_timers();
    _precompute(ago,nlocal,nall,host_x,host_type,success,host_q,boxlo,delxinv,
                delyinv,delzinv);
  }

  device->stop_host_timer();
  
  if (!success || nlocal==0)
    return 0;
    
  double t=MPI_Wtime();
  time_out.sync_stop();
  _cpu_idle_time+=MPI_Wtime()-t;

  _precompute_done=false;

  if (error_flag[0]==2) {
    // Not enough storage for atoms on the brick
    _max_brick_atoms*=2;
    error_flag.device.zero();
    d_brick_atoms.resize(_atom_stride*_max_brick_atoms);
    _max_bytes+=d_brick_atoms.row_bytes();
    return spread(ago,nlocal,nall,host_x,host_type,success,host_q,boxlo, 
                  delxinv,delyinv,delzinv);
  }
  
  return error_flag[0];
}

// ---------------------------------------------------------------------------
// Charge spreading stuff
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp, class grdtyp, class grdtyp4>
void PPPMT::interp(const grdtyp qqrd2e_scale) {
  time_in.start();
  vd_brick.update_device(true);
  time_in.stop();
  
  time_interp.start();
  // Compute the block size and grid size to keep all cores busy
  int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/BX));

  int ainum=this->ans->inum();
  
  k_interp.set_size(GX,BX);
  k_interp.run(&atom->x, &atom->q, &ainum, &vd_brick, &d_rho_coeff,
               &_npts_x, &_npts_yx, &_brick_x, &_brick_y, &_brick_z, &_delxinv,
               &_delyinv, &_delzinv, &_order, &_order2, &qqrd2e_scale, 
               &ans->force);
  time_interp.stop();

  ans->copy_answers(false,false,false,false);
  if (_kspace_split==false)
    device->add_ans_object(ans);
}

template <class numtyp, class acctyp, class grdtyp, class grdtyp4>
double PPPMT::host_memory_usage() const {
  return device->atom.host_memory_usage()+
         sizeof(PPPM<numtyp,acctyp,grdtyp,grdtyp4>);
}

template <class numtyp, class acctyp, class grdtyp, class grdtyp4>
void PPPMT::compile_kernels(UCL_Device &dev) {
  if (_compiled)
    return;

  if (sizeof(grdtyp)==sizeof(double) && ucl_device->double_precision()==false)
    return;

  std::string flags=device->compile_string();
  #ifdef USE_OPENCL
  flags+=std::string(" -Dgrdtyp=")+ucl_template_name<grdtyp>()+" -Dgrdtyp4="+
         ucl_template_name<grdtyp>()+"4";
  #endif

  pppm_program=new UCL_Program(dev);
  
  #ifdef USE_OPENCL
  pppm_program->load_string(pppm,flags.c_str());
  #else
  if (sizeof(grdtyp)==sizeof(float))
    pppm_program->load_string(pppm_f,flags.c_str());
  else
    pppm_program->load_string(pppm_d,flags.c_str());
  #endif

  k_particle_map.set_function(*pppm_program,"particle_map");
  k_make_rho.set_function(*pppm_program,"make_rho");
  k_interp.set_function(*pppm_program,"interp");
  pos_tex.get_texture(*pppm_program,"pos_tex");
  q_tex.get_texture(*pppm_program,"q_tex");

  _compiled=true;
}

template class PPPM<PRECISION,ACC_PRECISION,float,_lgpu_float4>;
template class PPPM<PRECISION,ACC_PRECISION,double,_lgpu_double4>;
