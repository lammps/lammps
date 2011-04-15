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

#ifdef USE_OPENCL
#include "gb_gpu_cl.h"
#include "gb_gpu_nbor_cl.h"
#else
#include "gb_gpu_ptx.h"
#endif

#include "gb_gpu_memory.h"
#include <cassert>
#define GB_GPU_MemoryT GB_GPU_Memory<numtyp, acctyp>

extern PairGPUDevice<PRECISION,ACC_PRECISION> pair_gpu_device;

template <class numtyp, class acctyp>
GB_GPU_MemoryT::GB_GPU_Memory() : _allocated(false), _compiled(false),
                                  _max_bytes(0.0) {
  device=&pair_gpu_device;
  ans=new PairGPUAns<numtyp,acctyp>();
  nbor=new PairGPUNbor;
}

template <class numtyp, class acctyp>
GB_GPU_MemoryT::~GB_GPU_Memory() { 
  clear();
  delete ans;
  delete nbor;
}
 
template <class numtyp, class acctyp>
int GB_GPU_MemoryT::bytes_per_atom(const int max_nbors) const {
  return device->atom.bytes_per_atom()+ans->bytes_per_atom()+
         nbor->bytes_per_atom(max_nbors);
}

template <class numtyp, class acctyp>
int GB_GPU_MemoryT::init(const int ntypes, const double gamma, 
                         const double upsilon, const double mu, 
                         double **host_shape, double **host_well, 
                         double **host_cutsq, double **host_sigma, 
                         double **host_epsilon, double *host_lshape, 
                         int **h_form, double **host_lj1, double **host_lj2,
                         double **host_lj3, double **host_lj4,
                         double **host_offset, const double *host_special_lj,
                         const int nlocal, const int nall,
                         const int max_nbors, const double cell_size,
                         const double gpu_split, FILE *_screen) {
  nbor_time_avail=false;
  screen=_screen;

  bool gpu_nbor=false;
  if (device->gpu_mode()==PairGPUDevice<numtyp,acctyp>::GPU_NEIGH)
    gpu_nbor=true;

  int _gpu_host=0;
  int host_nlocal=hd_balancer.first_host_count(nlocal,gpu_split,gpu_nbor);
  if (host_nlocal>0)
    _gpu_host=1;
  
  _threads_per_atom=device->threads_per_atom();
  int success=device->init(*ans,false,true,nlocal,host_nlocal,nall,nbor,0,
                           _gpu_host,max_nbors,cell_size,true);
  if (success!=0)
    return success;
    
  ucl_device=device->gpu;
  atom=&device->atom;

  _block_size=device->pair_block_size();
  compile_kernels(*ucl_device);

  // Initialize host-device load balancer
  hd_balancer.init(device,gpu_nbor,gpu_split);

  // Initialize timers for the selected GPU
  time_pair.init(*ucl_device);
  time_pair.zero();
  
  // If atom type constants fit in shared memory use fast kernel
  int lj_types=ntypes;
  shared_types=false;
  int max_shared_types=device->max_shared_types();
  if (lj_types<=max_shared_types && _block_size>=max_shared_types) {
    lj_types=max_shared_types;
    shared_types=true;
  }
  _lj_types=lj_types;

  // Allocate a host write buffer for copying type data
  UCL_H_Vec<numtyp> host_write(lj_types*lj_types*32,*ucl_device,
                               UCL_WRITE_OPTIMIZED);

  for (int i=0; i<lj_types*lj_types; i++)
    host_write[i]=0.0;

  sigma_epsilon.alloc(lj_types*lj_types,*ucl_device,UCL_READ_ONLY);
  this->atom->type_pack2(ntypes,lj_types,sigma_epsilon,host_write,
			 host_sigma,host_epsilon);

  cut_form.alloc(lj_types*lj_types,*ucl_device,UCL_READ_ONLY);
  this->atom->type_pack2(ntypes,lj_types,cut_form,host_write,
			 host_cutsq,h_form);

  lj1.alloc(lj_types*lj_types,*ucl_device,UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,lj1,host_write,host_lj1,host_lj2,
			 host_cutsq,h_form);

  lj3.alloc(lj_types*lj_types,*ucl_device,UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,lj3,host_write,host_lj3,host_lj4,
		         host_offset);

  dev_error.alloc(1,*ucl_device);
  dev_error.zero();
    
  _allocated=true;
    
  host_form=h_form;
    
  // Initialize timers for the selected GPU
  time_kernel.init(*ucl_device);
  time_gayberne.init(*ucl_device);
  time_kernel2.init(*ucl_device);
  time_gayberne2.init(*ucl_device);
  time_kernel.zero();
  time_gayberne.zero();
  time_kernel2.zero();
  time_gayberne2.zero();
    
  // Allocate, cast and asynchronous memcpy of constant data
  // Copy data for bonded interactions
  gamma_upsilon_mu.alloc(7,*ucl_device,UCL_READ_ONLY);
  host_write[0]=static_cast<numtyp>(gamma); 
  host_write[1]=static_cast<numtyp>(upsilon);
  host_write[2]=static_cast<numtyp>(mu);
  host_write[3]=static_cast<numtyp>(host_special_lj[0]);
  host_write[4]=static_cast<numtyp>(host_special_lj[1]);
  host_write[5]=static_cast<numtyp>(host_special_lj[2]);
  host_write[6]=static_cast<numtyp>(host_special_lj[3]);
  ucl_copy(gamma_upsilon_mu,host_write,7,false);

  lshape.alloc(ntypes,*ucl_device,UCL_READ_ONLY);
  UCL_H_Vec<double> d_view;
  d_view.view(host_lshape,lshape.numel(),*ucl_device);
  ucl_copy(lshape,d_view,false);
    
  // Copy shape, well, sigma, epsilon, and cutsq onto GPU
  // - cast if necessary
  shape.alloc(ntypes,*ucl_device,UCL_READ_ONLY);
  for (int i=0; i<ntypes; i++) {
    host_write[i*4]=host_shape[i][0];
    host_write[i*4+1]=host_shape[i][1];
    host_write[i*4+2]=host_shape[i][2];
  }
  UCL_H_Vec<numtyp4> view4;
  view4.view((numtyp4*)host_write.begin(),shape.numel(),*ucl_device);
  ucl_copy(shape,view4,false);

  well.alloc(ntypes,*ucl_device,UCL_READ_ONLY);
  for (int i=0; i<ntypes; i++) {
    host_write[i*4]=host_well[i][0];
    host_write[i*4+1]=host_well[i][1];
    host_write[i*4+2]=host_well[i][2];
  }
  view4.view((numtyp4*)host_write.begin(),well.numel(),*ucl_device);
  ucl_copy(well,view4,false);
  
  // See if we want fast GB-sphere or sphere-sphere calculations
  multiple_forms=false;
  for (int i=1; i<ntypes; i++)
    for (int j=i; j<ntypes; j++) 
      if (host_form[i][j]!=ELLIPSE_ELLIPSE)
        multiple_forms=true;
  if (multiple_forms && host_nlocal>0) {
    std::cerr << "Cannot use Gayberne with multiple forms and GPU neighbor.\n";
    exit(1);
  }
  
  if (multiple_forms)
    ans->dev_ans.zero();

  _max_bytes=ans->gpu_bytes()+nbor->gpu_bytes();

  // Memory for ilist ordered by particle type
  if (host_olist.alloc(nbor->max_atoms(),*ucl_device)==UCL_SUCCESS)
    return 0;
  else return -3;
}

template <class numtyp, class acctyp>
void GB_GPU_MemoryT::estimate_gpu_overhead() {
  device->estimate_gpu_overhead(2,_gpu_overhead,_driver_overhead);
}

template <class numtyp, class acctyp>
void GB_GPU_MemoryT::clear() {
  if (!_allocated)
    return;

  UCL_H_Vec<int> err_flag(1,*ucl_device);
  ucl_copy(err_flag,dev_error,false);
  if (err_flag[0] == 2)
    std::cerr << "BAD MATRIX INVERSION IN FORCE COMPUTATION.\n";  
  err_flag.clear();

  _allocated=false;

  // Output any timing information
  acc_timers();
  double single[9], times[9];

  single[0]=atom->transfer_time()+ans->transfer_time();
  single[1]=nbor->time_nbor.total_seconds();
  single[2]=time_kernel.total_seconds()+time_kernel2.total_seconds()+
            nbor->time_kernel.total_seconds();
  single[3]=time_gayberne.total_seconds()+time_gayberne2.total_seconds();
  if (multiple_forms)
    single[4]=time_pair.total_seconds();
  else
    single[4]=0;
  single[5]=atom->cast_time()+ans->cast_time();
  single[6]=_gpu_overhead;
  single[7]=_driver_overhead;
  single[8]=ans->cpu_idle_time();

  MPI_Reduce(single,times,9,MPI_DOUBLE,MPI_SUM,0,device->replica());
  double avg_split=hd_balancer.all_avg_split();

  _max_bytes+=dev_error.row_bytes()+lj1.row_bytes()+lj3.row_bytes()+
              sigma_epsilon.row_bytes()+cut_form.row_bytes()+
              shape.row_bytes()+well.row_bytes()+lshape.row_bytes()+
              gamma_upsilon_mu.row_bytes()+atom->max_gpu_bytes();
  double mpi_max_bytes;
  MPI_Reduce(&_max_bytes,&mpi_max_bytes,1,MPI_DOUBLE,MPI_MAX,0,
             device->replica());
  double max_mb=mpi_max_bytes/(1024*1024);

  if (device->replica_me()==0)
    if (screen && times[3]>0.0) {
      int replica_size=device->replica_size();

      fprintf(screen,"\n\n-------------------------------------");
      fprintf(screen,"--------------------------------\n");
      fprintf(screen,"      GPU Time Info (average): ");
      fprintf(screen,"\n-------------------------------------");
      fprintf(screen,"--------------------------------\n");

      if (device->procs_per_gpu()==1) {
        fprintf(screen,"Data Transfer:   %.4f s.\n",times[0]/replica_size);
        fprintf(screen,"Data Cast/Pack:  %.4f s.\n",times[5]/replica_size);
        fprintf(screen,"Neighbor copy:   %.4f s.\n",times[1]/replica_size);
        if (nbor->gpu_nbor())
          fprintf(screen,"Neighbor build:  %.4f s.\n",times[2]/replica_size);
        else
          fprintf(screen,"Neighbor unpack: %.4f s.\n",times[2]/replica_size);
        fprintf(screen,"Force calc:      %.4f s.\n",times[3]/replica_size);
        fprintf(screen,"LJ calc:         %.4f s.\n",times[4]/replica_size);
      }
      fprintf(screen,"GPU Overhead:    %.4f s.\n",times[6]/replica_size);
      fprintf(screen,"Average split:   %.4f.\n",avg_split);
      fprintf(screen,"Max Mem / Proc:  %.2f MB.\n",max_mb);
      fprintf(screen,"CPU Driver_Time: %.4f s.\n",times[7]/replica_size);
      fprintf(screen,"CPU Idle_Time:   %.4f s.\n",times[8]/replica_size);
      fprintf(screen,"-------------------------------------");
      fprintf(screen,"--------------------------------\n\n");


      fprintf(screen,"Average split:   %.4f.\n",avg_split);
      fprintf(screen,"Max Mem / Proc:  %.2f MB.\n",max_mb);


    }
  _max_bytes=0.0;

  dev_error.clear();
  lj1.clear();
  lj3.clear();
  sigma_epsilon.clear();
  cut_form.clear();

  shape.clear();
  well.clear();
  lshape.clear();
  gamma_upsilon_mu.clear();
  host_olist.clear();

  time_kernel.clear();
  time_gayberne.clear();
  time_kernel2.clear();
  time_gayberne2.clear();
  time_pair.clear();
  hd_balancer.clear();

  if (_compiled) {
    k_gb_nbor_fast.clear();
    k_gb_nbor.clear();
    k_gayberne.clear();
    k_sphere_gb.clear();
    k_lj_fast.clear();
    k_lj.clear();
    delete pair_program;
    delete gb_program;
    delete gb_lj_program;
    _compiled=false;
  }

  device->clear();
}

template <class numtyp, class acctyp>
double GB_GPU_MemoryT::host_memory_usage() const {
  return device->atom.host_memory_usage()+nbor->host_memory_usage()+
         4*sizeof(numtyp)+sizeof(GB_GPU_Memory<numtyp,acctyp>)+
         nbor->max_atoms()*sizeof(int);
}

template <class numtyp, class acctyp>
void GB_GPU_MemoryT::compile_kernels(UCL_Device &dev) {
  if (_compiled)
    return;

  std::string flags="-cl-fast-relaxed-math -cl-mad-enable "+
                    std::string(OCL_PRECISION_COMPILE);

  pair_program=new UCL_Program(dev);
  pair_program->load_string(gb_gpu_kernel_nbor,flags.c_str());
  k_gb_nbor_fast.set_function(*pair_program,"kernel_gb_nbor_fast");
  k_gb_nbor.set_function(*pair_program,"kernel_gb_nbor");

  gb_program=new UCL_Program(dev);
  gb_program->load_string(gb_gpu_kernel,flags.c_str());
  k_gayberne.set_function(*gb_program,"kernel_gayberne");

  gb_lj_program=new UCL_Program(dev);
  gb_lj_program->load_string(gb_gpu_kernel_lj,flags.c_str());
  k_sphere_gb.set_function(*gb_lj_program,"kernel_sphere_gb");
  k_lj_fast.set_function(*gb_lj_program,"kernel_lj_fast");
  k_lj.set_function(*gb_lj_program,"kernel_lj");

  _compiled=true;
}

template class GB_GPU_Memory<PRECISION,ACC_PRECISION>;

