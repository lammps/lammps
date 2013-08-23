/***************************************************************************
                                  device.cpp
                             -------------------
                            W. Michael Brown (ORNL)

  Class for management of the device where the computations are performed

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : brownw@ornl.gov
 ***************************************************************************/

#include "lal_device.h"
#include "lal_precision.h"
#include <map>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(USE_OPENCL)
#include "device_cl.h"
#elif defined(USE_CUDART)
const char *device=0;
#else
#include "device_cubin.h"
#endif

using namespace LAMMPS_AL;
#define DeviceT Device<numtyp, acctyp>

template <class numtyp, class acctyp>
DeviceT::Device() : _init_count(0), _device_init(false),
                                  _gpu_mode(GPU_FORCE), _first_device(0),
                                  _last_device(0), _compiled(false) {
}

template <class numtyp, class acctyp>
DeviceT::~Device() {
  clear_device();
}

template <class numtyp, class acctyp>
int DeviceT::init_device(MPI_Comm world, MPI_Comm replica, const int first_gpu,
                         const int last_gpu, const int gpu_mode, 
                         const double p_split, const int nthreads, 
                         const int t_per_atom, const double cell_size,
                         char *ocl_vendor) {
  _nthreads=nthreads;
  #ifdef _OPENMP
  omp_set_num_threads(nthreads);
  #endif
  _threads_per_atom=t_per_atom;
  _threads_per_charge=t_per_atom;

  if (_device_init)
    return 0;
  _device_init=true;
  _comm_world=replica; //world;
  _comm_replica=replica;
  _first_device=first_gpu;
  _last_device=last_gpu;
  _gpu_mode=gpu_mode;
  _particle_split=p_split;
  _cell_size=cell_size;

  // Get the rank/size within the world
  MPI_Comm_rank(_comm_world,&_world_me);
  MPI_Comm_size(_comm_world,&_world_size);
  // Get the rank/size within the replica
  MPI_Comm_rank(_comm_replica,&_replica_me);
  MPI_Comm_size(_comm_replica,&_replica_size);

  // Get the names of all nodes
  int name_length;
  char node_name[MPI_MAX_PROCESSOR_NAME];
  char node_names[MPI_MAX_PROCESSOR_NAME*_world_size];
  MPI_Get_processor_name(node_name,&name_length);
  MPI_Allgather(&node_name,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,&node_names,
                MPI_MAX_PROCESSOR_NAME,MPI_CHAR,_comm_world);
  std::string node_string=std::string(node_name);
  
  // Get the number of procs per node                
  std::map<std::string,int> name_map;
  std::map<std::string,int>::iterator np;
  for (int i=0; i<_world_size; i++) {
    std::string i_string=std::string(&node_names[i*MPI_MAX_PROCESSOR_NAME]);
    np=name_map.find(i_string);
    if (np==name_map.end())
      name_map[i_string]=1;
    else
      np->second++;
  }
  int procs_per_node=name_map.begin()->second;

  // Assign a unique id to each node
  int split_num=0, split_id=0;
  for (np=name_map.begin(); np!=name_map.end(); ++np) {
    if (np->first==node_string)
      split_id=split_num;
    split_num++;
  }
  
  // Set up a per node communicator and find rank within
  MPI_Comm node_comm;
  MPI_Comm_split(_comm_world, split_id, 0, &node_comm);  
  int node_rank;
  MPI_Comm_rank(node_comm,&node_rank);                  

  // set the device ID
  _procs_per_gpu=static_cast<int>(ceil(static_cast<double>(procs_per_node)/
                                       (last_gpu-first_gpu+1)));
  int my_gpu=node_rank/_procs_per_gpu+first_gpu;

  // Time on the device only if 1 proc per gpu
  _time_device=true;
  if (_procs_per_gpu>1)
    _time_device=false;
  
  // Set up a per device communicator
  MPI_Comm_split(node_comm,my_gpu,0,&_comm_gpu);
  MPI_Comm_rank(_comm_gpu,&_gpu_rank);

  gpu=new UCL_Device();
  if (my_gpu>=gpu->num_devices())
    return -2;
    
  #ifndef CUDA_PROXY
  if (_procs_per_gpu>1 && gpu->sharing_supported(my_gpu)==false)
    return -7;
  #endif
  
  if (gpu->set(my_gpu)!=UCL_SUCCESS)
    return -6;

  gpu->push_command_queue();
  gpu->set_command_queue(1);

  _long_range_precompute=0;

  if (set_ocl_params(ocl_vendor)!=0)
    return -11;
  
  int flag=0;
  for (int i=0; i<_procs_per_gpu; i++) {
    if (_gpu_rank==i)
      flag=compile_kernels();
    gpu_barrier();
  }
  return flag;
}

template <class numtyp, class acctyp>
int DeviceT::set_ocl_params(char *ocl_vendor) {
  #ifdef USE_OPENCL
  std::string s_vendor=OCL_DEFAULT_VENDOR;
  if (ocl_vendor!=NULL)
    s_vendor=ocl_vendor;
  if (s_vendor=="none")
    s_vendor="generic";
  
  if (s_vendor=="kepler") {
    _ocl_vendor_name="NVIDIA Kepler";
    #if defined (__APPLE__) || defined(MACOSX)
    _ocl_vendor_string="-DKEPLER_OCL -DNO_OCL_PTX";
    #else
    _ocl_vendor_string="-DKEPLER_OCL";
    #endif
  } else if (s_vendor=="fermi") {    
    _ocl_vendor_name="NVIDIA Fermi";
    _ocl_vendor_string="-DFERMI_OCL";
  } else if (s_vendor=="cypress") {    
    _ocl_vendor_name="AMD Cypress";
    _ocl_vendor_string="-DCYPRESS_OCL";
  } else if (s_vendor=="generic") {    
    _ocl_vendor_name="GENERIC";
    _ocl_vendor_string="-DGENERIC_OCL";
  } else {
    _ocl_vendor_name="CUSTOM";
    _ocl_vendor_string="-DUSE_OPENCL";
    int token_count=0;
    std::string params[13];
    char *pch = strtok(ocl_vendor,"\" ");
    while (pch != NULL) {
      if (token_count==13)
        return -11;
      params[token_count]=pch;
      token_count++;
      pch = strtok(NULL,"\" ");
    }
    _ocl_vendor_string+=" -DMEM_THREADS="+params[0]+
                        " -DTHREADS_PER_ATOM="+params[1]+
                        " -DTHREADS_PER_CHARGE="+params[2]+
                        " -DBLOCK_PAIR="+params[3]+
                        " -DMAX_SHARED_TYPES="+params[4]+
                        " -DBLOCK_NBOR_BUILD="+params[5]+
                        " -DBLOCK_BIO_PAIR="+params[6]+
                        " -DBLOCK_ELLIPSE="+params[7]+
                        " -DWARP_SIZE="+params[8]+
                        " -DPPPM_BLOCK_1D="+params[9]+
                        " -DBLOCK_CELL_2D="+params[10]+
                        " -DBLOCK_CELL_ID="+params[11]+
                        " -DMAX_BIO_SHARED_TYPES="+params[12];
  }
  _ocl_compile_string="-cl-fast-relaxed-math -cl-mad-enable "+
                      std::string(OCL_PRECISION_COMPILE)+" "+_ocl_vendor_string;
  #endif
  return 0;
}

template <class numtyp, class acctyp>
int DeviceT::init(Answer<numtyp,acctyp> &ans, const bool charge,
                         const bool rot, const int nlocal, 
                         const int host_nlocal, const int nall,
                         Neighbor *nbor, const int maxspecial,
                         const int gpu_host, const int max_nbors, 
                         const double cell_size, const bool pre_cut,
                         const int threads_per_atom) {
  if (!_device_init)
    return -1;
  if (sizeof(acctyp)==sizeof(double) && gpu->double_precision()==false)
    return -5;

  // Counts of data transfers for timing overhead estimates
  _data_in_estimate=0;
  _data_out_estimate=1;

  // Initial number of local particles
  int ef_nlocal=nlocal;
  if (_particle_split<1.0 && _particle_split>0.0)
    ef_nlocal=static_cast<int>(_particle_split*nlocal);

  int gpu_nbor=0;
  if (_gpu_mode==Device<numtyp,acctyp>::GPU_NEIGH)
    gpu_nbor=1;
  else if (_gpu_mode==Device<numtyp,acctyp>::GPU_HYB_NEIGH)
    gpu_nbor=2;
  #ifndef USE_CUDPP
  if (gpu_nbor==1)
    gpu_nbor=2;
  #endif

  if (_init_count==0) {
    // Initialize atom and nbor data
    if (!atom.init(nall,charge,rot,*gpu,gpu_nbor,gpu_nbor>0 && maxspecial>0))
      return -3;
      
    _data_in_estimate++;
    if (charge)
      _data_in_estimate++;
    if (rot)
      _data_in_estimate++;
  } else {
    if (atom.charge()==false && charge)
      _data_in_estimate++;
    if (atom.quaternion()==false && rot)
      _data_in_estimate++;
    if (!atom.add_fields(charge,rot,gpu_nbor,gpu_nbor>0 && maxspecial))
      return -3;
  }
  
  if (!ans.init(ef_nlocal,charge,rot,*gpu))
    return -3;

  if (!nbor->init(&_neighbor_shared,ef_nlocal,host_nlocal,max_nbors,maxspecial,
                  *gpu,gpu_nbor,gpu_host,pre_cut, _block_cell_2d, 
                  _block_cell_id, _block_nbor_build, threads_per_atom,
                  _warp_size, _time_device, compile_string()))
    return -3;
  if (_cell_size<0.0)
    nbor->cell_size(cell_size,cell_size);
  else
    nbor->cell_size(_cell_size,cell_size);

  _init_count++;
  return 0;
}

template <class numtyp, class acctyp>
int DeviceT::init(Answer<numtyp,acctyp> &ans, const int nlocal,
                         const int nall) {
  if (!_device_init)
    return -1;                          
  if (sizeof(acctyp)==sizeof(double) && gpu->double_precision()==false)
    return -5;

  if (_init_count==0) {
    // Initialize atom and nbor data
    if (!atom.init(nall,true,false,*gpu,false,false))
      return -3;
  } else
    if (!atom.add_fields(true,false,false,false))
      return -3;

  if (!ans.init(nlocal,true,false,*gpu))
    return -3;

  _init_count++;
  return 0;
}

template <class numtyp, class acctyp>
void DeviceT::set_single_precompute
                     (PPPM<numtyp,acctyp,float,_lgpu_float4> *pppm) {
  _long_range_precompute=1;
  pppm_single=pppm;
}

template <class numtyp, class acctyp>
void DeviceT::set_double_precompute
                     (PPPM<numtyp,acctyp,double,_lgpu_double4> *pppm) {
  _long_range_precompute=2;
  pppm_double=pppm;
}

template <class numtyp, class acctyp>
void DeviceT::init_message(FILE *screen, const char *name,
                                  const int first_gpu, const int last_gpu) {
  #if defined(USE_OPENCL)
  std::string fs="";
  #elif defined(USE_CUDART)
  std::string fs="";
  #else
  std::string fs=toa(gpu->free_gigabytes())+"/";
  #endif
  
  if (_replica_me == 0 && screen) {
    fprintf(screen,"\n-------------------------------------");
    fprintf(screen,"-------------------------------------\n");
    fprintf(screen,"- Using GPGPU acceleration for %s:\n",name);
    fprintf(screen,"-  with %d proc(s) per device.\n",_procs_per_gpu);
    #ifdef _OPENMP
    fprintf(screen,"-  with %d thread(s) per proc.\n",_nthreads);
    #endif
    #ifdef USE_OPENCL
    fprintf(screen,"-  with OpenCL Parameters for: %s\n",
            _ocl_vendor_name.c_str());
    #endif
    fprintf(screen,"-------------------------------------");
    fprintf(screen,"-------------------------------------\n");

    int last=last_gpu+1;
    if (last>gpu->num_devices())
      last=gpu->num_devices();
    for (int i=first_gpu; i<last; i++) {
      std::string sname;
      if (i==first_gpu)
        sname=gpu->name(i)+", "+toa(gpu->cores(i))+" cores, "+fs+
              toa(gpu->gigabytes(i))+" GB, "+toa(gpu->clock_rate(i))+" GHZ (";
      else              
        sname=gpu->name(i)+", "+toa(gpu->cores(i))+" cores, "+fs+
              toa(gpu->clock_rate(i))+" GHZ (";
      if (sizeof(PRECISION)==4) {
        if (sizeof(ACC_PRECISION)==4)
          sname+="Single Precision)";
        else
          sname+="Mixed Precision)";
      } else
        sname+="Double Precision)";

      fprintf(screen,"GPU %d: %s\n",i,sname.c_str());
    }

    fprintf(screen,"-------------------------------------");
    fprintf(screen,"-------------------------------------\n\n");
  }
}

template <class numtyp, class acctyp>
void DeviceT::estimate_gpu_overhead(const int kernel_calls, 
                                           double &gpu_overhead,
                                           double &gpu_driver_overhead) {
  UCL_H_Vec<int> *host_data_in=NULL, *host_data_out=NULL;
  UCL_D_Vec<int> *dev_data_in=NULL, *dev_data_out=NULL, *kernel_data=NULL;
  UCL_Timer *timers_in=NULL, *timers_out=NULL, *timers_kernel=NULL;
  UCL_Timer over_timer(*gpu);

  if (_data_in_estimate>0) {
    host_data_in=new UCL_H_Vec<int>[_data_in_estimate];
    dev_data_in=new UCL_D_Vec<int>[_data_in_estimate];
    timers_in=new UCL_Timer[_data_in_estimate];
  }
  
  if (_data_out_estimate>0) {
    host_data_out=new UCL_H_Vec<int>[_data_out_estimate];
    dev_data_out=new UCL_D_Vec<int>[_data_out_estimate];
    timers_out=new UCL_Timer[_data_out_estimate];
  }
  
  if (kernel_calls>0) {
    kernel_data=new UCL_D_Vec<int>[kernel_calls];
    timers_kernel=new UCL_Timer[kernel_calls];
  }
  
  for (int i=0; i<_data_in_estimate; i++) {
    host_data_in[i].alloc(1,*gpu);
    dev_data_in[i].alloc(1,*gpu);
    timers_in[i].init(*gpu);
  }  
  
  for (int i=0; i<_data_out_estimate; i++) {
    host_data_out[i].alloc(1,*gpu);
    dev_data_out[i].alloc(1,*gpu);
    timers_out[i].init(*gpu);
  }  
  
  for (int i=0; i<kernel_calls; i++) {
    kernel_data[i].alloc(1,*gpu);
    timers_kernel[i].init(*gpu);
  }  
  
  gpu_overhead=0.0;
  gpu_driver_overhead=0.0;
  
  for (int i=0; i<10; i++) {
    gpu->sync();
    gpu_barrier();
    over_timer.start();
    gpu->sync();
    gpu_barrier();

    double driver_time=MPI_Wtime();
    for (int i=0; i<_data_in_estimate; i++) {
      timers_in[i].start();
      ucl_copy(dev_data_in[i],host_data_in[i],true);
      timers_in[i].stop();
    }
    
    for (int i=0; i<kernel_calls; i++) {
      timers_kernel[i].start();
      zero(kernel_data[i],1);
      timers_kernel[i].stop();
    }

    for (int i=0; i<_data_out_estimate; i++) {
      timers_out[i].start();
      ucl_copy(host_data_out[i],dev_data_out[i],true);
      timers_out[i].stop();
    }
    over_timer.stop();

    double time=over_timer.seconds();
    driver_time=MPI_Wtime()-driver_time;
     
    if (time_device()) {
      for (int i=0; i<_data_in_estimate; i++)
        timers_in[i].add_to_total();
      for (int i=0; i<kernel_calls; i++)
        timers_kernel[i].add_to_total();
      for (int i=0; i<_data_out_estimate; i++)
        timers_out[i].add_to_total();
    }
    
    double mpi_time, mpi_driver_time;
    MPI_Allreduce(&time,&mpi_time,1,MPI_DOUBLE,MPI_MAX,gpu_comm());
    MPI_Allreduce(&driver_time,&mpi_driver_time,1,MPI_DOUBLE,MPI_MAX,gpu_comm());
    gpu_overhead+=mpi_time;
    gpu_driver_overhead+=mpi_driver_time;
  }
  gpu_overhead/=10.0;
  gpu_driver_overhead/=10.0;

  if (_data_in_estimate>0) {
    delete [] host_data_in;
    delete [] dev_data_in;
    delete [] timers_in;
  }
  
  if (_data_out_estimate>0) {
    delete [] host_data_out;
    delete [] dev_data_out;
    delete [] timers_out;
  }
  
  if (kernel_calls>0) {
    delete [] kernel_data;
    delete [] timers_kernel;
  }
}              

template <class numtyp, class acctyp>
void DeviceT::output_times(UCL_Timer &time_pair, Answer<numtyp,acctyp> &ans, 
                           Neighbor &nbor, const double avg_split, 
                           const double max_bytes, const double gpu_overhead,
                           const double driver_overhead, 
                           const int threads_per_atom, FILE *screen) {
  double single[9], times[9];

  single[0]=atom.transfer_time()+ans.transfer_time();
  single[1]=nbor.time_nbor.total_seconds()+nbor.time_hybrid1.total_seconds()+
            nbor.time_hybrid2.total_seconds();
  single[2]=nbor.time_kernel.total_seconds();
  single[3]=time_pair.total_seconds();
  single[4]=atom.cast_time()+ans.cast_time();
  single[5]=gpu_overhead;
  single[6]=driver_overhead;
  single[7]=ans.cpu_idle_time();
  single[8]=nbor.bin_time();

  MPI_Reduce(single,times,9,MPI_DOUBLE,MPI_SUM,0,_comm_replica);

  double my_max_bytes=max_bytes+atom.max_gpu_bytes();
  double mpi_max_bytes;
  MPI_Reduce(&my_max_bytes,&mpi_max_bytes,1,MPI_DOUBLE,MPI_MAX,0,_comm_replica);
  double max_mb=mpi_max_bytes/(1024.0*1024.0);

  if (replica_me()==0)
    if (screen && times[5]>0.0) {
      fprintf(screen,"\n\n-------------------------------------");
      fprintf(screen,"--------------------------------\n");
      fprintf(screen,"      GPU Time Info (average): ");
      fprintf(screen,"\n-------------------------------------");
      fprintf(screen,"--------------------------------\n");

      if (time_device()) {
        fprintf(screen,"Data Transfer:   %.4f s.\n",times[0]/_replica_size);
        fprintf(screen,"Data Cast/Pack:  %.4f s.\n",times[4]/_replica_size);
        fprintf(screen,"Neighbor copy:   %.4f s.\n",times[1]/_replica_size);
        if (nbor.gpu_nbor()>0)
          fprintf(screen,"Neighbor build:  %.4f s.\n",times[2]/_replica_size);
        else
          fprintf(screen,"Neighbor unpack: %.4f s.\n",times[2]/_replica_size);
        fprintf(screen,"Force calc:      %.4f s.\n",times[3]/_replica_size);
      }
      if (nbor.gpu_nbor()==2)
        fprintf(screen,"Neighbor (CPU):  %.4f s.\n",times[8]/_replica_size);
      fprintf(screen,"GPU Overhead:    %.4f s.\n",times[5]/_replica_size);
      fprintf(screen,"Average split:   %.4f.\n",avg_split);
      fprintf(screen,"Threads / atom:  %d.\n",threads_per_atom);
      fprintf(screen,"Max Mem / Proc:  %.2f MB.\n",max_mb);
      fprintf(screen,"CPU Driver_Time: %.4f s.\n",times[6]/_replica_size);
      fprintf(screen,"CPU Idle_Time:   %.4f s.\n",times[7]/_replica_size);

      fprintf(screen,"-------------------------------------");
      fprintf(screen,"--------------------------------\n\n");
    }
}

template <class numtyp, class acctyp>
void DeviceT::output_kspace_times(UCL_Timer &time_in, 
                                         UCL_Timer &time_out,
                                         UCL_Timer &time_map,
                                         UCL_Timer &time_rho,
                                         UCL_Timer &time_interp,
                                         Answer<numtyp,acctyp> &ans, 
                                         const double max_bytes, 
                                         const double cpu_time, 
                                         const double idle_time, FILE *screen) {
  double single[8], times[8];

  single[0]=time_out.total_seconds();
  single[1]=time_in.total_seconds()+atom.transfer_time()+atom.cast_time();
  single[2]=time_map.total_seconds();
  single[3]=time_rho.total_seconds();
  single[4]=time_interp.total_seconds();
  single[5]=ans.transfer_time()+ans.cast_time();
  single[6]=cpu_time;
  single[7]=idle_time;

  MPI_Reduce(single,times,8,MPI_DOUBLE,MPI_SUM,0,_comm_replica);

  double my_max_bytes=max_bytes+atom.max_gpu_bytes();
  double mpi_max_bytes;
  MPI_Reduce(&my_max_bytes,&mpi_max_bytes,1,MPI_DOUBLE,MPI_MAX,0,_comm_replica);
  double max_mb=mpi_max_bytes/(1024.0*1024.0);

  if (replica_me()==0)
    if (screen && times[6]>0.0) {
      fprintf(screen,"\n\n-------------------------------------");
      fprintf(screen,"--------------------------------\n");
      fprintf(screen,"      GPU Time Info (average): ");
      fprintf(screen,"\n-------------------------------------");
      fprintf(screen,"--------------------------------\n");

      if (time_device()) {
        fprintf(screen,"Data Out:        %.4f s.\n",times[0]/_replica_size);
        fprintf(screen,"Data In:         %.4f s.\n",times[1]/_replica_size);
        fprintf(screen,"Kernel (map):    %.4f s.\n",times[2]/_replica_size);
        fprintf(screen,"Kernel (rho):    %.4f s.\n",times[3]/_replica_size);
        fprintf(screen,"Force interp:    %.4f s.\n",times[4]/_replica_size);
        fprintf(screen,"Total rho:       %.4f s.\n",
                (times[0]+times[2]+times[3])/_replica_size);
        fprintf(screen,"Total interp:    %.4f s.\n",
                (times[1]+times[4])/_replica_size);
        fprintf(screen,"Force copy/cast: %.4f s.\n",times[5]/_replica_size);
        fprintf(screen,"Total:           %.4f s.\n",
                (times[0]+times[1]+times[2]+times[3]+times[4]+times[5])/
                _replica_size);
      }
      fprintf(screen,"CPU Poisson:     %.4f s.\n",times[6]/_replica_size);
      fprintf(screen,"CPU Idle Time:   %.4f s.\n",times[7]/_replica_size);
      fprintf(screen,"Max Mem / Proc:  %.2f MB.\n",max_mb);

      fprintf(screen,"-------------------------------------");
      fprintf(screen,"--------------------------------\n\n");
    }
}

template <class numtyp, class acctyp>
void DeviceT::clear() {
  if (_init_count>0) {
    _long_range_precompute=0;
    _init_count--;
    if (_init_count==0) {
      atom.clear();
      _neighbor_shared.clear();
    }
  }
}

template <class numtyp, class acctyp>
void DeviceT::clear_device() {
  while (_init_count>0)
    clear();
  if (_compiled) {
    k_zero.clear();
    k_info.clear();
    delete dev_program;
    _compiled=false;
  }
  if (_device_init) {
    delete gpu;
    _device_init=false;
  }
}

template <class numtyp, class acctyp>
int DeviceT::compile_kernels() {
  int flag=0;

  if (_compiled)
  	return flag;
  	
  dev_program=new UCL_Program(*gpu);
  int success=dev_program->load_string(device,compile_string().c_str());
  if (success!=UCL_SUCCESS)
    return -4;
  k_zero.set_function(*dev_program,"kernel_zero");
  k_info.set_function(*dev_program,"kernel_info");
  _compiled=true;

  UCL_Vector<int,int> gpu_lib_data(15,*gpu,UCL_NOT_PINNED);
  k_info.set_size(1,1);
  k_info.run(&gpu_lib_data);
  gpu_lib_data.update_host(false);
  
  _ptx_arch=static_cast<double>(gpu_lib_data[0])/100.0;
  #ifndef USE_OPENCL
  if (_ptx_arch>gpu->arch() || floor(_ptx_arch)<floor(gpu->arch()))
    return -4;
  #endif

  _num_mem_threads=gpu_lib_data[1];
  _warp_size=gpu_lib_data[2];
  if (_threads_per_atom<1)
    _threads_per_atom=gpu_lib_data[3];
  if (_threads_per_charge<1)
    _threads_per_charge=gpu_lib_data[13];
  _pppm_max_spline=gpu_lib_data[4];
  _pppm_block=gpu_lib_data[5];
  _block_pair=gpu_lib_data[6];
  _max_shared_types=gpu_lib_data[7];
  _block_cell_2d=gpu_lib_data[8];
  _block_cell_id=gpu_lib_data[9];
  _block_nbor_build=gpu_lib_data[10];
  _block_bio_pair=gpu_lib_data[11];
  _max_bio_shared_types=gpu_lib_data[12];
  _block_ellipse=gpu_lib_data[14];

  if (static_cast<size_t>(_block_pair)>gpu->group_size())
    _block_pair=gpu->group_size();
  if (static_cast<size_t>(_block_bio_pair)>gpu->group_size())
    _block_bio_pair=gpu->group_size();
  if (_threads_per_atom>_warp_size)
    _threads_per_atom=_warp_size;
  if (_warp_size%_threads_per_atom!=0)
    _threads_per_atom=1;
  if (_threads_per_atom & (_threads_per_atom - 1))
    _threads_per_atom=1;
  if (_threads_per_charge>_warp_size)
    _threads_per_charge=_warp_size;
  if (_warp_size%_threads_per_charge!=0)
    _threads_per_charge=1;
  if (_threads_per_charge & (_threads_per_charge - 1))
    _threads_per_charge=1;

  return flag;    
}

template <class numtyp, class acctyp>
double DeviceT::host_memory_usage() const {
  return atom.host_memory_usage()+4*sizeof(numtyp)+
         sizeof(Device<numtyp,acctyp>);
}

template class Device<PRECISION,ACC_PRECISION>;
Device<PRECISION,ACC_PRECISION> global_device;

int lmp_init_device(MPI_Comm world, MPI_Comm replica, const int first_gpu,
                    const int last_gpu, const int gpu_mode, 
                    const double particle_split, const int nthreads,
                    const int t_per_atom, const double cell_size, 
                    char *opencl_vendor) {
  return global_device.init_device(world,replica,first_gpu,last_gpu,gpu_mode,
                                   particle_split,nthreads,t_per_atom, 
                                   cell_size,opencl_vendor);
}

void lmp_clear_device() {
  global_device.clear_device();
}

double lmp_gpu_forces(double **f, double **tor, double *eatom,
                      double **vatom, double *virial, double &ecoul) {
  return global_device.fix_gpu(f,tor,eatom,vatom,virial,ecoul);
}

