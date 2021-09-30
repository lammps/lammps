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
#include <cmath>
#include <cstdlib>
#include <iostream>
#if (LAL_USE_OMP == 1)
#include <omp.h>
#endif

#if defined(USE_OPENCL)
#include "device_cl.h"

#ifdef LAL_OCL_EXTRA_ARGS
#define LAL_DM_STRINGIFY(x) #x
#define LAL_PRE_STRINGIFY(x) LAL_DM_STRINGIFY(x)
#endif

#elif defined(USE_CUDART)
const char *device=0;
#else
#include "device_cubin.h"
#endif

namespace LAMMPS_AL {
#define DeviceT Device<numtyp, acctyp>

template <class numtyp, class acctyp>
DeviceT::Device() : _init_count(0), _device_init(false),
                    _gpu_mode(GPU_FORCE), _first_device(0),
                    _last_device(0), _platform_id(-1), _compiled(false) {
}

template <class numtyp, class acctyp>
DeviceT::~Device() {
  clear_device();
}

template <class numtyp, class acctyp>
int DeviceT::init_device(MPI_Comm world, MPI_Comm replica, const int ngpu,
                         const int first_gpu_id, const int gpu_mode,
                         const double p_split, const int t_per_atom,
                         const double user_cell_size, char *ocl_args,
                         const int ocl_platform, char *device_type_flags,
                         const int block_pair) {
  _threads_per_atom=t_per_atom;
  _threads_per_charge=t_per_atom;
  _threads_per_three=t_per_atom;

  if (_device_init)
    return 0;
  _device_init=true;
  _comm_world=replica; //world;
  _comm_replica=replica;
  int ndevices=ngpu;
  _first_device=first_gpu_id;
  _gpu_mode=gpu_mode;
  _particle_split=p_split;
  _user_cell_size=user_cell_size;
  _block_pair=block_pair;

  // support selecting OpenCL platform id with "package platform" keyword
  if (ocl_platform >= 0)
    _platform_id = ocl_platform;

  gpu=new UCL_Device();

  // ---------------------- OpenCL Compiler Args -------------------------
  std::string extra_args="";
  if (ocl_args) extra_args+=":"+std::string(ocl_args);
  #ifdef LAL_OCL_EXTRA_ARGS
  extra_args+=":" LAL_PRE_STRINGIFY(LAL_OCL_EXTRA_ARGS);
  #endif
  for (int i=0; i<extra_args.length(); i++)
    if (extra_args[i]==':') extra_args[i]=' ';

  // --------------------------- MPI setup -------------------------------

  // Get the rank/size within the world
  MPI_Comm_rank(_comm_world,&_world_me);
  MPI_Comm_size(_comm_world,&_world_size);
  // Get the rank/size within the replica
  MPI_Comm_rank(_comm_replica,&_replica_me);
  MPI_Comm_size(_comm_replica,&_replica_size);

  // Get the names of all nodes
  int name_length;
  char node_name[MPI_MAX_PROCESSOR_NAME];
  char *node_names = new char[MPI_MAX_PROCESSOR_NAME*_world_size];
  MPI_Get_processor_name(node_name,&name_length);
  MPI_Allgather(&node_name,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,&node_names[0],
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
  delete[] node_names;

  // Set up a per node communicator and find rank within
  MPI_Comm node_comm;
  MPI_Comm_split(_comm_world, split_id, 0, &node_comm);
  int node_rank;
  MPI_Comm_rank(node_comm,&node_rank);

  // ------------------- Device selection parameters----------------------

  if (ndevices > procs_per_node)
    ndevices = procs_per_node;

  // --------------------- OCL Platform Selection  -----------------------

  // Setup OpenCL platform and parameters based on platform
  // and device type specifications
  std::string ocl_vstring="";
  if (device_type_flags != nullptr) ocl_vstring=device_type_flags;

  // Setup the OpenCL platform
  // If multiple platforms and no user platform specified,
  // try to match platform from config matching any user specified
  // device type. Give preference to platforms with GPUs.
  // Priority under these conditions to platform with device with
  // highest compute unit count.
  int pres;
  enum UCL_DEVICE_TYPE type=UCL_GPU;
  #ifndef USE_OPENCL
  pres=gpu->set_platform(0);
  #else
  if (_platform_id>=0)
    pres=gpu->set_platform(_platform_id);
  else {
    std::string vendor="";
    if (device_type_flags!=nullptr) {
      if (ocl_vstring=="intelgpu")
        vendor="intel";
      else if (ocl_vstring=="intelcpu") {
        vendor="intel";
        type=UCL_CPU;
      } else if (ocl_vstring=="nvidiagpu")
        vendor="nvidia";
      else if (ocl_vstring=="amdgpu")
        vendor="amd";
      else if (ocl_vstring=="applegpu")
        vendor="apple";
    }
    pres=gpu->auto_set_platform(type,vendor,ndevices,_first_device);
  }
  #endif
  if (pres != UCL_SUCCESS)
    return -12;

  // ------------------------ Device Selection ---------------------------
  if (_first_device > -1 && _first_device >= gpu->num_devices())
    return -2;
  if (ndevices > gpu->num_devices())
    return -2;
  if (_first_device + ndevices > gpu->num_devices())
    return -2;
  if (gpu->num_devices()==0)
    return -2;

  // Fully specified deviceIDs
  if (_first_device > -1 && ndevices > 0)
    _last_device = _first_device + ndevices - 1;

  // Find deviceID with most CUs (priority given to the accelerator type)
  if (_first_device < 0) {
    int best_device = 0;
    int best_cus = gpu->cus(0);
    bool type_match = (gpu->device_type(0) == type);
    for (int i = 1; i < gpu->num_devices(); i++) {
      if (type_match==true && gpu->device_type(i)!=type)
        continue;
      if (type_match == false && gpu->device_type(i) == type) {
        type_match = true;
        best_cus = gpu->cus(i);
        best_device = i;
      }
      if (gpu->cus(i) > best_cus) {
        best_cus = gpu->cus(i);
        best_device = i;
      }
    }
    _first_device = _last_device = best_device;
    type = gpu->device_type(_first_device);

    if (ndevices > 0) {
      // Expand range to meet specified number of devices
      while (_last_device - _first_device < ndevices - 1) {
        if (_last_device + 1 == gpu->num_devices())
          _first_device--;
        else if (_first_device == 0)
          _last_device++;
        else {
          if (gpu->device_type(_last_device+1)==type &&
              gpu->device_type(_first_device-1)!=type)
            _last_device++;
          else if (gpu->device_type(_last_device+1)!=type &&
                   gpu->device_type(_first_device-1)==type)
            _first_device--;
          else if (gpu->cus(_last_device+1) > gpu->cus(_first_device-1))
            _last_device++;
          else
            _first_device--;
        }
      }
    }
  }

  // If ngpus not specified, expand range to include matching devices
  if (ndevices == 0) {
    for (int i = _first_device; i < gpu->num_devices(); i++) {
      if (gpu->device_type(i)==gpu->device_type(_first_device) &&
          gpu->cus(i)==gpu->cus(_first_device))
        _last_device = i;
      else
        break;
    }
    ndevices = _last_device - _first_device + 1;
    if (ndevices > procs_per_node) {
      ndevices = procs_per_node;
      _last_device=_first_device + ndevices - 1;
    }
  }

  // ------------------------ MPI Device ID Setup -----------------------

  // set the device ID
  _procs_per_gpu=static_cast<int>(ceil(static_cast<double>(procs_per_node)/
                                       ndevices));
  int my_gpu=node_rank/_procs_per_gpu+_first_device;

  // Time on the device only if 1 proc per gpu
  _time_device=true;

#if 0
  // XXX: the following setting triggers a memory leak with OpenCL and MPI
  //      setting _time_device=true for all processes doesn't seem to be a
  //      problem with either (no segfault, no (large) memory leak.
  //      thus keeping this disabled for now. may need to review later.
  //      2018-07-23 <akohlmey@gmail.com>
  if (_procs_per_gpu>1)
    _time_device=false;
#endif

  // Set up a per device communicator
  MPI_Comm_split(node_comm,my_gpu,0,&_comm_gpu);
  MPI_Comm_rank(_comm_gpu,&_gpu_rank);

  #if !defined(CUDA_PROXY) && !defined(CUDA_MPS_SUPPORT)
  if (_procs_per_gpu>1 && gpu->sharing_supported(my_gpu)==false)
    return -7;
  #endif

  // --------------- Device Configuration and Setup  -------------------------

  if (gpu->set(my_gpu)!=UCL_SUCCESS)
    return -6;

  #if !defined(USE_OPENCL) && !defined(USE_HIP)
  if (gpu->arch()<7.0) {
    gpu->push_command_queue();
    gpu->set_command_queue(1);
  }
  #endif

  _long_range_precompute=0;

  // If OpenCL parameters not specified by user, try to auto detect
  // best option from the platform config
  #ifdef USE_OPENCL
  if (device_type_flags==nullptr) {
    std::string pname = gpu->platform_name();
    for (int i=0; i<pname.length(); i++)
      if (pname[i]<='z' && pname[i]>='a')
        pname[i]=toupper(pname[i]);
    if (pname.find("NVIDIA")!=std::string::npos)
      ocl_vstring="nvidiagpu";
    else if (pname.find("INTEL")!=std::string::npos) {
      if (gpu->device_type()==UCL_GPU)
        ocl_vstring="intelgpu";
      else if (gpu->device_type()==UCL_CPU)
        ocl_vstring="intelcpu";
    } else if (pname.find("AMD")!=std::string::npos) {
      if (gpu->device_type()==UCL_GPU)
        ocl_vstring="amdgpu";
    } else if (pname.find("APPLE")!=std::string::npos) {
      if (gpu->device_type()==UCL_GPU)
        ocl_vstring="applegpu";
    }
  }
  #endif

  if (set_ocl_params(ocl_vstring, extra_args)!=0)
    return -11;

  int flag=0;
  for (int i=0; i<_procs_per_gpu; i++) {
    if (_gpu_rank==i)
      flag=compile_kernels();
    gpu_barrier();
  }

  // check if double precision support is available
  #if defined(_SINGLE_DOUBLE) || defined(_DOUBLE_DOUBLE)
  if (!gpu->double_precision())
    return -16;
  #endif

  // Setup auto bin size calculation for calls from atom::sort
  // - This is repeated in neighbor init with additional info
  if (_user_cell_size<0.0) {
    #ifndef LAL_USE_OLD_NEIGHBOR
    _neighbor_shared.setup_auto_cell_size(true,0,_simd_size);
    #else
    _neighbor_shared.setup_auto_cell_size(false,0,_simd_size);
    #endif
  } else
    _neighbor_shared.setup_auto_cell_size(false,_user_cell_size,_simd_size);

  return flag;
}

template <class numtyp, class acctyp>
int DeviceT::set_ocl_params(std::string s_config, const std::string &extra_args) {
  #ifdef USE_OPENCL

  #include "lal_pre_ocl_config.h"

  if (s_config=="" || s_config=="none")
    s_config="generic";

  int config_index=-1;
  for (int i=0; i<nconfigs; i++)
    if (s_config==std::string(ocl_config_names[i]))
      config_index=i;

  if (config_index != -1)
    s_config=ocl_config_strings[config_index];

  _ocl_config_name="CUSTOM";
  int token_count=0;
  std::string params[18];
  char ocl_config[2048];
  strncpy(ocl_config,s_config.c_str(),2047);
  char *pch = strtok(ocl_config,",");
  _ocl_config_name=pch;
  pch = strtok(nullptr,",");
  if (pch == nullptr) return -11;
  while (pch != nullptr) {
    if (token_count==18)
      return -11;
    params[token_count]=pch;
    token_count++;
    pch = strtok(nullptr,",");
  }

  _ocl_compile_string="-cl-mad-enable ";
  if (params[4]!="0") _ocl_compile_string+="-cl-fast-relaxed-math ";
  _ocl_compile_string+=std::string(OCL_INT_TYPE)+" "+
    std::string(OCL_PRECISION_COMPILE);
  if (gpu->has_subgroup_support())
    _ocl_compile_string+=" -DUSE_OPENCL_SUBGROUPS";
  #ifdef LAL_USE_OLD_NEIGHBOR
  _ocl_compile_string+=" -DLAL_USE_OLD_NEIGHBOR";
  #endif

  _ocl_compile_string += " -DCONFIG_ID="+params[0]+
                         " -DSIMD_SIZE="+params[1]+
                         " -DMEM_THREADS="+params[2];
  if (gpu->has_shuffle_support()==false)
    _ocl_compile_string+=" -DSHUFFLE_AVAIL=0";
  else
    _ocl_compile_string+=" -DSHUFFLE_AVAIL="+params[3];
  _ocl_compile_string += " -DFAST_MATH="+params[4]+

                         " -DTHREADS_PER_ATOM="+params[5]+
                         " -DTHREADS_PER_CHARGE="+params[6]+
                         " -DTHREADS_PER_THREE="+params[7]+

                         " -DBLOCK_PAIR="+params[8]+
                         " -DBLOCK_BIO_PAIR="+params[9]+
                         " -DBLOCK_ELLIPSE="+params[10]+
                         " -DPPPM_BLOCK_1D="+params[11]+
                         " -DBLOCK_NBOR_BUILD="+params[12]+
                         " -DBLOCK_CELL_2D="+params[13]+
                         " -DBLOCK_CELL_ID="+params[14]+

                         " -DMAX_SHARED_TYPES="+params[15]+
                         " -DMAX_BIO_SHARED_TYPES="+params[16]+
                         " -DPPPM_MAX_SPLINE="+params[17];
  _ocl_compile_string += extra_args;
  #endif
  return 0;
}

template <class numtyp, class acctyp>
std::string DeviceT::compile_string_nofast() {
  std::string no_fast = _ocl_compile_string;
  size_t p = no_fast.find("-cl-fast-relaxed-math ");
  if (p != std::string::npos) no_fast.erase(p,22);
  p = no_fast.find("-DFAST_MATH=");
  if (p != std::string::npos) no_fast[p + 12]='0';
  return no_fast;
}

template <class numtyp, class acctyp>
int DeviceT::init(Answer<numtyp,acctyp> &ans, const bool charge,
                  const bool rot, const int nlocal,
                  const int nall, const int maxspecial,
                  const bool vel) {
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
  #if !defined(USE_CUDPP) && !defined(USE_HIP_DEVICE_SORT)
  if (gpu_nbor==1) gpu_nbor=2;
  #endif
  #ifndef LAL_USE_OLD_NEIGHBOR
  if (gpu_nbor==1) gpu_nbor=2;
  #endif

  if (_init_count==0) {
    // Initialize atom and nbor data
    if (!atom.init(nall,charge,rot,*gpu,gpu_nbor,gpu_nbor>0 && maxspecial>0,vel))
      return -3;

    _data_in_estimate++;
    if (charge)
      _data_in_estimate++;
    if (rot)
      _data_in_estimate++;
    if (vel)
      _data_in_estimate++;
  } else {
    if (atom.charge()==false && charge)
      _data_in_estimate++;
    if (atom.quaternion()==false && rot)
      _data_in_estimate++;
    if (atom.velocity()==false && vel)
      _data_in_estimate++;
    if (!atom.add_fields(charge,rot,gpu_nbor,gpu_nbor>0 && maxspecial,vel))
      return -3;
  }

  if (!ans.init(ef_nlocal,charge,rot,*gpu))
    return -3;

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
int DeviceT::init_nbor(Neighbor *nbor, const int nlocal,
                       const int host_nlocal, const int nall,
                       const int maxspecial, const int gpu_host,
                       const int max_nbors, const double cutoff,
                       const bool pre_cut, const int threads_per_atom,
                       const bool ilist_map) {
  int ef_nlocal=nlocal;
  if (_particle_split<1.0 && _particle_split>0.0)
    ef_nlocal=static_cast<int>(_particle_split*nlocal);

  int gpu_nbor=0;
  if (_gpu_mode==Device<numtyp,acctyp>::GPU_NEIGH)
    gpu_nbor=1;
  else if (_gpu_mode==Device<numtyp,acctyp>::GPU_HYB_NEIGH)
    gpu_nbor=2;
  #if !defined(USE_CUDPP) && !defined(USE_HIP_DEVICE_SORT)
  if (gpu_nbor==1)
    gpu_nbor=2;
  #endif
  #ifndef LAL_USE_OLD_NEIGHBOR
  if (gpu_nbor==1)
    gpu_nbor=2;
  #endif

  if (!nbor->init(&_neighbor_shared,ef_nlocal,host_nlocal,max_nbors,maxspecial,
                  *gpu,gpu_nbor,gpu_host,pre_cut,_block_cell_2d,
                  _block_cell_id, _block_nbor_build, threads_per_atom,
                  _simd_size, _time_device, compile_string(), ilist_map))
    return -3;

  if (_user_cell_size<0.0) {
    _neighbor_shared.setup_auto_cell_size(false,cutoff,nbor->simd_size());
  } else
    _neighbor_shared.setup_auto_cell_size(false,_user_cell_size,nbor->simd_size());
  nbor->set_cutoff(cutoff);

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
    fprintf(screen,"- Using acceleration for %s:\n",name);
    fprintf(screen,"-  with %d proc(s) per device.\n",_procs_per_gpu);
    #if (LAL_USE_OMP == 1)
    fprintf(screen,"-  with %d thread(s) per proc.\n", omp_get_max_threads());
    #endif
    #ifdef USE_OPENCL
    fprintf(screen,"-  with OpenCL Parameters for: %s (%d)\n",
            _ocl_config_name.c_str(),_config_id);
    #endif
    if (shuffle_avail())
      fprintf(screen,"-  Horizontal vector operations: ENABLED\n");
    else
      fprintf(screen,"-  Horizontal vector operations: DISABLED\n");
    if (gpu->shared_memory(first_gpu))
      fprintf(screen,"-  Shared memory system: Yes\n");
    else
      fprintf(screen,"-  Shared memory system: No\n");
    fprintf(screen,"-------------------------------------");
    fprintf(screen,"-------------------------------------\n");

    int last=last_gpu+1;
    if (last>gpu->num_devices())
      last=gpu->num_devices();
    for (int i=first_gpu; i<last; i++) {
      std::string sname;
      if (i==first_gpu)
        sname=gpu->name(i)+", "+toa(gpu->cus(i))+" CUs, "+fs+
              toa(gpu->gigabytes(i))+" GB, "+toa(gpu->clock_rate(i))+" GHZ (";
      else
        sname=gpu->name(i)+", "+toa(gpu->cus(i))+" CUs, "+
              toa(gpu->clock_rate(i))+" GHZ (";
      if (sizeof(PRECISION)==4) {
        if (sizeof(ACC_PRECISION)==4)
          sname+="Single Precision)";
        else
          sname+="Mixed Precision)";
      } else
        sname+="Double Precision)";

      fprintf(screen,"Device %d: %s\n",i,sname.c_str());
    }

    fprintf(screen,"-------------------------------------");
    fprintf(screen,"-------------------------------------\n\n");
  }
}

template <class numtyp, class acctyp>
void DeviceT::estimate_gpu_overhead(const int kernel_calls,
                                    double &gpu_overhead,
                                    double &gpu_driver_overhead) {
  UCL_H_Vec<int> *host_data_in=nullptr, *host_data_out=nullptr;
  UCL_D_Vec<int> *dev_data_in=nullptr, *dev_data_out=nullptr,
    *kernel_data=nullptr;
  UCL_Timer *timers_in=nullptr, *timers_out=nullptr, *timers_kernel=nullptr;
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

  for (int z=0; z<11; z++) {
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

    const int numel=1;
    for (int i=0; i<kernel_calls; i++) {
      timers_kernel[i].start();
      k_zero.set_size(1,_block_pair);
      k_zero.run(&(kernel_data[i]),&numel);
      timers_kernel[i].stop();
    }

    for (int i=0; i<_data_out_estimate; i++) {
      timers_out[i].start();
      ucl_copy(host_data_out[i],dev_data_out[i],true);
      timers_out[i].stop();
    }
    over_timer.stop();
    #ifndef GERYON_OCL_FLUSH
    if (_data_out_estimate)
      dev_data_out[0].flush();
    #endif
    driver_time=MPI_Wtime()-driver_time;
    double time=over_timer.seconds();

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
    MPI_Allreduce(&driver_time,&mpi_driver_time,1,MPI_DOUBLE,MPI_MAX,
                  gpu_comm());
    if (z>0) {
      gpu_overhead+=mpi_time;
      gpu_driver_overhead+=mpi_driver_time;
    }
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
  int post_final=0;

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

  MPI_Finalized(&post_final);
  if (post_final) return;

  MPI_Reduce(single,times,9,MPI_DOUBLE,MPI_SUM,0,_comm_replica);

  double my_max_bytes=max_bytes+atom.max_gpu_bytes();
  double mpi_max_bytes;
  MPI_Reduce(&my_max_bytes,&mpi_max_bytes,1,MPI_DOUBLE,MPI_MAX,0,_comm_replica);
  double max_mb=mpi_max_bytes/(1024.0*1024.0);

  #ifdef USE_OPENCL
  // Workaround for timing issue on Intel OpenCL
  if (times[0] > 80e6) times[0]=0.0;
  if (times[3] > 80e6) times[3]=0.0;
  if (times[5] > 80e6) times[5]=0.0;
  #endif

  if (replica_me()==0)
    if (screen && (times[6] > 0.0)) {
      fprintf(screen,"\n\n-------------------------------------");
      fprintf(screen,"--------------------------------\n");
      fprintf(screen,"      Device Time Info (average): ");
      fprintf(screen,"\n-------------------------------------");
      fprintf(screen,"--------------------------------\n");

      if (time_device() && (times[3] > 0.0)) {
        if (times[0] > 0.0)
          fprintf(screen,"Data Transfer:   %.4f s.\n",times[0]/_replica_size);
        fprintf(screen,"Neighbor copy:   %.4f s.\n",times[1]/_replica_size);
        if (nbor.gpu_nbor() > 0.0)
          fprintf(screen,"Neighbor build:  %.4f s.\n",times[2]/_replica_size);
        else
          fprintf(screen,"Neighbor unpack: %.4f s.\n",times[2]/_replica_size);
        fprintf(screen,"Force calc:      %.4f s.\n",times[3]/_replica_size);
      }
      if (times[5] > 0.0)
        fprintf(screen,"Device Overhead: %.4f s.\n",times[5]/_replica_size);
      fprintf(screen,"Average split:   %.4f.\n",avg_split);
      fprintf(screen,"Lanes / atom:    %d.\n",threads_per_atom);
      fprintf(screen,"Vector width:    %d.\n", simd_size());
      fprintf(screen,"Max Mem / Proc:  %.2f MB.\n",max_mb);
      if (nbor.gpu_nbor()==2)
        fprintf(screen,"CPU Neighbor:    %.4f s.\n",times[8]/_replica_size);
      fprintf(screen,"CPU Cast/Pack:   %.4f s.\n",times[4]/_replica_size);
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
  double single[9], times[9];

  single[0]=time_out.total_seconds();
  single[1]=time_in.total_seconds()+atom.transfer_time()+atom.cast_time();
  single[2]=time_map.total_seconds();
  single[3]=time_rho.total_seconds();
  single[4]=time_interp.total_seconds();
  single[5]=ans.transfer_time();
  single[6]=cpu_time;
  single[7]=idle_time;
  single[8]=ans.cast_time();

  MPI_Reduce(single,times,9,MPI_DOUBLE,MPI_SUM,0,_comm_replica);

  double my_max_bytes=max_bytes+atom.max_gpu_bytes();
  double mpi_max_bytes;
  MPI_Reduce(&my_max_bytes,&mpi_max_bytes,1,MPI_DOUBLE,MPI_MAX,0,_comm_replica);
  double max_mb=mpi_max_bytes/(1024.0*1024.0);
  #ifdef USE_OPENCL
  // Workaround for timing issue on Intel OpenCL
  if (times[3] > 80e6) times[3]=0.0;
  #endif


  if (replica_me()==0)
    if (screen && times[6]>0.0) {
      fprintf(screen,"\n\n-------------------------------------");
      fprintf(screen,"--------------------------------\n");
      fprintf(screen,"    Device Time Info (average) for kspace: ");
      fprintf(screen,"\n-------------------------------------");
      fprintf(screen,"--------------------------------\n");

      if (time_device() && times[3]>0) {
        fprintf(screen,"Data Out:        %.4f s.\n",times[0]/_replica_size);
        fprintf(screen,"Data In:         %.4f s.\n",times[1]/_replica_size);
        fprintf(screen,"Kernel (map):    %.4f s.\n",times[2]/_replica_size);
        fprintf(screen,"Kernel (rho):    %.4f s.\n",times[3]/_replica_size);
        fprintf(screen,"Force interp:    %.4f s.\n",times[4]/_replica_size);
        fprintf(screen,"Total rho:       %.4f s.\n",
                (times[0]+times[2]+times[3])/_replica_size);
        fprintf(screen,"Total interp:    %.4f s.\n",
                (times[1]+times[4])/_replica_size);
        fprintf(screen,"Force copy:      %.4f s.\n",times[5]/_replica_size);
        fprintf(screen,"Total:           %.4f s.\n",
                (times[0]+times[1]+times[2]+times[3]+times[4]+times[5])/
                _replica_size);
      }
      fprintf(screen,"CPU Poisson:     %.4f s.\n",times[6]/_replica_size);
      fprintf(screen,"CPU Data Cast:   %.4f s.\n",times[8]/_replica_size);
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
  int success=dev_program->load_string(device,compile_string().c_str(),
                                       nullptr,stderr);
  if (success!=UCL_SUCCESS)
    return -6;
  k_zero.set_function(*dev_program,"kernel_zero");
  k_info.set_function(*dev_program,"kernel_info");
  _compiled=true;

  UCL_Vector<int,int> gpu_lib_data(19,*gpu,UCL_NOT_PINNED);
  k_info.set_size(1,1);
  k_info.run(&gpu_lib_data);
  gpu_lib_data.update_host(false);

  _ptx_arch=static_cast<double>(gpu_lib_data[0])/100.0;
  #if !(defined(USE_OPENCL) || defined(USE_HIP))
  if (_ptx_arch>gpu->arch() || floor(_ptx_arch)<floor(gpu->arch()))
    return -4;
  #endif

  _config_id=gpu_lib_data[1];

  if (sizeof(numtyp)==sizeof(float))
    _simd_size=std::max(gpu_lib_data[2],gpu->preferred_fp32_width());
  else
    _simd_size=std::max(gpu_lib_data[2],gpu->preferred_fp64_width());

  _num_mem_threads=gpu_lib_data[3];
  _shuffle_avail=gpu_lib_data[4];
  _fast_math=gpu_lib_data[5];

  if (_threads_per_atom<1)
    _threads_per_atom=gpu_lib_data[6];
  if (_threads_per_charge<1)
    _threads_per_charge=gpu_lib_data[7];
  if (_threads_per_three<1)
    _threads_per_three=gpu_lib_data[8];

  if (_block_pair == -1) {
    _block_pair=gpu_lib_data[9];
    _block_bio_pair=gpu_lib_data[10];
    _block_ellipse=gpu_lib_data[11];
  } else {
    _block_bio_pair=_block_pair;
    _block_ellipse=_block_pair;
  }
  _pppm_block=gpu_lib_data[12];
  _block_nbor_build=gpu_lib_data[13];
  _block_cell_2d=gpu_lib_data[14];
  _block_cell_id=gpu_lib_data[15];

  _max_shared_types=gpu_lib_data[16];
  _max_bio_shared_types=gpu_lib_data[17];
  _pppm_max_spline=gpu_lib_data[18];

  if (static_cast<size_t>(_block_pair) > gpu->group_size_dim(0) ||
      static_cast<size_t>(_block_bio_pair) > gpu->group_size_dim(0) ||
      static_cast<size_t>(_block_ellipse) > gpu->group_size_dim(0) ||
      static_cast<size_t>(_pppm_block) > gpu->group_size_dim(0) ||
      static_cast<size_t>(_block_nbor_build) > gpu->group_size_dim(0) ||
      static_cast<size_t>(_block_cell_2d) > gpu->group_size_dim(0) ||
      static_cast<size_t>(_block_cell_2d) > gpu->group_size_dim(1) ||
      static_cast<size_t>(_block_cell_id) > gpu->group_size_dim(0) ||
      static_cast<size_t>(_max_shared_types*_max_shared_types*sizeof(numtyp)*17 > gpu->slm_size()) ||
      static_cast<size_t>(_max_bio_shared_types*2*sizeof(numtyp) > gpu->slm_size()))
    return -13;

  if (_block_pair % _simd_size != 0 || _block_bio_pair % _simd_size != 0 ||
      _block_ellipse % _simd_size != 0 || _pppm_block % _simd_size != 0 ||
      _block_nbor_build % _simd_size != 0 ||
      _block_pair < _max_shared_types * _max_shared_types ||
      _block_bio_pair * 2 < _max_bio_shared_types ||
      _pppm_block < _pppm_max_spline * _pppm_max_spline)
    return -11;

  if (_threads_per_atom>_simd_size)
    _threads_per_atom=_simd_size;
  if (_simd_size%_threads_per_atom!=0)
    _threads_per_atom=1;
  if (_threads_per_atom & (_threads_per_atom - 1))
    _threads_per_atom=1;
  if (_threads_per_charge>_simd_size)
    _threads_per_charge=_simd_size;
  if (_simd_size%_threads_per_charge!=0)
    _threads_per_charge=1;
  if (_threads_per_charge & (_threads_per_charge - 1))
    _threads_per_charge=1;
  if (_threads_per_three>_simd_size)
    _threads_per_three=_simd_size;
  if (_simd_size%_threads_per_three!=0)
    _threads_per_three=1;
  if (_threads_per_three & (_threads_per_three - 1))
    _threads_per_three=1;

  return flag;
}

template <class numtyp, class acctyp>
double DeviceT::host_memory_usage() const {
  return atom.host_memory_usage()+4*sizeof(numtyp)+
         sizeof(Device<numtyp,acctyp>);
}

template class Device<PRECISION,ACC_PRECISION>;
Device<PRECISION,ACC_PRECISION> global_device;
}

using namespace LAMMPS_AL;

bool lmp_has_gpu_device()
{
  UCL_Device gpu;
  return (gpu.num_platforms() > 0);
}

std::string lmp_gpu_device_info()
{
  std::ostringstream out;
  UCL_Device gpu;
  if (gpu.num_platforms() > 0)
    gpu.print_all(out);
  return out.str();
}

int lmp_init_device(MPI_Comm world, MPI_Comm replica, const int ngpu,
                    const int first_gpu_id, const int gpu_mode,
                    const double particle_split, const int t_per_atom,
                    const double user_cell_size, char *opencl_config,
                    const int ocl_platform, char *device_type_flags,
                    const int block_pair) {
  return global_device.init_device(world,replica,ngpu,first_gpu_id,gpu_mode,
                                   particle_split,t_per_atom,user_cell_size,
                                   opencl_config,ocl_platform,
                                   device_type_flags,block_pair);
}

void lmp_clear_device() {
  global_device.clear_device();
}

double lmp_gpu_forces(double **f, double **tor, double *eatom, double **vatom,
                      double *virial, double &ecoul, int &error_flag) {
  return global_device.fix_gpu(f,tor,eatom,vatom,virial,ecoul,error_flag);
}

double lmp_gpu_update_bin_size(const double subx, const double suby,
                               const double subz, const int nlocal,
                               const double cut) {
  return global_device._neighbor_shared.update_cell_size(subx, suby,
                                                         subz, nlocal, cut);
}

bool lmp_gpu_config(const std::string &category, const std::string &setting)
{
  if (category == "api") {
#if defined(USE_OPENCL)
    return setting == "opencl";
#elif defined(USE_HIP)
    return setting == "hip";
#elif defined(USE_CUDA) || defined(USE_CUDART)
    return setting == "cuda";
#endif
    return false;
  }
  if (category == "precision") {
    if (setting == "single") {
#if defined(_SINGLE_SINGLE)
      return true;
#else
      return false;
#endif
    } else if (setting == "mixed") {
#if defined(_SINGLE_DOUBLE)
      return true;
#else
      return false;
#endif
    } else if (setting == "double") {
#if defined(_DOUBLE_DOUBLE)
      return true;
#else
      return false;
#endif
    } else return false;
  }
  return false;
}
