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

#include "pair_gpu_device.h"
#include "pair_gpu_precision.h"
#include <map>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define PairGPUDeviceT PairGPUDevice<numtyp, acctyp>

template <class numtyp, class acctyp>
PairGPUDeviceT::PairGPUDevice() : _init_count(0), _device_init(false),
                                  _gpu_mode(GPU_FORCE), _first_device(0),
                                  _last_device(0) {
}

template <class numtyp, class acctyp>
PairGPUDeviceT::~PairGPUDevice() {
  clear_device();
}

template <class numtyp, class acctyp>
bool PairGPUDeviceT::init_device(MPI_Comm world, MPI_Comm replica, 
                                 const int first_gpu, const int last_gpu,
                                 const int gpu_mode, const double p_split,
                                 const int nthreads) {
  _nthreads=nthreads;
  #ifdef _OPENMP
  omp_set_num_threads(nthreads);
  #endif

  if (_device_init)
    return true;
  _device_init=true;
  _comm_world=world;
  _comm_replica=replica;
  _first_device=first_gpu;
  _last_device=last_gpu;
  _gpu_mode=gpu_mode;
  _particle_split=p_split;

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
  int my_gpu=node_rank/_procs_per_gpu;
  
  // Set up a per device communicator
  MPI_Comm_split(node_comm,my_gpu,0,&_comm_gpu);
  MPI_Comm_rank(_comm_gpu,&_gpu_rank);

  gpu=new UCL_Device();
  if (my_gpu>=gpu->num_devices())
    return false;
  
  gpu->set(my_gpu);
  return true;
}

template <class numtyp, class acctyp>
bool PairGPUDeviceT::init(PairGPUAns<numtyp,acctyp> &ans, const bool charge,
                          const bool rot, const int nlocal, 
                          const int host_nlocal, const int nall,
                          PairGPUNbor *nbor, const int maxspecial,
                          const bool gpu_nbor, const int gpu_host,
                          const int max_nbors, const double cell_size,
                          const bool pre_cut) {
  if (!_device_init)
    return false;                          

  // Initial number of local particles
  int ef_nlocal=nlocal;
  if (_particle_split<1.0 && _particle_split>0.0)
    ef_nlocal=static_cast<int>(_particle_split*nlocal);

  if (_init_count==0) {
    // Initialize atom and nbor data
    if (!atom.init(nall,charge,rot,*gpu,gpu_nbor,gpu_nbor && maxspecial>0))
      return false;
  } else
    atom.add_fields(charge,rot);

  if (!ans.init(ef_nlocal,charge,rot,*gpu))
    return false;

  if (!nbor->init(&_nbor_shared,ef_nlocal,host_nlocal,max_nbors,maxspecial,
                  *gpu,gpu_nbor,gpu_host,pre_cut))
    return false;
  nbor->cell_size(cell_size);

  _init_count++;
  return true;
}

template <class numtyp, class acctyp>
bool PairGPUDeviceT::init(PairGPUAns<numtyp,acctyp> &ans, const bool charge,
                          const bool rot, const int nlocal, const int nall) {
  if (!_device_init)
    return false;                          

  if (_init_count==0) {
    // Initialize atom and nbor data
    if (!atom.init(nall,charge,rot,*gpu,false,false))
      return false;
  } else
    atom.add_fields(charge,rot);

  if (!ans.init(nlocal,charge,rot,*gpu))
    return false;

  _init_count++;
  return true;
}

template <class numtyp, class acctyp>
void PairGPUDeviceT::init_message(FILE *screen, const char *name,
                                  const int first_gpu, const int last_gpu) {
  #ifdef USE_OPENCL
  std::string fs="";
  #else
  std::string fs=toa(gpu->free_gigabytes())+"/";
  #endif
  
  if (_replica_me == 0 && screen) {
    fprintf(screen,"\n-------------------------------------");
    fprintf(screen,"-------------------------------------\n");
    fprintf(screen,"- Using GPGPU acceleration for %s:\n",name);
    fprintf(screen,"-  with %d procs per device.\n",_procs_per_gpu);
    #ifdef _OPENMP
    fprintf(screen,"-  with %d threads per proc.\n",_nthreads);
    #endif
    fprintf(screen,"-------------------------------------");
    fprintf(screen,"-------------------------------------\n");

    for (int i=first_gpu; i<=last_gpu; i++) {
      std::string sname=gpu->name(i)+", "+toa(gpu->cores(i))+" cores, "+fs+
                        toa(gpu->gigabytes(i))+" GB, "+toa(gpu->clock_rate(i))+
                        " GHZ (";
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
void PairGPUDeviceT::output_times(UCL_Timer &time_pair, 
                                  PairGPUAns<numtyp,acctyp> &ans, 
                                  PairGPUNbor &nbor, const double avg_split, 
                                  const double max_bytes, FILE *screen) {
  double single[5], times[5];

  single[0]=atom.transfer_time()+ans.transfer_time();
  single[1]=nbor.time_nbor.total_seconds();
  single[2]=nbor.time_kernel.total_seconds();
  single[3]=time_pair.total_seconds();
  single[4]=atom.cast_time()+ans.cast_time();

  MPI_Reduce(single,times,5,MPI_DOUBLE,MPI_SUM,0,_comm_replica);

  double my_max_bytes=max_bytes+atom.max_gpu_bytes();
  double mpi_max_bytes;
  MPI_Reduce(&my_max_bytes,&mpi_max_bytes,1,MPI_DOUBLE,MPI_MAX,0,_comm_replica);
  double max_mb=mpi_max_bytes/(1024.0*1024.0);

  if (replica_me()==0)
    if (screen && times[3]>0.0) {
      fprintf(screen,"\n\n-------------------------------------");
      fprintf(screen,"--------------------------------\n");
      fprintf(screen,"      GPU Time Info (average): ");
      fprintf(screen,"\n-------------------------------------");
      fprintf(screen,"--------------------------------\n");

      if (procs_per_gpu()==1) {
        fprintf(screen,"Data Transfer:   %.4f s.\n",times[0]/_replica_size);
        fprintf(screen,"Data Cast/Pack:  %.4f s.\n",times[4]/_replica_size);
        fprintf(screen,"Neighbor copy:   %.4f s.\n",times[1]/_replica_size);
        if (nbor.gpu_nbor())
          fprintf(screen,"Neighbor build:  %.4f s.\n",times[2]/_replica_size);
        else
          fprintf(screen,"Neighbor unpack: %.4f s.\n",times[2]/_replica_size);
        fprintf(screen,"Force calc:      %.4f s.\n",times[3]/_replica_size);
      }
      fprintf(screen,"Average split:   %.4f.\n",avg_split);
      fprintf(screen,"Max Mem / Proc:  %.2f MB.\n",max_mb);

      fprintf(screen,"-------------------------------------");
      fprintf(screen,"--------------------------------\n\n");
    }
}

template <class numtyp, class acctyp>
void PairGPUDeviceT::output_kspace_times(UCL_Timer &time_in, 
                                         UCL_Timer &time_out,
                                         UCL_Timer &time_map,
                                         UCL_Timer &time_rho,
                                         PairGPUAns<numtyp,acctyp> &ans, 
                                         const double max_bytes, FILE *screen) {
  double single[4], times[4];

  single[0]=atom.transfer_time()+ans.transfer_time()+time_in.total_seconds()+
            time_out.total_seconds();
  single[1]=atom.cast_time()+ans.cast_time();
  single[2]=time_map.total_seconds();
  single[3]=time_rho.total_seconds();

  MPI_Reduce(single,times,4,MPI_DOUBLE,MPI_SUM,0,_comm_replica);

  double my_max_bytes=max_bytes+atom.max_gpu_bytes();
  double mpi_max_bytes;
  MPI_Reduce(&my_max_bytes,&mpi_max_bytes,1,MPI_DOUBLE,MPI_MAX,0,_comm_replica);
  double max_mb=mpi_max_bytes/(1024.0*1024.0);

  if (replica_me()==0)
    if (screen && times[3]>0.0) {
      fprintf(screen,"\n\n-------------------------------------");
      fprintf(screen,"--------------------------------\n");
      fprintf(screen,"      GPU Time Info (average): ");
      fprintf(screen,"\n-------------------------------------");
      fprintf(screen,"--------------------------------\n");

      if (procs_per_gpu()==1) {
        fprintf(screen,"Data Transfer:   %.4f s.\n",times[0]/_replica_size);
        fprintf(screen,"Data Cast/Pack:  %.4f s.\n",times[1]/_replica_size);
        fprintf(screen,"Kernel (map):    %.4f s.\n",times[2]/_replica_size);
        fprintf(screen,"Kernel (rho):    %.4f s.\n",times[3]/_replica_size);
//        fprintf(screen,"Force calc:      %.4f s.\n",times[3]/_replica_size);
        fprintf(screen,"Total:           %.4f s.\n",(times[0]+times[1]+times[2]+times[3])/_replica_size);
      }
      fprintf(screen,"Max Mem / Proc:  %.2f MB.\n",max_mb);

      fprintf(screen,"-------------------------------------");
      fprintf(screen,"--------------------------------\n\n");
    }
}

template <class numtyp, class acctyp>
void PairGPUDeviceT::clear() {
  if (_init_count>0) {
    _init_count--;
    if (_init_count==0) {
      atom.clear();
      _nbor_shared.clear();
    }
  }
}

template <class numtyp, class acctyp>
void PairGPUDeviceT::clear_device() {
  while (_init_count>0)
    clear();
  if (_device_init) {
    delete gpu;
    _device_init=false;
  }
}

template <class numtyp, class acctyp>
double PairGPUDeviceT::host_memory_usage() const {
  return atom.host_memory_usage()+4*sizeof(numtyp)+
         sizeof(PairGPUDevice<numtyp,acctyp>);
}

template class PairGPUDevice<PRECISION,ACC_PRECISION>;
PairGPUDevice<PRECISION,ACC_PRECISION> pair_gpu_device;

bool lmp_init_device(MPI_Comm world, MPI_Comm replica, const int first_gpu,
                     const int last_gpu, const int gpu_mode, 
                     const double particle_split, const int nthreads) {
  return pair_gpu_device.init_device(world,replica,first_gpu,last_gpu,gpu_mode,
                                     particle_split,nthreads);
}

void lmp_clear_device() {
  pair_gpu_device.clear_device();
}

double lmp_gpu_forces(double **f, double **tor, double *eatom,
                      double **vatom, double *virial, double &ecoul) {
  return pair_gpu_device.fix_gpu(f,tor,eatom,vatom,virial,ecoul);
}
