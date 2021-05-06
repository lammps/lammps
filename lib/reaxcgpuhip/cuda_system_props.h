
#ifndef __CUDA_SYSTEM_PROPS_H__
#define __CUDA_SYSTEM_PROPS_H__

#include "reaxc_types.h"


#ifdef __cplusplus
extern "C"  {
#endif

void Cuda_Compute_Total_Mass( reax_system *, control_params *,
        storage *, simulation_data *, MPI_Comm );

void Cuda_Sync_Simulation_Data( simulation_data * );

void Cuda_Generate_Initial_Velocities( reax_system *, real );

void Cuda_Compute_Kinetic_Energy( reax_system *, control_params *,
        storage *, simulation_data *, MPI_Comm );

void Cuda_Compute_Center_of_Mass( reax_system *, control_params *,
        storage *, simulation_data *, mpi_datatypes *, MPI_Comm );

void Cuda_Compute_Pressure( reax_system *, control_params *,
        storage *, simulation_data *, mpi_datatypes * );

#ifdef __cplusplus
}
#endif


#endif
