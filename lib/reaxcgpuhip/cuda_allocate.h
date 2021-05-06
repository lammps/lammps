#ifndef __CUDA_ALLOCATE_H_
#define __CUDA_ALLOCATE_H_

#include "reaxc_types.h"

#ifdef __cplusplus
extern "C"  {
#endif


void Cuda_Allocate_System( reax_system * );

void Cuda_Allocate_Grid( reax_system * );

void Cuda_Allocate_Simulation_Data( simulation_data * );

void Cuda_Allocate_Workspace( reax_system *, control_params *,storage *, int, int );


void Cuda_Allocate_Control( control_params * );

void Cuda_Deallocate_Grid_Cell_Atoms( reax_system * );

void Cuda_Allocate_Grid_Cell_Atoms( reax_system *, int );

void Cuda_Deallocate_Workspace( control_params *, storage * );


void Cuda_ReAllocate( reax_system*, control_params*, simulation_data*, storage*,
        reax_list**);
void Cuda_Allocate_Simulation_Data( simulation_data *data );



#ifdef __cplusplus
}
#endif

#endif
