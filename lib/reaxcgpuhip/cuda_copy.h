#ifndef __CUDA_COPY_H_
#define __CUDA_COPY_H_

#include "reaxc_types.h"


#ifdef __cplusplus
extern "C"  {
#endif

void Sync_Atoms( reax_system * );

void Sync_Grid( grid *, grid * );

void Sync_System( reax_system * );

void Prep_Device_For_Output( reax_system *, simulation_data * );

void Output_Sync_Lists( reax_list *host, reax_list *device, int type );

void Output_Sync_Atoms( reax_system * );

void Output_Sync_Forces(storage *workspace, int total_cap);


void Output_Sync_Simulation_Data( simulation_data *, simulation_data * );

#ifdef __cplusplus
}
#endif


#endif
