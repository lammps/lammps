
#include "cuda_reset_tools.h"

#include "cuda_list.h"
#include "cuda_utils.h"
#include "cuda_reduction.h"

#include "reset_tools.h"


extern "C"
{

void Cuda_Reset_Workspace( reax_system *system , storage *workspace)
{
    cuda_memset( workspace->d_workspace->total_bond_order, 0,
            system->total_cap * sizeof(real), "total_bond_order" );
    cuda_memset( workspace->d_workspace->dDeltap_self, 0,
            system->total_cap * sizeof(rvec), "dDeltap_self" );
    cuda_memset( workspace->d_workspace->CdDelta, 0,
            system->total_cap * sizeof(real), "CdDelta" );
    cuda_memset(workspace->d_workspace->f, 0,
            system->total_cap * sizeof(rvec), "f" );
}


CUDA_GLOBAL void k_reset_hindex( reax_atom *my_atoms, single_body_parameters *sbp,
        int * hindex, int N )
{
    int i;

    i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i >= N )
    {
        return;
    }

    if ( sbp[ my_atoms[i].type ].p_hbond == H_ATOM ||
      sbp[ my_atoms[i].type ].p_hbond == H_BONDING_ATOM )
    {
        hindex[i] = 1;
    }
    else
    {
        hindex[i] = 0;
    }

//    my_atoms[i].Hindex = hindex[i];
    my_atoms[i].Hindex = i;
}


void Cuda_Reset_Atoms( reax_system* system, control_params *control,
        storage *workspace )
{
    int blocks;
    int *hindex;

    hindex = (int *) workspace->scratch;

    blocks = system->n / DEF_BLOCK_SIZE
        + ((system->n % DEF_BLOCK_SIZE == 0 ) ? 0 : 1);

    hipLaunchKernelGGL(k_reset_hindex, dim3(blocks), dim3(DEF_BLOCK_SIZE ), 0, 0,  system->d_my_atoms, system->reax_param.d_sbp, hindex, system->N );
    hipDeviceSynchronize( );
    cudaCheckError( );

    Cuda_Reduction_Sum( hindex, system->d_numH, system->n);

    copy_host_device( &system->numH, system->d_numH, sizeof(int), 
            hipMemcpyDeviceToHost, "Cuda_Reset_Atoms::d_numH" );

    system->Hcap = MAX( (int)(system->numH * SAFER_ZONE), MIN_CAP );

    //printf("Num  H %d, %d \n", system->numH, system->Hcap);

}


void Cuda_Reset( reax_system *system, control_params *control,
        simulation_data *data, storage *workspace, reax_list **lists )
{
    Cuda_Reset_Atoms( system, control, workspace );

    Reset_Simulation_Data_Host( data );

    if ( control->virial )
    {
        Reset_Pressures_Host( data );
    }

    Cuda_Reset_Workspace( system, workspace );

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "p%d @ step%d: reset done\n", system->my_rank, data->step );
    MPI_Barrier( MPI_COMM_WORLD );
#endif

}


}
