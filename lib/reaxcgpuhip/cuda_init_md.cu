
#include "cuda_init_md.h"

#include "cuda_allocate.h"
#include "cuda_list.h"
#include "cuda_copy.h"
#include "cuda_forces.h"
#include "cuda_integrate.h"
#include "cuda_neighbors.h"
#include "cuda_reset_tools.h"
#include "cuda_system_props.h"
#include "cuda_utils.h"

#include "init_md.h"


#include "random.h"
#include "reset_tools.h"
#include "tool_box.h"
#include "vector.h"

static void Cuda_Init_Scratch_Space( storage *workspace )
{
	cuda_malloc( (void **)&workspace->scratch, DEVICE_SCRATCH_SIZE, TRUE,
			"Cuda_Init_Scratch_Space::workspace->scratch" );

	//printf("Scratch size %d \n", DEVICE_SCRATCH_SIZE);

	workspace->host_scratch = (void *) smalloc( HOST_SCRATCH_SIZE,
			"Cuda_Init_Scratch_Space::workspace->host_scratch" );
}


int Cuda_Init_System( reax_system *system, control_params *control,
		simulation_data *data, storage *workspace,
		mpi_datatypes *mpi_data, char *msg )
{
	int i;
	int nrecv[MAX_NBRS];


	reax_atom *atom;

	int mincap = system->mincap;
	double safezone = system->safezone;
	double saferzone = system->saferzone;


	system->local_cap = MAX( (int)(system->n * safezone), mincap);
	system->total_cap = MAX( (int)(system->N * safezone), mincap);

	/* estimate numH and Hcap */
	system->numH = 0;
	if (control->hbond_cut > 0)
		for( i = 0; i < system->n; ++i ) {
			atom = &(system->my_atoms[i]);
			if (system->reax_param.sbp[ atom->type ].p_hbond == 1 && atom->type >= 0)
				atom->Hindex = system->numH++;
			else atom->Hindex = -1;
		}
	system->Hcap = (int)(MAX( system->numH * saferzone, mincap ));

	Cuda_Allocate_System(system);
	Sync_System( system );

	/* estimate numH and Hcap */
	Cuda_Reset_Atoms( system, control, workspace );


	return SUCCESS;
}


void Cuda_Init_Simulation_Data( reax_system *system, control_params *control,
		simulation_data *data )
{
	Cuda_Allocate_Simulation_Data( data );

	Reset_Simulation_Data_Host( data );

	if ( !control->restart )
	{
		data->step = data->prev_steps = 0;
	}



}


void Cuda_Init_Workspace( reax_system *system, control_params *control,
		storage *workspace )
{
	Cuda_Allocate_Workspace( system, control, workspace->d_workspace,
			system->local_cap, system->total_cap );

	memset( &workspace->realloc, 0, sizeof(reallocate_data) );


	Cuda_Reset_Workspace( system, workspace );
	Init_Host_Taper( control, workspace->d_workspace );


}


void Cuda_Init_Lists( reax_system *system, control_params *control,
		simulation_data *data, storage *workspace, reax_list **lists,
		reax_list *cpu_lists,mpi_datatypes *mpi_data )
{

	Cuda_Estimate_Storages( system, control, lists,
			TRUE, TRUE, TRUE, data->step );


	if (control->hbond_cut > 0) {
		Cuda_Make_List( system->total_cap, system->total_hbonds, TYP_HBOND, lists[HBONDS]);
		Cuda_Init_HBond_Indices( system, workspace, lists );
	}


	/* bonds list */

	//printf("total bonds %d, %d \n",system->total_bonds,system->total_hbonds);
	Cuda_Make_List( system->total_cap, system->total_bonds, TYP_BOND, lists[BONDS]);
	Cuda_Init_Bond_Indices( system, lists );


	/* 3bodies list: since a more accurate estimate of the num.
	 * three body interactions requires that bond orders have
	 * been computed, delay estimation until computation */
}



void Cuda_Write_Reax_Lists(reax_system *system, reax_list **gpu_lists, reax_list *cpu_lists) {

	copy_host_device( (cpu_lists+FAR_NBRS)->index, gpu_lists[FAR_NBRS]->index,
			system->total_cap * sizeof(int),
			hipMemcpyHostToDevice, "Output_Sync_Lists::far_neighbor_list" );


	copy_host_device( (cpu_lists+FAR_NBRS)->end_index, gpu_lists[FAR_NBRS]->end_index,
			system->total_cap * sizeof(int),
			hipMemcpyHostToDevice, "Output_Sync_Lists::far_neighbor_list" );



	copy_host_device( (cpu_lists+FAR_NBRS)->select.far_nbr_list, gpu_lists[FAR_NBRS]->select.far_nbr_list,
			system->total_far_nbrs * sizeof(far_neighbor_data),
			hipMemcpyHostToDevice, "Output_Sync_Lists::far_neighbor_list" );


}



void Cuda_Initialize( reax_system *system, control_params *control,
		simulation_data *data, storage *workspace,
		reax_list **lists,reax_list *cpu_lists, output_controls *out_control,
		mpi_datatypes *mpi_data )
{
	char msg[MAX_STR];

	Cuda_Init_Scratch_Space( workspace );


	if ( Cuda_Init_System( system, control, data, workspace, mpi_data, msg ) == FAILURE )
	{
		fprintf( stderr, "[ERROR] p%d: %s\n", system->my_rank, msg );
		fprintf( stderr, "[ERROR] p%d: system could not be initialized! terminating.\n",
				system->my_rank );
		MPI_Abort( MPI_COMM_WORLD, CANNOT_INITIALIZE );
	}

	Cuda_Init_Simulation_Data( system, control, data );

	Cuda_Init_Workspace( system, control, workspace );

	Cuda_Allocate_Control( control );

	Cuda_Init_Lists( system, control, data, workspace, lists,cpu_lists, mpi_data );


	int deviceID;
	hipGetDevice(&deviceID);

	//printf("Device id %d \n", deviceID)
}
