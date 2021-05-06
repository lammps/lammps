
#include "cuda_copy.h"

#include "cuda_utils.h"

#include "list.h"
#include "vector.h"
#include "tool_box.h"


/* Copy grid info from host to device */
void Sync_Grid( grid *host, grid *device )
{
	int total;

	total = host->ncells[0] * host->ncells[1] * host->ncells[2];

	ivec_Copy( device->ncells, host->ncells);
	rvec_Copy( device->cell_len, host->cell_len);
	rvec_Copy( device->inv_len, host->inv_len);

	ivec_Copy( device->bond_span, host->bond_span );
	ivec_Copy( device->nonb_span, host->nonb_span );
	ivec_Copy( device->vlist_span, host->vlist_span );

	ivec_Copy( device->native_cells, host->native_cells );
	ivec_Copy( device->native_str, host->native_str );
	ivec_Copy( device->native_end, host->native_end );

	device->ghost_cut = host->ghost_cut;
	ivec_Copy( device->ghost_span, host->ghost_span );
	ivec_Copy( device->ghost_nonb_span, host->ghost_nonb_span );
	ivec_Copy( device->ghost_hbond_span, host->ghost_hbond_span );
	ivec_Copy( device->ghost_bond_span, host->ghost_bond_span );

	copy_host_device( host->str, device->str, sizeof(int) * total,
			hipMemcpyHostToDevice, "grid:str" );
	copy_host_device( host->end, device->end, sizeof(int) * total,
			hipMemcpyHostToDevice, "grid:end" );
	copy_host_device( host->cutoff, device->cutoff, sizeof(real) * total,
			hipMemcpyHostToDevice, "grid:cutoff" );
	copy_host_device( host->nbrs_x, device->nbrs_x, sizeof(ivec) * total *
			host->max_nbrs, hipMemcpyHostToDevice, "grid:nbrs_x" );
	copy_host_device( host->nbrs_cp, device->nbrs_cp, sizeof(rvec) * total *
			host->max_nbrs, hipMemcpyHostToDevice, "grid:nbrs_cp" );

	copy_host_device( host->rel_box, device->rel_box, sizeof(ivec) * total,
			hipMemcpyHostToDevice, "grid:rel_box" );

	device->max_nbrs = host->max_nbrs;
}


void Output_Sync_Forces(storage *workspace, int total_cap)
{
	copy_host_device(workspace->f, workspace->d_workspace->f,
			total_cap * sizeof(rvec), hipMemcpyDeviceToHost,
			"Output_Sync_Atoms::my_atoms" );

}

/* Copy atom info from host to device */
void Sync_Atoms( reax_system *sys )
{
	//TODO METIN FIX, coredump on his machine
	//    copy_host_device( sys->my_atoms, sys->d_my_atoms, sizeof(reax_atom) * sys->total_cap,
	//            hipMemcpyHostToDevice, "Sync_Atoms::system->my_atoms" );

#if defined(__CUDA_DEBUG_LOG__)
	fprintf( stderr, "p:%d - Synching atoms: n: %d N: %d, total_cap: %d \n",
			sys->my_rank, sys->n, sys->N, sys->total_cap );
#endif

	copy_host_device( sys->my_atoms, sys->d_my_atoms, sizeof(reax_atom) * sys->N,
			hipMemcpyHostToDevice, "Sync_Atoms::system->my_atoms" );
	//TODO METIN FIX, coredump on his machine
}


/* Copy atomic system info from host to device */
void Sync_System( reax_system *sys )
{

	Sync_Atoms(sys);
	copy_host_device( sys->reax_param.sbp, sys->reax_param.d_sbp,
			sizeof(single_body_parameters) * sys->reax_param.num_atom_types,
			hipMemcpyHostToDevice, "Sync_System::system->sbp" );
	copy_host_device( sys->reax_param.tbp, sys->reax_param.d_tbp,
			sizeof(two_body_parameters) * POW(sys->reax_param.num_atom_types, 2),
			hipMemcpyHostToDevice, "Sync_System::system->tbp" );
	copy_host_device( sys->reax_param.thbp, sys->reax_param.d_thbp,
			sizeof(three_body_header) * POW(sys->reax_param.num_atom_types, 3),
			hipMemcpyHostToDevice, "Sync_System::system->thbh" );
	copy_host_device( sys->reax_param.hbp, sys->reax_param.d_hbp,
			sizeof(hbond_parameters) * POW(sys->reax_param.num_atom_types, 3),
			hipMemcpyHostToDevice, "Sync_System::system->hbond" );
	copy_host_device( sys->reax_param.fbp, sys->reax_param.d_fbp,
			sizeof(four_body_header) * POW(sys->reax_param.num_atom_types, 4),
			hipMemcpyHostToDevice, "Sync_System::system->four_header" );

	copy_host_device( sys->reax_param.gp.l, sys->reax_param.d_gp.l,
			sizeof(real) * sys->reax_param.gp.n_global, hipMemcpyHostToDevice,
			"Sync_System::system->global_parameters" );

	sys->reax_param.d_gp.n_global = sys->reax_param.gp.n_global;
	sys->reax_param.d_gp.vdw_type = sys->reax_param.gp.vdw_type;
}


/* Copy atom info from device to host */
void Output_Sync_Atoms( reax_system *sys )
{
	copy_host_device( sys->my_atoms, sys->d_my_atoms,
			sizeof(reax_atom) * sys->total_cap, hipMemcpyDeviceToHost,
			"Output_Sync_Atoms::my_atoms" );
}


/* Copy simulation data from device to host */
void Output_Sync_Simulation_Data( simulation_data *host, simulation_data *dev )
{
	copy_host_device( &host->my_en, &dev->my_en, sizeof(energy_data),
			hipMemcpyDeviceToHost, "simulation_data:energy_data" );
	copy_host_device( &host->kin_press, &dev->kin_press, sizeof(real),
			hipMemcpyDeviceToHost, "simulation_data:kin_press" );
	copy_host_device(host->int_press, dev->int_press, sizeof(rvec),
			hipMemcpyDeviceToHost, "simulation_data:int_press" );
	copy_host_device( host->ext_press, dev->ext_press, sizeof(rvec),
			hipMemcpyDeviceToHost, "simulation_data:ext_press" );
}


/* Copy interaction lists from device to host,
 * with allocation for the host list */
void Output_Sync_Lists( reax_list *l, reax_list *device_list, int type )
{

   if ( l->allocated == TRUE )
	{
		fprintf( stderr, "[WARNING] attempted to allocate list which was already allocated."
				" Returning without allocation...\n" );
		return;
	}

	int n = device_list->n;
	int num_intrs = device_list->num_intrs;


	l->allocated = TRUE;
	l->n = n;
	l->num_intrs = num_intrs;
	l->type = type;

	l->index = (int*)malloc( sizeof(int) * n);
	l->end_index = (int*)malloc( sizeof(int) * n);


	switch ( l->type )
	{
	case TYP_VOID:
		l->select.v = smalloc( sizeof(void*) * l->num_intrs, "Make_List::v" );
		break;

	case TYP_BOND:
		l->select.bond_list = (bond_data*)smalloc( sizeof(bond_data) * l->num_intrs, "Make_List::bonds" );
		break;

	case TYP_THREE_BODY:
		l->select.three_body_list = (three_body_interaction_data*)smalloc( sizeof(three_body_interaction_data) * l->num_intrs,
				"Make_List::three_bodies" );
		break;

	case TYP_HBOND:
		l->select.hbond_list = (hbond_data*)smalloc( sizeof(hbond_data) * l->num_intrs, "Make_List::hbonds" );
		break;

	case TYP_FAR_NEIGHBOR:
		l->select.far_nbr_list = (far_neighbor_data*)smalloc( sizeof(far_neighbor_data) * l->num_intrs,
				"Make_List::far_nbrs" );
		break;

	case TYP_DBO:
		l->select.dbo_list = (dbond_data*)smalloc( sizeof(dbond_data) * l->num_intrs, "Make_List::dbonds" );
		break;

	case TYP_DDELTA:
		l->select.dDelta_list = (dDelta_data*)smalloc( sizeof(dDelta_data) * l->num_intrs, "Make_List::dDeltas" );
		break;

	default:
		fprintf( stderr, "[ERROR] unknown list type (%d)\n", l->type );
		MPI_Abort( MPI_COMM_WORLD, INVALID_INPUT );
		break;
	}


	//printf("List allocated\n");


	copy_host_device( l->index, device_list->index, sizeof(int) *  device_list->n,
	       hipMemcpyDeviceToHost, "Output_Sync_Lists::list->index" );
    copy_host_device( l->end_index, device_list->end_index, sizeof(int) *
	       device_list->n, hipMemcpyDeviceToHost, "Output_Sync_Lists::list->end_index" );


	 switch ( type )
    {   
        case TYP_FAR_NEIGHBOR:
            copy_host_device(l->select.far_nbr_list, device_list->select.far_nbr_list,
                    sizeof(far_neighbor_data) * device_list->num_intrs,
                    hipMemcpyDeviceToHost, "Output_Sync_Lists::far_neighbor_list" );
            break;

        case TYP_BOND:
            copy_host_device( l->select.bond_list, device_list->select.bond_list,
                    sizeof(bond_data) * device_list->num_intrs,
                    hipMemcpyDeviceToHost, "Output_Sync_Lists::bond_list" );
            break;

        case TYP_HBOND:
            copy_host_device( l->select.hbond_list, device_list->select.hbond_list,
                    sizeof(hbond_data) * device_list->num_intrs,
                    hipMemcpyDeviceToHost, "Output_Sync_Lists::hbond_list" );
            break;

        case TYP_THREE_BODY:
            copy_host_device(l->select.three_body_list,
                    device_list->select.three_body_list,
                    sizeof(three_body_interaction_data ) * device_list->num_intrs,
                    hipMemcpyDeviceToHost, "Output_Sync_Lists::three_body_list" );
            break;

        default:
            fprintf( stderr, "[ERROR] Unknown list synching from device to host (%d)\n",
                    type );
            exit( INVALID_INPUT );
            break;
    }

	//printf("Copy finish\n");

}
