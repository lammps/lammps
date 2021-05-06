/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, haktulga@cs.purdue.edu
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "reaxc_types.h"


#include "reset_tools.h"
#include "list.h"
#include "tool_box.h"
#include "vector.h"
#include "index_utils.h"


static void Reset_Atoms_Host( reax_system * const system, control_params * const control )
{
    int i;
    reax_atom *atom;

    system->numH = 0;
    if ( control->hbond_cut > 0.0 )
    {
        for ( i = 0; i < system->N; ++i )
        {
            atom = &system->my_atoms[i];

            if( system->reax_param.sbp[ atom->type ].p_hbond == H_ATOM )
            {
                atom->Hindex = system->numH++;
            }
            else
            {
                atom->Hindex = -1;
            }
        }
    }
}


static void Reset_Energies_Host( energy_data * const en )
{
    en->e_bond = 0.0;
    en->e_ov = 0.0;
    en->e_un = 0.0;
    en->e_lp = 0.0;
    en->e_ang = 0.0;
    en->e_pen = 0.0;
    en->e_coa = 0.0;
    en->e_hb = 0.0;
    en->e_tor = 0.0;
    en->e_con = 0.0;
    en->e_vdW = 0.0;
    en->e_ele = 0.0;

    en->e_pot = 0.0;
    en->e_kin = 0.0;
    en->e_tot = 0.0;
}


static void Reset_Temperatures_Host( simulation_data * const data )
{
    data->therm.T = 0.0;
}


void Reset_Pressures_Host( simulation_data * const data )
{
    data->flex_bar.P_scalar = 0.0;
    rtensor_MakeZero( data->flex_bar.P );

    data->iso_bar.P = 0.0;
    rvec_MakeZero( data->int_press );
    rvec_MakeZero( data->my_ext_press );
    rvec_MakeZero( data->ext_press );
}


void Reset_Simulation_Data_Host( simulation_data * const data )
{
    Reset_Energies_Host( &data->my_en );
    Reset_Energies_Host( &data->sys_en );
    Reset_Temperatures_Host( data );
    Reset_Pressures_Host( data );
}


void Reset_Timing_Host( reax_timing * const rt )
{
   // rt->total = Get_Time( );
    rt->comm = 0.0;
    rt->nbrs = 0.0;
    rt->init_forces = 0.0;
    rt->bonded = 0.0;
    rt->nonb = 0.0;
    rt->cm = 0.0;
    rt->cm_solver_iters = 0;
    rt->num_retries = 0;
}


#ifdef TEST_FORCES
void Reset_Test_Forces( reax_system * const system, storage * const workspace )
{
    memset( workspace->f_ele, 0, system->total_cap * sizeof(rvec) );
    memset( workspace->f_vdw, 0, system->total_cap * sizeof(rvec) );
    memset( workspace->f_bo, 0, system->total_cap * sizeof(rvec) );
    memset( workspace->f_be, 0, system->total_cap * sizeof(rvec) );
    memset( workspace->f_lp, 0, system->total_cap * sizeof(rvec) );
    memset( workspace->f_ov, 0, system->total_cap * sizeof(rvec) );
    memset( workspace->f_un, 0, system->total_cap * sizeof(rvec) );
    memset( workspace->f_ang, 0, system->total_cap * sizeof(rvec) );
    memset( workspace->f_coa, 0, system->total_cap * sizeof(rvec) );
    memset( workspace->f_pen, 0, system->total_cap * sizeof(rvec) );
    memset( workspace->f_hb, 0, system->total_cap * sizeof(rvec) );
    memset( workspace->f_tor, 0, system->total_cap * sizeof(rvec) );
    memset( workspace->f_con, 0, system->total_cap * sizeof(rvec) );
}
#endif


void Reset_Workspace_Host( reax_system * const system, storage * const workspace )
{
    memset( workspace->total_bond_order, 0, system->total_cap * sizeof( real ) );
    memset( workspace->dDeltap_self, 0, system->total_cap * sizeof( rvec ) );
    memset( workspace->CdDelta, 0, system->total_cap * sizeof( real ) );
    memset( workspace->f, 0, system->total_cap * sizeof( rvec ) );

#ifdef TEST_FORCES
    memset( workspace->dDelta, 0, sizeof(rvec) * system->total_cap );
    Reset_Test_Forces( system, workspace );
#endif
}


void Reset_Grid_Host( grid * const g )
{
    int i, j, k;

    for ( i = 0; i < g->ncells[0]; i++ )
    {
        for ( j = 0; j < g->ncells[1]; j++ )
        {
            for ( k = 0; k < g->ncells[2]; k++ )
            {
                g->cells[ index_grid_3d(i, j, k, g) ].top = 0;
                //g->cells[ index_grid_3d(i, j, k, g) ].str = 0;
                //g->cells[ index_grid_3d(i, j, k, g) ].end = 0;
            }
        }
    }
}


void Reset_Out_Buffers_Host( mpi_out_data * const out_buf, int n )
{
    int i;

    for ( i = 0; i < n; ++i )
    {
        out_buf[i].cnt = 0;
    }
}


void Reset_Lists_Host (reax_system * const system, control_params * const control,
        storage * const workspace, reax_list ** const lists )
{
    int i;
    reax_list * const bond_list = lists[BONDS];
    reax_list * const hbond_list = lists[HBONDS];

    if ( system->N > 0 )
    {
        for ( i = 0; i < system->total_cap; ++i )
        {
            Set_End_Index( i, Start_Index( i, bond_list ), bond_list );
        }

        if ( control->hbond_cut > 0.0 && system->numH > 0 )
        {
            for ( i = 0; i < system->total_cap; ++i )
            {
                /* do not use Hindex, unconditionally reset end indices */
                Set_End_Index( i, Start_Index( i, hbond_list ), hbond_list );
            }
        }

        for ( i = 0; i < system->local_cap; ++i )
        {
            workspace->H.end[i] = workspace->H.start[i];
        }
    }
}


void Reset_Host( reax_system * const system, control_params * const control,
        simulation_data * const data, storage * const workspace,
        reax_list ** const lists )
{
    Reset_Atoms_Host( system, control );

    Reset_Simulation_Data_Host( data );

    if ( control->virial )
    {
        Reset_Pressures_Host( data );
    }

    Reset_Workspace_Host( system, workspace );

    Reset_Lists_Host( system, control, workspace, lists );

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "p%d @ step%d: reset done\n", system->my_rank, data->step );
    MPI_Barrier( MPI_COMM_WORLD );
#endif
}
