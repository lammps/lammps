/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

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

#include "pair_reax_c.h"
#if defined(PURE_REAX)
#include "init_md.h"
#include "allocate.h"
#include "box.h"
#include "comm_tools.h"
#include "forces.h"
#include "grid.h"
#include "integrate.h"
#include "io_tools.h"
#include "list.h"
#include "lookup.h"
#include "neighbors.h"
#include "random.h"
#include "reset_tools.h"
#include "system_props.h"
#include "tool_box.h"
#include "vector.h"
#elif defined(LAMMPS_REAX)
#include "reaxc_init_md.h"
#include "reaxc_allocate.h"
#include "reaxc_forces.h"
#include "reaxc_io_tools.h"
#include "reaxc_list.h"
#include "reaxc_lookup.h"
#include "reaxc_reset_tools.h"
#include "reaxc_system_props.h"
#include "reaxc_tool_box.h"
#include "reaxc_vector.h"
#endif


#if defined(PURE_REAX)
/************************ initialize system ************************/
int Reposition_Atoms( reax_system *system, control_params *control,
                      simulation_data *data, mpi_datatypes *mpi_data,
                      char *msg )
{
  int   i;
  rvec  dx;

  /* reposition atoms */
  if( control->reposition_atoms == 0 ) { //fit atoms to periodic box
    rvec_MakeZero( dx );
  }
  else if( control->reposition_atoms == 1 ) { //put center of mass to center
    rvec_Scale( dx, 0.5, system->big_box.box_norms );
    rvec_ScaledAdd( dx, -1., data->xcm );
  }
  else if( control->reposition_atoms == 2 ) { //put center of mass to origin
    rvec_Scale( dx, -1., data->xcm );
  }
  else {
    strcpy( msg, "reposition_atoms: invalid option" );
    return FAILURE;
  }

  for( i = 0; i < system->n; ++i )
    // Inc_on_T3_Gen( system->my_atoms[i].x, dx, &(system->big_box) );
    rvec_Add( system->my_atoms[i].x, dx );

  return SUCCESS;
}



void Generate_Initial_Velocities( reax_system *system, real T )
{
  int i;
  real m, scale, norm;


  if( T <= 0.1 ) {
    for( i = 0; i < system->n; i++ )
      rvec_MakeZero( system->my_atoms[i].v );
  }
  else {
    Randomize();

    for( i = 0; i < system->n; i++ ) {
      rvec_Random( system->my_atoms[i].v );

      norm = rvec_Norm_Sqr( system->my_atoms[i].v );
      m = system->reax_param.sbp[ system->my_atoms[i].type ].mass;
      scale = sqrt( m * norm / (3.0 * K_B * T) );

      rvec_Scale( system->my_atoms[i].v, 1./scale, system->my_atoms[i].v );

      // fprintf( stderr, "v = %f %f %f\n",
      // system->my_atoms[i].v[0],
      // system->my_atoms[i].v[1],
      // system->my_atoms[i].v[2] );

      // fprintf( stderr, "scale = %f\n", scale );
      // fprintf( stderr, "v = %f %f %f\n",
      // system->my_atoms[i].v[0],
      // system->my_atoms[i].v[1],
      // system->my_atoms[i].v[2] );
    }
  }
}


int Init_System( reax_system *system, control_params *control,
                 simulation_data *data, storage *workspace,
                 mpi_datatypes *mpi_data, char *msg )
{
  int i;
  reax_atom *atom;
  int nrecv[MAX_NBRS];

  Setup_New_Grid( system, control, mpi_data->world );
#if defined(DEBUG_FOCUS)
  fprintf( stderr, "p%d GRID:\n", system->my_rank );
  Print_Grid( &(system->my_grid), stderr );
#endif
  Bin_My_Atoms( system, &(workspace->realloc) );
  Reorder_My_Atoms( system, workspace );

  /* estimate N and total capacity */
  for( i = 0; i < MAX_NBRS; ++i ) nrecv[i] = 0;
  system->max_recved = 0;
  system->N = SendRecv( system, mpi_data, mpi_data->boundary_atom_type, nrecv,
                        Estimate_Boundary_Atoms, Unpack_Estimate_Message, 1 );
  system->total_cap = MAX( (int)(system->N * SAFE_ZONE), MIN_CAP );
  Bin_Boundary_Atoms( system );

  /* estimate numH and Hcap */
  system->numH = 0;
  if( control->hbond_cut > 0 )
    for( i = 0; i < system->n; ++i ) {
      atom = &(system->my_atoms[i]);
      if( system->reax_param.sbp[ atom->type ].p_hbond == 1 )
        atom->Hindex = system->numH++;
      else atom->Hindex = -1;
    }
  system->Hcap = MAX( system->numH * SAFER_ZONE, MIN_CAP );

  //Allocate_System( system, system->local_cap, system->total_cap, msg );
#if defined(DEBUG_FOCUS)
  fprintf( stderr, "p%d: n=%d local_cap=%d\n",
           system->my_rank, system->n, system->local_cap );
  fprintf( stderr, "p%d: N=%d total_cap=%d\n",
           system->my_rank, system->N, system->total_cap );
  fprintf( stderr, "p%d: numH=%d H_cap=%d\n",
           system->my_rank, system->numH, system->Hcap );
#endif

  // if( Reposition_Atoms( system, control, data, mpi_data, msg ) == FAILURE )
  //   return FAILURE;

  /* initialize velocities so that desired init T can be attained */
  if( !control->restart || (control->restart && control->random_vel) )
    Generate_Initial_Velocities( system, control->T_init );

  return SUCCESS;
}


/************************ initialize simulation data ************************/
int Init_Simulation_Data( reax_system *system, control_params *control,
                          simulation_data *data, mpi_datatypes *mpi_data,
                          char *msg )
{
  Reset_Simulation_Data( data, control->virial );

  if( !control->restart )
    data->step = data->prev_steps = 0;

  Compute_Total_Mass( system, data, mpi_data->comm_mesh3D );
  Compute_Center_of_Mass( system, data, mpi_data, mpi_data->comm_mesh3D );
  Compute_Kinetic_Energy( system, data, mpi_data->comm_mesh3D );

  switch( control->ensemble ){
  case NVE:
    data->N_f = 3 * system->bigN;
    Evolve = Velocity_Verlet_NVE;
    break;

  case bNVT:
    data->N_f = 3 * system->bigN + 1;
    Evolve = Velocity_Verlet_Berendsen_NVT;
    break;

  case nhNVT:
    fprintf( stderr, "WARNING: Nose-Hoover NVT is still under testing.\n" );
    //return FAILURE;
    data->N_f = 3 * system->bigN + 1;
    Evolve = Velocity_Verlet_Nose_Hoover_NVT_Klein;
    if( !control->restart || (control->restart && control->random_vel) ) {
      data->therm.G_xi = control->Tau_T *
        (2.0 * data->sys_en.e_kin - data->N_f * K_B * control->T );
      data->therm.v_xi = data->therm.G_xi * control->dt;
      data->therm.v_xi_old = 0;
      data->therm.xi = 0;
    }
    break;

  case sNPT: /* Semi-Isotropic NPT */
    data->N_f = 3 * system->bigN + 4;
    Evolve = Velocity_Verlet_Berendsen_NPT;
    if( !control->restart )
      Reset_Pressures( data );
    break;

  case iNPT: /* Isotropic NPT */
    data->N_f = 3 * system->bigN + 2;
    Evolve = Velocity_Verlet_Berendsen_NPT;
    if( !control->restart )
      Reset_Pressures( data );
    break;

  case NPT: /* Anisotropic NPT */
    strcpy( msg, "init_simulation_data: option not yet implemented" );
    return FAILURE;

    data->N_f = 3 * system->bigN + 9;
    Evolve = Velocity_Verlet_Berendsen_NPT;
    /*if( !control->restart ) {
      data->therm.G_xi = control->Tau_T *
      (2.0 * data->my_en.e_Kin - data->N_f * K_B * control->T );
      data->therm.v_xi = data->therm.G_xi * control->dt;
      data->iso_bar.eps = 0.33333 * log(system->box.volume);
      data->inv_W = 1.0 /
      ( data->N_f * K_B * control->T * SQR(control->Tau_P) );
      Compute_Pressure( system, control, data, out_control );
      }*/
    break;

  default:
    strcpy( msg, "init_simulation_data: ensemble not recognized" );
    return FAILURE;
  }

  /* initialize the timer(s) */
  MPI_Barrier( mpi_data->world );  // wait for everyone to come here
  if( system->my_rank == MASTER_NODE ) {
    data->timing.start = Get_Time( );
#if defined(LOG_PERFORMANCE)
    Reset_Timing( &data->timing );
#endif
  }


#if defined(DEBUG)
      fprintf( stderr, "data->N_f: %8.3f\n", data->N_f );
#endif
  return SUCCESS;
}

#elif defined(LAMMPS_REAX)
int Init_System( reax_system *system, control_params *control, char *msg )
{
  int i;
  reax_atom *atom;

  int mincap = system->mincap;
  double safezone = system->safezone;
  double saferzone = system->saferzone;

  // determine the local and total capacity

  system->local_cap = MAX( (int)(system->n * safezone), mincap);
  system->total_cap = MAX( (int)(system->N * safezone), mincap);

  /* estimate numH and Hcap */
  system->numH = 0;
  if( control->hbond_cut > 0 )
    for( i = 0; i < system->n; ++i ) {
      atom = &(system->my_atoms[i]);
      if( system->reax_param.sbp[ atom->type ].p_hbond == 1 )
        atom->Hindex = system->numH++;
      else atom->Hindex = -1;
    }
  system->Hcap = (int)(MAX( system->numH * saferzone, mincap ));

#if defined(DEBUG_FOCUS)
  fprintf( stderr, "p%d: n=%d local_cap=%d\n",
           system->my_rank, system->n, system->local_cap );
  fprintf( stderr, "p%d: N=%d total_cap=%d\n",
           system->my_rank, system->N, system->total_cap );
  fprintf( stderr, "p%d: numH=%d H_cap=%d\n",
           system->my_rank, system->numH, system->Hcap );
#endif

  return SUCCESS;
}


int Init_Simulation_Data( reax_system *system, control_params *control,
                          simulation_data *data, char *msg )
{
  Reset_Simulation_Data( data, control->virial );

  /* initialize the timer(s) */
  if( system->my_rank == MASTER_NODE ) {
    data->timing.start = Get_Time( );
#if defined(LOG_PERFORMANCE)
    Reset_Timing( &data->timing );
#endif
  }

  //if( !control->restart )
  data->step = data->prev_steps = 0;

  return SUCCESS;
}
#endif



/************************ initialize workspace ************************/
/* Initialize Taper params */
void Init_Taper( control_params *control,  storage *workspace, MPI_Comm comm )
{
  real d1, d7;
  real swa, swa2, swa3;
  real swb, swb2, swb3;

  swa = control->nonb_low;
  swb = control->nonb_cut;

  if( fabs( swa ) > 0.01 )
    fprintf( stderr, "Warning: non-zero lower Taper-radius cutoff\n" );

  if( swb < 0 ) {
    fprintf( stderr, "Negative upper Taper-radius cutoff\n" );
    MPI_Abort( comm,  INVALID_INPUT );
  }
  else if( swb < 5 )
    fprintf( stderr, "Warning: very low Taper-radius cutoff: %f\n", swb );

  d1 = swb - swa;
  d7 = pow( d1, 7.0 );
  swa2 = SQR( swa );
  swa3 = CUBE( swa );
  swb2 = SQR( swb );
  swb3 = CUBE( swb );

  workspace->Tap[7] =  20.0 / d7;
  workspace->Tap[6] = -70.0 * (swa + swb) / d7;
  workspace->Tap[5] =  84.0 * (swa2 + 3.0*swa*swb + swb2) / d7;
  workspace->Tap[4] = -35.0 * (swa3 + 9.0*swa2*swb + 9.0*swa*swb2 + swb3 ) / d7;
  workspace->Tap[3] = 140.0 * (swa3*swb + 3.0*swa2*swb2 + swa*swb3 ) / d7;
  workspace->Tap[2] =-210.0 * (swa3*swb2 + swa2*swb3) / d7;
  workspace->Tap[1] = 140.0 * swa3 * swb3 / d7;
  workspace->Tap[0] = (-35.0*swa3*swb2*swb2 + 21.0*swa2*swb3*swb2 +
                     7.0*swa*swb3*swb3 + swb3*swb3*swb ) / d7;
}


int Init_Workspace( reax_system *system, control_params *control,
                    storage *workspace, MPI_Comm comm, char *msg )
{
  int ret;

  ret = Allocate_Workspace( system, control, workspace,
                            system->local_cap, system->total_cap, comm, msg );
  if( ret != SUCCESS )
    return ret;

  memset( &(workspace->realloc), 0, sizeof(reallocate_data) );
  Reset_Workspace( system, workspace );

  /* Initialize the Taper function */
  Init_Taper( control, workspace, comm );

  return SUCCESS;
}


/************** setup communication data structures  **************/
int Init_MPI_Datatypes( reax_system *system, storage *workspace,
                        mpi_datatypes *mpi_data, MPI_Comm comm, char *msg )
{
#if defined(PURE_REAX)
  int           i, block[11];
  MPI_Aint      base, disp[11];
  MPI_Datatype  type[11];
  mpi_atom      sample;
  boundary_atom b_sample;
  restart_atom  r_sample;
  rvec          rvec_sample;
  rvec2         rvec2_sample;
#endif

  /* setup the world */
  mpi_data->world = comm;
  MPI_Comm_size( comm, &(system->wsize) );

#if defined(PURE_REAX)
  /* init mpi buffers  */
  mpi_data->in1_buffer = NULL;
  mpi_data->in2_buffer = NULL;

  /* mpi_atom - [orig_id, imprt_id, type, num_bonds, num_hbonds, name,
                 x, v, f_old, s, t] */
  block[0] = block[1] = block[2] = block[3] = block[4] = 1;
  block[5] = 8;
  block[6] = block[7] = block[8] = 3;
  block[9] = block[10] = 4;

  MPI_Address( &(sample.orig_id),    disp + 0 );
  MPI_Address( &(sample.imprt_id),   disp + 1 );
  MPI_Address( &(sample.type),       disp + 2 );
  MPI_Address( &(sample.num_bonds),  disp + 3 );
  MPI_Address( &(sample.num_hbonds), disp + 4 );
  MPI_Address( &(sample.name),       disp + 5 );
  MPI_Address( &(sample.x[0]),       disp + 6 );
  MPI_Address( &(sample.v[0]),       disp + 7 );
  MPI_Address( &(sample.f_old[0]),   disp + 8 );
  MPI_Address( &(sample.s[0]),       disp + 9 );
  MPI_Address( &(sample.t[0]),       disp + 10 );

  base = (MPI_Aint)(&(sample));
  for( i = 0; i < 11; ++i ) disp[i] -= base;

  type[0] = type[1] = type[2] = type[3] = type[4] = MPI_INT;
  type[5] = MPI_CHAR;
  type[6] = type[7] = type[8] = type[9] = type[10] = MPI_DOUBLE;

  MPI_Type_struct( 11, block, disp, type, &(mpi_data->mpi_atom_type) );
  MPI_Type_commit( &(mpi_data->mpi_atom_type) );

  /* boundary_atom - [orig_id, imprt_id, type, num_bonds, num_hbonds, x] */
  block[0] = block[1] = block[2] = block[3] = block[4] = 1;
  block[5] = 3;

  MPI_Address( &(b_sample.orig_id),    disp + 0 );
  MPI_Address( &(b_sample.imprt_id),   disp + 1 );
  MPI_Address( &(b_sample.type),       disp + 2 );
  MPI_Address( &(b_sample.num_bonds),  disp + 3 );
  MPI_Address( &(b_sample.num_hbonds), disp + 4 );
  MPI_Address( &(b_sample.x[0]),       disp + 5 );

  base = (MPI_Aint)(&(b_sample));
  for( i = 0; i < 6; ++i ) disp[i] -= base;

  type[0] = type[1] = type[2] = type[3] = type[4] = MPI_INT;
  type[5] = MPI_DOUBLE;

  MPI_Type_struct( 6, block, disp, type, &(mpi_data->boundary_atom_type) );
  MPI_Type_commit( &(mpi_data->boundary_atom_type) );

  /* mpi_rvec */
  block[0] = 3;
  MPI_Address( &(rvec_sample[0]), disp + 0 );
  base = disp[0];
  for( i = 0; i < 1; ++i ) disp[i] -= base;
  type[0] = MPI_DOUBLE;
  MPI_Type_struct( 1, block, disp, type, &(mpi_data->mpi_rvec) );
  MPI_Type_commit( &(mpi_data->mpi_rvec) );

  /* mpi_rvec2 */
  block[0] = 2;
  MPI_Address( &(rvec2_sample[0]), disp + 0 );
  base = disp[0];
  for( i = 0; i < 1; ++i ) disp[i] -= base;
  type[0] = MPI_DOUBLE;
  MPI_Type_struct( 1, block, disp, type, &(mpi_data->mpi_rvec2) );
  MPI_Type_commit( &(mpi_data->mpi_rvec2) );

  /* restart_atom - [orig_id, type, name[8], x, v] */
  block[0] = block[1] = 1 ;
  block[2] = 8;
  block[3] = block[4] = 3;

  MPI_Address( &(r_sample.orig_id),    disp + 0 );
  MPI_Address( &(r_sample.type),       disp + 1 );
  MPI_Address( &(r_sample.name),       disp + 2 );
  MPI_Address( &(r_sample.x[0]),       disp + 3 );
  MPI_Address( &(r_sample.v[0]),       disp + 4 );

  base = (MPI_Aint)(&(r_sample));
  for( i = 0; i < 5; ++i ) disp[i] -= base;

  type[0] = type[1] = MPI_INT;
  type[2] = MPI_CHAR;
  type[3] = type[4] = MPI_DOUBLE;

  MPI_Type_struct( 5, block, disp, type, &(mpi_data->restart_atom_type) );
  MPI_Type_commit( &(mpi_data->restart_atom_type) );
#endif

  return SUCCESS;
}


/********************** allocate lists *************************/
#if defined(PURE_REAX)
int  Init_Lists( reax_system *system, control_params *control,
                 simulation_data *data, storage *workspace, reax_list **lists,
                 mpi_datatypes *mpi_data, char *msg )
{
  int i, num_nbrs;
  int total_hbonds, total_bonds, bond_cap, num_3body, cap_3body, Htop;
  int *hb_top, *bond_top;
  MPI_Comm comm;

  comm = mpi_data->world;
  //for( i = 0; i < MAX_NBRS; ++i ) nrecv[i] = system->my_nbrs[i].est_recv;
  //system->N = SendRecv( system, mpi_data, mpi_data->boundary_atom_type, nrecv,
  //                Sort_Boundary_Atoms, Unpack_Exchange_Message, 1 );
  num_nbrs = Estimate_NumNeighbors( system, lists );
  if(!Make_List( system->total_cap, num_nbrs, TYP_FAR_NEIGHBOR,
                 *lists+FAR_NBRS, comm )){
    fprintf(stderr, "Problem in initializing far nbrs list. Terminating!\n");
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }
#if defined(DEBUG_FOCUS)
  fprintf( stderr, "p%d: allocated far_nbrs: num_far=%d, space=%dMB\n",
           system->my_rank, num_nbrs,
           (int)(num_nbrs*sizeof(far_neighbor_data)/(1024*1024)) );
#endif

  Generate_Neighbor_Lists( system, data, workspace, lists );
  bond_top = (int*) calloc( system->total_cap, sizeof(int) );
  hb_top = (int*) calloc( system->local_cap, sizeof(int) );
  Estimate_Storages( system, control, lists,
                     &Htop, hb_top, bond_top, &num_3body, comm );

  Allocate_Matrix( &(workspace->H), system->local_cap, Htop, comm );
  workspace->L = NULL;
  workspace->U = NULL;
#if defined(DEBUG_FOCUS)
  fprintf( stderr, "p%d: allocated H matrix: Htop=%d, space=%dMB\n",
           system->my_rank, Htop,
           (int)(Htop * sizeof(sparse_matrix_entry) / (1024*1024)) );
#endif

  if( control->hbond_cut > 0 ) {
    /* init H indexes */
    total_hbonds = 0;
    for( i = 0; i < system->n; ++i ) {
      system->my_atoms[i].num_hbonds = hb_top[i];
      total_hbonds += hb_top[i];
    }
    total_hbonds = MAX( total_hbonds*SAFER_ZONE, MIN_CAP*MIN_HBONDS );

    if( !Make_List( system->Hcap, total_hbonds, TYP_HBOND,
                    *lists+HBONDS, comm ) ) {
      fprintf( stderr, "not enough space for hbonds list. terminating!\n" );
      MPI_Abort( comm, INSUFFICIENT_MEMORY );
    }
#if defined(DEBUG_FOCUS)
    fprintf( stderr, "p%d: allocated hbonds: total_hbonds=%d, space=%dMB\n",
             system->my_rank, total_hbonds,
             (int)(total_hbonds*sizeof(hbond_data)/(1024*1024)) );
#endif
  }

  /* bonds list */
  //Allocate_Bond_List( system->N, bond_top, (*lists)+BONDS );
  //num_bonds = bond_top[system->N-1];
  total_bonds = 0;
  for( i = 0; i < system->N; ++i ) {
    system->my_atoms[i].num_bonds = bond_top[i];
    total_bonds += bond_top[i];
  }
  bond_cap = MAX( total_bonds*SAFE_ZONE, MIN_CAP*MIN_BONDS );

  if( !Make_List( system->total_cap, bond_cap, TYP_BOND,
                  *lists+BONDS, comm ) ) {
    fprintf( stderr, "not enough space for bonds list. terminating!\n" );
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }
#if defined(DEBUG_FOCUS)
  fprintf( stderr, "p%d: allocated bonds: total_bonds=%d, space=%dMB\n",
           system->my_rank, bond_cap,
           (int)(bond_cap*sizeof(bond_data)/(1024*1024)) );
#endif

  /* 3bodies list */
  cap_3body = MAX( num_3body*SAFE_ZONE, MIN_3BODIES );
  if( !Make_List( bond_cap, cap_3body, TYP_THREE_BODY,
                  *lists+THREE_BODIES, comm ) ){
    fprintf( stderr, "Problem in initializing angles list. Terminating!\n" );
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }
#if defined(DEBUG_FOCUS)
  fprintf( stderr, "p%d: allocated 3-body list: num_3body=%d, space=%dMB\n",
           system->my_rank, cap_3body,
           (int)(cap_3body*sizeof(three_body_interaction_data)/(1024*1024)) );
#endif

#if defined(TEST_FORCES)
  if( !Make_List( system->total_cap, bond_cap*8, TYP_DDELTA,
                  *lists+DDELTAS, comm ) ) {
    fprintf( stderr, "Problem in initializing dDelta list. Terminating!\n" );
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }
  fprintf( stderr, "p%d: allocated dDelta list: num_ddelta=%d space=%ldMB\n",
           system->my_rank, bond_cap*30,
           bond_cap*8*sizeof(dDelta_data)/(1024*1024) );

  if( !Make_List( bond_cap, bond_cap*50, TYP_DBO, *lists+DBOS, comm ) ) {
    fprintf( stderr, "Problem in initializing dBO list. Terminating!\n" );
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }
  fprintf( stderr, "p%d: allocated dbond list: num_dbonds=%d space=%ldMB\n",
           system->my_rank, bond_cap*MAX_BONDS*3,
           bond_cap*MAX_BONDS*3*sizeof(dbond_data)/(1024*1024) );
#endif

  free( hb_top );
  free( bond_top );

  return SUCCESS;
}
#elif defined(LAMMPS_REAX)
int  Init_Lists( reax_system *system, control_params *control,
                 simulation_data *data, storage *workspace, reax_list **lists,
                 mpi_datatypes *mpi_data, char *msg )
{
  int i, num_nbrs;
  int total_hbonds, total_bonds, bond_cap, num_3body, cap_3body, Htop;
  int *hb_top, *bond_top;
  int nrecv[MAX_NBRS];
  MPI_Comm comm;

  int mincap = system->mincap;
  double safezone = system->safezone;
  double saferzone = system->saferzone;

  comm = mpi_data->world;
  bond_top = (int*) calloc( system->total_cap, sizeof(int) );
  hb_top = (int*) calloc( system->local_cap, sizeof(int) );
  Estimate_Storages( system, control, lists,
                     &Htop, hb_top, bond_top, &num_3body, comm );

  if( control->hbond_cut > 0 ) {
    /* init H indexes */
    total_hbonds = 0;
    for( i = 0; i < system->n; ++i ) {
      system->my_atoms[i].num_hbonds = hb_top[i];
      total_hbonds += hb_top[i];
    }
    total_hbonds = (int)(MAX( total_hbonds*saferzone, mincap*MIN_HBONDS ));

    if( !Make_List( system->Hcap, total_hbonds, TYP_HBOND,
                    *lists+HBONDS, comm ) ) {
      fprintf( stderr, "not enough space for hbonds list. terminating!\n" );
      MPI_Abort( comm, INSUFFICIENT_MEMORY );
    }
#if defined(DEBUG_FOCUS)
    fprintf( stderr, "p%d: allocated hbonds: total_hbonds=%d, space=%dMB\n",
             system->my_rank, total_hbonds,
             (int)(total_hbonds*sizeof(hbond_data)/(1024*1024)) );
#endif
  }

  /* bonds list */
  //Allocate_Bond_List( system->N, bond_top, (*lists)+BONDS );
  //num_bonds = bond_top[system->N-1];
  total_bonds = 0;
  for( i = 0; i < system->N; ++i ) {
    system->my_atoms[i].num_bonds = bond_top[i];
    total_bonds += bond_top[i];
  }
  bond_cap = (int)(MAX( total_bonds*safezone, mincap*MIN_BONDS ));

  if( !Make_List( system->total_cap, bond_cap, TYP_BOND,
                  *lists+BONDS, comm ) ) {
    fprintf( stderr, "not enough space for bonds list. terminating!\n" );
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }
#if defined(DEBUG_FOCUS)
  fprintf( stderr, "p%d: allocated bonds: total_bonds=%d, space=%dMB\n",
           system->my_rank, bond_cap,
           (int)(bond_cap*sizeof(bond_data)/(1024*1024)) );
#endif

  /* 3bodies list */
  cap_3body = (int)(MAX( num_3body*safezone, MIN_3BODIES ));
  if( !Make_List( bond_cap, cap_3body, TYP_THREE_BODY,
                  *lists+THREE_BODIES, comm ) ){
    fprintf( stderr, "Problem in initializing angles list. Terminating!\n" );
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }
#if defined(DEBUG_FOCUS)
  fprintf( stderr, "p%d: allocated 3-body list: num_3body=%d, space=%dMB\n",
           system->my_rank, cap_3body,
           (int)(cap_3body*sizeof(three_body_interaction_data)/(1024*1024)) );
#endif

#if defined(TEST_FORCES)
  if( !Make_List( system->total_cap, bond_cap*8, TYP_DDELTA,
                  *lists+DDELTAS, comm ) ) {
    fprintf( stderr, "Problem in initializing dDelta list. Terminating!\n" );
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }
  fprintf( stderr, "p%d: allocated dDelta list: num_ddelta=%d space=%ldMB\n",
           system->my_rank, bond_cap*30,
           bond_cap*8*sizeof(dDelta_data)/(1024*1024) );

  if( !Make_List( bond_cap, bond_cap*50, TYP_DBO, (*lists)+DBOS, comm ) ) {
    fprintf( stderr, "Problem in initializing dBO list. Terminating!\n" );
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }
  fprintf( stderr, "p%d: allocated dbond list: num_dbonds=%d space=%ldMB\n",
           system->my_rank, bond_cap*MAX_BONDS*3,
           bond_cap*MAX_BONDS*3*sizeof(dbond_data)/(1024*1024) );
#endif

  free( hb_top );
  free( bond_top );

  return SUCCESS;
}
#endif



#if defined(PURE_REAX)
void Initialize( reax_system *system, control_params *control,
                 simulation_data *data, storage *workspace,
                 reax_list **lists, output_controls *out_control,
                 mpi_datatypes *mpi_data )
{
  char msg[MAX_STR];

  if( Init_MPI_Datatypes( system, workspace, mpi_data, MPI_COMM_WORLD, msg ) ==
      FAILURE ) {
    fprintf( stderr, "p%d: init_mpi_datatypes: could not create datatypes\n",
             system->my_rank );
    fprintf( stderr, "p%d: mpi_data couldn't be initialized! terminating.\n",
             system->my_rank );
    MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
  }
#if defined(DEBUG)
  fprintf( stderr, "p%d: initialized mpi datatypes\n", system->my_rank );
#endif

  if( Init_System(system, control, data, workspace, mpi_data, msg) == FAILURE ){
    fprintf( stderr, "p%d: %s\n", system->my_rank, msg );
    fprintf( stderr, "p%d: system could not be initialized! terminating.\n",
             system->my_rank );
    MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
  }
#if defined(DEBUG)
  fprintf( stderr, "p%d: system initialized\n", system->my_rank );
#endif

  if( Init_Simulation_Data(system, control, data, mpi_data, msg) == FAILURE ) {
    fprintf( stderr, "p%d: %s\n", system->my_rank, msg );
    fprintf( stderr, "p%d: sim_data couldn't be initialized! terminating.\n",
             system->my_rank );
    MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
  }
#if defined(DEBUG)
  fprintf( stderr, "p%d: initialized simulation data\n", system->my_rank );
#endif

  if( Init_Workspace( system, control, workspace, mpi_data->world, msg ) ==
      FAILURE ) {
    fprintf( stderr, "p%d:init_workspace: not enough memory\n",
             system->my_rank );
    fprintf( stderr, "p%d:workspace couldn't be initialized! terminating.\n",
             system->my_rank );
    MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
  }
#if defined(DEBUG)
  fprintf( stderr, "p%d: initialized workspace\n", system->my_rank );
#endif

  if( Init_Lists( system, control, data, workspace, lists, mpi_data, msg ) ==
      FAILURE ) {
      fprintf( stderr, "p%d: %s\n", system->my_rank, msg );
      fprintf( stderr, "p%d: system could not be initialized! terminating.\n",
               system->my_rank );
      MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
    }
#if defined(DEBUG)
  fprintf( stderr, "p%d: initialized lists\n", system->my_rank );
#endif

  if(Init_Output_Files(system,control,out_control,mpi_data,msg) == FAILURE) {
    fprintf( stderr, "p%d: %s\n", system->my_rank, msg );
    fprintf( stderr, "p%d: could not open output files! terminating...\n",
             system->my_rank );
    MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
  }
#if defined(DEBUG)
  fprintf( stderr, "p%d: output files opened\n", system->my_rank );
#endif

  if( control->tabulate ) {
    if( Init_Lookup_Tables(system,control,workspace,mpi_data,msg) == FAILURE ) {
      fprintf( stderr, "p%d: %s\n", system->my_rank, msg );
      fprintf( stderr, "p%d: couldn't create lookup table! terminating.\n",
               system->my_rank );
      MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
    }
#if defined(DEBUG)
    fprintf( stderr, "p%d: initialized lookup tables\n", system->my_rank );
#endif
  }

  Init_Force_Functions( control );
#if defined(DEBUG)
  fprintf( stderr, "p%d: initialized force functions\n", system->my_rank );
#endif
  /*#ifdef TEST_FORCES
    Init_Force_Test_Functions();
    fprintf(stderr,"p%d: initialized force test functions\n",system->my_rank);
    #endif */
}

#elif defined(LAMMPS_REAX)
void Initialize( reax_system *system, control_params *control,
                 simulation_data *data, storage *workspace,
                 reax_list **lists, output_controls *out_control,
                 mpi_datatypes *mpi_data, MPI_Comm comm )
{
  char msg[MAX_STR];


  if( Init_MPI_Datatypes(system, workspace, mpi_data, comm, msg) == FAILURE ) {
    fprintf( stderr, "p%d: init_mpi_datatypes: could not create datatypes\n",
             system->my_rank );
    fprintf( stderr, "p%d: mpi_data couldn't be initialized! terminating.\n",
             system->my_rank );
    MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
  }
#if defined(DEBUG)
  fprintf( stderr, "p%d: initialized mpi datatypes\n", system->my_rank );
#endif

  if( Init_System(system, control, msg) == FAILURE ){
    fprintf( stderr, "p%d: %s\n", system->my_rank, msg );
    fprintf( stderr, "p%d: system could not be initialized! terminating.\n",
             system->my_rank );
    MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
  }
#if defined(DEBUG)
  fprintf( stderr, "p%d: system initialized\n", system->my_rank );
#endif

  if( Init_Simulation_Data( system, control, data, msg ) == FAILURE ) {
    fprintf( stderr, "p%d: %s\n", system->my_rank, msg );
    fprintf( stderr, "p%d: sim_data couldn't be initialized! terminating.\n",
             system->my_rank );
    MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
  }
#if defined(DEBUG)
  fprintf( stderr, "p%d: initialized simulation data\n", system->my_rank );
#endif

  if( Init_Workspace( system, control, workspace, mpi_data->world, msg ) ==
      FAILURE ) {
    fprintf( stderr, "p%d:init_workspace: not enough memory\n",
             system->my_rank );
    fprintf( stderr, "p%d:workspace couldn't be initialized! terminating.\n",
             system->my_rank );
    MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
  }
#if defined(DEBUG)
  fprintf( stderr, "p%d: initialized workspace\n", system->my_rank );
#endif

  if( Init_Lists( system, control, data, workspace, lists, mpi_data, msg ) ==
      FAILURE ) {
      fprintf( stderr, "p%d: %s\n", system->my_rank, msg );
      fprintf( stderr, "p%d: system could not be initialized! terminating.\n",
               system->my_rank );
      MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
    }
#if defined(DEBUG)
  fprintf( stderr, "p%d: initialized lists\n", system->my_rank );
#endif

  if( Init_Output_Files(system,control,out_control,mpi_data,msg)== FAILURE) {
    fprintf( stderr, "p%d: %s\n", system->my_rank, msg );
    fprintf( stderr, "p%d: could not open output files! terminating...\n",
             system->my_rank );
    MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
  }
#if defined(DEBUG)
  fprintf( stderr, "p%d: output files opened\n", system->my_rank );
#endif

  if( control->tabulate ) {
    if( Init_Lookup_Tables( system, control, workspace, mpi_data, msg ) == FAILURE ) {
      fprintf( stderr, "p%d: %s\n", system->my_rank, msg );
      fprintf( stderr, "p%d: couldn't create lookup table! terminating.\n",
               system->my_rank );
      MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
    }
#if defined(DEBUG)
    fprintf( stderr, "p%d: initialized lookup tables\n", system->my_rank );
#endif
  }


  Init_Force_Functions( control );
#if defined(DEBUG)
  fprintf( stderr, "p%d: initialized force functions\n", system->my_rank );
#endif
  /*#if defined(TEST_FORCES)
    Init_Force_Test_Functions();
    fprintf(stderr,"p%d: initialized force test functions\n",system->my_rank);
  #endif*/
}
#endif
