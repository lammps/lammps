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
#include "forces.h"
#include "bond_orders.h"
#include "bonds.h"
#include "basic_comm.h"
#include "hydrogen_bonds.h"
#include "io_tools.h"
#include "list.h"
#include "lookup.h"
#include "multi_body.h"
#include "nonbonded.h"
#include "qEq.h"
#include "tool_box.h"
#include "torsion_angles.h"
#include "valence_angles.h"
#include "vector.h"
#elif defined(LAMMPS_REAX)
#include "reaxc_forces.h"
#include "reaxc_bond_orders.h"
#include "reaxc_bonds.h"
#include "reaxc_basic_comm.h"
#include "reaxc_hydrogen_bonds.h"
#include "reaxc_io_tools.h"
#include "reaxc_list.h"
#include "reaxc_lookup.h"
#include "reaxc_multi_body.h"
#include "reaxc_nonbonded.h"
#include "reaxc_tool_box.h"
#include "reaxc_torsion_angles.h"
#include "reaxc_valence_angles.h"
#include "reaxc_vector.h"
#endif

interaction_function Interaction_Functions[NUM_INTRS];

void Dummy_Interaction( reax_system *system, control_params *control,
                        simulation_data *data, storage *workspace,
                        reax_list **lists, output_controls *out_control )
{
}


void Init_Force_Functions( control_params *control )
{
  Interaction_Functions[0] = BO;
  Interaction_Functions[1] = Bonds; //Dummy_Interaction;
  Interaction_Functions[2] = Atom_Energy; //Dummy_Interaction;
  Interaction_Functions[3] = Valence_Angles; //Dummy_Interaction;
  Interaction_Functions[4] = Torsion_Angles; //Dummy_Interaction;
  if( control->hbond_cut > 0 )
    Interaction_Functions[5] = Hydrogen_Bonds;
  else Interaction_Functions[5] = Dummy_Interaction;
  Interaction_Functions[6] = Dummy_Interaction; //empty
  Interaction_Functions[7] = Dummy_Interaction; //empty
  Interaction_Functions[8] = Dummy_Interaction; //empty
  Interaction_Functions[9] = Dummy_Interaction; //empty
}


void Compute_Bonded_Forces( reax_system *system, control_params *control,
                            simulation_data *data, storage *workspace,
                            reax_list **lists, output_controls *out_control,
                            MPI_Comm comm )
{
  int i;

  /* Mark beginning of a new timestep in bonded energy files */
#if defined(TEST_ENERGY)
  Debug_Marker_Bonded( out_control, data->step );
#endif

  /* Implement all force calls as function pointers */
  for( i = 0; i < NUM_INTRS; i++ ) {
#if defined(DEBUG)
    fprintf( stderr, "p%d: starting f%d\n", system->my_rank, i );
    MPI_Barrier( comm );
#endif
    (Interaction_Functions[i])( system, control, data, workspace,
                                lists, out_control );
#if defined(DEBUG)
    fprintf( stderr, "p%d: f%d done\n", system->my_rank, i );
    MPI_Barrier( comm );
#endif
  }
}


void Compute_NonBonded_Forces( reax_system *system, control_params *control,
                               simulation_data *data, storage *workspace,
                               reax_list **lists, output_controls *out_control,
                               MPI_Comm comm )
{
  /* Mark beginning of a new timestep in nonbonded energy files */
#if defined(TEST_ENERGY)
  Debug_Marker_Nonbonded( out_control, data->step );
#endif

  /* van der Waals and Coulomb interactions */
  if( control->tabulate == 0 )
    vdW_Coulomb_Energy( system, control, data, workspace,
                        lists, out_control );
  else
    Tabulated_vdW_Coulomb_Energy( system, control, data, workspace,
                                  lists, out_control );

#if defined(DEBUG)
  fprintf( stderr, "p%d: nonbonded forces done\n", system->my_rank );
  MPI_Barrier( comm );
#endif
}



/* this version of Compute_Total_Force computes forces from
   coefficients accumulated by all interaction functions.
   Saves enormous time & space! */
void Compute_Total_Force( reax_system *system, control_params *control,
                          simulation_data *data, storage *workspace,
                          reax_list **lists, mpi_datatypes *mpi_data )
{
  int i, pj;
  reax_list *bonds = (*lists) + BONDS;

  for( i = 0; i < system->N; ++i )
    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj )
      if( i < bonds->select.bond_list[pj].nbr ) {
        if( control->virial == 0 )
          Add_dBond_to_Forces( system, i, pj, workspace, lists );
        else
          Add_dBond_to_Forces_NPT( i, pj, data, workspace, lists );
      }

  //Print_Total_Force( system, data, workspace );
#if defined(PURE_REAX)
  /* now all forces are computed to their partially-final values
     based on the neighbors information each processor has had.
     final values of force on each atom needs to be computed by adding up
     all partially-final pieces */
  Coll( system, mpi_data, workspace->f, mpi_data->mpi_rvec,
        sizeof(rvec)/sizeof(void), rvec_unpacker );
  for( i = 0; i < system->n; ++i )
    rvec_Copy( system->my_atoms[i].f, workspace->f[i] );

#if defined(TEST_FORCES)
  Coll( system, mpi_data, workspace->f_ele, mpi_data->mpi_rvec, rvec_unpacker);
  Coll( system, mpi_data, workspace->f_vdw, mpi_data->mpi_rvec, rvec_unpacker);
  Coll( system, mpi_data, workspace->f_be, mpi_data->mpi_rvec, rvec_unpacker );
  Coll( system, mpi_data, workspace->f_lp, mpi_data->mpi_rvec, rvec_unpacker );
  Coll( system, mpi_data, workspace->f_ov, mpi_data->mpi_rvec, rvec_unpacker );
  Coll( system, mpi_data, workspace->f_un, mpi_data->mpi_rvec, rvec_unpacker );
  Coll( system, mpi_data, workspace->f_ang, mpi_data->mpi_rvec, rvec_unpacker);
  Coll( system, mpi_data, workspace->f_coa, mpi_data->mpi_rvec, rvec_unpacker);
  Coll( system, mpi_data, workspace->f_pen, mpi_data->mpi_rvec, rvec_unpacker);
  Coll( system, mpi_data, workspace->f_hb, mpi_data->mpi_rvec, rvec_unpacker );
  Coll( system, mpi_data, workspace->f_tor, mpi_data->mpi_rvec, rvec_unpacker);
  Coll( system, mpi_data, workspace->f_con, mpi_data->mpi_rvec, rvec_unpacker);
#endif

#endif
}

void Validate_Lists( reax_system *system, storage *workspace, reax_list **lists,
                     int step, int n, int N, int numH, MPI_Comm comm )
{
  int i, comp, Hindex;
  reax_list *bonds, *hbonds;
  reallocate_data *realloc;
  realloc = &(workspace->realloc);

  double saferzone = system->saferzone;

  /* bond list */
  if( N > 0 ) {
    bonds = *lists + BONDS;

    for( i = 0; i < N; ++i ) {
      // if( i < n ) - we need to update ghost estimates for delayed nbrings
      system->my_atoms[i].num_bonds = MAX(Num_Entries(i,bonds)*2, MIN_BONDS);

      //if( End_Index(i, bonds) >= Start_Index(i+1, bonds)-2 )
      //workspace->realloc.bonds = 1;

      if( i < N-1 )
        comp = Start_Index(i+1, bonds);
      else comp = bonds->num_intrs;

      if( End_Index(i, bonds) > comp ) {
        fprintf( stderr, "step%d-bondchk failed: i=%d end(i)=%d str(i+1)=%d\n",
                 step, i, End_Index(i,bonds), comp );
        MPI_Abort( comm, INSUFFICIENT_MEMORY );
      }
    }
  }


  /* hbonds list */
  if( numH > 0 ) {
    hbonds = *lists + HBONDS;

    for( i = 0; i < n; ++i ) {
      Hindex = system->my_atoms[i].Hindex;
      if( Hindex > -1 ) {
        system->my_atoms[i].num_hbonds =
          (int)(MAX( Num_Entries(Hindex, hbonds)*saferzone, MIN_HBONDS ));

        //if( Num_Entries(i, hbonds) >=
        //(Start_Index(i+1,hbonds)-Start_Index(i,hbonds))*0.90/*DANGER_ZONE*/){
        //  workspace->realloc.hbonds = 1;

        if( Hindex < numH-1 )
          comp = Start_Index(Hindex+1, hbonds);
        else comp = hbonds->num_intrs;

        if( End_Index(Hindex, hbonds) > comp ) {
          fprintf(stderr,"step%d-hbondchk failed: H=%d end(H)=%d str(H+1)=%d\n",
                  step, Hindex, End_Index(Hindex,hbonds), comp );
          MPI_Abort( comm, INSUFFICIENT_MEMORY );
        }
      }
    }
  }
}


real Compute_H( real r, real gamma, real *ctap )
{
  real taper, dr3gamij_1, dr3gamij_3;

  taper = ctap[7] * r + ctap[6];
  taper = taper * r + ctap[5];
  taper = taper * r + ctap[4];
  taper = taper * r + ctap[3];
  taper = taper * r + ctap[2];
  taper = taper * r + ctap[1];
  taper = taper * r + ctap[0];

  dr3gamij_1 = ( r*r*r + gamma );
  dr3gamij_3 = pow( dr3gamij_1 , 0.33333333333333 );
  return taper * EV_to_KCALpMOL / dr3gamij_3;
}


real Compute_tabH( real r_ij, int ti, int tj )
{
  int r, tmin, tmax;
  real val, dif, base;
  LR_lookup_table *t;

  tmin  = MIN( ti, tj );
  tmax  = MAX( ti, tj );
  t = &( LR[tmin][tmax] );

  /* cubic spline interpolation */
  r = (int)(r_ij * t->inv_dx);
  if( r == 0 )  ++r;
  base = (real)(r+1) * t->dx;
  dif = r_ij - base;
  val = ((t->ele[r].d*dif + t->ele[r].c)*dif + t->ele[r].b)*dif +
    t->ele[r].a;
  val *= EV_to_KCALpMOL / C_ele;

  return val;
}


void Init_Forces( reax_system *system, control_params *control,
                  simulation_data *data, storage *workspace, reax_list **lists,
                  output_controls *out_control, MPI_Comm comm ) {
  int i, j, pj;
  int start_i, end_i;
  int type_i, type_j;
  int Htop, btop_i, btop_j, num_bonds, num_hbonds;
  int ihb, jhb, ihb_top, jhb_top;
  int local, flag, renbr;
  real r_ij, cutoff;
  sparse_matrix *H;
  reax_list *far_nbrs, *bonds, *hbonds;
  single_body_parameters *sbp_i, *sbp_j;
  two_body_parameters *twbp;
  far_neighbor_data *nbr_pj;
  reax_atom *atom_i, *atom_j;

  far_nbrs = *lists + FAR_NBRS;
  bonds = *lists + BONDS;
  hbonds = *lists + HBONDS;

  for( i = 0; i < system->n; ++i )
    workspace->bond_mark[i] = 0;
  for( i = system->n; i < system->N; ++i ) {
    workspace->bond_mark[i] = 1000; // put ghost atoms to an infinite distance
    //workspace->done_after[i] = Start_Index( i, far_nbrs );
  }

  H = workspace->H;
  H->n = system->n;
  Htop = 0;
  num_bonds = 0;
  num_hbonds = 0;
  btop_i = btop_j = 0;
  renbr = (data->step-data->prev_steps) % control->reneighbor == 0;

  for( i = 0; i < system->N; ++i ) {
    atom_i = &(system->my_atoms[i]);
    type_i  = atom_i->type;
    start_i = Start_Index(i, far_nbrs);
    end_i   = End_Index(i, far_nbrs);
    btop_i = End_Index( i, bonds );
    sbp_i = &(system->reax_param.sbp[type_i]);

    if( i < system->n ) {
      local = 1;
      cutoff = control->nonb_cut;
    }
    else {
      local = 0;
      cutoff = control->bond_cut;
    }

    ihb = -1;
    ihb_top = -1;
    if( local ) {
      H->start[i] = Htop;
      H->entries[Htop].j = i;
      H->entries[Htop].val = sbp_i->eta;
      ++Htop;

      if( control->hbond_cut > 0 ) {
        ihb = sbp_i->p_hbond;
        if( ihb == 1 )
          ihb_top = End_Index( atom_i->Hindex, hbonds );
        else ihb_top = -1;
      }
    }

    /* update i-j distance - check if j is within cutoff */
    for( pj = start_i; pj < end_i; ++pj ) {
      nbr_pj = &( far_nbrs->select.far_nbr_list[pj] );
      j = nbr_pj->nbr;
      atom_j = &(system->my_atoms[j]);
      //fprintf( stderr, "%d%d i=%d x_i: %f %f %f,j=%d x_j: %f %f %f, d=%f\n",
      //         MIN(atom_i->orig_id, atom_j->orig_id),
      //         MAX(atom_i->orig_id, atom_j->orig_id),
      //         i, atom_i->x[0], atom_i->x[1], atom_i->x[2],
      //         j, atom_j->x[0], atom_j->x[1], atom_j->x[2], nbr_pj->d );
      if( renbr ) {
        if(nbr_pj->d <= cutoff)
          flag = 1;
        else flag = 0;
      }
      else{
        nbr_pj->dvec[0] = atom_j->x[0] - atom_i->x[0];
        nbr_pj->dvec[1] = atom_j->x[1] - atom_i->x[1];
        nbr_pj->dvec[2] = atom_j->x[2] - atom_i->x[2];
        nbr_pj->d = rvec_Norm_Sqr( nbr_pj->dvec );
        if( nbr_pj->d <= SQR(cutoff) ) {
          nbr_pj->d = sqrt(nbr_pj->d);
          flag = 1;
        }
        else {
          flag = 0;
        }
      }

      if( flag ){
        type_j = atom_j->type;
        r_ij = nbr_pj->d;
        sbp_j = &(system->reax_param.sbp[type_j]);
        twbp = &(system->reax_param.tbp[type_i][type_j]);

        if( local ) {
          /* H matrix entry */
          if( j < system->n || atom_i->orig_id < atom_j->orig_id ) {//tryQEq||1
            H->entries[Htop].j = j;
            //fprintf( stdout, "%d%d %d %d\n",
            //     MIN(atom_i->orig_id, atom_j->orig_id),
            //     MAX(atom_i->orig_id, atom_j->orig_id),
            //     MIN(atom_i->orig_id, atom_j->orig_id),
            //     MAX(atom_i->orig_id, atom_j->orig_id) );
            if( control->tabulate == 0 )
              H->entries[Htop].val = Compute_H(r_ij,twbp->gamma,workspace->Tap);
            else H->entries[Htop].val = Compute_tabH(r_ij, type_i, type_j);
            ++Htop;
          }

          /* hydrogen bond lists */
          if( control->hbond_cut > 0 && (ihb==1 || ihb==2) &&
              nbr_pj->d <= control->hbond_cut ) {
            // fprintf( stderr, "%d %d\n", atom1, atom2 );
            jhb = sbp_j->p_hbond;
            if( ihb == 1 && jhb == 2 ) {
              hbonds->select.hbond_list[ihb_top].nbr = j;
              hbonds->select.hbond_list[ihb_top].scl = 1;
              hbonds->select.hbond_list[ihb_top].ptr = nbr_pj;
              ++ihb_top;
              ++num_hbonds;
            }
            else if( j < system->n && ihb == 2 && jhb == 1 ) {
              jhb_top = End_Index( atom_j->Hindex, hbonds );
              hbonds->select.hbond_list[jhb_top].nbr = i;
              hbonds->select.hbond_list[jhb_top].scl = -1;
              hbonds->select.hbond_list[jhb_top].ptr = nbr_pj;
              Set_End_Index( atom_j->Hindex, jhb_top+1, hbonds );
              ++num_hbonds;
            }
          }
        }

        /* uncorrected bond orders */
        if( //(workspace->bond_mark[i] < 3 || workspace->bond_mark[j] < 3) &&
            nbr_pj->d <= control->bond_cut &&
            BOp( workspace, bonds, control->bo_cut,
                 i , btop_i, nbr_pj, sbp_i, sbp_j, twbp ) ) {
          num_bonds += 2;
          ++btop_i;

          if( workspace->bond_mark[j] > workspace->bond_mark[i] + 1 )
            workspace->bond_mark[j] = workspace->bond_mark[i] + 1;
          else if( workspace->bond_mark[i] > workspace->bond_mark[j] + 1 ) {
            workspace->bond_mark[i] = workspace->bond_mark[j] + 1;
            //if( workspace->bond_mark[i] == 1000 )
            //  workspace->done_after[i] = pj;
          }
          //fprintf( stdout, "%d%d - %d(%d) %d(%d)\n",
          //   i , j, i, workspace->bond_mark[i], j, workspace->bond_mark[j] );
        }
      }
    }

    Set_End_Index( i, btop_i, bonds );
    if( local ) {
      H->end[i] = Htop;
      if( ihb == 1 )
        Set_End_Index( atom_i->Hindex, ihb_top, hbonds );
    }
  }

  //fprintf( stderr, "after the first init loop\n" );
  /*for( i = system->n; i < system->N; ++i )
    if( workspace->bond_mark[i] > 3 ) {
      start_i = Start_Index(i, bonds);
      end_i = End_Index(i, bonds);
      num_bonds -= (end_i - start_i);
      Set_End_Index(i, start_i, bonds );
      }*/

  /*for( i = system->n; i < system->N; ++i ) {
    start_i = Start_Index(i, far_nbrs);
    end_i = workspace->done_after[i];

    if( workspace->bond_mark[i] >= 2 && start_i < end_i ) {
      atom_i = &(system->my_atoms[i]);
      type_i = atom_i->type;
      btop_i = End_Index( i, bonds );
      sbp_i = &(system->reax_param.sbp[type_i]);

      for( pj = start_i; pj < end_i; ++pj ) {
        nbr_pj = &( far_nbrs->select.far_nbr_list[pj] );
        j = nbr_pj->nbr;

        if( workspace->bond_mark[j] >= 2 && nbr_pj->d <= control->bond_cut ) {
          atom_j = &(system->my_atoms[j]);
          type_j = atom_j->type;
          sbp_j = &(system->reax_param.sbp[type_j]);
          twbp = &(system->reax_param.tbp[type_i][type_j]);

          if( BOp( workspace, bonds, control->bo_cut,
                   i , btop_i, nbr_pj, sbp_i, sbp_j, twbp ) ) {
            num_bonds += 2;
            ++btop_i;

            if( workspace->bond_mark[j] > workspace->bond_mark[i] + 1 )
              workspace->bond_mark[j] = workspace->bond_mark[i] + 1;
            else if( workspace->bond_mark[i] > workspace->bond_mark[j] + 1 )
              workspace->bond_mark[i] = workspace->bond_mark[j] + 1;

            //fprintf( stdout, "%d%d - %d(%d) %d(%d) new\n",
            // i , j, i, workspace->bond_mark[i], j, workspace->bond_mark[j] );
          }
        }
      }
      Set_End_Index( i, btop_i, bonds );
    }
    }*/

  workspace->realloc.Htop = Htop;
  workspace->realloc.num_bonds = num_bonds;
  workspace->realloc.num_hbonds = num_hbonds;

#if defined(DEBUG_FOCUS)
  fprintf( stderr, "p%d @ step%d: Htop = %d num_bonds = %d num_hbonds = %d\n",
           system->my_rank, data->step, Htop, num_bonds, num_hbonds );
  MPI_Barrier( comm );
#endif
#if defined( DEBUG )
  Print_Bonds( system, bonds, "debugbonds.out" );
  Print_Bond_List2( system, bonds, "pbonds.out" );
  Print_Sparse_Matrix( system, H );
  for( i = 0; i < H->n; ++i )
    for( j = H->start[i]; j < H->end[i]; ++j )
      fprintf( stderr, "%d %d %.15e\n",
               MIN(system->my_atoms[i].orig_id,
                   system->my_atoms[H->entries[j].j].orig_id),
               MAX(system->my_atoms[i].orig_id,
                   system->my_atoms[H->entries[j].j].orig_id),
               H->entries[j].val );
#endif

  Validate_Lists( system, workspace, lists, data->step,
                  system->n, system->N, system->numH, comm );
}


void Init_Forces_noQEq( reax_system *system, control_params *control,
                        simulation_data *data, storage *workspace,
                        reax_list **lists, output_controls *out_control,
                        MPI_Comm comm ) {
  int i, j, pj;
  int start_i, end_i;
  int type_i, type_j;
  int btop_i, btop_j, num_bonds, num_hbonds;
  int ihb, jhb, ihb_top, jhb_top;
  int local, flag, renbr;
  real r_ij, cutoff;
  reax_list *far_nbrs, *bonds, *hbonds;
  single_body_parameters *sbp_i, *sbp_j;
  two_body_parameters *twbp;
  far_neighbor_data *nbr_pj;
  reax_atom *atom_i, *atom_j;

  far_nbrs = *lists + FAR_NBRS;
  bonds = *lists + BONDS;
  hbonds = *lists + HBONDS;

  for( i = 0; i < system->n; ++i )
    workspace->bond_mark[i] = 0;
  for( i = system->n; i < system->N; ++i ) {
    workspace->bond_mark[i] = 1000; // put ghost atoms to an infinite distance
    //workspace->done_after[i] = Start_Index( i, far_nbrs );
  }

  num_bonds = 0;
  num_hbonds = 0;
  btop_i = btop_j = 0;
  renbr = (data->step-data->prev_steps) % control->reneighbor == 0;

  for( i = 0; i < system->N; ++i ) {
    atom_i = &(system->my_atoms[i]);
    type_i  = atom_i->type;
    start_i = Start_Index(i, far_nbrs);
    end_i   = End_Index(i, far_nbrs);
    btop_i = End_Index( i, bonds );
    sbp_i = &(system->reax_param.sbp[type_i]);

    if( i < system->n ) {
      local = 1;
      cutoff = MAX( control->hbond_cut, control->bond_cut );
    }
    else {
      local = 0;
      cutoff = control->bond_cut;
    }

    ihb = -1;
    ihb_top = -1;
    if( local && control->hbond_cut > 0 ) {
      ihb = sbp_i->p_hbond;
      if( ihb == 1 )
        ihb_top = End_Index( atom_i->Hindex, hbonds );
      else ihb_top = -1;
    }

    /* update i-j distance - check if j is within cutoff */
    for( pj = start_i; pj < end_i; ++pj ) {
      nbr_pj = &( far_nbrs->select.far_nbr_list[pj] );
      j = nbr_pj->nbr;
      atom_j = &(system->my_atoms[j]);

      if( renbr ) {
        if( nbr_pj->d <= cutoff )
          flag = 1;
        else flag = 0;
      }
      else{
        nbr_pj->dvec[0] = atom_j->x[0] - atom_i->x[0];
        nbr_pj->dvec[1] = atom_j->x[1] - atom_i->x[1];
        nbr_pj->dvec[2] = atom_j->x[2] - atom_i->x[2];
        nbr_pj->d = rvec_Norm_Sqr( nbr_pj->dvec );
        if( nbr_pj->d <= SQR(cutoff) ) {
          nbr_pj->d = sqrt(nbr_pj->d);
          flag = 1;
        }
        else {
          flag = 0;
        }
      }

      if( flag ) {
        type_j = atom_j->type;
        r_ij = nbr_pj->d;
        sbp_j = &(system->reax_param.sbp[type_j]);
        twbp = &(system->reax_param.tbp[type_i][type_j]);

        if( local ) {
          /* hydrogen bond lists */
          if( control->hbond_cut > 0 && (ihb==1 || ihb==2) &&
              nbr_pj->d <= control->hbond_cut ) {
            // fprintf( stderr, "%d %d\n", atom1, atom2 );
            jhb = sbp_j->p_hbond;
            if( ihb == 1 && jhb == 2 ) {
              hbonds->select.hbond_list[ihb_top].nbr = j;
              hbonds->select.hbond_list[ihb_top].scl = 1;
              hbonds->select.hbond_list[ihb_top].ptr = nbr_pj;
              ++ihb_top;
              ++num_hbonds;
            }
            else if( j < system->n && ihb == 2 && jhb == 1 ) {
              jhb_top = End_Index( atom_j->Hindex, hbonds );
              hbonds->select.hbond_list[jhb_top].nbr = i;
              hbonds->select.hbond_list[jhb_top].scl = -1;
              hbonds->select.hbond_list[jhb_top].ptr = nbr_pj;
              Set_End_Index( atom_j->Hindex, jhb_top+1, hbonds );
              ++num_hbonds;
            }
          }
        }


        /* uncorrected bond orders */
        if( //(workspace->bond_mark[i] < 3 || workspace->bond_mark[j] < 3) &&
            nbr_pj->d <= control->bond_cut &&
            BOp( workspace, bonds, control->bo_cut,
                 i , btop_i, nbr_pj, sbp_i, sbp_j, twbp ) ) {
          num_bonds += 2;
          ++btop_i;

          if( workspace->bond_mark[j] > workspace->bond_mark[i] + 1 )
            workspace->bond_mark[j] = workspace->bond_mark[i] + 1;
          else if( workspace->bond_mark[i] > workspace->bond_mark[j] + 1 ) {
            workspace->bond_mark[i] = workspace->bond_mark[j] + 1;
            //if( workspace->bond_mark[i] == 1000 )
            //  workspace->done_after[i] = pj;
          }
          //fprintf( stdout, "%d%d - %d(%d) %d(%d)\n",
          //   i , j, i, workspace->bond_mark[i], j, workspace->bond_mark[j] );
        }
      }
    }

    Set_End_Index( i, btop_i, bonds );
    if( local && ihb == 1 )
      Set_End_Index( atom_i->Hindex, ihb_top, hbonds );
  }

  /*for( i = system->n; i < system->N; ++i )
    if( workspace->bond_mark[i] > 3 ) {
      start_i = Start_Index(i, bonds);
      end_i = End_Index(i, bonds);
      num_bonds -= (end_i - start_i);
      Set_End_Index(i, start_i, bonds );
      }*/

  workspace->realloc.num_bonds = num_bonds;
  workspace->realloc.num_hbonds = num_hbonds;

#if defined(DEBUG_FOCUS)
  fprintf( stderr, "p%d @ step%d: num_bonds = %d num_hbonds = %d\n",
           system->my_rank, data->step, num_bonds, num_hbonds );
  MPI_Barrier( comm );
#endif
#if defined( DEBUG )
  Print_Bonds( system, bonds, "debugbonds.out" );
  Print_Bond_List2( system, bonds, "pbonds.out" );
#endif

  Validate_Lists( system, workspace, lists, data->step,
                  system->n, system->N, system->numH, comm );
}


void Estimate_Storages( reax_system *system, control_params *control,
                        reax_list **lists, int *Htop, int *hb_top,
                        int *bond_top, int *num_3body, MPI_Comm comm )
{
  int i, j, pj;
  int start_i, end_i;
  int type_i, type_j;
  int ihb, jhb;
  int local;
  real cutoff;
  real r_ij, r2;
  real C12, C34, C56;
  real BO, BO_s, BO_pi, BO_pi2;
  reax_list *far_nbrs;
  single_body_parameters *sbp_i, *sbp_j;
  two_body_parameters *twbp;
  far_neighbor_data *nbr_pj;
  reax_atom *atom_i, *atom_j;

  int mincap = system->mincap;
  double safezone = system->safezone;
  double saferzone = system->saferzone;

  far_nbrs = *lists + FAR_NBRS;
  *Htop = 0;
  memset( hb_top, 0, sizeof(int) * system->local_cap );
  memset( bond_top, 0, sizeof(int) * system->total_cap );
  *num_3body = 0;

  for( i = 0; i < system->N; ++i ) {
    atom_i = &(system->my_atoms[i]);
    type_i  = atom_i->type;
    start_i = Start_Index(i, far_nbrs);
    end_i   = End_Index(i, far_nbrs);
    sbp_i = &(system->reax_param.sbp[type_i]);

    if( i < system->n ) {
      local = 1;
      cutoff = control->nonb_cut;
      ++(*Htop);
      ihb = sbp_i->p_hbond;
    }
    else {
      local = 0;
      cutoff = control->bond_cut;
      ihb = -1;
    }

    for( pj = start_i; pj < end_i; ++pj ) {
      nbr_pj = &( far_nbrs->select.far_nbr_list[pj] );
      j = nbr_pj->nbr;
      atom_j = &(system->my_atoms[j]);

      if(nbr_pj->d <= cutoff) {
        type_j = system->my_atoms[j].type;
        r_ij = nbr_pj->d;
        sbp_j = &(system->reax_param.sbp[type_j]);
        twbp = &(system->reax_param.tbp[type_i][type_j]);

        if( local ) {
          if( j < system->n || atom_i->orig_id < atom_j->orig_id ) //tryQEq ||1
            ++(*Htop);

          /* hydrogen bond lists */
          if( control->hbond_cut > 0.1 && (ihb==1 || ihb==2) &&
              nbr_pj->d <= control->hbond_cut ) {
            jhb = sbp_j->p_hbond;
            if( ihb == 1 && jhb == 2 )
              ++hb_top[i];
            else if( j < system->n && ihb == 2 && jhb == 1 )
              ++hb_top[j];
          }
        }

        /* uncorrected bond orders */
        if( nbr_pj->d <= control->bond_cut ) {
          r2 = SQR(r_ij);

          if( sbp_i->r_s > 0.0 && sbp_j->r_s > 0.0) {
            C12 = twbp->p_bo1 * pow( r_ij / twbp->r_s, twbp->p_bo2 );
            BO_s = (1.0 + control->bo_cut) * exp( C12 );
          }
          else BO_s = C12 = 0.0;

          if( sbp_i->r_pi > 0.0 && sbp_j->r_pi > 0.0) {
            C34 = twbp->p_bo3 * pow( r_ij / twbp->r_p, twbp->p_bo4 );
            BO_pi = exp( C34 );
          }
          else BO_pi = C34 = 0.0;

          if( sbp_i->r_pi_pi > 0.0 && sbp_j->r_pi_pi > 0.0) {
            C56 = twbp->p_bo5 * pow( r_ij / twbp->r_pp, twbp->p_bo6 );
            BO_pi2= exp( C56 );
          }
          else BO_pi2 = C56 = 0.0;

          /* Initially BO values are the uncorrected ones, page 1 */
          BO = BO_s + BO_pi + BO_pi2;

          if( BO >= control->bo_cut ) {
            ++bond_top[i];
            ++bond_top[j];
          }
        }
      }
    }
  }

  *Htop = (int)(MAX( *Htop * safezone, mincap * MIN_HENTRIES ));
  for( i = 0; i < system->n; ++i )
    hb_top[i] = (int)(MAX( hb_top[i] * saferzone, MIN_HBONDS ));

  for( i = 0; i < system->N; ++i ) {
    *num_3body += SQR(bond_top[i]);
    //if( i < system->n )
    bond_top[i] = MAX( bond_top[i] * 2, MIN_BONDS );
      //else bond_top[i] = MAX_BONDS;
  }

#if defined(DEBUG_FOCUS)
  fprintf( stderr, "p%d @ estimate storages: Htop = %d, num_3body = %d\n",
           system->my_rank, *Htop, *num_3body );
  MPI_Barrier( comm );
#endif
}


void Compute_Forces( reax_system *system, control_params *control,
                     simulation_data *data, storage *workspace,
                     reax_list **lists, output_controls *out_control,
                     mpi_datatypes *mpi_data )
{
  MPI_Comm comm;
  int qeq_flag;
#if defined(LOG_PERFORMANCE)
  real t_start = 0;

  //MPI_Barrier( mpi_data->world );
  if( system->my_rank == MASTER_NODE )
    t_start = Get_Time( );
#endif

  comm = mpi_data->world;
  /********* init forces ************/
#if defined(PURE_REAX)
  if( control->qeq_freq && (data->step-data->prev_steps)%control->qeq_freq==0 )
    qeq_flag = 1;
  else qeq_flag = 0;
#elif defined(LAMMPS_REAX)
  qeq_flag = 0;
#endif

  if( qeq_flag )
    Init_Forces( system, control, data, workspace, lists, out_control, comm );
  else
    Init_Forces_noQEq( system, control, data, workspace,
                       lists, out_control, comm );

#if defined(LOG_PERFORMANCE)
  //MPI_Barrier( mpi_data->world );
  if( system->my_rank == MASTER_NODE )
    Update_Timing_Info( &t_start, &(data->timing.init_forces) );
#endif


  /********* bonded interactions ************/
  Compute_Bonded_Forces( system, control, data, workspace,
                         lists, out_control, mpi_data->world );

#if defined(LOG_PERFORMANCE)
  //MPI_Barrier( mpi_data->world );
  if( system->my_rank == MASTER_NODE )
    Update_Timing_Info( &t_start, &(data->timing.bonded) );
#endif
#if defined(DEBUG_FOCUS)
  fprintf( stderr, "p%d @ step%d: completed bonded\n",
           system->my_rank, data->step );
  MPI_Barrier( mpi_data->world );
#endif


  /**************** qeq ************************/
#if defined(PURE_REAX)
  if( qeq_flag )
    QEq( system, control, data, workspace, out_control, mpi_data );

#if defined(LOG_PERFORMANCE)
  //MPI_Barrier( mpi_data->world );
  if( system->my_rank == MASTER_NODE )
    Update_Timing_Info( &t_start, &(data->timing.qEq) );
#endif
#if defined(DEBUG_FOCUS)
  fprintf(stderr, "p%d @ step%d: qeq completed\n", system->my_rank, data->step);
  MPI_Barrier( mpi_data->world );
#endif
#endif //PURE_REAX


  /********* nonbonded interactions ************/
  Compute_NonBonded_Forces( system, control, data, workspace,
                            lists, out_control, mpi_data->world );

#if defined(LOG_PERFORMANCE)
  //MPI_Barrier( mpi_data->world );
  if( system->my_rank == MASTER_NODE )
    Update_Timing_Info( &t_start, &(data->timing.nonb) );
#endif
#if defined(DEBUG_FOCUS)
  fprintf( stderr, "p%d @ step%d: nonbonded forces completed\n",
           system->my_rank, data->step );
  MPI_Barrier( mpi_data->world );
#endif


  /*********** total force ***************/
  Compute_Total_Force( system, control, data, workspace, lists, mpi_data );

#if defined(LOG_PERFORMANCE)
  //MPI_Barrier( mpi_data->world );
  if( system->my_rank == MASTER_NODE )
    Update_Timing_Info( &t_start, &(data->timing.bonded) );
#endif
#if defined(DEBUG_FOCUS)
  fprintf( stderr, "p%d @ step%d: total forces computed\n",
           system->my_rank, data->step );
  //Print_Total_Force( system, data, workspace );
  MPI_Barrier( mpi_data->world );
#endif

#if defined(TEST_FORCES)
  Print_Force_Files( system, control, data, workspace,
                     lists, out_control, mpi_data );
#endif
}
