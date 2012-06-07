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
#include "update.h"
#if defined(PURE_REAX)
#include "io_tools.h"
#include "basic_comm.h"
#include "list.h"
#include "reset_tools.h"
#include "system_props.h"
#include "tool_box.h"
#include "traj.h"
#include "vector.h"
#elif defined(LAMMPS_REAX)
#include "reaxc_io_tools.h"
#include "reaxc_basic_comm.h"
#include "reaxc_list.h"
#include "reaxc_reset_tools.h"
#include "reaxc_system_props.h"
#include "reaxc_tool_box.h"
#include "reaxc_traj.h"
#include "reaxc_vector.h"
#endif

print_interaction Print_Interactions[NUM_INTRS];

/************************ initialize output controls ************************/
int Init_Output_Files( reax_system *system, control_params *control,
                       output_controls *out_control, mpi_datatypes *mpi_data,
                       char *msg )
{
  char temp[MAX_STR];
  int ret;

  if( out_control->write_steps > 0 ){
    ret = Init_Traj( system, control, out_control, mpi_data, msg );
    if( ret == FAILURE )
      return ret;
  }

  if( system->my_rank == MASTER_NODE ) {
    /* These files are written only by the master node */
    if( out_control->energy_update_freq > 0 ) {
#if defined(PURE_REAX)
      /* init out file */
      sprintf( temp, "%s.out", control->sim_name );
      if( (out_control->out = fopen( temp, "w" )) != NULL ) {
#if !defined(DEBUG) && !defined(DEBUG_FOCUS)
        fprintf( out_control->out, "%-6s%14s%14s%14s%11s%13s%13s\n",
                 "step", "total energy", "potential", "kinetic",
                 "T(K)", "V(A^3)", "P(Gpa)" );
#else
        fprintf( out_control->out, "%-6s%24s%24s%24s%13s%16s%13s\n",
                 "step", "total energy", "potential", "kinetic",
                 "T(K)", "V(A^3)", "P(GPa)" );
#endif
        fflush( out_control->out );
      }
      else {
        strcpy( msg, "init_out_controls: .out file could not be opened\n" );
        return FAILURE;
      }
#endif

      /* init potentials file */
      sprintf( temp, "%s.pot", control->sim_name );
      if( (out_control->pot = fopen( temp, "w" )) != NULL ) {
#if !defined(DEBUG) && !defined(DEBUG_FOCUS)
        fprintf( out_control->pot,
                 "%-6s%14s%14s%14s%14s%14s%14s%14s%14s%14s%14s%14s\n",
                 "step", "ebond", "eatom", "elp",
                 "eang", "ecoa", "ehb", "etor", "econj",
                 "evdw","ecoul", "epol" );
#else
        fprintf( out_control->pot,
                 "%-6s%24s%24s%24s%24s%24s%24s%24s%24s%24s%24s%24s\n",
                 "step", "ebond", "eatom", "elp",
                 "eang", "ecoa", "ehb", "etor", "econj",
                 "evdw","ecoul", "epol" );
#endif
        fflush( out_control->pot );
      }
      else {
        strcpy( msg, "init_out_controls: .pot file could not be opened\n" );
        return FAILURE;
      }

      /* init log file */
#if defined(LOG_PERFORMANCE)
      sprintf( temp, "%s.log", control->sim_name );
      if( (out_control->log = fopen( temp, "w" )) != NULL ) {
        fprintf( out_control->log, "%6s%8s%8s%8s%8s%8s%8s%8s%8s\n",
                 "step", "total", "comm", "nbrs", "init", "bonded", "nonb",
                 "qeq", "matvecs" );
        fflush( out_control->log );
      }
      else {
        strcpy( msg, "init_out_controls: .log file could not be opened\n" );
        return FAILURE;
      }
#endif
    }

    /* init pressure file */
    if( control->ensemble == NPT  ||
        control->ensemble == iNPT ||
        control->ensemble == sNPT ) {
      sprintf( temp, "%s.prs", control->sim_name );
      if( (out_control->prs = fopen( temp, "w" )) != NULL ) {
        fprintf(out_control->prs,"%8s%13s%13s%13s%13s%13s%13s%13s\n",
                "step", "Pint/norm[x]", "Pint/norm[y]", "Pint/norm[z]",
                "Pext/Ptot[x]", "Pext/Ptot[y]", "Pext/Ptot[z]", "Pkin/V" );
        fflush( out_control->prs );
      }
      else {
        strcpy(msg,"init_out_controls: .prs file couldn't be opened\n");
        return FAILURE;
      }
    }

    /* init electric dipole moment analysis file */
    // not yet implemented in the parallel version!!!
    // if( control->dipole_anal ) {
    //   sprintf( temp, "%s.dpl", control->sim_name );
    //   if( (out_control->dpl = fopen( temp, "w" )) != NULL ) {
    //         fprintf( out_control->dpl, "%6s%20s%30s",
    //                  "step", "molecule count", "avg dipole moment norm" );
    //         fflush( out_control->dpl );
    //   }
    //   else {
    //         strcpy(msg, "init_out_controls: .dpl file couldn't be opened\n");
    //         return FAILURE;
    //   }
    // }

    /* init diffusion coef analysis file */
    // not yet implemented in the parallel version!!!
    // if( control->diffusion_coef ) {
    //   sprintf( temp, "%s.drft", control->sim_name );
    //   if( (out_control->drft = fopen( temp, "w" )) != NULL ) {
    //         fprintf( out_control->drft, "%7s%20s%20s\n",
    //                  "step", "type count", "avg disp^2" );
    //         fflush( out_control->drft );
    //   }
    //   else {
    //         strcpy(msg,"init_out_controls: .drft file couldn't be opened\n");
    //         return FAILURE;
    //   }
    // }
  }


  /* init molecular analysis file */
  /* proc0 opens this file and shares it with everyone.
     then all processors write into it in a round-robin
     fashion controlled by their rank */
  /*if( control->molecular_analysis ) {
    if( system->my_rank == MASTER_NODE ) {
      sprintf( temp, "%s.mol", control->sim_name );
      if( (out_control->mol = fopen( temp, "w" )) == NULL ) {
        strcpy(msg,"init_out_controls: .mol file could not be opened\n");
        return FAILURE;
      }
    }

    MPI_Bcast( &(out_control->mol), 1, MPI_LONG, 0, MPI_COMM_WORLD );
    }*/


#ifdef TEST_ENERGY
  /* open bond energy file */
  sprintf( temp, "%s.ebond.%d", control->sim_name, system->my_rank );
  if( (out_control->ebond = fopen( temp, "w" )) == NULL ) {
    strcpy(msg,"Init_Out_Files: .ebond file couldn't be opened\n");
    return FAILURE;
  }

  /* open lone-pair energy file */
  sprintf( temp, "%s.elp.%d", control->sim_name, system->my_rank );
  if( (out_control->elp = fopen( temp, "w" )) == NULL ) {
    strcpy(msg,"Init_Out_Files: .elp file couldn't be opened\n");
    return FAILURE;
  }

  /* open overcoordination energy file */
  sprintf( temp, "%s.eov.%d", control->sim_name, system->my_rank );
  if( (out_control->eov = fopen( temp, "w" )) == NULL )        {
    strcpy(msg,"Init_Out_Files: .eov file couldn't be opened\n");
    return FAILURE;
  }

  /* open undercoordination energy file */
  sprintf( temp, "%s.eun.%d", control->sim_name, system->my_rank );
  if( (out_control->eun = fopen( temp, "w" )) == NULL )        {
    strcpy(msg,"Init_Out_Files: .eun file couldn't be opened\n");
    return FAILURE;
  }

  /* open angle energy file */
  sprintf( temp, "%s.eval.%d", control->sim_name, system->my_rank );
  if( (out_control->eval = fopen( temp, "w" )) == NULL ) {
    strcpy(msg,"Init_Out_Files: .eval file couldn't be opened\n");
    return FAILURE;
  }

  /* open coalition energy file */
  sprintf( temp, "%s.ecoa.%d", control->sim_name, system->my_rank );
  if( (out_control->ecoa = fopen( temp, "w" )) == NULL ) {
    strcpy(msg,"Init_Out_Files: .ecoa file couldn't be opened\n");
    return FAILURE;
  }

  /* open penalty energy file */
  sprintf( temp, "%s.epen.%d", control->sim_name, system->my_rank );
  if( (out_control->epen = fopen( temp, "w" )) == NULL ) {
    strcpy(msg,"Init_Out_Files: .epen file couldn't be opened\n");
    return FAILURE;
  }

  /* open torsion energy file */
  sprintf( temp, "%s.etor.%d", control->sim_name, system->my_rank );
  if( (out_control->etor = fopen( temp, "w" )) == NULL ) {
    strcpy(msg,"Init_Out_Files: .etor file couldn't be opened\n");
    return FAILURE;
  }

  /* open conjugation energy file */
  sprintf( temp, "%s.econ.%d", control->sim_name, system->my_rank );
  if( (out_control->econ = fopen( temp, "w" )) == NULL ) {
    strcpy(msg,"Init_Out_Files: .econ file couldn't be opened\n");
    return FAILURE;
  }

  /* open hydrogen bond energy file */
  sprintf( temp, "%s.ehb.%d", control->sim_name, system->my_rank );
  if( (out_control->ehb = fopen( temp, "w" )) == NULL )        {
    strcpy(msg,"Init_Out_Files: .ehb file couldn't be opened\n");
    return FAILURE;
  }

  /* open vdWaals energy file */
  sprintf( temp, "%s.evdw.%d", control->sim_name, system->my_rank );
  if( (out_control->evdw = fopen( temp, "w" )) == NULL ) {
    strcpy(msg,"Init_Out_Files: .evdw file couldn't be opened\n");
    return FAILURE;
  }

  /* open coulomb energy file */
  sprintf( temp, "%s.ecou.%d", control->sim_name, system->my_rank );
  if( (out_control->ecou = fopen( temp, "w" )) == NULL ) {
    strcpy(msg,"Init_Out_Files: .ecou file couldn't be opened\n");
    return FAILURE;
  }
#endif


#ifdef TEST_FORCES
  /* open bond orders file */
  sprintf( temp, "%s.fbo.%d", control->sim_name, system->my_rank );
  if( (out_control->fbo = fopen( temp, "w" )) == NULL ) {
    strcpy(msg,"Init_Out_Files: .fbo file couldn't be opened\n");
    return FAILURE;
  }

  /* open bond orders derivatives file */
  sprintf( temp, "%s.fdbo.%d", control->sim_name, system->my_rank );
  if( (out_control->fdbo = fopen( temp, "w" )) == NULL ) {
    strcpy(msg,"Init_Out_Files: .fdbo file couldn't be opened\n");
    return FAILURE;
  }

  /* produce a single force file - to be written by p0 */
  if( system->my_rank == MASTER_NODE ) {
    /* open bond forces file */
    sprintf( temp, "%s.fbond", control->sim_name );
    if( (out_control->fbond = fopen( temp, "w" )) == NULL ) {
      strcpy(msg,"Init_Out_Files: .fbond file couldn't be opened\n");
      return FAILURE;
    }

    /* open lone-pair forces file */
    sprintf( temp, "%s.flp", control->sim_name );
    if( (out_control->flp = fopen( temp, "w" )) == NULL ) {
      strcpy(msg,"Init_Out_Files: .flp file couldn't be opened\n");
      return FAILURE;
    }

    /* open overcoordination forces file */
    sprintf( temp, "%s.fov", control->sim_name );
    if( (out_control->fov = fopen( temp, "w" )) == NULL ) {
      strcpy(msg,"Init_Out_Files: .fov file couldn't be opened\n");
      return FAILURE;
    }

    /* open undercoordination forces file */
    sprintf( temp, "%s.fun", control->sim_name );
    if( (out_control->fun = fopen( temp, "w" )) == NULL ) {
      strcpy(msg,"Init_Out_Files: .fun file couldn't be opened\n");
      return FAILURE;
    }

    /* open angle forces file */
    sprintf( temp, "%s.fang", control->sim_name );
    if( (out_control->fang = fopen( temp, "w" )) == NULL ) {
      strcpy(msg,"Init_Out_Files: .fang file couldn't be opened\n");
      return FAILURE;
    }

    /* open coalition forces file */
    sprintf( temp, "%s.fcoa", control->sim_name );
    if( (out_control->fcoa = fopen( temp, "w" )) == NULL ) {
      strcpy(msg,"Init_Out_Files: .fcoa file couldn't be opened\n");
      return FAILURE;
    }

    /* open penalty forces file */
    sprintf( temp, "%s.fpen", control->sim_name );
    if( (out_control->fpen = fopen( temp, "w" )) == NULL ) {
      strcpy(msg,"Init_Out_Files: .fpen file couldn't be opened\n");
      return FAILURE;
    }

    /* open torsion forces file */
    sprintf( temp, "%s.ftor", control->sim_name );
    if( (out_control->ftor = fopen( temp, "w" )) == NULL ) {
      strcpy(msg,"Init_Out_Files: .ftor file couldn't be opened\n");
      return FAILURE;
    }

    /* open conjugation forces file */
    sprintf( temp, "%s.fcon", control->sim_name );
    if( (out_control->fcon = fopen( temp, "w" )) == NULL ) {
      strcpy(msg,"Init_Out_Files: .fcon file couldn't be opened\n");
      return FAILURE;
    }

    /* open hydrogen bond forces file */
    sprintf( temp, "%s.fhb", control->sim_name );
    if( (out_control->fhb = fopen( temp, "w" )) == NULL )        {
      strcpy(msg,"Init_Out_Files: .fhb file couldn't be opened\n");
      return FAILURE;
    }

    /* open vdw forces file */
    sprintf( temp, "%s.fvdw", control->sim_name );
    if( (out_control->fvdw = fopen( temp, "w" )) == NULL ) {
      strcpy(msg,"Init_Out_Files: .fvdw file couldn't be opened\n");
      return FAILURE;
    }

    /* open nonbonded forces file */
    sprintf( temp, "%s.fele", control->sim_name );
    if( (out_control->fele = fopen( temp, "w" )) == NULL ) {
      strcpy(msg,"Init_Out_Files: .fele file couldn't be opened\n");
      return FAILURE;
    }

    /* open total force file */
    sprintf( temp, "%s.ftot", control->sim_name );
    if( (out_control->ftot = fopen( temp, "w" )) == NULL ) {
      strcpy(msg,"Init_Out_Files: .ftot file couldn't be opened\n");
      return FAILURE;
    }

    /* open force comprison file */
    sprintf( temp, "%s.fcomp", control->sim_name );
    if( (out_control->fcomp = fopen( temp, "w" )) == NULL ) {
      strcpy(msg,"Init_Out_Files: .fcomp file couldn't be opened\n");
      return FAILURE;
    }
  }
#endif

#if defined(PURE_REAX)
#if defined(TEST_FORCES) || defined(TEST_ENERGY)
    /* open far neighbor list file */
  sprintf( temp, "%s.far_nbrs_list.%d", control->sim_name, system->my_rank );
  if( (out_control->flist = fopen( temp, "w" )) == NULL ) {
    strcpy(msg,"Init_Out_Files: .far_nbrs_list file couldn't be opened\n");
    return FAILURE;
  }

  /* open bond list file */
  sprintf( temp, "%s.bond_list.%d", control->sim_name, system->my_rank );
  if( (out_control->blist = fopen( temp, "w" )) == NULL ) {
    strcpy(msg,"Init_Out_Files: .bond_list file couldn't be opened\n");
    return FAILURE;
  }

  /* open near neighbor list file */
  sprintf( temp, "%s.near_nbrs_list.%d", control->sim_name, system->my_rank );
  if( (out_control->nlist = fopen( temp, "w" )) == NULL ) {
    strcpy(msg,"Init_Out_Files: .near_nbrs_list file couldn't be opened\n");
    return FAILURE;
  }
#endif
#endif

  return SUCCESS;
}


/************************ close output files ************************/
int Close_Output_Files( reax_system *system, control_params *control,
                        output_controls *out_control, mpi_datatypes *mpi_data )
{
  if( out_control->write_steps > 0 )
    End_Traj( system->my_rank, out_control );

  if( system->my_rank == MASTER_NODE ) {
    if( out_control->energy_update_freq > 0 ) {
#if defined(PURE_REAX)
      fclose( out_control->out );
#endif
      fclose( out_control->pot );
#if defined(LOG_PERFORMANCE)
      fclose( out_control->log );
#endif
    }

    if( control->ensemble == NPT || control->ensemble == iNPT ||
        control->ensemble == sNPT )
      fclose( out_control->prs );

    // not yet implemented in the parallel version
    //if( control->dipole_anal ) fclose( out_control->dpl );
    //if( control->diffusion_coef ) fclose( out_control->drft );
    //if( control->molecular_analysis ) fclose( out_control->mol );
  }

#ifdef TEST_ENERGY
  fclose( out_control->ebond );
  fclose( out_control->elp );
  fclose( out_control->eov );
  fclose( out_control->eun );
  fclose( out_control->eval );
  fclose( out_control->epen );
  fclose( out_control->ecoa );
  fclose( out_control->ehb );
  fclose( out_control->etor );
  fclose( out_control->econ );
  fclose( out_control->evdw );
  fclose( out_control->ecou );
#endif

#ifdef TEST_FORCES
  fclose( out_control->fbo );
  fclose( out_control->fdbo );

  if( system->my_rank == MASTER_NODE ) {
    fclose( out_control->fbond );
    fclose( out_control->flp );
    fclose( out_control->fov );
    fclose( out_control->fun );
    fclose( out_control->fang );
    fclose( out_control->fcoa );
    fclose( out_control->fpen );
    fclose( out_control->ftor );
    fclose( out_control->fcon );
    fclose( out_control->fhb );
    fclose( out_control->fvdw );
    fclose( out_control->fele );
    fclose( out_control->ftot );
    fclose( out_control->fcomp );
  }
#endif

#if defined(PURE_REAX)
#if defined(TEST_FORCES) || defined(TEST_ENERGY)
  fclose( out_control->flist );
  fclose( out_control->blist );
  fclose( out_control->nlist );
#endif
#endif

  return SUCCESS;
}



void Print_Box( simulation_box* box, char *name, FILE *out )
{
  // int i, j;

  fprintf( out, "%s:\n", name );
  fprintf( out, "\tmin[%8.3f %8.3f %8.3f]\n",
           box->min[0], box->min[1], box->min[2] );
  fprintf( out, "\tmax[%8.3f %8.3f %8.3f]\n",
           box->max[0], box->max[1], box->max[2] );
  fprintf( out, "\tdims[%8.3f%8.3f%8.3f]\n",
           box->box_norms[0], box->box_norms[1], box->box_norms[2] );

  // fprintf( out, "box: {" );
  // for( i = 0; i < 3; ++i )
  //   {
  //     fprintf( out, "{" );
  //     for( j = 0; j < 3; ++j )
  //       fprintf( out, "%8.3f ", box->box[i][j] );
  //     fprintf( out, "}" );
  //   }
  // fprintf( out, "}\n" );

  // fprintf( out, "box_trans: {" );
  // for( i = 0; i < 3; ++i )
  //   {
  //     fprintf( out, "{" );
  //     for( j = 0; j < 3; ++j )
  //         fprintf( out, "%8.3f ", box->trans[i][j] );
  //     fprintf( out, "}" );
  //   }
  // fprintf( out, "}\n" );

  // fprintf( out, "box_trinv: {" );
  // for( i = 0; i < 3; ++i )
  //   {
  //     fprintf( out, "{" );
  //     for( j = 0; j < 3; ++j )
  //         fprintf( out, "%8.3f ", box->trans_inv[i][j] );
  //     fprintf( out, "}" );
  //   }
  // fprintf( out, "}\n" );
}



void Print_Grid( grid* g, FILE *out )
{
  int x, y, z, gc_type;
  ivec gc_str;
  char gcell_type_text[10][12] =
    { "NO_NBRS", "NEAR_ONLY", "HBOND_ONLY", "FAR_ONLY",
      "NEAR_HBOND", "NEAR_FAR", "HBOND_FAR", "FULL_NBRS", "NATIVE" };

  fprintf( out, "\tnumber of grid cells: %d %d %d\n",
           g->ncells[0], g->ncells[1], g->ncells[2] );
  fprintf( out, "\tgcell lengths: %8.3f %8.3f %8.3f\n",
           g->cell_len[0], g->cell_len[1], g->cell_len[2] );
  fprintf( out, "\tinverses of gcell lengths: %8.3f %8.3f %8.3f\n",
           g->inv_len[0], g->inv_len[1], g->inv_len[2] );
  fprintf( out, "\t---------------------------------\n" );
  fprintf( out, "\tnumber of native gcells: %d %d %d\n",
           g->native_cells[0], g->native_cells[1], g->native_cells[2] );
  fprintf( out, "\tnative gcell span: %d-%d  %d-%d  %d-%d\n",
           g->native_str[0], g->native_end[0],
           g->native_str[1], g->native_end[1],
           g->native_str[2], g->native_end[2] );
  fprintf( out, "\t---------------------------------\n" );
  fprintf( out, "\tvlist gcell stretch: %d %d %d\n",
           g->vlist_span[0], g->vlist_span[1], g->vlist_span[2] );
  fprintf( out, "\tnonbonded nbrs gcell stretch: %d %d %d\n",
           g->nonb_span[0], g->nonb_span[1], g->nonb_span[2] );
  fprintf( out, "\tbonded nbrs gcell stretch: %d %d %d\n",
           g->bond_span[0], g->bond_span[1], g->bond_span[2] );
  fprintf( out, "\t---------------------------------\n" );
  fprintf( out, "\tghost gcell span: %d %d %d\n",
           g->ghost_span[0], g->ghost_span[1], g->ghost_span[2] );
  fprintf( out, "\tnonbonded ghost gcell span: %d %d %d\n",
           g->ghost_nonb_span[0],g->ghost_nonb_span[1],g->ghost_nonb_span[2]);
  fprintf(out, "\thbonded ghost gcell span: %d %d %d\n",
          g->ghost_hbond_span[0],g->ghost_hbond_span[1],g->ghost_hbond_span[2]);
  fprintf( out, "\tbonded ghost gcell span: %d %d %d\n",
           g->ghost_bond_span[0],g->ghost_bond_span[1],g->ghost_bond_span[2]);
  //fprintf(out, "\t---------------------------------\n" );
  //fprintf(out, "\tmax number of gcells at the boundary: %d\n", g->gcell_cap);
  fprintf( out, "\t---------------------------------\n" );

  fprintf( stderr, "GCELL MARKS:\n" );
  gc_type = g->cells[0][0][0].type;
  ivec_MakeZero( gc_str );

  x = y = z = 0;
  for( x = 0; x < g->ncells[0]; ++x )
    for( y = 0; y < g->ncells[1]; ++y )
      for( z = 0; z < g->ncells[2]; ++z )
        if( g->cells[x][y][z].type != gc_type ){
          fprintf( stderr,
                   "\tgcells from(%2d %2d %2d) to (%2d %2d %2d): %d - %s\n",
                   gc_str[0], gc_str[1], gc_str[2], x, y, z,
                   gc_type, gcell_type_text[gc_type] );
          gc_type = g->cells[x][y][z].type;
          gc_str[0] = x;
          gc_str[1] = y;
          gc_str[2] = z;
        }
  fprintf( stderr, "\tgcells from(%2d %2d %2d) to (%2d %2d %2d): %d - %s\n",
           gc_str[0], gc_str[1], gc_str[2], x, y, z,
           gc_type, gcell_type_text[gc_type] );
  fprintf( out, "-------------------------------------\n" );
}



void Print_GCell_Exchange_Bounds( int my_rank, neighbor_proc *my_nbrs )
{
  ivec r;
  int nbr;
  neighbor_proc *nbr_pr;
  char fname[100];
  FILE *f;
  char exch[3][10] = { "NONE", "NEAR_EXCH", "FULL_EXCH" };

  sprintf( fname, "gcell_exchange_bounds%d", my_rank );
  f = fopen( fname, "w" );

  /* loop over neighbor processes */
  for( r[0] = -1; r[0] <= 1; ++r[0])
    for( r[1] = -1; r[1] <= 1; ++r[1] )
      for( r[2] = -1; r[2] <= 1; ++r[2] )
        if( r[0]!=0 || r[1]!=0 || r[2]!=0 ) {
          nbr_pr = &(my_nbrs[nbr]);

          fprintf( f, "p%-2d GCELL BOUNDARIES with r(%2d %2d %2d):\n",
                   my_rank, r[0], r[1], r[2] );

          fprintf( f, "\tsend_type %s: send(%d %d %d) to (%d %d %d)\n",
                   exch[nbr_pr->send_type],
                   nbr_pr->str_send[0], nbr_pr->str_send[1],
                   nbr_pr->str_send[2],
                   nbr_pr->end_send[0], nbr_pr->end_send[1],
                   nbr_pr->end_send[2] );

          fprintf( f, "\trecv_type %s: recv(%d %d %d) to (%d %d %d)\n",
                   exch[nbr_pr->recv_type],
                   nbr_pr->str_recv[0], nbr_pr->str_recv[1],
                   nbr_pr->str_recv[2],
                   nbr_pr->end_recv[0], nbr_pr->end_recv[1],
                   nbr_pr->end_recv[2] );
        }

  fclose(f);
}



void Print_Native_GCells( reax_system *system )
{
  int        i, j, k, l;
  char       fname[100];
  FILE      *f;
  grid      *g;
  grid_cell *gc;
  char gcell_type_text[10][12] =
    { "NO_NBRS", "NEAR_ONLY", "HBOND_ONLY", "FAR_ONLY",
      "NEAR_HBOND", "NEAR_FAR", "HBOND_FAR", "FULL_NBRS", "NATIVE" };

  sprintf( fname, "native_gcells.%d", system->my_rank );
  f = fopen( fname, "w" );
  g = &(system->my_grid);

  for( i = g->native_str[0]; i < g->native_end[0]; i++ )
    for( j = g->native_str[1]; j < g->native_end[1]; j++ )
      for( k = g->native_str[2]; k < g->native_end[2]; k++ )
          {
            gc = &( g->cells[i][j][k] );

            fprintf( f, "p%d gcell(%2d %2d %2d) of type %d(%s)\n",
                     system->my_rank, i, j, k,
                   gc->type, gcell_type_text[gc->type] );

            fprintf( f, "\tatom list start: %d, end: %d\n\t", gc->str, gc->end );

          for( l = gc->str; l < gc->end; ++l )
              fprintf( f, "%5d", system->my_atoms[l].orig_id );
            fprintf( f, "\n" );
          }

  fclose(f);
}



void Print_All_GCells( reax_system *system )
{
  int        i, j, k, l;
  char       fname[100];
  FILE      *f;
  grid      *g;
  grid_cell *gc;
  char gcell_type_text[10][12] =
    { "NO_NBRS", "NEAR_ONLY", "HBOND_ONLY", "FAR_ONLY",
      "NEAR_HBOND", "NEAR_FAR", "HBOND_FAR", "FULL_NBRS", "NATIVE" };

  sprintf( fname, "all_gcells.%d", system->my_rank );
  f = fopen( fname, "w" );
  g = &(system->my_grid);

  for( i = 0; i < g->ncells[0]; i++ )
    for( j = 0; j < g->ncells[1]; j++ )
      for( k = 0; k < g->ncells[2]; k++ )
          {
            gc = &( g->cells[i][j][k] );

            fprintf( f, "p%d gcell(%2d %2d %2d) of type %d(%s)\n",
                     system->my_rank, i, j, k,
                   gc->type, gcell_type_text[gc->type] );

            fprintf( f, "\tatom list start: %d, end: %d\n\t", gc->str, gc->end );

          for( l = gc->str; l < gc->end; ++l )
              fprintf( f, "%5d", system->my_atoms[l].orig_id );
            fprintf( f, "\n" );
          }

  fclose(f);
}



void Print_My_Atoms( reax_system *system )
{
  int   i;
  char  fname[100];
  FILE *fh;

  sprintf( fname, "my_atoms.%d", system->my_rank );
  if( (fh = fopen( fname, "w" )) == NULL )
    {
      fprintf( stderr, "error in opening my_atoms file" );
      MPI_Abort( MPI_COMM_WORLD, FILE_NOT_FOUND );
    }

  // fprintf( stderr, "p%d had %d atoms\n",
  //   system->my_rank, system->n );

  for( i = 0; i < system->n; ++i )
    fprintf( fh, "p%-2d %-5d %2d %24.15e%24.15e%24.15e\n",
             system->my_rank,
             system->my_atoms[i].orig_id, system->my_atoms[i].type,
             system->my_atoms[i].x[0],
             system->my_atoms[i].x[1],
             system->my_atoms[i].x[2] );

  fclose( fh );
}


void Print_My_Ext_Atoms( reax_system *system )
{
  int   i;
  char  fname[100];
  FILE *fh;

  sprintf( fname, "my_ext_atoms.%d", system->my_rank );
  if( (fh = fopen( fname, "w" )) == NULL )
    {
      fprintf( stderr, "error in opening my_ext_atoms file" );
      MPI_Abort( MPI_COMM_WORLD, FILE_NOT_FOUND );
    }

  // fprintf( stderr, "p%d had %d atoms\n",
  //   system->my_rank, system->n );

  for( i = 0; i < system->N; ++i )
    fprintf( fh, "p%-2d %-5d imprt%-5d %2d %24.15e%24.15e%24.15e\n",
             system->my_rank, system->my_atoms[i].orig_id,
             system->my_atoms[i].imprt_id, system->my_atoms[i].type,
             system->my_atoms[i].x[0],
             system->my_atoms[i].x[1],
             system->my_atoms[i].x[2] );

  fclose( fh );
}


void Print_Far_Neighbors( reax_system *system, reax_list **lists,
                          control_params *control )
{
  char  fname[100];
  int   i, j, id_i, id_j, nbr, natoms;
  FILE *fout;
  reax_list *far_nbrs;

  sprintf( fname, "%s.far_nbrs.%d", control->sim_name, system->my_rank );
  fout      = fopen( fname, "w" );
  far_nbrs = (*lists) + FAR_NBRS;
  natoms = system->N;

  for( i = 0; i < natoms; ++i ) {
    id_i = system->my_atoms[i].orig_id;

    for( j = Start_Index(i,far_nbrs); j < End_Index(i,far_nbrs); ++j ) {
      nbr = far_nbrs->select.far_nbr_list[j].nbr;
      id_j = system->my_atoms[nbr].orig_id;

      fprintf( fout, "%6d%6d%24.15e%24.15e%24.15e%24.15e\n",
               id_i, id_j, far_nbrs->select.far_nbr_list[j].d,
               far_nbrs->select.far_nbr_list[j].dvec[0],
               far_nbrs->select.far_nbr_list[j].dvec[1],
               far_nbrs->select.far_nbr_list[j].dvec[2] );

      fprintf( fout, "%6d%6d%24.15e%24.15e%24.15e%24.15e\n",
               id_j, id_i, far_nbrs->select.far_nbr_list[j].d,
               -far_nbrs->select.far_nbr_list[j].dvec[0],
               -far_nbrs->select.far_nbr_list[j].dvec[1],
               -far_nbrs->select.far_nbr_list[j].dvec[2] );
    }
  }

  fclose( fout );
}


void Print_Sparse_Matrix( reax_system *system, sparse_matrix *A )
{
  int i, j;

  for( i = 0; i < A->n; ++i )
    for( j = A->start[i]; j < A->end[i]; ++j )
      fprintf( stderr, "%d %d %.15e\n",
               system->my_atoms[i].orig_id,
               system->my_atoms[A->entries[j].j].orig_id,
               A->entries[j].val );
}


void Print_Sparse_Matrix2( reax_system *system, sparse_matrix *A, char *fname )
{
  int i, j;
  FILE *f = fopen( fname, "w" );

  for( i = 0; i < A->n; ++i )
    for( j = A->start[i]; j < A->end[i]; ++j )
      fprintf( f, "%d %d %.15e\n",
               system->my_atoms[i].orig_id,
               system->my_atoms[A->entries[j].j].orig_id,
               A->entries[j].val );

  fclose(f);
}


void Print_Symmetric_Sparse(reax_system *system, sparse_matrix *A, char *fname)
{
  int i, j;
  reax_atom *ai, *aj;
  FILE *f = fopen( fname, "w" );

  for( i = 0; i < A->n; ++i ) {
    ai = &(system->my_atoms[i]);
    for( j = A->start[i]; j < A->end[i]; ++j ) {
      aj = &(system->my_atoms[A->entries[j].j]);
      fprintf( f, "%d %d %.15e\n",
               ai->renumber, aj->renumber, A->entries[j].val );
      if( A->entries[j].j < system->n && ai->renumber != aj->renumber )
        fprintf( f, "%d %d %.15e\n",
                 aj->renumber, ai->renumber, A->entries[j].val );
    }
  }

  fclose(f);
}


void Print_Linear_System( reax_system *system, control_params *control,
                          storage *workspace, int step )
{
  int   i, j;
  char  fname[100];
  reax_atom *ai, *aj;
  sparse_matrix *H;
  FILE *out;

  // print rhs and init guesses for QEq
  sprintf( fname, "%s.p%dstate%d", control->sim_name, system->my_rank, step );
  out = fopen( fname, "w" );
  for( i = 0; i < system->n; i++ ) {
    ai = &(system->my_atoms[i]);
    fprintf( out, "%6d%2d%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e\n",
             ai->renumber, ai->type, ai->x[0], ai->x[1], ai->x[2],
             workspace->s[i], workspace->b_s[i],
             workspace->t[i], workspace->b_t[i] );
  }
  fclose( out );

  // print QEq coef matrix
  sprintf( fname, "%s.p%dH%d", control->sim_name, system->my_rank, step );
  Print_Symmetric_Sparse( system, workspace->H, fname );

  // print the incomplete H matrix
  /*sprintf( fname, "%s.p%dHinc%d", control->sim_name, system->my_rank, step );
  out = fopen( fname, "w" );
  H = workspace->H;
  for( i = 0; i < H->n; ++i ) {
    ai = &(system->my_atoms[i]);
    for( j = H->start[i]; j < H->end[i]; ++j )
      if( H->entries[j].j < system->n ) {
        aj = &(system->my_atoms[H->entries[j].j]);
        fprintf( out, "%d %d %.15e\n",
                 ai->orig_id, aj->orig_id, H->entries[j].val );
        if( ai->orig_id != aj->orig_id )
          fprintf( out, "%d %d %.15e\n",
                   aj->orig_id, ai->orig_id, H->entries[j].val );
      }
  }
  fclose( out );*/

  // print the L from incomplete cholesky decomposition
  /*sprintf( fname, "%s.p%dL%d", control->sim_name, system->my_rank, step );
    Print_Sparse_Matrix2( system, workspace->L, fname );*/
}


void Print_LinSys_Soln( reax_system *system, real *x, real *b_prm, real *b )
{
  int    i;
  char   fname[100];
  FILE  *fout;

  sprintf( fname, "qeq.%d.out", system->my_rank );
  fout = fopen( fname, "w" );

  for( i = 0; i < system->n; ++i )
    fprintf( fout, "%6d%10.4f%10.4f%10.4f\n",
             system->my_atoms[i].orig_id, x[i], b_prm[i], b[i] );

  fclose( fout );
}


void Print_Charges( reax_system *system )
{
  int    i;
  char   fname[100];
  FILE  *fout;

  sprintf( fname, "q.%d.out", system->my_rank );
  fout = fopen( fname, "w" );

  for( i = 0; i < system->n; ++i )
    fprintf( fout, "%6d %10.7f %10.7f %10.7f\n",
             system->my_atoms[i].orig_id,
             system->my_atoms[i].s[0],
             system->my_atoms[i].t[0],
             system->my_atoms[i].q );

  fclose( fout );
}


void Print_Bonds( reax_system *system, reax_list *bonds, char *fname )
{
  int i, j, pj;
  bond_data *pbond;
  bond_order_data *bo_ij;
  FILE *f = fopen( fname, "w" );

  for( i = 0; i < system->N; ++i )
    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj ) {
      pbond = &(bonds->select.bond_list[pj]);
      bo_ij = &(pbond->bo_data);
      j = pbond->nbr;
      //fprintf( f, "%6d%6d%23.15e%23.15e%23.15e%23.15e%23.15e\n",
      //       system->my_atoms[i].orig_id, system->my_atoms[j].orig_id,
      //       pbond->d, bo_ij->BO, bo_ij->BO_s, bo_ij->BO_pi, bo_ij->BO_pi2 );
      fprintf( f, "%8d%8d %24.15f %24.15f\n",
               i, j,//system->my_atoms[i].orig_id, system->my_atoms[j].orig_id,
               pbond->d, bo_ij->BO );
    }

  fclose(f);
}


int fn_qsort_intcmp( const void *a, const void *b )
{
  return( *(int *)a - *(int *)b );
}

void Print_Bond_List2( reax_system *system, reax_list *bonds, char *fname )
{
  int i,j, id_i, id_j, nbr, pj;
  FILE *f = fopen( fname, "w" );
  int temp[500];
  int num=0;

  for( i = 0; i < system->n; ++i ) {
    num=0;
    id_i = system->my_atoms[i].orig_id;
    fprintf( f, "%6d:", id_i);
    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj ) {
      nbr = bonds->select.bond_list[pj].nbr;
      id_j = system->my_atoms[nbr].orig_id;
      if( id_i < id_j )
        temp[num++] = id_j;
    }

    qsort(&temp, num, sizeof(int), fn_qsort_intcmp);
    for(j=0; j < num; j++)
      fprintf(f, "%6d", temp[j] );
    fprintf(f, "\n");
  }
}


void Print_Total_Force( reax_system *system, simulation_data *data,
                        storage *workspace )
{
  int    i;

  fprintf( stderr, "step: %d\n", data->step );
  fprintf( stderr, "%6s\t%-38s\n", "atom", "atom.f[0,1,2]");

  for( i = 0; i < system->N; ++i )
    fprintf( stderr, "%6d %f %f %f\n",
             //"%6d%24.15e%24.15e%24.15e\n",
             system->my_atoms[i].orig_id,
             workspace->f[i][0], workspace->f[i][1], workspace->f[i][2] );
}

void fixbond( reax_system *system, control_params *control,
                     simulation_data *data, reax_list **lists,
                     output_controls *out_control, mpi_datatypes *mpi_data )
{
  // count the number of bonds around each atom, for fix reax/c/bond
  int i, j, pj, my_bonds_0, i_id, j_id;
  int my_bonds_max = 0;
  double BO_tmp;

  control->bg_cut = 0.3;  // this values will not change with control file
  reax_list *bonds = (*lists) + BONDS;

  for( i=0; i < system->n; ++i ) {
    my_bonds_0 = 0;
    i_id = system->my_atoms[i].orig_id;  // orig_id is atom->tag
    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj ) {
      j = bonds->select.bond_list[pj].nbr;
      j_id = system->my_atoms[j].orig_id;
      BO_tmp = bonds->select.bond_list[pj].bo_data.BO;
      if( i_id != j_id && BO_tmp >= control->bg_cut ) {
        ++my_bonds_0;
        system->my_atoms[i].nbr_id[my_bonds_0] = j_id;
        system->my_atoms[i].nbr_bo[my_bonds_0] = BO_tmp;
      }
    }
    my_bonds_max = MAX(my_bonds_0, my_bonds_max);
    system->my_atoms[i].numbonds = my_bonds_0;
    system->my_bonds = my_bonds_max;
  }
}


void fixspecies( reax_system *system, control_params *control,
                     simulation_data *data, reax_list **lists,
                     output_controls *out_control, mpi_datatypes *mpi_data )
{
  // count the number of bonds around each atom, for fix reax/c/bond
  int i, j, pj, my_bonds_0, i_id, j_id;
  int my_bonds_max = 0;
  double BO_tmp;

  control->bg_cut = 0.3;  // this values will not change with control file
  reax_list *bonds = (*lists) + BONDS;

  for( i=0; i < system->n; ++i ) {
    my_bonds_0 = 0;
    i_id = system->my_atoms[i].orig_id;
    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj ) {
      j = bonds->select.bond_list[pj].nbr;
      j_id = system->my_atoms[j].orig_id;
      BO_tmp = bonds->select.bond_list[pj].bo_data.BO;
      if( i_id != j_id && BO_tmp >= control->bg_cut ) {
        ++my_bonds_0;
        system->my_atoms[i].nbr_id[my_bonds_0] = j_id;
        system->my_atoms[i].nbr_bo[my_bonds_0] = BO_tmp;
      }
    }
    my_bonds_max = MAX(my_bonds_0, my_bonds_max);
    system->my_atoms[i].numbonds = my_bonds_0;
    system->my_bonds = my_bonds_max;
  }
}


void Output_Results( reax_system *system, control_params *control,
                     simulation_data *data, reax_list **lists,
                     output_controls *out_control, mpi_datatypes *mpi_data )
{
#if defined(LOG_PERFORMANCE)
  real t_elapsed, denom;
#endif

  if((out_control->energy_update_freq > 0 &&
      data->step%out_control->energy_update_freq == 0) ||
     (out_control->write_steps > 0 &&
      data->step%out_control->write_steps == 0)){
    /* update system-wide energies */
    Compute_System_Energy( system, data, mpi_data->world );

    /* output energies */
    if( system->my_rank == MASTER_NODE &&
        out_control->energy_update_freq > 0 &&
        data->step % out_control->energy_update_freq == 0 ) {
#if !defined(DEBUG) && !defined(DEBUG_FOCUS)
#if defined(PURE_REAX)
      fprintf( out_control->out,
               "%-6d%14.2f%14.2f%14.2f%11.2f%13.2f%13.5f\n",
               data->step, data->sys_en.e_tot, data->sys_en.e_pot,
               E_CONV * data->sys_en.e_kin, data->therm.T,
               system->big_box.V, data->iso_bar.P );
      fflush( out_control->out );
#endif

      fprintf( out_control->pot,
               "%-6d%14.2f%14.2f%14.2f%14.2f%14.2f%14.2f%14.2f%14.2f%14.2f%14.2f%14.2f\n",
               data->step,
               data->sys_en.e_bond,
               data->sys_en.e_ov + data->sys_en.e_un,  data->sys_en.e_lp,
               data->sys_en.e_ang + data->sys_en.e_pen, data->sys_en.e_coa,
               data->sys_en.e_hb,
               data->sys_en.e_tor, data->sys_en.e_con,
               data->sys_en.e_vdW, data->sys_en.e_ele, data->sys_en.e_pol);
      fflush( out_control->pot );
#else
#if defined(PURE_REAX)
      fprintf( out_control->out,
               "%-6d%24.15e%24.15e%24.15e%13.5f%16.5f%13.5f\n",
               data->step, data->sys_en.e_tot, data->sys_en.e_pot,
               E_CONV * data->sys_en.e_kin, data->therm.T,
               system->big_box.V, data->iso_bar.P );
      fflush( out_control->out );
#endif

      fprintf( out_control->pot,
               "%-6d%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e\n",
               data->step,
               data->sys_en.e_bond,
               data->sys_en.e_ov + data->sys_en.e_un,  data->sys_en.e_lp,
               data->sys_en.e_ang + data->sys_en.e_pen, data->sys_en.e_coa,
               data->sys_en.e_hb,
               data->sys_en.e_tor, data->sys_en.e_con,
               data->sys_en.e_vdW, data->sys_en.e_ele, data->sys_en.e_pol);
      fflush( out_control->pot );
#endif //DEBUG

#if defined(LOG_PERFORMANCE)
      t_elapsed = Get_Timing_Info( data->timing.total );
      if( data->step - data->prev_steps > 0 )
        denom = 1.0 / out_control->energy_update_freq;
      else denom = 1;

      fprintf( out_control->log, "%6d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%6d\n",
               data->step,
               t_elapsed * denom,
               data->timing.comm * denom,
               data->timing.nbrs * denom,
               data->timing.init_forces * denom,
               data->timing.bonded * denom,
               data->timing.nonb * denom,
               data->timing.qEq * denom,
               (int)((data->timing.s_matvecs+data->timing.t_matvecs)*denom) );

      Reset_Timing( &(data->timing) );
      fflush( out_control->log );
#endif //LOG_PERFORMANCE

      if( control->virial ){
        fprintf( out_control->prs,
                 "%8d%13.6f%13.6f%13.6f%13.6f%13.6f%13.6f%13.6f\n",
                 data->step,
                 data->int_press[0], data->int_press[1], data->int_press[2],
                 data->ext_press[0], data->ext_press[1], data->ext_press[2],
                 data->kin_press );

        fprintf( out_control->prs,
                 "%8s%13.6f%13.6f%13.6f%13.6f%13.6f%13.6f%13.6f\n",
                 "",system->big_box.box_norms[0], system->big_box.box_norms[1],
                 system->big_box.box_norms[2],
                 data->tot_press[0], data->tot_press[1], data->tot_press[2],
                 system->big_box.V );

        fflush( out_control->prs);
      }
    }

    /* write current frame */
    if( out_control->write_steps > 0 &&
        (data->step-data->prev_steps) % out_control->write_steps == 0 ) {
      Append_Frame( system, control, data, lists, out_control, mpi_data );
    }
  }

#if defined(DEBUG)
  fprintf( stderr, "output_results... done\n" );
#endif
}


#ifdef TEST_ENERGY
void Debug_Marker_Bonded( output_controls *out_control, int step )
{
  fprintf( out_control->ebond, "step: %d\n%6s%6s%12s%12s%12s\n",
           step, "atom1", "atom2", "bo", "ebond", "total" );
  fprintf( out_control->elp, "step: %d\n%6s%12s%12s%12s\n",
           step, "atom", "nlp", "elp", "total" );
  fprintf( out_control->eov, "step: %d\n%6s%12s%12s\n",
           step, "atom", "eov", "total" );
  fprintf( out_control->eun, "step: %d\n%6s%12s%12s\n",
           step, "atom", "eun", "total" );
  fprintf( out_control->eval, "step: %d\n%6s%6s%6s%12s%12s%12s%12s%12s%12s\n",
           step, "atom1", "atom2", "atom3", "angle", "theta0",
           "bo(12)", "bo(23)", "eval", "total" );
  fprintf( out_control->epen, "step: %d\n%6s%6s%6s%12s%12s%12s%12s%12s\n",
           step, "atom1", "atom2", "atom3", "angle", "bo(12)", "bo(23)",
           "epen", "total" );
  fprintf( out_control->ecoa, "step: %d\n%6s%6s%6s%12s%12s%12s%12s%12s\n",
           step, "atom1", "atom2", "atom3", "angle", "bo(12)", "bo(23)",
           "ecoa", "total" );
  fprintf( out_control->ehb,  "step: %d\n%6s%6s%6s%12s%12s%12s%12s%12s\n",
           step, "atom1", "atom2", "atom3", "r(23)", "angle", "bo(12)",
           "ehb", "total" );
  fprintf( out_control->etor, "step: %d\n%6s%6s%6s%6s%12s%12s%12s%12s\n",
           step, "atom1", "atom2", "atom3", "atom4", "phi", "bo(23)",
           "etor", "total" );
  fprintf( out_control->econ,"step:%d\n%6s%6s%6s%6s%12s%12s%12s%12s%12s%12s\n",
           step, "atom1", "atom2", "atom3", "atom4",
           "phi", "bo(12)", "bo(23)", "bo(34)", "econ", "total" );
}

void Debug_Marker_Nonbonded( output_controls *out_control, int step )
{
  fprintf( out_control->evdw, "step: %d\n%6s%6s%12s%12s%12s\n",
           step, "atom1", "atom2", "r12", "evdw", "total" );
  fprintf( out_control->ecou, "step: %d\n%6s%6s%12s%12s%12s%12s%12s\n",
           step, "atom1", "atom2", "r12", "q1", "q2", "ecou", "total" );
}

#endif


#ifdef TEST_FORCES
void Dummy_Printer( reax_system *system, control_params *control,
                    simulation_data *data, storage *workspace,
                    reax_list **lists, output_controls *out_control )
{
}



void Print_Bond_Orders( reax_system *system, control_params *control,
                        simulation_data *data, storage *workspace,
                        reax_list **lists, output_controls *out_control )
{
  int i, pj, pk;
  bond_order_data *bo_ij;
  dbond_data *dbo_k;
  reax_list *bonds = (*lists) + BONDS;
  reax_list *dBOs  = (*lists) + DBOS;

  /* bond orders */
  fprintf( out_control->fbo, "step: %d\n", data->step );
  fprintf( out_control->fbo, "%6s%6s%12s%12s%12s%12s%12s\n",
           "atom1", "atom2", "r_ij", "total_bo", "bo_s", "bo_p", "bo_pp" );

  for( i = 0; i < system->N; ++i )
    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj ) {
        bo_ij = &(bonds->select.bond_list[pj].bo_data);
        fprintf( out_control->fbo,
                 "%6d%6d%24.15e%24.15e%24.15e%24.15e%24.15e\n",
                 system->my_atoms[i].orig_id,
                 system->my_atoms[bonds->select.bond_list[pj].nbr].orig_id,
                 bonds->select.bond_list[pj].d,
                 bo_ij->BO, bo_ij->BO_s, bo_ij->BO_pi, bo_ij->BO_pi2 );
    }


  /* derivatives of bond orders */
  fprintf( out_control->fdbo, "step: %d\n", data->step );
  fprintf( out_control->fdbo, "%6s%6s%6s%24s%24s%24s\n",
           "atom1", "atom2", "atom2", "dBO", "dBOpi", "dBOpi2" );
  for( i = 0; i < system->N; ++i )
    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj ) {
      /* fprintf( out_control->fdbo, "%6d %6d\tstart: %6d\tend: %6d\n",
         system->my_atoms[i].orig_id,
         system->my_atoms[bonds->select.bond_list[pj].nbr].orig_id,
         Start_Index( pj, dBOs ), End_Index( pj, dBOs ) ); */
      for( pk = Start_Index(pj, dBOs); pk < End_Index(pj, dBOs); ++pk ) {
        dbo_k = &(dBOs->select.dbo_list[pk]);
        fprintf( out_control->fdbo, "%6d%6d%6d%24.15e%24.15e%24.15e\n",
                 system->my_atoms[i].orig_id,
                 system->my_atoms[bonds->select.bond_list[pj].nbr].orig_id,
                 system->my_atoms[dbo_k->wrt].orig_id,
                 dbo_k->dBO[0], dbo_k->dBO[1], dbo_k->dBO[2] );

        fprintf( out_control->fdbo, "%6d%6d%6d%24.15e%24.15e%24.15e\n",
                 system->my_atoms[i].orig_id,
                 system->my_atoms[bonds->select.bond_list[pj].nbr].orig_id,
                 system->my_atoms[dbo_k->wrt].orig_id,
                 dbo_k->dBOpi[0], dbo_k->dBOpi[1], dbo_k->dBOpi[2] );

        fprintf( out_control->fdbo, "%6d%6d%6d%24.15e%24.15e%24.15e\n",
                 system->my_atoms[i].orig_id,
                 system->my_atoms[bonds->select.bond_list[pj].nbr].orig_id,
                 system->my_atoms[dbo_k->wrt].orig_id,
                 dbo_k->dBOpi2[0], dbo_k->dBOpi2[1], dbo_k->dBOpi2[2] );
      }
    }
}


void Print_Forces( FILE *f, storage *workspace, int N, int step )
{
  int i;

  fprintf( f, "step: %d\n", step );
  for( i = 0; i < N; ++i )
    //fprintf( f, "%6d %23.15e %23.15e %23.15e\n",
    //fprintf( f, "%6d%12.6f%12.6f%12.6f\n",
    fprintf( f, "%6d %19.9e %19.9e %19.9e\n",
             workspace->id_all[i], workspace->f_all[i][0],
             workspace->f_all[i][1], workspace->f_all[i][2] );
}


void Print_Force_Files( reax_system *system, control_params *control,
                        simulation_data *data, storage *workspace,
                        reax_list **lists, output_controls *out_control,
                        mpi_datatypes *mpi_data )
{
  int i, d;

  Coll_ids_at_Master( system, workspace, mpi_data );

  Print_Bond_Orders( system, control, data, workspace, lists, out_control );

  Coll_rvecs_at_Master( system, workspace, mpi_data, workspace->f_be );
  if( system->my_rank == MASTER_NODE )
    Print_Forces( out_control->fbond, workspace, system->bigN, data->step );

  Coll_rvecs_at_Master( system, workspace, mpi_data, workspace->f_lp );
  if( system->my_rank == MASTER_NODE )
    Print_Forces( out_control->flp, workspace, system->bigN, data->step );

  Coll_rvecs_at_Master( system, workspace, mpi_data, workspace->f_ov );
  if( system->my_rank == MASTER_NODE )
    Print_Forces( out_control->fov, workspace, system->bigN, data->step );

  Coll_rvecs_at_Master( system, workspace, mpi_data, workspace->f_un );
  if( system->my_rank == MASTER_NODE )
    Print_Forces( out_control->fun, workspace, system->bigN, data->step );

  Coll_rvecs_at_Master( system, workspace, mpi_data, workspace->f_ang );
  if( system->my_rank == MASTER_NODE )
    Print_Forces( out_control->fang, workspace, system->bigN, data->step );

  Coll_rvecs_at_Master( system, workspace, mpi_data, workspace->f_coa );
  if( system->my_rank == MASTER_NODE )
    Print_Forces( out_control->fcoa, workspace, system->bigN, data->step );

  Coll_rvecs_at_Master( system, workspace, mpi_data, workspace->f_pen );
  if( system->my_rank == MASTER_NODE )
    Print_Forces( out_control->fpen, workspace, system->bigN, data->step );

  Coll_rvecs_at_Master( system, workspace, mpi_data, workspace->f_tor );
  if( system->my_rank == MASTER_NODE )
    Print_Forces( out_control->ftor, workspace, system->bigN, data->step );

  Coll_rvecs_at_Master( system, workspace, mpi_data, workspace->f_con );
  if( system->my_rank == MASTER_NODE )
    Print_Forces( out_control->fcon, workspace, system->bigN, data->step );

  Coll_rvecs_at_Master( system, workspace, mpi_data, workspace->f_hb );
  if( system->my_rank == MASTER_NODE )
    Print_Forces( out_control->fhb, workspace, system->bigN, data->step );

  Coll_rvecs_at_Master( system, workspace, mpi_data, workspace->f_vdw );
  if( system->my_rank == MASTER_NODE )
    Print_Forces( out_control->fvdw, workspace, system->bigN, data->step );

  Coll_rvecs_at_Master( system, workspace, mpi_data, workspace->f_ele );
  if( system->my_rank == MASTER_NODE )
    Print_Forces( out_control->fele, workspace, system->bigN, data->step );

  Coll_rvecs_at_Master( system, workspace, mpi_data, workspace->f );
  if( system->my_rank == MASTER_NODE )
    Print_Forces( out_control->ftot, workspace, system->bigN, data->step );

  for( i = 0; i < system->n; ++i ) {
    for( d = 0; d < 3; ++d )
      workspace->f_tot[i][d] = workspace->f_be[i][d] +
        workspace->f_lp[i][d]+workspace->f_ov[i][d]+workspace->f_un[i][d] +
        workspace->f_ang[i][d]+workspace->f_pen[i][d]+workspace->f_coa[i][d] +
        workspace->f_tor[i][d]+workspace->f_con[i][d] +
        workspace->f_vdw[i][d]+workspace->f_ele[i][d] +
        workspace->f_hb[i][d];
  }

  Coll_rvecs_at_Master( system, workspace, mpi_data, workspace->f_tot );
  if( system->my_rank == MASTER_NODE )
    Print_Forces( out_control->fcomp, workspace, system->bigN, data->step );
}
#endif


#if defined(TEST_FORCES) || defined(TEST_ENERGY)

void Print_Far_Neighbors_List( reax_system *system, reax_list **lists,
                               control_params *control, simulation_data *data,
                               output_controls *out_control )
{
  int   i, j, id_i, id_j, nbr, natoms;
  int num=0;
  int temp[500];
  reax_list *far_nbrs;

  far_nbrs = (*lists) + FAR_NBRS;
  fprintf( out_control->flist, "step: %d\n", data->step );
  fprintf( out_control->flist, "%6s\t%-38s\n", "atom", "Far_nbrs_list");


  natoms = system->n;
  for( i = 0; i < natoms; ++i ) {
    id_i = system->my_atoms[i].orig_id;
    fprintf( out_control->flist, "%6d:",id_i);
    num=0;

    for( j = Start_Index(i,far_nbrs); j < End_Index(i,far_nbrs); ++j ) {
      nbr = far_nbrs->select.far_nbr_list[j].nbr;
      id_j = system->my_atoms[nbr].orig_id;
      temp[num++] = id_j;
    }

    qsort(&temp, num, sizeof(int), fn_qsort_intcmp);
    for(j=0; j < num; j++)
      fprintf(out_control->flist, "%6d",temp[j]);
    fprintf( out_control->flist, "\n");
  }
}

void Print_Bond_List( reax_system *system, control_params *control,
                      simulation_data *data, reax_list **lists,
                      output_controls *out_control)
{
  int i,j, id_i, id_j, nbr, pj;
  reax_list *bonds = (*lists) + BONDS;

  int temp[500];
  int num=0;

  fprintf( out_control->blist, "step: %d\n", data->step );
  fprintf( out_control->blist, "%6s\t%-38s\n", "atom", "Bond_list");

  /* bond list */
  for( i = 0; i < system->n; ++i ) {
    num=0;
    id_i = system->my_atoms[i].orig_id;
    fprintf( out_control->blist, "%6d:", id_i);
    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj ) {
      nbr = bonds->select.bond_list[pj].nbr;
      id_j = system->my_atoms[nbr].orig_id;
      if( id_i < id_j )
        temp[num++] = id_j;
    }

    qsort(&temp, num, sizeof(int), fn_qsort_intcmp);
    for(j=0; j < num; j++)
      fprintf(out_control->blist, "%6d",temp[j]);
    fprintf(out_control->blist, "\n");
  }
}


#endif


#ifdef OLD_VERSION
void Print_Init_Atoms( reax_system *system, storage *workspace )
{
  int i;

  fprintf( stderr, "p%d had %d atoms\n",
           system->my_rank, workspace->init_cnt );

  for( i = 0; i < workspace->init_cnt; ++i )
    fprintf( stderr, "p%d, atom%d: %d  %s  %8.3f %8.3f %8.3f\n",
             system->my_rank, i,
             workspace->init_atoms[i].type, workspace->init_atoms[i].name,
             workspace->init_atoms[i].x[0],
             workspace->init_atoms[i].x[1],
             workspace->init_atoms[i].x[2] );
}
#endif //OLD_VERSION


/*void Print_Bond_Forces( reax_system *system, control_params *control,
                        simulation_data *data, storage *workspace,
                        reax_list **lists, output_controls *out_control )
{
  int i;

  fprintf( out_control->fbond, "step: %d\n", data->step );
  fprintf( out_control->fbond, "%6s%24s%24s%24s\n",
           "atom", "f_be[0]", "f_be[1]", "f_be[2]" );

  for( i = 0; i < system->bigN; ++i )
    fprintf(out_control->fbond, "%6d%24.15e%24.15e%24.15e\n",
            system->my_atoms[i].orig_id,
            workspace->f_all[i][0], workspace->f_all[i][1],
            workspace->f_all[i][2]);
}

void Print_LonePair_Forces( reax_system *system, control_params *control,
                            simulation_data *data, storage *workspace,
                            reax_list **lists, output_controls *out_control )
{
  int i;

  fprintf( out_control->flp, "step: %d\n", data->step );
  fprintf( out_control->flp, "%6s%24s\n", "atom", "f_lonepair" );

  for( i = 0; i < system->bigN; ++i )
    fprintf(out_control->flp, "%6d%24.15e%24.15e%24.15e\n",
            system->my_atoms[i].orig_id,
            workspace->f_all[i][0], workspace->f_all[i][1],
            workspace->f_all[i][2]);
}


void Print_OverCoor_Forces( reax_system *system, control_params *control,
                            simulation_data *data, storage *workspace,
                            reax_list **lists, output_controls *out_control )
{
  int i;

  fprintf( out_control->fov, "step: %d\n", data->step );
  fprintf( out_control->fov, "%6s%-38s%-38s%-38s\n",
           "atom","f_over[0]", "f_over[1]", "f_over[2]" );

  for( i = 0; i < system->bigN; ++i )
    fprintf( out_control->fov,
             "%6d %24.15e%24.15e%24.15e 0 0 0\n",
             system->my_atoms[i].orig_id,
             workspace->f_all[i][0], workspace->f_all[i][1],
             workspace->f_all[i][2] );
}


void Print_UnderCoor_Forces( reax_system *system, control_params *control,
                             simulation_data *data, storage *workspace,
                             reax_list **lists, output_controls *out_control )
{
  int i;

  fprintf( out_control->fun, "step: %d\n", data->step );
  fprintf( out_control->fun, "%6s%-38s%-38s%-38s\n",
           "atom","f_under[0]", "f_under[1]", "f_under[2]" );

  for( i = 0; i < system->bigN; ++i )
    fprintf( out_control->fun,
             "%6d %24.15e%24.15e%24.15e 0 0 0\n",
             system->my_atoms[i].orig_id,
             workspace->f_all[i][0], workspace->f_all[i][1],
             workspace->f_all[i][2] );
}


void Print_ValAngle_Forces( reax_system *system, control_params *control,
                            simulation_data *data, storage *workspace,
                            reax_list **lists, output_controls *out_control )
{
  int j;

  fprintf( out_control->f3body, "step: %d\n", data->step );
  fprintf( out_control->f3body, "%6s%-37s%-37s%-37s%-38s\n",
           "atom", "3-body total", "f_ang", "f_pen", "f_coa" );

  for( j = 0; j < system->N; ++j ){
    if( rvec_isZero(workspace->f_pen[j]) && rvec_isZero(workspace->f_coa[j]) )
      fprintf( out_control->f3body,
               "%6d %24.15e%24.15e%24.15e  0 0 0  0 0 0\n",
               system->my_atoms[j].orig_id,
               workspace->f_ang[j][0], workspace->f_ang[j][1],
               workspace->f_ang[j][2] );
    else if( rvec_isZero(workspace->f_coa[j]) )
      fprintf( out_control->f3body,
               "%6d %24.15e%24.15e%24.15e %24.15e%24.15e%24.15e "        \
               "%24.15e%24.15e%24.15e\n",
               system->my_atoms[j].orig_id,
               workspace->f_ang[j][0] + workspace->f_pen[j][0],
               workspace->f_ang[j][1] + workspace->f_pen[j][1],
               workspace->f_ang[j][2] + workspace->f_pen[j][2],
               workspace->f_ang[j][0], workspace->f_ang[j][1],
               workspace->f_ang[j][2],
               workspace->f_pen[j][0], workspace->f_pen[j][1],
               workspace->f_pen[j][2] );
    else{
      fprintf( out_control->f3body, "%6d %24.15e%24.15e%24.15e ",
             system->my_atoms[j].orig_id,
               workspace->f_ang[j][0] + workspace->f_pen[j][0] +
               workspace->f_coa[j][0],
               workspace->f_ang[j][1] + workspace->f_pen[j][1] +
               workspace->f_coa[j][1],
               workspace->f_ang[j][2] + workspace->f_pen[j][2] +
               workspace->f_coa[j][2] );

      fprintf( out_control->f3body,
               "%24.15e%24.15e%24.15e %24.15e%24.15e%24.15e "\
               "%24.15e%24.15e%24.15e\n",
               workspace->f_ang[j][0], workspace->f_ang[j][1],
               workspace->f_ang[j][2],
               workspace->f_pen[j][0], workspace->f_pen[j][1],
               workspace->f_pen[j][2],
               workspace->f_coa[j][0], workspace->f_coa[j][1],
               workspace->f_coa[j][2] );
        }
  }
}


void Print_Hydrogen_Bond_Forces( reax_system *system, control_params *control,
                                 simulation_data *data, storage *workspace,
                                 reax_list **lists, output_controls *out_control)
{
  int j;

  fprintf( out_control->fhb, "step: %d\n", data->step );
  fprintf( out_control->fhb, "%6s\t%-38s\n", "atom", "f_hb[0,1,2]" );

  for( j = 0; j < system->N; ++j )
    fprintf(out_control->fhb, "%6d%24.15e%24.15e%24.15e\n",
             system->my_atoms[j].orig_id,
             workspace->f_hb[j][0],
             workspace->f_hb[j][1],
             workspace->f_hb[j][2] );
}


void Print_Four_Body_Forces( reax_system *system, control_params *control,
                             simulation_data *data, storage *workspace,
                             reax_list **lists, output_controls *out_control )
{
  int j;

  fprintf( out_control->f4body, "step: %d\n", data->step );
  fprintf( out_control->f4body, "%6s\t%-38s%-38s%-38s\n",
           "atom", "4-body total", "f_tor", "f_con" );

  for( j = 0; j < system->N; ++j ){
    if( !rvec_isZero( workspace->f_con[j] ) )
      fprintf( out_control->f4body,
               "%6d %24.15e%24.15e%24.15e %24.15e%24.15e%24.15e "\
               "%24.15e%24.15e%24.15e\n",
             system->my_atoms[j].orig_id,
               workspace->f_tor[j][0] + workspace->f_con[j][0],
               workspace->f_tor[j][1] + workspace->f_con[j][1],
               workspace->f_tor[j][2] + workspace->f_con[j][2],
               workspace->f_tor[j][0], workspace->f_tor[j][1],
               workspace->f_tor[j][2],
               workspace->f_con[j][0], workspace->f_con[j][1],
               workspace->f_con[j][2] );
    else
      fprintf( out_control->f4body,
               "%6d %24.15e%24.15e%24.15e  0 0 0\n",
               system->my_atoms[j].orig_id, workspace->f_tor[j][0],
               workspace->f_tor[j][1], workspace->f_tor[j][2] );
  }

}


void Print_vdW_Coulomb_Forces( reax_system *system, control_params *control,
                               simulation_data *data, storage *workspace,
                               reax_list **lists, output_controls *out_control )
{
  int  i;

  return;

  fprintf( out_control->fnonb, "step: %d\n", data->step );
  fprintf( out_control->fnonb, "%6s\t%-38s%-38s%-38s\n",
           "atom", "nonbonded_total[0,1,2]", "f_vdw[0,1,2]", "f_ele[0,1,2]" );

  for( i = 0; i < system->N; ++i )
    fprintf( out_control->fnonb,
             "%6d%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e\n",
             system->my_atoms[i].orig_id,
             workspace->f_vdw[i][0] + workspace->f_ele[i][0],
             workspace->f_vdw[i][1] + workspace->f_ele[i][1],
             workspace->f_vdw[i][2] + workspace->f_ele[i][2],
             workspace->f_vdw[i][0],
             workspace->f_vdw[i][1],
             workspace->f_vdw[i][2],
             workspace->f_ele[i][0],
             workspace->f_ele[i][1],
             workspace->f_ele[i][2] );
}


void Print_Total_Force( reax_system *system, control_params *control,
                        simulation_data *data, storage *workspace,
                        reax_list **lists, output_controls *out_control )
{
  int    i;

  return;

  fprintf( out_control->ftot, "step: %d\n", data->step );
  fprintf( out_control->ftot, "%6s\t%-38s\n", "atom", "atom.f[0,1,2]");

  for( i = 0; i < system->n; ++i )
    fprintf( out_control->ftot, "%6d%24.15e%24.15e%24.15e\n",
             system->my_atoms[i].orig_id,
             system->my_atoms[i].f[0],
             system->my_atoms[i].f[1],
             system->my_atoms[i].f[2] );
}


void Compare_Total_Forces( reax_system *system, control_params *control,
                           simulation_data *data, storage *workspace,
                           reax_list **lists, output_controls *out_control )
{
  int i;

  return;

  fprintf( out_control->ftot2, "step: %d\n", data->step );
  fprintf( out_control->ftot2, "%6s\t%-38s%-38s\n",
           "atom", "f_total[0,1,2]", "test_force_total[0,1,2]" );

  for( i = 0; i < system->N; ++i )
    fprintf( out_control->ftot2, "%6d%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e\n",
             system->my_atoms[i].orig_id,
             system->my_atoms[i].f[0],
             system->my_atoms[i].f[1],
             system->my_atoms[i].f[2],
             workspace->f_be[i][0] + workspace->f_lp[i][0] +
             workspace->f_ov[i][0] + workspace->f_un[i][0] +
             workspace->f_ang[i][0]+ workspace->f_pen[i][0]+
             workspace->f_coa[i][0]+ + workspace->f_hb[i][0] +
             workspace->f_tor[i][0] + workspace->f_con[i][0] +
             workspace->f_vdw[i][0] + workspace->f_ele[i][0],
             workspace->f_be[i][1] + workspace->f_lp[i][1] +
             workspace->f_ov[i][1] + workspace->f_un[i][1] +
             workspace->f_ang[i][1]+ workspace->f_pen[i][1]+
             workspace->f_coa[i][1]+ + workspace->f_hb[i][1] +
             workspace->f_tor[i][1] + workspace->f_con[i][1] +
             workspace->f_vdw[i][1] + workspace->f_ele[i][1],
             workspace->f_be[i][2] + workspace->f_lp[i][2] +
             workspace->f_ov[i][2] + workspace->f_un[i][2] +
             workspace->f_ang[i][2]+ workspace->f_pen[i][2] +
             workspace->f_coa[i][2]+ + workspace->f_hb[i][2] +
             workspace->f_tor[i][2] + workspace->f_con[i][2] +
             workspace->f_vdw[i][2] + workspace->f_ele[i][2] );
}*/

/*void Init_Force_Test_Functions( )
{
  Print_Interactions[0] = Print_Bond_Orders;
  Print_Interactions[1] = Print_Bond_Forces;
  Print_Interactions[2] = Print_LonePair_Forces;
  Print_Interactions[3] = Print_OverUnderCoor_Forces;
  Print_Interactions[4] = Print_Three_Body_Forces;
  Print_Interactions[5] = Print_Four_Body_Forces;
  Print_Interactions[6] = Print_Hydrogen_Bond_Forces;
  Print_Interactions[7] = Print_vdW_Coulomb_Forces;
  Print_Interactions[8] = Print_Total_Force;
  Print_Interactions[9] = Compare_Total_Forces;
  }*/
