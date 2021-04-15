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
  <https://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "reaxc_traj.h"
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "reaxc_defs.h"
#include "reaxc_list.h"

#include "error.h"

int Reallocate_Output_Buffer( LAMMPS_NS::Error *error_ptr, output_controls *out_control, int req_space )
{
  if (out_control->buffer_len > 0)
    free( out_control->buffer );

  out_control->buffer_len = (int)(req_space*REAX_SAFE_ZONE);
  out_control->buffer = (char*) malloc(out_control->buffer_len*sizeof(char));
  if (out_control->buffer == nullptr) {
    char errmsg[256];
    snprintf(errmsg, 256, "Insufficient memory for required buffer size %d", (int) (req_space*REAX_SAFE_ZONE));
    error_ptr->one(FLERR,errmsg);
  }

  return SUCCESS;
}

void Write_Skip_Line(output_controls *out_control,
                     int my_rank, int skip, int num_section)
{
  if (my_rank == MASTER_NODE)
    fprintf( out_control->strj, INT2_LINE,
             "chars_to_skip_section:", skip, num_section );
}


int Write_Header(reax_system *system, control_params *control,
                 output_controls *out_control)
{
  int  num_hdr_lines, my_hdr_lines, buffer_req;

  char traj_methods[TF_N][10] = { "custom", "xyz" };
  char atom_formats[8][40] =  { "none", "invalid", "invalid", "invalid",
                                "xyz_q",
                                "xyz_q_fxfyfz",
                                "xyz_q_vxvyvz",
                                "detailed_atom_info" };
  char bond_formats[3][30] = { "none",
                               "basic_bond_info",
                               "detailed_bond_info" };
  char angle_formats[2][30] = { "none", "basic_angle_info" };

  /* set header lengths */
  num_hdr_lines = NUM_HEADER_LINES;
  my_hdr_lines = num_hdr_lines * ( system->my_rank == MASTER_NODE );
  buffer_req = my_hdr_lines * HEADER_LINE_LEN;
  if (buffer_req > out_control->buffer_len * DANGER_ZONE)
    Reallocate_Output_Buffer( control->error_ptr, out_control, buffer_req );

  /* only the master node writes into trajectory header */
  if (system->my_rank == MASTER_NODE) {
    /* clear the contents of line & buffer */
    out_control->line[0] = 0;
    out_control->buffer[0] = 0;

    /* to skip the header */
    sprintf( out_control->line, INT_LINE, "chars_to_skip_header:",
             (num_hdr_lines-1) * HEADER_LINE_LEN );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    /* general simulation info */
    sprintf( out_control->line, STR_LINE, "simulation_name:",
             out_control->traj_title );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, BIGINT_LINE, "number_of_atoms:", system->bigN );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, STR_LINE, "ensemble_type:", "(unknown)");
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, INT_LINE, "number_of_steps:", 0);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "timestep_length_(in_fs):", 0.0);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    /* restart info */
    sprintf( out_control->line, STR_LINE, "is_this_a_restart?:", "no");
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, STR_LINE, "write_restart_files?:", "no");
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, INT_LINE, "frequency_to_write_restarts:", 0);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    /* preferences */
    sprintf( out_control->line, STR_LINE, "tabulate_long_range_intrs?:",
             (control->tabulate ? "yes" : "no") );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, INT_LINE, "table_size:", control->tabulate );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, STR_LINE, "restrict_bonds?:", "no");
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, INT_LINE, "bond_restriction_length:", 0);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, STR_LINE, "reposition_atoms?:", "fit to periodic box");
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, INT_LINE, "remove_CoM_velocity?:", 0);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    /* cut-off values */
    sprintf( out_control->line, REAL_LINE, "bonded_intr_dist_cutoff:",
             control->bond_cut );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "nonbonded_intr_dist_cutoff:",
             control->nonb_cut );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "hbond_dist_cutoff:",
             control->hbond_cut );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "reax_bond_threshold:",
             control->bo_cut );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "physical_bond_threshold:",
             control->bg_cut );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "valence_angle_threshold:",
             control->thb_cut );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, SCI_LINE, "QEq_tolerance:", 0.0 );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    /* temperature controls */
    sprintf( out_control->line, REAL_LINE, "initial_temperature:", 0.0);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "target_temperature:", 0.0);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "thermal_inertia:", 0.0);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, STR_LINE, "temperature_regime:", "(unknown)");
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "temperature_change_rate_(K/ps):", 0.0);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    /* pressure controls */
    sprintf( out_control->line, REAL3_LINE, "target_pressure_(GPa):", 0.0, 0.0, 0.0);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL3_LINE, "virial_inertia:", 0.0, 0.0, 0.0);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    /* trajectory */
    sprintf( out_control->line, INT_LINE, "energy_dumping_freq:",
             out_control->energy_update_freq );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, INT_LINE, "trajectory_dumping_freq:",
             out_control->write_steps );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, STR_LINE, "compress_trajectory_output?:",
             (out_control->traj_compress ? "yes" : "no") );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, STR_LINE, "trajectory_format:",
             traj_methods[ out_control->traj_method ] );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, STR_LINE, "atom_info:",
             atom_formats[ out_control->atom_info ] );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, STR_LINE, "bond_info:",
             bond_formats[ out_control->bond_info ] );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, STR_LINE, "angle_info:",
             angle_formats[ out_control->angle_info ] );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    /* analysis */
    sprintf( out_control->line, STR_LINE, "molecular_analysis:", "no");
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );
    sprintf( out_control->line, INT_LINE, "molecular_analysis_frequency:", 0);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );
  }

  /* dump out the buffer */
  if (system->my_rank == MASTER_NODE)
    fprintf( out_control->strj, "%s", out_control->buffer );

  return SUCCESS;
}


int Write_Init_Desc(reax_system *system, output_controls *out_control,
                    MPI_Comm world)
{
  int i, me, np, cnt, buffer_len, buffer_req;
  reax_atom *p_atom;
  MPI_Status status;

  me = system->my_rank;
  np = system->wsize;

  /* skip info */
  Write_Skip_Line(out_control, me, system->bigN*INIT_DESC_LEN, system->bigN);

  if (out_control->traj_method == REG_TRAJ && me == MASTER_NODE)
    buffer_req = system->bigN * INIT_DESC_LEN + 1;
  else buffer_req = system->n * INIT_DESC_LEN + 1;

  if (buffer_req > out_control->buffer_len * DANGER_ZONE)
    Reallocate_Output_Buffer( system->error_ptr, out_control, buffer_req );

  out_control->line[0] = 0;
  out_control->buffer[0] = 0;
  for (i = 0; i < system->n; ++i) {
    p_atom = &( system->my_atoms[i] );
    sprintf( out_control->line, INIT_DESC,
             p_atom->orig_id, p_atom->type, p_atom->name,
             system->reax_param.sbp[ p_atom->type ].mass );
    strncpy( out_control->buffer + i*INIT_DESC_LEN,
             out_control->line, INIT_DESC_LEN+1 );
  }

  if (me != MASTER_NODE) {
    MPI_Send(out_control->buffer, buffer_req-1, MPI_CHAR, MASTER_NODE,
              np * INIT_DESCS + me, world);
  } else {
    buffer_len = system->n * INIT_DESC_LEN;
    for (i = 0; i < np; ++i)
      if (i != MASTER_NODE) {
        MPI_Recv( out_control->buffer + buffer_len, buffer_req - buffer_len,
                  MPI_CHAR, i, np*INIT_DESCS+i, world, &status );
        MPI_Get_count( &status, MPI_CHAR, &cnt );
        buffer_len += cnt;
      }
    out_control->buffer[buffer_len] = 0;
    fprintf( out_control->strj, "%s", out_control->buffer );
  }

  return SUCCESS;
}


int Init_Traj(reax_system *system, control_params *control,
              output_controls *out_control, MPI_Comm world, char *msg)
{
  char fname[MAX_STR+8];
  int  atom_line_len[ NR_OPT_ATOM ] = { 0, 0, 0, 0,
                                        ATOM_BASIC_LEN, ATOM_wV_LEN,
                                        ATOM_wF_LEN, ATOM_FULL_LEN };
  int  bond_line_len[ NR_OPT_BOND ] = { 0, BOND_BASIC_LEN, BOND_FULL_LEN };
  int  angle_line_len[ NR_OPT_ANGLE ] = { 0, ANGLE_BASIC_LEN };

  /* generate trajectory name */
  sprintf( fname, "%s.trj", control->sim_name );

  /* how should I write atoms? */
  out_control->atom_line_len = atom_line_len[ out_control->atom_info ];
  out_control->write_atoms = ( out_control->atom_line_len ? 1 : 0 );
  /* bonds? */
  out_control->bond_line_len = bond_line_len[ out_control->bond_info ];
  out_control->write_bonds = ( out_control->bond_line_len ? 1 : 0 );
  /* angles? */
  out_control->angle_line_len = angle_line_len[ out_control->angle_info ];
  out_control->write_angles = ( out_control->angle_line_len ? 1 : 0 );

  /* allocate line & buffer space */
  out_control->line = (char*) calloc( MAX_TRJ_LINE_LEN + 1, sizeof(char) );
  out_control->buffer_len = 0;
  out_control->buffer = nullptr;

  /* write trajectory header and atom info, if applicable */
  if (out_control->traj_method == REG_TRAJ) {
    if (system->my_rank == MASTER_NODE)
      out_control->strj = fopen(fname, "w");
  } else {
    strcpy( msg, "init_traj: unknown trajectory option" );
    return FAILURE;
  }
  Write_Header(system, control, out_control);
  Write_Init_Desc(system, out_control, world);

  return SUCCESS;
}


int Write_Frame_Header(reax_system *system, control_params *control,
                       simulation_data *data, output_controls *out_control)
{
  int me, num_frm_hdr_lines, my_frm_hdr_lines, buffer_req;

  me = system->my_rank;
  /* frame header lengths */
  num_frm_hdr_lines = 22;
  my_frm_hdr_lines = num_frm_hdr_lines * ( me == MASTER_NODE );
  buffer_req = my_frm_hdr_lines * HEADER_LINE_LEN;
  if (buffer_req > out_control->buffer_len * DANGER_ZONE)
    Reallocate_Output_Buffer( control->error_ptr, out_control, buffer_req );

  /* only the master node writes into trajectory header */
  if (me == MASTER_NODE) {
    /* clear the contents of line & buffer */
    out_control->line[0] = 0;
    out_control->buffer[0] = 0;

    /* skip info */
    sprintf( out_control->line, INT_LINE, "chars_to_skip_frame_header:",
             (num_frm_hdr_lines - 1) * HEADER_LINE_LEN );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    /* step & time */
    sprintf( out_control->line, INT_LINE, "step:", data->step );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "time_in_ps:", 0.0);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    /* box info */
    sprintf( out_control->line, REAL_LINE, "volume:", 0.0 );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL3_LINE, "box_dimensions:", 0.0, 0.0, 0.0);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL3_LINE,
             "coordinate_angles:", 90.0, 90.0, 90.0 );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    /* system T and P */
    sprintf( out_control->line, REAL_LINE, "temperature:", 0.0);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "pressure:", 0.0);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    /* energies */
    double epot = data->sys_en.e_bond + data->sys_en.e_ov + data->sys_en.e_un
      + data->sys_en.e_lp + data->sys_en.e_ang + data->sys_en.e_pen
      + data->sys_en.e_coa + data->sys_en.e_hb + data->sys_en.e_tor
      + data->sys_en.e_con + data->sys_en.e_vdW
      + data->sys_en.e_ele + data->my_en.e_pol;

    sprintf( out_control->line, REAL_LINE, "total_energy:", epot);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "total_kinetic:", 0.0);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "total_potential:", epot);
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "bond_energy:",
             data->sys_en.e_bond );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "atom_energy:",
             data->sys_en.e_ov + data->sys_en.e_un  );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "lone_pair_energy:",
             data->sys_en.e_lp  );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "valence_angle_energy:",
             data->sys_en.e_ang + data->sys_en.e_pen );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "3-body_conjugation:",
             data->sys_en.e_coa );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "hydrogen_bond_energy:",
             data->sys_en.e_hb );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "torsion_angle_energy:",
             data->sys_en.e_tor );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "4-body_conjugation:",
             data->sys_en.e_con );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "vdWaals_energy:",
             data->sys_en.e_vdW );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "electrostatics_energy:",
             data->sys_en.e_ele );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );

    sprintf( out_control->line, REAL_LINE, "polarization_energy:",
             data->sys_en.e_pol );
    strncat( out_control->buffer, out_control->line, HEADER_LINE_LEN+1 );
  }

  /* dump out the buffer */
  if (system->my_rank == MASTER_NODE)
    fprintf( out_control->strj, "%s", out_control->buffer );

  return SUCCESS;
}



int Write_Atoms(reax_system *system, output_controls *out_control,
                MPI_Comm world)
{
  int i, me, np, line_len, buffer_len, buffer_req, cnt;
  MPI_Status  status;
  reax_atom  *p_atom;

  me = system->my_rank;
  np = system->wsize;
  line_len = out_control->atom_line_len;

  Write_Skip_Line(out_control, me, system->bigN*line_len, system->bigN);

  if (out_control->traj_method == REG_TRAJ && me == MASTER_NODE)
    buffer_req = system->bigN * line_len + 1;
  else buffer_req = system->n * line_len + 1;

  if (buffer_req > out_control->buffer_len * DANGER_ZONE)
    Reallocate_Output_Buffer( system->error_ptr, out_control, buffer_req );

  /* fill in buffer */
  out_control->line[0] = 0;
  out_control->buffer[0] = 0;
  for (i = 0; i < system->n; ++i) {
    p_atom = &( system->my_atoms[i] );

    switch (out_control->atom_info) {
    case OPT_ATOM_BASIC:
      sprintf( out_control->line, ATOM_BASIC,
               p_atom->orig_id, p_atom->x[0], p_atom->x[1], p_atom->x[2],
               p_atom->q );
      break;
    case OPT_ATOM_wF:
      sprintf( out_control->line, ATOM_wF,
               p_atom->orig_id, p_atom->x[0], p_atom->x[1], p_atom->x[2],
               p_atom->f[0], p_atom->f[1], p_atom->f[2], p_atom->q );
      break;
    case OPT_ATOM_wV:
      sprintf( out_control->line, ATOM_wV,
               p_atom->orig_id, p_atom->x[0], p_atom->x[1], p_atom->x[2],
               p_atom->v[0], p_atom->v[1], p_atom->v[2], p_atom->q );
      break;
    case OPT_ATOM_FULL:
      sprintf( out_control->line, ATOM_FULL,
               p_atom->orig_id, p_atom->x[0], p_atom->x[1], p_atom->x[2],
               p_atom->v[0], p_atom->v[1], p_atom->v[2],
               p_atom->f[0], p_atom->f[1], p_atom->f[2], p_atom->q );
      break;
    default:
      system->error_ptr->one(FLERR,"Write_traj_atoms: unknown atom trajectory format");
    }

    strncpy( out_control->buffer + i*line_len, out_control->line, line_len+1 );
  }

  if (me != MASTER_NODE) {
    MPI_Send( out_control->buffer, buffer_req-1, MPI_CHAR, MASTER_NODE,
              np*ATOM_LINES+me, world);
  } else {
    buffer_len = system->n * line_len;
    for (i = 0; i < np; ++i)
      if (i != MASTER_NODE) {
        MPI_Recv( out_control->buffer + buffer_len, buffer_req - buffer_len,
                  MPI_CHAR, i, np*ATOM_LINES+i, world, &status );
        MPI_Get_count( &status, MPI_CHAR, &cnt );
        buffer_len += cnt;
      }
    out_control->buffer[buffer_len] = 0;
    fprintf( out_control->strj, "%s", out_control->buffer );
  }

  return SUCCESS;
}


int Write_Bonds(reax_system *system, control_params *control, reax_list *bonds,
                output_controls *out_control, MPI_Comm world)
{
  int i, j, pj, me, np;
  int my_bonds, num_bonds;
  int line_len, buffer_len, buffer_req, cnt;
  MPI_Status  status;
  bond_data  *bo_ij;

  me = system->my_rank;
  np = system->wsize;
  line_len = out_control->bond_line_len;

  /* count the number of bonds I will write */
  my_bonds = 0;
  for (i=0; i < system->n; ++i)
    for (pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj) {
      j = bonds->select.bond_list[pj].nbr;
      if ( system->my_atoms[i].orig_id <= system->my_atoms[j].orig_id &&
          bonds->select.bond_list[pj].bo_data.BO >= control->bg_cut )
        ++my_bonds;
    }

  /* allreduce - total number of bonds */
  MPI_Allreduce( &my_bonds, &num_bonds, 1, MPI_INT, MPI_SUM, world );

  Write_Skip_Line( out_control, me, num_bonds*line_len, num_bonds );

  if (out_control->traj_method == REG_TRAJ && me == MASTER_NODE)
    buffer_req = num_bonds * line_len + 1;
  else buffer_req = my_bonds * line_len + 1;

  if (buffer_req > out_control->buffer_len * DANGER_ZONE)
    Reallocate_Output_Buffer( system->error_ptr, out_control, buffer_req );

  /* fill in the buffer */
  out_control->line[0] = 0;
  out_control->buffer[0] = 0;

  my_bonds = 0;
  for (i=0; i < system->n; ++i) {
    for (pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj) {
      bo_ij = &( bonds->select.bond_list[pj] );
      j = bo_ij->nbr;

      if ( system->my_atoms[i].orig_id <= system->my_atoms[j].orig_id &&
          bo_ij->bo_data.BO >= control->bg_cut) {
        switch (out_control->bond_info) {
        case OPT_BOND_BASIC:
          sprintf( out_control->line, BOND_BASIC,
                   system->my_atoms[i].orig_id, system->my_atoms[j].orig_id,
                   bo_ij->d, bo_ij->bo_data.BO );
          break;
        case OPT_BOND_FULL:
          sprintf( out_control->line, BOND_FULL,
                   system->my_atoms[i].orig_id, system->my_atoms[j].orig_id,
                   bo_ij->d, bo_ij->bo_data.BO, bo_ij->bo_data.BO_s,
                   bo_ij->bo_data.BO_pi, bo_ij->bo_data.BO_pi2 );
          break;
        default:
          system->error_ptr->one(FLERR, "Write_traj_bonds: FATAL! invalid bond_info option");
        }
        strncpy( out_control->buffer + my_bonds*line_len,
                 out_control->line, line_len+1 );
        ++my_bonds;
      }
    }
  }

  if (me != MASTER_NODE) {
    MPI_Send( out_control->buffer, buffer_req-1, MPI_CHAR, MASTER_NODE,
              np*BOND_LINES+me, world );
  } else {
    buffer_len = my_bonds * line_len;
    for (i = 0; i < np; ++i)
      if (i != MASTER_NODE) {
        MPI_Recv( out_control->buffer + buffer_len, buffer_req - buffer_len,
                  MPI_CHAR, i, np*BOND_LINES+i, world, &status );
        MPI_Get_count( &status, MPI_CHAR, &cnt );
        buffer_len += cnt;
      }
    out_control->buffer[buffer_len] = 0;
    fprintf( out_control->strj, "%s", out_control->buffer );
  }

  return SUCCESS;
}


int Write_Angles( reax_system *system, control_params *control,
                  reax_list *bonds, reax_list *thb_intrs,
                  output_controls *out_control, MPI_Comm world)
{
  int i, j, k, pi, pk, me, np;
  int my_angles, num_angles;
  int line_len, buffer_len, buffer_req, cnt;
  bond_data  *bo_ij, *bo_jk;
  three_body_interaction_data *angle_ijk;
  MPI_Status  status;

  me = system->my_rank;
  np = system->wsize;
  line_len = out_control->angle_line_len;

  /* count the number of valence angles I will output */
  my_angles = 0;
  for (j = 0; j < system->n; ++j)
    for (pi = Start_Index(j, bonds); pi < End_Index(j, bonds); ++pi) {
      bo_ij = &(bonds->select.bond_list[pi]);
      i     = bo_ij->nbr;

      if (bo_ij->bo_data.BO >= control->bg_cut) // physical j&i bond
        for (pk = Start_Index( pi, thb_intrs);
             pk < End_Index( pi, thb_intrs ); ++pk) {
          angle_ijk = &(thb_intrs->select.three_body_list[pk]);
          k       = angle_ijk->thb;
          bo_jk   = &(bonds->select.bond_list[ angle_ijk->pthb ]);

          if ( system->my_atoms[i].orig_id < system->my_atoms[k].orig_id &&
              bo_jk->bo_data.BO >= control->bg_cut ) // physical j&k bond
            ++my_angles;
        }
    }
  /* total number of valences */
  MPI_Allreduce(&my_angles, &num_angles, 1, MPI_INT, MPI_SUM,  world);

  Write_Skip_Line( out_control, me, num_angles*line_len, num_angles );

  if (out_control->traj_method == REG_TRAJ && me == MASTER_NODE)
    buffer_req = num_angles * line_len + 1;
  else buffer_req = my_angles * line_len + 1;

  if (buffer_req > out_control->buffer_len * DANGER_ZONE)
    Reallocate_Output_Buffer( system->error_ptr, out_control, buffer_req );

  /* fill in the buffer */
  my_angles = 0;
  out_control->line[0] = 0;
  out_control->buffer[0] = 0;
  for (j = 0; j < system->n; ++j)
    for (pi = Start_Index(j, bonds); pi < End_Index(j, bonds); ++pi) {
      bo_ij = &(bonds->select.bond_list[pi]);
      i     = bo_ij->nbr;

      if (bo_ij->bo_data.BO >= control->bg_cut) // physical j&i bond
        for (pk = Start_Index( pi, thb_intrs);
             pk < End_Index( pi, thb_intrs ); ++pk) {
          angle_ijk = &(thb_intrs->select.three_body_list[pk]);
          k       = angle_ijk->thb;
          bo_jk   = &(bonds->select.bond_list[ angle_ijk->pthb ]);

          if ( system->my_atoms[i].orig_id < system->my_atoms[k].orig_id &&
              bo_jk->bo_data.BO >= control->bg_cut) { // physical j&k bond
            sprintf( out_control->line, ANGLE_BASIC,
                     system->my_atoms[i].orig_id, system->my_atoms[j].orig_id,
                     system->my_atoms[k].orig_id, RAD2DEG( angle_ijk->theta ) );

            strncpy( out_control->buffer + my_angles*line_len,
                     out_control->line, line_len+1 );
            ++my_angles;
          }
        }
    }

  if (me != MASTER_NODE) {
    MPI_Send( out_control->buffer, buffer_req-1, MPI_CHAR, MASTER_NODE,
              np*ANGLE_LINES+me, world );
  } else {
    buffer_len = my_angles * line_len;
    for (i = 0; i < np; ++i)
      if (i != MASTER_NODE) {
        MPI_Recv( out_control->buffer + buffer_len, buffer_req - buffer_len,
                  MPI_CHAR, i, np*ANGLE_LINES+i, world, &status );
        MPI_Get_count( &status, MPI_CHAR, &cnt );
        buffer_len += cnt;
      }
    out_control->buffer[buffer_len] = 0;
    fprintf( out_control->strj, "%s", out_control->buffer );
  }

  return SUCCESS;
}


int Append_Frame(reax_system *system, control_params *control,
                 simulation_data *data, reax_list **lists,
                 output_controls *out_control, MPI_Comm world)
{
  Write_Frame_Header(system, control, data, out_control);

  if (out_control->write_atoms)
    Write_Atoms(system, out_control, world);

  if (out_control->write_bonds)
    Write_Bonds(system, control, (*lists + BONDS), out_control, world);

  if (out_control->write_angles)
    Write_Angles(system, control, (*lists + BONDS), (*lists + THREE_BODIES),
                  out_control, world);

  return SUCCESS;
}


int End_Traj(int my_rank, output_controls *out_control)
{
  if (my_rank == MASTER_NODE)
    fclose(out_control->strj);

  free(out_control->buffer);
  free(out_control->line);

  return SUCCESS;
}
