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

#define MAX_TRJ_LINE_LEN     120
#define MAX_TRJ_BUFFER_SIZE  (MAX_TRJ_LINE_LEN * 100)

#define NUM_HEADER_LINES 37
#define HEADER_LINE_LEN 62
#define INIT_DESC_LEN 32

#define ATOM_BASIC_LEN 50
#define ATOM_wV_LEN 80
#define ATOM_wF_LEN 80
#define ATOM_FULL_LEN 110

#define BOND_BASIC_LEN 39
#define BOND_FULL_LEN 69
#define ANGLE_BASIC_LEN 38

enum ATOM_LINE_OPTS  { OPT_NOATOM = 0, OPT_ATOM_BASIC = 4, OPT_ATOM_wF = 5, OPT_ATOM_wV = 6, OPT_ATOM_FULL = 7, NR_OPT_ATOM = 8 };
enum BOND_LINE_OPTS  { OPT_NOBOND, OPT_BOND_BASIC, OPT_BOND_FULL, NR_OPT_BOND };
enum ANGLE_LINE_OPTS { OPT_NOANGLE, OPT_ANGLE_BASIC, NR_OPT_ANGLE };

std::string fmtline(const char *key, double val)
{
  return fmt::format("{:<37}{:<24.3f}\n",key,val);
}

std::string fmtline(const char *key, double val1, double val2, double val3)
{
  return fmt::format("{:<32}{:>9.3f},{:>9.3f},{:>9.3f}\n",key,val1,val2,val3);
}

template <typename T>
std::string fmtline(const char *key, T val)
{
  return fmt::format("{:<37}{:<24}\n",key,val);
}

template <typename T>
std::string fmtline(const char *key, T val1, T val2)
{
  return fmt::format("{:<36}{:<12},{:<12}\n",key,val1,val2);
}

void Reallocate_Output_Buffer(LAMMPS_NS::Error *error_ptr, output_controls *out_control, int req_space)
{
  if (out_control->buffer_len > 0)
    free(out_control->buffer);

  out_control->buffer_len = (int)(req_space*REAX_SAFE_ZONE);
  out_control->buffer = (char*) malloc(out_control->buffer_len*sizeof(char));
  if (!out_control->buffer)
    error_ptr->one(FLERR,fmt::format("Insufficient memory for required buffer "
                                     "size {}", (req_space*REAX_SAFE_ZONE)));
}

void Write_Skip_Line(output_controls *out_control, int my_rank, int skip, int num_section)
{
  if (my_rank == MASTER_NODE)
    fmt::print(out_control->strj,fmtline("chars_to_skip_section:",skip,num_section));
}

void Write_Header(reax_system *system, control_params *control, output_controls *out_control)
{
  char atom_formats[8][40] =  { "none", "invalid", "invalid", "invalid",
                                "xyz_q",
                                "xyz_q_fxfyfz",
                                "xyz_q_vxvyvz",
                                "detailed_atom_info" };
  char bond_formats[3][30] = { "none",
                               "basic_bond_info",
                               "detailed_bond_info" };
  char angle_formats[2][30] = { "none", "basic_angle_info" };

  /* only the master node writes into trajectory header */
  if (system->my_rank == MASTER_NODE) {
    std::string buffer("");

    /* to skip the header */
    buffer += fmtline("chars_to_skip_header:",(NUM_HEADER_LINES-1) * HEADER_LINE_LEN);

    /* general simulation info */
    buffer += fmtline("simulation_name:", out_control->traj_title);
    buffer += fmtline("number_of_atoms:", system->bigN);
    buffer += fmtline("ensemble_type:", "(unknown)");
    buffer += fmtline("number_of_steps:", 0);
    buffer += fmtline("timestep_length_(in_fs):", 0.0);

    /* restart info */
    buffer += fmtline("is_this_a_restart?:", "no");
    buffer += fmtline("write_restart_files?:", "no");
    buffer += fmtline("frequency_to_write_restarts:", 0);

    /* preferences */
    buffer += fmtline("tabulate_long_range_intrs?:",(control->tabulate ? "yes" : "no"));
    buffer += fmtline("table_size:", control->tabulate);
    buffer += fmtline("restrict_bonds?:", "no");
    buffer += fmtline("bond_restriction_length:", 0);
    buffer += fmtline("reposition_atoms?:", "fit to periodic box");
    buffer += fmtline("remove_CoM_velocity?:", 0);
    buffer += fmtline("bonded_intr_dist_cutoff:",control->bond_cut);
    buffer += fmtline("nonbonded_intr_dist_cutoff:",control->nonb_cut);
    buffer += fmtline("hbond_dist_cutoff:",control->hbond_cut);
    buffer += fmtline("reax_bond_threshold:",control->bo_cut);
    buffer += fmtline("physical_bond_threshold:", control->bg_cut);
    buffer += fmtline("valence_angle_threshold:",control->thb_cut);
    buffer += fmtline("QEq_tolerance:", 0.0);

    /* temperature controls */
    buffer += fmtline("initial_temperature:", 0.0);
    buffer += fmtline("target_temperature:", 0.0);
    buffer += fmtline("thermal_inertia:", 0.0);
    buffer += fmtline("temperature_regime:", "(unknown)");
    buffer += fmtline("temperature_change_rate_(K/ps):", 0.0);

    /* pressure controls */
    buffer += fmtline("target_pressure_(GPa):", 0.0, 0.0, 0.0);
    buffer += fmtline("virial_inertia:", 0.0, 0.0, 0.0);

    /* trajectory */
    buffer += fmtline("energy_dumping_freq:", out_control->energy_update_freq);
    buffer += fmtline("trajectory_dumping_freq:",out_control->write_steps);
    buffer += fmtline("compress_trajectory_output?:", "no");
    buffer += fmtline("trajectory_format:", "regular");
    buffer += fmtline("atom_info:", atom_formats[out_control->atom_info]);
    buffer += fmtline("bond_info:", bond_formats[out_control->bond_info]);
    buffer += fmtline("angle_info:", angle_formats[out_control->angle_info]);

    /* analysis */
    buffer += fmtline("molecular_analysis:", "no");
    buffer += fmtline("molecular_analysis_frequency:", 0);

    /* dump out the buffer */

    fputs(buffer.c_str(),out_control->strj);
  }
}

void Write_Init_Desc(reax_system *system, output_controls *out_control, MPI_Comm world)
{
  int i, me, np, cnt, buffer_len, buffer_req;
  reax_atom *p_atom;
  MPI_Status status;

  me = system->my_rank;
  np = system->wsize;

  /* skip info */
  Write_Skip_Line(out_control, me, system->bigN*INIT_DESC_LEN, system->bigN);

  if (me == MASTER_NODE)
    buffer_req = system->bigN * INIT_DESC_LEN + 1;
  else buffer_req = system->n * INIT_DESC_LEN + 1;

  if (buffer_req > out_control->buffer_len * DANGER_ZONE)
    Reallocate_Output_Buffer(system->error_ptr, out_control, buffer_req);

  out_control->buffer[0] = 0;
  for (i = 0; i < system->n; ++i) {
    p_atom = &(system->my_atoms[i]);
    auto buffer = fmt::format("{:9}{:3}{:9}{:10.3f}\n",
                              p_atom->orig_id, p_atom->type, p_atom->name,
                              system->reax_param.sbp[p_atom->type].mass);
    strncpy(out_control->buffer + i*INIT_DESC_LEN, buffer.c_str(), INIT_DESC_LEN+1);
  }

  if (me != MASTER_NODE) {
    MPI_Send(out_control->buffer, buffer_req-1, MPI_CHAR, MASTER_NODE,
             np * INIT_DESCS + me, world);
  } else {
    buffer_len = system->n * INIT_DESC_LEN;
    for (i = 0; i < np; ++i)
      if (i != MASTER_NODE) {
        MPI_Recv(out_control->buffer + buffer_len, buffer_req - buffer_len,
                 MPI_CHAR, i, np*INIT_DESCS+i, world, &status);
        MPI_Get_count(&status, MPI_CHAR, &cnt);
        buffer_len += cnt;
      }
    out_control->buffer[buffer_len] = 0;
    fputs(out_control->buffer,out_control->strj);
  }
}

void Init_Traj(reax_system *system, control_params *control,
               output_controls *out_control, MPI_Comm world)
{
  int atom_line_len[NR_OPT_ATOM] = { 0, 0, 0, 0, ATOM_BASIC_LEN, ATOM_wV_LEN, ATOM_wF_LEN, ATOM_FULL_LEN};
  int bond_line_len[NR_OPT_BOND] = { 0, BOND_BASIC_LEN, BOND_FULL_LEN };
  int angle_line_len[NR_OPT_ANGLE] = { 0, ANGLE_BASIC_LEN };

  /* generate trajectory name */
  auto fname = std::string(control->sim_name) + ".trj";

  /* how should I write atoms? */
  out_control->atom_line_len = atom_line_len[out_control->atom_info];
  out_control->write_atoms = (out_control->atom_line_len ? 1 : 0);
  /* bonds? */
  out_control->bond_line_len = bond_line_len[out_control->bond_info];
  out_control->write_bonds = (out_control->bond_line_len ? 1 : 0);
  /* angles? */
  out_control->angle_line_len = angle_line_len[out_control->angle_info];
  out_control->write_angles = (out_control->angle_line_len ? 1 : 0);

  /* allocate line & buffer space */
  out_control->buffer_len = 0;
  out_control->buffer = nullptr;

  /* write trajectory header and atom info, if applicable */
  if (system->my_rank == MASTER_NODE)
    out_control->strj = fopen(fname.c_str(), "w");

  Write_Header(system, control, out_control);
  Write_Init_Desc(system, out_control, world);
}

void Write_Frame_Header(reax_system *system, simulation_data *data, output_controls *out_control)
{
  const int me = system->my_rank;
  /* frame header lengths */
  constexpr int num_frm_hdr_lines = 22;

  /* only the master node writes into trajectory header */
  if (me == MASTER_NODE) {
    std::string buffer("");

    /* skip info */
    buffer += fmtline("chars_to_skip_frame_header:", (num_frm_hdr_lines-1)*HEADER_LINE_LEN);

    /* step & time */
    buffer += fmtline("step:", data->step);
    buffer += fmtline("time_in_ps:", 0.0);

    /* box info */
    buffer += fmtline("volume:", 0.0);
    buffer += fmtline("box_dimensions:", 0.0, 0.0, 0.0);
    buffer += fmtline("coordinate_angles:", 90.0, 90.0, 90.0);

    /* system T and P */
    buffer += fmtline("temperature:", 0.0);
    buffer += fmtline("pressure:", 0.0);

    /* energies */
    double epot = data->sys_en.e_bond + data->sys_en.e_ov + data->sys_en.e_un
      + data->sys_en.e_lp + data->sys_en.e_ang + data->sys_en.e_pen
      + data->sys_en.e_coa + data->sys_en.e_hb + data->sys_en.e_tor
      + data->sys_en.e_con + data->sys_en.e_vdW
      + data->sys_en.e_ele + data->my_en.e_pol;

    buffer += fmtline("total_energy:", epot);
    buffer += fmtline("total_kinetic:", 0.0);
    buffer += fmtline("total_potential:", epot);
    buffer += fmtline("bond_energy:", data->sys_en.e_bond);
    buffer += fmtline("atom_energy:", data->sys_en.e_ov + data->sys_en.e_un);
    buffer += fmtline("lone_pair_energy:", data->sys_en.e_lp);
    buffer += fmtline("valence_angle_energy:", data->sys_en.e_ang + data->sys_en.e_pen);
    buffer += fmtline("3-body_conjugation:", data->sys_en.e_coa);
    buffer += fmtline("hydrogen_bond_energy:", data->sys_en.e_hb);
    buffer += fmtline("torsion_angle_energy:", data->sys_en.e_tor);
    buffer += fmtline("4-body_conjugation:", data->sys_en.e_con);
    buffer += fmtline("vdWaals_energy:", data->sys_en.e_vdW);
    buffer += fmtline("electrostatics_energy:", data->sys_en.e_ele);
    buffer += fmtline("polarization_energy:", data->sys_en.e_pol);

    /* dump out the buffer */
    fputs(buffer.c_str(),out_control->strj);
  }
}

void Write_Atoms(reax_system *system, output_controls *out_control,
                 MPI_Comm world)
{
  int i, me, np, line_len, buffer_len, buffer_req, cnt;
  MPI_Status  status;
  reax_atom  *p_atom;

  me = system->my_rank;
  np = system->wsize;
  line_len = out_control->atom_line_len;

  Write_Skip_Line(out_control, me, system->bigN*line_len, system->bigN);

  if (me == MASTER_NODE)
    buffer_req = system->bigN * line_len + 1;
  else buffer_req = system->n * line_len + 1;

  if (buffer_req > out_control->buffer_len * DANGER_ZONE)
    Reallocate_Output_Buffer(system->error_ptr, out_control, buffer_req);

  /* fill in buffer */
  out_control->buffer[0] = 0;

  for (i = 0; i < system->n; ++i) {
    std::string buffer("");
    p_atom = &(system->my_atoms[i]);
    buffer += fmt::format("{:9}{:10.3f}{:10.3f}{:10.3f}",p_atom->orig_id,
                          p_atom->x[0], p_atom->x[1],p_atom->x[2]);

    switch (out_control->atom_info) {
    case OPT_ATOM_BASIC:
      buffer += fmt::format("{:10.3f}\n",p_atom->q);
      break;
    case OPT_ATOM_wF:
      buffer += fmt::format("{:10.3f}{:10.3f}{:10.3f}{:10.3f}\n",
                            p_atom->f[0], p_atom->f[1], p_atom->f[2], p_atom->q);
      break;
    case OPT_ATOM_wV:
      buffer += fmt::format("{:10.3f}{:10.3f}{:10.3f}{:10.3f}\n",
                            p_atom->v[0], p_atom->v[1], p_atom->v[2], p_atom->q);
      break;
    case OPT_ATOM_FULL:
      buffer += fmt::format("{:10.3f}{:10.3f}{:10.3f}{:10.3f}{:10.3f}{:10.3f}\n",
                            p_atom->v[0], p_atom->v[1], p_atom->v[2],
                            p_atom->f[0], p_atom->f[1], p_atom->f[2], p_atom->q);
      break;
    default:
      system->error_ptr->one(FLERR,"Write_traj_atoms: unknown atom trajectory format");
    }

    strncpy(out_control->buffer + i*line_len, buffer.c_str(), line_len+1);
  }

  if (me != MASTER_NODE) {
    MPI_Send(out_control->buffer, buffer_req-1, MPI_CHAR, MASTER_NODE,
             np*ATOM_LINES+me, world);
  } else {
    buffer_len = system->n * line_len;
    for (i = 0; i < np; ++i)
      if (i != MASTER_NODE) {
        MPI_Recv(out_control->buffer + buffer_len, buffer_req - buffer_len,
                 MPI_CHAR, i, np*ATOM_LINES+i, world, &status);
        MPI_Get_count(&status, MPI_CHAR, &cnt);
        buffer_len += cnt;
      }
    out_control->buffer[buffer_len] = 0;
    fputs(out_control->buffer, out_control->strj);
  }
}

void Write_Bonds(reax_system *system, control_params *control, reax_list *bonds,
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
      if (system->my_atoms[i].orig_id <= system->my_atoms[j].orig_id &&
          bonds->select.bond_list[pj].bo_data.BO >= control->bg_cut)
        ++my_bonds;
    }

  /* allreduce - total number of bonds */
  MPI_Allreduce(&my_bonds, &num_bonds, 1, MPI_INT, MPI_SUM, world);

  Write_Skip_Line(out_control, me, num_bonds*line_len, num_bonds);

  if (me == MASTER_NODE)
    buffer_req = num_bonds * line_len + 1;
  else buffer_req = my_bonds * line_len + 1;

  if (buffer_req > out_control->buffer_len * DANGER_ZONE)
    Reallocate_Output_Buffer(system->error_ptr, out_control, buffer_req);

  /* fill in the buffer */
  out_control->buffer[0] = 0;

  my_bonds = 0;
  for (i=0; i < system->n; ++i) {
    for (pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj) {
      bo_ij = &(bonds->select.bond_list[pj]);
      j = bo_ij->nbr;

      if (system->my_atoms[i].orig_id <= system->my_atoms[j].orig_id &&
          bo_ij->bo_data.BO >= control->bg_cut) {
        auto buffer = fmt::format("{:9}{:9}{:10.3f}{:10.3f}", system->my_atoms[i].orig_id,
                                  system->my_atoms[j].orig_id,bo_ij->d,bo_ij->bo_data.BO);

        switch (out_control->bond_info) {
        case OPT_BOND_BASIC:
          buffer += "\n";
          break;
        case OPT_BOND_FULL:
          buffer += fmt::format("{:10.3f}{:10.3f}{:10.3f}\n", bo_ij->bo_data.BO_s,
                                bo_ij->bo_data.BO_pi, bo_ij->bo_data.BO_pi2);
          break;
        default:
          system->error_ptr->one(FLERR, "Write_traj_bonds: FATAL! invalid bond_info option");
        }
        strncpy(out_control->buffer + my_bonds*line_len, buffer.c_str(), line_len+1);
        ++my_bonds;
      }
    }
  }

  if (me != MASTER_NODE) {
    MPI_Send(out_control->buffer, buffer_req-1, MPI_CHAR, MASTER_NODE, np*BOND_LINES+me, world);
  } else {
    buffer_len = my_bonds * line_len;
    for (i = 0; i < np; ++i)
      if (i != MASTER_NODE) {
        MPI_Recv(out_control->buffer + buffer_len, buffer_req - buffer_len,
                 MPI_CHAR, i, np*BOND_LINES+i, world, &status);
        MPI_Get_count(&status, MPI_CHAR, &cnt);
        buffer_len += cnt;
      }
    out_control->buffer[buffer_len] = 0;
    fputs(out_control->buffer,out_control->strj);
  }
}

void Write_Angles(reax_system *system, control_params *control,
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
        for (pk = Start_Index(pi, thb_intrs);
             pk < End_Index(pi, thb_intrs); ++pk) {
          angle_ijk = &(thb_intrs->select.three_body_list[pk]);
          k       = angle_ijk->thb;
          bo_jk   = &(bonds->select.bond_list[ angle_ijk->pthb ]);

          if (system->my_atoms[i].orig_id < system->my_atoms[k].orig_id &&
              bo_jk->bo_data.BO >= control->bg_cut) // physical j&k bond
            ++my_angles;
        }
    }
  /* total number of valences */
  MPI_Allreduce(&my_angles, &num_angles, 1, MPI_INT, MPI_SUM,  world);

  Write_Skip_Line(out_control, me, num_angles*line_len, num_angles);

  if (me == MASTER_NODE)
    buffer_req = num_angles * line_len + 1;
  else buffer_req = my_angles * line_len + 1;

  if (buffer_req > out_control->buffer_len * DANGER_ZONE)
    Reallocate_Output_Buffer(system->error_ptr, out_control, buffer_req);

  /* fill in the buffer */
  my_angles = 0;
  out_control->buffer[0] = 0;
  for (j = 0; j < system->n; ++j)
    for (pi = Start_Index(j, bonds); pi < End_Index(j, bonds); ++pi) {
      bo_ij = &(bonds->select.bond_list[pi]);
      i     = bo_ij->nbr;

      if (bo_ij->bo_data.BO >= control->bg_cut) // physical j&i bond
        for (pk = Start_Index(pi, thb_intrs);
             pk < End_Index(pi, thb_intrs); ++pk) {
          angle_ijk = &(thb_intrs->select.three_body_list[pk]);
          k       = angle_ijk->thb;
          bo_jk   = &(bonds->select.bond_list[ angle_ijk->pthb ]);

          if (system->my_atoms[i].orig_id < system->my_atoms[k].orig_id &&
              bo_jk->bo_data.BO >= control->bg_cut) { // physical j&k bond
            auto buffer = fmt::format("{:9}{:9}{:9}{:10.3f}\n",
                                      system->my_atoms[i].orig_id, system->my_atoms[j].orig_id,
                                      system->my_atoms[k].orig_id, RAD2DEG(angle_ijk->theta));

            strncpy(out_control->buffer + my_angles*line_len, buffer.c_str(), line_len+1);
            ++my_angles;
          }
        }
    }

  if (me != MASTER_NODE) {
    MPI_Send(out_control->buffer, buffer_req-1, MPI_CHAR, MASTER_NODE,
             np*ANGLE_LINES+me, world);
  } else {
    buffer_len = my_angles * line_len;
    for (i = 0; i < np; ++i)
      if (i != MASTER_NODE) {
        MPI_Recv(out_control->buffer + buffer_len, buffer_req - buffer_len,
                 MPI_CHAR, i, np*ANGLE_LINES+i, world, &status);
        MPI_Get_count(&status, MPI_CHAR, &cnt);
        buffer_len += cnt;
      }
    out_control->buffer[buffer_len] = 0;
    fputs(out_control->buffer,out_control->strj);
  }
}

void Append_Frame(reax_system *system, control_params *control,
                  simulation_data *data, reax_list **lists,
                  output_controls *out_control, MPI_Comm world)
{
  Write_Frame_Header(system, data, out_control);

  if (out_control->write_atoms)
    Write_Atoms(system, out_control, world);

  if (out_control->write_bonds)
    Write_Bonds(system, control, (*lists + BONDS), out_control, world);

  if (out_control->write_angles)
    Write_Angles(system, control, (*lists + BONDS), (*lists + THREE_BODIES),
                 out_control, world);
}

void End_Traj(int my_rank, output_controls *out_control)
{
  if (my_rank == MASTER_NODE)
    fclose(out_control->strj);

  free(out_control->buffer);
}
