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

#include "reaxff_api.h"

#include "error.h"
#include "utils.h"
#include "tokenizer.h"

#include <cstring>
#include <exception>
#include <string>
#include <unordered_set>

using LAMMPS_NS::utils::getsyserror;
using LAMMPS_NS::utils::sfgets;
using LAMMPS_NS::utils::logmesg;
using LAMMPS_NS::ValueTokenizer;

namespace ReaxFF {
  static std::unordered_set<std::string> ignored_keywords = {
    "ensemble_type", "nsteps", "dt", "proc_by_dim", "random_vel",
    "restart_format", "restart_freq", "reposition_atoms",
    "restrict_bonds", "remove_CoM_vel", "debug_level", "reneighbor",
    "vlist_buffer", "ghost_cutoff", "qeq_freq", "q_err", "ilu_refactor",
    "ilu_droptol", "temp_init", "temp_final", "t_mass", "t_mode", "t_rate",
    "t_freq", "pressure", "p_mass", "pt_mass", "compress", "press_mode",
    "geo_format", "traj_compress", "traj_method", "molecular_analysis",
    "ignore", "dipole_anal", "freq_dipole_anal", "diffusion_coef",
    "freq_diffusion_coef", "restrict_type"
  };

  class parser_error : public std::exception {
    std::string message;
  public:
    parser_error(const std::string &mesg) { message = mesg; }
    const char *what() const noexcept { return message.c_str(); }
  };

  void Read_Control_File(const char *control_file, control_params *control,
                         output_controls *out_control)
  {
    FILE *fp;
    char line[MAX_LINE];
    auto error = control->error_ptr;
    auto lmp = control->lmp_ptr;

    /* open control file */
    fp = fopen(control_file, "r");
    if (!fp)
      error->one(FLERR,fmt::format("The control file {} cannot be opened: {}",
                                   control_file, getsyserror()));
    /* assign default values */
    strcpy(control->sim_name, "simulate");
    control->nthreads = 1;
    control->tabulate = 0;
    control->virial = 0;
    control->bond_cut = 5.0;
    control->bg_cut = 0.3;
    control->thb_cut = 0.001;
    control->thb_cutsq = 0.00001;
    control->hbond_cut = 7.5;

    out_control->write_steps = 0;
    out_control->energy_update_freq = 0;
    strcpy(out_control->traj_title, "default_title");
    out_control->atom_info = 0;
    out_control->bond_info = 0;
    out_control->angle_info = 0;

    /* read control parameters file */
    while (fgets(line, MAX_LINE, fp)) {
      ValueTokenizer values(line);

      // empty line
      if (values.count() == 0) continue;

      try {
        auto keyword = values.next_string();

        if (!values.has_next())
          throw parser_error(fmt::format("No value(s) for control parameter: {}\n",keyword));

        if (ignored_keywords.find(keyword) != ignored_keywords.end()) {
          logmesg(lmp,fmt::format("Ignoring inactive control parameter: {}\n",keyword));
        } else if (keyword == "simulation_name") {
          strcpy(control->sim_name, values.next_string().c_str());
        } else if (keyword == "energy_update_freq") {
          out_control->energy_update_freq = values.next_int();
        } else if (keyword == "nbrhood_cutoff") {
          control->bond_cut = values.next_double();
        } else if (keyword == "bond_graph_cutoff") {
          control->bg_cut = values.next_double();
        } else if (keyword == "thb_cutoff") {
          control->thb_cut = values.next_double();
        } else if (keyword == "thb_cutoff_sq") {
          control->thb_cutsq = values.next_double();
        } else if (keyword == "hbond_cutoff") {
          control->hbond_cut = values.next_double();
        } else if (keyword == "tabulate_long_range") {
          control->tabulate = values.next_int();
        } else if (keyword == "write_freq") {
          out_control->write_steps = values.next_int();
        } else if (keyword == "traj_title") {
          strcpy(out_control->traj_title, values.next_string().c_str());
        } else if (keyword == "atom_info") {
          out_control->atom_info += values.next_int() * 4;
        } else if (keyword == "atom_velocities") {
          out_control->atom_info += values.next_int() * 2;
        } else if (keyword == "atom_forces") {
          out_control->atom_info += values.next_int() * 1;
        } else if (keyword == "bond_info") {
          out_control->bond_info = values.next_int();
        } else if (keyword == "angle_info") {
          out_control->angle_info = values.next_int();
        } else {
          throw parser_error(fmt::format("Unknown parameter {} in control file", keyword));
        }
      } catch (std::exception &e) {
        error->one(FLERR, e.what());
      }
    }
    fclose(fp);
  }
}
