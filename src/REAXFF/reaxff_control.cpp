// clang-format off
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
#include "text_file_reader.h"

#include <cstring>
#include <exception>
#include <string>
#include <unordered_set>

using LAMMPS_NS::utils::getsyserror;
using LAMMPS_NS::utils::sfgets;
using LAMMPS_NS::utils::logmesg;
using LAMMPS_NS::ValueTokenizer;

namespace ReaxFF {
  static std::unordered_set<std::string> inactive_keywords = {
    "ensemble_type", "nsteps", "dt", "proc_by_dim", "random_vel",
    "restart_format", "restart_freq", "reposition_atoms",
    "restrict_bonds", "remove_CoM_vel", "debug_level", "reneighbor",
    "vlist_buffer", "ghost_cutoff", "qeq_freq", "q_err", "ilu_refactor",
    "ilu_droptol", "temp_init", "temp_final", "t_mass", "t_mode", "t_rate",
    "t_freq", "pressure", "p_mass", "pt_mass", "compress", "press_mode",
    "geo_format", "traj_compress", "traj_method", "molecular_analysis",
    "ignore", "dipole_anal", "freq_dipole_anal", "diffusion_coef",
    "freq_diffusion_coef", "restrict_type", "traj_title", "simulation_name",
    "energy_update_freq", "atom_info", "atom_velocities", "atom_forces",
    "bond_info", "angle_info" };

  class parser_error : public std::exception {
    std::string message;
  public:
    parser_error(const std::string &mesg) { message = mesg; }
    const char *what() const noexcept { return message.c_str(); }
  };

  // NOTE: this function is run on MPI rank 0 only

  void Read_Control_File(const char *control_file, control_params *control)
  {
    auto error = control->error_ptr;

    /* assign default values */
    control->nthreads = 1;
    control->tabulate = 0;
    control->bond_cut = 5.0;
    control->bg_cut = 0.3;
    control->thb_cut = 0.001;
    control->thb_cutsq = 0.00001;
    control->hbond_cut = 7.5;

    /* read control parameters file */
    try {
      LAMMPS_NS::TextFileReader reader(control_file, "ReaxFF control");
      reader.ignore_comments = false;

      while (1) {
        auto values = reader.next_values(0);

        // empty line
        if (values.count() == 0) continue;

        auto keyword = values.next_string();

        if (!values.has_next())
          throw parser_error(fmt::format("No value(s) for control parameter: {}\n",keyword));

        if (inactive_keywords.find(keyword) != inactive_keywords.end()) {
          error->warning(FLERR,fmt::format("Ignoring inactive control "
                                           "parameter: {}",keyword));
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
          if (values.next_int() > 0)
            error->warning(FLERR,"Support for writing native trajectories has "
                           "been removed after LAMMPS version 8 April 2021");
        } else {
          throw parser_error(fmt::format("Unknown parameter {} in "
                                         "control file", keyword));
        }
      }
    } catch (LAMMPS_NS::EOFException &) {
      ; // catch and ignore
    } catch (std::exception &e) {
      error->one(FLERR, e.what());
    }
  }
}
