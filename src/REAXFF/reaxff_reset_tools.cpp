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

#include <cstring>

namespace ReaxFF {

  static void Reset_Atoms(reax_system* system, control_params *control)
  {
    int i;
    reax_atom *atom;

    system->numH = 0;
    if (control->hbond_cut > 0)
      for (i = 0; i < system->n; ++i) {
        atom = &(system->my_atoms[i]);
        if (atom->type < 0) continue;
        if (system->reax_param.sbp[atom->type].p_hbond == 1)
          atom->Hindex = system->numH++;
        else atom->Hindex = -1;
      }
  }

  void Reset_Simulation_Data(simulation_data* data)
  {
    memset(&data->my_en,0,sizeof(energy_data));
  }

  void Reset_Workspace(reax_system *system, storage *workspace)
  {
    memset(workspace->total_bond_order, 0, system->total_cap * sizeof(double));
    memset(workspace->dDeltap_self, 0, system->total_cap * sizeof(rvec));
    memset(workspace->CdDelta, 0, system->total_cap * sizeof(double));
    memset(workspace->f, 0, system->total_cap * sizeof(rvec));

  }

  static void Reset_Neighbor_Lists(reax_system *system, control_params *control,
                             storage *workspace, reax_list **lists)
  {
    int i, total_bonds, Hindex, total_hbonds;
    reax_list *bonds, *hbonds;

    /* bonds list */
    if (system->N > 0) {
      bonds = (*lists) + BONDS;
      total_bonds = 0;

      /* reset start-end indexes */
      for (i = 0; i < system->N; ++i) {
        Set_Start_Index(i, total_bonds, bonds);
        Set_End_Index(i, total_bonds, bonds);
        total_bonds += system->my_atoms[i].num_bonds;
      }

      /* is reallocation needed? */
      if (total_bonds >= bonds->num_intrs * DANGER_ZONE) {
        workspace->realloc.bonds = 1;
        if (total_bonds >= bonds->num_intrs)
          control->error_ptr->one(FLERR,fmt::format("Not enough space for bonds! "
                                                    "total={} allocated={}\n",
                                                    total_bonds, bonds->num_intrs));
      }
    }

    if (control->hbond_cut > 0 && system->numH > 0) {
      hbonds = (*lists) + HBONDS;
      total_hbonds = 0;

      /* reset start-end indexes */
      for (i = 0; i < system->n; ++i) {
        Hindex = system->my_atoms[i].Hindex;
        if (Hindex > -1) {
          Set_Start_Index(Hindex, total_hbonds, hbonds);
          Set_End_Index(Hindex, total_hbonds, hbonds);
          total_hbonds += system->my_atoms[i].num_hbonds;
        }
      }

      /* is reallocation needed? */
      if (total_hbonds >= hbonds->num_intrs * 0.90/*DANGER_ZONE*/) {
        workspace->realloc.hbonds = 1;
        if (total_hbonds >= hbonds->num_intrs)
          control->error_ptr->one(FLERR,fmt::format("Not enough space for hbonds! "
                                                    "total={} allocated={}\n",
                                                    total_hbonds, hbonds->num_intrs));
      }
    }
  }


  void Reset(reax_system *system, control_params *control,
             simulation_data *data, storage *workspace, reax_list **lists)
  {
    Reset_Atoms(system, control);
    Reset_Simulation_Data(data);
    Reset_Workspace(system, workspace);
    Reset_Neighbor_Lists(system, control, workspace, lists);
  }
}
