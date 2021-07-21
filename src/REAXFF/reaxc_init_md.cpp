// clang-format off
/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University

  Contributing authors:
  H. M. Aktulga, J. Fogarty, S. Pandit, A. Grama
  Corresponding author:
  Hasan Metin Aktulga, Michigan State University, hma@cse.msu.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, 38 (4-5), 245-259

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

#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>

namespace ReaxFF {

  void Init_System(reax_system *system, control_params *control)
  {
    int i;
    reax_atom *atom;

    int mincap = system->mincap;
    double safezone = system->safezone;
    double saferzone = system->saferzone;

    // determine the local and total capacity

    system->local_cap = MAX((int)(system->n * safezone), mincap);
    system->total_cap = MAX((int)(system->N * safezone), mincap);

    /* estimate numH and Hcap */
    system->numH = 0;
    if (control->hbond_cut > 0)
      for (i = 0; i < system->n; ++i) {
        atom = &(system->my_atoms[i]);
        if (system->reax_param.sbp[atom->type].p_hbond == 1 && atom->type >= 0)
          atom->Hindex = system->numH++;
        else atom->Hindex = -1;
      }
    system->Hcap = (int)(MAX(system->numH * saferzone, mincap));
  }

  void Init_Simulation_Data(simulation_data *data)
  {
    Reset_Simulation_Data(data);
    data->step = 0;
  }

  static void Init_Taper(control_params *control, storage *workspace)
  {
    double d1, d7;
    double swa, swa2, swa3;
    double swb, swb2, swb3;
    LAMMPS_NS::Error *error = control->error_ptr;

    swa = control->nonb_low;
    swb = control->nonb_cut;

    if (fabs(swa) > 0.01 && control->me == 0)
      error->warning(FLERR, "Non-zero lower Taper-radius cutoff");

    if (swb < 0) {
      error->all(FLERR,"Negative upper Taper-radius cutoff");
    }
    else if (swb < 5 && control->me == 0)
      error->warning(FLERR,fmt::format("Warning: very low Taper-radius cutoff: "
                                       "{}\n", swb));
    d1 = swb - swa;
    d7 = pow(d1, 7.0);
    swa2 = SQR(swa);
    swa3 = CUBE(swa);
    swb2 = SQR(swb);
    swb3 = CUBE(swb);

    workspace->Tap[7] =  20.0 / d7;
    workspace->Tap[6] = -70.0 * (swa + swb) / d7;
    workspace->Tap[5] =  84.0 * (swa2 + 3.0*swa*swb + swb2) / d7;
    workspace->Tap[4] = -35.0 * (swa3 + 9.0*swa2*swb + 9.0*swa*swb2 + swb3) / d7;
    workspace->Tap[3] = 140.0 * (swa3*swb + 3.0*swa2*swb2 + swa*swb3) / d7;
    workspace->Tap[2] =-210.0 * (swa3*swb2 + swa2*swb3) / d7;
    workspace->Tap[1] = 140.0 * swa3 * swb3 / d7;
    workspace->Tap[0] = (-35.0*swa3*swb2*swb2 + 21.0*swa2*swb3*swb2 -
                         7.0*swa*swb3*swb3 + swb3*swb3*swb) / d7;
  }

  void Init_Workspace(reax_system *system, control_params *control, storage *workspace)
  {
    Allocate_Workspace(control, workspace,system->total_cap);

    memset(&(workspace->realloc), 0, sizeof(reallocate_data));
    Reset_Workspace(system, workspace);

    /* Initialize the Taper function */
    Init_Taper(control, workspace);
  }

  static void Init_Lists(reax_system *system, control_params *control, reax_list **lists)
  {
    int i, total_hbonds, total_bonds, bond_cap, num_3body, cap_3body, Htop;
    int *hb_top, *bond_top;

    int mincap = system->mincap;
    double safezone = system->safezone;
    double saferzone = system->saferzone;

    bond_top = (int*) calloc(system->total_cap, sizeof(int));
    hb_top = (int*) calloc(system->local_cap, sizeof(int));
    Estimate_Storages(system, control, lists,
                      &Htop, hb_top, bond_top, &num_3body);

    if (control->hbond_cut > 0) {
      /* init H indexes */
      total_hbonds = 0;
      for (i = 0; i < system->n; ++i) {
        system->my_atoms[i].num_hbonds = hb_top[i];
        total_hbonds += hb_top[i];
      }
      total_hbonds = (int)(MAX(total_hbonds*saferzone,mincap*system->minhbonds));

      Make_List(system->Hcap, total_hbonds, TYP_HBOND,*lists+HBONDS);
      (*lists+HBONDS)->error_ptr = system->error_ptr;
    }

    total_bonds = 0;
    for (i = 0; i < system->N; ++i) {
      system->my_atoms[i].num_bonds = bond_top[i];
      total_bonds += bond_top[i];
    }
    bond_cap = (int)(MAX(total_bonds*safezone, mincap*MIN_BONDS));

    Make_List(system->total_cap, bond_cap, TYP_BOND,*lists+BONDS);
    (*lists+BONDS)->error_ptr = system->error_ptr;

    /* 3bodies list */
    cap_3body = (int)(MAX(num_3body*safezone, MIN_3BODIES));
    Make_List(bond_cap, cap_3body, TYP_THREE_BODY,*lists+THREE_BODIES);
    (*lists+THREE_BODIES)->error_ptr = system->error_ptr;

    free(hb_top);
    free(bond_top);
  }

  void Initialize(reax_system *system, control_params *control,
                  simulation_data *data, storage *workspace,
                  reax_list **lists, MPI_Comm world)
  {
    Init_System(system,control);
    Init_Simulation_Data(data);
    Init_Workspace(system,control,workspace);
    Init_Lists(system,control,lists);
    if (control->tabulate)
      Init_Lookup_Tables(system,control,workspace,world);
  }
}
