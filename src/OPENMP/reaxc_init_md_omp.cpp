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

#include "reaxff_omp.h"
#include "reaxff_api.h"

#include "error.h"

#include <mpi.h>
#include <cstdlib>

namespace ReaxFF {
  static void Init_ListsOMP(reax_system *system, control_params *control,
                            reax_list **lists)
  {
    int i, total_hbonds, total_bonds, bond_cap, num_3body, cap_3body, Htop;
    int *hb_top, *bond_top;

    int mincap = system->mincap;
    double safezone = system->safezone;
    double saferzone = system->saferzone;
    auto error = system->error_ptr;

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

    int nthreads = control->nthreads;
    reax_list *bonds = (*lists)+BONDS;

    for (i = 0; i < bonds->num_intrs; ++i)
      bonds->select.bond_list[i].bo_data.CdboReduction =
        (double*) smalloc(error, sizeof(double)*nthreads, "CdboReduction");

    /* 3bodies list */
    cap_3body = (int)(MAX(num_3body*safezone, MIN_3BODIES));
    Make_List(bond_cap, cap_3body, TYP_THREE_BODY,*lists+THREE_BODIES);
    (*lists+THREE_BODIES)->error_ptr = system->error_ptr;

    free(hb_top);
    free(bond_top);
  }

  void InitializeOMP(reax_system *system, control_params *control,
                     simulation_data *data, storage *workspace,
                     reax_list **lists, MPI_Comm world)
  {
    Init_System(system,control);
    Init_Simulation_Data(data);
    Init_Workspace(system,control,workspace);
    Init_ListsOMP(system,control,lists);
    if (control->tabulate)
      Init_Lookup_Tables(system,control,workspace,world);
  }
}
