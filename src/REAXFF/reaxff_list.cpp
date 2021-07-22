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

namespace ReaxFF {

  /************* allocate list space ******************/
  void Make_List(int n, int num_intrs, int type, reax_list *l)
  {
    l->allocated = 1;

    l->n = n;
    l->num_intrs = num_intrs;

    if (l->index) sfree(l->error_ptr, l->index, "list:index");
    if (l->end_index) sfree(l->error_ptr, l->end_index, "list:end_index");
    l->index = (int*) smalloc(l->error_ptr, n * sizeof(int), "list:index");
    l->end_index = (int*) smalloc(l->error_ptr, n * sizeof(int), "list:end_index");

    l->type = type;

    switch(l->type) {
    case TYP_THREE_BODY:
      if (l->select.three_body_list) sfree(l->error_ptr, l->select.three_body_list,"list:three_bodies");
      l->select.three_body_list = (three_body_interaction_data*)
        smalloc(l->error_ptr, (rc_bigint) num_intrs * sizeof(three_body_interaction_data),
                "list:three_bodies");
      break;

    case TYP_BOND:
      if (l->select.bond_list) sfree(l->error_ptr, l->select.bond_list,"list:bonds");
      l->select.bond_list = (bond_data*)
        smalloc(l->error_ptr, (rc_bigint) num_intrs * sizeof(bond_data), "list:bonds");
      break;

    case TYP_FAR_NEIGHBOR:
      if (l->select.far_nbr_list) sfree(l->error_ptr, l->select.far_nbr_list,"list:far_nbrs");
      l->select.far_nbr_list = (far_neighbor_data*)
        smalloc(l->error_ptr, (rc_bigint) num_intrs * sizeof(far_neighbor_data), "list:far_nbrs");
      break;

    case TYP_HBOND:
      if (l->select.hbond_list) sfree(l->error_ptr, l->select.hbond_list,"list:hbonds");
      l->select.hbond_list = (hbond_data*)
        smalloc(l->error_ptr, (rc_bigint) num_intrs * sizeof(hbond_data), "list:hbonds");
      break;

    default:
      l->error_ptr->all(FLERR,fmt::format("No list type {} defined", l->type));
    }
  }

  void Delete_List(reax_list *l)
  {
    if (l->allocated == 0)
      return;
    l->allocated = 0;

    sfree(l->error_ptr, l->index, "list:index");
    sfree(l->error_ptr, l->end_index, "list:end_index");
    l->index = nullptr;
    l->end_index = nullptr;

    switch(l->type) {
    case TYP_HBOND:
      sfree(l->error_ptr, l->select.hbond_list, "list:hbonds");
      l->select.hbond_list = nullptr;
      break;
    case TYP_FAR_NEIGHBOR:
      sfree(l->error_ptr, l->select.far_nbr_list, "list:far_nbrs");
      l->select.far_nbr_list = nullptr;
      break;
    case TYP_BOND:
      sfree(l->error_ptr, l->select.bond_list, "list:bonds");
      l->select.bond_list = nullptr;
      break;
    case TYP_THREE_BODY:
      sfree(l->error_ptr, l->select.three_body_list, "list:three_bodies");
      l->select.three_body_list = nullptr;
      break;

    default:
      l->error_ptr->all(FLERR,fmt::format("No list type {} defined", l->type));
    }
  }
}
