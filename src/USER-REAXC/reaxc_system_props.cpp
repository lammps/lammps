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

#include "reaxc_system_props.h"
#include <mpi.h>
#include "reaxc_defs.h"

void Compute_System_Energy( reax_system *system, simulation_data *data,
                            MPI_Comm comm )
{
  double my_en[15], sys_en[15];

  my_en[0] = data->my_en.e_bond;
  my_en[1] = data->my_en.e_ov;
  my_en[2] = data->my_en.e_un;
  my_en[3] = data->my_en.e_lp;
  my_en[4] = data->my_en.e_ang;
  my_en[5] = data->my_en.e_pen;
  my_en[6] = data->my_en.e_coa;
  my_en[7] = data->my_en.e_hb;
  my_en[8] = data->my_en.e_tor;
  my_en[9] = data->my_en.e_con;
  my_en[10] = data->my_en.e_vdW;
  my_en[11] = data->my_en.e_ele;
  my_en[12] = data->my_en.e_pol;
  my_en[13] = data->my_en.e_kin;
  MPI_Reduce( my_en, sys_en, 14, MPI_DOUBLE, MPI_SUM, MASTER_NODE, comm );

  data->my_en.e_pot = data->my_en.e_bond +
    data->my_en.e_ov + data->my_en.e_un  + data->my_en.e_lp +
    data->my_en.e_ang + data->my_en.e_pen + data->my_en.e_coa +
    data->my_en.e_hb +
    data->my_en.e_tor + data->my_en.e_con +
    data->my_en.e_vdW + data->my_en.e_ele + data->my_en.e_pol;

  data->my_en.e_tot = data->my_en.e_pot + E_CONV * data->my_en.e_kin;

  if (system->my_rank == MASTER_NODE) {
    data->sys_en.e_bond = sys_en[0];
    data->sys_en.e_ov = sys_en[1];
    data->sys_en.e_un = sys_en[2];
    data->sys_en.e_lp = sys_en[3];
    data->sys_en.e_ang = sys_en[4];
    data->sys_en.e_pen = sys_en[5];
    data->sys_en.e_coa = sys_en[6];
    data->sys_en.e_hb = sys_en[7];
    data->sys_en.e_tor = sys_en[8];
    data->sys_en.e_con = sys_en[9];
    data->sys_en.e_vdW = sys_en[10];
    data->sys_en.e_ele = sys_en[11];
    data->sys_en.e_pol = sys_en[12];
    data->sys_en.e_kin = sys_en[13];

    data->sys_en.e_pot = data->sys_en.e_bond +
      data->sys_en.e_ov + data->sys_en.e_un  + data->sys_en.e_lp +
      data->sys_en.e_ang + data->sys_en.e_pen + data->sys_en.e_coa +
      data->sys_en.e_hb +
      data->sys_en.e_tor + data->sys_en.e_con +
      data->sys_en.e_vdW + data->sys_en.e_ele + data->sys_en.e_pol;

    data->sys_en.e_tot = data->sys_en.e_pot + E_CONV * data->sys_en.e_kin;
  }
}
