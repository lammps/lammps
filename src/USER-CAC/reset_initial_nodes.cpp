/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "reset_initial_nodes.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "comm.h"
#include "special.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ResetInitialNodes::ResetInitialNodes(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void ResetInitialNodes::command(int narg, char **/*arg*/)
{
  //check if simulation box has been defined
  if (domain->box_exist == 0)
    error->all(FLERR,"Reset_initial_nodes command before simulation box is defined");
  //check if atom and element data has been read in already
  if (atom->natoms == 0)
    error->all(FLERR,"Reset_initial_nodes command before material data has been defined");
  //check if CAC atom style is defined
  if(!atom->CAC_flag)
  error->all(FLERR, "CAC dump styles require a CAC atom style");

  //reset initial nodal positions with the current nodal positions
  int *type = atom->type;
  int *mask = atom->mask;
  double ****nodal_positions = atom->nodal_positions;
  double ****initial_nodal_positions = atom->initial_nodal_positions;
  int *nodes_per_element_list = atom->nodes_per_element_list;
  //double ****initial_nodal_positions = atom->initial_nodal_positions;
  int nlocal = atom->nlocal;
  int *poly_count = atom->poly_count;
  int *element_type = atom->element_type;
  int **node_types = atom->node_types;
  int **element_scale = atom->element_scale;

  for (int i = 0; i < nlocal; i++) {
    for (int k = 0; k < poly_count[i]; k++) {
      for (int j = 0; j < nodes_per_element_list[element_type[i]]; j++) {
        initial_nodal_positions[i][k][j][0] = nodal_positions[i][k][j][0];
        initial_nodal_positions[i][k][j][1] = nodal_positions[i][k][j][1];
        initial_nodal_positions[i][k][j][2] = nodal_positions[i][k][j][2];
      }
      }
  }
}
