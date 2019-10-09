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

#include <mpi.h>
#include "compute_cac_ke.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCACKE::ComputeCACKE(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute ke command");

  scalar_flag = 1;
  extscalar = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeCACKE::init()
{
  if (!atom->CAC_flag) error->all(FLERR,"compute cac/nodal/temp requires a CAC atom style");
  pfactor = 0.5 * force->mvv2e;
}

/* ---------------------------------------------------------------------- */

double ComputeCACKE::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double ****nodal_velocities = atom->nodal_velocities;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int **node_types = atom->node_types;
  int **element_scale = atom->element_scale;
  int *nodes_count_list = atom->nodes_per_element_list;	
  int nodes_per_element;

  double ke = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit){
         nodes_per_element = nodes_count_list[element_type[i]];
				  for (int ipoly = 0; ipoly < poly_count[i]; ipoly++)
					  for (int n = 0; n < nodes_per_element; n++)
						  ke += (nodal_velocities[i][n][ipoly][0] * nodal_velocities[i][n][ipoly][0]
							  + nodal_velocities[i][n][ipoly][1] * nodal_velocities[i][n][ipoly][1]
							  + nodal_velocities[i][n][ipoly][2] * nodal_velocities[i][n][ipoly][2])
						  * rmass[i]*element_scale[i][0]*element_scale[i][1]*element_scale[i][2]/nodes_per_element;
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit){
        nodes_per_element = nodes_count_list[element_type[i]];
				  for (int ipoly = 0; ipoly < poly_count[i]; ipoly++)
					  for (int n = 0; n < nodes_per_element; n++)
						  ke += (nodal_velocities[i][n][ipoly][0] * nodal_velocities[i][n][ipoly][0]
							  + nodal_velocities[i][n][ipoly][1] * nodal_velocities[i][n][ipoly][1]
							  + nodal_velocities[i][n][ipoly][2] * nodal_velocities[i][n][ipoly][2])*
						  mass[node_types[i][ipoly]]*element_scale[i][0]*element_scale[i][1]*element_scale[i][2]/nodes_per_element;
      }
    }
  }

  MPI_Allreduce(&ke,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar *= pfactor;
  return scalar;
}
