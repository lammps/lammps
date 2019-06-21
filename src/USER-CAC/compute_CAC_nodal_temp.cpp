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

#include <cstring>
#include <mpi.h>
#include "compute_cac_nodal_temp.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "comm.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNodalTemp::ComputeNodalTemp(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute temp command");

  scalar_flag = 1;
  //size_vector = 6;
  extscalar = 0;
  //extvector = 1;
  tempflag = 1;

  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeNodalTemp::~ComputeNodalTemp()
{
  if (!copymode)
    delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeNodalTemp::setup()
{
  if (!atom->CAC_flag) error->all(FLERR,"compute cac/nodal/temp requires a CAC atom style");
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeNodalTemp::dof_compute()
{
  adjust_dof_fix();
  natoms_temp = group->count(igroup);
  dof = domain->dimension * natoms_temp;
  dof -= extra_dof + fix_dof;
  if (dof > 0.0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeNodalTemp::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **v = atom->v;
  double ****nodal_positions = atom->nodal_positions;
  double ****nodal_velocities = atom->nodal_velocities;
  double ****nodal_forces = atom->nodal_forces;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int **node_types = atom->node_types;
  int **element_scale = atom->element_scale;
  int *nodes_count_list = atom->nodes_per_element_list;	
  int nodes_per_element;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
 
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;

  if (rmass) {

	  for (int i = 0; i < nlocal; i++) {
		  
			  if (mask[i] & groupbit)
           nodes_per_element = nodes_count_list[element_type[i]];
				  for (int ipoly = 0; ipoly < poly_count[i]; ipoly++)
					  for (int n = 0; n < nodes_per_element; n++)
						  t += (nodal_velocities[i][n][ipoly][0] * nodal_velocities[i][n][ipoly][0]
							  + nodal_velocities[i][n][ipoly][1] * nodal_velocities[i][n][ipoly][1]
							  + nodal_velocities[i][n][ipoly][2] * nodal_velocities[i][n][ipoly][2])
						  * rmass[i] / nodes_per_element / poly_count[i];
		  
	  }
  }
  else {
	  
	  for (int i = 0; i < nlocal; i++) {

			  if (mask[i] & groupbit)
         nodes_per_element = nodes_count_list[element_type[i]];
				  for (int ipoly = 0; ipoly < poly_count[i]; ipoly++)
					  for (int n = 0; n < nodes_per_element; n++)
						  t += (nodal_velocities[i][n][ipoly][0] * nodal_velocities[i][n][ipoly][0]
							  + nodal_velocities[i][n][ipoly][1] * nodal_velocities[i][n][ipoly][1]
							  + nodal_velocities[i][n][ipoly][2] * nodal_velocities[i][n][ipoly][2])*
						  mass[node_types[i][ipoly]] / nodes_per_element / poly_count[i];
		  
	  }
  }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}


