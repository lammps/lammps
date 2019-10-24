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
#include "compute_cac_ke_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCACKEAtom::ComputeCACKEAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  ke(NULL)
{
  if (narg != 3) error->all(FLERR,"Illegal compute ke/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeCACKEAtom::~ComputeCACKEAtom()
{
  memory->destroy(ke);
}

/* ---------------------------------------------------------------------- */

void ComputeCACKEAtom::init()
{
  int count = 0;
  if (!atom->CAC_flag) error->all(FLERR,"compute cac/ke/atom requires a CAC atom style");
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"cac/ke/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute cac/ke/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeCACKEAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow ke array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(ke);
    nmax = atom->nmax;
    memory->create(ke,nmax,"ke/atom:ke");
    vector_atom = ke;
  }

  // compute kinetic energy for each atom in group

  double mvv2e = force->mvv2e;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
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

  if (rmass) {
    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit){
         nodes_per_element = nodes_count_list[element_type[i]];
         ke[i]=0;
				  for (int ipoly = 0; ipoly < poly_count[i]; ipoly++)
					  for (int n = 0; n < nodes_per_element; n++)
						  ke[i] = 0.5 * mvv2e *(nodal_velocities[i][ipoly][n][0] * nodal_velocities[i][ipoly][n][0]
							  + nodal_velocities[i][ipoly][n][1] * nodal_velocities[i][ipoly][n][1]
							  + nodal_velocities[i][ipoly][n][2] * nodal_velocities[i][ipoly][n][2])
						  * rmass[i]*element_scale[i][0]*element_scale[i][1]*element_scale[i][2]/nodes_per_element;
      }
      else ke[i] = 0;
    }
  } else {
    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit){
        nodes_per_element = nodes_count_list[element_type[i]];
        ke[i]=0;
				  for (int ipoly = 0; ipoly < poly_count[i]; ipoly++)
					  for (int n = 0; n < nodes_per_element; n++)
						  ke[i] += 0.5 * mvv2e *(nodal_velocities[i][ipoly][n][0] * nodal_velocities[i][ipoly][n][0]
							  + nodal_velocities[i][ipoly][n][1] * nodal_velocities[i][ipoly][n][1]
							  + nodal_velocities[i][ipoly][n][2] * nodal_velocities[i][ipoly][n][2])*
						  mass[node_types[i][ipoly]]*element_scale[i][0]*element_scale[i][1]*element_scale[i][2]/nodes_per_element;
      }
      else ke[i] = 0;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeCACKEAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
