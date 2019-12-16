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

#include <cstdio>
#include <cstring>
#include "fix_cac_nve.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVECAC::FixNVECAC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"nve/sphere") != 0 && narg < 3)
    error->all(FLERR,"Illegal fix cac/nve command");

  dynamic_group_allow = 1;
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixNVECAC::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVECAC::init()
{
  if (!atom->CAC_flag) error->all(FLERR,"CAC fix styles require a CAC atom style");
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  if (strstr(update->integrate_style,"respa"))
    step_respa = ((Respa *) update->integrate)->step;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVECAC::initial_integrate(int vflag)
{
  double dtfm;

  // update v and x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double ****nodal_positions=atom->nodal_positions;
  double ****nodal_velocities=atom->nodal_velocities;
  double ****nodal_forces=atom->nodal_forces;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int **node_types = atom->node_types;
  int *nodes_count_list = atom->nodes_per_element_list;	

  int nodes_per_element;

  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
  for (int i = 0; i < nlocal; i++){
  nodes_per_element = nodes_count_list[element_type[i]];

  if (mask[i] & groupbit) {
    x[i][0] = 0;
    x[i][1] = 0;
    x[i][2] = 0;
    v[i][0] = 0;
    v[i][1] = 0;
    v[i][2] = 0;
  for (int poly_counter = 0; poly_counter < poly_count[i];poly_counter++) {	
      for(int k=0; k<nodes_per_element; k++){	
        dtfm = dtf / rmass[i];
        nodal_velocities[i][poly_counter][k][0] += dtfm * nodal_forces[i][poly_counter][k][0];
        nodal_velocities[i][poly_counter][k][1] += dtfm * nodal_forces[i][poly_counter][k][1];
        nodal_velocities[i][poly_counter][k][2] += dtfm * nodal_forces[i][poly_counter][k][2];
        nodal_positions[i][poly_counter][k][0] += dtv * nodal_velocities[i][poly_counter][k][0];
        nodal_positions[i][poly_counter][k][1] += dtv * nodal_velocities[i][poly_counter][k][1];
        nodal_positions[i][poly_counter][k][2] += dtv * nodal_velocities[i][poly_counter][k][2];

        x[i][0] += nodal_positions[i][poly_counter][k][0];
        x[i][1] += nodal_positions[i][poly_counter][k][1];
        x[i][2] += nodal_positions[i][poly_counter][k][2];
        v[i][0] += nodal_velocities[i][poly_counter][k][0];
        v[i][1] += nodal_velocities[i][poly_counter][k][1];
        v[i][2] += nodal_velocities[i][poly_counter][k][2];

      }
    }
  x[i][0] = x[i][0] / nodes_per_element / poly_count[i];
  x[i][1] = x[i][1] / nodes_per_element / poly_count[i];
  x[i][2] = x[i][2] / nodes_per_element / poly_count[i];
  v[i][0] = v[i][0] / nodes_per_element / poly_count[i];
  v[i][1] = v[i][1] / nodes_per_element / poly_count[i];
  v[i][2] = v[i][2] / nodes_per_element / poly_count[i];
  }



  }

  } else {
    for (int i = 0; i < nlocal; i++){
    
    nodes_per_element = nodes_count_list[element_type[i]];

      if (mask[i] & groupbit) {
        x[i][0] = 0;
        x[i][1] = 0;
        x[i][2] = 0;
        v[i][0] = 0;
        v[i][1] = 0;
        v[i][2] = 0;
      for (int poly_counter = 0; poly_counter < poly_count[i];poly_counter++) {	
        for(int k=0; k<nodes_per_element; k++){	
          
          dtfm = dtf / mass[node_types[i][poly_counter]];
          nodal_velocities[i][poly_counter][k][0] += dtfm * nodal_forces[i][poly_counter][k][0];
          nodal_velocities[i][poly_counter][k][1] += dtfm * nodal_forces[i][poly_counter][k][1];
          nodal_velocities[i][poly_counter][k][2] += dtfm * nodal_forces[i][poly_counter][k][2];
          nodal_positions[i][poly_counter][k][0] += dtv * nodal_velocities[i][poly_counter][k][0];
          nodal_positions[i][poly_counter][k][1] += dtv * nodal_velocities[i][poly_counter][k][1];
          nodal_positions[i][poly_counter][k][2] += dtv * nodal_velocities[i][poly_counter][k][2];

          x[i][0] += nodal_positions[i][poly_counter][k][0];
          x[i][1] += nodal_positions[i][poly_counter][k][1];
          x[i][2] += nodal_positions[i][poly_counter][k][2];
          v[i][0] += nodal_velocities[i][poly_counter][k][0];
          v[i][1] += nodal_velocities[i][poly_counter][k][1];
          v[i][2] += nodal_velocities[i][poly_counter][k][2];

        }
      }
      x[i][0] = x[i][0] / nodes_per_element / poly_count[i];
      x[i][1] = x[i][1] / nodes_per_element / poly_count[i];
      x[i][2] = x[i][2] / nodes_per_element / poly_count[i];
      v[i][0] = v[i][0] / nodes_per_element / poly_count[i];
      v[i][1] = v[i][1] / nodes_per_element / poly_count[i];
      v[i][2] = v[i][2] / nodes_per_element / poly_count[i];
      }

    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVECAC::final_integrate()
{
  double dtfm;

  // update v of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double ****nodal_positions = atom->nodal_positions;
  double ****nodal_velocities = atom->nodal_velocities;
  double ****nodal_forces = atom->nodal_forces;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int **node_types = atom->node_types;
  int *nodes_count_list = atom->nodes_per_element_list;	

  int nodes_per_element;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++){
    nodes_per_element = nodes_count_list[element_type[i]];
         
      if (mask[i] & groupbit) {
        v[i][0] = 0;
        v[i][1] = 0;
        v[i][2] = 0;
      for (int poly_counter = 0; poly_counter < poly_count[i];poly_counter++) {	
        for(int k=0; k<nodes_per_element; k++){	
          
            dtfm = dtf / rmass[i];
            nodal_velocities[i][poly_counter][k][0] += dtfm * nodal_forces[i][poly_counter][k][0];
            nodal_velocities[i][poly_counter][k][1] += dtfm * nodal_forces[i][poly_counter][k][1];
            nodal_velocities[i][poly_counter][k][2] += dtfm * nodal_forces[i][poly_counter][k][2];

            v[i][0] += nodal_velocities[i][poly_counter][k][0];
            v[i][1] += nodal_velocities[i][poly_counter][k][1];
            v[i][2] += nodal_velocities[i][poly_counter][k][2];

          }
        }
      v[i][0] = v[i][0] / nodes_per_element / poly_count[i];
      v[i][1] = v[i][1] / nodes_per_element / poly_count[i];
      v[i][2] = v[i][2] / nodes_per_element / poly_count[i];
      }

    }

  } else {
    for (int i = 0; i < nlocal; i++){
    nodes_per_element = nodes_count_list[element_type[i]];

      if (mask[i] & groupbit) {
        v[i][0] = 0;
        v[i][1] = 0;
        v[i][2] = 0;
      for (int poly_counter = 0; poly_counter < poly_count[i];poly_counter++) {	
        for(int k=0; k<nodes_per_element; k++){	
          
            dtfm = dtf / mass[node_types[i][poly_counter]];
            nodal_velocities[i][poly_counter][k][0] += dtfm * nodal_forces[i][poly_counter][k][0];
            nodal_velocities[i][poly_counter][k][1] += dtfm * nodal_forces[i][poly_counter][k][1];
            nodal_velocities[i][poly_counter][k][2] += dtfm * nodal_forces[i][poly_counter][k][2];

            v[i][0] += nodal_velocities[i][poly_counter][k][0];
            v[i][1] += nodal_velocities[i][poly_counter][k][1];
            v[i][2] += nodal_velocities[i][poly_counter][k][2];

          }
        }
      v[i][0] = v[i][0] / nodes_per_element / poly_count[i];
      v[i][1] = v[i][1] / nodes_per_element / poly_count[i];
      v[i][2] = v[i][2] / nodes_per_element / poly_count[i];
      }

    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVECAC::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  // innermost level - NVE update of v and x
  // all other levels - NVE update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVECAC::final_integrate_respa(int ilevel, int iloop)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVECAC::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
