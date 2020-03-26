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

#include <cstdlib>
#include <cstring>
#include "atom_vec_cac_charge.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "asa_data.h"

#define MAX_ELEMENT_NAME 256

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecCAC_Charge::AtomVecCAC_Charge(LAMMPS *lmp) : AtomVecCAC(lmp)
{
  molecular = 0;
  mass_type = 1;
  element_type_count = 0;
  scale_count=0;
  comm_x_only = comm_f_only = 0;
  size_forward = 3;
  size_reverse = 3;
  size_border = 6;
  size_velocity = 3;
  size_data_atom = 5;
  size_data_vel = 4;
  xcol_data = 4;
  forceclearflag = 1;
  atom->CAC_flag=1;
  atom->q_flag = 1;
  search_range_max = 0;
  initial_size=0;
  check_distance_flag=1;
}

//--------------------------------------------------------------------------

AtomVecCAC_Charge::~AtomVecCAC_Charge() {
for(int element_index=0; element_index < alloc_counter; element_index++){
  memory->destroy(node_charges[element_index]);
}
}

/* ----------------------------------------------------------------------
   process user input
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecCAC_Charge::process_args(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Invalid atom_style cac command");
  nodes_per_element=force->inumeric(FLERR,arg[0]);
  maxpoly = force->inumeric(FLERR, arg[1]);
  atom->nodes_per_element=nodes_per_element;
  atom-> words_per_node = 7;
  atom->maxpoly = maxpoly;

  size_forward = 9*nodes_per_element*maxpoly +8+ 2*maxpoly;
  size_reverse = 3; // 3 + drho + de
  size_border = 9*nodes_per_element*maxpoly +11+ 2 * maxpoly;
  size_velocity = 9*nodes_per_element*maxpoly +11+ 2 * maxpoly;
  size_data_atom = 3*nodes_per_element*maxpoly +10+ 2 * maxpoly;
  size_data_vel = 9*nodes_per_element*maxpoly +9+ 2 * maxpoly;
  xcol_data = 4;

  maxexchange=size_border;

  //initialize node counts associated with each element type
  //call setup for element types
  define_elements();

  //create array that tests in data_atom for odd node to iDod counts
  memory->create(node_count_per_poly, maxpoly, "AtomVecCAC: node_count_per_poly");
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecCAC_Charge::grow(int n)
{
  if (n == 0) grow_nmax();
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
  error->one(FLERR, "Per-processor system is too big");

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f, nmax*comm->nthreads, 3, "atom:f");
  poly_count = memory->grow(atom->poly_count, nmax, "atom:type_count");
  element_type= memory->grow(atom->element_type, nmax, "atom:element_type");
  element_scale = memory->grow(atom->element_scale, nmax,3, "atom:element_scales");

  //grow pointers for a ragged allocation strategy since atoms allocate far less memory
  if(CAC_nmax==0){
  atom->node_types = node_types = (int **) memory->smalloc(sizeof(int *)*nmax, "atom:node_types");
  atom->node_charges = node_charges = (double **) memory->smalloc(sizeof(double *)*nmax, "atom:node_charges");
  atom->nodal_positions = nodal_positions =
    (double ****) memory->smalloc(sizeof(double ***)*nmax, "atom:nodal_positions");
  hold_nodal_positions = (double ****) memory->smalloc(sizeof(double ***)*nmax, "atom:hold_nodal_positions");
  atom->initial_nodal_positions = initial_nodal_positions =
    (double ****) memory->smalloc(sizeof(double ***)*nmax, "atom:initial_nodal_positions");
  atom->nodal_velocities = nodal_velocities =
    (double ****) memory->smalloc(sizeof(double ***)*nmax, "atom:nodal_velocities");
  atom->nodal_forces = nodal_forces =
    (double ****) memory->smalloc(sizeof(double ***)*nmax, "atom:nodal_forces");
  atom->nodal_virial = nodal_virial =
    (double ****) memory->smalloc(sizeof(double ***)*nmax, "atom:nodal_virial");
  CAC_nmax = nmax;
  }
  else{
  atom->node_types = node_types = (int **) memory->srealloc(node_types,sizeof(int *)*nmax, "atom:node_types");
  atom->node_charges = node_charges = (double **) memory->srealloc(node_charges,sizeof(double *)*nmax, "atom:node_charges");
  atom->nodal_positions = nodal_positions =
    (double ****) memory->srealloc(nodal_positions,sizeof(double ***)*nmax, "atom:nodal_positions");
  hold_nodal_positions =
  (double ****) memory->srealloc(hold_nodal_positions,sizeof(double ***)*nmax, "atom:hold_nodal_positions");
  atom->initial_nodal_positions = initial_nodal_positions =
    (double ****) memory->srealloc(initial_nodal_positions,sizeof(double ***)*nmax, "atom:initial_nodal_positions");
  atom->nodal_velocities = nodal_velocities =
    (double ****) memory->srealloc(nodal_velocities,sizeof(double ***)*nmax, "atom:nodal_velocities");
  atom->nodal_forces = nodal_forces =
    (double ****) memory->srealloc(nodal_forces,sizeof(double ***)*nmax, "atom:nodal_forces");
  atom->nodal_virial = nodal_virial =
    (double ****) memory->srealloc(nodal_virial,sizeof(double ***)*nmax, "atom:nodal_virial");
  CAC_nmax = nmax;
  }
  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   resize atom arrays
------------------------------------------------------------------------- */

void AtomVecCAC_Charge::shrink_array(int n)
{
  if(n>nmax)
  error->one(FLERR, "resize function is called to shrink atom arrays; use grow instead");
  atom->nmax=nmax=n;
  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f, nmax*comm->nthreads, 3, "atom:f");
  poly_count = memory->grow(atom->poly_count, nmax, "atom:type_count");
  element_type= memory->grow(atom->element_type, nmax, "atom:element_type");
  element_scale = memory->grow(atom->element_scale, nmax,3, "atom:element_scales");

  //deallocate element contents if n is smaller than the alloc counter for elements
  for(int element_index=alloc_counter-1; element_index >= n; element_index--){
  memory->destroy(node_types[element_index]);
  memory->destroy(node_charges[element_index]);
  memory->destroy(nodal_positions[element_index]);
  memory->destroy(hold_nodal_positions[element_index]);
  memory->destroy(initial_nodal_positions[element_index]);
  memory->destroy(nodal_velocities[element_index]);
  memory->destroy(nodal_forces[element_index]);
  memory->destroy(nodal_virial[element_index]);
  }
  if(alloc_counter>n)
  alloc_counter = n;

  //shrink pointer arrays
  atom->node_types = node_types = (int **) memory->srealloc(node_types,sizeof(int *)*nmax, "atom:node_types");
  atom->node_charges = node_charges = (double **) memory->srealloc(node_charges,sizeof(double *)*nmax, "atom:node_charges");
  atom->nodal_positions = nodal_positions =
    (double ****) memory->srealloc(nodal_positions,sizeof(double ***)*nmax, "atom:nodal_positions");
  hold_nodal_positions =
  (double ****) memory->srealloc(hold_nodal_positions,sizeof(double ***)*nmax, "atom:hold_nodal_positions");
  atom->initial_nodal_positions = initial_nodal_positions =
    (double ****) memory->srealloc(initial_nodal_positions,sizeof(double ***)*nmax, "atom:initial_nodal_positions");
  atom->nodal_velocities = nodal_velocities =
    (double ****) memory->srealloc(nodal_velocities,sizeof(double ***)*nmax, "atom:nodal_velocities");
  atom->nodal_forces = nodal_forces =
    (double ****) memory->srealloc(nodal_forces,sizeof(double ***)*nmax, "atom:nodal_forces");
  atom->nodal_virial = nodal_virial =
    (double ****) memory->srealloc(nodal_virial,sizeof(double ***)*nmax, "atom:nodal_virial");
  CAC_nmax = n;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   allocate finite element data for nodal arrays
------------------------------------------------------------------------- */

void AtomVecCAC_Charge::allocate_element(int element_index, int node_count, int poly_count)
{
  //destroy previous contents at that element index if present
  if(element_index<alloc_counter){
  memory->destroy(node_types[element_index]);
  memory->destroy(node_charges[element_index]);
  memory->destroy(nodal_positions[element_index]);
  memory->destroy(hold_nodal_positions[element_index]);
  memory->destroy(initial_nodal_positions[element_index]);
  memory->destroy(nodal_velocities[element_index]);
  memory->destroy(nodal_forces[element_index]);
  memory->destroy(nodal_virial[element_index]);
  }
  else
  alloc_counter++;
  //create new allocation for this element
  memory->create(node_types[element_index], poly_count,  "atom:node_types");
  memory->create(node_charges[element_index], poly_count,  "atom:node_charges");
  memory->create(nodal_positions[element_index], poly_count, node_count, 3, "atom:nodal_positions");
  memory->create(hold_nodal_positions[element_index], poly_count, node_count, 3, "atom:hold_nodal_positions");
  memory->create(initial_nodal_positions[element_index], poly_count, node_count, 3, "atom:initial_nodal_positions");
  memory->create(nodal_velocities[element_index], poly_count, node_count, 3, "atom:nodal_velocities");
  memory->create(nodal_forces[element_index], poly_count, node_count, 3, "atom:nodal_forces");
  memory->create(nodal_virial[element_index], poly_count, node_count, 6, "atom:nodal_virial");

}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecCAC_Charge::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
    nodal_positions = atom->nodal_positions;
  initial_nodal_positions = atom->initial_nodal_positions;
  nodal_velocities = atom->nodal_velocities;
  nodal_forces = atom->nodal_forces;
  nodal_virial = atom->nodal_virial;
  poly_count = atom->poly_count;
  element_type = atom->element_type;
  element_scale = atom->element_scale;
  node_types = atom->node_types;
  node_charges = atom->node_charges;
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecCAC_Charge::copy(int i, int j, int delflag)
{
  int *nodes_count_list = atom->nodes_per_element_list;
  int node_count;
  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];
  element_type[j] = element_type[i];
  element_scale[j][0] = element_scale[i][0];
  element_scale[j][1] = element_scale[i][1];
  element_scale[j][2] = element_scale[i][2];
  poly_count[j] = poly_count[i];
  node_count = nodes_count_list[element_type[j]];
  //copy nodal information; requires resizing since copy might be for
  //an element of a different size
    //destroy previous contents at that element index if present
  if(j<alloc_counter){
  memory->destroy(node_types[j]);
  memory->destroy(node_charges[j]);
  memory->destroy(nodal_positions[j]);
  memory->destroy(hold_nodal_positions[j]);
  memory->destroy(initial_nodal_positions[j]);
  memory->destroy(nodal_velocities[j]);
  memory->destroy(nodal_forces[j]);
  memory->destroy(nodal_virial[j]);
  }
  else
  alloc_counter++;

  //create new allocation for this element
  memory->create(node_types[j], poly_count[j],  "atom:node_types");
  memory->create(node_charges[j], poly_count[j],  "atom:node_charges");
  memory->create(nodal_positions[j], poly_count[j], node_count, 3, "atom:nodal_positions");
  memory->create(hold_nodal_positions[j], poly_count[j], node_count, 3, "atom:hold_nodal_positions");
  memory->create(initial_nodal_positions[j], poly_count[j], node_count, 3, "atom:initial_nodal_positions");
  memory->create(nodal_velocities[j], poly_count[j], node_count, 3, "atom:nodal_velocities");
  memory->create(nodal_forces[j], poly_count[j], node_count, 3, "atom:nodal_forces");
  memory->create(nodal_virial[j], poly_count[j], node_count, 6, "atom:nodal_virial");

  for (int type_map = 0; type_map < poly_count[j]; type_map++) {
    node_types[j][type_map] = node_types[i][type_map];
    node_charges[j][type_map] = node_charges[i][type_map];
  }

   for (int poly_index = 0; poly_index < poly_count[j]; poly_index++){
    for(int nodecount=0; nodecount< nodes_count_list[element_type[j]]; nodecount++ )
      {
        nodal_positions[j][poly_index][nodecount][0] = nodal_positions[i][poly_index][nodecount][0];
        nodal_positions[j][poly_index][nodecount][1] = nodal_positions[i][poly_index][nodecount][1];
        nodal_positions[j][poly_index][nodecount][2] = nodal_positions[i][poly_index][nodecount][2];
        initial_nodal_positions[j][poly_index][nodecount][0] = initial_nodal_positions[i][poly_index][nodecount][0];
        initial_nodal_positions[j][poly_index][nodecount][1] = initial_nodal_positions[i][poly_index][nodecount][1];
        initial_nodal_positions[j][poly_index][nodecount][2] = initial_nodal_positions[i][poly_index][nodecount][2];
        nodal_velocities[j][poly_index][nodecount][0] = nodal_velocities[i][poly_index][nodecount][0];
        nodal_velocities[j][poly_index][nodecount][1] = nodal_velocities[i][poly_index][nodecount][1];
        nodal_velocities[j][poly_index][nodecount][2] = nodal_velocities[i][poly_index][nodecount][2];
      }
    }

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVecCAC_Charge::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;
  int *nodes_count_list = atom->nodes_per_element_list;
  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];

    buf[m++] = ubuf(element_type[j]).d;

    buf[m++] = ubuf(element_scale[j][0]).d;
    buf[m++] = ubuf(element_scale[j][1]).d;
    buf[m++] = ubuf(element_scale[j][2]).d;

    buf[m++] = ubuf(poly_count[j]).d;
    for (int type_map = 0; type_map < poly_count[j]; type_map++) {
      buf[m++] = ubuf(node_types[j][type_map]).d;
      buf[m++] = node_charges[j][type_map];
    }

    for (int poly_index = 0; poly_index < poly_count[j]; poly_index++){
      for (int nodecount = 0; nodecount< nodes_count_list[element_type[j]]; nodecount++)
        {
          buf[m++] = nodal_positions[j][poly_index][nodecount][0];
          buf[m++] = nodal_positions[j][poly_index][nodecount][1];
          buf[m++] = nodal_positions[j][poly_index][nodecount][2];
          buf[m++] = initial_nodal_positions[j][poly_index][nodecount][0];
          buf[m++] = initial_nodal_positions[j][poly_index][nodecount][1];
          buf[m++] = initial_nodal_positions[j][poly_index][nodecount][2];
          buf[m++] = nodal_velocities[j][poly_index][nodecount][0];
          buf[m++] = nodal_velocities[j][poly_index][nodecount][1];
          buf[m++] = nodal_velocities[j][poly_index][nodecount][2];
        }
      }
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;

    buf[m++] = ubuf(element_type[j]).d;

    buf[m++] = ubuf(element_scale[j][0]).d;
    buf[m++] = ubuf(element_scale[j][1]).d;
    buf[m++] = ubuf(element_scale[j][2]).d;

    buf[m++] = ubuf(poly_count[j]).d;
    for (int type_map = 0; type_map < poly_count[j]; type_map++) {
      buf[m++] = ubuf(node_types[j][type_map]).d;
      buf[m++] = node_charges[j][type_map];
    }

    for (int poly_index = 0; poly_index < poly_count[j]; poly_index++){
      for (int nodecount = 0; nodecount< nodes_count_list[element_type[j]]; nodecount++)
      {
        buf[m++] = nodal_positions[j][poly_index][nodecount][0]+dx;
        buf[m++] = nodal_positions[j][poly_index][nodecount][1] + dy;
        buf[m++] = nodal_positions[j][poly_index][nodecount][2] + dz;
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][0] + dx;
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][1] + dy;
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][2] + dz;
        buf[m++] = nodal_velocities[j][poly_index][nodecount][0];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][1];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][2];
      }
    }

    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecCAC_Charge::pack_comm_vel(int n, int *list, double *buf,
                                 int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;
  int *nodes_count_list = atom->nodes_per_element_list;
  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = ubuf(element_type[j]).d;

      buf[m++] = ubuf(element_scale[j][0]).d;
      buf[m++] = ubuf(element_scale[j][1]).d;
      buf[m++] = ubuf(element_scale[j][2]).d;

      buf[m++] = ubuf(poly_count[j]).d;
      for (int type_map = 0; type_map < poly_count[j]; type_map++) {
        buf[m++] = ubuf(node_types[j][type_map]).d;
        buf[m++] = node_charges[j][type_map];
      }

      for (int poly_index = 0; poly_index < poly_count[j]; poly_index++){
        for (int nodecount = 0; nodecount < nodes_count_list[element_type[j]]; nodecount++)
        {
          buf[m++] = nodal_positions[j][poly_index][nodecount][0];
          buf[m++] = nodal_positions[j][poly_index][nodecount][1];
          buf[m++] = nodal_positions[j][poly_index][nodecount][2];
          buf[m++] = initial_nodal_positions[j][poly_index][nodecount][0];
          buf[m++] = initial_nodal_positions[j][poly_index][nodecount][1];
          buf[m++] = initial_nodal_positions[j][poly_index][nodecount][2];
          buf[m++] = nodal_velocities[j][poly_index][nodecount][0];
          buf[m++] = nodal_velocities[j][poly_index][nodecount][1];
          buf[m++] = nodal_velocities[j][poly_index][nodecount][2];
        }
      }
    }
  }
  else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
    buf[m++] = ubuf(element_type[j]).d;

    buf[m++] = ubuf(element_scale[j][0]).d;
    buf[m++] = ubuf(element_scale[j][1]).d;
    buf[m++] = ubuf(element_scale[j][2]).d;

    buf[m++] = ubuf(poly_count[j]).d;
    for (int type_map = 0; type_map < poly_count[j]; type_map++) {
      buf[m++] = ubuf(node_types[j][type_map]).d;
      buf[m++] = node_charges[j][type_map];
    }

    for (int poly_index = 0; poly_index < poly_count[j]; poly_index++){
      for (int nodecount = 0; nodecount< nodes_count_list[element_type[j]]; nodecount++)
      {
        buf[m++] = nodal_positions[j][poly_index][nodecount][0] + dx;
        buf[m++] = nodal_positions[j][poly_index][nodecount][1] + dy;
        buf[m++] = nodal_positions[j][poly_index][nodecount][2] + dz;
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][0] + dx;
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][1] + dy;
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][2] + dz;
        buf[m++] = nodal_velocities[j][poly_index][nodecount][0];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][1];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][2];
      }
    }
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;

        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
      buf[m++] = ubuf(element_type[j]).d;

      buf[m++] = ubuf(element_scale[j][0]).d;
      buf[m++] = ubuf(element_scale[j][1]).d;
      buf[m++] = ubuf(element_scale[j][2]).d;

      buf[m++] = ubuf(poly_count[j]).d;
      for (int type_map = 0; type_map < poly_count[j]; type_map++) {
        buf[m++] = ubuf(node_types[j][type_map]).d;
        buf[m++] = node_charges[j][type_map];
      }

      for (int poly_index = 0; poly_index < poly_count[j]; poly_index++){
        for (int nodecount = 0; nodecount< nodes_count_list[element_type[j]]; nodecount++)
        {
        buf[m++] = nodal_positions[j][poly_index][nodecount][0] + dx;
        buf[m++] = nodal_positions[j][poly_index][nodecount][1] + dy;
        buf[m++] = nodal_positions[j][poly_index][nodecount][2] + dz;
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][0] + dx;
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][1] + dy;
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][2] + dz;
        buf[m++] = nodal_velocities[j][poly_index][nodecount][0];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][1];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][2];
        }
      }
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
      buf[m++] = ubuf(element_type[j]).d;

      buf[m++] = ubuf(element_scale[j][0]).d;
      buf[m++] = ubuf(element_scale[j][1]).d;
      buf[m++] = ubuf(element_scale[j][2]).d;

      buf[m++] = ubuf(poly_count[j]).d;
      for (int type_map = 0; type_map < poly_count[j]; type_map++) {
        buf[m++] = ubuf(node_types[j][type_map]).d;
        buf[m++] = node_charges[j][type_map];
      }

      for (int poly_index = 0; poly_index < poly_count[j]; poly_index++){
        for (int nodecount = 0; nodecount< nodes_count_list[element_type[j]]; nodecount++)
        {
        buf[m++] = nodal_positions[j][poly_index][nodecount][0] + dx;
        buf[m++] = nodal_positions[j][poly_index][nodecount][1] + dy;
        buf[m++] = nodal_positions[j][poly_index][nodecount][2] + dz;
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][0] + dx;
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][1] + dy;
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][2] + dz;
        buf[m++] = nodal_velocities[j][poly_index][nodecount][0];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][1];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][2];
        }
      }
        }
      }
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecCAC_Charge::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;
  int *nodes_count_list = atom->nodes_per_element_list;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
   element_type[i]= (int) ubuf(buf[m++]).i;
   element_scale[i][0] = (int) ubuf(buf[m++]).i;
   element_scale[i][1] = (int) ubuf(buf[m++]).i;
   element_scale[i][2] = (int) ubuf(buf[m++]).i;
   poly_count[i] = (int)ubuf(buf[m++]).i;
  for (int type_map = 0; type_map < poly_count[i]; type_map++) {
     node_types[i][type_map]= (int) ubuf(buf[m++]).i;
     node_charges[i][type_map] = buf[m++];
  }
  for (int poly_index = 0; poly_index < poly_count[i]; poly_index++){
    for (int nodecount = 0; nodecount < nodes_count_list[element_type[i]]; nodecount++)
    {
      nodal_positions[i][poly_index][nodecount][0] = buf[m++];
      nodal_positions[i][poly_index][nodecount][1] = buf[m++];
      nodal_positions[i][poly_index][nodecount][2] = buf[m++];
      initial_nodal_positions[i][poly_index][nodecount][0] = buf[m++];
      initial_nodal_positions[i][poly_index][nodecount][1] = buf[m++];
      initial_nodal_positions[i][poly_index][nodecount][2] = buf[m++];
      nodal_velocities[i][poly_index][nodecount][0] = buf[m++];
      nodal_velocities[i][poly_index][nodecount][1] = buf[m++];
      nodal_velocities[i][poly_index][nodecount][2] = buf[m++];
    }
  }
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecCAC_Charge::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;
  int *nodes_count_list = atom->nodes_per_element_list;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  element_type[i] = (int)ubuf(buf[m++]).i;
  element_scale[i][0] = (int)ubuf(buf[m++]).i;
  element_scale[i][1] = (int)ubuf(buf[m++]).i;
  element_scale[i][2] = (int)ubuf(buf[m++]).i;
  poly_count[i] = (int)ubuf(buf[m++]).i;
  for (int type_map = 0; type_map < poly_count[i]; type_map++) {
    node_types[i][type_map] = (int)ubuf(buf[m++]).i;
    node_charges[i][type_map] = buf[m++];
  }
  for (int poly_index = 0; poly_index < poly_count[i]; poly_index++){
    for (int nodecount = 0; nodecount < nodes_count_list[element_type[i]]; nodecount++)
    {
      nodal_positions[i][poly_index][nodecount][0] = buf[m++];
      nodal_positions[i][poly_index][nodecount][1] = buf[m++];
      nodal_positions[i][poly_index][nodecount][2] = buf[m++];
      initial_nodal_positions[i][poly_index][nodecount][0] = buf[m++];
      initial_nodal_positions[i][poly_index][nodecount][1] = buf[m++];
      initial_nodal_positions[i][poly_index][nodecount][2] = buf[m++];
      nodal_velocities[i][poly_index][nodecount][0] = buf[m++];
      nodal_velocities[i][poly_index][nodecount][1] = buf[m++];
      nodal_velocities[i][poly_index][nodecount][2] = buf[m++];
    }
  }
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecCAC_Charge::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecCAC_Charge::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecCAC_Charge::pack_border(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;
  double lamda_temp[3];
  double nodal_temp[3];
  int *nodes_count_list = atom->nodes_per_element_list;
  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
    buf[m++] = ubuf(element_type[j]).d;

    buf[m++] = ubuf(element_scale[j][0]).d;
    buf[m++] = ubuf(element_scale[j][1]).d;
    buf[m++] = ubuf(element_scale[j][2]).d;

    buf[m++] = ubuf(poly_count[j]).d;
    for (int type_map = 0; type_map < poly_count[j]; type_map++) {
      buf[m++] = ubuf(node_types[j][type_map]).d;
      buf[m++] = node_charges[j][type_map];
    }

    for (int poly_index = 0; poly_index < poly_count[j]; poly_index++){
      for (int nodecount = 0; nodecount< nodes_count_list[element_type[j]]; nodecount++)
      {
        buf[m++] = nodal_positions[j][poly_index][nodecount][0];
        buf[m++] = nodal_positions[j][poly_index][nodecount][1];
        buf[m++] = nodal_positions[j][poly_index][nodecount][2];
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][0];
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][1];
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][2];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][0];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][1];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][2];
      }
    }
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
    buf[m++] = ubuf(element_type[j]).d;

    buf[m++] = ubuf(element_scale[j][0]).d;
    buf[m++] = ubuf(element_scale[j][1]).d;
    buf[m++] = ubuf(element_scale[j][2]).d;

    buf[m++] = ubuf(poly_count[j]).d;
    for (int type_map = 0; type_map < poly_count[j]; type_map++) {
      buf[m++] = ubuf(node_types[j][type_map]).d;
      buf[m++] = node_charges[j][type_map];
    }

    for (int poly_index = 0; poly_index < poly_count[j]; poly_index++){
      for (int nodecount = 0; nodecount< nodes_count_list[element_type[j]]; nodecount++)
      {
        nodal_temp[0] = nodal_positions[j][poly_index][nodecount][0];
        nodal_temp[1] = nodal_positions[j][poly_index][nodecount][1];
        nodal_temp[2] = nodal_positions[j][poly_index][nodecount][2];
        if (domain->triclinic != 0) {
          domain->x2lamda(nodal_temp, lamda_temp);
          lamda_temp[0] += dx;
          lamda_temp[1] += dy;
          lamda_temp[2] += dz;
          domain->lamda2x(lamda_temp, nodal_temp);
          buf[m++] = nodal_temp[0];
          buf[m++] = nodal_temp[1];
          buf[m++] = nodal_temp[2];
        }
        else {
          buf[m++] = nodal_temp[0]+dx;
          buf[m++] = nodal_temp[1]+dy;
          buf[m++] = nodal_temp[2]+dz;
        }

        nodal_temp[0] = initial_nodal_positions[j][poly_index][nodecount][0];
        nodal_temp[1] = initial_nodal_positions[j][poly_index][nodecount][1];
        nodal_temp[2] = initial_nodal_positions[j][poly_index][nodecount][2];
        if (domain->triclinic != 0) {
          domain->x2lamda(nodal_temp, lamda_temp);
          lamda_temp[0] += dx;
          lamda_temp[1] += dy;
          lamda_temp[2] += dz;
          domain->lamda2x(lamda_temp, nodal_temp);
          buf[m++] = nodal_temp[0];
          buf[m++] = nodal_temp[1];
          buf[m++] = nodal_temp[2];
        }
        else {
          buf[m++] = nodal_temp[0] + dx;
          buf[m++] = nodal_temp[1] + dy;
          buf[m++] = nodal_temp[2] + dz;
        }
        buf[m++] = nodal_velocities[j][poly_index][nodecount][0];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][1];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][2];
      }
    }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecCAC_Charge::pack_border_vel(int n, int *list, double *buf,
                                   int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;
  double lamda_temp[3];
  double nodal_temp[3];
  int *nodes_count_list = atom->nodes_per_element_list;
  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    buf[m++] = ubuf(element_type[j]).d;

    buf[m++] = ubuf(element_scale[j][0]).d;
    buf[m++] = ubuf(element_scale[j][1]).d;
    buf[m++] = ubuf(element_scale[j][2]).d;

    buf[m++] = ubuf(poly_count[j]).d;
    for (int type_map = 0; type_map < poly_count[j]; type_map++) {
      buf[m++] = ubuf(node_types[j][type_map]).d;
      buf[m++] = node_charges[j][type_map];
    }

    for (int poly_index = 0; poly_index < poly_count[j]; poly_index++){
      for (int nodecount = 0; nodecount< nodes_count_list[element_type[j]]; nodecount++)
      {
        buf[m++] = nodal_positions[j][poly_index][nodecount][0];
        buf[m++] = nodal_positions[j][poly_index][nodecount][1];
        buf[m++] = nodal_positions[j][poly_index][nodecount][2];
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][0];
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][1];
        buf[m++] = initial_nodal_positions[j][poly_index][nodecount][2];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][0];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][1];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][2];
      }
    }
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
    buf[m++] = ubuf(element_type[j]).d;

    buf[m++] = ubuf(element_scale[j][0]).d;
    buf[m++] = ubuf(element_scale[j][1]).d;
    buf[m++] = ubuf(element_scale[j][2]).d;

    buf[m++] = ubuf(poly_count[j]).d;
    for (int type_map = 0; type_map < poly_count[j]; type_map++) {
      buf[m++] = ubuf(node_types[j][type_map]).d;
      buf[m++] = node_charges[j][type_map];
    }

    for (int poly_index = 0; poly_index < poly_count[j]; poly_index++){
      for (int nodecount = 0; nodecount< nodes_count_list[element_type[j]]; nodecount++)
      {
        nodal_temp[0] = nodal_positions[j][poly_index][nodecount][0];
        nodal_temp[1] = nodal_positions[j][poly_index][nodecount][1];
        nodal_temp[2] = nodal_positions[j][poly_index][nodecount][2];
        if (domain->triclinic != 0) {
          domain->x2lamda(nodal_temp, lamda_temp);
          lamda_temp[0] += dx;
          lamda_temp[1] += dy;
          lamda_temp[2] += dz;
          domain->lamda2x(lamda_temp, nodal_temp);
          buf[m++] = nodal_temp[0];
          buf[m++] = nodal_temp[1];
          buf[m++] = nodal_temp[2];
        }
        else {
          buf[m++] = nodal_temp[0] + dx;
          buf[m++] = nodal_temp[1] + dy;
          buf[m++] = nodal_temp[2] + dz;
        }

        nodal_temp[0] = initial_nodal_positions[j][poly_index][nodecount][0];
        nodal_temp[1] = initial_nodal_positions[j][poly_index][nodecount][1];
        nodal_temp[2] = initial_nodal_positions[j][poly_index][nodecount][2];
        if (domain->triclinic != 0) {
          domain->x2lamda(nodal_temp, lamda_temp);
          lamda_temp[0] += dx;
          lamda_temp[1] += dy;
          lamda_temp[2] += dz;
          domain->lamda2x(lamda_temp, nodal_temp);
          buf[m++] = nodal_temp[0];
          buf[m++] = nodal_temp[1];
          buf[m++] = nodal_temp[2];
        }
        else {
          buf[m++] = nodal_temp[0] + dx;
          buf[m++] = nodal_temp[1] + dy;
          buf[m++] = nodal_temp[2] + dz;
        }
        buf[m++] = nodal_velocities[j][poly_index][nodecount][0];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][1];
        buf[m++] = nodal_velocities[j][poly_index][nodecount][2];
      }
    }
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
      buf[m++] = ubuf(element_type[j]).d;

      buf[m++] = ubuf(element_scale[j][0]).d;
      buf[m++] = ubuf(element_scale[j][1]).d;
      buf[m++] = ubuf(element_scale[j][2]).d;

      buf[m++] = ubuf(poly_count[j]).d;
      for (int type_map = 0; type_map < poly_count[j]; type_map++) {
        buf[m++] = ubuf(node_types[j][type_map]).d;
        buf[m++] = node_charges[j][type_map];
      }

      for (int poly_index = 0; poly_index < poly_count[j]; poly_index++){
        for (int nodecount = 0; nodecount< nodes_count_list[element_type[j]]; nodecount++)
        {
          nodal_temp[0] = nodal_positions[j][poly_index][nodecount][0];
          nodal_temp[1] = nodal_positions[j][poly_index][nodecount][1];
          nodal_temp[2] = nodal_positions[j][poly_index][nodecount][2];
          if (domain->triclinic != 0) {
            domain->x2lamda(nodal_temp, lamda_temp);
            lamda_temp[0] += dx;
            lamda_temp[1] += dy;
            lamda_temp[2] += dz;
            domain->lamda2x(lamda_temp, nodal_temp);
            buf[m++] = nodal_temp[0];
            buf[m++] = nodal_temp[1];
            buf[m++] = nodal_temp[2];
          }
          else {
            buf[m++] = nodal_temp[0] + dx;
            buf[m++] = nodal_temp[1] + dy;
            buf[m++] = nodal_temp[2] + dz;
          }

          nodal_temp[0] = initial_nodal_positions[j][poly_index][nodecount][0];
          nodal_temp[1] = initial_nodal_positions[j][poly_index][nodecount][1];
          nodal_temp[2] = initial_nodal_positions[j][poly_index][nodecount][2];
          if (domain->triclinic != 0) {
            domain->x2lamda(nodal_temp, lamda_temp);
            lamda_temp[0] += dx;
            lamda_temp[1] += dy;
            lamda_temp[2] += dz;
            domain->lamda2x(lamda_temp, nodal_temp);
            buf[m++] = nodal_temp[0];
            buf[m++] = nodal_temp[1];
            buf[m++] = nodal_temp[2];
          }
          else {
            buf[m++] = nodal_temp[0] + dx;
            buf[m++] = nodal_temp[1] + dy;
            buf[m++] = nodal_temp[2] + dz;
          }
          buf[m++] = nodal_velocities[j][poly_index][nodecount][0] + dvx;
          buf[m++] = nodal_velocities[j][poly_index][nodecount][1] + dvy;
          buf[m++] = nodal_velocities[j][poly_index][nodecount][2] + dvz;
        }
      }
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
          buf[m++] = ubuf(element_type[j]).d;

          buf[m++] = ubuf(element_scale[j][0]).d;
          buf[m++] = ubuf(element_scale[j][1]).d;
          buf[m++] = ubuf(element_scale[j][2]).d;

          buf[m++] = ubuf(poly_count[j]).d;
      for (int type_map = 0; type_map < poly_count[j]; type_map++) {
        buf[m++] = ubuf(node_types[j][type_map]).d;
        buf[m++] = node_charges[j][type_map];
      }

      for (int poly_index = 0; poly_index < poly_count[j]; poly_index++){
        for (int nodecount = 0; nodecount< nodes_count_list[element_type[j]]; nodecount++)
        {
          nodal_temp[0] = nodal_positions[j][poly_index][nodecount][0];
          nodal_temp[1] = nodal_positions[j][poly_index][nodecount][1];
          nodal_temp[2] = nodal_positions[j][poly_index][nodecount][2];
          if (domain->triclinic != 0) {
            domain->x2lamda(nodal_temp, lamda_temp);
            lamda_temp[0] += dx;
            lamda_temp[1] += dy;
            lamda_temp[2] += dz;
            domain->lamda2x(lamda_temp, nodal_temp);
            buf[m++] = nodal_temp[0];
            buf[m++] = nodal_temp[1];
            buf[m++] = nodal_temp[2];
          }
          else {
            buf[m++] = nodal_temp[0] + dx;
            buf[m++] = nodal_temp[1] + dy;
            buf[m++] = nodal_temp[2] + dz;
          }

          nodal_temp[0] = initial_nodal_positions[j][poly_index][nodecount][0];
          nodal_temp[1] = initial_nodal_positions[j][poly_index][nodecount][1];
          nodal_temp[2] = initial_nodal_positions[j][poly_index][nodecount][2];
          if (domain->triclinic != 0) {
            domain->x2lamda(nodal_temp, lamda_temp);
            lamda_temp[0] += dx;
            lamda_temp[1] += dy;
            lamda_temp[2] += dz;
            domain->lamda2x(lamda_temp, nodal_temp);
            buf[m++] = nodal_temp[0];
            buf[m++] = nodal_temp[1];
            buf[m++] = nodal_temp[2];
          }
          else {
            buf[m++] = nodal_temp[0] + dx;
            buf[m++] = nodal_temp[1] + dy;
            buf[m++] = nodal_temp[2] + dz;
          }
          buf[m++] = nodal_velocities[j][poly_index][nodecount][0];
          buf[m++] = nodal_velocities[j][poly_index][nodecount][1];
          buf[m++] = nodal_velocities[j][poly_index][nodecount][2];
        }
      }
        }
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecCAC_Charge::unpack_border(int n, int first, double *buf)
{
  int i,m,last;
  int *nodes_count_list = atom->nodes_per_element_list;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
  element_type[i] = (int)ubuf(buf[m++]).i;
  element_scale[i][0] = (int)ubuf(buf[m++]).i;
  element_scale[i][1] = (int)ubuf(buf[m++]).i;
  element_scale[i][2] = (int)ubuf(buf[m++]).i;
  poly_count[i] = (int)ubuf(buf[m++]).i;
  allocate_element(i,nodes_count_list[element_type[i]],poly_count[i]);
  for (int type_map = 0; type_map < poly_count[i]; type_map++) {
    node_types[i][type_map] = (int)ubuf(buf[m++]).i;
    node_charges[i][type_map] = buf[m++];
  }
  for (int poly_index = 0; poly_index < poly_count[i]; poly_index++){
    for (int nodecount = 0; nodecount < nodes_count_list[element_type[i]]; nodecount++)
    {
      nodal_positions[i][poly_index][nodecount][0] = buf[m++];
      nodal_positions[i][poly_index][nodecount][1] = buf[m++];
      nodal_positions[i][poly_index][nodecount][2] = buf[m++];
      initial_nodal_positions[i][poly_index][nodecount][0] = buf[m++];
      initial_nodal_positions[i][poly_index][nodecount][1] = buf[m++];
      initial_nodal_positions[i][poly_index][nodecount][2] = buf[m++];
      nodal_velocities[i][poly_index][nodecount][0] = buf[m++];
      nodal_velocities[i][poly_index][nodecount][1] = buf[m++];
      nodal_velocities[i][poly_index][nodecount][2] = buf[m++];
    }
  }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecCAC_Charge::unpack_border_vel(int n, int first, double *buf)
{
  int i,m,last;
  int *nodes_count_list = atom->nodes_per_element_list;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  element_type[i] = (int)ubuf(buf[m++]).i;
  element_scale[i][0] = (int)ubuf(buf[m++]).i;
  element_scale[i][1] = (int)ubuf(buf[m++]).i;
  element_scale[i][2] = (int)ubuf(buf[m++]).i;
  poly_count[i] = (int)ubuf(buf[m++]).i;
  allocate_element(i,nodes_count_list[element_type[i]],poly_count[i]);
  for (int type_map = 0; type_map < poly_count[i]; type_map++) {
    node_types[i][type_map] = (int)ubuf(buf[m++]).i;
    node_charges[i][type_map] = buf[m++];
  }
  for (int poly_index = 0; poly_index < poly_count[i]; poly_index++){
    for (int nodecount = 0; nodecount < nodes_count_list[element_type[i]]; nodecount++)
    {
      nodal_positions[i][poly_index][nodecount][0] = buf[m++];
      nodal_positions[i][poly_index][nodecount][1] = buf[m++];
      nodal_positions[i][poly_index][nodecount][2] = buf[m++];
      initial_nodal_positions[i][poly_index][nodecount][0] = buf[m++];
      initial_nodal_positions[i][poly_index][nodecount][1] = buf[m++];
      initial_nodal_positions[i][poly_index][nodecount][2] = buf[m++];
      nodal_velocities[i][poly_index][nodecount][0] = buf[m++];
      nodal_velocities[i][poly_index][nodecount][1] = buf[m++];
      nodal_velocities[i][poly_index][nodecount][2] = buf[m++];
    }
  }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecCAC_Charge::pack_exchange(int i, double *buf)
{
  int m = 1;
  int *nodes_count_list = atom->nodes_per_element_list;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = ubuf(element_type[i]).d;

  buf[m++] = ubuf(element_scale[i][0]).d;
  buf[m++] = ubuf(element_scale[i][1]).d;
  buf[m++] = ubuf(element_scale[i][2]).d;

  buf[m++] = ubuf(poly_count[i]).d;
  for (int type_map = 0; type_map < poly_count[i]; type_map++) {
    buf[m++] = ubuf(node_types[i][type_map]).d;
    buf[m++] = node_charges[i][type_map];
  }

  for (int poly_index = 0; poly_index < poly_count[i]; poly_index++){
    for (int nodecount = 0; nodecount< nodes_count_list[element_type[i]]; nodecount++)
    {
      buf[m++] = nodal_positions[i][poly_index][nodecount][0];
      buf[m++] = nodal_positions[i][poly_index][nodecount][1];
      buf[m++] = nodal_positions[i][poly_index][nodecount][2];
      buf[m++] = initial_nodal_positions[i][poly_index][nodecount][0];
      buf[m++] = initial_nodal_positions[i][poly_index][nodecount][1];
      buf[m++] = initial_nodal_positions[i][poly_index][nodecount][2];
      buf[m++] = nodal_velocities[i][poly_index][nodecount][0];
      buf[m++] = nodal_velocities[i][poly_index][nodecount][1];
      buf[m++] = nodal_velocities[i][poly_index][nodecount][2];
    }
  }

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecCAC_Charge::unpack_exchange(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);
  int *nodes_count_list = atom->nodes_per_element_list;
  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  element_type[nlocal] = (int)ubuf(buf[m++]).i;
  element_scale[nlocal][0] = (int)ubuf(buf[m++]).i;
  element_scale[nlocal][1] = (int)ubuf(buf[m++]).i;
  element_scale[nlocal][2] = (int)ubuf(buf[m++]).i;
  poly_count[nlocal] = (int)ubuf(buf[m++]).i;
  allocate_element(nlocal,nodes_count_list[element_type[nlocal]],poly_count[nlocal]);
  for (int type_map = 0; type_map < poly_count[nlocal]; type_map++) {
    node_types[nlocal][type_map] = (int)ubuf(buf[m++]).i;
    node_charges[nlocal][type_map] = buf[m++];
  }
  for (int poly_index = 0; poly_index < poly_count[nlocal]; poly_index++){
    for (int nodecount = 0; nodecount < nodes_count_list[element_type[nlocal]]; nodecount++)
    {
      nodal_positions[nlocal][poly_index][nodecount][0] = buf[m++];
      nodal_positions[nlocal][poly_index][nodecount][1] = buf[m++];
      nodal_positions[nlocal][poly_index][nodecount][2] = buf[m++];
      initial_nodal_positions[nlocal][poly_index][nodecount][0] = buf[m++];
      initial_nodal_positions[nlocal][poly_index][nodecount][1] = buf[m++];
      initial_nodal_positions[nlocal][poly_index][nodecount][2] = buf[m++];
      nodal_velocities[nlocal][poly_index][nodecount][0] = buf[m++];
      nodal_velocities[nlocal][poly_index][nodecount][1] = buf[m++];
      nodal_velocities[nlocal][poly_index][nodecount][2] = buf[m++];
    }
  }

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecCAC_Charge::size_restart()
{
  int i;
  int current_node_count;
  int *nodes_count_list = atom->nodes_per_element_list;
  int nlocal = atom->nlocal;
  int n=0;
  for (i=0; i < nlocal; i++){
  current_node_count=nodes_count_list[element_type[i]];
   n += (16+9*current_node_count*poly_count[i]+2*poly_count[i]);
  }

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecCAC_Charge::pack_restart(int i, double *buf)
{
  int m = 1;
  int current_node_count;
  int *nodes_count_list = atom->nodes_per_element_list;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = ubuf(element_type[i]).d;
  buf[m++] = ubuf(element_scale[i][0]).d;
  buf[m++] = ubuf(element_scale[i][1]).d;
  buf[m++] = ubuf(element_scale[i][2]).d;
  buf[m++] = ubuf(poly_count[i]).d;
  current_node_count=nodes_count_list[element_type[i]];

  for (int type_map = 0; type_map < poly_count[i]; type_map++) {
    buf[m++] = ubuf(node_types[i][type_map]).d;
    buf[m++] = node_charges[i][type_map];
  }

  for (int poly_index = 0; poly_index < poly_count[i]; poly_index++) {
    for (int nodecount = 0; nodecount< current_node_count; nodecount++)
    {
      buf[m++] = nodal_positions[i][poly_index][nodecount][0];
      buf[m++] = nodal_positions[i][poly_index][nodecount][1];
      buf[m++] = nodal_positions[i][poly_index][nodecount][2];
      buf[m++] = initial_nodal_positions[i][poly_index][nodecount][0];
      buf[m++] = initial_nodal_positions[i][poly_index][nodecount][1];
      buf[m++] = initial_nodal_positions[i][poly_index][nodecount][2];
      buf[m++] = nodal_velocities[i][poly_index][nodecount][0];
      buf[m++] = nodal_velocities[i][poly_index][nodecount][1];
      buf[m++] = nodal_velocities[i][poly_index][nodecount][2];
    }
  }
  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecCAC_Charge::unpack_restart(double *buf)
{
  int current_node_count;
  int nlocal = atom->nlocal;

  int *nodes_count_list = atom->nodes_per_element_list;
  initial_size=atom->initial_size;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  element_type[nlocal] = (int) ubuf(buf[m++]).i;
  element_scale[nlocal][0] = (int) ubuf(buf[m++]).i;
  element_scale[nlocal][1] = (int) ubuf(buf[m++]).i;
  element_scale[nlocal][2] = (int) ubuf(buf[m++]).i;
  poly_count[nlocal] = (int) ubuf(buf[m++]).i;
  current_node_count=nodes_count_list[element_type[nlocal]];
  allocate_element(nlocal,current_node_count,poly_count[nlocal]);

  for (int type_map = 0; type_map < poly_count[nlocal]; type_map++) {
    node_types[nlocal][type_map] = (int) ubuf(buf[m++]).i;
    node_charges[nlocal][type_map] = buf[m++];
  }
  for (int poly_index = 0; poly_index < poly_count[nlocal]; poly_index++){
    for (int nodecount = 0; nodecount < current_node_count; nodecount++)
    {
      nodal_positions[nlocal][poly_index][nodecount][0] = buf[m++];
      nodal_positions[nlocal][poly_index][nodecount][1] = buf[m++];
      nodal_positions[nlocal][poly_index][nodecount][2] = buf[m++];
      initial_nodal_positions[nlocal][poly_index][nodecount][0] = buf[m++];
      initial_nodal_positions[nlocal][poly_index][nodecount][1] = buf[m++];
      initial_nodal_positions[nlocal][poly_index][nodecount][2] = buf[m++];
      nodal_velocities[nlocal][poly_index][nodecount][0] = buf[m++];
      nodal_velocities[nlocal][poly_index][nodecount][1] = buf[m++];
      nodal_velocities[nlocal][poly_index][nodecount][2] = buf[m++];
    }
  }


  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecCAC_Charge::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) |
    ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  element_type[nlocal] = 0;

  poly_count[nlocal] =1;
  allocate_element(nlocal,1,poly_count[nlocal]);
  for (int type_map = 0; type_map < poly_count[nlocal]; type_map++) {
    node_types[nlocal][type_map] = itype;
    node_charges[nlocal][type_map] = 0;
  }
  for (int poly_index = 0; poly_index < 1; poly_index++){
    for (int nodecount = 0; nodecount < 1; nodecount++)
    {
      nodal_positions[nlocal][poly_index][nodecount][0] = coord[0];
      nodal_positions[nlocal][poly_index][nodecount][1] = coord[1];
      nodal_positions[nlocal][poly_index][nodecount][2] = coord[2];
      initial_nodal_positions[nlocal][poly_index][nodecount][0] = coord[0];
      initial_nodal_positions[nlocal][poly_index][nodecount][1] = coord[1];
      initial_nodal_positions[nlocal][poly_index][nodecount][2] = coord[2];
      nodal_velocities[nlocal][poly_index][nodecount][0] = 0;
      nodal_velocities[nlocal][poly_index][nodecount][1] = 0;
      nodal_velocities[nlocal][poly_index][nodecount][2] = 0;
    }
  }

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one element from CAC_Elements section of data file
   initialize other Element quantities
------------------------------------------------------------------------- */

void AtomVecCAC_Charge::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  int node_index,node_type,poly_index;
  int typefound;
  double node_charge;
  if (nlocal == nmax) grow(0);
  int nodetotal, npoly;
  int tmp;

  int *nodes_count_list = atom->nodes_per_element_list;

  int types_filled = 0;
    initial_size=atom->initial_size;

  poly_index = 0;
  tag[nlocal] = ATOTAGINT(values[0]);
  char* element_type_read;
  element_type_read = values[1];
  type[nlocal] = 1;

  npoly = force->inumeric(FLERR,values[2]);
  if (npoly > maxpoly)
    error->one(FLERR, "poly count declared in data file was greater than maxpoly in input file");

//loop through defined element types
  int type_found=0;
  for(int string_check=1; string_check < element_type_count; string_check++){
    if (strcmp(element_type_read, element_names[string_check]) == 0){
    type_found=1;
    element_type[nlocal] = string_check;
    nodetotal = nodes_count_list[string_check];
    poly_count[nlocal] = npoly;
    element_scale[nlocal][0] = force->inumeric(FLERR,values[3]);
    element_scale[nlocal][1] = force->inumeric(FLERR,values[4]);
    element_scale[nlocal][2] = force->inumeric(FLERR,values[5]);
    }
  }
  //if (strcmp(element_type_read, "Eight_Node") == 0) {//add a control block for new types of elements

  //}
  //set atom type explicitly in case default values werent set to convention
  if (strcmp(element_type_read, "Atom") == 0) {
    type_found=1;
    element_type[nlocal] = 0;
    nodetotal = nodes_count_list[element_type[nlocal]];
    npoly = 1;
    poly_count[nlocal] = npoly;
    element_scale[nlocal][0] = 1;
    element_scale[nlocal][1] = 1;
    element_scale[nlocal][2] = 1;

  }

  if(!type_found) {
    error->one(FLERR, "element type not yet defined, add definition in process_args function of atom_vec_CAC.cpp style");
  }
  if (nodetotal > nodes_per_element)
    error->one(FLERR, "element type requires a greater number of nodes than the specified maximum nodes per element passed to atom style cac/charge");

  allocate_element(nlocal,nodetotal,poly_count[nlocal]);
  for (int polycount = 0; polycount < npoly; polycount++) {
    node_types[nlocal][polycount] = 0; //initialize
    node_charges[nlocal][polycount] = 0; //initialize
    node_count_per_poly[polycount]=0;
  }


  int m = 6;

  for (int polycount = 0; polycount < npoly; polycount++) {
    for (int nodecount = 0; nodecount < nodetotal; nodecount++)
    {


    node_index = force->inumeric(FLERR,values[m++]);
    if (node_index < 1 ||node_index > nodetotal)
      error->one(FLERR, "Invalid node index in CAC_Elements section of data file");
    poly_index = force->inumeric(FLERR,values[m++]);
    if (poly_index < 1 || poly_index > npoly)
      error->one(FLERR, "Invalid poly index in CAC_Elements section of data file");
    node_index = node_index - 1;
    poly_index = poly_index - 1;
    node_type = force->inumeric(FLERR,values[m++]);
    node_charge = force->numeric(FLERR,values[m++]);
    node_count_per_poly[poly_index]++;
    if (node_type <= 0 || node_type > atom->ntypes)
      error->one(FLERR, "Invalid atom type in CAC_Elements section of data file");

    if (node_types[nlocal][poly_index] == 0 || node_types[nlocal][poly_index] == node_type) {
      node_types[nlocal][poly_index] = node_type;
      node_charges[nlocal][poly_index] = node_charge;
    }
    else {
      error->one(FLERR, "more than one type assigned to the same poly index in an element");
    }

    if(node_count_per_poly[poly_index]>nodetotal)
    error->one(FLERR, "there are more nodes for one internal DOF than the element type admits");

    nodal_positions[nlocal][poly_index][node_index][0] = force->numeric(FLERR,values[m++]);
    nodal_positions[nlocal][poly_index][node_index][1] = force->numeric(FLERR,values[m++]);
    nodal_positions[nlocal][poly_index][node_index][2] = force->numeric(FLERR,values[m++]);
    initial_nodal_positions[nlocal][poly_index][node_index][0] = nodal_positions[nlocal][poly_index][node_index][0];
    initial_nodal_positions[nlocal][poly_index][node_index][1] = nodal_positions[nlocal][poly_index][node_index][1];
    initial_nodal_positions[nlocal][poly_index][node_index][2] = nodal_positions[nlocal][poly_index][node_index][2];
    nodal_velocities[nlocal][poly_index][node_index][0] = 0;
    nodal_velocities[nlocal][poly_index][node_index][1] = 0;
    nodal_velocities[nlocal][poly_index][node_index][2] = 0;
    nodal_forces[nlocal][poly_index][node_index][0] = 0;
    nodal_forces[nlocal][poly_index][node_index][1] = 0;
    nodal_forces[nlocal][poly_index][node_index][2] = 0;
  }
  }




  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  atom->nlocal++;

}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecCAC_Charge::pack_data(double **buf)
{
  error->all(FLERR,"cac atom style does not yet support writing data files");
  int nlocal = atom->nlocal;
  int *nodes_count_list = atom->nodes_per_element_list;
  for (int i = 0; i < nlocal; i++) {
       int m=0;
    buf[i][m++] = ubuf(tag[i]).d;
    buf[i][m++] = ubuf(type[i]).d;
  buf[i][m++] = ubuf(element_type[i]).d;

  buf[i][m++] = ubuf(element_scale[i][0]).d;
  buf[i][m++] = ubuf(element_scale[i][1]).d;
  buf[i][m++] = ubuf(element_scale[i][2]).d;

  buf[i][m++] = ubuf(poly_count[i]).d;
  for (int type_map = 0; type_map < poly_count[i]; type_map++) {
    buf[i][m++] = ubuf(node_types[i][type_map]).d;
    buf[i][m++] = node_charges[i][type_map];
  }

  for (int poly_index = 0; poly_index < poly_count[i]; poly_index++){
    for (int nodecount = 0; nodecount< nodes_count_list[element_type[i]]; nodecount++)
    {
      buf[i][m++] = nodal_positions[i][poly_index][nodecount][0];
      buf[i][m++] = nodal_positions[i][poly_index][nodecount][1];
      buf[i][m++] = nodal_positions[i][poly_index][nodecount][2];
      buf[i][m++] = initial_nodal_positions[i][poly_index][nodecount][0];
      buf[i][m++] = initial_nodal_positions[i][poly_index][nodecount][1];
      buf[i][m++] = initial_nodal_positions[i][poly_index][nodecount][2];
      buf[i][m++] = nodal_velocities[i][poly_index][nodecount][0];
      buf[i][m++] = nodal_velocities[i][poly_index][nodecount][1];
      buf[i][m++] = nodal_velocities[i][poly_index][nodecount][2];
    }
  }
    buf[i][m++] = x[i][0];
    buf[i][m++] = x[i][1];
    buf[i][m++] = x[i][2];
    buf[i][m++] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][m++] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][m++] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecCAC_Charge::write_data(FILE *fp, int n, double **buf)
{
  error->all(FLERR,"cac atom style does not yet support writing data files");
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            buf[i][2],buf[i][3],buf[i][4],
            (int) ubuf(buf[i][5]).i,(int) ubuf(buf[i][6]).i,
            (int) ubuf(buf[i][7]).i);
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecCAC_Charge::memory_usage()
{
  bigint bytes = 0;
  int *nodes_count_list = atom->nodes_per_element_list;
  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);
  if (atom->memcheck("element_types")) bytes += memory->usage(element_type, nmax);
  if (atom->memcheck("poly_counts")) bytes += memory->usage(poly_count, nmax);
  if (atom->memcheck("element_scale")) bytes += memory->usage(element_scale, nmax, 3);
  for(int usage_index=0; usage_index < alloc_counter; usage_index++){
  int current_poly_count = poly_count[usage_index];
  int node_count = nodes_count_list[element_type[usage_index]];
  if (atom->memcheck("node_types")) bytes += memory->usage(node_types[usage_index],current_poly_count);
  if (atom->memcheck("node_charges")) bytes += memory->usage(node_charges[usage_index], current_poly_count);
  if (atom->memcheck("nodal_positions")) bytes += memory->usage(nodal_positions[usage_index],current_poly_count, node_count,3);
  if (atom->memcheck("hold_nodal_positions")) bytes += memory->usage(hold_nodal_positions[usage_index],current_poly_count, node_count,3);
  if (atom->memcheck("initial_nodal_positions")) bytes += memory->usage(initial_nodal_positions[usage_index], current_poly_count, node_count, 3);
  if (atom->memcheck("nodal_velocities")) bytes += memory->usage(nodal_velocities[usage_index], current_poly_count, node_count,3);
  if (atom->memcheck("nodal_forces")) bytes += memory->usage(nodal_forces[usage_index], current_poly_count, node_count,3);
  if (atom->memcheck("nodal_virial")) bytes += memory->usage(nodal_virial[usage_index], current_poly_count, node_count,6);
  }

  return bytes;
}