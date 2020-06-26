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
#include "fix_cac_minimize.h"
#include "atom.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCACMinimize::FixCACMinimize(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nvector(0), peratom(NULL), nodal_vectors(NULL)
{
  // register callback to this fix from Atom class
  // don't perform initial allocation here, must wait until add_vector()

  atom->add_callback(0);
  alloc_counter=0;
}

/* ---------------------------------------------------------------------- */

FixCACMinimize::~FixCACMinimize()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  // delete locally stored data

  memory->destroy(peratom);
  for (int m = 0; m < nvector; m++){
  for(int element_index=0; element_index < alloc_counter; element_index++)
    memory->destroy(nodal_vectors[m][element_index]);
  memory->sfree(nodal_vectors[m]);
  }
  memory->sfree(nodal_vectors);
}

/* ---------------------------------------------------------------------- */

int FixCACMinimize::setmask()
{
  return 0;
}

/* ----------------------------------------------------------------------
   allocate/initialize memory for a new vector with N elements per atom
------------------------------------------------------------------------- */

void FixCACMinimize::add_vector(int n)
{
  int nlocal = atom->nlocal;
  int *nodes_per_element_list = atom->nodes_per_element_list;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int node_count;
  memory->grow(peratom,nvector+1,"minimize:peratom");
  peratom[nvector] = n;

  nodal_vectors = (double *****)
    memory->srealloc(nodal_vectors,(nvector+1)*sizeof(double ****),"minimize:vectors");
  nodal_vectors[nvector] = (double ****)
    memory->smalloc(atom->nmax*sizeof(double ***),"minimize:vector");

  //allocate element sizing for the additional vector
  for(int ecount=0; ecount < nlocal; ecount++){
  node_count = nodes_per_element_list[element_type[ecount]];
  memory->create(nodal_vectors[nvector][ecount], poly_count[ecount], node_count, peratom[nvector], "minimize:vectors");
  }
  if(nlocal>alloc_counter) alloc_counter = nlocal;
  //initialize to zero
  for(int ecount=0; ecount < nlocal; ecount++){
  node_count = nodes_per_element_list[element_type[ecount]];
    for(int ipoly = 0; ipoly < poly_count[ecount]; ipoly++)
      for(int inode = 0; inode < node_count; inode++)
        for(int ix = 0; ix < peratom[nvector]; ix++)
        nodal_vectors[nvector][ecount][ipoly][inode][ix] = 0;
  }
  nvector++;
}

/* ----------------------------------------------------------------------
   return a pointer to the Mth vector
------------------------------------------------------------------------- */

double ****FixCACMinimize::request_vector(int m)
{
  return nodal_vectors[m];
}

/* ----------------------------------------------------------------------
   store box size at beginning of line search
------------------------------------------------------------------------- */

void FixCACMinimize::store_box()
{
  boxlo[0] = domain->boxlo[0];
  boxlo[1] = domain->boxlo[1];
  boxlo[2] = domain->boxlo[2];
  boxhi[0] = domain->boxhi[0];
  boxhi[1] = domain->boxhi[1];
  boxhi[2] = domain->boxhi[2];
}

/* ----------------------------------------------------------------------
   reset x0 for atoms that moved across PBC via reneighboring in line search
   x0 = 1st vector
   must do minimum_image using original box stored at beginning of line search
   swap & set_global_box() change to original box, then restore current box
------------------------------------------------------------------------- */

void FixCACMinimize::reset_coords()
{
  int *poly_count = atom->poly_count;
  int *nodes_count_list = atom->nodes_per_element_list;
  int *element_type = atom->element_type;
  int node_count;

  box_swap();
  domain->set_global_box();

  double ****nodal_x = atom->nodal_positions;
  double ****nodal_x0 = nodal_vectors[0];
  int nlocal = atom->nlocal;
  double dx,dy,dz,dx0,dy0,dz0;

  for (int i = 0; i < nlocal; i++) {
    dx = dx0 = nodal_x[i][0][0][0] - nodal_x0[i][0][0][0];
    dy = dy0 = nodal_x[i][0][0][1] - nodal_x0[i][0][0][1];
    dz = dz0 = nodal_x[i][0][0][2] - nodal_x0[i][0][0][2];
    domain->minimum_image(dx,dy,dz);
    if (dx != dx0){
      //remap nodal information
      node_count = nodes_count_list[element_type[i]];
      for(int ipoly=0; ipoly < poly_count[i]; ipoly++){
        for(int inode=0; inode < node_count; inode++){
          nodal_x0[i][ipoly][inode][0] = nodal_x[i][ipoly][inode][0] - dx;
        }
      }
    }
    if (dy != dy0) {
      //remap nodal information
      node_count = nodes_count_list[element_type[i]];
      for(int ipoly=0; ipoly < poly_count[i]; ipoly++){
        for(int inode=0; inode < node_count; inode++){
          nodal_x0[i][ipoly][inode][1] = nodal_x[i][ipoly][inode][1] - dy;
        }
      }
    }
    if (dz != dz0) {
      //remap nodal information
      node_count = nodes_count_list[element_type[i]];
      for(int ipoly=0; ipoly < poly_count[i]; ipoly++){
        for(int inode=0; inode < node_count; inode++){
          nodal_x0[i][ipoly][inode][2] = nodal_x[i][ipoly][inode][2] - dz;
        }
      }
    }
  }

  box_swap();
  domain->set_global_box();
}

/* ----------------------------------------------------------------------
   swap current box size with stored box size
------------------------------------------------------------------------- */

void FixCACMinimize::box_swap()
{
  double tmp;

  tmp = boxlo[0];
  boxlo[0] = domain->boxlo[0];
  domain->boxlo[0] = tmp;
  tmp = boxlo[1];
  boxlo[1] = domain->boxlo[1];
  domain->boxlo[1] = tmp;
  tmp = boxlo[2];
  boxlo[2] = domain->boxlo[2];
  domain->boxlo[2] = tmp;

  tmp = boxhi[0];
  boxhi[0] = domain->boxhi[0];
  domain->boxhi[0] = tmp;
  tmp = boxhi[1];
  boxhi[1] = domain->boxhi[1];
  domain->boxhi[1] = tmp;
  tmp = boxhi[2];
  boxhi[2] = domain->boxhi[2];
  domain->boxhi[2] = tmp;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixCACMinimize::memory_usage()
{
  double bytes = 0.0;
  for (int m = 0; m < nvector; m++)
    bytes += atom->nmax*peratom[m]*sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixCACMinimize::grow_arrays(int nmax)
{
  for (int m = 0; m < nvector; m++)
  nodal_vectors[m] =
    (double ****) memory->srealloc(nodal_vectors[m],sizeof(double ***)*nmax, "atom:nodal_positions");
}

/* ----------------------------------------------------------------------
   shrink nodal arrays
------------------------------------------------------------------------- */

void FixCACMinimize::shrink_arrays(int n)
{
  //deallocate element contents if n is smaller than the alloc counter for elements
  for(int element_index=alloc_counter-1; element_index >= n; element_index--){
  for (int m = 0; m < nvector; m++)
  memory->destroy(nodal_vectors[m][element_index]);
  }
  if(alloc_counter>n)
  alloc_counter = n;

  //shrink pointer arrays
  for (int m = 0; m < nvector; m++)
  nodal_vectors[m] =
    (double ****) memory->srealloc(nodal_vectors[m],sizeof(double ***)*n, "atom:nodal_positions");
}

/* ----------------------------------------------------------------------
   allocate finite element data for nodal arrays
------------------------------------------------------------------------- */

void FixCACMinimize::allocate_element(int element_index, int node_count, int poly_count)
{
  //destroy previous contents at that element index if present
  if(element_index<alloc_counter){
  for (int m = 0; m < nvector; m++)
  memory->destroy(nodal_vectors[m][element_index]);
  }
  else
  alloc_counter++;
  //create new allocation for this element
  for (int m = 0; m < nvector; m++)
  memory->create(nodal_vectors[m][element_index], poly_count, node_count, peratom[m], "atom:nodal_positions");

}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixCACMinimize::copy_arrays(int i, int j, int /*delflag*/)
{
  int m,iper,nper,ni,nj;
  int *poly_count = atom->poly_count;
  int *nodes_count_list = atom->nodes_per_element_list;
  int *element_type = atom->element_type;
  int node_count;

  node_count = nodes_count_list[element_type[j]];
  //copy nodal information; requires resizing since copy might be for
  //an element of a different size
    //destroy previous contents at that element index if present
  if(j<alloc_counter){
  for (m = 0; m < nvector; m++)
  memory->destroy(nodal_vectors[m][j]);
  }
  else
  alloc_counter++;

  //create new allocation for this element
  for (m = 0; m < nvector; m++)
  memory->create(nodal_vectors[m][j], poly_count[j], node_count, peratom[m], "atom:nodal_positions");
  
  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    for (int poly_index = 0; poly_index < poly_count[j]; poly_index++){
      for(int nodecount=0; nodecount< nodes_count_list[element_type[j]]; nodecount++ )
      {
        for (iper = 0; iper < nper; iper++)
        nodal_vectors[m][j][poly_index][nodecount][iper] = nodal_vectors[m][i][poly_index][nodecount][iper];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixCACMinimize::pack_exchange(int i, double *buf)
{
  int m,iper,nper,ni;
  int *nodes_count_list = atom->nodes_per_element_list;
  int *poly_count = atom->poly_count;
  int *element_type = atom->element_type;
  int n = 0;

  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    for (int poly_index = 0; poly_index < poly_count[i]; poly_index++){
      for (int nodecount = 0; nodecount < nodes_count_list[element_type[i]]; nodecount++)
      {
        for (iper = 0; iper < nper; iper++) buf[n++] = nodal_vectors[m][i][poly_index][nodecount][iper];
      }
    }
  }

  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixCACMinimize::unpack_exchange(int nlocal, double *buf)
{
  int m,iper,nper,ni;
  int *nodes_count_list = atom->nodes_per_element_list;
  int *poly_count = atom->poly_count;
  int *element_type = atom->element_type;
  int n = 0;
  allocate_element(nlocal,nodes_count_list[element_type[nlocal]],poly_count[nlocal]);
  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    for (int poly_index = 0; poly_index < poly_count[nlocal]; poly_index++){
      for (int nodecount = 0; nodecount < nodes_count_list[element_type[nlocal]]; nodecount++)
      {
        for (iper = 0; iper < nper; iper++) nodal_vectors[m][nlocal][poly_index][nodecount][iper] = buf[n++];
      }
    }
  }

  return n;
}
