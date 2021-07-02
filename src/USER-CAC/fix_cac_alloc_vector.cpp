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
#include "fix_cac_alloc_vector.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCACAllocVector::FixCACAllocVector(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nvector(0), peratom(NULL), nodal_vectors(NULL)
{
  if(narg==4)
    callback = utils::inumeric(FLERR,arg[3],false,lmp);
  else
    callback = 0;
  // register callback to this fix from Atom class
  // don't perform initial allocation here, must wait until add_vector()
  if(callback==0)
    atom->add_callback(0);
  else if(callback==3){
    atom->add_callback(0);
    atom->add_callback(3);
  }
  else{
    error->all(FLERR,"fix CAC/ALLOCVECTOR was invoked with an unimplemented callback flag");
  }
  alloc_counter=0;
}

/* ---------------------------------------------------------------------- */

FixCACAllocVector::~FixCACAllocVector()
{
  // unregister callbacks to this fix from Atom class
  if(callback==0)
    atom->delete_callback(id,0);
  else if(callback==3){
    atom->delete_callback(id,0);
    atom->delete_callback(id,3);
  }

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

int FixCACAllocVector::setmask()
{
  return 0;
}

/* ----------------------------------------------------------------------
   allocate/initialize memory for a new vector with N elements per atom
------------------------------------------------------------------------- */

void FixCACAllocVector::add_vector(int n)
{
  int nlocal = atom->nlocal;
  int *nodes_per_element_list = atom->nodes_per_element_list;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int node_count;
  memory->grow(peratom,nvector+1,"CACAllocVector:peratom");
  peratom[nvector] = n;

  nodal_vectors = (double *****)
    memory->srealloc(nodal_vectors,(nvector+1)*sizeof(double ****),"CACAllocVector:nodal_vectors");
  nodal_vectors[nvector] = (double ****)
    memory->smalloc(atom->nmax*sizeof(double ***),"CACAllocVector:nodal_vector");

  //allocate element sizing for the additional vector
  for(int ecount=0; ecount < nlocal; ecount++){
  node_count = nodes_per_element_list[element_type[ecount]];
  memory->create(nodal_vectors[nvector][ecount], poly_count[ecount], node_count, peratom[nvector], "CACAllocVector:vectors");
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

  //boost max exchange variables for communication buffer
  atom->avec->maxexchange += n*atom->nodes_per_element*atom->maxpoly;

  nvector++;
}

/* ----------------------------------------------------------------------
   return a pointer to the Mth vector
------------------------------------------------------------------------- */

double ****FixCACAllocVector::request_vector(int m)
{
  return nodal_vectors[m];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixCACAllocVector::memory_usage()
{
  double bytes = 0.0;
  for (int m = 0; m < nvector; m++)
    bytes += atom->nmax*peratom[m]*sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixCACAllocVector::grow_arrays(int nmax)
{
  for (int m = 0; m < nvector; m++)
  nodal_vectors[m] =
    (double ****) memory->srealloc(nodal_vectors[m],sizeof(double ***)*nmax, "CACAllocVector:nodal_vectors");
}

/* ----------------------------------------------------------------------
   shrink nodal arrays
------------------------------------------------------------------------- */

void FixCACAllocVector::shrink_arrays(int n)
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
    (double ****) memory->srealloc(nodal_vectors[m],sizeof(double ***)*n, "CACAllocVector:nodal_vectors");
}

/* ----------------------------------------------------------------------
   allocate finite element data for nodal arrays
------------------------------------------------------------------------- */

void FixCACAllocVector::allocate_element(int element_index, int node_count, int poly_count)
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
  memory->create(nodal_vectors[m][element_index], poly_count, node_count, peratom[m], "CACAllocVector:nodal_vectors");

}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixCACAllocVector::copy_arrays(int i, int j, int /*delflag*/)
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
  memory->create(nodal_vectors[m][j], poly_count[j], node_count, peratom[m], "CACAllocVector:nodal_vectors");
  
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

int FixCACAllocVector::pack_exchange(int i, double *buf)
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

int FixCACAllocVector::unpack_exchange(int nlocal, double *buf)
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

/* ----------------------------------------------------------------------
   clear nodal array values such as forces and fluxes every timestep
------------------------------------------------------------------------- */

void FixCACAllocVector::clear_arrays(int a, size_t b)
{
  int m,iper,nper,ni;
  int nlocal = atom->nlocal;
  int *nodes_count_list = atom->nodes_per_element_list;
  int *poly_count = atom->poly_count;
  int *element_type = atom->element_type;

  for (m = 0; m < nvector; m++) {
    nper = peratom[m];
    for (int i = 0; i < atom->nlocal; i++)
      for (int poly_index = 0; poly_index < poly_count[i]; poly_index++)
        for (int nodecount = 0; nodecount < nodes_count_list[element_type[i]]; nodecount++)
          for (iper = 0; iper < nper; iper++) nodal_vectors[m][i][poly_index][nodecount][iper] = 0;
  }
}