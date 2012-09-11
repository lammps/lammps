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

// C or Fortran style library interface to LAMMPS
// customize by adding new LAMMPS-specific functions

#include "lmptype.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "library.h"
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   create an instance of LAMMPS and return pointer to it
   pass in command-line args and MPI communicator to run on
------------------------------------------------------------------------- */

void lammps_open(int argc, char **argv, MPI_Comm communicator, void **ptr)
{
  LAMMPS *lmp = new LAMMPS(argc,argv,communicator);
  *ptr = (void *) lmp;
}

/* ----------------------------------------------------------------------
   create an instance of LAMMPS and return pointer to it
   caller doesn't know MPI communicator, so use MPI_COMM_WORLD
   intialize MPI if needed
------------------------------------------------------------------------- */

void lammps_open_no_mpi(int argc, char **argv, void **ptr)
{
  int flag;
  MPI_Initialized(&flag);

  if (!flag) {
    int argc = 0;
    char **argv = NULL;
    MPI_Init(&argc,&argv);
  }

  MPI_Comm communicator = MPI_COMM_WORLD;

  LAMMPS *lmp = new LAMMPS(argc,argv,communicator);
  *ptr = (void *) lmp;
}

/* ----------------------------------------------------------------------
   destruct an instance of LAMMPS
------------------------------------------------------------------------- */

void lammps_close(void *ptr)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  delete lmp;
}

/* ----------------------------------------------------------------------
   process an input script in filename str
------------------------------------------------------------------------- */

void lammps_file(void *ptr, char *str)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  lmp->input->file(str);
}

/* ----------------------------------------------------------------------
   process a single input command in str
------------------------------------------------------------------------- */

char *lammps_command(void *ptr, char *str)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  return lmp->input->one(str);
}

/* ----------------------------------------------------------------------
   clean-up function to free memory allocated by lib and returned to caller
------------------------------------------------------------------------- */

void lammps_free(void *ptr)
{
  free(ptr);
}

/* ----------------------------------------------------------------------
   add LAMMPS-specific library functions
   all must receive LAMMPS pointer as argument
   customize by adding a function here and in library.h header file
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   extract a pointer to an internal LAMMPS global entity
   name = desired quantity, e.g. dt or boxyhi or natoms
   returns a void pointer to the entity
     which the caller can cast to the proper data type
   returns a NULL if name not listed below
   customize by adding names
------------------------------------------------------------------------- */

void *lammps_extract_global(void *ptr, char *name)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  if (strcmp(name,"dt") == 0) return (void *) &lmp->update->dt;
  if (strcmp(name,"boxxlo") == 0) return (void *) &lmp->domain->boxlo[0];
  if (strcmp(name,"boxxhi") == 0) return (void *) &lmp->domain->boxhi[0];
  if (strcmp(name,"boxylo") == 0) return (void *) &lmp->domain->boxlo[1];
  if (strcmp(name,"boxyhi") == 0) return (void *) &lmp->domain->boxhi[1];
  if (strcmp(name,"boxzlo") == 0) return (void *) &lmp->domain->boxlo[2];
  if (strcmp(name,"boxzhi") == 0) return (void *) &lmp->domain->boxhi[2];
  if (strcmp(name,"natoms") == 0) return (void *) &lmp->atom->natoms;
  if (strcmp(name,"nlocal") == 0) return (void *) &lmp->atom->nlocal;
  return NULL;
}

/* ----------------------------------------------------------------------
   extract a pointer to an internal LAMMPS atom-based entity
   name = desired quantity, e.g. x or mass
   returns a void pointer to the entity
     which the caller can cast to the proper data type
   returns a NULL if Atom::extract() does not recognize the name
   customize by adding names to Atom::extract()
------------------------------------------------------------------------- */

void *lammps_extract_atom(void *ptr, char *name)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  return lmp->atom->extract(name);
}

/* ----------------------------------------------------------------------
   extract a pointer to an internal LAMMPS compute-based entity
   id = compute ID
   style = 0 for global data, 1 for per-atom data, 2 for local data
   type = 0 for scalar, 1 for vector, 2 for array
   for global data, returns a pointer to the
     compute's internal data structure for the entity
     caller should cast it to (double *) for a scalar or vector
     caller should cast it to (double **) for an array
   for per-atom or local data, returns a pointer to the
     compute's internal data structure for the entity
     caller should cast it to (double *) for a vector
     caller should cast it to (double **) for an array
   returns a void pointer to the compute's internal data structure
     for the entity which the caller can cast to the proper data type
   returns a NULL if id is not recognized or style/type not supported
   IMPORTANT: if the compute is not current it will be invoked
     LAMMPS cannot easily check here if it is valid to invoke the compute,
     so caller must insure that it is OK
------------------------------------------------------------------------- */

void *lammps_extract_compute(void *ptr, char *id, int style, int type)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  int icompute = lmp->modify->find_compute(id);
  if (icompute < 0) return NULL;
  Compute *compute = lmp->modify->compute[icompute];

  if (style == 0) {
    if (type == 0) {
      if (!compute->scalar_flag) return NULL;
      if (compute->invoked_scalar != lmp->update->ntimestep)
        compute->compute_scalar();
      return (void *) &compute->scalar;
    }
    if (type == 1) {
      if (!compute->vector_flag) return NULL;
      if (compute->invoked_vector != lmp->update->ntimestep)
        compute->compute_vector();
      return (void *) compute->vector;
    }
    if (type == 2) {
      if (!compute->array_flag) return NULL;
      if (compute->invoked_array != lmp->update->ntimestep)
        compute->compute_array();
      return (void *) compute->array;
    }
  }

  if (style == 1) {
    if (!compute->peratom_flag) return NULL;
    if (type == 1) {
      if (compute->invoked_peratom != lmp->update->ntimestep)
        compute->compute_peratom();
      return (void *) compute->vector_atom;
    }
    if (type == 2) {
      if (compute->invoked_peratom != lmp->update->ntimestep)
        compute->compute_peratom();
      return (void *) compute->array_atom;
    }
  }

  if (style == 2) {
    if (!compute->local_flag) return NULL;
    if (type == 1) {
      if (compute->invoked_local != lmp->update->ntimestep)
        compute->compute_local();
      return (void *) compute->vector_local;
    }
    if (type == 2) {
      if (compute->invoked_local != lmp->update->ntimestep)
        compute->compute_local();
      return (void *) compute->array_local;
    }
  }

  return NULL;
}

/* ----------------------------------------------------------------------
   extract a pointer to an internal LAMMPS fix-based entity
   id = fix ID
   style = 0 for global data, 1 for per-atom data, 2 for local data
   type = 0 for scalar, 1 for vector, 2 for array
   i,j = indices needed only to specify which global vector or array value
   for global data, returns a pointer to a memory location
     which is allocated by this function
     which the caller can cast to a (double *) which points to the value
   for per-atom or local data, returns a pointer to the
     fix's internal data structure for the entity
     caller should cast it to (double *) for a vector
     caller should cast it to (double **) for an array
   returns a NULL if id is not recognized or style/type not supported
   IMPORTANT: for global data,
     this function allocates a double to store the value in,
     so the caller must free this memory to avoid a leak, e.g.
       double *dptr = (double *) lammps_extract_fix();
       double value = *dptr;
       free(dptr);
   IMPORTANT: LAMMPS cannot easily check here when info extracted from
     the fix is valid, so caller must insure that it is OK
------------------------------------------------------------------------- */

void *lammps_extract_fix(void *ptr, char *id, int style, int type,
                         int i, int j)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  int ifix = lmp->modify->find_fix(id);
  if (ifix < 0) return NULL;
  Fix *fix = lmp->modify->fix[ifix];

  if (style == 0) {
    double *dptr = (double *) malloc(sizeof(double));
    if (type == 0) {
      if (!fix->scalar_flag) return NULL;
      *dptr = fix->compute_scalar();
      return (void *) dptr;
    }
    if (type == 1) {
      if (!fix->vector_flag) return NULL;
      *dptr = fix->compute_vector(i);
      return (void *) dptr;
    }
    if (type == 2) {
      if (!fix->array_flag) return NULL;
      *dptr = fix->compute_array(i,j);
      return (void *) dptr;
    }
  }

  if (style == 1) {
    if (!fix->peratom_flag) return NULL;
    if (type == 1) return (void *) fix->vector_atom;
    if (type == 2) return (void *) fix->array_atom;
  }

  if (style == 2) {
    if (!fix->local_flag) return NULL;
    if (type == 1) return (void *) fix->vector_local;
    if (type == 2) return (void *) fix->array_local;
  }

  return NULL;
}

/* ----------------------------------------------------------------------
   extract a pointer to an internal LAMMPS evaluated variable
   name = variable name, must be equal-style or atom-style variable
   group = group ID for evaluating an atom-style variable, else NULL
   for equal-style variable, returns a pointer to a memory location
     which is allocated by this function
     which the caller can cast to a (double *) which points to the value
   for atom-style variable, returns a pointer to the
     vector of per-atom values on each processor,
     which the caller can cast to a (double *) which points to the values
   returns a NULL if name is not recognized or not equal-style or atom-style
   IMPORTANT: for both equal-style and atom-style variables,
     this function allocates memory to store the variable data in
     so the caller must free this memory to avoid a leak
     e.g. for equal-style variables
       double *dptr = (double *) lammps_extract_variable();
       double value = *dptr;
       free(dptr);
     e.g. for atom-style variables
       double *vector = (double *) lammps_extract_variable();
       use the vector values
       free(vector);
   IMPORTANT: LAMMPS cannot easily check here when it is valid to evaluate
     the variable or any fixes or computes or thermodynamic info it references,
     so caller must insure that it is OK
------------------------------------------------------------------------- */

void *lammps_extract_variable(void *ptr, char *name, char *group)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  int ivar = lmp->input->variable->find(name);
  if (ivar < 0) return NULL;

  if (lmp->input->variable->equalstyle(ivar)) {
    double *dptr = (double *) malloc(sizeof(double));
    *dptr = lmp->input->variable->compute_equal(ivar);
    return (void *) dptr;
  }

  if (lmp->input->variable->atomstyle(ivar)) {
    int igroup = lmp->group->find(group);
    if (igroup < 0) return NULL;
    int nlocal = lmp->atom->nlocal;
    double *vector = (double *) malloc(nlocal*sizeof(double));
    lmp->input->variable->compute_atom(ivar,igroup,vector,1,0);
    return (void *) vector;
  }

  return NULL;
}

/* ----------------------------------------------------------------------
   return the total number of atoms in the system
   useful before call to lammps_get_atoms() so can pre-allocate vector
------------------------------------------------------------------------- */

int lammps_get_natoms(void *ptr)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  if (lmp->atom->natoms > MAXSMALLINT) return 0;
  int natoms = static_cast<int> (lmp->atom->natoms);
  return natoms;
}

/* ----------------------------------------------------------------------
   gather the named atom-based entity across all processors
   name = desired quantity, e.g. x or charge
   type = 0 for integer values, 1 for double values
   count = # of per-atom values, e.g. 1 for type or charge, 3 for x or f
   return atom-based values in data, ordered by count, then by atom ID
     e.g. x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...
   data must be pre-allocated by caller to correct length
------------------------------------------------------------------------- */

void lammps_gather_atoms(void *ptr, char *name, 
                         int type, int count, void *data)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  // error if tags are not defined or not consecutive

  int flag = 0;
  if (lmp->atom->tag_enable == 0 || lmp->atom->tag_consecutive() == 0) flag = 1;
  if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
  if (flag && lmp->comm->me == 0) {
    lmp->error->warning(FLERR,"Library error in lammps_gather_atoms");
    return;
  }

  int natoms = static_cast<int> (lmp->atom->natoms);

  int i,j,offset;
  void *vptr = lmp->atom->extract(name);

  // copy = Natom length vector of per-atom values
  // use atom ID to insert each atom's values into copy
  // MPI_Allreduce with MPI_SUM to merge into data, ordered by atom ID

  if (type == 0) {
    int *vector = NULL;
    int **array = NULL;
    if (count == 1) vector = (int *) vptr;
    else array = (int **) vptr;

    int *copy;
    lmp->memory->create(copy,count*natoms,"lib/gather:copy");
    for (i = 0; i < count*natoms; i++) copy[i] = 0;

    int *tag = lmp->atom->tag;
    int nlocal = lmp->atom->nlocal;

    if (count == 1)
      for (i = 0; i < nlocal; i++)
        copy[tag[i]-1] = vector[i];
    else
      for (i = 0; i < nlocal; i++) {
        offset = count*(tag[i]-1);
        for (j = 0; j < count; j++)
          copy[offset++] = array[i][0];
      }

    MPI_Allreduce(copy,data,count*natoms,MPI_INT,MPI_SUM,lmp->world);
    lmp->memory->destroy(copy);

  } else {
    double *vector = NULL;
    double **array = NULL;
    if (count == 1) vector = (double *) vptr;
    else array = (double **) vptr;

    double *copy;
    lmp->memory->create(copy,count*natoms,"lib/gather:copy");
    for (i = 0; i < count*natoms; i++) copy[i] = 0.0;

    int *tag = lmp->atom->tag;
    int nlocal = lmp->atom->nlocal;

    if (count == 1) {
      for (i = 0; i < nlocal; i++)
        copy[tag[i]-1] = vector[i];
    } else {
      for (i = 0; i < nlocal; i++) {
        offset = count*(tag[i]-1);
        for (j = 0; j < count; j++)
          copy[offset++] = array[i][j];
      }
    }

    MPI_Allreduce(copy,data,count*natoms,MPI_DOUBLE,MPI_SUM,lmp->world);
    lmp->memory->destroy(copy);
  }
}

/* ----------------------------------------------------------------------
   scatter the named atom-based entity across all processors
   name = desired quantity, e.g. x or charge
   type = 0 for integer values, 1 for double values
   count = # of per-atom values, e.g. 1 for type or charge, 3 for x or f
   data = atom-based values in data, ordered by count, then by atom ID
     e.g. x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...
------------------------------------------------------------------------- */

void lammps_scatter_atoms(void *ptr, char *name,
                          int type, int count, void *data)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  // error if tags are not defined or not consecutive or no atom map

  int flag = 0;
  if (lmp->atom->tag_enable == 0 || lmp->atom->tag_consecutive() == 0) flag = 1;
  if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
  if (lmp->atom->map_style == 0) flag = 1;
  if (flag && lmp->comm->me == 0) {
    lmp->error->warning(FLERR,"Library error in lammps_scatter_atoms");
    return;
  }

  int natoms = static_cast<int> (lmp->atom->natoms);

  int i,j,m,offset;
  void *vptr = lmp->atom->extract(name);

  // copy = Natom length vector of per-atom values
  // use atom ID to insert each atom's values into copy
  // MPI_Allreduce with MPI_SUM to merge into data, ordered by atom ID

  if (type == 0) {
    int *vector = NULL;
    int **array = NULL;
    if (count == 1) vector = (int *) vptr;
    else array = (int **) vptr;
    int *dptr = (int *) data;

    if (count == 1) {
      for (i = 0; i < natoms; i++)
        if ((m = lmp->atom->map(i+1)) >= 0)
          vector[m] = dptr[i];
    } else {
      for (i = 0; i < natoms; i++)
        if ((m = lmp->atom->map(i+1)) >= 0) {
          offset = count*i;
          for (j = 0; j < count; j++)
            array[m][j] = dptr[offset++];
        }
    }
  } else {
    double *vector = NULL;
    double **array = NULL;
    if (count == 1) vector = (double *) vptr;
    else array = (double **) vptr;
    double *dptr = (double *) data;

    if (count == 1) {
      for (i = 0; i < natoms; i++)
        if ((m = lmp->atom->map(i+1)) >= 0)
          vector[m] = dptr[i];
    } else {
      for (i = 0; i < natoms; i++) {
        if ((m = lmp->atom->map(i+1)) >= 0) {
          offset = count*i;
          for (j = 0; j < count; j++)
            array[m][j] = dptr[offset++];
        }
      }
    }
  }
}
