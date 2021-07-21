/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* libwrapper = fortran wrappers for LAMMPS library functions.
   See README for compilation instructions */

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "stdint.h"
#include "library.h"        /* this is a LAMMPS include file */

/* wrapper for creating a lammps instance from fortran.
   since fortran has no simple way to emit a C-compatible
   argument array, we don't support it. for simplicity,
   the address of the pointer to the lammps object is
   stored in a 64-bit integer on all platforms. */

void lammps_open_(MPI_Fint *comm, int64_t *ptr)
{
    *ptr = (int64_t) lammps_open_fortran(0,NULL,*comm);
}

/* no-MPI version of the wrapper from above. */

void lammps_open_no_mpi_(int64_t *ptr)
{
    void *obj;

    lammps_open_no_mpi(0,NULL,&obj);
    *ptr = (int64_t) obj;
}

/* wrapper for shutting down a lammps instance from fortran. */

void lammps_close_(int64_t *ptr)
{
    void *obj;
    obj = (void *) *ptr;

    lammps_close(obj);
}

/* wrapper for passing an input file to lammps from fortran.
   since fortran strings are not zero terminated, we have
   to pass the length explicitly and make a copy that is. */

void lammps_file_(int64_t *ptr, char *fname, MPI_Fint *len)
{
    void *obj;
    char *cpy;

    obj = (void *) *ptr;

    cpy = (char *)calloc(*len + 1,sizeof(char));
    memcpy(cpy,fname,*len);

    lammps_file(obj,cpy);
    free(cpy);
}

/* wrapper for passing a line input to lammps from fortran.
   since fortran strings are not zero terminated, we have
   to pass the length explicitly and make a copy that is. */

void lammps_command_(int64_t *ptr, char *line, MPI_Fint *len)
{
    void *obj;
    char *cpy;

    obj = (void *) *ptr;

    cpy = (char *)calloc(*len + 1,sizeof(char));
    memcpy(cpy,line,*len);

    lammps_command(obj,cpy);
    free(cpy);
}

/* fortran wrapper to get the number of atoms from lammps.
   return values require an interface in fortran, so we 
   make the wrapper into a procedure. */

void lammps_get_natoms_(int64_t *ptr, MPI_Fint *natoms)
{
    void *obj;
    obj = (void *) *ptr;

    *natoms = lammps_get_natoms(obj);
}

