/* -----------------------------------------------------------------------
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    www.cs.sandia.gov/~sjplimp/lammps.html
    Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
 
    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under 
    the GNU General Public License.
 
    See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
    Contributing author:  Karl D. Hammond <hammondkd@missouri.edu>
                          University of Missouri (USA), 2012
------------------------------------------------------------------------- */

/* This is set of "wrapper" functions to assist LAMMPS.F90, which itself
   provides a (I hope) robust Fortran interface to library.cpp and
   library.h.  All prototypes herein COULD be added to library.h instead of
   including this as a separate file. See the README for instructions. */
#ifdef __cplusplus
extern "C" {
#endif

/* Prototypes for auxiliary functions */
void lammps_open_fortran_wrapper (int, char**, MPI_Fint, void**);
int lammps_get_ntypes (void*);
int lammps_extract_compute_vectorsize (void*, char*, int);
void lammps_extract_compute_arraysize (void*, char*, int, int*, int*);
int lammps_extract_fix_vectorsize (void*, char*, int);
void lammps_extract_fix_arraysize (void*, char*, int, int*, int*);
void lammps_error_all (void*, const char*, int, const char*);

#ifdef __cplusplus
}
#endif

/* vim: set ts=3 sts=3 expandtab: */
