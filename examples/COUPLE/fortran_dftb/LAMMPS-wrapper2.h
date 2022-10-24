/* -----------------------------------------------------------------------
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    https://www.lammps.org/, Sandia National Laboratories
    LAMMPS development team: developers@lammps.org
 
    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under 
    the GNU General Public License.
 
    See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
    Contributing author:  Nir Goldman, LLNL <ngoldman@llnl.gov>, 2016
------------------------------------------------------------------------- */

/* This is set of "wrapper" functions to assist LAMMPS.F90, which itself
   provides a (I hope) robust Fortran interface to library.cpp and
   library.h.  All prototypes herein COULD be added to library.h instead of
   including this as a separate file. See the README for instructions. */
#ifdef __cplusplus
extern "C" {
#endif

/* Prototypes for auxiliary functions */
void lammps_set_callback (void *); 
void lammps_set_user_energy (void*, double); 
void lammps_set_user_virial (void*, double*); 
void lammps_set_external_vector_length (void*, int); 
void lammps_set_external_vector (void*, int, double); 

#ifdef __cplusplus
}
#endif

/* vim: set ts=3 sts=3 expandtab: */
