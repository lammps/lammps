/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*
   C or Fortran style library interface to LAMMPS
   new LAMMPS-specific functions can be added
*/

/*
 * Follow the behavior of regular LAMMPS compilation and assume
 * -DLAMMPS_SMALLBIG when no define is set.
 */
#if !defined(LAMMPS_BIGBIG) && !defined(LAMMPS_SMALLBIG) && !defined(LAMMPS_SMALLSMALL)
#define LAMMPS_SMALLBIG
#endif

#include <mpi.h>
#if defined(LAMMPS_BIGBIG) || defined(LAMMPS_SMALLBIG)
#include <inttypes.h>  /* for int64_t */
#endif

/* ifdefs allow this file to be included in a C program */

#ifdef __cplusplus
extern "C" {
#endif

void lammps_open(int, char **, MPI_Comm, void **);
void lammps_open_no_mpi(int, char **, void **);
void lammps_close(void *);
int  lammps_version(void *);
void lammps_file(void *, char *);
char *lammps_command(void *, char *);
void lammps_commands_list(void *, int, char **);
void lammps_commands_string(void *, char *);
void lammps_free(void *);

int lammps_extract_setting(void *, char *);
void *lammps_extract_global(void *, char *);
void lammps_extract_box(void *, double *, double *,
                        double *, double *, double *, int *, int *);
void *lammps_extract_atom(void *, char *);
void *lammps_extract_compute(void *, char *, int, int);
void *lammps_extract_fix(void *, char *, int, int, int, int);
void *lammps_extract_variable(void *, char *, char *);

double lammps_get_thermo(void *, char *);
int lammps_get_natoms(void *);

int lammps_set_variable(void *, char *, char *);
void lammps_reset_box(void *, double *, double *, double, double, double);

void lammps_gather_atoms(void *, char *, int, int, void *);
void lammps_gather_atoms_concat(void *, char *, int, int, void *);
void lammps_gather_atoms_subset(void *, char *, int, int, int, int *, void *);
void lammps_scatter_atoms(void *, char *, int, int, void *);
void lammps_scatter_atoms_subset(void *, char *, int, int, int, int *, void *);

#if defined(LAMMPS_BIGBIG)
typedef void (*FixExternalFnPtr)(void *, int64_t, int, int64_t *, double **, double **);
void lammps_set_fix_external_callback(void *, char *, FixExternalFnPtr, void*);
#elif defined(LAMMPS_SMALLBIG)
typedef void (*FixExternalFnPtr)(void *, int64_t, int, int *, double **, double **);
void lammps_set_fix_external_callback(void *, char *, FixExternalFnPtr, void*);
#else
typedef void (*FixExternalFnPtr)(void *, int, int, int *, double **, double **);
void lammps_set_fix_external_callback(void *, char *, FixExternalFnPtr, void*);
#endif

int lammps_config_has_package(char * package_name);
int lammps_config_package_count();
int lammps_config_package_name(int index, char * buffer, int max_size);
int lammps_config_has_gzip_support();
int lammps_config_has_png_support();
int lammps_config_has_jpeg_support();
int lammps_config_has_ffmpeg_support();
int lammps_config_has_exceptions();

// lammps_create_atoms() takes tagint and imageint as args
// ifdef insures they are compatible with rest of LAMMPS
// caller must match to how LAMMPS library is built

#ifdef LAMMPS_BIGBIG
void lammps_create_atoms(void *, int, int64_t *, int *,
                         double *, double *, int64_t *, int);
#else
void lammps_create_atoms(void *, int, int *, int *,
                         double *, double *, int *, int);
#endif

#ifdef LAMMPS_EXCEPTIONS
int lammps_has_error(void *);
int lammps_get_last_error_message(void *, char *, int);
#endif

#undef LAMMPS
#ifdef __cplusplus
}
#endif

/* ERROR/WARNING messages:

E: Library error: issuing LAMMPS command during run

UNDOCUMENTED

W: Library error in lammps_gather_atoms

This library function cannot be used if atom IDs are not defined
or are not consecutively numbered.

W: lammps_gather_atoms: unknown property name

UNDOCUMENTED

W: Library error in lammps_gather_atoms_subset

UNDOCUMENTED

W: lammps_gather_atoms_subset: unknown property name

UNDOCUMENTED

W: Library error in lammps_scatter_atoms

This library function cannot be used if atom IDs are not defined or
are not consecutively numbered, or if no atom map is defined.  See the
atom_modify command for details about atom maps.

W: lammps_scatter_atoms: unknown property name

UNDOCUMENTED

W: Library error in lammps_scatter_atoms_subset

UNDOCUMENTED

W: lammps_scatter_atoms_subset: unknown property name

UNDOCUMENTED

W: Library error in lammps_create_atoms

UNDOCUMENTED

W: Library warning in lammps_create_atoms, invalid total atoms %ld %ld

UNDOCUMENTED

*/
