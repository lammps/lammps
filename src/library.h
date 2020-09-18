/* -*- c -*- ------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LAMMPS_LIBRARY_H
#define LAMMPS_LIBRARY_H

/* C style library interface to LAMMPS which allows to create and
 * control instances of the LAMMPS C++ class and exchange data with it.
 * The C bindings are the basis for the Python and Fortran modules.
 *
 * If needed, new LAMMPS-specific functions can be added to expose
 * additional LAMMPS functionality to this library interface. */

/* We follow the behavior of regular LAMMPS compilation and assume
 * -DLAMMPS_SMALLBIG when no define is set. */

#if  !defined(LAMMPS_BIGBIG)     \
  && !defined(LAMMPS_SMALLBIG)   \
  && !defined(LAMMPS_SMALLSMALL)
#define LAMMPS_SMALLBIG
#endif

/* To allow including the library interface without MPI */

#if !defined(LAMMPS_LIB_NO_MPI)
#include <mpi.h>
#endif

#if defined(LAMMPS_BIGBIG) || defined(LAMMPS_SMALLBIG)
#include <inttypes.h>  /* for int64_t */
#endif

/** Data type constants for extracting data from atoms, computes and fixes
 *
 * Must be kept in sync with the equivalent constants in lammps.py */

enum _LMP_DATATYPE_CONST {
  LAMMPS_INT    = 0,       /*!< 32-bit integer (array) */
  LAMMPS_INT_2D  = 1,      /*!< two-dimensional 32-bit integer array */
  LAMMPS_DOUBLE = 2,       /*!< 64-bit double (array) */
  LAMMPS_DOUBLE_2D = 3,    /*!< two-dimensional 64-bit double array */
  LAMMPS_INT64 = 4,        /*!< 64-bit integer (array) */
  LAMMPS_INT64_2D = 5,     /*!< two-dimensional 64-bit integer array */
  LAMMPS_STRING = 6        /*!< C-String */
};

/** Style constants for extracting data from computes and fixes.
 *
 * Must be kept in sync with the equivalent constants in lammps.py */

enum _LMP_STYLE_CONST {
  LMP_STYLE_GLOBAL=0,           /*!< return global data */
  LMP_STYLE_ATOM  =1,           /*!< return per-atom data */
  LMP_STYLE_LOCAL =2            /*!< return local data */
};

/** Type and size constants for extracting data from computes and fixes.
 *
 * Must be kept in sync with the equivalent constants in lammps.py */

enum _LMP_TYPE_CONST {
  LMP_TYPE_SCALAR=0,            /*!< return scalar */
  LMP_TYPE_VECTOR=1,            /*!< return vector */
  LMP_TYPE_ARRAY =2,            /*!< return array */
  LMP_SIZE_VECTOR=3,            /*!< return length of vector */
  LMP_SIZE_ROWS  =4,            /*!< return number of rows */
  LMP_SIZE_COLS  =5             /*!< return number of columns */
};

/* Ifdefs to allow this file to be included in C and C++ programs */

#ifdef __cplusplus
extern "C" {
#endif

/* ----------------------------------------------------------------------
 * Library functions to create/destroy an instance of LAMMPS
 * ---------------------------------------------------------------------- */

#if !defined(LAMMPS_LIB_NO_MPI)
void *lammps_open(int argc, char **argv, MPI_Comm comm, void **ptr);
#endif
void *lammps_open_no_mpi(int argc, char **argv, void **ptr);
void *lammps_open_fortran(int argc, char **argv, int f_comm);
void  lammps_close(void *handle);
void  lammps_mpi_init();
void  lammps_mpi_finalize();
void  lammps_free(void *ptr);

/* ----------------------------------------------------------------------
 * Library functions to process commands
 * ---------------------------------------------------------------------- */

void  lammps_file(void *handle, const char *file);

char *lammps_command(void *handle, const char *cmd);
void  lammps_commands_list(void *handle, int ncmd, const char **cmds);
void  lammps_commands_string(void *handle, const char *str);

/* -----------------------------------------------------------------------
 * Library functions to extract info from LAMMPS or set data in LAMMPS
 * ----------------------------------------------------------------------- */

int    lammps_version(void *handle);
void   lammps_memory_usage(void *handle, double *meminfo);
int    lammps_get_mpi_comm(void *handle);
double lammps_get_natoms(void *handle);
double lammps_get_thermo(void *handle, const char *keyword);

void   lammps_extract_box(void *handle, double *boxlo, double *boxhi,
                          double *xy, double *yz, double *xz,
                          int *pflags, int *boxflag);
void   lammps_reset_box(void *handle, double *boxlo, double *boxhi,
                        double xy, double yz, double xz);

int    lammps_extract_setting(void *handle, const char *keyword);
void  *lammps_extract_global(void *handle, const char *name);
void  *lammps_extract_atom(void *handle, const char *name);

int lammps_extract_global_datatype(void *handle, const char *name);
int lammps_extract_atom_datatype(void *handle, const char *name);

#if !defined(LAMMPS_BIGBIG)
int    lammps_create_atoms(void *handle, int n, int *id, int *type,
                           double *x, double *v, int *image, int bexpand);
#else
int    lammps_create_atoms(void *handle, int n, int64_t *id, int *type,
                           double *x, double *v, int64_t* image, int bexpand);
#endif

/* ----------------------------------------------------------------------
 * Library functions to access data from computes, fixes, variables in LAMMPS
 * ---------------------------------------------------------------------- */

void *lammps_extract_compute(void *handle, char *id, int, int);
void *lammps_extract_fix(void *handle, char *, int, int, int, int);
void *lammps_extract_variable(void *handle, char *, char *);
int   lammps_set_variable(void *, char *, char *);

/* ----------------------------------------------------------------------
 * Library functions for scatter/gather operations of data
 * ---------------------------------------------------------------------- */


void lammps_gather(void *, char *, int, int, void *);
void lammps_gather_concat(void *, char *, int, int, void *);
void lammps_gather_subset(void *, char *, int, int, int, int *, void *);
void lammps_scatter(void *, char *, int, int, void *);
void lammps_scatter_subset(void *, char *, int, int, int, int *, void *);


void lammps_gather_atoms(void *, char *, int, int, void *);
void lammps_gather_atoms_concat(void *, char *, int, int, void *);
void lammps_gather_atoms_subset(void *, char *, int, int, int, int *, void *);
void lammps_scatter_atoms(void *, char *, int, int, void *);
void lammps_scatter_atoms_subset(void *, char *, int, int, int, int *, void *);

/* ----------------------------------------------------------------------
 * Library functions for retrieving configuration information
 * ---------------------------------------------------------------------- */

int lammps_config_has_mpi_support();
int lammps_config_has_package(const char *);
int lammps_config_package_count();
int lammps_config_package_name(int, char *, int);
int lammps_config_has_gzip_support();
int lammps_config_has_png_support();
int lammps_config_has_jpeg_support();
int lammps_config_has_ffmpeg_support();
int lammps_config_has_exceptions();

int lammps_has_style(void *, const char *, const char *);
int lammps_style_count(void *, const char *);
int lammps_style_name(void *, const char *, int, char *, int);

/* ----------------------------------------------------------------------
 * Library functions for accessing neighbor lists
 * ---------------------------------------------------------------------- */

int lammps_find_pair_neighlist(void*, char *, int, int, int);
int lammps_find_fix_neighlist(void*, char *, int);
int lammps_find_compute_neighlist(void*, char *, int);
int lammps_neighlist_num_elements(void*, int);
void lammps_neighlist_element_neighbors(void *, int, int, int *, int *, int ** );

/* ----------------------------------------------------------------------
 * Utility functions
 * ---------------------------------------------------------------------- */

#if !defined(LAMMPS_BIGBIG)
int lammps_encode_image_flags(int ix, int iy, int iz);
void lammps_decode_image_flags(int image, int *flags);
#else
int64_t lammps_encode_image_flags(int ix, int iy, int iz);
void lammps_decode_image_flags(int64_t image, int *flags);
#endif

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
void lammps_fix_external_set_energy_global(void *, char *, double);
void lammps_fix_external_set_virial_global(void *, char *, double *);

int lammps_has_error(void *handle);
int lammps_get_last_error_message(void *handle, char *buffer, int buf_size);

#ifdef __cplusplus
}
#endif

#endif /* LAMMPS_LIBRARY_H */

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

/* Local Variables:
 * fill-column: 72
 * End: */
