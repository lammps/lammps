/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LIBLAMMPSPLUGIN_H
#define LIBLAMMPSPLUGIN_H
/*
   Variant of the C style library interface to LAMMPS
   that uses a shared library and dynamically opens it,
   so this can be used as a prototype code to integrate
   a LAMMPS plugin to some other software.
*/

/*
 * Follow the behavior of regular LAMMPS compilation and assume
 * -DLAMMPS_SMALLBIG when no define is set.
 */
#if !defined(LAMMPS_BIGBIG) && !defined(LAMMPS_SMALLBIG) && !defined(LAMMPS_SMALLSMALL)
#define LAMMPS_SMALLBIG
#endif

#if defined(LAMMPS_LIB_MPI)
#include <mpi.h>
#endif

#if defined(LAMMPS_BIGBIG) || defined(LAMMPS_SMALLBIG)
#include <stdint.h>  /* for int64_t */
#endif

/* The following enums must be kept in sync with the equivalent enums
 * or constants in python/lammps/constants.py, fortran/lammps.f90,
 * tools/swig/lammps.i, and examples/COUPLE/plugin/liblammpsplugin.h */

/* Data type constants for extracting data from atoms, computes and fixes */

enum _LMP_DATATYPE_CONST {
  LAMMPS_INT = 0,       /*!< 32-bit integer (array) */
  LAMMPS_INT_2D = 1,    /*!< two-dimensional 32-bit integer array */
  LAMMPS_DOUBLE = 2,    /*!< 64-bit double (array) */
  LAMMPS_DOUBLE_2D = 3, /*!< two-dimensional 64-bit double array */
  LAMMPS_INT64 = 4,     /*!< 64-bit integer (array) */
  LAMMPS_INT64_2D = 5,  /*!< two-dimensional 64-bit integer array */
  LAMMPS_STRING = 6     /*!< C-String */
};

/* Style constants for extracting data from computes and fixes. */

enum _LMP_STYLE_CONST {
  LMP_STYLE_GLOBAL = 0, /*!< return global data */
  LMP_STYLE_ATOM = 1,   /*!< return per-atom data */
  LMP_STYLE_LOCAL = 2   /*!< return local data */
};

/* Type and size constants for extracting data from computes and fixes. */

enum _LMP_TYPE_CONST {
  LMP_TYPE_SCALAR = 0, /*!< return scalar */
  LMP_TYPE_VECTOR = 1, /*!< return vector */
  LMP_TYPE_ARRAY = 2,  /*!< return array */
  LMP_SIZE_VECTOR = 3, /*!< return length of vector */
  LMP_SIZE_ROWS = 4,   /*!< return number of rows */
  LMP_SIZE_COLS = 5    /*!< return number of columns */
};

/* Error codes to select the suitable function in the Error class */

enum _LMP_ERROR_CONST {
  LMP_ERROR_WARNING = 0, /*!< call Error::warning() */
  LMP_ERROR_ONE = 1,     /*!< called from one MPI rank */
  LMP_ERROR_ALL = 2,     /*!< called from all MPI ranks */
  LMP_ERROR_WORLD = 4,   /*!< error on Comm::world */
  LMP_ERROR_UNIVERSE = 8 /*!< error on Comm::universe */
};

/** Variable style constants for extracting data from variables.
 *
 * Must be kept in sync with the equivalent constants in python/lammps/constants.py,
 * fortran/lammps.f90, and tools/swig/lammps.i */

enum _LMP_VAR_CONST {
  LMP_VAR_EQUAL = 0,  /*!< compatible with equal-style variables */
  LMP_VAR_ATOM = 1,   /*!< compatible with atom-style variables */
  LMP_VAR_VECTOR = 2, /*!< compatible with vector-style variables */
  LMP_VAR_STRING = 3  /*!< return value will be a string (catch-all) */
};

#ifdef __cplusplus
extern "C" {
#endif

#if defined(LAMMPS_BIGBIG)
typedef void (*FixExternalFnPtr)(void *, int64_t, int, int64_t *, double **, double **);
#elif defined(LAMMPS_SMALLSMALL)
typedef void (*FixExternalFnPtr)(void *, int, int, int *, double **, double **);
#else
typedef void (*FixExternalFnPtr)(void *, int64_t, int, int *, double **, double **);
#endif

#define LAMMPSPLUGIN_ABI_VERSION 1
struct _liblammpsplugin {
  int abiversion;
  int has_exceptions;
  void *handle;
#if defined(LAMMPS_LIB_MPI)
  void *(*open)(int, char **, MPI_Comm, void **);
#else
  void *open;
#endif
  void *(*open_no_mpi)(int, char **, void **);
  void *(*open_fortran)(int, char **, void **, int);
  void (*close)(void *);

  void (*mpi_init)();
  void (*mpi_finalize)();
  void (*kokkos_finalize)();
  void (*python_finalize)();

  void (*error)(void *, int, const char *);

  void (*file)(void *, char *);
  char *(*command)(void *, const char *);
  void (*commands_list)(void *, int, const char **);
  void (*commands_string)(void *, const char *);

  double (*get_natoms)(void *);
  double (*get_thermo)(void *, const char *);
  void *(*last_thermo)(void *, const char *, int);

  void (*extract_box)(void *, double *, double *,
                      double *, double *, double *, int *, int *);
  void (*reset_box)(void *, double *, double *, double, double, double);

  void (*memory_usage)(void *, double *);
  int (*get_mpi_comm)(void *);

  int (*extract_setting)(void *, const char *);
  int *(*extract_global_datatype)(void *, const char *);
  void *(*extract_global)(void *, const char *);

  int *(*extract_atom_datatype)(void *, const char *);
  void *(*extract_atom)(void *, const char *);

  void *(*extract_compute)(void *, const char *, int, int);
  void *(*extract_fix)(void *, const char *, int, int, int, int);
  void *(*extract_variable)(void *, const char *, char *);
  int (*extract_variable_datatype)(void *, const char *);
  int (*set_variable)(void *, char *, char *);

  void (*gather_atoms)(void *, const char *, int, int, void *);
  void (*gather_atoms_concat)(void *, const char *, int, int, void *);
  void (*gather_atoms_subset)(void *, const char *, int, int, int, int *, void *);
  void (*scatter_atoms)(void *, const char *, int, int, void *);
  void (*scatter_atoms_subset)(void *, const char *, int, int, int, int *, void *);

  void (*gather_bonds)(void *, void *);
  void (*gather_angles)(void *, void *);
  void (*gather_dihedrals)(void *, void *);
  void (*gather_impropers)(void *, void *);

  void (*gather)(void *, const char *, int, int, void *);
  void (*gather_concat)(void *, const char *, int, int, void *);
  void (*gather_subset)(void *, const char *, int, int, int, int *,void *);
  void (*scatter)(void *, const char *, int, int, void *);
  void (*scatter_subset)(void *, const char *, int, int, int, int *, void *);

/* lammps_create_atoms() takes tagint and imageint as args
 * the ifdef ensures they are compatible with rest of LAMMPS
 * caller must match to how LAMMPS library is built */

#ifndef LAMMPS_BIGBIG
 int (*create_atoms)(void *, int, int *, int *, double *, double *, int *, int);
#else
  int (*create_atoms)(void *, int, int64_t *, int *, double *, double *, int64_t *, int);
#endif

  int (*find_pair_neighlist)(void *, const char *, int, int, int);
  int (*find_fix_neighlist)(void *, const char *, int);
  int (*find_compute_neighlist)(void *, const char *, int);
  int (*neighlist_num_elements)(void *, int);
  void (*neighlist_element_neighbors)(void *, int, int, int *, int *, int **);

  int (*version)(void *);
  void (*get_os_info)(char *, int);

  int (*config_has_mpi_support)();
  int (*config_has_gzip_support)();
  int (*config_has_png_support)();
  int (*config_has_jpeg_support)();
  int (*config_has_ffmpeg_support)();
  int (*config_has_exceptions)();

  int (*config_has_package)(const char *);
  int (*config_package_count)();
  int (*config_package_name)(int, char *, int);

  int (*config_accelerator)(const char *, const char *, const char *);
  int (*has_gpu_device)();
  void (*get_gpu_device_info)(char *, int);

  int (*has_style)(void *, const char *, const char *);
  int (*style_count)(void *, const char *);
  int (*style_name)(void *, const char *, int, char *, int);

  int (*has_id)(void *, const char *, const char *);
  int (*id_count)(void *, const char *);
  int (*id_name)(void *, const char *, int, char *, int);

  int (*plugin_count)();
  int (*plugin_name)(int, char *, char *, int);

#if !defined(LAMMPS_BIGBIG)
  int (*encode_image_flags)(int, int, int);
  void (*decode_image_flags)(int, int *);
#else
  int64_t (*encode_image_flags)(int, int, int);
  void (*decode_image_flags)(int64_t, int *);
#endif

  void (*set_fix_external_callback)(void *, const char *, FixExternalFnPtr, void *);
  double **(*fix_external_get_force)(void *, const char *);
  void (*fix_external_set_energy_global)(void *, const char *, double);
  void (*fix_external_set_energy_peratom)(void *, const char *, double *);
  void (*fix_external_set_virial_global)(void *, const char *, double *);
  void (*fix_external_set_virial_peratom)(void *, const char *, double **);
  void (*fix_external_set_vector_length)(void *, const char *, int);
  void (*fix_external_set_vector)(void *, const char *, int, double);

  void (*flush_buffers)(void *);

  void (*free)(void *);

  void (*is_running)(void *);
  void (*force_timeout)(void *);

  int (*has_error)(void *);
  int (*get_last_error_message)(void *, char *, int);

  int (*python_api_version)();
};

typedef struct _liblammpsplugin liblammpsplugin_t;

liblammpsplugin_t *liblammpsplugin_load(const char *);
int liblammpsplugin_release(liblammpsplugin_t *);

#undef LAMMPS
#ifdef __cplusplus
}
#endif

#endif
