/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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

#include <mpi.h>
#if defined(LAMMPS_BIGBIG) || defined(LAMMPS_SMALLBIG)
#include <inttypes.h>  /* for int64_t */
#endif

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

struct _liblammpsplugin {
  int abiversion;
  int has_exceptions;
  void *handle;
  void *(*open)(int, char **, MPI_Comm, void **);
  void *(*open_no_mpi)(int, char **, void **);
  void *(*open_fortran)(int, char **, void **, int);
  void (*close)(void *);

  void (*mpi_init)();
  void (*mpi_finalize)();
  void (*kokkos_finalize)();
  void (*python_finalize)();

  void (*file)(void *, char *);
  char *(*command)(void *, const char *);
  void (*commands_list)(void *, int, const char **);
  void (*commands_string)(void *, const char *);

  double (*get_natoms)(void *);
  double (*get_thermo)(void *, char *);

  void (*extract_box)(void *, double *, double *,
                      double *, double *, double *, int *, int *);
  void (*reset_box)(void *, double *, double *, double, double, double);

  void (*memory_usage)(void *, double *);
  int (*get_mpi_comm)(void *);

  int (*extract_setting)(void *, const char *);
  int *(*extract_global_datatype)(void *, const char *);
  void *(*extract_global)(void *, const char *);

  void *(*extract_atom_datatype)(void *, const char *);
  void *(*extract_atom)(void *, const char *);

  void *(*extract_compute)(void *, const char *, int, int);
  void *(*extract_fix)(void *, const char *, int, int, int, int);
  void *(*extract_variable)(void *, const char *, char *);
  int (*set_variable)(void *, char *, char *);

  void (*gather_atoms)(void *, char *, int, int, void *);
  void (*gather_atoms_concat)(void *, char *, int, int, void *);
  void (*gather_atoms_subset)(void *, char *, int, int, int, int *, void *);
  void (*scatter_atoms)(void *, char *, int, int, void *);
  void (*scatter_atoms_subset)(void *, char *, int, int, int, int *, void *);

  void (*gather_bonds)(void *, void *);
  
// lammps_create_atoms() takes tagint and imageint as args
// ifdef insures they are compatible with rest of LAMMPS
// caller must match to how LAMMPS library is built

#ifndef LAMMPS_BIGBIG
 void (*create_atoms)(void *, int, int *, int *, double *,
                      double *, int *, int);
#else
  void (*create_atoms)(void *, int, int64_t *, int *, double *,
                       double *, int64_t *, int);
#endif

  int (*find_pair_neighlist)(void *, const char *, int, int, int);
  int (*find_fix_neighlist)(void *, const char *, int);
  int (*find_compute_neighlist)(void *, char *, int);
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

  void (*set_fix_external_callback)(void *, const char *, FixExternalFnPtr, void*);
  void (*fix_external_get_force)(void *, const char *);
  void (*fix_external_set_energy_global)(void *, const char *, double);
  void (*fix_external_set_energy_peratom)(void *, const char *, double *);
  void (*fix_external_set_virial_global)(void *, const char *, double *);
  void (*fix_external_set_virial_peratom)(void *, const char *, double **);
  void (*fix_external_set_vector_length)(void *, const char *, int);
  void (*fix_external_set_vector)(void *, const char *, int, double);

  void (*free)(void *);

  void (*is_running)(void *);
  void (*force_timeout)(void *);

  int (*has_error)(void *);
  int (*get_last_error_message)(void *, char *, int);
};

typedef struct _liblammpsplugin liblammpsplugin_t;

liblammpsplugin_t *liblammpsplugin_load(const char *);
int liblammpsplugin_release(liblammpsplugin_t *);

#undef LAMMPS
#ifdef __cplusplus
}
#endif

#endif
