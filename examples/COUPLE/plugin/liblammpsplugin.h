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
#elif defined(LAMMPS_SMALLBIG)
typedef void (*FixExternalFnPtr)(void *, int64_t, int, int *, double **, double **);
#else
typedef void (*FixExternalFnPtr)(void *, int, int, int *, double **, double **);
#endif

  
struct _liblammpsplugin {
  int abiversion;
  int has_exceptions;
  void *handle;
  void (*open)(int, char **, MPI_Comm, void **);
  void (*open_no_mpi)(int, char **, void **);
  void (*close)(void *);
  int  (*version)(void *);
  void (*file)(void *, char *);
  char *(*command)(void *, char *);
  void (*commands_list)(void *, int, char **);
  void (*commands_string)(void *, char *);
  void (*free)(void *);
  int (*extract_setting)(void *, char *);
  void *(*extract_global)(void *, char *);
  void (*extract_box)(void *, double *, double *,
		      double *, double *, double *, int *, int *);
  void *(*extract_atom)(void *, char *);
  void *(*extract_compute)(void *, char *, int, int);
  void *(*extract_fix)(void *, char *, int, int, int, int);
  void *(*extract_variable)(void *, char *, char *);

  double (*get_thermo)(void *, char *);
  int (*get_natoms)(void *);

  int (*set_variable)(void *, char *, char *);
  void (*reset_box)(void *, double *, double *, double, double, double);

  void (*gather_atoms)(void *, char *, int, int, void *);
  void (*gather_atoms_concat)(void *, char *, int, int, void *);
  void (*gather_atoms_subset)(void *, char *, int, int, int, int *, void *);
  void (*scatter_atoms)(void *, char *, int, int, void *);
  void (*scatter_atoms_subset)(void *, char *, int, int, int, int *, void *);

  void (*set_fix_external_callback)(void *, char *, FixExternalFnPtr, void*);

  int (*config_has_package)(char * package_name);
  int (*config_package_count)();
  int (*config_package_name)(int index, char * buffer, int max_size);
  int (*config_has_gzip_support)();
  int (*config_has_png_support)();
  int (*config_has_jpeg_support)();
  int (*config_has_ffmpeg_support)();
  int (*config_has_exceptions)();

  int (*find_pair_neighlist)(void* ptr, char * style, int exact, int nsub, int request);
  int (*find_fix_neighlist)(void* ptr, char * id, int request);
  int (*find_compute_neighlist)(void* ptr, char * id, int request);
  int (*neighlist_num_elements)(void* ptr, int idx);
  void (*neighlist_element_neighbors)(void * ptr, int idx, int element, int * iatom, int * numneigh, int ** neighbors);

// lammps_create_atoms() takes tagint and imageint as args
// ifdef insures they are compatible with rest of LAMMPS
// caller must match to how LAMMPS library is built

#ifdef LAMMPS_BIGBIG
  void (*create_atoms)(void *, int, int64_t *, int *,
                         double *, double *, int64_t *, int);
#else
  void (*create_atoms)(void *, int, int *, int *,
                         double *, double *, int *, int);
#endif

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
