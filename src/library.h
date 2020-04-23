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

#ifndef LAMMPS_LIBRARY_H
#define LAMMPS_LIBRARY_H

/*
 * C style library interface to LAMMPS which allows to create
 * and control instances of the LAMMPS C++ class and exchange
 * data with it.  The C bindings are then used as basis for
 * the ctypes based Python wrapper in lammps.py and also there
 * are Fortran interfaces that expose the functionality by
 * using the ISO_C_BINDINGS module in modern Fortran.
 * If needed, new LAMMPS-specific functions can be added to
 * expose additional LAMMPS functionality to this library interface.
 */

/*
 * Follow the behavior of regular LAMMPS compilation and assume
 * -DLAMMPS_SMALLBIG when no define is set.
 */
#if !defined(LAMMPS_BIGBIG) && !defined(LAMMPS_SMALLBIG) && !defined(LAMMPS_SMALLSMALL)
#define LAMMPS_SMALLBIG
#endif

#if !defined(LAMMPS_LIB_NO_MPI)
#include <mpi.h>
#endif

#if defined(LAMMPS_BIGBIG) || defined(LAMMPS_SMALLBIG)
#include <inttypes.h>  /* for int64_t */
#endif

/*! \brief Style constants for extracting data from computes
 * and fixes.
 *
 * Must be kept in sync with constants in lammps.py
 */
enum _LMP_STYLE_CONST {
  LMP_STYLE_GLOBAL=0,           /*!< return global data */
  LMP_STYLE_ATOM  =1,           /*!< return per-atom data */
  LMP_STYLE_LOCAL =2            /*!< return local data */
};

/*! \brief Type and size constants for extracting data from
 * computes and fixes.
 *
 * Must be kept in sync with constants in lammps.py
 */
enum _LMP_TYPE_CONST {
  LMP_TYPE_SCALAR=0,            /*!< return scalar */
  LMP_TYPE_VECTOR=1,            /*!< return vector */
  LMP_TYPE_ARRAY =2,            /*!< return array */
  LMP_SIZE_VECTOR=3,            /*!< return length of vector */
  LMP_SIZE_ROWS  =4,            /*!< return number of rows */
  LMP_SIZE_COLS  =5             /*!< return number of columns */
};

/* ifdefs allow this file to be included in a C program */

#ifdef __cplusplus
extern "C" {
#endif

#if !defined(LAMMPS_LIB_NO_MPI)
void lammps_open(int, char **, MPI_Comm, void **);
#endif
void lammps_open_no_mpi(int, char **, void **);
void lammps_close(void *);
void lammps_finalize();
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

int lammps_config_has_package(char *);
int lammps_config_package_count();
int lammps_config_package_name(int, char *, int);
int lammps_config_has_gzip_support();
int lammps_config_has_png_support();
int lammps_config_has_jpeg_support();
int lammps_config_has_ffmpeg_support();
int lammps_config_has_exceptions();

int lammps_has_style(void *, char *, char *);
int lammps_style_count(void *, char *);
int lammps_style_name(void *, char *, int, char *, int);

int lammps_find_pair_neighlist(void*, char *, int, int, int);
int lammps_find_fix_neighlist(void*, char *, int);
int lammps_find_compute_neighlist(void*, char *, int);
int lammps_neighlist_num_elements(void*, int);
void lammps_neighlist_element_neighbors(void *, int, int, int *, int *, int ** );

/*
 * lammps_create_atoms() takes tagint and imageint as args
 * ifdef insures they are compatible with rest of LAMMPS
 * caller must match to how LAMMPS library is built
 */

#if !defined(LAMMPS_BIGBIG)
/** Create N atoms from list of coordinates
 *
\verbatim embed:rst
This function creates additional atoms from a given list of coordinates
and a list of atom types.  Additionally the atom-IDs, velocities, and
image flags may be provided.  If atom-IDs are not provided, they will
be automatically created as a sequence following the largest existing
atom-ID.

This function is useful to add atoms to a simulation or - in tandem
with :cpp:func:`lammps_reset_box` - to restore a previously extracted
and saved state of a simulation.  Additional properties for the new
atoms can then be assigned via the :cpp:func:`lammps_scatter_atoms`
:cpp:func:`lammps_extract_atom` functions.

For non-periodic boundaries, atoms will **not** be created that have
coordinates outside the box unless it is a shrink-wrap boundary and
the shrinkexceed flag has been set to a non-zero value.  For periodic
boundaries atoms will be wrapped back into the simulation cell and
its image flags adjusted accordingly, unless explicit image flags are
provided.

The function returns the number of atoms created or -1 on failure, e.g.
when called before as box has been created.

Coordinates and velocities have to be given in a 1d-array in the order
X(1),Y(1),Z(1),X(2),Y(2),Z(2),...,X(N),Y(N),Z(N).
\endverbatim
 *
 * \param handle pointer to a previously created LAMMPS instance cast to ``void *``.
 * \param n number of atoms, N, to be added to the system
 * \param id pointer to N atom IDs; NULL will generate IDs 1 to N
 * \param type pointer to N atom types (required)
 * \param x pointer to 3N doubles with x-,y-,z- positions of new atoms (required)
 * \param v pointer to 3N doubles with x-,y-,z- velocities of new atoms (set to 0.0 if NULL)
 * \param image pointer to N imageint sets of image flags, or NULL
 * \param shrinkexceed if 1 atoms outside shrink-wrap boundaries will still be created and not dropped.
 * \return number of atoms created on success, -1 on failure (no box, no atom IDs, etc.)
 */
int lammps_create_atoms(void *handle, int n, int *id, int *type,
                        double *x, double *v, int *image, int shrinkexceed);
#else
/** \brief Create N atoms from list of coordinates
\verbatim embed:rst
This is the interface of the :cpp:func:`lammps_create_atoms`
function if LAMMPS has been compiled with the -DLAMMPS_BIGBIG setting.
\endverbatim
 */
int lammps_create_atoms(void *, int, int64_t *, int *,
                        double *, double *, int64_t *, int);
#endif

#if !defined(LAMMPS_BIGBIG)
/** \brief Encode three integer image flags into a single imageint
 *
 * This function performs the bit-shift, addition, and bit-wise OR
 * operations necessary to combine the values of three integers
 * representing the image flags in x-, y-, and z-direction.  Unless
 * LAMMPS is compiled with -DLAMMPS_BIGBIG, those integers are
 * limited 10-bit signed integers [-512, 511].  Otherwise the return
 * type changes from ``int`` to ``int64_t`` and the valid range for
 * the individual image flags becomes [-1048576,1048575],
 * i.e. that of a 21-bit signed integer.  There is no check on whether
 * the arguments conform to these requirements.
 *
 * \param ix  image flag value in x
 * \param iy  image flag value in y
 * \param iz  image flag value in z
 * \return encoded image flags
 */
int lammps_encode_image_flags(int ix, int iy, int iz);

/** \brief Decode a single image flag integer into three regular integers
 *
\verbatim embed:rst
This function does the reverse operation of :cpp:func:`lammps_encode_image_flags`
and takes an image flag integer does the bit-shift and bit-masking operations to
decode it and stores the resulting three regular integers into the buffer pointed
to by *flags*.
\endverbatim
 * \param image encoded image flag integer
 * \param [out] flags pointer to storage where the decoded image flags are stored.
 */
void lammps_decode_image_flags(int image, int *flags);
#else
int64_t lammps_encode_image_flags(int, int, int);
/** \brief Decode a single image flag integer into three regular integers
\verbatim embed:rst
This is the interface of the :cpp:func:`lammps_decode_image_flags`
function if LAMMPS has been compiled with the -DLAMMPS_BIGBIG setting.
\endverbatim
 * \param image encoded image flag integer
 * \param [out] flags pointer to storage where the decoded image flags are stored.
 */
void lammps_decode_image_flags(int64_t image, int *flags);
#endif

#ifdef LAMMPS_EXCEPTIONS
int lammps_has_error(void *);
int lammps_get_last_error_message(void *, char *, int);
#endif

#undef LAMMPS
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
