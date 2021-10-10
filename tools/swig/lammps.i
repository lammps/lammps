%include "cpointer.i"
%include "carrays.i"
%include "cdata.i"
%include "cstring.i"

%pointer_functions(int,    int_p);
%pointer_functions(int,    int64_p);
%pointer_functions(double, double_p);

%array_functions(char, char_1d);
%array_functions(int, int_1d);
%array_functions(double, double_1d);

%pointer_cast(void *, int *,     void_p_to_int_p);
%pointer_cast(void *, int **,    void_p_to_int2d_p);
%pointer_cast(void *, int *,     void_p_to_int64_p);
%pointer_cast(void *, int **,    void_p_to_int64_2d_p);
%pointer_cast(void *, double *,  void_p_to_double_p);
%pointer_cast(void *, double **, void_p_to_double_2d_p);

%cstring_output_maxsize(char *buffer, int buf_size);

%{

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
 * Must be kept in sync with the equivalent constants in lammps/constants.py */

enum _LMP_STYLE_CONST {
  LMP_STYLE_GLOBAL=0,           /*!< return global data */
  LMP_STYLE_ATOM  =1,           /*!< return per-atom data */
  LMP_STYLE_LOCAL =2            /*!< return local data */
};

/** Type and size constants for extracting data from computes and fixes.
 *
 * Must be kept in sync with the equivalent constants in lammps/constants.py */

enum _LMP_TYPE_CONST {
  LMP_TYPE_SCALAR=0,            /*!< return scalar */
  LMP_TYPE_VECTOR=1,            /*!< return vector */
  LMP_TYPE_ARRAY =2,            /*!< return array */
  LMP_SIZE_VECTOR=3,            /*!< return length of vector */
  LMP_SIZE_ROWS  =4,            /*!< return number of rows */
  LMP_SIZE_COLS  =5             /*!< return number of columns */
};

/*
 extern void *lammps_open(int argc, char **argv, MPI_Comm comm, void **ptr);
*/
extern void  *lammps_open_no_mpi(int argc, char **argv, void **ptr);
extern void  *lammps_open_fortran(int argc, char **argv, int f_comm);
extern void   lammps_close(void *handle);
extern void   lammps_mpi_init();
extern void   lammps_mpi_finalize();
extern void   lammps_kokkos_finalize();
extern void   lammps_python_finalize();
extern void   lammps_file(void *handle, const char *file);
extern char  *lammps_command(void *handle, const char *cmd);
extern void   lammps_commands_list(void *handle, int ncmd, const char **cmds);
extern void   lammps_commands_string(void *handle, const char *str);
extern double lammps_get_natoms(void *handle);
extern double lammps_get_thermo(void *handle, const char *keyword);
extern void   lammps_extract_box(void *handle, double *boxlo, double *boxhi,
                          double *xy, double *yz, double *xz,
                          int *pflags, int *boxflag);
extern void   lammps_reset_box(void *handle, double *boxlo, double *boxhi,
                        double xy, double yz, double xz);
extern void   lammps_memory_usage(void *handle, double *meminfo);
extern int    lammps_get_mpi_comm(void *handle);
extern int    lammps_extract_setting(void *handle, const char *keyword);
extern int    lammps_extract_global_datatype(void *handle, const char *name);
extern void  *lammps_extract_global(void *handle, const char *name);
extern int    lammps_extract_atom_datatype(void *handle, const char *name);
extern void  *lammps_extract_atom(void *handle, const char *name);
extern void  *lammps_extract_compute(void *handle, char *id, int, int);
extern void  *lammps_extract_fix(void *handle, char *, int, int, int, int);
extern void  *lammps_extract_variable(void *handle, char *, char *);
extern int    lammps_set_variable(void *, char *, char *);
extern void   lammps_gather_atoms(void *, char *, int, int, void *);
extern void   lammps_gather_atoms_concat(void *, char *, int, int, void *);
extern void   lammps_gather_atoms_subset(void *, char *, int, int, int, int *, void *);
extern void   lammps_scatter_atoms(void *, char *, int, int, void *);
extern void   lammps_scatter_atoms_subset(void *, char *, int, int, int, int *, void *);
extern void   lammps_gather_bonds(void *handle, void *data);
extern void   lammps_gather(void *, char *, int, int, void *);
extern void   lammps_gather_concat(void *, char *, int, int, void *);
extern void   lammps_gather_subset(void *, char *, int, int, int, int *, void *);
extern void   lammps_scatter(void *, char *, int, int, void *);
extern void   lammps_scatter_subset(void *, char *, int, int, int, int *, void *);
extern int    lammps_create_atoms(void *handle, int n, int *id, int *type,
                           double *x, double *v, int *image, int bexpand);
/*
 extern int    lammps_create_atoms(void *handle, int n, int64_t *id, int *type, */
extern int    lammps_find_pair_neighlist(void*, char *, int, int, int);
extern int    lammps_find_fix_neighlist(void*, char *, int);
extern int    lammps_find_compute_neighlist(void*, char *, int);
extern int    lammps_neighlist_num_elements(void*, int);
extern void   lammps_neighlist_element_neighbors(void *, int, int, int *, int *, int ** );
extern int    lammps_version(void *handle);
extern void   lammps_get_os_info(char *buffer, int buf_size);
extern int    lammps_config_has_mpi_support();
extern int    lammps_config_has_gzip_support();
extern int    lammps_config_has_png_support();
extern int    lammps_config_has_jpeg_support();
extern int    lammps_config_has_ffmpeg_support();
extern int    lammps_config_has_exceptions();
extern int    lammps_config_has_package(const char *);
extern int    lammps_config_package_count();
extern int    lammps_config_package_name(int, char *, int);
extern int    lammps_config_accelerator(const char *, const char *, const char *);
extern int    lammps_has_gpu_device();
extern void   lammps_get_gpu_device_info(char *buffer, int buf_size);
extern int    lammps_has_style(void *, const char *, const char *);
extern int    lammps_style_count(void *, const char *);
extern int    lammps_style_name(void *, const char *, int, char *buffer, int buf_size);
extern int    lammps_has_id(void *, const char *, const char *);
extern int    lammps_id_count(void *, const char *);
extern int    lammps_id_name(void *, const char *, int, char *buffer, int buf_size);
extern int    lammps_plugin_count();
extern int    lammps_plugin_name(int, char *, char *, int);
/*
 * Have not found a good way to map these functions in a general way.
 * So some individual customization for the specific use case and compilation is needed.
 *
  extern int lammps_encode_image_flags(int ix, int iy, int iz);
  extern void lammps_decode_image_flags(int image, int *flags);
  extern int64_t lammps_encode_image_flags(int ix, int iy, int iz);
  extern void lammps_decode_image_flags(int64_t image, int *flags);

 * Supporting the fix external callback mechanism will require extra code specific to the application.
  typedef void (*FixExternalFnPtr)(void *, int64_t, int, int64_t *, double **, double **);
  extern  void lammps_set_fix_external_callback(void *handle, const char *id, FixExternalFnPtr funcptr, void *ptr);
 * these two functions can only be used from the callback, so we don't support them either
  extern  void lammps_fix_external_set_energy_peratom(void *handle, const char *id, double *eng);
  extern void lammps_fix_external_set_virial_peratom(void *handle, const char *id, double **virial);
*/
extern double **lammps_fix_external_get_force(void *handle, const char *id);
extern void   lammps_fix_external_set_energy_global(void *handle, const char *id, double eng);
extern void   lammps_fix_external_set_virial_global(void *handle, const char *id, double *virial);
extern void   lammps_fix_external_set_vector_length(void *handle, const char *id, int len);
extern void   lammps_fix_external_set_vector(void *handle, const char *id, int idx, double val);

extern void   lammps_free(void *ptr);
extern int    lammps_is_running(void *handle);
extern void   lammps_force_timeout(void *handle);
extern int    lammps_has_error(void *handle);
extern int    lammps_get_last_error_message(void *handle, char *buffer, int buf_size);
%}

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
 * Must be kept in sync with the equivalent constants in lammps/constants.py */

enum _LMP_STYLE_CONST {
  LMP_STYLE_GLOBAL=0,           /*!< return global data */
  LMP_STYLE_ATOM  =1,           /*!< return per-atom data */
  LMP_STYLE_LOCAL =2            /*!< return local data */
};

/** Type and size constants for extracting data from computes and fixes.
 *
 * Must be kept in sync with the equivalent constants in lammps/constants.py */

enum _LMP_TYPE_CONST {
  LMP_TYPE_SCALAR=0,            /*!< return scalar */
  LMP_TYPE_VECTOR=1,            /*!< return vector */
  LMP_TYPE_ARRAY =2,            /*!< return array */
  LMP_SIZE_VECTOR=3,            /*!< return length of vector */
  LMP_SIZE_ROWS  =4,            /*!< return number of rows */
  LMP_SIZE_COLS  =5             /*!< return number of columns */
};

/*
 extern void *lammps_open(int argc, char **argv, MPI_Comm comm, void **ptr);
*/
extern void  *lammps_open_no_mpi(int argc, char **argv, void **ptr);
extern void  *lammps_open_fortran(int argc, char **argv, int f_comm);
extern void   lammps_close(void *handle);
extern void   lammps_mpi_init();
extern void   lammps_mpi_finalize();
extern void   lammps_kokkos_finalize();
extern void   lammps_python_finalize();
extern void   lammps_file(void *handle, const char *file);
extern char  *lammps_command(void *handle, const char *cmd);
extern void   lammps_commands_list(void *handle, int ncmd, const char **cmds);
extern void   lammps_commands_string(void *handle, const char *str);
extern double lammps_get_natoms(void *handle);
extern double lammps_get_thermo(void *handle, const char *keyword);
extern void   lammps_extract_box(void *handle, double *boxlo, double *boxhi,
                          double *xy, double *yz, double *xz,
                          int *pflags, int *boxflag);
extern void   lammps_reset_box(void *handle, double *boxlo, double *boxhi,
                        double xy, double yz, double xz);
extern void   lammps_memory_usage(void *handle, double *meminfo);
extern int    lammps_get_mpi_comm(void *handle);
extern int    lammps_extract_setting(void *handle, const char *keyword);
extern int    lammps_extract_global_datatype(void *handle, const char *name);
extern void  *lammps_extract_global(void *handle, const char *name);
extern int    lammps_extract_atom_datatype(void *handle, const char *name);
extern void  *lammps_extract_atom(void *handle, const char *name);
extern void  *lammps_extract_compute(void *handle, char *id, int, int);
extern void  *lammps_extract_fix(void *handle, char *, int, int, int, int);
extern void  *lammps_extract_variable(void *handle, char *, char *);
extern int    lammps_set_variable(void *, char *, char *);
extern void   lammps_gather_atoms(void *, char *, int, int, void *);
extern void   lammps_gather_atoms_concat(void *, char *, int, int, void *);
extern void   lammps_gather_atoms_subset(void *, char *, int, int, int, int *, void *);
extern void   lammps_scatter_atoms(void *, char *, int, int, void *);
extern void   lammps_scatter_atoms_subset(void *, char *, int, int, int, int *, void *);
extern void   lammps_gather_bonds(void *handle, void *data);
extern void   lammps_gather(void *, char *, int, int, void *);
extern void   lammps_gather_concat(void *, char *, int, int, void *);
extern void   lammps_gather_subset(void *, char *, int, int, int, int *, void *);
extern void   lammps_scatter(void *, char *, int, int, void *);
extern void   lammps_scatter_subset(void *, char *, int, int, int, int *, void *);
extern int    lammps_create_atoms(void *handle, int n, int *id, int *type,
                           double *x, double *v, int *image, int bexpand);
/*
 extern int    lammps_create_atoms(void *handle, int n, int64_t *id, int *type, */
extern int    lammps_find_pair_neighlist(void*, char *, int, int, int);
extern int    lammps_find_fix_neighlist(void*, char *, int);
extern int    lammps_find_compute_neighlist(void*, char *, int);
extern int    lammps_neighlist_num_elements(void*, int);
extern void   lammps_neighlist_element_neighbors(void *, int, int, int *, int *, int ** );
extern int    lammps_version(void *handle);
extern void   lammps_get_os_info(char *buffer, int buf_size);
extern int    lammps_config_has_mpi_support();
extern int    lammps_config_has_gzip_support();
extern int    lammps_config_has_png_support();
extern int    lammps_config_has_jpeg_support();
extern int    lammps_config_has_ffmpeg_support();
extern int    lammps_config_has_exceptions();
extern int    lammps_config_has_package(const char *);
extern int    lammps_config_package_count();
extern int    lammps_config_package_name(int, char *, int);
extern int    lammps_config_accelerator(const char *, const char *, const char *);
extern int    lammps_has_gpu_device();
extern void   lammps_get_gpu_device_info(char *buffer, int buf_size);
extern int    lammps_has_style(void *, const char *, const char *);
extern int    lammps_style_count(void *, const char *);
extern int    lammps_style_name(void *, const char *, int, char *buffer, int buf_size);
extern int    lammps_has_id(void *, const char *, const char *);
extern int    lammps_id_count(void *, const char *);
extern int    lammps_id_name(void *, const char *, int, char *buffer, int buf_size);
extern int    lammps_plugin_count();
extern int    lammps_plugin_name(int, char *, char *, int);
/*
extern int lammps_encode_image_flags(int ix, int iy, int iz);
extern void lammps_decode_image_flags(int image, int *flags);
extern int64_t lammps_encode_image_flags(int ix, int iy, int iz);
extern void lammps_decode_image_flags(int64_t image, int *flags);
typedef void (*FixExternalFnPtr)(void *, int64_t, int, int64_t *, double **, double **);
extern  void lammps_set_fix_external_callback(void *handle, const char *id, FixExternalFnPtr funcptr, void *ptr);
extern  void lammps_fix_external_set_energy_peratom(void *handle, const char *id, double *eng);
extern void lammps_fix_external_set_virial_peratom(void *handle, const char *id, double **virial);
*/
extern double **lammps_fix_external_get_force(void *handle, const char *id);
extern void   lammps_fix_external_set_energy_global(void *handle, const char *id, double eng);
extern void   lammps_fix_external_set_virial_global(void *handle, const char *id, double *virial);
extern void   lammps_fix_external_set_vector_length(void *handle, const char *id, int len);
extern void   lammps_fix_external_set_vector(void *handle, const char *id, int idx, double val);

extern void   lammps_free(void *ptr);
extern int    lammps_is_running(void *handle);
extern void   lammps_force_timeout(void *handle);
extern int    lammps_has_error(void *handle);
extern int    lammps_get_last_error_message(void *handle, char *buffer, int buf_size);

/* last revised on 21 July 2021 */
